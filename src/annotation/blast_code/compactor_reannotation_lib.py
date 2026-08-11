#!/usr/bin/env python3
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from shutil import which

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    limit = sys.maxsize
    while True:
        limit = int(limit / 10)
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            continue

DEFAULT_INPUT_DIR = "/scratch/users/jiamuyu/proj_botryllus/splash2/260713_01_3ants_challenge"
DEFAULT_OUTPUT_DIR = (
    "/scratch/users/jiamuyu/proj_botryllus/flash/results/"
    "260714-00-3ants-challenge/filter1/noCluster/hyena/normalized"
)
DEFAULT_RESULTS_DIR = "/scratch/users/jiamuyu/proj_botryllus/flash/results/260714-00-3ants-challenge"
DEFAULT_PLOT_PREFIX = (
    DEFAULT_OUTPUT_DIR
    + "/260714-00-3ants-challenge_hyena_adelie_results_top2000_target1_k41_s41_trainProp0.8"
    + "_nonzero_coefficients"
)
DEFAULT_BLASTP_ANNOTATED = DEFAULT_PLOT_PREFIX + "_blastp_annotated.tsv"
DEFAULT_BLAST_ANNOTATED = DEFAULT_PLOT_PREFIX + "_blast_annotated.tsv"
DEFAULT_PLOT_PDF = DEFAULT_PLOT_PREFIX + "_blast_annotated_plots.pdf"
DEFAULT_PLOT_SUMMARY = DEFAULT_PLOT_PREFIX + "_blast_annotated_plots_summary.tsv"
DEFAULT_CLUSTERS = (
    DEFAULT_RESULTS_DIR
    + "/filter1/noCluster/260714-00-3ants-challenge_sequences_per_cluster_top2000-clusters_target1_k41_s41.tsv"
)
DEFAULT_FEATHER = (
    DEFAULT_RESULTS_DIR
    + "/260714-00-3ants-challenge_hyena_top_variance_features_for_glmnet_filter1_noCluster_top2000_target1_k41_s41_normalized.feather"
)
DEFAULT_SAMPLE_SEQS = (
    DEFAULT_RESULTS_DIR
    + "/260714-00-3ants-challenge_prepared_sequences_filter1_noCluster_top2000_target1_k41_s41_sample_sequences.tsv"
)
DEFAULT_METADATA = "/scratch/users/jiamuyu/proj_botryllus/splash2/260713_01_3ants_challenge/metadata.tsv"
DEFAULT_TAXIDS = "300111;102681;104421"
DEFAULT_BLAST_DB = "/scratch/users/jiamuyu/dabs_ref/blast/"
DEFAULT_PROTEIN_DB = "refseq_protein"
DEFAULT_OUTPUT_STEM = "resfungi_compactors"
REPO_ROOT = Path(__file__).resolve().parent


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Pick representative resfungi compactors, BLAST/BLASTP them, and annotate "
            "seed extendors by their anchor-matched compactors."
        )
    )
    parser.add_argument("--input_dir", default=DEFAULT_INPUT_DIR)
    parser.add_argument(
        "--compactor_table",
        default=None,
        help=(
            "Optional TSV with anchor/compactor rows to blast directly. When omitted, "
            "the legacy mode scans INPUT_DIR/*resfungi.tsv and picks one representative "
            "compactor per file."
        ),
    )
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--output_stem",
        default=DEFAULT_OUTPUT_STEM,
        help=(
            "Stem for FASTA, BLAST, and selected-compactor outputs. The default keeps "
            "legacy filenames such as resfungi_compactors_blast.tsv."
        ),
    )
    parser.add_argument(
        "--compactor_source_stems",
        default=None,
        help=(
            "Comma/semicolon-separated compactor output stems to combine in --noblastp "
            "plot-rescue mode. Defaults to --output_stem."
        ),
    )
    parser.add_argument("--seeds", default=None, help="Defaults to INPUT_DIR/seeds.resfungi.raw")
    parser.add_argument("--taxids", default=DEFAULT_TAXIDS)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--translation_table", default="1")
    parser.add_argument("--entrez_email", default=os.environ.get("ENTREZ_EMAIL", "v8514616@outlook.com"))
    parser.add_argument("--blast_db", default=DEFAULT_BLAST_DB)
    parser.add_argument(
        "--protein_db",
        default=DEFAULT_PROTEIN_DB,
        help=(
            "BLAST-formatted protein database name for the BLASTX/BLASTP helper. "
            "Examples: refseq_protein, swissprot, clusteredNR."
        ),
    )
    parser.add_argument("--anchor_len", type=int, default=31)
    parser.add_argument("--thresholds", default="1000,100,5")
    parser.add_argument("--skip_blast", action="store_true")
    parser.add_argument(
        "--reblast_mode",
        choices=["missing", "none"],
        default="missing",
        help="Run unrestricted reblast for restricted-taxid misses, or skip it.",
    )
    parser.add_argument(
        "--blast_modes",
        choices=["both", "blastn", "blastp"],
        default="both",
        help="Which BLAST searches to run or reuse.",
    )
    parser.add_argument(
        "--overwrite_blast",
        action="store_true",
        help="Re-run BLAST/BLASTP even if the expected output TSVs already exist.",
    )
    parser.add_argument(
        "--noblastp",
        action="store_true",
        help=(
            "Skip compactor selection/BLAST and only reuse existing compactor annotations "
            "to fill NO BLAST/NO MATCH rows in the FLASH blast plot outputs."
        ),
    )
    parser.add_argument("--plot_blastp_annotated", default=DEFAULT_BLASTP_ANNOTATED)
    parser.add_argument("--plot_blast_annotated", default=DEFAULT_BLAST_ANNOTATED)
    parser.add_argument("--plot_pdf", default=DEFAULT_PLOT_PDF)
    parser.add_argument("--plot_summary", default=DEFAULT_PLOT_SUMMARY)
    parser.add_argument("--plot_clusters", default=DEFAULT_CLUSTERS)
    parser.add_argument("--plot_feather", default=DEFAULT_FEATHER)
    parser.add_argument("--plot_sample_seqs", default=DEFAULT_SAMPLE_SEQS)
    parser.add_argument("--plot_metadata", default=DEFAULT_METADATA)
    parser.add_argument("--plot_target_vars", default="")
    parser.add_argument("--plot_confound_vars", default="")
    parser.add_argument(
        "--fungus_output",
        default=str(REPO_ROOT / "fungus.tsv"),
        help=(
            "Output path for rows from the generated compactor plot summary containing "
            "fungus_species_residual. Defaults to fungus.tsv in the repository root."
        ),
    )
    parser.add_argument(
        "--plot_rscript",
        default=os.environ.get("RESFUNGI_PLOT_RSCRIPT", "Rscript"),
        help="Rscript executable for regenerating the compactor plot PDF.",
    )
    return parser.parse_args()


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_id(value):
    value = re.sub(r"\.resfungi\.tsv$", "", Path(value).name)
    value = re.sub(r"\.tsv$", "", value)
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_") or "resfungi"


def choose_compactor_from_rows(rows, thresholds):
    if not rows:
        return None

    best_by_length = {}
    for row in rows:
        current = best_by_length.get(row["length"])
        if current is None:
            best_by_length[row["length"]] = row
        elif (
            row["exact_support"] > current["exact_support"]
            or (
                row["exact_support"] == current["exact_support"]
                and row["row_index"] < current["row_index"]
            )
        ):
            best_by_length[row["length"]] = row

    representatives = sorted(best_by_length.values(), key=lambda item: item["row_index"])
    for threshold in thresholds:
        candidates = [row for row in representatives if row["exact_support"] > threshold]
        if len(candidates) == 1:
            chosen = candidates[0]
            chosen["support_threshold"] = threshold
            chosen["selection_reason"] = "only_representative_above_threshold"
            return chosen
        for prev, current in zip(candidates, candidates[1:]):
            if current["exact_support"] > prev["exact_support"]:
                current["support_threshold"] = threshold
                current["selection_reason"] = "support_increase_vs_previous_representative"
                return current
    return None


def choose_compactor(path, thresholds):
    rows = []
    for idx, row in enumerate(read_tsv(path), start=1):
        compactor = (row.get("compactor") or "").strip()
        anchor = (row.get("anchor") or "").strip()
        exact_support = numeric(row.get("exact_support"))
        if not compactor or not anchor or exact_support is None:
            continue
        rows.append(
            {
                "row_index": idx,
                "anchor": anchor,
                "compactor": compactor,
                "length": len(compactor),
                "exact_support": exact_support,
                "source_file": str(path),
            }
        )
    return choose_compactor_from_rows(rows, thresholds)


def read_compactor_table(path, thresholds):
    rows_by_anchor = defaultdict(list)
    records = []
    for idx, row in enumerate(read_tsv(path), start=1):
        compactor = (row.get("compactor") or "").strip()
        anchor = (row.get("anchor") or "").strip()
        exact_support = numeric(row.get("exact_support"))
        if not compactor or not anchor or exact_support is None:
            continue
        length_value = int(numeric(row.get("total_length"), len(compactor)) or len(compactor))
        rows_by_anchor[anchor].append(
            {
                "anchor": anchor,
                "compactor": compactor,
                "length": length_value,
                "exact_support": exact_support,
                "row_index": idx,
                "row_id": row.get("id", idx),
                "source_file": str(path),
                "support": row.get("support", "NA"),
                "expected_read_count": row.get("expected_read_count", "NA"),
                "extender_specificity": row.get("extender_specificity", "NA"),
                "num_extended": row.get("num_extended", "NA"),
            }
        )

    skipped = 0
    for anchor, rows in sorted(rows_by_anchor.items(), key=lambda item: min(row["row_index"] for row in item[1])):
        chosen = choose_compactor_from_rows(rows, thresholds)
        if chosen is None:
            skipped += 1
            continue
        row_id = chosen.get("row_id")
        row_id = row_id if has_text(row_id) else chosen["row_index"]
        chosen = dict(chosen)
        chosen["query"] = (
            f"{safe_id(path)}__anchor{len(records) + 1}__row{chosen['row_index']}"
            f"__id{row_id}__len{chosen['length']}"
        )
        chosen["selection_reason"] = f"per_anchor_{chosen['selection_reason']}"
        records.append(chosen)

    print(
        f"Selected {len(records)}/{len(rows_by_anchor)} anchor-level compactors from {path}; "
        f"skipped {skipped} anchors with no representative passing thresholds {thresholds}."
    )
    return records


def output_paths(output_dir, stem):
    return {
        "fasta": output_dir / f"{stem}.fasta",
        "selected": output_dir / f"{stem}_selected.tsv",
        "blast": output_dir / f"{stem}_blast.tsv",
        "reblast": output_dir / f"{stem}_reblast.tsv",
        "blastp": output_dir / f"{stem}_blastp.tsv",
        "reblastp": output_dir / f"{stem}_reblastp.tsv",
        "tmp": output_dir / f"{stem}_tmp",
        "compactor_annotations": compactor_annotations_path(output_dir, stem),
        "seed_annotations": seed_annotations_path(output_dir, stem),
    }


def compactor_annotations_path(output_dir, stem):
    if stem == DEFAULT_OUTPUT_STEM:
        return output_dir / "resfungi_compactor_annotations.tsv"
    return output_dir / f"{stem}_annotations.tsv"


def seed_annotations_path(output_dir, stem):
    if stem == DEFAULT_OUTPUT_STEM:
        return output_dir / "seeds_resfungi_compactor_annotations.tsv"
    return output_dir / f"seeds_{stem}_annotations.tsv"


def write_fasta(records, output_fasta):
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w", newline="") as handle:
        for record in records:
            header = (
                f"{record['query']} anchor={record['anchor']} length={record['length']}"
                f" exact_support={record['exact_support']:g} row={record['row_index']}"
                f" source={Path(record['source_file']).name}"
            )
            handle.write(f">{header}\n{record['compactor']}\n")


def write_selected_compactors(records, output_path):
    columns = [
        "query",
        "anchor",
        "compactor",
        "length",
        "exact_support",
        "row_index",
        "support_threshold",
        "selection_reason",
        "source_file",
        "support",
        "expected_read_count",
        "extender_specificity",
        "num_extended",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for record in records:
            writer.writerow({col: record.get(col, "NA") for col in columns})


def read_selected_compactors(path):
    records = []
    if not path.exists() or path.stat().st_size == 0:
        return records
    for row in read_tsv(path):
        record = dict(row)
        record["length"] = int(numeric(record.get("length"), 0) or 0)
        record["exact_support"] = numeric(record.get("exact_support"), "NA")
        record["row_index"] = int(numeric(record.get("row_index"), 0) or 0)
        records.append(record)
    return records


def make_python_shim(output_dir):
    shim_dir = output_dir / "resfungi_compactors_python_shim"
    shim_dir.mkdir(parents=True, exist_ok=True)
    python_link = shim_dir / "python"
    python3_link = shim_dir / "python3"
    for link in (python_link, python3_link):
        if link.exists() or link.is_symlink():
            link.unlink()
        try:
            link.symlink_to(sys.executable)
        except OSError:
            # Some filesystems disallow symlinks; a tiny shell shim is enough.
            link.write_text(f"#!/bin/sh\nexec {sys.executable} \"$@\"\n")
            link.chmod(0o755)
    env = os.environ.copy()
    env["PATH"] = f"{shim_dir}{os.pathsep}{env.get('PATH', '')}"
    return env


def run_command(command, env=None, cwd=None):
    print("+ " + " ".join(map(str, command)), flush=True)
    subprocess.run(command, check=True, env=env, cwd=cwd)


def require_output(path, label):
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(
            f"{label} did not create a non-empty output file: {path}. "
            "Check the module-loaded python/Rscript versions in the job log."
        )


def has_output(path):
    return path.exists() and path.stat().st_size > 0


def write_empty_blastn(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "query\tsubject\tidentity\talignment_length\tmismatches\tgap_opens\tq_start\tq_end\t"
        "s_start\ts_end\tsstrand\tevalue\tqcovs\tsgi\tsacc\tslen\tstaxids\tstitle\tspecies_origin\tfeatures\tfeatures_10000_window\n"
    )


def write_empty_blastp(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "query\tidentity\tevalue\tqcovs\tqframe\tstaxids\tstitle\t"
        "NCBI_protein_accession\tUniProt_accession\tmethod\tGO\n"
    )


def count_tsv_rows(path):
    if not path.exists() or path.stat().st_size == 0:
        return 0
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def run_blasts(args, output_fasta, output_dir):
    paths = output_paths(output_dir, args.output_stem)
    blast = paths["blast"]
    reblast = paths["reblast"]
    blastp = paths["blastp"]
    reblastp = paths["reblastp"]
    temp_root = paths["tmp"]
    command_env = make_python_shim(temp_root)
    print(f"Wrapper subprocess python: {which('python', path=command_env['PATH'])}", flush=True)
    run_command(["python", "--version"], env=command_env)

    should_run_blastn = args.blast_modes in {"both", "blastn"}
    should_run_blastp = args.blast_modes in {"both", "blastp"}
    reblast_arg = str(reblast) if args.reblast_mode == "missing" else ""
    reblastp_arg = str(reblastp) if args.reblast_mode == "missing" else ""

    need_blastn = args.overwrite_blast or not has_output(blast) or (
        args.reblast_mode == "missing" and not has_output(reblast)
    )
    need_blastp = args.overwrite_blast or not has_output(blastp) or (
        args.reblast_mode == "missing" and not has_output(reblastp)
    )

    if should_run_blastn and need_blastn:
        run_command(
            [
                "bash",
                str(REPO_ROOT / "src/annotation/blast_code/run_blast.sh"),
                str(output_fasta),
                str(temp_root / "split_fasta"),
                str(temp_root / "blast"),
                str(blast),
                str(args.threads),
                args.taxids,
                args.entrez_email,
                str(temp_root),
                args.blast_db,
                "0",
                "all",
                "",
                "10",
                "",
                "0",
                reblast_arg,
            ],
            env=command_env,
            cwd=REPO_ROOT,
        )
    else:
        print(f"Reusing existing restricted BLASTN output: {blast}")
    if should_run_blastn or has_output(blast):
        require_output(blast, "Restricted BLASTN")
    else:
        write_empty_blastn(blast)
    if args.reblast_mode == "missing":
        if should_run_blastn:
            require_output(reblast, "Unrestricted reBLASTN")
        elif not has_output(reblast):
            write_empty_blastn(reblast)
    else:
        write_empty_blastn(reblast)

    if should_run_blastp and need_blastp:
        run_command(
            [
                "bash",
                str(REPO_ROOT / "src/annotation/blast_code/run_blastp.sh"),
                str(output_fasta),
                str(temp_root / "split_fasta_blastp"),
                str(temp_root / "blastp"),
                str(blastp),
                str(args.threads),
                args.taxids,
                str(args.translation_table),
                args.blast_db,
                args.protein_db,
                "0",
                "all",
                "",
                "10",
                "",
                "0",
                reblastp_arg,
            ],
            env=command_env,
            cwd=REPO_ROOT,
        )
    else:
        print(f"Reusing existing restricted BLASTP output: {blastp}")
    if should_run_blastp or has_output(blastp):
        require_output(blastp, "Restricted BLASTP")
    else:
        write_empty_blastp(blastp)
    if args.reblast_mode == "missing":
        if should_run_blastp:
            require_output(reblastp, "Unrestricted reBLASTP")
        elif not has_output(reblastp):
            write_empty_blastp(reblastp)
    else:
        write_empty_blastp(reblastp)

    print(
        "BLAST row counts: "
        f"restricted_blastn={count_tsv_rows(blast)}, "
        f"unrestricted_blastn={count_tsv_rows(reblast)}, "
        f"restricted_blastp={count_tsv_rows(blastp)}, "
        f"unrestricted_blastp={count_tsv_rows(reblastp)}"
    )
    return blast, reblast, blastp, reblastp


def has_text(value):
    value = "" if value is None else str(value).strip()
    return value not in {"", "NA", "NaN", "None", "none", "[]"}


def clean_label(value):
    value = "" if value is None else str(value)
    value = re.sub(r"LOC\d+[- ]*", "", value)
    value = re.sub(r"\s+isoform\s+X\d+\b", "", value)
    value = re.sub(r"\s+transcript\s+variant\s+X?\d+\b", "", value)
    value = re.sub(r"\s+variant\s+X?\d+\b", "", value)
    value = re.sub(r"\s+", " ", value).strip(" ,;")
    return value or None


def extract_qualifier(features, qualifier):
    if not has_text(features):
        return []
    pattern = rf"['\"]{qualifier}['\"]:\s*(?:\[([^\]]*)\]|([^,}}\]]+))"
    values = []
    for match in re.finditer(pattern, str(features)):
        text = match.group(1) or match.group(2) or ""
        values.extend(re.findall(r"(?<=[\'\"])[^\'\"]+(?=[\'\"])", text))
    cleaned = []
    for value in values:
        label = clean_label(value)
        if label and label not in {"None", "NA"} and label not in cleaned:
            cleaned.append(label)
    return cleaned


def blastn_label(row):
    feature_text = ";".join(
        str(row.get(col, ""))
        for col in row
        if col == "features" or col.startswith("features_")
    )
    products = extract_qualifier(feature_text, "product")
    genes = extract_qualifier(feature_text, "gene")
    labels = products or genes
    if labels:
        return ";".join(labels)
    if has_text(row.get("identity")) or has_text(row.get("qcovs")):
        return "UNANNOTATED"
    return None


def blastp_label(row):
    title = clean_label(re.sub(r"\s*\[[^\]]+\]\s*$", "", str(row.get("stitle", ""))))
    if title:
        return title
    for col in ("NCBI_protein_accession", "UniProt_accession", "GO"):
        if has_text(row.get(col)):
            return str(row.get(col)).strip()
    if has_text(row.get("identity")) or has_text(row.get("qcovs")):
        return "UNANNOTATED"
    return None


def species_from_stitle(value):
    value = "" if value is None else str(value)
    matches = re.findall(r"\[([^\]]+)\]", value)
    if not matches:
        return "NA"
    species = clean_label(matches[-1])
    return species if species else "NA"


def raw_blastn_annotation(row):
    values = []
    for col in ("stitle", "features", "features_10000_window"):
        if has_text(row.get(col)):
            values.append(f"{col}={str(row.get(col)).replace(chr(9), ' ')}")
    return " | ".join(values) if values else "NA"


def raw_blastp_annotation(row):
    values = []
    for col in ("stitle", "NCBI_protein_accession", "UniProt_accession", "GO"):
        if has_text(row.get(col)):
            values.append(f"{col}={str(row.get(col)).replace(chr(9), ' ')}")
    return " | ".join(values) if values else "NA"


def normalize_query_id(query):
    if query is None:
        return ""
    query = str(query).strip()
    if not query:
        return ""
    query = query.split()[0]
    return query.split("|")[0]


def read_annotation_table(path, mode, source):
    annotations = defaultdict(list)
    if not path.exists() or path.stat().st_size == 0:
        return annotations
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "query" not in reader.fieldnames:
            return annotations
        for row in reader:
            query = normalize_query_id(row.get("query"))
            if not query:
                continue
            label = blastp_label(row) if mode == "blastp" else blastn_label(row)
            raw_annotation = raw_blastp_annotation(row) if mode == "blastp" else raw_blastn_annotation(row)
            if not label:
                label = "UNANNOTATED" if raw_annotation != "NA" else "NO MATCH"
            identity = row.get("identity", "NA")
            qcovs = row.get("qcovs", "NA")
            species_origin = species_from_stitle(row.get("stitle"))
            staxids = row.get("staxids", "NA")
            annotations[query].append(
                {
                    "annotation_label": label,
                    "annotation_source": (
                        "restricted_taxid" if source == "restricted" else "outside_taxid"
                    ),
                    "blast_mode": mode,
                    "blast_scope": source,
                    "identity": identity,
                    "qcovs": qcovs,
                    "species_origin": species_origin,
                    "staxids": staxids,
                    "raw_annotation": raw_annotation,
                    f"{source}_{mode}_label": label,
                    f"{source}_{mode}_identity": identity,
                    f"{source}_{mode}_qcovs": qcovs,
                    f"{source}_{mode}_species": species_origin,
                    f"{source}_{mode}_staxids": staxids,
                }
            )
    return annotations


def merge_annotation_maps(*maps):
    merged = defaultdict(list)
    for annotation_map in maps:
        for query, entries in annotation_map.items():
            merged[query].extend(entries)
    return dict(merged)


def no_hit_entry():
    return {
        "annotation_label": "NO MATCH",
        "annotation_source": "no_hit",
        "blast_mode": "NA",
        "blast_scope": "NA",
        "identity": "NA",
        "qcovs": "NA",
        "species_origin": "NA",
        "staxids": "NA",
        "raw_annotation": "NA",
    }


def write_compactor_annotation_summary(records, annotations, output_path):
    columns = [
        "query",
        "anchor",
        "compactor",
        "length",
        "exact_support",
        "source_file",
        "support",
        "expected_read_count",
        "extender_specificity",
        "num_extended",
        "annotation_label",
        "annotation_source",
        "blast_mode",
        "blast_scope",
        "identity",
        "qcovs",
        "species_origin",
        "staxids",
        "raw_annotation",
        "restricted_blastp_label",
        "restricted_blast_label",
        "unrestricted_blastp_label",
        "unrestricted_blast_label",
        "restricted_blastp_species",
        "restricted_blast_species",
        "unrestricted_blastp_species",
        "unrestricted_blast_species",
        "restricted_blastp_staxids",
        "restricted_blast_staxids",
        "unrestricted_blastp_staxids",
        "unrestricted_blast_staxids",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_counts = defaultdict(int)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for record in records:
            entries = annotations.get(record["query"], []) or [no_hit_entry()]
            for entry in entries:
                row = {**record, **entry}
                source_counts[row["annotation_source"]] += 1
                writer.writerow({col: row.get(col, "NA") for col in columns})
    print(
        "Compactor annotation source counts: "
        + ", ".join(f"{key}={source_counts[key]}" for key in sorted(source_counts))
    )


def read_seed_rows(seeds_path, anchor_len):
    rows = []
    if not seeds_path.exists():
        return rows
    with open(seeds_path, newline="") as handle:
        for idx, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            extendor = fields[3] if len(fields) >= 4 else ""
            rows.append(
                {
                    "seed_row": idx,
                    "seed_extendor": extendor,
                    "seed_anchor": extendor[-anchor_len:] if len(extendor) >= anchor_len else extendor,
                    "seed_source": str(seeds_path),
                    "raw_seed_row": line.replace("\t", "\\t"),
                }
            )
    return rows


def clean_sequence_candidate(value):
    value = "" if value is None else str(value).strip()
    value = re.sub(r"^cluster_\d+_", "", value)
    value = value.replace("-", "")
    value = re.sub(r"[^ACGTNacgtn]", "", value)
    return value.upper()


def read_plot_summary_seed_rows(summary_path, anchor_len):
    rows = []
    summary_path = Path(summary_path)
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        return rows
    preferred_columns = [
        "seed_extendor",
        "sequence",
        "query",
        "extendor",
        "anchor_target",
        "anchor_target_sequence",
    ]
    with open(summary_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, row in enumerate(reader, start=1):
            extendor = ""
            for col in preferred_columns:
                candidate = clean_sequence_candidate(row.get(col))
                if len(candidate) >= anchor_len:
                    extendor = candidate
                    break
            if not extendor:
                for value in row.values():
                    candidate = clean_sequence_candidate(value)
                    if len(candidate) >= anchor_len:
                        extendor = candidate
                        break
            if not extendor:
                continue
            rows.append(
                {
                    "seed_row": idx,
                    "seed_extendor": extendor,
                    "seed_anchor": extendor[-anchor_len:],
                    "seed_source": str(summary_path),
                    "raw_seed_row": "\t".join(str(row.get(col, "")) for col in (reader.fieldnames or [])).replace(
                        "\t", "\\t"
                    ),
                }
            )
    return rows


def combine_seed_rows(*seed_groups):
    combined = []
    seen = set()
    for rows in seed_groups:
        for row in rows:
            key = row.get("seed_extendor")
            if not key or key in seen:
                continue
            seen.add(key)
            row = dict(row)
            row["seed_row"] = len(combined) + 1
            combined.append(row)
    return combined


def records_by_anchor(records):
    by_anchor = defaultdict(list)
    for record in records:
        by_anchor[record["anchor"]].append(record)
    return by_anchor


def write_seed_annotation_summary(seeds, records_by_anchor, annotations, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "seed_row",
        "seed_anchor",
        "seed_extendor",
        "compactor_query",
        "compactor",
        "compactor_length",
        "compactor_exact_support",
        "compactor_anchor",
        "compactor_source_file",
        "annotation_label",
        "annotation_source",
        "blast_mode",
        "blast_scope",
        "identity",
        "qcovs",
        "species_origin",
        "staxids",
        "raw_annotation",
        "restricted_blastp_label",
        "restricted_blast_label",
        "unrestricted_blastp_label",
        "unrestricted_blast_label",
        "restricted_blastp_species",
        "restricted_blast_species",
        "unrestricted_blastp_species",
        "unrestricted_blast_species",
        "restricted_blastp_staxids",
        "restricted_blast_staxids",
        "unrestricted_blastp_staxids",
        "unrestricted_blast_staxids",
        "seed_source",
        "raw_seed_row",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for seed in seeds:
            matched_records = records_by_anchor.get(seed["seed_anchor"], [])
            if not matched_records:
                matched_records = [{}]
            for record in matched_records:
                query = record.get("query", "")
                row = {**seed}
                if record:
                    row.update(
                        {
                            "compactor_query": query,
                            "compactor": record["compactor"],
                            "compactor_length": record["length"],
                            "compactor_exact_support": record["exact_support"],
                            "compactor_anchor": record["anchor"],
                            "compactor_source_file": record["source_file"],
                        }
                    )
                    entries = annotations.get(query, []) or [no_hit_entry()]
                else:
                    entries = [
                        {
                            **no_hit_entry(),
                            "annotation_source": "no_anchor_matched_compactor",
                        }
                    ]
                for entry in entries:
                    out_row = {**row, **entry}
                    writer.writerow({col: out_row.get(col, "NA") for col in columns})


UNRESOLVED_LABELS = {
    "",
    "NA",
    "NAN",
    "NONE",
    "NO MATCH",
    "NO BLAST",
    "UNANNOTATED",
    "NO PROTEIN/GENE HIT",
    "BLAST",
    "BLASTP",
}


def is_real_annotation(label):
    label = "" if label is None else str(label).strip()
    return label.upper() not in UNRESOLVED_LABELS


def sequence_from_plot_query(query):
    query = normalize_query_id(query)
    return re.sub(r"^cluster_\d+_", "", query)


def tsv_output_path(path, suffix):
    path = Path(path)
    if path.name.endswith(".tsv"):
        return path.with_name(path.name[:-4] + suffix + ".tsv")
    return path.with_name(path.name + suffix + ".tsv")


def pdf_output_path(path, suffix):
    path = Path(path)
    if path.name.endswith(".pdf"):
        return path.with_name(path.name[:-4] + suffix + ".pdf")
    return path.with_name(path.name + suffix + ".pdf")


def read_seed_compactor_annotation_map(path):
    annotation_map = {}
    if not path.exists() or path.stat().st_size == 0:
        return annotation_map
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            extendor = (row.get("seed_extendor") or "").strip()
            label = clean_label(row.get("annotation_label")) or ""
            if not extendor or not is_real_annotation(label):
                continue
            candidate = {
                "label": f"{label} (COMPACTOR)",
                "identity": row.get("identity", "NA"),
                "qcovs": row.get("qcovs", "NA"),
                "raw_annotation": row.get("raw_annotation", "NA"),
                "species_origin": row.get("species_origin", "NA"),
                "staxids": row.get("staxids", "NA"),
                "compactor_query": row.get("compactor_query", "NA"),
                "compactor": row.get("compactor", "NA"),
                "compactor_length": row.get("compactor_length", "NA"),
                "compactor_exact_support": row.get("compactor_exact_support", "NA"),
                "annotation_source": row.get("annotation_source", "NA"),
                "blast_mode": row.get("blast_mode", "NA"),
            }
            current = annotation_map.get(extendor)
            if current is None:
                annotation_map[extendor] = candidate
                continue
            current_priority = 0 if current.get("annotation_source") == "restricted_taxid" else 1
            candidate_priority = 0 if candidate.get("annotation_source") == "restricted_taxid" else 1
            if candidate_priority < current_priority:
                annotation_map[extendor] = candidate
    return annotation_map


def annotation_priority(hit):
    source_rank = 0 if hit.get("annotation_source") == "restricted_taxid" else 1
    mode_rank = 0 if hit.get("blast_mode") == "blastp" else 1
    label_rank = len(str(hit.get("label", "")))
    return (source_rank, mode_rank, label_rank)


def compactor_hit_from_row(row):
    label = clean_label(row.get("annotation_label")) or ""
    if not is_real_annotation(label):
        return None
    return {
        "label": f"{label} (COMPACTOR)",
        "identity": row.get("identity", "NA"),
        "qcovs": row.get("qcovs", "NA"),
        "raw_annotation": row.get("raw_annotation", "NA"),
        "species_origin": row.get("species_origin", "NA"),
        "staxids": row.get("staxids", "NA"),
        "compactor_query": row.get("query") or row.get("compactor_query", "NA"),
        "compactor": row.get("compactor", "NA"),
        "compactor_length": row.get("length") or row.get("compactor_length", "NA"),
        "compactor_exact_support": row.get("exact_support") or row.get("compactor_exact_support", "NA"),
        "annotation_source": row.get("annotation_source", "NA"),
        "blast_mode": row.get("blast_mode", "NA"),
        "compactor_anchor": row.get("anchor") or row.get("compactor_anchor", "NA"),
    }


def read_compactor_anchor_annotation_map(path):
    anchor_map = defaultdict(list)
    if not path.exists() or path.stat().st_size == 0:
        return anchor_map
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            anchor = clean_sequence_candidate(row.get("anchor") or row.get("compactor_anchor"))
            hit = compactor_hit_from_row(row)
            if not anchor or hit is None:
                continue
            anchor_map[anchor].append(hit)
    return dict(anchor_map)


def select_compactor_hit(candidates):
    candidates = [hit for hit in candidates if hit and is_real_annotation(hit.get("label", "").replace("(COMPACTOR)", ""))]
    if not candidates:
        return None
    return sorted(candidates, key=annotation_priority)[0]


def merge_compactor_maps(target, source):
    for sequence, hit in source.items():
        current = target.get(sequence)
        if current is None or annotation_priority(hit) < annotation_priority(current):
            target[sequence] = hit


def add_anchor_hits(anchor_map, source):
    for anchor, hits in source.items():
        anchor_map.setdefault(anchor, []).extend(hits)


def parse_stem_list(value, default_stem):
    value = value or default_stem
    stems = [item.strip() for item in re.split(r"[;,]", value) if item.strip()]
    return stems or [default_stem]


def write_fungus_subset(summary_path, output_path):
    summary_path = Path(summary_path)
    output_path = Path(output_path)
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        print(f"Skipping fungus subset because generated summary does not exist: {summary_path}")
        return

    kept = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, newline="") as in_handle, open(output_path, "w", newline="") as out_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return
        writer.writerow(header)
        for row in reader:
            if "fungus_species_residual" in "\t".join(row):
                writer.writerow(row)
                kept += 1
    print(f"Wrote {kept} fungus_species_residual rows to {output_path}")


def build_anchor_annotation_map(compactor_map, anchor_len):
    anchor_map = {}
    for sequence, hit in compactor_map.items():
        sequence = clean_sequence_candidate(sequence)
        if len(sequence) < anchor_len:
            continue
        anchor = sequence[-anchor_len:]
        current = anchor_map.get(anchor)
        if current is None:
            anchor_map[anchor] = hit
            continue
        current_priority = 0 if current.get("annotation_source") == "restricted_taxid" else 1
        hit_priority = 0 if hit.get("annotation_source") == "restricted_taxid" else 1
        if hit_priority < current_priority:
            anchor_map[anchor] = hit
    return {anchor: [hit] for anchor, hit in anchor_map.items()}


def lookup_compactor_hit(sequence, compactor_map, anchor_map, anchor_len):
    sequence = clean_sequence_candidate(sequence)
    if not sequence:
        return None
    if sequence in compactor_map:
        return compactor_map[sequence]
    if len(sequence) >= anchor_len:
        return select_compactor_hit(anchor_map.get(sequence[-anchor_len:], []))
    return None


def row_has_real_plot_annotation(row, mode):
    if mode == "blastp":
        return is_real_annotation(row.get("annotation")) or is_real_annotation(row.get("stitle"))
    feature_text = " ".join(str(row.get(col, "")) for col in ("features", "features_all", "features_10000_window"))
    return is_real_annotation(blastn_label(row)) or bool(extract_qualifier(feature_text, "product")) or bool(
        extract_qualifier(feature_text, "gene")
    )


def fake_feature_annotation(label):
    safe_label = str(label).strip()
    safe_label = safe_label.replace("'", "").replace('"', "")
    return (
        "[{'type': 'compactor', 'start': '0', 'end': '0', "
        f"'gene': ['{safe_label}'], 'product': ['{safe_label}'], "
        "'protein_seq': None, 'protein_id': None, 'note': ['COMPACTOR']}]"
    )


def synthetic_query(cluster, sequence):
    cluster = "" if cluster is None else str(cluster).strip()
    sequence = clean_sequence_candidate(sequence)
    if cluster and sequence:
        return f"{cluster}_{sequence}"
    return sequence


def summary_compactor_hit(row):
    label = clean_label(row.get("compactor_annotation")) or ""
    if not is_real_annotation(label):
        return None
    return {
        "label": label,
        "identity": row.get("identity", "NA"),
        "qcovs": row.get("qcovs", "NA"),
        "raw_annotation": row.get("compactor_raw_annotation", "NA"),
        "species_origin": row.get("compactor_species", row.get("species_origin", "NA")),
        "staxids": row.get("compactor_staxids", row.get("staxids", "NA")),
        "compactor_query": row.get("compactor_query", "NA"),
        "compactor": row.get("compactor", "NA"),
        "compactor_length": row.get("compactor_length", "NA"),
        "compactor_exact_support": row.get("compactor_exact_support", "NA"),
        "annotation_source": row.get("annotation_source", "NA"),
        "blast_mode": row.get("blast_mode", "NA"),
        "match_source": "summary",
        "metadata_category": row.get("metadata_category", ""),
        "feature": row.get("feature", ""),
        "cluster": row.get("cluster", ""),
        "sequence": clean_sequence_candidate(row.get("sequence")),
    }


def plot_match_key(row, sequence):
    return (
        row.get("metadata_category", ""),
        row.get("feature", ""),
        row.get("cluster", ""),
        clean_sequence_candidate(sequence),
    )


def build_summary_compactor_maps(summary_path):
    key_map = {}
    sequence_hits = defaultdict(list)
    summary_rows = []
    summary_path = Path(summary_path)
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        return key_map, {}, summary_rows

    with open(summary_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hit = summary_compactor_hit(row)
            if hit is None:
                continue
            sequence = ""
            for col in ("sequence", "query", "extendor", "anchor_target", "anchor_target_sequence"):
                sequence = clean_sequence_candidate(row.get(col))
                if sequence:
                    break
            if not sequence:
                continue
            key = plot_match_key(row, sequence)
            key_map[key] = hit
            sequence_hits[sequence].append(hit)
            summary_rows.append({**row, "_sequence": sequence, "_key": key, "_hit": hit})

    sequence_map = {}
    for sequence, hits in sequence_hits.items():
        sequence_map[sequence] = select_compactor_hit(hits)
    return key_map, sequence_map, summary_rows


def find_row_compactor_hit(row, compactor_map, anchor_map, anchor_len, summary_key_map=None, summary_sequence_map=None):
    sequences = []
    query_sequence = sequence_from_plot_query(row.get("query"))
    if query_sequence:
        sequences.append(query_sequence)
    for col in ("sequence", "extendor", "anchor_target", "anchor_target_sequence"):
        sequence = clean_sequence_candidate(row.get(col))
        if sequence and sequence not in sequences:
            sequences.append(sequence)

    summary_key_map = summary_key_map or {}
    summary_sequence_map = summary_sequence_map or {}
    for sequence in sequences:
        hit = summary_key_map.get(plot_match_key(row, sequence))
        if hit:
            return hit
    for sequence in sequences:
        hit = summary_sequence_map.get(sequence)
        if hit:
            return hit
    for sequence in sequences:
        hit = lookup_compactor_hit(sequence, compactor_map, anchor_map, anchor_len)
        if hit:
            return hit
    return None


def add_sequence_candidate(sequences, value):
    sequence = clean_sequence_candidate(value)
    if sequence and sequence not in sequences:
        sequences.append(sequence)


def plot_annotation_sequences(row, fieldnames):
    sequences = []
    add_sequence_candidate(sequences, sequence_from_plot_query(row.get("query")))
    for col in ("sequence", "extendor", "anchor_target", "anchor_target_sequence"):
        add_sequence_candidate(sequences, row.get(col))
    if len(fieldnames) >= 11:
        add_sequence_candidate(sequences, row.get(fieldnames[10]))
    return sequences


def find_plot_annotation_compactor_hit(
    row,
    fieldnames,
    compactor_map,
    anchor_map,
    anchor_len,
    summary_key_map=None,
    summary_sequence_map=None,
):
    sequences = plot_annotation_sequences(row, fieldnames)

    summary_key_map = summary_key_map or {}
    summary_sequence_map = summary_sequence_map or {}
    for sequence in sequences:
        hit = summary_key_map.get(plot_match_key(row, sequence))
        if hit:
            return hit
    for sequence in sequences:
        hit = summary_sequence_map.get(sequence)
        if hit:
            return hit
    for sequence in sequences:
        hit = lookup_compactor_hit(sequence, compactor_map, anchor_map, anchor_len)
        if hit:
            return hit
    return None


def template_keys(row):
    metadata_category = row.get("metadata_category", "")
    feature = row.get("feature", "")
    cluster = row.get("cluster", "")
    return [
        (metadata_category, feature, cluster),
        (metadata_category, feature, ""),
        ("", feature, cluster),
        ("", feature, ""),
    ]


def apply_compactor_hit_to_annotation_row(row, compactor_hit, mode, force=False):
    if not compactor_hit:
        return False
    should_fill = force or compactor_hit.get("match_source") == "summary" or not row_has_real_plot_annotation(row, mode)
    if not should_fill:
        return False

    row["identity"] = row.get("identity") if has_text(row.get("identity")) else compactor_hit["identity"]
    row["qcovs"] = row.get("qcovs") if has_text(row.get("qcovs")) else compactor_hit["qcovs"]
    if mode == "blastp":
        row["annotation"] = compactor_hit["label"]
        row["stitle"] = row.get("stitle") if has_text(row.get("stitle")) else compactor_hit["label"]
    else:
        feature_annotation = fake_feature_annotation(compactor_hit["label"])
        row["features"] = row.get("features") if has_text(row.get("features")) else feature_annotation
        row["features_all"] = row.get("features_all") if has_text(row.get("features_all")) else feature_annotation
        if "features_10000_window" in row:
            row["features_10000_window"] = (
                row.get("features_10000_window")
                if has_text(row.get("features_10000_window"))
                else feature_annotation
            )
    row["compactor_annotation"] = compactor_hit["label"]
    row["compactor_query"] = compactor_hit["compactor_query"]
    row["compactor_length"] = compactor_hit["compactor_length"]
    row["compactor_exact_support"] = compactor_hit["compactor_exact_support"]
    row["compactor_raw_annotation"] = compactor_hit["raw_annotation"]
    row["compactor_species"] = compactor_hit.get("species_origin", "NA")
    row["compactor_staxids"] = compactor_hit.get("staxids", "NA")
    return True


def make_synthetic_annotation_row(summary_row, template, fieldnames, mode):
    row = {col: template.get(col, "NA") for col in fieldnames}
    for col in ("metadata_category", "feature", "cluster"):
        if col in fieldnames:
            row[col] = summary_row.get(col, "NA")
    if "query" in fieldnames:
        row["query"] = synthetic_query(summary_row.get("cluster"), summary_row.get("_sequence"))
    if "evalue" in fieldnames:
        row["evalue"] = "NA"
    hit = summary_row["_hit"]
    apply_compactor_hit_to_annotation_row(row, hit, mode, force=True)
    return row


def fill_plot_annotation_tsv(
    input_path,
    output_path,
    compactor_map,
    anchor_map,
    anchor_len,
    mode,
    summary_key_map=None,
    summary_sequence_map=None,
    summary_rows=None,
):
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    filled = 0
    summary_matches = 0
    appended = 0
    skipped_no_template = 0
    if not input_path.exists():
        raise FileNotFoundError(f"Missing plot annotation TSV: {input_path}")
    with open(input_path, newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        original_fieldnames = list(reader.fieldnames or [])
        fieldnames = list(original_fieldnames)
        rows = list(reader)
        for extra_col in (
            "compactor_annotation",
            "compactor_query",
            "compactor_length",
            "compactor_exact_support",
            "compactor_raw_annotation",
            "compactor_species",
            "compactor_staxids",
        ):
            if extra_col not in fieldnames:
                fieldnames.append(extra_col)
        if mode == "blastp" and "annotation" not in fieldnames:
            fieldnames.append("annotation")
        if mode == "blast" and "features" not in fieldnames:
            fieldnames.append("features")
        if mode == "blast" and "features_all" not in fieldnames:
            fieldnames.append("features_all")

        templates = {}
        existing_summary_keys = set()
        for row in rows:
            for key in template_keys(row):
                templates.setdefault(key, row)
            for sequence in plot_annotation_sequences(row, original_fieldnames):
                existing_summary_keys.add(plot_match_key(row, sequence))

        with open(output_path, "w", newline="") as out_handle:
            writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                compactor_hit = find_plot_annotation_compactor_hit(
                    row,
                    original_fieldnames,
                    compactor_map,
                    anchor_map,
                    anchor_len,
                    summary_key_map,
                    summary_sequence_map,
                )
                if compactor_hit and compactor_hit.get("match_source") == "summary":
                    summary_matches += 1
                if apply_compactor_hit_to_annotation_row(row, compactor_hit, mode):
                    filled += 1
                else:
                    row.setdefault("compactor_annotation", "NA")
                    row.setdefault("compactor_query", "NA")
                    row.setdefault("compactor_length", "NA")
                    row.setdefault("compactor_exact_support", "NA")
                    row.setdefault("compactor_raw_annotation", "NA")
                writer.writerow({col: row.get(col, "NA") for col in fieldnames})

            for summary_row in summary_rows or []:
                key = summary_row.get("_key")
                if key in existing_summary_keys:
                    continue
                template = None
                for template_key in template_keys(summary_row):
                    template = templates.get(template_key)
                    if template is not None:
                        break
                if template is None:
                    skipped_no_template += 1
                    continue
                synthetic_row = make_synthetic_annotation_row(summary_row, template, fieldnames, mode)
                writer.writerow({col: synthetic_row.get(col, "NA") for col in fieldnames})
                existing_summary_keys.add(key)
                appended += 1
    print(
        f"Filled {filled} {mode} rows with compactor annotations in {output_path} "
        f"({summary_matches} exact-summary matches found, {appended} summary rows appended, "
        f"{skipped_no_template} summary rows skipped without a template)."
    )
    return output_path, filled


def fill_plot_summary_tsv(input_path, output_path, compactor_map, anchor_map, anchor_len):
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        print(f"Skipping compactor summary fill because {input_path} does not exist.")
        return output_path, 0

    filled = 0
    with open(input_path, newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for col in (
            "compactor_annotation",
            "compactor_query",
            "compactor_length",
            "compactor_exact_support",
            "compactor_raw_annotation",
            "compactor_species",
            "compactor_staxids",
        ):
            if col not in fieldnames:
                fieldnames.append(col)
        with open(output_path, "w", newline="") as out_handle:
            writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            for row in reader:
                sequence = ""
                for col in ("sequence", "query", "extendor", "anchor_target", "anchor_target_sequence"):
                    sequence = clean_sequence_candidate(row.get(col))
                    if sequence:
                        break
                compactor_hit = lookup_compactor_hit(sequence, compactor_map, anchor_map, anchor_len)
                label = row.get("Blast Label") or row.get("label") or row.get("annotation")
                if compactor_hit and not is_real_annotation(label):
                    if "Blast Label" in fieldnames:
                        row["Blast Label"] = compactor_hit["label"]
                    elif "label" in fieldnames:
                        row["label"] = compactor_hit["label"]
                    elif "annotation" in fieldnames:
                        row["annotation"] = compactor_hit["label"]
                    row["compactor_annotation"] = compactor_hit["label"]
                    row["compactor_query"] = compactor_hit["compactor_query"]
                    row["compactor_length"] = compactor_hit["compactor_length"]
                    row["compactor_exact_support"] = compactor_hit["compactor_exact_support"]
                    row["compactor_raw_annotation"] = compactor_hit["raw_annotation"]
                    row["compactor_species"] = compactor_hit.get("species_origin", "NA")
                    row["compactor_staxids"] = compactor_hit.get("staxids", "NA")
                    filled += 1
                else:
                    row.setdefault("compactor_annotation", "NA")
                    row.setdefault("compactor_query", "NA")
                    row.setdefault("compactor_length", "NA")
                    row.setdefault("compactor_exact_support", "NA")
                    row.setdefault("compactor_raw_annotation", "NA")
                    row.setdefault("compactor_species", "NA")
                    row.setdefault("compactor_staxids", "NA")
                writer.writerow({col: row.get(col, "NA") for col in fieldnames})
    print(f"Filled {filled} unresolved plot summary rows with compactor annotations in {output_path}")
    return output_path, filled


