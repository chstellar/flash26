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


DEFAULT_INPUT_DIR = "/scratch/users/jiamuyu/proj_botryllus/splash2/260713_01_3ants_challenge"
DEFAULT_OUTPUT_DIR = (
    "/scratch/users/jiamuyu/proj_botryllus/flash/results/"
    "260714-00-3ants-challenge/filter1/noCluster/hyena/normalized"
)
DEFAULT_TAXIDS = "300111;102681;104421"
DEFAULT_BLAST_DB = "/scratch/users/jiamuyu/dabs_ref/blast/"
REPO_ROOT = Path(__file__).resolve().parent


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Pick representative resfungi compactors, BLAST/BLASTP them, and annotate "
            "seed extendors by their anchor-matched compactors."
        )
    )
    parser.add_argument("--input_dir", default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seeds", default=None, help="Defaults to INPUT_DIR/seeds.resfungi.raw")
    parser.add_argument("--taxids", default=DEFAULT_TAXIDS)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--translation_table", default="1")
    parser.add_argument("--entrez_email", default=os.environ.get("ENTREZ_EMAIL", "v8514616@outlook.com"))
    parser.add_argument("--blast_db", default=DEFAULT_BLAST_DB)
    parser.add_argument("--anchor_len", type=int, default=31)
    parser.add_argument("--thresholds", default="1000,100,5")
    parser.add_argument("--skip_blast", action="store_true")
    parser.add_argument(
        "--overwrite_blast",
        action="store_true",
        help="Re-run BLAST/BLASTP even if the expected output TSVs already exist.",
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
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_") or "resfungi"


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


def write_fasta(records, output_fasta):
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w", newline="") as handle:
        for record in records:
            header = (
                f"{record['query']}|anchor={record['anchor']}|length={record['length']}"
                f"|exact_support={record['exact_support']:g}|row={record['row_index']}"
                f"|source={Path(record['source_file']).name}"
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
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for record in records:
            writer.writerow({col: record.get(col, "NA") for col in columns})


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


def run_blasts(args, output_fasta, output_dir):
    blast = output_dir / "resfungi_compactors_blast.tsv"
    reblast = output_dir / "resfungi_compactors_reblast.tsv"
    blastp = output_dir / "resfungi_compactors_blastp.tsv"
    reblastp = output_dir / "resfungi_compactors_reblastp.tsv"
    temp_root = output_dir / "resfungi_compactors_tmp"
    command_env = make_python_shim(temp_root)
    print(f"Wrapper subprocess python: {which('python', path=command_env['PATH'])}", flush=True)
    run_command(["python", "--version"], env=command_env)

    if args.overwrite_blast or not (has_output(blast) and has_output(reblast)):
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
                str(reblast),
            ],
            env=command_env,
            cwd=REPO_ROOT,
        )
    else:
        print(f"Reusing existing BLASTN outputs: {blast} and {reblast}")
    require_output(blast, "Restricted BLASTN")
    require_output(reblast, "Unrestricted reBLASTN")
    if args.overwrite_blast or not (has_output(blastp) and has_output(reblastp)):
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
                "0",
                "all",
                "",
                "10",
                "",
                "0",
                str(reblastp),
            ],
            env=command_env,
            cwd=REPO_ROOT,
        )
    else:
        print(f"Reusing existing BLASTP outputs: {blastp} and {reblastp}")
    require_output(blastp, "Restricted BLASTP")
    require_output(reblastp, "Unrestricted reBLASTP")
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


def read_annotation_table(path, mode, source):
    annotations = {}
    if not path.exists() or path.stat().st_size == 0:
        return annotations
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "query" not in reader.fieldnames:
            return annotations
        for row in reader:
            query = row.get("query")
            if not query:
                continue
            label = blastp_label(row) if mode == "blastp" else blastn_label(row)
            if not label:
                continue
            current = annotations.get(query)
            identity = row.get("identity", "NA")
            qcovs = row.get("qcovs", "NA")
            if current is None:
                annotations[query] = {
                    f"{source}_{mode}_label": label,
                    f"{source}_{mode}_identity": identity,
                    f"{source}_{mode}_qcovs": qcovs,
                }
    return annotations


def merge_annotation_maps(*maps):
    merged = defaultdict(dict)
    for annotation_map in maps:
        for query, values in annotation_map.items():
            merged[query].update(values)
    return dict(merged)


def choose_final_annotation(row):
    restricted = [
        row.get("restricted_blastp_label"),
        row.get("restricted_blast_label"),
    ]
    outside = [
        row.get("unrestricted_blastp_label"),
        row.get("unrestricted_blast_label"),
    ]
    restricted = [label for label in restricted if has_text(label) and label != "UNANNOTATED"]
    outside = [label for label in outside if has_text(label) and label != "UNANNOTATED"]
    if restricted:
        return ";".join(dict.fromkeys(restricted)), "restricted_taxid"
    if outside:
        return ";".join(dict.fromkeys(outside)), "outside_taxid"
    if any(row.get(col) == "UNANNOTATED" for col in row):
        return "UNANNOTATED", "hit_without_parsed_annotation"
    return "NO MATCH", "no_hit"


def write_compactor_annotation_summary(records, annotations, output_path):
    columns = [
        "query",
        "anchor",
        "compactor",
        "length",
        "exact_support",
        "source_file",
        "annotation_label",
        "annotation_source",
        "restricted_blastp_label",
        "restricted_blast_label",
        "unrestricted_blastp_label",
        "unrestricted_blast_label",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = {**record, **annotations.get(record["query"], {})}
            label, source = choose_final_annotation(row)
            row["annotation_label"] = label
            row["annotation_source"] = source
            writer.writerow({col: row.get(col, "NA") for col in columns})


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
                    "raw_seed_row": line,
                }
            )
    return rows


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
        "restricted_blastp_label",
        "restricted_blast_label",
        "unrestricted_blastp_label",
        "unrestricted_blast_label",
        "raw_seed_row",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        for seed in seeds:
            record = records_by_anchor.get(seed["seed_anchor"], {})
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
                row.update(annotations.get(query, {}))
            label, source = choose_final_annotation(row)
            row["annotation_label"] = label
            row["annotation_source"] = source if record else "no_anchor_matched_compactor"
            writer.writerow({col: row.get(col, "NA") for col in columns})


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    seeds_path = Path(args.seeds) if args.seeds else input_dir / "seeds.resfungi.raw"
    thresholds = [float(item) for item in args.thresholds.split(",") if item.strip()]

    records = []
    for path in sorted(input_dir.glob("*resfungi.tsv")):
        if path.stat().st_size == 0:
            print(f"Skipping empty file: {path}")
            continue
        chosen = choose_compactor(path, thresholds)
        if chosen is None:
            print(f"Skipping {path}: no representative passed thresholds {thresholds}")
            continue
        chosen["query"] = f"{safe_id(path)}__row{chosen['row_index']}__len{chosen['length']}"
        records.append(chosen)

    output_fasta = output_dir / "resfungi_compactors.fasta"
    write_fasta(records, output_fasta)
    selected_path = output_dir / "resfungi_compactors_selected.tsv"
    write_selected_compactors(records, selected_path)
    print(f"Wrote {len(records)} representative compactors to {output_fasta}")
    print(f"Wrote selected compactor table to {selected_path}")

    blast = output_dir / "resfungi_compactors_blast.tsv"
    reblast = output_dir / "resfungi_compactors_reblast.tsv"
    blastp = output_dir / "resfungi_compactors_blastp.tsv"
    reblastp = output_dir / "resfungi_compactors_reblastp.tsv"
    if records and not args.skip_blast:
        blast, reblast, blastp, reblastp = run_blasts(args, output_fasta, output_dir)

    annotations = merge_annotation_maps(
        read_annotation_table(blastp, "blastp", "restricted"),
        read_annotation_table(blast, "blast", "restricted"),
        read_annotation_table(reblastp, "blastp", "unrestricted"),
        read_annotation_table(reblast, "blast", "unrestricted"),
    )
    compactor_summary_path = output_dir / "resfungi_compactor_annotations.tsv"
    write_compactor_annotation_summary(records, annotations, compactor_summary_path)
    print(f"Wrote compactor annotation summary to {compactor_summary_path}")

    records_by_anchor = {record["anchor"]: record for record in records}
    summary_path = output_dir / "seeds_resfungi_compactor_annotations.tsv"
    write_seed_annotation_summary(
        read_seed_rows(seeds_path, args.anchor_len),
        records_by_anchor,
        annotations,
        summary_path,
    )
    print(f"Wrote seed extendor compactor annotation summary to {summary_path}")


if __name__ == "__main__":
    main()
