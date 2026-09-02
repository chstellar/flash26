#!/usr/bin/env python3
"""
Temporary helper for replotting final compactor BLAST plots split by metadata.

Expected input in RESULTS_DIR:
  *_nonzero_coefficients_blast_annotated_plots_compactor.pdf

For each matching final PDF, this derives the companion compactor TSVs, filters
them to the selected metadata_category values, then reruns the normal compactor
plot script once per class in SPLIT_METADATA_COL.
"""

import argparse
import csv
import os
import re
import shutil
import shlex
import subprocess
import sys
from pathlib import Path


FINAL_PDF_GLOB = "*_nonzero_coefficients_blast_annotated_plots_compactor.pdf"
PLOT_SCRIPT = "src/annotation/blast_code/plot_blast_annotations_each_feature.R"
NUM_HITS = 10

FINAL_PDF_SUFFIX = "_blast_annotated_plots_compactor.pdf"
PREFIX_SUFFIX = "_nonzero_coefficients"
SUMMARY_SUFFIX = "_blast_annotated_plots_summary_compactor.tsv"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--metadata_file", required=True)
    parser.add_argument("--metadata_categories", default="infectant,infection_status")
    parser.add_argument("--split_metadata_col", default="sra_study")
    parser.add_argument("--dry_run", action="store_true")
    return parser.parse_args()


def split_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def sanitize(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "value"


def sniff_delimiter(path):
    with open(path, newline="") as handle:
        sample = handle.read(4096)
    return "," if sample.count(",") > sample.count("\t") else "\t"


def read_dicts(path):
    delimiter = sniff_delimiter(path)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader), list(reader.fieldnames or []), delimiter


def write_dicts(path, rows, fieldnames, delimiter="\t"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in fieldnames})


def filter_category_table(src, dst, categories):
    rows, fieldnames, _ = read_dicts(src)
    if "metadata_category" not in fieldnames:
        raise ValueError(f"{src} lacks metadata_category")
    kept = [row for row in rows if row.get("metadata_category") in categories]
    write_dicts(dst, kept, fieldnames)
    return len(kept)


def filter_metadata_by_split(metadata_file, split_col, split_value, dst):
    rows, fieldnames, delimiter = read_dicts(metadata_file)
    if split_col not in fieldnames:
        raise ValueError(f"{metadata_file} lacks split column {split_col!r}")
    kept = [row for row in rows if row.get(split_col) == split_value]
    write_dicts(dst, kept, fieldnames, delimiter)
    return len(kept)


def split_values(metadata_file, split_col):
    rows, fieldnames, _ = read_dicts(metadata_file)
    if split_col not in fieldnames:
        raise ValueError(f"{metadata_file} lacks split column {split_col!r}")
    values = sorted({row.get(split_col, "").strip() for row in rows if row.get(split_col, "").strip()})
    if not values:
        raise ValueError(f"{metadata_file} has no non-empty values in {split_col!r}")
    return values


def find_dataset_parts(path):
    parts = path.parts
    if "results" in parts:
        idx = len(parts) - 1 - list(reversed(parts)).index("results")
        if len(parts) >= idx + 6:
            return {
                "results_root": Path(*parts[: idx + 1]),
                "dataset": parts[idx + 1],
                "select_type": parts[idx + 2],
                "cluster_type": parts[idx + 3],
                "model": parts[idx + 4],
                "normalize": parts[idx + 5],
            }
    parent = path.parent.parts
    if len(parent) < 5:
        raise ValueError(f"Could not infer result path structure from {path}")
    return {
        "results_root": Path("results"),
        "dataset": parent[-5],
        "select_type": parent[-4],
        "cluster_type": parent[-3],
        "model": parent[-2],
        "normalize": parent[-1],
    }


def parse_prefix(prefix_name, parts):
    pattern = (
        rf"^{re.escape(parts['dataset'])}_{re.escape(parts['model'])}_(?P<prediction_task>.+?)_results_"
        rf"top(?P<num_clusters>\d+)_target(?P<target_rank>\d+)_"
        rf"k(?P<kmer_width>\d+)_s(?P<kmer_step>\d+)_"
        rf"trainProp(?P<train_proportion>[^_]+){re.escape(PREFIX_SUFFIX)}$"
    )
    match = re.match(pattern, prefix_name)
    if not match:
        raise ValueError(f"Could not parse final plot prefix: {prefix_name}")
    return match.groupdict()


def companion_paths(final_pdf):
    if not final_pdf.name.endswith(FINAL_PDF_SUFFIX):
        raise ValueError(f"Not a final compactor BLAST plot: {final_pdf}")
    prefix = final_pdf.name[: -len(FINAL_PDF_SUFFIX)]
    parts = find_dataset_parts(final_pdf)
    params = parse_prefix(prefix, parts)
    base = final_pdf.parent / prefix
    dataset = parts["dataset"]

    clusters = (
        parts["results_root"]
        / dataset
        / parts["select_type"]
        / parts["cluster_type"]
        / (
            f"{dataset}_sequences_per_cluster_top{params['num_clusters']}-clusters_"
            f"target{params['target_rank']}_k{params['kmer_width']}_s{params['kmer_step']}.tsv"
        )
    )
    sample_seqs = (
        parts["results_root"]
        / dataset
        / (
            f"{dataset}_prepared_sequences_{parts['select_type']}_{parts['cluster_type']}_"
            f"top{params['num_clusters']}_target{params['target_rank']}_"
            f"k{params['kmer_width']}_s{params['kmer_step']}_sample_sequences.tsv"
        )
    )
    feather = (
        parts["results_root"]
        / dataset
        / (
            f"{dataset}_{parts['model']}_top_variance_features_for_glmnet_"
            f"{parts['select_type']}_{parts['cluster_type']}_top{params['num_clusters']}_"
            f"target{params['target_rank']}_k{params['kmer_width']}_s{params['kmer_step']}_"
            f"{parts['normalize']}.feather"
        )
    )

    return {
        "prefix": prefix,
        "blastp": Path(str(base) + "_blastp_annotated_compactor.tsv"),
        "blast": Path(str(base) + "_blast_annotated_compactor.tsv"),
        "summary": Path(str(base) + SUMMARY_SUFFIX),
        "clusters": clusters,
        "sample_seqs": sample_seqs,
        "feather": feather,
        "cluster_length": params["kmer_width"],
    }


def require(path, label):
    if not Path(path).exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def stage_compactor_aux(summary_path, work_dir):
    prefix = str(summary_path)[: -len(SUMMARY_SUFFIX)]
    for pattern in ("_compactor_*_selected.tsv", "_compactor_*_seed_annotations.tsv"):
        matches = sorted(Path(summary_path).parent.glob(Path(prefix).name + pattern))
        matches = [path for path in matches if path.exists() and path.stat().st_size > 0]
        for src in matches:
            dst = work_dir / src.name
            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)


def main():
    args = parse_args()
    project_dir = Path(__file__).resolve().parent
    plot_script = project_dir / PLOT_SCRIPT
    results_dir = Path(args.results_dir)
    metadata_file = Path(args.metadata_file)
    categories = set(split_csv(args.metadata_categories))

    require(results_dir, "results_dir")
    require(metadata_file, "metadata_file")
    require(plot_script, "plot script")
    if not categories:
        raise SystemExit("metadata_categories is empty")

    final_pdfs = sorted(results_dir.glob(FINAL_PDF_GLOB))
    if not final_pdfs:
        raise SystemExit(f"No final compactor PDFs matched {results_dir / FINAL_PDF_GLOB}")

    values = split_values(metadata_file, args.split_metadata_col)
    failures = 0
    for final_pdf in final_pdfs:
        try:
            paths = companion_paths(final_pdf)
            for label in ("blastp", "blast", "summary", "clusters", "sample_seqs", "feather"):
                require(paths[label], label)

            category_tag = sanitize("-".join(sorted(categories)))
            for value in values:
                value_tag = sanitize(value)
                work_dir = final_pdf.parent / "split_replot_inputs" / value_tag
                filtered_blastp = work_dir / f"{paths['prefix']}_{category_tag}_blastp_annotated_compactor.tsv"
                filtered_blast = work_dir / f"{paths['prefix']}_{category_tag}_blast_annotated_compactor.tsv"
                filtered_summary = work_dir / f"{paths['prefix']}{SUMMARY_SUFFIX}"
                split_metadata = work_dir / f"metadata_{sanitize(args.split_metadata_col)}_{value_tag}.tsv"

                if filter_category_table(paths["blastp"], filtered_blastp, categories) == 0:
                    print(f"Skipping {final_pdf}: no matching metadata_category rows", file=sys.stderr)
                    break
                filter_category_table(paths["blast"], filtered_blast, categories)
                filter_category_table(paths["summary"], filtered_summary, categories)
                stage_compactor_aux(paths["summary"], work_dir)

                n_samples = filter_metadata_by_split(
                    metadata_file, args.split_metadata_col, value, split_metadata
                )
                if n_samples == 0:
                    continue

                output_pdf = (
                    final_pdf.parent
                    / f"{paths['prefix']}_blast_annotated_plots_compactor_split-by-"
                    f"{sanitize(args.split_metadata_col)}_{value_tag}.pdf"
                )
                cmd = [
                    "Rscript",
                    "--vanilla",
                    str(plot_script),
                    "--nonzero_annotations",
                    str(filtered_blastp),
                    "--clusters",
                    str(paths["clusters"]),
                    "--feather_file",
                    str(paths["feather"]),
                    "--sample_seqs",
                    str(paths["sample_seqs"]),
                    "--metadata",
                    str(split_metadata),
                    "--compactor_summary",
                    str(filtered_summary),
                    "--output",
                    str(output_pdf),
                    "--num_hits",
                    str(NUM_HITS),
                    "--cluster_length",
                    str(paths["cluster_length"]),
                ]

                print(" ".join(shlex.quote(part) for part in cmd))
                if not args.dry_run:
                    env = os.environ.copy()
                    env.setdefault("MPLBACKEND", "Agg")
                    subprocess.run(cmd, check=True, env=env)
                    print(f"Wrote {output_pdf}")
        except Exception as exc:
            failures += 1
            print(f"Error processing {final_pdf}: {exc}", file=sys.stderr)

    if failures:
        raise SystemExit(f"{failures} final plot(s) failed")


if __name__ == "__main__":
    main()
