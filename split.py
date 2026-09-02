#!/usr/bin/env python3
"""
Temporary helper to replot BLAST annotation PDFs split by one metadata column.

This is intentionally not wired into Snakemake. It scans a results folder for
BLASTP annotation TSVs, keeps selected metadata_category rows, and calls the
existing split-aware R plotting script once per matching TSV.
"""

import argparse
import csv
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path


BLASTP_SUFFIX = "_nonzero_coefficients_blastp_annotated.tsv"
BLAST_SUFFIX = "_nonzero_coefficients_blast_annotated.tsv"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results_dir", required=True, help="Dataset results folder to scan.")
    parser.add_argument(
        "--metadata_categories",
        default="infectant,infection_status",
        help="Comma-separated metadata_category values to replot.",
    )
    parser.add_argument(
        "--split_metadata_col",
        default="sra_study",
        help="Metadata column used to split each detailed BLAST plot.",
    )
    parser.add_argument(
        "--metadata_file",
        default="",
        help="Metadata TSV/CSV. If blank, infer from --dataset_table.",
    )
    parser.add_argument(
        "--dataset_table",
        default="dataset_table.csv",
        help="CSV with dataset_short_name and metadata_file columns.",
    )
    parser.add_argument(
        "--temp_dir",
        default="results",
        help="Root containing prepared sample sequence and feather files.",
    )
    parser.add_argument(
        "--plot_script",
        default="src/annotation/blast_code/plot_blast_annotations_each_feature_split_species.R",
        help="R plotting script to call.",
    )
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable.")
    parser.add_argument("--num_hits", type=int, default=10, help="Top coefficients to plot.")
    parser.add_argument(
        "--output_dir",
        default="",
        help="Output folder. Default: <results_dir>/split_replots.",
    )
    parser.add_argument(
        "--work_dir",
        default="",
        help="Folder for filtered temporary TSVs. Default: <output_dir>/filtered_inputs.",
    )
    parser.add_argument(
        "--file_glob",
        default=f"*{BLASTP_SUFFIX}",
        help="Glob, relative to --results_dir, selecting BLASTP annotation TSVs.",
    )
    parser.add_argument(
        "--products",
        action="store_true",
        help="Pass --products to the R plotting script.",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print commands without running R.",
    )
    return parser.parse_args()


def split_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def sanitize(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "value"


def read_table(path):
    with open(path, newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "," if sample.count(",") > sample.count("\t") else "\t"
        yield from csv.DictReader(handle, delimiter=delimiter)


def write_filtered_table(src, dst, categories):
    rows_written = 0
    with open(src, newline="") as in_handle:
        sample = in_handle.read(4096)
        in_handle.seek(0)
        delimiter = "," if sample.count(",") > sample.count("\t") else "\t"
        reader = csv.DictReader(in_handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"{src} has no header")
        if "metadata_category" not in reader.fieldnames:
            raise ValueError(f"{src} lacks required column metadata_category")
        dst.parent.mkdir(parents=True, exist_ok=True)
        with open(dst, "w", newline="") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                if row.get("metadata_category") in categories:
                    writer.writerow(row)
                    rows_written += 1
    return rows_written


def find_dataset_parts(blastp_path):
    parts = blastp_path.parts
    if "results" in parts:
        idx = len(parts) - 1 - list(reversed(parts)).index("results")
        if len(parts) >= idx + 6:
            return {
                "dataset": parts[idx + 1],
                "select_type": parts[idx + 2],
                "cluster_type": parts[idx + 3],
                "model": parts[idx + 4],
                "normalize": parts[idx + 5],
            }
    parent_parts = blastp_path.parent.parts
    if len(parent_parts) >= 5:
        return {
            "dataset": parent_parts[-5],
            "select_type": parent_parts[-4],
            "cluster_type": parent_parts[-3],
            "model": parent_parts[-2],
            "normalize": parent_parts[-1],
        }
    raise ValueError(f"Could not infer dataset path parts from {blastp_path}")


def infer_results_root(path):
    parts = path.parts
    if "results" in parts:
        idx = len(parts) - 1 - list(reversed(parts)).index("results")
        return Path(*parts[: idx + 1])
    return Path("results")


def parse_filename(path, parts):
    name = path.name
    escaped_dataset = re.escape(parts["dataset"])
    escaped_model = re.escape(parts["model"])
    pattern = (
        rf"^{escaped_dataset}_{escaped_model}_(?P<prediction_task>.+?)_results_"
        rf"top(?P<num_clusters>\d+)_target(?P<target_rank>\d+)_"
        rf"k(?P<kmer_width>\d+)_s(?P<kmer_step>\d+)_"
        rf"trainProp(?P<train_proportion>[^_]+){re.escape(BLASTP_SUFFIX)}$"
    )
    match = re.match(pattern, name)
    if not match:
        raise ValueError(f"Could not parse FLASH plot filename: {path}")
    out = match.groupdict()
    out["cluster_length"] = int(out["kmer_width"])
    return out


def infer_metadata_file(dataset, dataset_table):
    if not dataset_table or not dataset_table.exists():
        return ""
    rows = list(read_table(dataset_table))
    exact = [row for row in rows if row.get("dataset_short_name") == dataset]
    fuzzy = [
        row
        for row in rows
        if row.get("dataset_short_name", "").startswith(dataset)
        or dataset.startswith(row.get("dataset_short_name", ""))
    ]
    matches = exact or fuzzy
    if len(matches) == 1:
        return matches[0].get("metadata_file", "")
    return ""


def require_path(path, label):
    if not path or not Path(path).exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def main():
    args = parse_args()
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir) if args.output_dir else results_dir / "split_replots"
    work_dir = Path(args.work_dir) if args.work_dir else output_dir / "filtered_inputs"
    categories = set(split_csv(args.metadata_categories))
    if not categories:
        raise SystemExit("--metadata_categories did not contain any values")

    require_path(results_dir, "results_dir")
    require_path(args.plot_script, "plot_script")

    blastp_files = sorted(results_dir.glob(args.file_glob))
    blastp_files = [path for path in blastp_files if path.name.endswith(BLASTP_SUFFIX)]
    if not blastp_files:
        raise SystemExit(f"No BLASTP annotation files matched {results_dir / args.file_glob}")

    failures = 0
    for blastp_path in blastp_files:
        try:
            parts = find_dataset_parts(blastp_path)
            params = parse_filename(blastp_path, parts)
            dataset = parts["dataset"]
            results_root = infer_results_root(blastp_path)
            metadata_file = args.metadata_file or infer_metadata_file(
                dataset, Path(args.dataset_table)
            )

            blast_path = Path(str(blastp_path).replace("blastp_annotated", "blast_annotated"))
            cluster_path = (
                results_root
                / dataset
                / parts["select_type"]
                / parts["cluster_type"]
                / (
                    f"{dataset}_sequences_per_cluster_top{params['num_clusters']}-clusters_"
                    f"target{params['target_rank']}_k{params['kmer_width']}_s{params['kmer_step']}.tsv"
                )
            )
            sample_seqs = (
                Path(args.temp_dir)
                / dataset
                / (
                    f"{dataset}_prepared_sequences_{parts['select_type']}_{parts['cluster_type']}_"
                    f"top{params['num_clusters']}_target{params['target_rank']}_"
                    f"k{params['kmer_width']}_s{params['kmer_step']}_sample_sequences.tsv"
                )
            )
            feather_file = (
                Path(args.temp_dir)
                / dataset
                / (
                    f"{dataset}_{parts['model']}_top_variance_features_for_glmnet_"
                    f"{parts['select_type']}_{parts['cluster_type']}_top{params['num_clusters']}_"
                    f"target{params['target_rank']}_k{params['kmer_width']}_s{params['kmer_step']}_"
                    f"{parts['normalize']}.feather"
                )
            )

            require_path(blast_path, "matching BLAST annotation TSV")
            require_path(cluster_path, "clusters TSV")
            require_path(sample_seqs, "sample sequences TSV")
            require_path(feather_file, "feather file")
            require_path(metadata_file, "metadata file")

            stem = blastp_path.name.replace(BLASTP_SUFFIX, "")
            category_tag = sanitize("-".join(sorted(categories)))
            split_tag = sanitize(args.split_metadata_col)
            filtered_blastp = work_dir / f"{stem}_{category_tag}_blastp_annotated.tsv"
            filtered_blast = work_dir / f"{stem}_{category_tag}_blast_annotated.tsv"
            rows_blastp = write_filtered_table(blastp_path, filtered_blastp, categories)
            rows_blast = write_filtered_table(blast_path, filtered_blast, categories)
            if rows_blastp == 0:
                print(f"Skipping {blastp_path}: no rows matched {sorted(categories)}", file=sys.stderr)
                continue
            if rows_blast == 0:
                print(f"Warning: {blast_path} had no matching rows", file=sys.stderr)

            run_out_dir = output_dir / blastp_path.parent.relative_to(results_dir)
            run_out_dir.mkdir(parents=True, exist_ok=True)
            output_pdf = run_out_dir / f"{stem}_{category_tag}_split-by-{split_tag}.pdf"

            cmd = [
                args.rscript,
                "--vanilla",
                args.plot_script,
                "--nonzero_annotations",
                str(filtered_blastp),
                "--clusters",
                str(cluster_path),
                "--feather_file",
                str(feather_file),
                "--sample_seqs",
                str(sample_seqs),
                "--metadata",
                str(metadata_file),
                "--output",
                str(output_pdf),
                "--num_hits",
                str(args.num_hits),
                "--cluster_length",
                str(params["cluster_length"]),
                "--split_metadata_col",
                args.split_metadata_col,
            ]
            if args.products:
                cmd.append("--products")

            print(" ".join(shlex.quote(str(part)) for part in cmd))
            if not args.dry_run:
                env = os.environ.copy()
                env.setdefault("FLASH_SPLIT_METADATA_COL", args.split_metadata_col)
                subprocess.run(cmd, check=True, env=env)
                print(f"Wrote {output_pdf}")
        except Exception as exc:
            failures += 1
            print(f"Error processing {blastp_path}: {exc}", file=sys.stderr)

    if failures:
        raise SystemExit(f"{failures} file(s) failed")


if __name__ == "__main__":
    main()
