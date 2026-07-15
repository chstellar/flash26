#!/usr/bin/env python3
import argparse
import csv
import glob
import math
from pathlib import Path


def resolve_results_dir(arg):
    path = Path(arg)
    if path.exists():
        return path

    run_prefix = arg.replace("_", "-")
    matches = sorted(Path("results").glob(f"{run_prefix}*"))
    matches = [match for match in matches if match.is_dir()]
    if not matches:
        raise FileNotFoundError(
            f"No results directory found for {arg!r}; tried results/{run_prefix}*"
        )
    if len(matches) > 1:
        print(
            "Multiple matching results directories found; using "
            f"{matches[0]}: {', '.join(str(match) for match in matches)}"
        )
    return matches[0]


def parse_metadata_counts(metadata):
    counts = []
    for item in str(metadata).split("/"):
        item = item.strip()
        if not item or ":" not in item:
            continue
        _, value = item.rsplit(":", 1)
        try:
            count = float(value)
        except ValueError:
            continue
        if count > 0:
            counts.append(count)
    return counts


def shannon_entropy(counts):
    total = sum(counts)
    if total <= 0:
        return math.inf
    entropy = 0.0
    for count in counts:
        p = count / total
        entropy -= p * math.log2(p)
    return entropy


def is_unannotated(label):
    label = str(label).strip()
    return label == "" or label.upper() == "NA" or label.upper() == "NAN"


def output_path_for(summary_path):
    suffix = "blast_annotated_plots_summary.tsv"
    path = str(summary_path)
    if path.endswith(suffix):
        return Path(path[: -len(suffix)] + "blast_annotated_plots_summary_unanno.tsv")
    return summary_path.with_name(summary_path.stem + "_unanno.tsv")


def process_summary(summary_path):
    with open(summary_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{summary_path} is empty or missing a header")
        if "label" not in reader.fieldnames:
            raise ValueError(f"{summary_path} does not contain a 'label' column")
        if "metadata" not in reader.fieldnames:
            raise ValueError(f"{summary_path} does not contain a 'metadata' column")

        rows = []
        for row in reader:
            if not is_unannotated(row.get("label", "")):
                continue
            counts = parse_metadata_counts(row.get("metadata", ""))
            total = sum(counts)
            entropy = shannon_entropy(counts)
            max_fraction = max(counts) / total if total > 0 else 0.0
            row["_metadata_entropy"] = f"{entropy:.8g}" if math.isfinite(entropy) else "NA"
            row["_metadata_total"] = f"{total:.8g}"
            row["_metadata_max_fraction"] = f"{max_fraction:.8g}"
            rows.append((entropy, -max_fraction, -total, row))

    rows.sort(key=lambda item: (item[0], item[1], item[2]))
    output_path = output_path_for(summary_path)
    fieldnames = reader.fieldnames + [
        "_metadata_entropy",
        "_metadata_total",
        "_metadata_max_fraction",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for _, _, _, row in rows:
            writer.writerow(row)
    return output_path, len(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Collect unannotated BLAST plot summary rows and rank them by metadata entropy."
    )
    parser.add_argument(
        "run",
        help="Run shorthand such as 260622_02, or a path such as results/260622-02-temnothorax-challenge",
    )
    args = parser.parse_args()

    results_dir = resolve_results_dir(args.run)
    pattern = str(results_dir / "**" / "*blast_annotated_plots_summary.tsv")
    summaries = [Path(path) for path in sorted(glob.glob(pattern, recursive=True))]
    summaries = [
        path for path in summaries if not str(path).endswith("_unanno.tsv")
    ]
    if not summaries:
        raise FileNotFoundError(
            f"No blast_annotated_plots_summary.tsv files found under {results_dir}"
        )

    for summary_path in summaries:
        output_path, n_rows = process_summary(summary_path)
        print(f"Wrote {n_rows} unannotated rows to {output_path}")


if __name__ == "__main__":
    main()
