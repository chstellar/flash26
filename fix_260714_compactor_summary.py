#!/usr/bin/env python3
"""
One-off fixer for the 260714 compactor plot summary.

It updates legacy placeholder labels and backfills class-specific coefficient
metadata from sibling annotated coefficient TSVs when possible.
"""

import csv
import math
import re
import shutil
from pathlib import Path


SUMMARY = Path(
    "/scratch/users/jiamuyu/proj_botryllus/flash/results/"
    "260714-00-3ants-challenge/filter1/noCluster/hyena/normalized/"
    "260714-00-3ants-challenge_hyena_adelie_results_top2000_target1_k41_s41_"
    "trainProp0.8_nonzero_coefficients_blast_annotated_plots_summary_compactor.tsv"
)


def split_values(value):
    value = "" if value is None else str(value).strip()
    value = re.sub(r"^\[|\]$", "", value)
    return [item.strip() for item in value.split(",") if item.strip()]


def coef_values(value):
    value = "" if value is None else str(value)
    return [
        float(item)
        for item in re.findall(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?", value)
    ]


def normalize_text(value):
    if value is None:
        return value
    text = str(value)
    if not text:
        return text
    text = text.replace("UNCHARACTERISED", "UNCHARACTERIZED")
    text = re.sub(r"\b([A-Za-z]+)RISED\b", r"\1RIZED", text)
    text = re.sub(r"\b([A-Za-z]+)rised\b", r"\1rized", text)
    if re.search(
        r"uncharacteri[sz]ed|hypothetical protein|predicted protein|unnamed protein",
        text,
        flags=re.IGNORECASE,
    ):
        return "UNANNOTATED"
    if text == "UNCHARACTERIZED":
        return "UNANNOTATED"
    return text


def class_info(row):
    classes = split_values(row.get("classes", ""))
    coefs = coef_values(row.get("coefficients", ""))
    if classes and coefs:
        n = min(len(classes), len(coefs))
        idx = max(range(n), key=lambda i: abs(coefs[i]))
        return {
            "classes": ",".join(classes),
            "max_coefficient_class": classes[idx],
            "max_coefficient_signed": f"{coefs[idx]:.10g}",
            "significant_class": classes[idx],
            "significant_coefficient": f"{coefs[idx]:.10g}",
        }

    first_class = row.get("first_class", "")
    first_coef = row.get("first_coef", "")
    if first_class and first_coef:
        return {
            "max_coefficient_class": first_class,
            "max_coefficient_signed": first_coef,
            "significant_class": first_class,
            "significant_coefficient": first_coef,
        }
    return {}


def row_key(row):
    return (
        row.get("metadata_category", ""),
        row.get("feature", ""),
        row.get("cluster", ""),
    )


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sibling_sources(summary_path):
    name = summary_path.name.replace(
        "_nonzero_coefficients_blast_annotated_plots_summary_compactor.tsv", ""
    )
    suffixes = [
        "_nonzero_coefficients_blastp_annotated_compactor.tsv",
        "_nonzero_coefficients_blast_annotated_compactor.tsv",
        "_nonzero_coefficients_blastp_annotated.tsv",
        "_nonzero_coefficients_blast_annotated.tsv",
        "_nonzero_coefficients_annotated.tsv",
        "_nonzero_coefficients.tsv",
    ]
    return [summary_path.with_name(name + suffix) for suffix in suffixes]


def build_class_lookup(summary_path):
    lookup = {}
    for path in sibling_sources(summary_path):
        if not path.exists() or path.stat().st_size == 0:
            continue
        for row in read_tsv(path):
            info = class_info(row)
            if info:
                lookup.setdefault(row_key(row), info)
    return lookup


def main():
    summary_path = SUMMARY
    if not summary_path.exists():
        raise SystemExit(f"Missing summary TSV: {summary_path}")

    backup = summary_path.with_suffix(summary_path.suffix + ".bak")
    if not backup.exists():
        shutil.copy2(summary_path, backup)

    rows = read_tsv(summary_path)
    class_lookup = build_class_lookup(summary_path)

    extra_cols = [
        "max_coefficient_class",
        "max_coefficient_signed",
        "significant_class",
        "significant_coefficient",
    ]
    fieldnames = list(rows[0].keys()) if rows else []
    for col in extra_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    changed_labels = 0
    filled_classes = 0
    for row in rows:
        for col in list(row.keys()):
            old = row[col]
            new = normalize_text(old)
            if new != old:
                row[col] = new
                changed_labels += 1

        info = class_info(row) or class_lookup.get(row_key(row), {})
        if info:
            for col in extra_cols:
                if not row.get(col):
                    row[col] = info.get(col, "")
            filled_classes += 1

    write_tsv(summary_path, rows, fieldnames)
    print(f"Wrote {summary_path}")
    print(f"Backup: {backup}")
    print(f"Normalized label cells: {changed_labels}")
    print(f"Rows with class info filled: {filled_classes}")


if __name__ == "__main__":
    main()
