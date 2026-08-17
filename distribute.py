#!/usr/bin/env python3
"""
Compute per-partition count distributions for extendor sequences.

Inputs:
  1. A TSV containing extendor/seed sequences in one selected column.
  2. all_satc_merged.txt, with rows like: sample, anchor, target, count.
  3. A sample partition TSV/CSV with: sample, major_partition, minor_partition.

The output TSV preserves the input rows and appends distributions across major
partitions and across every major.minor intersection.
"""

import argparse
import csv
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path


MISSING_PARTITION = "NA"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize each extendor's SATC counts by major partitions and by "
            "major.minor partition intersections."
        )
    )
    parser.add_argument("--input", required=True, help="Input TSV containing extendors.")
    parser.add_argument(
        "--output",
        required=True,
        help="Output summary TSV. Input columns are preserved and distribution columns are appended.",
    )
    parser.add_argument(
        "--satc",
        default="all_satc_merged.txt",
        help="sample-anchor-target-count TSV used by expression.py.",
    )
    parser.add_argument(
        "--partition_tsv",
        required=True,
        help="TSV/CSV with sample name, major partition, and minor partition columns.",
    )
    parser.add_argument(
        "--partition_delimiter",
        choices=("auto", "tab", "comma"),
        default="auto",
        help="Delimiter for --partition_tsv. Default: auto from file content/extension.",
    )
    parser.add_argument(
        "--extendor_col",
        default="1",
        help="Extendor column in --input. Use a header name or 1-based index. Default: 1.",
    )
    parser.add_argument(
        "--sample_col",
        default="1",
        help="Sample column in --partition_tsv. Use a header name or 1-based index. Default: 1.",
    )
    parser.add_argument(
        "--major_col",
        default="2",
        help="Major partition column in --partition_tsv. Use a header name or 1-based index. Default: 2.",
    )
    parser.add_argument(
        "--minor_col",
        default="3",
        help="Minor partition column in --partition_tsv. Use a header name or 1-based index. Default: 3.",
    )
    parser.add_argument(
        "--anchor_len",
        type=int,
        default=None,
        help="Anchor length. If omitted, inferred from --anchor_satc, then --satc.",
    )
    parser.add_argument(
        "--target_len",
        type=int,
        default=None,
        help="Target length. If omitted, inferred from --satc or from extendor length minus anchor length.",
    )
    parser.add_argument(
        "--anchor_satc",
        default="anchor_satc_merged.txt",
        help="Optional sample-anchor-count TSV used to infer anchor length.",
    )
    parser.add_argument(
        "--extendor_order",
        choices=("target-anchor", "anchor-target"),
        default="target-anchor",
        help=(
            "How the extendor sequence is arranged. Existing helper scripts in this "
            "repo use target-anchor, where the anchor is the suffix. Default: target-anchor."
        ),
    )
    parser.add_argument(
        "--input_has_header",
        choices=("auto", "yes", "no"),
        default="auto",
        help="Whether --input has a header. Default: auto.",
    )
    parser.add_argument(
        "--partition_has_header",
        choices=("auto", "yes", "no"),
        default="auto",
        help="Whether --partition_tsv has a header. Default: auto.",
    )
    parser.add_argument(
        "--satc_has_header",
        choices=("auto", "yes", "no"),
        default="auto",
        help="Whether --satc has a header. Default: auto.",
    )
    parser.add_argument(
        "--case_sensitive",
        action="store_true",
        help="Keep sequence case when matching. By default sequences are uppercased.",
    )
    parser.add_argument(
        "--distribution_sep",
        default="/",
        help="Separator between entries in distribution strings. Default: '/'.",
    )
    parser.add_argument(
        "--pair_sep",
        default=":",
        help="Separator between partition labels and counts. Default: ':'.",
    )
    parser.add_argument(
        "--intersection_sep",
        default=".",
        help="Separator between major and minor labels. Default: '.'.",
    )
    parser.add_argument(
        "--long_output",
        default="",
        help="Optional long-format TSV with one row per extendor x partition.",
    )
    return parser.parse_args()


def natural_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def looks_like_header(fields, requested_cols):
    lower = [field.strip().lower() for field in fields]
    if lower[:4] == ["sample", "anchor", "target", "count"]:
        return True
    common_header_names = {
        "sample",
        "sample_id",
        "anchor",
        "target",
        "count",
        "major",
        "minor",
        "major_partition",
        "minor_partition",
        "partition",
        "extendor",
        "seed_extendor",
        "seed",
        "sequence",
        "query",
    }
    if any(field in common_header_names for field in lower):
        return True
    for col in requested_cols:
        if not col.isdigit() and col.strip().lower() in lower:
            return True
    return False


def resolve_col(fields, spec, path):
    if spec.isdigit():
        index = int(spec) - 1
        if index < 0 or index >= len(fields):
            raise ValueError(f"Column index {spec} is out of range for {path}")
        return index
    try:
        return fields.index(spec)
    except ValueError as exc:
        raise ValueError(f"Column {spec!r} was not found in {path}") from exc


def normalize_seq(value, case_sensitive):
    value = value.strip()
    return value if case_sensitive else value.upper()


def infer_delimiter(path, mode):
    if mode == "tab":
        return "\t"
    if mode == "comma":
        return ","
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return ","
    if path.suffix.lower() in {".tsv", ".txt"}:
        return "\t"
    with open(path, newline="") as handle:
        for line in handle:
            if line.strip():
                return "," if line.count(",") > line.count("\t") else "\t"
    return "\t"


def read_rows(path, delimiter="\t"):
    with open(path, newline="") as handle:
        yield from csv.reader(handle, delimiter=delimiter)


def read_tsv_rows(path):
    yield from read_rows(path, "\t")


def detect_first_data_row(path, requested_cols, header_mode, delimiter="\t"):
    for fields in read_rows(path, delimiter):
        if not fields or all(not field.strip() for field in fields):
            continue
        if header_mode == "yes":
            return fields, True
        if header_mode == "no":
            return fields, False
        return fields, looks_like_header(fields, requested_cols)
    raise ValueError(f"{path} is empty")


def read_input(path, extendor_col, header_mode, case_sensitive):
    first, has_header = detect_first_data_row(path, [extendor_col], header_mode)
    extendor_idx = resolve_col(first, extendor_col, path) if has_header else resolve_col(first, extendor_col, path)
    header = list(first) if has_header else [f"col{i}" for i in range(1, len(first) + 1)]
    rows = []
    extendors = OrderedDict()

    def add_row(fields, line_no):
        if len(fields) <= extendor_idx:
            raise ValueError(f"{path}:{line_no} does not contain extendor column {extendor_col}")
        extendor = normalize_seq(fields[extendor_idx], case_sensitive)
        row = dict(zip(header, fields + [""] * (len(header) - len(fields))))
        row["_extendor_value"] = extendor
        rows.append(row)
        extendors.setdefault(extendor, None)

    line_no = 0
    for fields in read_tsv_rows(path):
        if not fields or all(not field.strip() for field in fields):
            continue
        line_no += 1
        if line_no == 1 and has_header:
            continue
        add_row(fields, line_no)

    return header, rows, list(extendors.keys())


def infer_anchor_len(anchor_satc, satc, satc_header_mode, case_sensitive):
    anchor_path = Path(anchor_satc)
    if anchor_path.exists() and anchor_path.stat().st_size > 0:
        for fields in read_tsv_rows(anchor_path):
            if len(fields) < 2:
                continue
            if fields[0].strip().lower() in {"sample", "sample_id"}:
                continue
            anchor = normalize_seq(fields[1], case_sensitive)
            if anchor:
                return len(anchor), f"{anchor_satc} column 2"

    for fields in iter_satc_rows(satc, satc_header_mode):
        _, anchor, _, _ = fields
        if anchor:
            return len(normalize_seq(anchor, case_sensitive)), f"{satc} anchor column"
    raise ValueError("Could not infer anchor length; pass --anchor_len")


def infer_target_len(satc, satc_header_mode, case_sensitive):
    for fields in iter_satc_rows(satc, satc_header_mode):
        _, _, target, _ = fields
        target = normalize_seq(target, case_sensitive)
        if target:
            return len(target), f"{satc} target column"
    raise ValueError("Could not infer target length; pass --target_len")


def iter_satc_rows(path, header_mode):
    first, has_header = detect_first_data_row(path, [], header_mode)
    del first
    line_no = 0
    for fields in read_tsv_rows(path):
        if not fields or all(not field.strip() for field in fields):
            continue
        line_no += 1
        if line_no == 1 and has_header:
            continue
        if len(fields) < 4:
            continue
        yield fields[:4]


def split_extendor(extendor, anchor_len, target_len, order):
    expected_len = anchor_len + target_len
    if len(extendor) < expected_len:
        raise ValueError(
            f"Extendor {extendor!r} is shorter than anchor_len + target_len ({expected_len})"
        )
    if len(extendor) > expected_len:
        extendor = extendor[:expected_len]
    if order == "target-anchor":
        return extendor[-anchor_len:], extendor[:target_len]
    return extendor[:anchor_len], extendor[anchor_len : anchor_len + target_len]


def add_ordered(value, values, seen):
    if value not in seen:
        values.append(value)
        seen.add(value)


def read_partitions(path, sample_col, major_col, minor_col, header_mode, delimiter_mode):
    delimiter = infer_delimiter(path, delimiter_mode)
    first, has_header = detect_first_data_row(
        path, [sample_col, major_col, minor_col], header_mode, delimiter
    )
    sample_idx = resolve_col(first, sample_col, path)
    major_idx = resolve_col(first, major_col, path)
    minor_idx = resolve_col(first, minor_col, path)

    sample_to_partition = {}
    major_values, minor_values = [], []
    seen_major, seen_minor = set(), set()

    line_no = 0
    for fields in read_rows(path, delimiter):
        if not fields or all(not field.strip() for field in fields):
            continue
        line_no += 1
        if line_no == 1 and has_header:
            continue
        max_idx = max(sample_idx, major_idx, minor_idx)
        if len(fields) <= max_idx:
            raise ValueError(f"{path}:{line_no} does not contain all requested partition columns")
        sample = fields[sample_idx].strip()
        major = fields[major_idx].strip() or MISSING_PARTITION
        minor = fields[minor_idx].strip() or MISSING_PARTITION
        if not sample:
            continue
        sample_to_partition[sample] = (major, minor)
        add_ordered(major, major_values, seen_major)
        add_ordered(minor, minor_values, seen_minor)

    if not sample_to_partition:
        raise ValueError(f"No sample partitions were read from {path}")
    return sample_to_partition, major_values, minor_values


def format_count(value):
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.10g}"


def distribution_string(labels, counts, sep, pair_sep):
    return sep.join(f"{label}{pair_sep}{format_count(counts.get(label, 0.0))}" for label in labels)


def write_summary(path, header, rows, major_labels, combo_labels, totals, major_counts, combo_counts, args):
    out_header = header + [
        "anchor",
        "target",
        "total_count",
        "major_distribution",
        "major_minor_distribution",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(out_header)
        for row in rows:
            extendor = row["_extendor_value"]
            writer.writerow(
                [row.get(col, "") for col in header]
                + [
                    row["_anchor"],
                    row["_target"],
                    format_count(totals.get(extendor, 0.0)),
                    distribution_string(
                        major_labels,
                        major_counts.get(extendor, {}),
                        args.distribution_sep,
                        args.pair_sep,
                    ),
                    distribution_string(
                        combo_labels,
                        combo_counts.get(extendor, {}),
                        args.distribution_sep,
                        args.pair_sep,
                    ),
                ]
            )


def write_long(path, rows, major_labels, combo_labels, totals, major_counts, combo_counts):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["extendor", "anchor", "target", "partition_type", "partition", "count", "total_count"])
        for row in rows:
            extendor = row["_extendor_value"]
            for major in major_labels:
                writer.writerow(
                    [
                        extendor,
                        row["_anchor"],
                        row["_target"],
                        "major",
                        major,
                        format_count(major_counts.get(extendor, {}).get(major, 0.0)),
                        format_count(totals.get(extendor, 0.0)),
                    ]
                )
            for combo in combo_labels:
                writer.writerow(
                    [
                        extendor,
                        row["_anchor"],
                        row["_target"],
                        "major_minor",
                        combo,
                        format_count(combo_counts.get(extendor, {}).get(combo, 0.0)),
                        format_count(totals.get(extendor, 0.0)),
                    ]
                )


def main():
    args = parse_args()
    header, input_rows, extendors = read_input(
        args.input, args.extendor_col, args.input_has_header, args.case_sensitive
    )
    if not input_rows:
        raise SystemExit("No extendors were read from --input")

    anchor_len = args.anchor_len
    target_len = args.target_len
    if anchor_len is None:
        anchor_len, source = infer_anchor_len(
            args.anchor_satc, args.satc, args.satc_has_header, args.case_sensitive
        )
        print(f"Inferred anchor_len={anchor_len} from {source}", file=sys.stderr)
    if target_len is None:
        try:
            target_len, source = infer_target_len(args.satc, args.satc_has_header, args.case_sensitive)
            print(f"Inferred target_len={target_len} from {source}", file=sys.stderr)
        except ValueError:
            lengths = sorted({len(extendor) for extendor in extendors})
            if len(lengths) == 1 and lengths[0] > anchor_len:
                target_len = lengths[0] - anchor_len
                print(f"Inferred target_len={target_len} from extendor length", file=sys.stderr)
            else:
                raise

    requested_pairs = {}
    for row in input_rows:
        anchor, target = split_extendor(row["_extendor_value"], anchor_len, target_len, args.extendor_order)
        row["_anchor"] = anchor
        row["_target"] = target
        requested_pairs[(anchor, target)] = row["_extendor_value"]

    sample_to_partition, major_labels, minor_labels = read_partitions(
        args.partition_tsv,
        args.sample_col,
        args.major_col,
        args.minor_col,
        args.partition_has_header,
        args.partition_delimiter,
    )
    def make_combo_labels():
        return [
            f"{major}{args.intersection_sep}{minor}"
            for major in major_labels
            for minor in minor_labels
        ]

    combo_labels = make_combo_labels()

    totals = defaultdict(float)
    major_counts = defaultdict(lambda: defaultdict(float))
    combo_counts = defaultdict(lambda: defaultdict(float))
    missing_samples = set()

    for sample, anchor, target, count_text in iter_satc_rows(args.satc, args.satc_has_header):
        anchor = normalize_seq(anchor, args.case_sensitive)
        target = normalize_seq(target, args.case_sensitive)
        extendor = requested_pairs.get((anchor, target))
        if extendor is None:
            continue
        try:
            count = float(count_text)
        except ValueError:
            continue
        major, minor = sample_to_partition.get(sample, (MISSING_PARTITION, MISSING_PARTITION))
        if sample not in sample_to_partition:
            missing_samples.add(sample)
        labels_changed = False
        if major == MISSING_PARTITION and major not in major_labels:
            major_labels.append(major)
            labels_changed = True
        if minor == MISSING_PARTITION and minor not in minor_labels:
            minor_labels.append(minor)
            labels_changed = True
        if labels_changed:
            combo_labels = make_combo_labels()
        combo = f"{major}{args.intersection_sep}{minor}"
        totals[extendor] += count
        major_counts[extendor][major] += count
        combo_counts[extendor][combo] += count

    write_summary(
        args.output,
        header,
        input_rows,
        major_labels,
        combo_labels,
        totals,
        major_counts,
        combo_counts,
        args,
    )
    if args.long_output:
        write_long(args.long_output, input_rows, major_labels, combo_labels, totals, major_counts, combo_counts)

    print(f"Wrote {args.output}")
    if args.long_output:
        print(f"Wrote {args.long_output}")
    print(f"Processed {len(input_rows)} input row(s), {len(requested_pairs)} unique anchor-target pair(s).")
    if missing_samples:
        print(
            f"Warning: {len(missing_samples)} SATC sample(s) were missing from {args.partition_tsv}; "
            f"their counts were assigned to {MISSING_PARTITION}.{MISSING_PARTITION}.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
