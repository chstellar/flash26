#!/usr/bin/env python3
"""
Plot target occurrence across samples for one SPLASH anchor.

The expected SATC dump format is the FLASH intermediate produced by
satc_dump after sample names are restored:

    sample  anchor  target  count

The script writes a tidy TSV with one row per sample-target pair, a grouped
TSV for time-course plotting, and a PDF with one panel per target.
"""

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create target occurrence line charts for one anchor."
    )
    parser.add_argument("--anchor", required=True, help="Anchor sequence to plot.")
    parser.add_argument(
        "--satc",
        required=True,
        help="Path to all_satc.filtered.dump or all_satc_merged.txt.",
    )
    parser.add_argument(
        "--sample_order",
        required=True,
        help=(
            "Sample order string. Use comma/semicolon/whitespace-separated sample "
            "names, or 'auto' for natural observed-sample order. Unlisted samples "
            "are appended naturally unless --strict_order is set."
        ),
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Output prefix for .target_counts.tsv and .target_expression.pdf.",
    )
    parser.add_argument(
        "--top_targets",
        type=int,
        default=0,
        help="Limit plotting to the top N targets by total count. 0 means all targets.",
    )
    parser.add_argument(
        "--min_total_count",
        type=float,
        default=0,
        help="Only plot targets with total count at least this value.",
    )
    parser.add_argument(
        "--strict_order",
        action="store_true",
        help="Alias for --only_ordered_samples.",
    )
    parser.add_argument(
        "--only_ordered_samples",
        action="store_true",
        help="Plot only samples listed in --sample_order.",
    )
    parser.add_argument(
        "--plot_mode",
        choices=("plain", "tc7", "manipulation"),
        default="plain",
        help=(
            "plain: one line over samples. tc7: one line each for A/B/C over "
            "TC7 time indices. manipulation: six lines for control/O/B x H/T."
        ),
    )
    parser.add_argument(
        "--aggregate",
        choices=("mean", "sum"),
        default="mean",
        help="How to combine replicate samples within a plotted time/group point.",
    )
    parser.add_argument(
        "--grid_cols",
        type=int,
        default=2,
        help="Number of target panels per PDF row.",
    )
    parser.add_argument(
        "--grid_rows",
        type=int,
        default=3,
        help="Number of target panels per PDF page column.",
    )
    parser.add_argument(
        "--no_individual_pages",
        action="store_true",
        help="Deprecated; retained for compatibility.",
    )
    return parser.parse_args()


def natural_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def split_order(order_string):
    if order_string.strip().lower() in {"auto", "natural"}:
        return []
    tokens = [
        token.strip()
        for token in re.split(r"[,;\s]+", order_string.strip())
        if token.strip()
    ]
    seen = set()
    ordered = []
    for token in tokens:
        if token not in seen:
            ordered.append(token)
            seen.add(token)
    return ordered


def read_anchor_rows(satc_path, anchor):
    rows = []
    all_samples = set()
    with open(satc_path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, fields in enumerate(reader, start=1):
            if not fields or len(fields) < 4:
                continue
            if line_no == 1 and fields[0].lower() in {"sample", "sample_id"}:
                continue
            sample, row_anchor, target, count = fields[:4]
            sample = sample.strip()
            row_anchor = row_anchor.strip()
            target = target.strip()
            all_samples.add(sample)
            if row_anchor != anchor:
                continue
            try:
                count_value = float(count)
            except ValueError:
                continue
            rows.append((sample, target, count_value))
    return rows, all_samples


def collapse_counts(rows):
    counts = defaultdict(float)
    for sample, target, count in rows:
        counts[(sample, target)] += count
    return counts


def target_totals(counts):
    totals = defaultdict(float)
    for (_, target), count in counts.items():
        totals[target] += count
    return totals


def write_tidy_tsv(path, anchor, samples, targets, counts, totals):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "sample_index", "anchor", "target", "count", "target_total_count"])
        for sample_index, sample in enumerate(samples, start=1):
            for target in targets:
                writer.writerow(
                    [
                        sample,
                        sample_index,
                        anchor,
                        target,
                        format_count(counts.get((sample, target), 0.0)),
                        format_count(totals[target]),
                    ]
                )


def tc7_metadata(sample):
    match = re.match(r"^TC7-(\d+)([ABC])_", sample)
    if not match:
        return None
    return {
        "x": int(match.group(1)),
        "x_label": str(int(match.group(1))),
        "line": match.group(2),
        "line_label": match.group(2),
    }


def manipulation_metadata(sample):
    if match := re.match(r"^0_C([HT])(\d+)$", sample):
        tissue = "Head" if match.group(1) == "H" else "Thorax"
        return {
            "x": 0,
            "x_label": "Control",
            "line": f"Control {tissue}",
            "line_label": f"Control {tissue}",
        }
    if match := re.match(r"^([123])_([BO])L?([A-Z0-9]+)([HT])(\d+)$", sample):
        stage_code = match.group(3)
        stage_order = {"T25": 1, "M": 2, "DM": 3}
        stage_label = {"T25": "LT25", "M": "LM", "DM": "DM"}
        infection = {"B": "B. bassiana", "O": "Ophiocordyceps"}[match.group(2)]
        tissue = "Head" if match.group(4) == "H" else "Thorax"
        if stage_code not in stage_order:
            return None
        return {
            "x": stage_order[stage_code],
            "x_label": stage_label[stage_code],
            "line": f"{infection} {tissue}",
            "line_label": f"{infection} {tissue}",
        }
    return None


def plain_metadata(sample, index):
    return {
        "x": index,
        "x_label": sample,
        "line": "count",
        "line_label": "count",
    }


def build_sample_metadata(samples, mode):
    metadata = {}
    for index, sample in enumerate(samples, start=1):
        if mode == "tc7":
            item = tc7_metadata(sample)
        elif mode == "manipulation":
            item = manipulation_metadata(sample)
        else:
            item = plain_metadata(sample, index)
        if item is not None:
            metadata[sample] = item
    return metadata


def aggregate_values(values, method):
    if not values:
        return 0.0
    if method == "sum":
        return sum(values)
    return sum(values) / len(values)


def grouped_counts(samples, targets, counts, sample_metadata, aggregate):
    grouped = defaultdict(list)
    for sample in samples:
        meta = sample_metadata.get(sample)
        if meta is None:
            continue
        for target in targets:
            grouped[(target, meta["x"], meta["x_label"], meta["line"], meta["line_label"])].append(
                counts.get((sample, target), 0.0)
            )
    return {
        key: (aggregate_values(values, aggregate), len(values))
        for key, values in grouped.items()
    }


def write_grouped_tsv(path, anchor, grouped):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["anchor", "target", "x", "x_label", "line", "line_label", "count", "n_samples"]
        )
        for key in sorted(grouped, key=lambda item: (natural_key(item[0]), item[1], natural_key(item[3]))):
            target, x_value, x_label, line, line_label = key
            count, n_samples = grouped[key]
            writer.writerow(
                [
                    anchor,
                    target,
                    x_value,
                    x_label,
                    line,
                    line_label,
                    format_count(count),
                    n_samples,
                ]
            )


def format_count(value):
    if math.isfinite(value) and float(value).is_integer():
        return str(int(value))
    return f"{value:.6g}"


def shorten_sequence(sequence, max_len=42):
    if len(sequence) <= max_len:
        return sequence
    left = max_len // 2 - 2
    right = max_len - left - 3
    return f"{sequence[:left]}...{sequence[-right:]}"


def plot_panel_pdf(path, anchor, targets, totals, grouped, grid_cols, grid_rows):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    line_colors = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#9467bd",
        "#ff7f0e",
        "#17becf",
        "#8c564b",
        "#e377c2",
    ]
    panels_per_page = max(1, grid_cols * grid_rows)
    with PdfPages(path) as pdf:
        for page_start in range(0, len(targets), panels_per_page):
            page_targets = targets[page_start : page_start + panels_per_page]
            fig, axes = plt.subplots(
                grid_rows,
                grid_cols,
                figsize=(5.8 * grid_cols, 3.9 * grid_rows),
                squeeze=False,
            )
            fig.suptitle(f"Anchor {anchor}: target occurrences", fontsize=13)
            for ax, target in zip(axes.ravel(), page_targets):
                target_rows = {
                    key: value
                    for key, value in grouped.items()
                    if key[0] == target
                }
                lines = sorted({key[3] for key in target_rows}, key=natural_key)
                x_labels = {
                    key[1]: key[2]
                    for key in target_rows
                }
                x_values = sorted(x_labels)
                for color_index, line in enumerate(lines):
                    y_values = []
                    for x_value in x_values:
                        matches = [
                            value[0]
                            for key, value in target_rows.items()
                            if key[1] == x_value and key[3] == line
                        ]
                        y_values.append(matches[0] if matches else math.nan)
                    ax.plot(
                        x_values,
                        y_values,
                        marker="o",
                        linewidth=1.5,
                        markersize=3.5,
                        label=line,
                        color=line_colors[color_index % len(line_colors)],
                    )
                ax.set_title(
                    f"{shorten_sequence(target, 46)}\ntotal={format_count(totals[target])}",
                    fontsize=9,
                )
                ax.set_ylabel("Occurrence")
                ax.set_xticks(x_values)
                ax.set_xticklabels([x_labels[x] for x in x_values], rotation=45, ha="right", fontsize=8)
                ax.grid(axis="y", alpha=0.25)
                if lines:
                    ax.legend(fontsize=6, frameon=False)
            for ax in axes.ravel()[len(page_targets):]:
                ax.axis("off")
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig)
            plt.close(fig)


def main():
    args = parse_args()
    rows, all_samples = read_anchor_rows(args.satc, args.anchor)
    if not rows:
        raise SystemExit(f"No rows found for anchor {args.anchor} in {args.satc}")

    counts = collapse_counts(rows)
    totals = target_totals(counts)

    ordered_samples = split_order(args.sample_order)
    observed_samples = {sample for sample, _, _ in rows}
    only_ordered_samples = args.strict_order or args.only_ordered_samples
    if only_ordered_samples:
        samples = [sample for sample in ordered_samples if sample in all_samples or sample in observed_samples]
    else:
        samples = ordered_samples + sorted(
            observed_samples.difference(ordered_samples), key=natural_key
        )
    if not samples:
        raise SystemExit("No samples selected for plotting.")
    rows = [row for row in rows if row[0] in set(samples)]
    counts = collapse_counts(rows)
    totals = target_totals(counts)

    targets = [
        target
        for target, total in sorted(totals.items(), key=lambda item: (-item[1], natural_key(item[0])))
        if total >= args.min_total_count
    ]
    if args.top_targets > 0:
        targets = targets[: args.top_targets]
    if not targets:
        raise SystemExit("No targets passed --min_total_count/--top_targets filters.")

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    tidy_tsv = output_prefix.with_suffix(".target_counts.tsv")
    grouped_tsv = output_prefix.with_suffix(".target_grouped_counts.tsv")
    plot_path = output_prefix.with_suffix(".target_expression.pdf")
    write_tidy_tsv(tidy_tsv, args.anchor, samples, targets, counts, totals)
    sample_metadata = build_sample_metadata(samples, args.plot_mode)
    if args.plot_mode != "plain":
        missing_metadata = [sample for sample in samples if sample not in sample_metadata]
        if missing_metadata:
            print(
                "Skipping samples not matching "
                f"{args.plot_mode} naming pattern: {', '.join(missing_metadata)}"
            )
    grouped = grouped_counts(samples, targets, counts, sample_metadata, args.aggregate)
    write_grouped_tsv(grouped_tsv, args.anchor, grouped)
    plot_panel_pdf(
        plot_path,
        args.anchor,
        targets,
        totals,
        grouped,
        grid_cols=args.grid_cols,
        grid_rows=args.grid_rows,
    )

    print(f"Anchor: {args.anchor}")
    print(f"Observed samples with this anchor: {len(observed_samples)}")
    print(f"Targets plotted: {len(targets)}")
    print(f"Wrote counts: {tidy_tsv}")
    print(f"Wrote grouped counts: {grouped_tsv}")
    print(f"Wrote plot: {plot_path}")


if __name__ == "__main__":
    main()
