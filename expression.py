#!/usr/bin/env python3
"""
Plot target occurrence across samples for one SPLASH anchor.

The expected SATC dump format is the FLASH intermediate produced by
satc_dump after sample names are restored:

    sample  anchor  target  count

The script writes a tidy TSV with one row per sample-target pair and a PDF
with one combined panel plus one page per target.
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
        help="Keep only samples listed in --sample_order.",
    )
    parser.add_argument(
        "--no_individual_pages",
        action="store_true",
        help="Only write the combined target plot, not one page per target.",
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


def plot_pdf(path, anchor, samples, targets, counts, totals, individual_pages):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    x = list(range(len(samples)))
    with PdfPages(path) as pdf:
        fig_width = max(10, min(28, 0.35 * len(samples) + 4))
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        for target in targets:
            y = [counts.get((sample, target), 0.0) for sample in samples]
            ax.plot(x, y, marker="o", linewidth=1.4, markersize=3, label=shorten_sequence(target))
        ax.set_title(f"Anchor {anchor}: target occurrences")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Occurrence count")
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=90, fontsize=7)
        ax.grid(axis="y", alpha=0.25)
        if len(targets) <= 20:
            ax.legend(title="Target", fontsize=7, title_fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        if individual_pages:
            for target in targets:
                y = [counts.get((sample, target), 0.0) for sample in samples]
                fig, ax = plt.subplots(figsize=(fig_width, 5))
                ax.plot(x, y, marker="o", linewidth=1.6, markersize=4, color="#3366AA")
                ax.set_title(
                    f"Anchor {anchor}\nTarget {shorten_sequence(target, 80)} | total={format_count(totals[target])}"
                )
                ax.set_xlabel("Sample")
                ax.set_ylabel("Occurrence count")
                ax.set_xticks(x)
                ax.set_xticklabels(samples, rotation=90, fontsize=7)
                ax.grid(axis="y", alpha=0.25)
                fig.tight_layout()
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
    if args.strict_order:
        samples = [sample for sample in ordered_samples if sample in all_samples or sample in observed_samples]
    else:
        samples = ordered_samples + sorted(
            observed_samples.difference(ordered_samples), key=natural_key
        )

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
    plot_path = output_prefix.with_suffix(".target_expression.pdf")
    write_tidy_tsv(tidy_tsv, args.anchor, samples, targets, counts, totals)
    plot_pdf(
        plot_path,
        args.anchor,
        samples,
        targets,
        counts,
        totals,
        individual_pages=not args.no_individual_pages,
    )

    print(f"Anchor: {args.anchor}")
    print(f"Observed samples with this anchor: {len(observed_samples)}")
    print(f"Targets plotted: {len(targets)}")
    print(f"Wrote counts: {tidy_tsv}")
    print(f"Wrote plot: {plot_path}")


if __name__ == "__main__":
    main()
