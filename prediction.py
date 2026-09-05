#!/usr/bin/env python3

import argparse
import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D


REQUIRED_COLUMNS = {
    "row_type",
    "metadata_category",
    "matrix",
    "true_label",
    "predicted_label",
    "n_samples",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Replot selected classification or regression predictions from an Adelie sidecar table."
    )
    parser.add_argument("--input", required=True, help="Confusion-matrix sidecar TSV/CSV.")
    parser.add_argument(
        "--metadata-category", required=True, help="Exact metadata_category to plot."
    )
    parser.add_argument("--output", required=True, help="Output multipage PDF.")
    parser.add_argument(
        "--matrix",
        choices=("train", "test", "both"),
        default="both",
        help="Dataset split to plot. Default: both.",
    )
    parser.add_argument(
        "--plot-type",
        choices=("auto", "regression", "confusion"),
        default="auto",
        help="Plot type. Auto detects it from row_type. Default: auto.",
    )
    parser.add_argument(
        "--color-column",
        default="",
        help="Sidecar column used to color regression points.",
    )
    parser.add_argument(
        "--dot-colors",
        default="",
        help="Comma-separated value:color mapping for a categorical regression color column.",
    )
    parser.add_argument(
        "--class-colors",
        default="",
        help="Comma-separated class:color mapping for confusion-matrix tick labels.",
    )
    parser.add_argument(
        "--class-order",
        default="",
        help="Optional comma-separated confusion-matrix class order.",
    )
    parser.add_argument("--point-color", default="#2F6F9F")
    parser.add_argument("--missing-color", default="#A6A6A6")
    parser.add_argument("--regression-cmap", default="viridis")
    parser.add_argument("--confusion-cmap", default="viridis")
    parser.add_argument("--point-size", type=float, default=42)
    parser.add_argument("--point-alpha", type=float, default=0.82)
    parser.add_argument("--cell-font-size", type=float, default=12)
    parser.add_argument("--axis-title-font-size", type=float, default=13)
    parser.add_argument("--axis-text-font-size", type=float, default=11)
    parser.add_argument("--title-font-size", type=float, default=14)
    parser.add_argument("--other-font-size", type=float, default=10)
    parser.add_argument("--legend-font-size", type=float, default=9)
    parser.add_argument("--width", type=float, default=7.2)
    parser.add_argument("--height", type=float, default=6.6)
    return parser.parse_args()


def delimiter_for(path):
    suffixes = [suffix.lower() for suffix in Path(path).suffixes]
    if suffixes and suffixes[-1] in {".gz", ".bz2", ".xz", ".zip"}:
        suffixes.pop()
    return "," if suffixes and suffixes[-1] == ".csv" else "\t"


def resolve_input_path(path_pattern):
    matches = sorted(Path(match) for match in glob.glob(str(path_pattern)))
    matches = [path for path in matches if path.is_file()]
    if not matches:
        raise ValueError(f"No sidecar file matched: {path_pattern}")
    if len(matches) > 1:
        formatted = "\n  ".join(str(path) for path in matches)
        raise ValueError(
            f"Sidecar wildcard matched {len(matches)} files; "
            f"make SIDECAR more specific:\n  {formatted}"
        )
    return matches[0]


def read_sidecar(path):
    path = resolve_input_path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Input sidecar does not exist or is empty: {path}")
    data = pd.read_csv(path, sep=delimiter_for(path), dtype=str, keep_default_na=False)
    missing = sorted(REQUIRED_COLUMNS - set(data.columns))
    if missing:
        raise ValueError(f"Sidecar is missing required columns: {', '.join(missing)}")
    return data


def parse_list(value):
    return [item.strip() for item in str(value).split(",") if item.strip()]


def parse_color_map(value, option_name):
    result = {}
    for item in parse_list(value):
        if ":" not in item:
            raise ValueError(f"{option_name} entry must use value:color format: {item}")
        label, color = (part.strip() for part in item.split(":", 1))
        if not label or not mcolors.is_color_like(color):
            raise ValueError(f"Invalid {option_name} entry: {item}")
        result[label] = color
    return result


def requested_matrices(value):
    return ["train", "test"] if value == "both" else [value]


def detect_plot_type(data, requested):
    row_types = set(data["row_type"])
    has_regression = "prediction" in row_types
    has_confusion = "entry" in row_types
    if requested != "auto":
        expected = "prediction" if requested == "regression" else "entry"
        if expected not in row_types:
            raise ValueError(f"Selected category does not contain {requested} rows.")
        return requested
    if has_regression and not has_confusion:
        return "regression"
    if has_confusion and not has_regression:
        return "confusion"
    raise ValueError("Could not uniquely infer plot type from the selected category.")


def ordered_labels(entries, requested_order):
    observed = set(entries["true_label"]) | set(entries["predicted_label"])
    ordered = [label for label in requested_order if label in observed]
    return ordered + sorted(observed - set(ordered))


def text_color_for_background(rgba):
    red, green, blue = rgba[:3]
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "white" if luminance < 0.48 else "black"


def plot_confusion(data, matrix_name, args, class_colors, class_order):
    entries = data[(data["matrix"] == matrix_name) & (data["row_type"] == "entry")].copy()
    if entries.empty:
        raise ValueError(f"No confusion-table entries found for matrix={matrix_name}.")
    entries["n_samples"] = pd.to_numeric(entries["n_samples"], errors="coerce").fillna(0)
    labels = ordered_labels(entries, class_order)
    counts = np.zeros((len(labels), len(labels)), dtype=float)
    label_index = {label: index for index, label in enumerate(labels)}
    for row in entries.itertuples(index=False):
        counts[label_index[row.true_label], label_index[row.predicted_label]] += row.n_samples

    row_sums = counts.sum(axis=1, keepdims=True)
    proportions = np.divide(
        counts, row_sums, out=np.zeros_like(counts), where=row_sums != 0
    )
    total = counts.sum()
    accuracy = np.trace(counts) / total if total else np.nan

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    image = ax.imshow(proportions, cmap=args.confusion_cmap, vmin=0, vmax=1)
    for row_index in range(len(labels)):
        for column_index in range(len(labels)):
            rgba = image.cmap(image.norm(proportions[row_index, column_index]))
            ax.text(
                column_index,
                row_index,
                f"{counts[row_index, column_index]:g}",
                ha="center",
                va="center",
                fontsize=args.cell_font_size,
                color=text_color_for_background(rgba),
            )

    ax.set_xticks(range(len(labels)), labels=labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels=labels)
    ax.tick_params(axis="both", labelsize=args.axis_text_font_size, length=0)
    for tick in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        tick.set_color(class_colors.get(tick.get_text(), "black"))
        if tick.get_text() in class_colors:
            tick.set_fontweight("bold")

    ax.set_xlabel("Predicted label", fontsize=args.axis_title_font_size)
    ax.set_ylabel("True label", fontsize=args.axis_title_font_size)
    metric = f"accuracy = {accuracy:.3f}" if np.isfinite(accuracy) else "accuracy = NA"
    ax.set_title(
        f"{args.metadata_category}\n{matrix_name.capitalize()}, n = {int(total)}, {metric}",
        fontsize=args.title_font_size,
        pad=14,
    )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.047, pad=0.04)
    colorbar.set_label("Row proportion", fontsize=args.axis_title_font_size)
    colorbar.ax.tick_params(labelsize=args.other_font_size)
    fig.tight_layout()
    return fig


def regression_r2(observed, predicted):
    denominator = np.sum((observed - np.mean(observed)) ** 2)
    if len(observed) < 2 or denominator == 0:
        return np.nan
    return 1 - np.sum((observed - predicted) ** 2) / denominator


def regression_limits(observed, predicted):
    minimum = min(np.min(observed), np.min(predicted))
    maximum = max(np.max(observed), np.max(predicted))
    padding = (maximum - minimum) * 0.06 or 0.5
    return minimum - padding, maximum + padding


def categorical_point_colors(values, supplied_colors):
    categories = sorted(pd.unique(values))
    fallback = plt.get_cmap("tab10")
    colors = {
        category: supplied_colors.get(category, fallback(index % 10))
        for index, category in enumerate(categories)
    }
    return categories, colors


def plot_regression(data, matrix_name, args, dot_colors):
    predictions = data[
        (data["matrix"] == matrix_name) & (data["row_type"] == "prediction")
    ].copy()
    if predictions.empty:
        raise ValueError(f"No regression predictions found for matrix={matrix_name}.")
    predictions["observed"] = pd.to_numeric(predictions["true_label"], errors="coerce")
    predictions["predicted"] = pd.to_numeric(
        predictions["predicted_label"], errors="coerce"
    )
    predictions = predictions.dropna(subset=["observed", "predicted"])
    if predictions.empty:
        raise ValueError(f"No finite regression predictions found for matrix={matrix_name}.")

    observed = predictions["observed"].to_numpy(dtype=float)
    predicted = predictions["predicted"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(args.width, args.height))

    if args.color_column:
        if args.color_column not in predictions.columns:
            raise ValueError(f"Color column not found in sidecar: {args.color_column}")
        raw_colors = predictions[args.color_column].astype(str).str.strip()
        missing_tokens = raw_colors.str.lower().isin(
            {"", "nan", "na", "none", "null", "n/a"}
        )
        raw_colors = raw_colors.mask(missing_tokens, "Missing")
        numeric_colors = pd.to_numeric(raw_colors.replace("Missing", np.nan), errors="coerce")
        nonmissing = raw_colors != "Missing"
        is_numeric = nonmissing.any() and numeric_colors[nonmissing].notna().all()
        if is_numeric:
            missing = ~np.isfinite(numeric_colors.to_numpy(dtype=float))
            if missing.any():
                ax.scatter(
                    observed[missing],
                    predicted[missing],
                    s=args.point_size,
                    alpha=args.point_alpha,
                    color=args.missing_color,
                    edgecolors="white",
                    linewidths=0.45,
                    label="Missing",
                )
            finite = ~missing
            points = ax.scatter(
                observed[finite],
                predicted[finite],
                c=numeric_colors.to_numpy(dtype=float)[finite],
                cmap=args.regression_cmap,
                s=args.point_size,
                alpha=args.point_alpha,
                edgecolors="white",
                linewidths=0.45,
            )
            colorbar = fig.colorbar(points, ax=ax, fraction=0.047, pad=0.04)
            colorbar.set_label(args.color_column, fontsize=args.axis_title_font_size)
            colorbar.ax.tick_params(labelsize=args.other_font_size)
            if missing.any():
                ax.legend(frameon=False, fontsize=args.legend_font_size)
        else:
            values = raw_colors.astype(str).to_numpy()
            categories, colors = categorical_point_colors(values, dot_colors)
            for category in categories:
                selected = values == category
                ax.scatter(
                    observed[selected],
                    predicted[selected],
                    s=args.point_size,
                    alpha=args.point_alpha,
                    color=colors[category],
                    edgecolors="white",
                    linewidths=0.45,
                )
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markerfacecolor=colors[category],
                    markeredgecolor="white",
                    markersize=7,
                    label=category,
                )
                for category in categories
            ]
            ax.legend(
                handles=handles,
                title=args.color_column,
                frameon=False,
                fontsize=args.legend_font_size,
                title_fontsize=args.other_font_size,
            )
    else:
        ax.scatter(
            observed,
            predicted,
            s=args.point_size,
            alpha=args.point_alpha,
            color=args.point_color,
            edgecolors="white",
            linewidths=0.45,
        )

    lower, upper = regression_limits(observed, predicted)
    ax.plot([lower, upper], [lower, upper], color="#222222", linewidth=1.2, linestyle="--")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#D9D9D9", linewidth=0.6, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=args.axis_text_font_size)
    ax.set_xlabel("Observed value", fontsize=args.axis_title_font_size)
    ax.set_ylabel("Predicted value", fontsize=args.axis_title_font_size)
    r2 = regression_r2(observed, predicted)
    metric = f"R2 = {r2:.3f}" if np.isfinite(r2) else "R2 = NA"
    ax.set_title(
        f"{args.metadata_category}\n{matrix_name.capitalize()}, n = {len(observed)}, {metric}",
        fontsize=args.title_font_size,
        pad=14,
    )
    fig.tight_layout()
    return fig


def main():
    args = parse_args()
    sidecar = read_sidecar(args.input)
    selected = sidecar[sidecar["metadata_category"] == args.metadata_category].copy()
    if selected.empty:
        available = sorted(sidecar["metadata_category"].drop_duplicates())
        raise ValueError(
            f"metadata_category not found: {args.metadata_category}. Available: {', '.join(available)}"
        )

    plot_type = detect_plot_type(selected, args.plot_type)
    matrices = requested_matrices(args.matrix)
    available_matrices = set(selected["matrix"])
    missing_matrices = [name for name in matrices if name not in available_matrices]
    if missing_matrices:
        raise ValueError(
            f"Requested matrix data not found: {', '.join(missing_matrices)}. "
            f"Available: {', '.join(sorted(available_matrices))}"
        )

    class_colors = parse_color_map(args.class_colors, "--class-colors")
    dot_colors = parse_color_map(args.dot_colors, "--dot-colors")
    class_order = parse_list(args.class_order)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(output) as pdf:
        for matrix_name in matrices:
            if plot_type == "regression":
                figure = plot_regression(selected, matrix_name, args, dot_colors)
            else:
                figure = plot_confusion(
                    selected, matrix_name, args, class_colors, class_order
                )
            pdf.savefig(figure, bbox_inches="tight")
            plt.close(figure)

    print(f"Wrote {len(matrices)} {plot_type} plot(s) to {output}")


if __name__ == "__main__":
    main()
