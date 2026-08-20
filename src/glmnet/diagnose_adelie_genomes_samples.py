import argparse
from pathlib import Path

import pandas as pd
import pyarrow.feather as feather


def parse_args():
    parser = argparse.ArgumentParser(
        description="Diagnose sample/class counts used by run_adelie_genomes.py"
    )
    parser.add_argument("--train_features", required=True)
    parser.add_argument("--train_metadata", required=True)
    parser.add_argument("--test_features", required=True)
    parser.add_argument("--test_metadata", required=True)
    parser.add_argument("--min_samples", type=int, default=28)
    parser.add_argument(
        "--output",
        default="",
        help="Optional TSV path for the per-metadata-column diagnostic table.",
    )
    return parser.parse_args()


def read_features(path):
    data = feather.read_feather(path)
    if "sample_name" not in data.columns:
        raise ValueError(f"{path} does not contain a sample_name column")
    data["sample_name"] = data["sample_name"].astype(str)
    return data


def read_metadata(path):
    data = pd.read_table(path, dtype=str)
    if "sample_name" not in data.columns:
        raise ValueError(f"{path} does not contain a sample_name column")
    data["sample_name"] = data["sample_name"].astype(str)
    return data


def clean_label_series(series):
    series = series.replace({"nan": pd.NA, "NA": pd.NA, "": pd.NA})
    return series.dropna()


def format_counts(counts):
    if counts.empty:
        return ""
    return ";".join(f"{idx}:{val}" for idx, val in counts.items())


def count_after_feature_merge(features, metadata, col):
    merged = features[["sample_name"]].merge(
        metadata[["sample_name", col]], on="sample_name", how="left"
    )
    return clean_label_series(merged[col]).value_counts()


def main():
    args = parse_args()

    train_features = read_features(args.train_features)
    test_features = read_features(args.test_features)
    train_metadata = read_metadata(args.train_metadata)
    test_metadata = read_metadata(args.test_metadata)

    train_feature_samples = set(train_features["sample_name"])
    test_feature_samples = set(test_features["sample_name"])
    train_metadata_samples = set(train_metadata["sample_name"])
    test_metadata_samples = set(test_metadata["sample_name"])

    print("=== Sample-name overlap ===")
    print(f"train feature rows:          {len(train_features):,}")
    print(f"train feature samples:       {len(train_feature_samples):,}")
    print(f"train metadata rows:         {len(train_metadata):,}")
    print(f"train metadata samples:      {len(train_metadata_samples):,}")
    print(
        f"train feature/metadata overlap: "
        f"{len(train_feature_samples & train_metadata_samples):,}"
    )
    print()
    print(f"test feature rows:           {len(test_features):,}")
    print(f"test feature samples:        {len(test_feature_samples):,}")
    print(f"test metadata rows:          {len(test_metadata):,}")
    print(f"test metadata samples:       {len(test_metadata_samples):,}")
    print(
        f"test feature/metadata overlap:  "
        f"{len(test_feature_samples & test_metadata_samples):,}"
    )
    print()

    common_metadata_cols = sorted(
        (set(train_metadata.columns) & set(test_metadata.columns)) - {"sample_name"}
    )
    print(f"metadata columns shared by train/test metadata: {len(common_metadata_cols):,}")
    print(f"min_samples threshold for train classes: {args.min_samples:,}")
    print()

    rows = []
    for col in common_metadata_cols:
        train_counts = count_after_feature_merge(train_features, train_metadata, col)
        test_counts = count_after_feature_merge(test_features, test_metadata, col)

        train_classes_ge_min = train_counts[train_counts >= args.min_samples]
        train_classes_to_keep = set(train_classes_ge_min.index.astype(str))
        test_counts_in_train = test_counts[test_counts.index.astype(str).isin(train_classes_to_keep)]

        if len(train_classes_ge_min) < 2:
            status = "FAIL_train_lt_2_classes_ge_min"
        elif test_counts.empty:
            status = "FAIL_no_test_metadata_after_merge"
        elif test_counts_in_train.empty:
            status = "FAIL_no_test_samples_in_train_classes"
        else:
            status = "OK"

        rows.append(
            {
                "metadata_column": col,
                "status": status,
                "train_n_with_metadata": int(train_counts.sum()),
                "test_n_with_metadata": int(test_counts.sum()),
                "train_classes_total": int(len(train_counts)),
                "test_classes_total": int(len(test_counts)),
                "train_classes_ge_min": int(len(train_classes_ge_min)),
                "test_n_in_train_classes": int(test_counts_in_train.sum()),
                "train_counts": format_counts(train_counts),
                "test_counts": format_counts(test_counts),
                "test_counts_in_train_classes": format_counts(test_counts_in_train),
            }
        )

    report = pd.DataFrame(rows)
    if report.empty:
        print("No shared metadata columns were found.")
    else:
        print("=== Columns run_adelie_genomes.py can process ===")
        print(report["status"].value_counts().to_string())
        print()
        printable_cols = [
            "metadata_column",
            "status",
            "train_n_with_metadata",
            "test_n_with_metadata",
            "train_classes_ge_min",
            "test_n_in_train_classes",
        ]
        print(report[printable_cols].to_string(index=False))

    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        report.to_csv(output, sep="\t", index=False)
        print()
        print(f"Wrote diagnostic table: {output}")


if __name__ == "__main__":
    main()
