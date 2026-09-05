from sklearn.metrics import confusion_matrix, r2_score
from sklearn.preprocessing import OneHotEncoder

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np
import scipy.stats as st

import pyarrow.feather as feather
import pandas as pd

import argparse
from pathlib import Path

np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a model to predict antibiotic resistance"
    )
    # Keep original separate train/test file structure
    parser.add_argument(
        "--train_features",
        type=str,
        help="Path to the training features file",
        required=True,
    )
    parser.add_argument(
        "--train_metadata",
        type=str,
        help="Path to the training metadata file",
        required=True,
    )
    parser.add_argument(
        "--test_features",
        type=str,
        help="Path to the test features file",
        required=True,
    )
    parser.add_argument(
        "--test_metadata",
        type=str,
        help="Path to the test metadata file",
        required=True,
    )
    parser.add_argument(
        "--output_prefix", type=str, help="Prefix for the output files", required=True
    )
    parser.add_argument(
        "--min_samples",
        type=int,
        default=28,
        help="Minimum number of samples per category to keep",
    )
    parser.add_argument(
        "--n_threads",
        type=int,
        default=1,
        help="Number of threads to use for training the model",
    )

    # New parameters from the generic version
    parser.add_argument(
        "--even_samples",
        action="store_true",
        help="Keep the same number of samples per class in the training set",
    )
    parser.add_argument(
        "--grouped",
        action="store_true",
        default=False,
        help="Use grouped elastic net based on feature name prefixes",
    )
    parser.add_argument(
        "--max_iters",
        type=float,
        default=1e5,
        help="Maximum number of iterations for the Adelie model training",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-7,
        help="Tolerance for the Adelie model training convergence",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1,
        help="Alpha parameter for elastic net regularization (0 = ridge, 1 = lasso)",
    )
    return parser.parse_args()


def read_feather_data(file_path):
    return feather.read_feather(file_path)


def get_metadata_delimiter(file_path):
    suffixes = [suffix.lower() for suffix in Path(file_path).suffixes]
    if suffixes and suffixes[-1] in {".gz", ".bz2", ".xz", ".zip"}:
        suffixes = suffixes[:-1]
    suffix = suffixes[-1] if suffixes else ""
    if suffix == ".csv":
        return ","
    if suffix in {".tsv", ".tab", ".txt"}:
        return "\t"
    raise ValueError(
        f"Unsupported metadata suffix for {file_path}. "
        "Use .csv, .tsv, .tab, or .txt, optionally with compression."
    )


def read_metadata(file_path):
    # Read all metadata columns as strings to avoid dtype mismatches later
    metadata = pd.read_csv(file_path, sep=get_metadata_delimiter(file_path), dtype=str)
    if "sample_name" not in metadata.columns:
        raise ValueError(
            f"Metadata file {file_path} must contain a sample_name column; "
            f"detected columns: {list(metadata.columns)}"
        )
    # Ensure all columns are strings (defensive)
    metadata = metadata.astype(str)
    return metadata


def get_metadata_columns(metadata, min_samples=50):
    """
    Returns the columns of the metadata file except for the sample_name column
    Filters the columns so it will only return metadata that have more than two
    discrete values with greater than min_samples per category
    """
    filtered_metadata = metadata.loc[:, metadata.columns != "sample_name"]
    if min_samples > 0:
        filtered_metadata = filtered_metadata.loc[
            :, filtered_metadata.apply(lambda x: len(x.unique()) >= 2, axis=0)
        ]
        filtered_metadata = filtered_metadata.loc[
            :,
            filtered_metadata.apply(
                lambda x: sum(x.value_counts() >= min_samples) > 1, axis=0
            ),
        ]
    else:
        filtered_metadata = filtered_metadata.loc[
            :, filtered_metadata.apply(lambda x: len(x.unique()) >= 1, axis=0)
        ]
        filtered_metadata = filtered_metadata.loc[
            :,
            filtered_metadata.apply(
                lambda x: sum(x.value_counts() >= min_samples) >= 1, axis=0
            ),
        ]
    return filtered_metadata.columns


def clean_metadata_series(series):
    series = series.copy()
    if series.dtype == object:
        series = series.astype(str).str.strip()
        series = series.replace(
            {"": np.nan, "nan": np.nan, "NaN": np.nan, "NA": np.nan, "None": np.nan}
        )
    return series


def get_numeric_metadata_columns(metadata, min_samples=50):
    """Return metadata columns whose non-missing values are all numeric."""
    columns = []
    for column in metadata.columns:
        if column == "sample_name":
            continue
        values = clean_metadata_series(metadata[column])
        numeric = pd.to_numeric(values, errors="coerce")
        if values.notna().sum() != numeric.notna().sum():
            continue
        if numeric.notna().sum() >= max(2, min_samples) and numeric.nunique() >= 2:
            columns.append(column)
    return columns


def merge_data(data, metadata, metadata_col, min_samples=50, even_samples=False):
    metadata = metadata[["sample_name", metadata_col]]
    metadata[metadata_col] = metadata[metadata_col].replace("nan", pd.NA)

    merged_data = pd.merge(data, metadata, on="sample_name", how="left")
    merged_data = merged_data.dropna(subset=[metadata_col])

    class_counts = merged_data[metadata_col].value_counts()
    class_counts = class_counts[class_counts >= min_samples]
    classes_to_keep = class_counts.index
    classes_to_keep = classes_to_keep[~pd.isna(classes_to_keep)]
    classes_to_keep = classes_to_keep[classes_to_keep != "nan"]

    if len(classes_to_keep) == 0:
        return None, None, None
    if len(classes_to_keep) < 2 and min_samples != 0:
        print("This logic is true")
        return None, None, None

    merged_data = merged_data[merged_data[metadata_col].isin(classes_to_keep)]

    if even_samples:
        num_to_keep = class_counts.min()
        indices_to_keep = (
            merged_data.groupby(metadata_col)
            .apply(
                lambda x: x.sample(n=num_to_keep, replace=False).index,
                include_groups=False,
            )
            .explode()
        )
        merged_data = merged_data.loc[indices_to_keep]

    X = merged_data.drop(["sample_name", metadata_col], axis=1)
    y = merged_data[metadata_col].to_numpy()
    return np.asfortranarray(np.asarray(X, dtype=np.float64)), y, X.columns


def merge_continuous_data(data, metadata, metadata_col, min_samples=50):
    metadata = metadata[["sample_name", metadata_col]].copy()
    metadata[metadata_col] = pd.to_numeric(
        clean_metadata_series(metadata[metadata_col]), errors="coerce"
    )
    merged_data = pd.merge(data, metadata, on="sample_name", how="left")
    merged_data = merged_data.dropna(subset=[metadata_col])

    minimum = max(2, min_samples)
    if len(merged_data) < minimum or (
        min_samples != 0 and merged_data[metadata_col].nunique() < 2
    ):
        return None, None, None

    X = merged_data.drop(["sample_name", metadata_col], axis=1)
    y = merged_data[metadata_col].to_numpy(dtype=np.float64)
    return np.asfortranarray(np.asarray(X, dtype=np.float64)), y, X.columns


def get_group_ids(column_names):
    """
    Given a list of the column names for X, return a list of the starting
    index of each group based on the number following the first underscore.
    The column names are expected to be in the format [cluster|kmer]_<group>_<feature>_NUM
    """
    group_ids = []
    current_group = None
    for i, col in enumerate(column_names):
        parts = col.split("_")
        if len(parts) < 3:
            raise ValueError(
                f"Column name {col} does not have the expected format [cluster|kmer]_<group>_<feature>_..."
            )
        group = parts[1]
        if group != current_group:
            group_ids.append(i)
            current_group = group

    return np.array(group_ids, dtype=np.int32)


def remove_zero_variance_groups(X_train, X_test, column_names):
    """Remove feature groups that are constant across the training genomes."""
    X_train = np.asarray(X_train, dtype=np.float64)
    column_names = pd.Index(column_names)
    group_starts = get_group_ids(column_names)
    group_ends = np.append(group_starts[1:], X_train.shape[1])
    keep_columns = np.zeros(X_train.shape[1], dtype=bool)

    for start, end in zip(group_starts, group_ends):
        group = X_train[:, start:end]
        if not np.isfinite(group).all() or np.any(np.var(group, axis=0) > 0):
            keep_columns[start:end] = True

    if not keep_columns.any():
        raise ValueError("all feature groups have zero variance in the training data")

    removed_groups = sum(
        not keep_columns[start:end].any() for start, end in zip(group_starts, group_ends)
    )
    if removed_groups:
        print(f"Removed {removed_groups} zero-variance groups from the training genomes.")

    filtered_test = np.asfortranarray(
        np.asarray(X_test, dtype=np.float64)[:, keep_columns]
    )
    return (
        np.asfortranarray(X_train[:, keep_columns]),
        filtered_test,
        column_names[keep_columns],
    )


def train_adelie_model(
    X_train, y_train, n_threads=1, group_ids=None, max_iters=1e5, tol=1e-7, alpha=0.5
):
    oh = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    y_train2 = oh.fit_transform(y_train[:, np.newaxis])

    X_train_wrap = ad.matrix.dense(
        np.asarray(X_train, dtype=np.float64), method="naive", n_threads=n_threads
    )

    max_iters = int(max_iters)

    model = ad.GroupElasticNet(solver="cv_grpnet", family="multinomial")
    if group_ids is not None:
        model.fit(
            X_train_wrap,
            y_train2.astype(np.float64),
            n_threads=n_threads,
            groups=group_ids,
            max_iters=max_iters,
            tol=tol,
            alpha=alpha,
        )
    else:
        model.fit(
            X_train_wrap,
            y_train2.astype(np.float64),
            n_threads=n_threads,
            max_iters=max_iters,
            tol=tol,
            alpha=alpha,
        )

    return model, oh


def train_adelie_regression_model(
    X_train, y_train, n_threads=1, group_ids=None, max_iters=1e5, tol=1e-7, alpha=0.5
):
    X_train_wrap = ad.matrix.dense(
        np.asarray(X_train, dtype=np.float64), method="naive", n_threads=n_threads
    )
    fit_kwargs = {
        "n_threads": n_threads,
        "max_iters": int(max_iters),
        "tol": tol,
        "alpha": alpha,
    }
    if group_ids is not None:
        fit_kwargs["groups"] = group_ids

    model = ad.GroupElasticNet(solver="cv_grpnet", family="gaussian")
    model.fit(X_train_wrap, np.asarray(y_train, dtype=np.float64), **fit_kwargs)
    return model


def flatten_coefficients(coef):
    if hasattr(coef, "toarray"):
        return coef.toarray().flatten()
    return np.asarray(coef).flatten()


def plot_regression_predictions(y_true, y_pred, metadata_col, train_r2, test_r2):
    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    minimum = min(np.min(y_true), np.min(y_pred))
    maximum = max(np.max(y_true), np.max(y_pred))
    padding = (maximum - minimum) * 0.06 or 0.5
    limits = [minimum - padding, maximum + padding]

    ax.scatter(
        y_true,
        y_pred,
        alpha=0.78,
        s=34,
        color="#2f6f9f",
        edgecolors="white",
        linewidths=0.45,
    )
    ax.plot(limits, limits, color="#222222", linewidth=1.2, linestyle="--")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(
        f"{metadata_col}\nTest R2: {test_r2:.2f} | Train R2: {train_r2:.2f}"
    )
    ax.set_xlabel("Observed value")
    ax.set_ylabel("Predicted value")
    fig.tight_layout()
    return fig


def main():
    args = parse_args()
    output_prefix = args.output_prefix
    output_pdf = output_prefix + "_confusion_matrices.pdf"
    output_coef = output_prefix + "_nonzero_coefficients.tsv"

    train_features = read_feather_data(args.train_features)
    train_metadata = read_metadata(args.train_metadata)
    test_features = read_feather_data(args.test_features)
    test_metadata = read_metadata(args.test_metadata)

    train_metadata_columns = get_metadata_columns(
        train_metadata, min_samples=args.min_samples
    )
    test_metadata_columns = get_metadata_columns(test_metadata, min_samples=0)
    train_numeric_columns = set(
        get_numeric_metadata_columns(train_metadata, min_samples=args.min_samples)
    )
    test_numeric_columns = set(
        get_numeric_metadata_columns(test_metadata, min_samples=0)
    )
    continuous_columns = train_numeric_columns & test_numeric_columns
    categorical_columns = (
        set(train_metadata_columns) & set(test_metadata_columns)
    ) - train_numeric_columns
    metadata_columns = [
        col
        for col in train_metadata.columns
        if col in categorical_columns or col in continuous_columns
    ]

    all_model_features = None

    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing metadata column: {metadata_col}")
            print()
            continuous_target = metadata_col in continuous_columns

            if continuous_target:
                X_train, y_train, model_features = merge_continuous_data(
                    train_features,
                    train_metadata,
                    metadata_col,
                    min_samples=args.min_samples,
                )
                X_test, y_test, _ = merge_continuous_data(
                    test_features, test_metadata, metadata_col, min_samples=0
                )
            else:
                X_train, y_train, model_features = merge_data(
                    train_features,
                    train_metadata,
                    metadata_col,
                    min_samples=args.min_samples,
                    even_samples=args.even_samples,
                )
                X_test, y_test, _ = merge_data(
                    test_features, test_metadata, metadata_col, min_samples=0
                )

            if X_test is not None and not continuous_target:
                test_classes_to_keep = np.isin(y_test, np.unique(y_train))
                X_test = X_test[test_classes_to_keep]
                y_test = y_test[test_classes_to_keep]

            if X_train is None or X_test is None:
                print(
                    f"Skipping {metadata_col} as there are not enough samples after merging and filtering..."
                )
                print()
                continue

            if continuous_target:
                print(f"Fitting continuous target for {metadata_col}.")
                if args.grouped:
                    try:
                        X_train, X_test, model_features = remove_zero_variance_groups(
                            X_train, X_test, model_features
                        )
                    except ValueError as e:
                        print(f"Skipping {metadata_col}: {e}")
                        print()
                        continue
                    group_ids = get_group_ids(model_features)
                    print(f"Using grouped elastic net with {len(group_ids)} groups.")
                else:
                    print("Not using grouped elastic net.")
                    group_ids = None

                try:
                    model = train_adelie_regression_model(
                        X_train,
                        y_train,
                        n_threads=args.n_threads,
                        group_ids=group_ids,
                        tol=args.tol,
                        max_iters=args.max_iters,
                        alpha=args.alpha,
                    )
                except Exception as e:
                    print(f"Failed to train model for {metadata_col}: {e}")
                    print()
                    continue

                try:
                    yhat_train = np.asarray(
                        model.predict(X_train.astype(np.float64))
                    ).flatten()
                    yhat = np.asarray(
                        model.predict(X_test.astype(np.float64))
                    ).flatten()
                    if not np.isfinite(yhat_train).all() or not np.isfinite(yhat).all():
                        raise ValueError("model predictions contain non-finite values")
                    train_r2 = r2_score(y_train, yhat_train)
                    test_r2 = r2_score(y_test, yhat)
                    if not np.isfinite(train_r2) or not np.isfinite(test_r2):
                        raise ValueError("R2 is not finite")
                    print(f"Train R2 for {metadata_col}: {train_r2:.4f}")
                    print(f"Test R2 for {metadata_col}: {test_r2:.4f}")
                except Exception as e:
                    print(f"Failed to evaluate model for {metadata_col}: {e}")
                    print()
                    continue

                coefficient_values = flatten_coefficients(model.coef_)
                model_features_df = pd.DataFrame(
                    {"feature": model_features, "coefficient": coefficient_values}
                )
                model_features_df = model_features_df[
                    model_features_df["coefficient"] != 0
                ]
                model_features_df["metadata_category"] = metadata_col
                model_features_df["accuracy"] = test_r2
                model_features_df["train_accuracy"] = train_r2
                model_features_df["sensitivity"] = "NA"
                model_features_df["specificity"] = "NA"
                model_features_df["confusion_matrix"] = "NA"
                model_features_df["classes"] = "[residual]"
                model_features_df["coefficients"] = model_features_df[
                    "coefficient"
                ].apply(lambda value: f"[{value}]")
                model_features_df = model_features_df[
                    [
                        "metadata_category",
                        "feature",
                        "accuracy",
                        "train_accuracy",
                        "sensitivity",
                        "specificity",
                        "confusion_matrix",
                        "classes",
                        "coefficients",
                    ]
                ]

                if all_model_features is None:
                    all_model_features = model_features_df
                else:
                    all_model_features = pd.concat(
                        [all_model_features, model_features_df], axis=0
                    )

                fig = plot_regression_predictions(
                    y_test, yhat, metadata_col, train_r2, test_r2
                )
                pdf.savefig(fig)
                plt.close(fig)
                print()
                continue

            num_classes = len(np.unique(y_train))
            print(f"Number of classes for {metadata_col}: {num_classes}")

            # Set group ids based on feature names if --grouped is supplied
            if args.grouped and num_classes < 4:
                group_ids = get_group_ids(model_features)
                print(f"Using grouped elastic net with {len(group_ids)} groups.")
            else:
                print("Not using grouped elastic net.")
                if args.grouped:
                    print(
                        f"Skipping grouped elastic net for {metadata_col} as it has {num_classes} classes (must be less than 4)."
                    )
                group_ids = None

            try:
                model, oh = train_adelie_model(
                    X_train,
                    y_train,
                    n_threads=args.n_threads,
                    group_ids=group_ids,
                    tol=args.tol,
                    max_iters=args.max_iters,
                    alpha=args.alpha,
                )
            except Exception as e:
                print(f"Failed to train model for {metadata_col}: {e}")
                continue

            # Test predictions
            print(X_test)
            yhat = model.predict(X_test.astype(np.float64))
            if len(np.unique(yhat)) < 2:
                print(f"Test predictions for {metadata_col} are all of one class.")
                unique_class = np.unique(yhat)[0]
                yhat_2d = np.zeros((y_test.size, len(oh.categories_[0])))
                yhat_2d[:, unique_class] = 1
                y_pred = oh.inverse_transform(yhat_2d).flatten()
            else:
                yhat_2d = np.zeros((yhat.size, len(oh.categories_[0])))
                yhat_2d[np.arange(yhat.size), yhat] = 1
                yhat = yhat_2d
                y_pred = oh.inverse_transform(yhat).flatten()

            cm = confusion_matrix(y_test, y_pred, labels=oh.categories_[0])
            print(f"Test confusion matrix for {metadata_col}")
            print(cm)

            # Train predictions for accuracy tracking
            yhat_train = model.predict(X_train.astype(np.float64))
            if len(np.unique(yhat_train)) < 2:
                print(f"Train predictions for {metadata_col} are all of one class.")
                unique_class_train = np.unique(yhat_train)[0]
                yhat_train_2d = np.zeros((y_train.size, len(oh.categories_[0])))
                yhat_train_2d[:, unique_class_train] = 1
                y_train_pred = oh.inverse_transform(yhat_train_2d).flatten()
            else:
                yhat_train_2d = np.zeros((yhat_train.size, len(oh.categories_[0])))
                yhat_train_2d[np.arange(yhat_train.size), yhat_train] = 1
                yhat_train = yhat_train_2d
                y_train_pred = oh.inverse_transform(yhat_train).flatten()

            cm_train = confusion_matrix(y_train, y_train_pred, labels=oh.categories_[0])
            print(f"Train confusion matrix for {metadata_col}")
            print(cm_train)

            # Extract coefficients
            coef = model.coef_
            metadata_categories = oh.categories_[0]
            model_features_expanded = [
                f"{feature}+{category}"
                for feature in model_features
                for category in metadata_categories
            ]

            model_features_df = pd.DataFrame(
                model_features_expanded, columns=["feature"]
            )
            model_features_df["coefficient"] = coef.toarray().flatten()
            model_features_df = model_features_df[model_features_df["coefficient"] != 0]

            model_features_df["feature"] = model_features_df["feature"].str.split("+")
            model_features_df["category"] = model_features_df["feature"].str[1]
            model_features_df["feature"] = model_features_df["feature"].str[0]

            model_features_df = (
                model_features_df.groupby("feature")
                .agg({"coefficient": list, "category": list})
                .reset_index()
            )
            model_features_df["classes"] = model_features_df["category"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )
            model_features_df["coefficients"] = model_features_df["coefficient"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )

            model_features_df["metadata_category"] = metadata_col

            # Calculate test metrics
            if cm.shape[0] > 2:
                accuracy = np.trace(cm) / np.sum(cm)
                specificity = None
                sensitivity = None
            else:
                accuracy = (cm[0][0] + cm[1][1]) / sum(sum(row) for row in cm)
                specificity = cm[0][0] / (cm[0][0] + cm[0][1])
                sensitivity = cm[1][1] / (cm[1][0] + cm[1][1])

            # Calculate train metrics
            train_accuracy = np.trace(cm_train) / np.sum(cm_train)

            if specificity is None:
                print(f"Test Accuracy: {accuracy:.2f}")
            else:
                print(
                    f"Test Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                )

            print(f"Train Accuracy: {train_accuracy:.2f}")

            model_features_df["accuracy"] = accuracy
            model_features_df["train_accuracy"] = train_accuracy
            model_features_df["specificity"] = (
                specificity if specificity is not None else "NA"
            )
            model_features_df["sensitivity"] = (
                sensitivity if sensitivity is not None else "NA"
            )

            # Format confusion matrix for output
            out_cm = [map(str, row) for row in cm]
            out_cm = "{" + ";".join([",".join(row) for row in out_cm]) + "}"
            model_features_df["confusion_matrix"] = out_cm

            model_features_df = model_features_df[
                [
                    "metadata_category",
                    "feature",
                    "accuracy",
                    "train_accuracy",
                    "sensitivity",
                    "specificity",
                    "confusion_matrix",
                    "classes",
                    "coefficients",
                ]
            ]

            if all_model_features is None:
                all_model_features = model_features_df
            else:
                all_model_features = pd.concat(
                    [all_model_features, model_features_df], axis=0
                )

            # Plot confusion matrix
            plt.figure()
            plt.imshow(
                cm / cm.sum(axis=1)[:, np.newaxis], cmap="viridis", vmin=0, vmax=1
            )
            plt.colorbar()
            for i in range(cm.shape[0]):
                for j in range(cm.shape[1]):
                    plt.text(j, i, f"{cm[i, j]}", ha="center", va="center")
            if cm.shape[0] > 2:
                plt.title(
                    f"{metadata_col}\nTest Accuracy: {accuracy:.2f}\nTrain Accuracy: {train_accuracy:.2f}"
                )
            else:
                plt.title(
                    f"{metadata_col}\nTest Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}\nTrain Accuracy: {train_accuracy:.2f}"
                )
            plt.xlabel("Predicted")
            plt.ylabel("True")
            plt.xticks(range(cm.shape[1]), metadata_categories, rotation=45)
            plt.yticks(range(cm.shape[0]), metadata_categories, rotation=45)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

            print()

        if all_model_features is None:
            plt.figure()
            plt.text(
                0.5, 0.5, "No metadata columns were processed", ha="center", va="center"
            )
            plt.axis("off")
            pdf.savefig()
            plt.close()

            columns = [
                "metadata_category",
                "feature",
                "accuracy",
                "train_accuracy",
                "sensitivity",
                "specificity",
                "confusion_matrix",
                "classes",
                "coefficients",
            ]
            pd.DataFrame(columns=columns).to_csv(output_coef, sep="\t", index=False)
        else:
            all_model_features.to_csv(
                output_coef, sep="\t", index=False, float_format="%.4f"
            )


if __name__ == "__main__":
    main()
