# from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import OneHotEncoder

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np

# import scipy.stats as st

import pyarrow.feather as feather
import pandas as pd
from math import floor

# from os.path import basename
import argparse
from pathlib import Path
import textwrap

np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a model to predict antibiotic resistance"
    )
    parser.add_argument("--data", type=str, help="Path to the data file", required=True)
    parser.add_argument(
        "--metadata", type=str, help="Path to the metadata file", required=True
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
    parser.add_argument(
        "--balanced_test",
        action="store_true",
        help="Keep the same number of samples per class in the test set",
    )
    parser.add_argument(
        "--train_prop",
        type=float,
        default=0.5,
        help="Proportion of the data to use for training."
        "Grabs this proportion from the smallest class and then evenly samples that number from all other classes.",
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
    parser.add_argument(
        "--target_vars",
        type=str,
        default="",
        help="Semicolon-delimited target variables to residualize. Use 'all' as the last field to also fit other metadata columns normally.",
    )
    parser.add_argument(
        "--confound_vars",
        type=str,
        default="",
        help="Semicolon-delimited confounder groups corresponding to --target_vars; confounders within a group are comma-delimited.",
    )
    return parser.parse_args()


def read_feather_data(file_path):
    return feather.read_feather(file_path)


def get_metadata_delimiter(file_path):
    suffix = Path(file_path).suffix.lower()
    if suffix == ".csv":
        return ","
    return "\t"


def get_confusion_log_path(output_pdf, metadata_file):
    suffix = ".csv" if Path(metadata_file).suffix.lower() == ".csv" else ".tsv"
    return str(Path(output_pdf).with_suffix(suffix))


def read_metadata(file_path):
    metadata = pd.read_csv(file_path, sep=get_metadata_delimiter(file_path))
    if "sample_name" not in metadata.columns:
        raise ValueError("Metadata file must contain a sample_name column")
    metadata = metadata.copy()
    metadata["sample_name"] = metadata["sample_name"].astype(str)
    return metadata


def append_confusion_log_rows(
    rows,
    raw_metadata,
    metadata_col,
    matrix_name,
    y_true,
    y_pred,
    sample_names,
):
    if sample_names is None:
        return

    raw_lookup = raw_metadata.set_index("sample_name", drop=False)
    table_row = {
        "row_type": "confusion_table",
        "metadata_category": metadata_col,
        "matrix": matrix_name,
        "true_label": "",
        "predicted_label": "",
        "n_samples": len(sample_names),
    }
    rows.append(table_row)

    labels = sorted(set(map(str, y_true)) | set(map(str, y_pred)))
    for true_label in labels:
        for predicted_label in labels:
            entry_mask = (np.asarray(y_true).astype(str) == true_label) & (
                np.asarray(y_pred).astype(str) == predicted_label
            )
            entry_samples = np.asarray(sample_names).astype(str)[entry_mask]
            entry_row = {
                "row_type": "entry",
                "metadata_category": metadata_col,
                "matrix": matrix_name,
                "true_label": true_label,
                "predicted_label": predicted_label,
                "n_samples": len(entry_samples),
            }
            rows.append(entry_row)

            for sample_name in entry_samples:
                sample_row = {
                    "row_type": "sample",
                    "metadata_category": metadata_col,
                    "matrix": matrix_name,
                    "true_label": true_label,
                    "predicted_label": predicted_label,
                    "n_samples": "",
                }
                if sample_name in raw_lookup.index:
                    raw_row = raw_lookup.loc[sample_name]
                    if isinstance(raw_row, pd.DataFrame):
                        raw_row = raw_row.iloc[0]
                    sample_row.update(raw_row.to_dict())
                else:
                    sample_row["sample_name"] = sample_name
                rows.append(sample_row)


def append_regression_log_rows(
    rows,
    raw_metadata,
    metadata_col,
    matrix_name,
    y_true,
    y_pred,
    sample_names,
):
    if sample_names is None:
        return

    raw_lookup = raw_metadata.set_index("sample_name", drop=False)
    rows.append(
        {
            "row_type": "regression_table",
            "metadata_category": metadata_col,
            "matrix": matrix_name,
            "true_label": "",
            "predicted_label": "",
            "n_samples": len(sample_names),
        }
    )

    for sample_name, observed, predicted in zip(sample_names, y_true, y_pred):
        sample_row = {
            "row_type": "prediction",
            "metadata_category": metadata_col,
            "matrix": matrix_name,
            "true_label": observed,
            "predicted_label": predicted,
            "n_samples": "",
        }
        sample_name = str(sample_name)
        if sample_name in raw_lookup.index:
            raw_row = raw_lookup.loc[sample_name]
            if isinstance(raw_row, pd.DataFrame):
                raw_row = raw_row.iloc[0]
            sample_row.update(raw_row.to_dict())
        else:
            sample_row["sample_name"] = sample_name
        rows.append(sample_row)


def clean_metadata_series(series):
    series = series.copy()
    if series.dtype == object:
        series = series.astype(str).str.strip()
        series = series.replace(
            {"": np.nan, "nan": np.nan, "NaN": np.nan, "NA": np.nan, "None": np.nan}
        )
    return series


def numericize_metadata_column(metadata, column, prefix=None, drop_first=True):
    values = clean_metadata_series(metadata[column])
    numeric = pd.to_numeric(values, errors="coerce")
    nonmissing = values.notna()

    if nonmissing.sum() == numeric.notna().sum():
        return pd.DataFrame({prefix or column: numeric.astype(float)})

    dummies = pd.get_dummies(values, prefix=prefix or column, dummy_na=False)
    if drop_first and dummies.shape[1] > 1:
        dummies = dummies.iloc[:, 1:]
    dummies = dummies.astype(float)
    dummies.loc[~nonmissing, :] = np.nan
    return dummies


def numericize_target_column(metadata, column):
    values = clean_metadata_series(metadata[column])
    numeric = pd.to_numeric(values, errors="coerce")
    nonmissing = values.notna()

    if nonmissing.sum() == numeric.notna().sum():
        return {f"{column}_residual": numeric.astype(float)}

    categories = sorted(values.dropna().unique())
    if len(categories) == 2:
        mapping = {categories[0]: 0.0, categories[1]: 1.0}
        print(f"Residual target {column}: encoding {mapping}")
        return {f"{column}_residual": values.map(mapping).astype(float)}

    targets = {}
    for category in categories:
        safe_category = make_safe_label(category)
        target_name = f"{column}_residual_{safe_category}"
        targets[target_name] = (values == category).where(nonmissing, np.nan).astype(float)
    print(
        f"Residual target {column}: created {len(targets)} one-vs-rest residual targets for multiclass metadata."
    )
    return targets


def make_safe_label(value):
    label = str(value)
    label = label.replace(" ", "_").replace("/", "_").replace("\\", "_")
    label = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in label)
    label = "_".join(part for part in label.split("_") if part)
    return label or "value"


def make_unique_column_name(base_name, existing_names):
    if base_name not in existing_names:
        existing_names.add(base_name)
        return base_name
    suffix = 2
    while f"{base_name}_{suffix}" in existing_names:
        suffix += 1
    unique_name = f"{base_name}_{suffix}"
    existing_names.add(unique_name)
    return unique_name


def make_display_metadata_name(metadata_col):
    return metadata_col.replace("__", "_")


def residual_title_font_size(confound_label):
    label_length = len(confound_label or "")
    if label_length > 140:
        return 6
    if label_length > 90:
        return 7
    if label_length > 50:
        return 8
    return 9


def add_regression_plot_title(metadata_col, metric_label, confound_label=""):
    ax = plt.gca()
    display_name = make_display_metadata_name(metadata_col)
    if not confound_label:
        ax.set_title(f"{display_name}\n{metric_label}")
        return

    ax.set_title(f"{display_name}\n{metric_label}", fontsize=11, pad=30)
    wrapped_label = textwrap.fill(
        f"Residualized against: {confound_label}",
        width=115,
        break_long_words=False,
    )
    ax.text(
        0.5,
        1.01,
        wrapped_label,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=residual_title_font_size(confound_label),
    )


def plot_regression_predictions(
    y_true,
    y_pred,
    metadata_col,
    metric_label,
    confound_label="",
):
    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    min_value = min(np.nanmin(y_true), np.nanmin(y_pred))
    max_value = max(np.nanmax(y_true), np.nanmax(y_pred))
    padding = (max_value - min_value) * 0.06
    if padding == 0:
        padding = 0.5
    axis_min = min_value - padding
    axis_max = max_value + padding

    ax.scatter(
        y_true,
        y_pred,
        alpha=0.78,
        s=34,
        color="#2f6f9f",
        edgecolors="white",
        linewidths=0.45,
    )
    ax.plot(
        [axis_min, axis_max],
        [axis_min, axis_max],
        color="#222222",
        linewidth=1.2,
        linestyle="--",
    )
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.04,
        0.96,
        f"n = {len(y_true)}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#bdbdbd", alpha=0.9),
    )
    add_regression_plot_title(metadata_col, metric_label, confound_label)
    x_label = "Observed residual" if confound_label else "Observed value"
    y_label = "Predicted residual" if confound_label else "Predicted value"
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    fig.tight_layout()
    return fig


def parse_residual_options(target_vars, confound_vars):
    target_vars = str(target_vars).strip().strip("\"'").strip()
    confound_vars = str(confound_vars).strip().strip("\"'").strip()
    if not target_vars or not confound_vars:
        return [], [], False

    targets = [item.strip() for item in target_vars.split(";") if item.strip()]
    confound_groups = [item.strip() for item in confound_vars.split(";")]
    if len(confound_groups) < len(targets):
        confound_groups.extend([""] * (len(targets) - len(confound_groups)))

    residual_specs = []
    raw_targets = []
    include_all = False
    for target, confound_group in zip(targets, confound_groups):
        if target.lower() == "all":
            include_all = True
            continue
        confounds = [item.strip() for item in confound_group.split(",") if item.strip()]
        if confounds:
            residual_specs.append((target, confounds))
        else:
            raw_targets.append(target)
    return residual_specs, raw_targets, include_all


def residualize_series(target, confound_matrix):
    residual = target.astype(float).copy()
    for confound_name, confound_df in confound_matrix:
        valid = residual.notna() & confound_df.notna().all(axis=1)
        if valid.sum() < 2:
            print(f"Skipping residual adjustment for {confound_name}: fewer than 2 complete samples.")
            continue

        x = confound_df.loc[valid].to_numpy(dtype=np.float64)
        x = np.column_stack([np.ones(x.shape[0]), x])
        y = residual.loc[valid].to_numpy(dtype=np.float64)
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        fitted = x @ beta
        residual.loc[valid] = y - fitted
        residual.loc[~valid] = np.nan
    return residual


def resolve_confound_columns(metadata, target_col, confound_cols):
    if any(confound_col.lower() == "all" for confound_col in confound_cols):
        resolved = [
            col
            for col in metadata.columns
            if col != "sample_name"
            and col != target_col
            and "_residual" not in col.replace("__", "_")
        ]
        print(
            f"Residual target {target_col}: expanding confounder 'all' to {len(resolved)} metadata columns."
        )
        return resolved
    return confound_cols


def add_residual_targets(metadata, target_vars, confound_vars):
    residual_specs, raw_targets, include_all = parse_residual_options(
        target_vars, confound_vars
    )
    if not residual_specs:
        return metadata, {}, raw_targets, include_all, {}, {}

    metadata = metadata.copy()
    residual_columns_by_target = {}
    residual_confound_labels = {}
    residual_split_columns = {}
    existing_residual_names = set(metadata.columns)
    for spec_index, (target_col, confound_cols) in enumerate(residual_specs, start=1):
        if target_col not in metadata.columns:
            print(f"Skipping residual target {target_col}: column not found in metadata.")
            continue
        if not confound_cols:
            print(f"Skipping residual target {target_col}: no confounders were provided.")
            continue

        confound_cols = resolve_confound_columns(metadata, target_col, confound_cols)
        confound_label = ", ".join(confound_cols)
        confound_matrix = []
        missing_confounds = []
        for confound_col in confound_cols:
            if confound_col not in metadata.columns:
                missing_confounds.append(confound_col)
                continue
            confound_matrix.append(
                (
                    confound_col,
                    numericize_metadata_column(
                        metadata, confound_col, prefix=f"confound__{confound_col}"
                    ),
                )
            )

        if missing_confounds:
            print(
                f"Residual target {target_col}: missing confounder columns ignored: {', '.join(missing_confounds)}"
            )
        if not confound_matrix:
            print(f"Skipping residual target {target_col}: no valid confounders remain.")
            continue

        residual_columns_by_target.setdefault(target_col, [])
        for residual_col, numeric_target in numericize_target_column(metadata, target_col).items():
            base_residual_col = residual_col
            if len([spec for spec in residual_specs if spec[0] == target_col]) > 1:
                base_residual_col = f"{residual_col}_adjustment{spec_index}"
            unique_residual_col = make_unique_column_name(
                base_residual_col, existing_residual_names
            )
            metadata[unique_residual_col] = residualize_series(
                numeric_target, confound_matrix
            )
            residual_columns_by_target[target_col].append(unique_residual_col)
            residual_confound_labels[unique_residual_col] = confound_label
            residual_split_columns[unique_residual_col] = target_col
            complete = metadata[unique_residual_col].notna().sum()
            print(
                f"Created residual target {unique_residual_col} from {target_col} after adjusting for {', '.join(confound_cols)} ({complete} complete samples)."
            )

    return (
        metadata,
        residual_columns_by_target,
        raw_targets,
        include_all,
        residual_confound_labels,
        residual_split_columns,
    )


def get_metadata_columns(metadata, min_samples=50):
    """
    Returns the columns of the metadata file except for the sample_name column
    Filters the columns so it will only return metadata that have more than two
    discrete values with greater than min_samples per category
    """
    filtered_metadata = metadata.loc[:, metadata.columns != "sample_name"].copy()
    filtered_metadata = filtered_metadata.apply(clean_metadata_series)
    # filter out columns with less than 2 unique values
    filtered_metadata = filtered_metadata.loc[
        :, filtered_metadata.apply(lambda x: len(x.unique()) >= 2, axis=0)
    ]
    # only grab columns with two or more categories that have more than min_samples
    filtered_metadata = filtered_metadata.loc[
        :,
        filtered_metadata.apply(
            lambda x: sum(x.value_counts() > min_samples) > 1, axis=0
        ),
    ]
    return filtered_metadata.columns


def get_numeric_metadata_columns(metadata, min_samples=50):
    columns = []
    for column in metadata.columns:
        if column == "sample_name":
            continue
        values = clean_metadata_series(metadata[column])
        numeric = pd.to_numeric(values, errors="coerce")
        if values.notna().sum() != numeric.notna().sum():
            continue
        if numeric.notna().sum() >= max(2, min_samples) and numeric.nunique(dropna=True) >= 2:
            columns.append(column)
    return columns


def merge_and_split_data(
    data,
    metadata,
    metadata_col,
    min_samples=50,
    train_prop=0.5,
    balanced_test=False,
    continuous=False,
    continuous_split_col=None,
):
    metadata_cols = ["sample_name", metadata_col]
    if continuous_split_col and continuous_split_col in metadata.columns:
        metadata_cols.append(continuous_split_col)
    metadata = metadata[metadata_cols]
    merged_data = pd.merge(data, metadata, on="sample_name", how="left")
    merged_data = merged_data.dropna(subset=[metadata_col])

    if continuous:
        split_col = continuous_split_col if continuous_split_col in merged_data.columns else None
        if split_col:
            merged_data = merged_data.dropna(subset=[split_col])
            split_values = clean_metadata_series(merged_data[split_col])
            split_numeric = pd.to_numeric(split_values, errors="coerce")
            use_categorical_split = split_values.notna().sum() != split_numeric.notna().sum()
            if not use_categorical_split and split_numeric.nunique(dropna=True) <= 10:
                use_categorical_split = True

            if use_categorical_split:
                class_counts = split_values.value_counts()
                class_counts = class_counts[class_counts >= min_samples]
                classes_to_keep = class_counts.index
                classes_to_keep = classes_to_keep[~pd.isna(classes_to_keep)]
                classes_to_keep = classes_to_keep[classes_to_keep != "nan"]
                if len(classes_to_keep) < 2:
                    return None, None, None, None, None, None, None
                merged_data = merged_data[split_values.isin(classes_to_keep)].copy()
                split_values = clean_metadata_series(merged_data[split_col])
                num_to_keep = floor(class_counts.min() * train_prop)
                if train_prop < 1 and num_to_keep < 1:
                    return None, None, None, None, None, None, None
                if train_prop == 1:
                    indices_to_keep = merged_data.index
                else:
                    indices_to_keep = (
                        merged_data.assign(__split_label=split_values)
                        .groupby("__split_label")
                        .apply(
                            lambda x: x.sample(n=num_to_keep, replace=False).index,
                            include_groups=False,
                        )
                        .explode()
                    )

                drop_cols = ["sample_name", metadata_col]
                if split_col not in drop_cols:
                    drop_cols.append(split_col)
                X_train = merged_data.drop(drop_cols, axis=1).loc[indices_to_keep]
                model_features = X_train.columns
                y_train = merged_data[metadata_col].loc[indices_to_keep].astype(float).to_numpy()
                train_sample_names = merged_data["sample_name"].loc[indices_to_keep].to_numpy()
                if train_prop == 1:
                    return (
                        np.asfortranarray(np.asarray(X_train, dtype=np.float64)),
                        None,
                        y_train,
                        None,
                        model_features,
                        train_sample_names,
                        None,
                    )
                X_test = merged_data.drop(drop_cols, axis=1).drop(indices_to_keep)
                y_test = merged_data[metadata_col].drop(indices_to_keep).astype(float).to_numpy()
                test_sample_names = merged_data["sample_name"].drop(indices_to_keep).to_numpy()
                return (
                    np.asfortranarray(np.asarray(X_train, dtype=np.float64)),
                    np.asfortranarray(np.asarray(X_test, dtype=np.float64)),
                    y_train,
                    y_test,
                    model_features,
                    train_sample_names,
                    test_sample_names,
                )

        if len(merged_data) < max(2, min_samples):
            return None, None, None, None, None, None, None
        train_size = len(merged_data) if train_prop == 1 else floor(len(merged_data) * train_prop)
        if train_prop < 1:
            train_size = min(max(train_size, 1), len(merged_data) - 1)
        indices_to_keep = merged_data.sample(n=train_size, replace=False).index
        drop_cols = ["sample_name", metadata_col]
        if split_col and split_col not in drop_cols:
            drop_cols.append(split_col)
        X_train = merged_data.drop(drop_cols, axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].astype(float).to_numpy()
        train_sample_names = merged_data["sample_name"].loc[indices_to_keep].to_numpy()
        if train_prop == 1:
            return (
                np.asfortranarray(X_train),
                None,
                y_train,
                None,
                model_features,
                train_sample_names,
                None,
            )
        X_test = merged_data.drop(drop_cols, axis=1).drop(
            indices_to_keep
        )
        y_test = merged_data[metadata_col].drop(indices_to_keep).astype(float).to_numpy()
        test_sample_names = merged_data["sample_name"].drop(indices_to_keep).to_numpy()
        return (
            np.asfortranarray(np.asarray(X_train, dtype=np.float64)),
            np.asfortranarray(np.asarray(X_test, dtype=np.float64)),
            y_train,
            y_test,
            model_features,
            train_sample_names,
            test_sample_names,
        )

    # Check the distribution of classes for this metadata category
    class_counts = merged_data[metadata_col].value_counts()

    # Drop any classes with less than min_samples
    class_counts = class_counts[class_counts >= min_samples]
    classes_to_keep = class_counts.index
    classes_to_keep = classes_to_keep[~pd.isna(classes_to_keep)]
    classes_to_keep = classes_to_keep[classes_to_keep != "nan"]
    # if there are not >= 2 classes with >= min_samples samples, return None
    if len(classes_to_keep) < 2:
        return None, None, None, None, None, None, None
    merged_data = merged_data[merged_data[metadata_col].isin(classes_to_keep)]

    # Get the minimum number of samples per class
    # keep exactly half of the samples for each class for the training set
    # and keep the rest of the samples for the test set
    num_to_keep = class_counts.min()
    if pd.isna(num_to_keep):
        return None, None, None, None, None, None, None
    num_to_keep = floor(num_to_keep * train_prop)
    indices_to_keep = (
        merged_data.groupby(metadata_col)
        .apply(
            lambda x: x.sample(n=num_to_keep, replace=False).index, include_groups=False
        )
        .explode()
    )

    if train_prop == 1:
        # Split the data into training and test sets
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()
        train_sample_names = merged_data["sample_name"].loc[indices_to_keep].to_numpy()

        return (
            np.asfortranarray(X_train),
            None,
            y_train,
            None,
            model_features,
            train_sample_names,
            None,
        )

    # If we want a balanced test set, keep the same number of samples per class in the test set
    if balanced_test:
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()
        train_sample_names = merged_data["sample_name"].loc[indices_to_keep].to_numpy()

        test_indices = (
            merged_data.drop(indices_to_keep)
            .groupby(metadata_col)
            .apply(
                lambda x: x.sample(n=num_to_keep, replace=False).index,
                include_groups=False,
            )
            .explode()
        )
        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            test_indices
        ]
        y_test = merged_data[metadata_col].loc[test_indices].to_numpy()
        test_sample_names = merged_data["sample_name"].loc[test_indices].to_numpy()
    else:
        # Split the data into training and test sets
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()
        train_sample_names = merged_data["sample_name"].loc[indices_to_keep].to_numpy()

        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).drop(
            indices_to_keep
        )
        y_test = merged_data[metadata_col].drop(indices_to_keep).to_numpy()
        test_sample_names = merged_data["sample_name"].drop(indices_to_keep).to_numpy()

    return (
        np.asfortranarray(np.asarray(X_train, dtype=np.float64)),
        np.asfortranarray(np.asarray(X_test, dtype=np.float64)),
        y_train,
        y_test,
        model_features,
        train_sample_names,
        test_sample_names,
    )


def get_group_ids(column_names):
    """
    Given a list of the column names for X, return a list of the starting
    index of each group based on the number following the first underscore.
    The column names are expected to be in the format [cluster|kmer]_<group>_<feature>_NUM

    Note that the column names must be sorted such that all features from the same group
    are together. This is the case for the current implementation of feature generation.

    Should return an ndarry of these starting indices.
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
    """Remove feature groups that are constant across the training samples."""
    X_train = np.asarray(X_train, dtype=np.float64)
    column_names = pd.Index(column_names)
    group_starts = get_group_ids(column_names)
    group_ends = np.append(group_starts[1:], X_train.shape[1])
    keep_columns = np.zeros(X_train.shape[1], dtype=bool)
    removed_groups = 0

    for start, end in zip(group_starts, group_ends):
        group = X_train[:, start:end]
        if not np.isfinite(group).all():
            keep_columns[start:end] = True
            continue
        if np.any(np.var(group, axis=0) > 0):
            keep_columns[start:end] = True
        else:
            removed_groups += 1

    if not keep_columns.any():
        raise ValueError("all feature groups have zero variance in the training data")

    removed_features = int((~keep_columns).sum())
    if removed_groups:
        print(
            f"Removed {removed_groups} zero-variance groups "
            f"({removed_features} features) from the training subset."
        )

    filtered_train = np.asfortranarray(X_train[:, keep_columns])
    filtered_test = None
    if X_test is not None:
        filtered_test = np.asfortranarray(
            np.asarray(X_test, dtype=np.float64)[:, keep_columns]
        )
    return filtered_train, filtered_test, column_names[keep_columns]


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

    max_iters = int(max_iters)
    model = ad.GroupElasticNet(solver="cv_grpnet", family="gaussian")
    fit_kwargs = dict(
        n_threads=n_threads,
        max_iters=max_iters,
        tol=tol,
        alpha=alpha,
    )
    if group_ids is not None:
        fit_kwargs["groups"] = group_ids
    model.fit(X_train_wrap, np.asarray(y_train, dtype=np.float64), **fit_kwargs)
    return model


def flatten_coefficients(coef):
    if hasattr(coef, "toarray"):
        return coef.toarray().flatten()
    return np.asarray(coef).flatten()


def main():
    args = parse_args()
    output_prefix = args.output_prefix
    output_pdf = output_prefix + "_confusion_matrices.pdf"
    output_coef = output_prefix + "_nonzero_coefficients.tsv"
    output_confusion_log = get_confusion_log_path(output_pdf, args.metadata)

    # Load teh data and metadata
    data = read_feather_data(args.data)
    metadata = read_metadata(args.metadata)
    raw_metadata = metadata.copy()
    (
        metadata,
        residual_columns_by_target,
        raw_target_columns,
        include_all_metadata,
        residual_confound_labels,
        residual_split_columns,
    ) = add_residual_targets(metadata, args.target_vars, args.confound_vars)
    # Get the metadata columns that have more than 2 unique values
    # and more than 50 samples per category
    original_metadata = metadata.drop(
        columns=[
            col
            for residual_cols in residual_columns_by_target.values()
            for col in residual_cols
        ],
        errors="ignore",
    )
    metadata_columns = list(
        get_metadata_columns(original_metadata, min_samples=args.min_samples)
    )
    numeric_metadata_columns = get_numeric_metadata_columns(
        original_metadata, min_samples=args.min_samples
    )
    numeric_metadata_column_set = set(numeric_metadata_columns)
    metadata_columns = [
        col for col in metadata_columns if col not in numeric_metadata_column_set
    ]
    continuous_metadata_columns = {
        col for residual_cols in residual_columns_by_target.values() for col in residual_cols
    }
    continuous_metadata_columns.update(numeric_metadata_column_set)
    if residual_columns_by_target or raw_target_columns or include_all_metadata:
        residual_targets = set(residual_columns_by_target)
        raw_targets_seen = set()
        selected_columns = []
        for target_col, residual_cols in residual_columns_by_target.items():
            selected_columns.extend(residual_cols)
        for target_col in raw_target_columns:
            if target_col not in metadata.columns:
                print(f"Skipping raw target {target_col}: column not found in metadata.")
                continue
            if target_col in raw_targets_seen:
                continue
            selected_columns.append(target_col)
            raw_targets_seen.add(target_col)
        if include_all_metadata:
            selected_columns.extend(
                [
                    col
                    for col in metadata_columns + numeric_metadata_columns
                    if col not in residual_targets and col not in raw_targets_seen
                ]
            )
        metadata_columns = selected_columns
    else:
        metadata_columns = metadata_columns + numeric_metadata_columns

    all_model_features = None
    confusion_log_rows = []

    # Iterate over the metadata columns
    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing metadata column: {metadata_col}")
            print()
            continuous_target = metadata_col in continuous_metadata_columns

            (
                X_train,
                X_test,
                y_train,
                y_test,
                model_features,
                train_sample_names,
                test_sample_names,
            ) = merge_and_split_data(
                data,
                metadata,
                metadata_col,
                min_samples=args.min_samples,
                balanced_test=args.balanced_test,
                train_prop=args.train_prop,
                continuous=continuous_target,
                continuous_split_col=residual_split_columns.get(metadata_col),
            )

            # skip the column if the merge and split function returns None
            if X_train is None:
                print(
                    f"Skipping {metadata_col} as there are not enough samples after merging and filtering..."
                )
                print()
                continue

            if continuous_target:
                print(f"Fitting continuous target for {metadata_col}.")
                confound_label = residual_confound_labels.get(metadata_col, "")
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
                    continue

                try:
                    yhat_train = np.asarray(
                        model.predict(X_train.astype(np.float64))
                    ).flatten()
                    train_r2 = r2_score(y_train, yhat_train)
                    print(f"Train R2 for {metadata_col}: {train_r2:.4f}")
                    append_regression_log_rows(
                        confusion_log_rows,
                        raw_metadata,
                        metadata_col,
                        "train",
                        y_train,
                        yhat_train,
                        train_sample_names,
                    )

                    if args.train_prop == 1:
                        test_r2 = None
                        print(
                            "Not calculating test R2 as we are not using test data..."
                        )
                    else:
                        yhat = np.asarray(
                            model.predict(X_test.astype(np.float64))
                        ).flatten()
                        test_r2 = r2_score(y_test, yhat)
                        print(f"Test R2 for {metadata_col}: {test_r2:.4f}")
                        append_regression_log_rows(
                            confusion_log_rows,
                            raw_metadata,
                            metadata_col,
                            "test",
                            y_test,
                            yhat,
                            test_sample_names,
                        )
                except Exception as e:
                    print(f"Failed to evaluate model for {metadata_col}: {e}")
                    print()
                    continue

                coef = flatten_coefficients(model.coef_)
                model_features = pd.DataFrame(
                    {"feature": model_features, "coefficient": coef}
                )
                model_features = model_features[model_features["coefficient"] != 0]
                model_features["metadata_category"] = metadata_col
                model_features["accuracy"] = test_r2 if test_r2 is not None else "NA"
                model_features["train_accuracy"] = train_r2
                model_features["sensitivity"] = "NA"
                model_features["specificity"] = "NA"
                model_features["confusion_matrix"] = "NA"
                model_features["classes"] = "[residual]"
                model_features["coefficients"] = model_features["coefficient"].apply(
                    lambda x: f"[{x}]"
                )
                model_features = model_features[
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
                    all_model_features = model_features
                else:
                    all_model_features = pd.concat(
                        [all_model_features, model_features], axis=0
                    )

                if args.train_prop < 1:
                    fig = plot_regression_predictions(
                        y_test,
                        yhat,
                        metadata_col,
                        f"Test R2: {test_r2:.2f}",
                        confound_label,
                    )
                    pdf.savefig(fig)
                    plt.close(fig)

                print()
                continue

            num_classes = len(np.unique(y_train))
            print(f"Number of classes for {metadata_col}: {num_classes}")

            # set group ids based on feature names if --grouped is supplied
            if args.grouped and num_classes < 4:
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

            if args.train_prop == 1:
                cm = []
                print("Not calculating test accuracy as we are not using test data...")
            else:
                # add a check to make sure there are more than 2 unique values in the predictions
                # an error can be thrown if inverse_transform gets the wrong number of columns
                # Handle cases where predictions are all of one class
                yhat = model.predict(X_test.astype(np.float64))
                if len(np.unique(yhat)) < 2:
                    print(f"Predictions for {metadata_col} are all of one class.")
                    unique_class = np.unique(yhat)[0]
                    yhat_2d = np.zeros((y_test.size, len(oh.categories_[0])))
                    yhat_2d[:, unique_class] = 1
                    y_pred = oh.inverse_transform(yhat_2d).flatten()
                    cm = confusion_matrix(y_test, y_pred)
                    append_confusion_log_rows(
                        confusion_log_rows,
                        raw_metadata,
                        metadata_col,
                        "test",
                        y_test,
                        y_pred,
                        test_sample_names,
                    )
                    accuracy = np.trace(cm) / np.sum(cm)
                    print(f"Test confusion matrix for {metadata_col}")
                    print(cm)
                    print(f"Accuracy: {accuracy:.2f}")
                else:
                    yhat_2d = np.zeros((yhat.size, len(oh.categories_[0])))
                    yhat_2d[np.arange(yhat.size), yhat] = 1
                    yhat = yhat_2d

                    try:
                        y_pred = oh.inverse_transform(yhat).flatten()
                        cm = confusion_matrix(y_test, y_pred, labels=oh.categories_[0])
                        append_confusion_log_rows(
                            confusion_log_rows,
                            raw_metadata,
                            metadata_col,
                            "test",
                            y_test,
                            y_pred,
                            test_sample_names,
                        )
                        print(f"Test confusion matrix for {metadata_col}")
                        print(cm)
                        accuracy = np.trace(cm) / np.sum(cm)
                        print(f"Accuracy: {accuracy:.2f}")
                    except Exception as e:
                        print(
                            f"Failed to transform predictions for {metadata_col}: {e}"
                        )
                        continue

            yhat_train = model.predict(X_train.astype(np.float64))
            if len(np.unique(yhat_train)) < 2:
                print(f"Train predictions for {metadata_col} are all of one class.")
                unique_class_train = np.unique(yhat_train)[0]
                yhat_train_2d = np.zeros((y_train.size, len(oh.categories_[0])))
                yhat_train_2d[:, unique_class_train] = 1
                y_train_pred = oh.inverse_transform(yhat_train_2d).flatten()
                cm_train = confusion_matrix(
                    y_train, y_train_pred, labels=oh.categories_[0]
                )
                append_confusion_log_rows(
                    confusion_log_rows,
                    raw_metadata,
                    metadata_col,
                    "train",
                    y_train,
                    y_train_pred,
                    train_sample_names,
                )
                train_accuracy = np.trace(cm_train) / np.sum(cm_train)
                print(f"Train confusion matrix for {metadata_col}")
                print(cm_train)
                print(f"Train accuracy: {train_accuracy:.2f}\n")
            else:
                yhat_train_2d = np.zeros((yhat_train.size, len(oh.categories_[0])))
                yhat_train_2d[np.arange(yhat_train.size), yhat_train] = 1
                yhat_train = yhat_train_2d

                try:
                    y_train_pred = oh.inverse_transform(yhat_train).flatten()
                    cm_train = confusion_matrix(y_train, y_train_pred)
                    append_confusion_log_rows(
                        confusion_log_rows,
                        raw_metadata,
                        metadata_col,
                        "train",
                        y_train,
                        y_train_pred,
                        train_sample_names,
                    )
                    print(f"Train confusion matrix for {metadata_col}")
                    print(cm_train)
                    train_accuracy = np.trace(cm_train) / np.sum(cm_train)
                    print(f"Train accuracy: {train_accuracy:.2f}\n")
                except Exception as e:
                    cm_train = []
                    print(
                        f"Failed to transform train predictions for {metadata_col}: {e}"
                    )
                    train_accuracy = 0
                    continue

            # extract the nonzero coefficients
            coef = model.coef_
            # get the feature names for the nonzero coefficients
            metadata_categories = oh.categories_[0]
            model_features = [
                f"{feature}+{category}"
                for feature in model_features
                for category in metadata_categories
            ]

            # get the names and values of the nonzero coefficients
            model_features = pd.DataFrame(model_features, columns=["feature"])
            model_features["coefficient"] = flatten_coefficients(coef)
            model_features = model_features[model_features["coefficient"] != 0]

            # separate the feature names into the feature and category
            model_features["feature"] = model_features["feature"].str.split("+")
            model_features["category"] = model_features["feature"].str[1]
            model_features["feature"] = model_features["feature"].str[0]

            # gather the coefficients and categorie names into two columns grouping by feature
            model_features = (
                model_features.groupby("feature")
                .agg({"coefficient": list, "category": list})
                .reset_index()
            )
            model_features["classes"] = model_features["category"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )
            model_features["coefficients"] = model_features["coefficient"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )

            # add column for metadata category
            model_features["metadata_category"] = metadata_col

            # assuming cm can be larger than 2x2
            if args.train_prop < 1:
                if cm.shape[0] > 2:
                    accuracy = np.trace(cm) / np.sum(cm)
                    specificity = None
                    sensitivity = None
                else:
                    accuracy = (cm[0][0] + cm[1][1]) / sum(sum(row) for row in cm)
                    specificity = cm[0][0] / (cm[0][0] + cm[0][1])
                    sensitivity = cm[1][1] / (cm[1][0] + cm[1][1])
            else:
                if cm_train.shape[0] > 2:
                    accuracy = None
                    specificity = None
                    sensitivity = None
                else:
                    accuracy = None
                    specificity = cm_train[0][0] / (cm_train[0][0] + cm_train[0][1])
                    sensitivity = cm_train[1][1] / (cm_train[1][0] + cm_train[1][1])

            if args.train_prop < 1:
                if specificity is None:
                    print(f"Accuracy: {accuracy:.2f}")
                else:
                    print(
                        f"Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                    )
            else:
                if specificity is None:
                    print(f"Train accuracy: {train_accuracy:.2f}")
                else:
                    print(
                        f"Train Accuracy: {train_accuracy:.2f}, Train Specificity: {specificity:.2f}, Train Sensitivity: {sensitivity:.2f}"
                    )

            # add the accuracy, specificity, and sensitivity to the model features
            model_features["accuracy"] = accuracy if accuracy is not None else "NA"
            model_features["train_accuracy"] = (
                train_accuracy if train_accuracy is not None else "NA"
            )
            model_features["specificity"] = (
                specificity if specificity is not None else "NA"
            )
            model_features["sensitivity"] = (
                sensitivity if sensitivity is not None else "NA"
            )

            if args.train_prop < 1:
                out_cm = [map(str, row) for row in cm]
                out_cm = "{" + ";".join([",".join(row) for row in out_cm]) + "}"
            else:
                out_cm = None
            model_features["confusion_matrix"] = out_cm if out_cm is not None else "NA"

            model_features = model_features[
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

            # join with the larger set of model features
            if all_model_features is None:
                all_model_features = model_features
            else:
                all_model_features = pd.concat(
                    [all_model_features, model_features], axis=0
                )

            # plot the confusion matrix and save to the pdf
            # color by relative frequency
            # add numbers to the cells of the confusion matrix
            # add accuracy, specificity, and sensitivity to the title
            if args.train_prop < 1:
                plt.figure()
                plt.imshow(
                    cm / cm.sum(axis=1)[:, np.newaxis], cmap="viridis", vmin=0, vmax=1
                )
                plt.colorbar()
                for i in range(cm.shape[0]):
                    for j in range(cm.shape[1]):
                        plt.text(j, i, f"{cm[i, j]}", ha="center", va="center")
                if cm.shape[0] > 2:
                    plt.title(f"{metadata_col}\nAccuracy: {accuracy:.2f}")
                else:
                    plt.title(
                        f"{metadata_col}\nAccuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                    )
                plt.xlabel("Predicted")
                plt.ylabel("True")
                plt.xticks(range(cm.shape[1]), metadata_categories, rotation=45)
                # these need to start from the bottom
                plt.yticks(range(cm.shape[0]), metadata_categories, rotation=45)
                plt.tight_layout()
                pdf.savefig()
                plt.close()

            # add a blank line
            print()

        if all_model_features is None:
            # Write a blank page to the PDF
            plt.figure()
            plt.text(
                0.5, 0.5, "No metadata columns were processed", ha="center", va="center"
            )
            plt.axis("off")
            pdf.savefig()
            plt.close()

            # Write an empty output TSV with just column names
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
            # output the nonzero coefficients to a tsv file
            all_model_features.to_csv(
                output_coef, sep="\t", index=False, float_format="%.4f"
            )

    confusion_log_columns = [
        "row_type",
        "metadata_category",
        "matrix",
        "true_label",
        "predicted_label",
        "n_samples",
    ] + [col for col in raw_metadata.columns if col not in {
        "row_type",
        "metadata_category",
        "matrix",
        "true_label",
        "predicted_label",
        "n_samples",
    }]
    confusion_log = pd.DataFrame(confusion_log_rows)
    confusion_log = confusion_log.reindex(columns=confusion_log_columns)
    confusion_log.to_csv(
        output_confusion_log,
        sep=get_metadata_delimiter(args.metadata),
        index=False,
    )


if __name__ == "__main__":
    main()
