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


def read_metadata(file_path):
    metadata = pd.read_table(file_path)
    if "sample_name" not in metadata.columns:
        raise ValueError("Metadata file must contain a sample_name column")
    metadata = metadata.copy()
    metadata["sample_name"] = metadata["sample_name"].astype(str)
    return metadata


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
        return {f"{column}__residual": numeric.astype(float)}

    categories = sorted(values.dropna().unique())
    if len(categories) == 2:
        mapping = {categories[0]: 0.0, categories[1]: 1.0}
        print(f"Residual target {column}: encoding {mapping}")
        return {f"{column}__residual": values.map(mapping).astype(float)}

    targets = {}
    for category in categories:
        safe_category = str(category).replace(" ", "_").replace("/", "_")
        target_name = f"{column}__residual__{safe_category}"
        targets[target_name] = (values == category).where(nonmissing, np.nan).astype(float)
    print(
        f"Residual target {column}: created {len(targets)} one-vs-rest residual targets for multiclass metadata."
    )
    return targets


def parse_residual_options(target_vars, confound_vars):
    target_vars = str(target_vars).strip().strip("\"'").strip()
    confound_vars = str(confound_vars).strip().strip("\"'").strip()
    if not target_vars or not confound_vars:
        return [], False

    targets = [item.strip() for item in target_vars.split(";") if item.strip()]
    confound_groups = [item.strip() for item in confound_vars.split(";")]
    if len(confound_groups) < len(targets):
        confound_groups.extend([""] * (len(targets) - len(confound_groups)))

    residual_specs = []
    include_all = False
    for target, confound_group in zip(targets, confound_groups):
        if target.lower() == "all":
            include_all = True
            continue
        confounds = [item.strip() for item in confound_group.split(",") if item.strip()]
        residual_specs.append((target, confounds))
    return residual_specs, include_all


def residualize_series(target, confound_matrix):
    residual = target.astype(float).copy()
    for confound_name, confound_df in confound_matrix:
        valid = residual.notna() & confound_df.notna().all(axis=1)
        if valid.sum() < 2:
            print(f"Skipping residual adjustment for {confound_name}: fewer than 2 complete samples.")
            residual.loc[~valid] = np.nan
            continue

        x = confound_df.loc[valid].to_numpy(dtype=np.float64)
        x = np.column_stack([np.ones(x.shape[0]), x])
        y = residual.loc[valid].to_numpy(dtype=np.float64)
        beta, *_ = np.linalg.lstsq(x, y, rcond=None)
        fitted = x @ beta
        residual.loc[valid] = y - fitted
        residual.loc[~valid] = np.nan
    return residual


def add_residual_targets(metadata, target_vars, confound_vars):
    residual_specs, include_all = parse_residual_options(target_vars, confound_vars)
    if not residual_specs:
        return metadata, {}, include_all

    metadata = metadata.copy()
    residual_columns_by_target = {}
    for target_col, confound_cols in residual_specs:
        if target_col not in metadata.columns:
            print(f"Skipping residual target {target_col}: column not found in metadata.")
            continue
        if not confound_cols:
            print(f"Skipping residual target {target_col}: no confounders were provided.")
            continue

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

        residual_columns_by_target[target_col] = []
        for residual_col, numeric_target in numericize_target_column(metadata, target_col).items():
            metadata[residual_col] = residualize_series(numeric_target, confound_matrix)
            residual_columns_by_target[target_col].append(residual_col)
            complete = metadata[residual_col].notna().sum()
            print(
                f"Created residual target {residual_col} from {target_col} after adjusting for {', '.join(confound_cols)} ({complete} complete samples)."
            )

    return metadata, residual_columns_by_target, include_all


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


def merge_and_split_data(
    data,
    metadata,
    metadata_col,
    min_samples=50,
    train_prop=0.5,
    balanced_test=False,
    continuous=False,
):
    metadata = metadata[["sample_name", metadata_col]]
    merged_data = pd.merge(data, metadata, on="sample_name", how="left")
    merged_data = merged_data.dropna(subset=[metadata_col])

    if continuous:
        if len(merged_data) < max(2, min_samples):
            return None, None, None, None, None
        train_size = len(merged_data) if train_prop == 1 else floor(len(merged_data) * train_prop)
        if train_prop < 1:
            train_size = min(max(train_size, 1), len(merged_data) - 1)
        indices_to_keep = merged_data.sample(n=train_size, replace=False).index
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].astype(float).to_numpy()
        if train_prop == 1:
            return np.asfortranarray(X_train), None, y_train, None, model_features
        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).drop(
            indices_to_keep
        )
        y_test = merged_data[metadata_col].drop(indices_to_keep).astype(float).to_numpy()
        return (
            np.asfortranarray(np.asarray(X_train, dtype=np.float64)),
            np.asfortranarray(np.asarray(X_test, dtype=np.float64)),
            y_train,
            y_test,
            model_features,
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
        return None, None, None, None, None
    merged_data = merged_data[merged_data[metadata_col].isin(classes_to_keep)]

    # Get the minimum number of samples per class
    # keep exactly half of the samples for each class for the training set
    # and keep the rest of the samples for the test set
    num_to_keep = class_counts.min()
    if pd.isna(num_to_keep):
        return None, None, None, None, None
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

        return np.asfortranarray(X_train), None, y_train, None, model_features

    # If we want a balanced test set, keep the same number of samples per class in the test set
    if balanced_test:
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()

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
    else:
        # Split the data into training and test sets
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[
            indices_to_keep
        ]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()

        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).drop(
            indices_to_keep
        )
        y_test = merged_data[metadata_col].drop(indices_to_keep).to_numpy()

    return (
        np.asfortranarray(np.asarray(X_train, dtype=np.float64)),
        np.asfortranarray(np.asarray(X_test, dtype=np.float64)),
        y_train,
        y_test,
        model_features,
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

    # Load teh data and metadata
    data = read_feather_data(args.data)
    metadata = read_metadata(args.metadata)
    metadata, residual_columns_by_target, include_all_metadata = add_residual_targets(
        metadata, args.target_vars, args.confound_vars
    )
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
    continuous_metadata_columns = {
        col for residual_cols in residual_columns_by_target.values() for col in residual_cols
    }
    if residual_columns_by_target:
        residual_targets = set(residual_columns_by_target)
        selected_columns = []
        for target_col, residual_cols in residual_columns_by_target.items():
            selected_columns.extend(residual_cols)
        if include_all_metadata:
            selected_columns.extend(
                [col for col in metadata_columns if col not in residual_targets]
            )
        metadata_columns = selected_columns

    all_model_features = None

    # Iterate over the metadata columns
    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing metadata column: {metadata_col}")
            print()
            continuous_target = metadata_col in continuous_metadata_columns

            X_train, X_test, y_train, y_test, model_features = merge_and_split_data(
                data,
                metadata,
                metadata_col,
                min_samples=args.min_samples,
                balanced_test=args.balanced_test,
                train_prop=args.train_prop,
                continuous=continuous_target,
            )

            # skip the column if the merge and split function returns None
            if X_train is None:
                print(
                    f"Skipping {metadata_col} as there are not enough samples after merging and filtering..."
                )
                print()
                continue

            if continuous_target:
                print(f"Fitting continuous residual target for {metadata_col}.")
                if args.grouped:
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

                yhat_train = np.asarray(model.predict(X_train.astype(np.float64))).flatten()
                train_r2 = r2_score(y_train, yhat_train)
                print(f"Train R2 for {metadata_col}: {train_r2:.4f}")

                if args.train_prop == 1:
                    test_r2 = None
                    print("Not calculating test R2 as we are not using test data...")
                else:
                    yhat = np.asarray(model.predict(X_test.astype(np.float64))).flatten()
                    test_r2 = r2_score(y_test, yhat)
                    print(f"Test R2 for {metadata_col}: {test_r2:.4f}")

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

                plt.figure()
                plt.scatter(y_train, yhat_train, alpha=0.8)
                min_value = min(np.min(y_train), np.min(yhat_train))
                max_value = max(np.max(y_train), np.max(yhat_train))
                plt.plot([min_value, max_value], [min_value, max_value], color="black")
                plt.title(f"{metadata_col}\nTrain R2: {train_r2:.2f}")
                plt.xlabel("Observed residual")
                plt.ylabel("Predicted residual")
                plt.tight_layout()
                pdf.savefig()
                plt.close()

                if args.train_prop < 1:
                    plt.figure()
                    plt.scatter(y_test, yhat, alpha=0.8)
                    min_value = min(np.min(y_test), np.min(yhat))
                    max_value = max(np.max(y_test), np.max(yhat))
                    plt.plot(
                        [min_value, max_value],
                        [min_value, max_value],
                        color="black",
                    )
                    plt.title(f"{metadata_col}\nTest R2: {test_r2:.2f}")
                    plt.xlabel("Observed residual")
                    plt.ylabel("Predicted residual")
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()

                print()
                continue

            num_classes = len(np.unique(y_train))
            print(f"Number of classes for {metadata_col}: {num_classes}")

            # set group ids based on feature names if --grouped is supplied
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
            else:
                plt.figure()
                plt.imshow(
                    cm_train / cm_train.sum(axis=1)[:, np.newaxis],
                    cmap="viridis",
                    vmin=0,
                    vmax=1,
                )
                plt.colorbar()
                for i in range(cm_train.shape[0]):
                    for j in range(cm_train.shape[1]):
                        plt.text(j, i, f"{cm_train[i, j]}", ha="center", va="center")
                if cm_train.shape[0] > 2:
                    plt.title(f"{metadata_col}\nTrain Accuracy: {train_accuracy:.2f}")
                else:
                    plt.title(
                        f"{metadata_col}\nTrain Accuracy: {train_accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}"
                    )
                plt.xlabel("Predicted")
                plt.ylabel("True")
                plt.xticks(range(cm_train.shape[1]), metadata_categories, rotation=45)
                # these need to start from the bottom
                plt.yticks(range(cm_train.shape[0]), metadata_categories, rotation=45)
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


if __name__ == "__main__":
    main()
