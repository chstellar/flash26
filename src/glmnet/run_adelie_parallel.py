# from sklearn.model_selection import train_test_split
# from sklearn.metrics import r2_score
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

# for parallelization
import concurrent.futures
import multiprocessing
import sys
import os

np.random.seed(42)

# Add globals to hold main process-loaded data
PROCESS_DATA = None
PROCESS_METADATA = None


# def worker_init(data_path, metadata_path):
#     """
#     Initializer run once per worker process to load data and metadata into
#     module-level globals so worker tasks don't reload files repeatedly.
#     """
#     global PROCESS_DATA, PROCESS_METADATA
#     try:
#         PROCESS_DATA = read_feather_data(data_path)
#         PROCESS_METADATA = read_metadata(metadata_path)
#         print(
#             f"Worker init complete (pid={os.getpid()}), data shape: {PROCESS_DATA.shape}, metadata shape: {PROCESS_METADATA.shape}",
#             file=sys.stderr,
#         )
#     except Exception:
#         tb = traceback.format_exc()
#         print(f"Worker init failed (pid={os.getpid()})", file=sys.stderr)
#         print(tb, file=sys.stderr)
#         raise
#     return None


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
    return parser.parse_args()


def read_feather_data(file_path):
    return feather.read_feather(file_path)


def read_metadata(file_path):
    metadata = pd.read_table(file_path)
    if "sample_name" not in metadata.columns:
        raise ValueError("Metadata file must contain a sample_name column")
    # mutate all columns to strings for categorical analysis
    metadata = metadata.apply(lambda x: x.astype(str))
    return metadata


def get_metadata_columns(metadata, min_samples=50):
    """
    Returns the columns of the metadata file except for the sample_name column
    Filters the columns so it will only return metadata that have more than two
    discrete values with greater than min_samples per category
    """
    filtered_metadata = metadata.loc[:, metadata.columns != "sample_name"]
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
    data, metadata, metadata_col, min_samples=50, train_prop=0.5, balanced_test=False
):
    metadata = metadata[["sample_name", metadata_col]]
    merged_data = pd.merge(data, metadata, on="sample_name", how="left")
    merged_data = merged_data.dropna(subset=[metadata_col])

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


def process_column(args_tuple):
    """
    Worker executed in separate process: loads data+metadata from paths and
    runs merge_and_split_data + train_adelie_model for one metadata column.
    Returns a serializable dict with results (dataframe as dict records + cms + metrics).
    """
    (
        data_path,
        metadata_path,
        metadata_col,
        min_samples,
        train_prop,
        balanced_test,
        grouped,
        tol,
        max_iters,
        alpha,
        per_worker_threads,
    ) = args_tuple

    # use global data loaded in main process
    global PROCESS_DATA, PROCESS_METADATA
    data = PROCESS_DATA
    metadata = PROCESS_METADATA

    if data is None or metadata is None:
        raise ValueError("Worker data or metadata not loaded")

    # call existing function
    X_train, X_test, y_train, y_test, model_features = merge_and_split_data(
        data,
        metadata,
        metadata_col,
        min_samples=min_samples,
        train_prop=train_prop,
        balanced_test=balanced_test,
    )
    if X_train is None:
        return {"skipped": True, "metadata_col": metadata_col}

    num_classes = len(np.unique(y_train))
    group_ids = None
    if grouped and num_classes < 4:
        group_ids = get_group_ids(model_features)

    try:
        print(
            f"Training model for metadata column {metadata_col} (pid={os.getpid()})",
            file=sys.stderr,
        )
        model, oh = train_adelie_model(
            X_train,
            y_train,
            n_threads=per_worker_threads,
            group_ids=group_ids,
            tol=tol,
            max_iters=max_iters,
            alpha=alpha,
        )
    except Exception as e:
        return {"skipped": True, "metadata_col": metadata_col, "error": str(e)}

    # predictions (same logic as main, but keep minimal serializable results)
    result = {"skipped": False, "metadata_col": metadata_col}
    try:
        if train_prop < 1:
            yhat = model.predict(X_test.astype(np.float64))
            if len(np.unique(yhat)) < 2:
                unique_class = np.unique(yhat)[0]
                yhat_2d = np.zeros((y_test.size, len(oh.categories_[0])))
                yhat_2d[:, unique_class] = 1
                y_pred = oh.inverse_transform(yhat_2d).flatten()
            else:
                yhat_2d = np.zeros((yhat.size, len(oh.categories_[0])))
                yhat_2d[np.arange(yhat.size), yhat] = 1
                y_pred = oh.inverse_transform(yhat_2d).flatten()
            cm = confusion_matrix(y_test, y_pred, labels=oh.categories_[0])
            result["cm"] = cm.tolist()
        else:
            result["cm"] = None
    except Exception as e:
        result["cm_error"] = str(e)
        result["cm"] = None

    # train confusion
    yhat_train = model.predict(X_train.astype(np.float64))
    if len(np.unique(yhat_train)) < 2:
        unique_class_train = np.unique(yhat_train)[0]
        yhat_train_2d = np.zeros((y_train.size, len(oh.categories_[0])))
        yhat_train_2d[:, unique_class_train] = 1
        y_train_pred = oh.inverse_transform(yhat_train_2d).flatten()
    else:
        yhat_train_2d = np.zeros((yhat_train.size, len(oh.categories_[0])))
        yhat_train_2d[np.arange(yhat_train.size), yhat_train] = 1
        y_train_pred = oh.inverse_transform(yhat_train_2d).flatten()
    cm_train = confusion_matrix(y_train, y_train_pred, labels=oh.categories_[0])
    result["cm_train"] = cm_train.tolist()

    # coefficients (serialize)
    coef = model.coef_.toarray().flatten().tolist()
    metadata_categories = list(oh.categories_[0])
    # expand feature names like in main
    expanded_feat = [
        f"{feature}+{cat}" for feature in model_features for cat in metadata_categories
    ]
    result["coef_features"] = expanded_feat
    result["coef_values"] = coef
    result["metadata_categories"] = metadata_categories
    result["model_features"] = model_features.tolist()

    return result


def main():
    args = parse_args()
    output_prefix = args.output_prefix
    output_pdf = output_prefix + "_confusion_matrices.pdf"
    output_coef = output_prefix + "_nonzero_coefficients.tsv"

    # Load the metadata to determine which columns to process
    metadata = read_metadata(args.metadata)
    metadata_columns = get_metadata_columns(metadata, min_samples=args.min_samples)

    # create the empty all_model_features dataframe
    all_model_features = None

    # Prepare worker args for running in parallel
    worker_args = []
    # for memory management we only want to use 4-8 workers at a time and use multiple threads per worker
    total_workers = min(
        4, len(metadata_columns) // 6
    )  # at max 4 workers, at least 6 columns per worker
    per_worker_threads = (args.n_threads - 2) // total_workers
    if per_worker_threads < 1:
        per_worker_threads = 1
    print(
        f"Preparing to process {len(metadata_columns)} metadata columns...",
        file=sys.stderr,
    )
    print(
        f"Using {total_workers} parallel workers with {per_worker_threads} threads each",
        file=sys.stderr,
    )
    for col in metadata_columns:
        worker_args.append(
            (
                args.data,
                args.metadata,
                col,
                args.min_samples,
                args.train_prop,
                args.balanced_test,
                args.grouped,
                args.tol,
                args.max_iters,
                args.alpha,
                per_worker_threads,  # tune as needed using args.n_threads
            )
        )

    # Load data and metadata once in main process for workers to inherit
    global PROCESS_DATA, PROCESS_METADATA
    PROCESS_DATA = read_feather_data(args.data)
    PROCESS_METADATA = read_metadata(args.metadata)
    print(
        f"Main process loaded data shape: {PROCESS_DATA.shape}, metadata shape: {PROCESS_METADATA.shape}",
        file=sys.stderr,
    )
    # Run processing in parallel
    # choose executor: prefer forked processes to inherit memory, otherwise use threads
    results = []  # hold completed results
    if "fork" in multiprocessing.get_all_start_methods():
        print("Using forked processes for parallelization", file=sys.stderr)
        ctx = multiprocessing.get_context("fork")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=total_workers, mp_context=ctx
        ) as ex:
            futs = [ex.submit(process_column, a) for a in worker_args]
            for fut in concurrent.futures.as_completed(futs):
                try:
                    results.append(fut.result())
                except Exception as e:
                    print(f"Error processing column: {e}", file=sys.stderr)
    else:
        # fall back to threads to avoid reloading huge file per process
        print("Using threads for parallelization", file=sys.stderr)
        with concurrent.futures.ThreadPoolExecutor(max_workers=total_workers) as ex:
            futs = [ex.submit(process_column, a) for a in worker_args]
            for fut in concurrent.futures.as_completed(futs):
                try:
                    results.append(fut.result())
                except Exception as e:
                    print(f"Error processing column: {e}", file=sys.stderr)

    # use `results` below instead of `futures`
    with PdfPages(output_pdf) as pdf:
        for res in results:
            if res.get("skipped"):
                continue
            # reconstruct DataFrame of nonzero coeffs similar to main
            feat = res["coef_features"]
            vals = res["coef_values"]
            df = pd.DataFrame({"feature": feat, "coefficient": vals})
            df = df[df["coefficient"] != 0]
            # split feature+category
            df["feature"] = df["feature"].str.split("+")
            df["category"] = df["feature"].str[1]
            df["feature"] = df["feature"].str[0]
            df = (
                df.groupby("feature")
                .agg({"coefficient": list, "category": list})
                .reset_index()
            )
            df["classes"] = df["category"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )
            df["coefficients"] = df["coefficient"].apply(
                lambda x: "[" + ",".join(map(str, x)) + "]"
            )
            df["metadata_category"] = res["metadata_col"]
            # append to all_model_features
            if all_model_features is None:
                all_model_features = df
            else:
                all_model_features = pd.concat([all_model_features, df], axis=0)

            # plot the confusion matrix and save to the pdf
            # color by relative frequency
            # add numbers to the cells of the confusion matrix
            # add accuracy, specificity, and sensitivity to the title
            if args.train_prop < 1:
                cm = np.array(res["cm"])
                metadata_col = res["metadata_col"]
                accuracy = res["accuracy"]
                specificity = res["specificity"]
                sensitivity = res["sensitivity"]
                metadata_categories = res["metadata_categories"]
                cm_train = np.array(res["cm_train"])
                train_accuracy = res["train_accuracy"]
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

        # write TSV as before
        if all_model_features is None:
            # Write a blank page to the PDF
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
            all_model_features.to_csv(
                output_coef, sep="\t", index=False, float_format="%.4f"
            )


# Run the main function
if __name__ == "__main__":
    main()
