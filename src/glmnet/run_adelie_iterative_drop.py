from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import OneHotEncoder

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np
import scipy.stats as st

import pyarrow.feather as feather
import pandas as pd
from math import floor

from os.path import basename
import argparse

import re

np.random.seed(42)

def parse_args():
    parser = argparse.ArgumentParser(description="Train a model to predict antibiotic resistance")
    parser.add_argument("--data", type=str, help="Path to the data file", required=True)
    parser.add_argument("--metadata", type=str, help="Path to the metadata file", required=True)
    parser.add_argument("--metadata_col", type=str, help="Metadata column to use for training the model", required=True)
    parser.add_argument("--output_prefix", type=str, help="Prefix for the output files", required=True)
    parser.add_argument("--min_samples", type=int, default=30, help="Minimum number of samples per category to keep")
    parser.add_argument("--n_threads", type=int, default=1, help="Number of threads to use for training the model")
    parser.add_argument("--balanced_test", action="store_true", help="Keep the same number of samples per class in the test set")
    parser.add_argument("--n_iter", type=int, default=100, help="Number of iterations to drop features")
    parser.add_argument("--drop_full_cluster", action="store_true", help="Drop all features that belong to the same cluster")
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
    filtered_metadata = filtered_metadata.loc[:, filtered_metadata.apply(lambda x: len(x.unique()) >= 2, axis=0)]
    # only grab columns with two or more categories that have more than min_samples
    filtered_metadata = filtered_metadata.loc[:, filtered_metadata.apply(lambda x: sum(x.value_counts() > min_samples) > 1, axis=0)]
    return filtered_metadata.columns


def merge_and_split_data(data, metadata, metadata_col, min_samples=50, train_prop=0.5, balanced_test=False):
    metadata = metadata[["sample_name", metadata_col]]
    merged_data = pd.merge(data, metadata, on='sample_name', how='left')
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
    indices_to_keep = merged_data.groupby(metadata_col).apply(lambda x: x.sample(n=num_to_keep, replace=False).index, include_groups=False).explode()
    
    # If we want a balanced test set, keep the same number of samples per class in the test set
    if balanced_test:
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[indices_to_keep]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()

        test_indices = merged_data.drop(indices_to_keep).groupby(metadata_col).apply(lambda x: x.sample(n=num_to_keep, replace=False).index, include_groups=False).explode()
        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).loc[test_indices]
        y_test = merged_data[metadata_col].loc[test_indices].to_numpy()
    else:
        # Split the data into training and test sets
        X_train = merged_data.drop(["sample_name", metadata_col], axis=1).loc[indices_to_keep]
        model_features = X_train.columns
        y_train = merged_data[metadata_col].loc[indices_to_keep].to_numpy()
        
        X_test = merged_data.drop(["sample_name", metadata_col], axis=1).drop(indices_to_keep)
        y_test = merged_data[metadata_col].drop(indices_to_keep).to_numpy()
    
    return np.asfortranarray(X_train), np.asfortranarray(X_test), y_train, y_test, model_features

def train_adelie_model(X_train, y_train,n_threads=1):
    oh = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    y_train2 = oh.fit_transform(y_train[:, np.newaxis])
    
    model = ad.GroupElasticNet(solver="cv_grpnet", family="multinomial")
    model.fit(X_train.astype(np.float64), y_train2.astype(np.float64), n_threads=n_threads)
    
    return model, oh

def main():
    args = parse_args()
    output_prefix = args.output_prefix
    output_pdf = output_prefix + "_confusion_matrices.pdf"
    output_coef = output_prefix + "_nonzero_coefficients.tsv"
    
    # Load teh data and metadata
    data = read_feather_data(args.data)
    metadata = read_metadata(args.metadata)
    # Get the metadata columns that have more than 2 unique values
    # and more than 50 samples per category
    metadata_columns = get_metadata_columns(metadata, min_samples=args.min_samples)

    N_ITER = args.n_iter
    drop_full_cluster = args.drop_full_cluster

    metadata_col = args.metadata_col
    if metadata_col not in metadata_columns:
        raise ValueError(f"Metadata column {metadata_col} not found in metadata file")
    
    all_model_features = None

    print(f"Processing metadata column: {metadata_col}")
    print()
    # initial loop use the full X matrix
    X_train, X_test, y_train, y_test, model_features = merge_and_split_data(data, metadata, metadata_col, min_samples=args.min_samples, balanced_test=args.balanced_test)
    
    original_X_train = X_train.copy()
    original_X_test = X_test.copy()
    original_model_features = model_features.copy()

    if X_train is None:
        print(f"Skipping {metadata_col} as there are not enough samples after merging and filtering...")
        print()
        return

    try:
        model, oh = train_adelie_model(X_train, y_train, n_threads=args.n_threads)
    except Exception as e:
        print(f"Failed to train model for {metadata_col}: {e}")
        return

    # add a check to make sure there are more than 2 unique values in the predictions
    # an error can be thrown if inverse_transform gets the wrong number of columns
    yhat = model.predict(X_test.astype(np.float64))
    if len(np.unique(yhat)) < 2:
        print(f"Skipping {metadata_col} as there are not enough unique predictions...")
        return
    yhat_2d = np.zeros((yhat.size, yhat.max() + 1))
    yhat_2d[np.arange(yhat.size), yhat] = 1
    yhat = yhat_2d
    
    try:
        y_pred = oh.inverse_transform(yhat).flatten()
        cm = confusion_matrix(y_test, y_pred)
        print(f"Confusion matrix for {metadata_col}")
        print(cm)
    except Exception as e:
        print(f"Failed to transform predictions for {metadata_col}: {e}")
        return

    # extract the nonzero coefficients
    coef = model.coef_
    # get the feature names for the nonzero coefficients
    metadata_categories = oh.categories_[0]
    model_features = [f"{feature}+{category}" for feature in model_features for category in metadata_categories]
    
    # get the names and values of the nonzero coefficients
    model_features = pd.DataFrame(model_features, columns=["feature"])
    model_features["coefficient"] = coef.toarray().flatten()
    model_features = model_features[model_features["coefficient"] != 0]

    # separate the feature names into the feature and category
    model_features["feature"] = model_features["feature"].str.split("+")
    model_features["category"] = model_features["feature"].str[1]
    model_features["feature"] = model_features["feature"].str[0]
    
    # gather the coefficients and categorie names into two columns grouping by feature
    model_features = model_features.groupby("feature").agg({"coefficient": list, "category": list}).reset_index()
    model_features["classes"] = model_features["category"].apply(lambda x: "[" + ",".join(map(str, x)) + "]")
    model_features["coefficients"] = model_features["coefficient"].apply(lambda x: "[" + ",".join(map(str, x)) + "]")

    # add column for metadata category
    model_features["metadata_category"] = metadata_col

    # assuming cm can be larger than 2x2
    if cm.shape[0] > 2:
        accuracy = np.trace(cm) / np.sum(cm)
        specificity = None
        sensitivity = None
    else:
        accuracy = (cm[0][0] + cm[1][1]) / sum(sum(row) for row in cm)
        specificity = cm[0][0] / (cm[0][0] + cm[0][1])
        sensitivity = cm[1][1] / (cm[1][0] + cm[1][1])
    
    if specificity is None:
        print(f"Accuracy: {accuracy:.2f}")
    else:
        print(f"Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}")

    # add the accuracy, specificity, and sensitivity to the model features
    model_features["accuracy"] = accuracy
    model_features["specificity"] = specificity if specificity is not None else "NA"
    model_features["sensitivity"] = sensitivity if sensitivity is not None else "NA"
    model_features = model_features[["metadata_category", "feature", "accuracy", "sensitivity", "specificity", "classes", "coefficients"]]
    model_features["iteration"] = 0

    # join with the larger set of model features
    all_model_features = model_features
    # empty numpy ndarray to store the features to drop
    features_to_drop = []
    features_to_drop = np.array(features_to_drop)
    

    # now remove the features with nonzero coefficients from the X matrix and repeat N_ITER times
    for i in range(1,N_ITER):
        if drop_full_cluster:
            # drop all features that belong to the same cluster
            clusters = np.unique([re.sub(r"_embedding_\d+$", "", feature) for feature in model_features["feature"].values])
            features_to_drop = np.concatenate((features_to_drop, clusters))
        else:
            features_to_drop = np.concatenate((features_to_drop, model_features["feature"].values))
        X_train = original_X_train.copy()
        X_test = original_X_test.copy()
        # use the model features list to get the column indices to drop and extend the cols_to_drop list
        if drop_full_cluster:
            cols_to_drop = [i for i, feature in enumerate(original_model_features) if any([re.sub(r"_embedding_\d+$", "", feature) in features_to_drop])]
        else:
            cols_to_drop = [i for i, feature in enumerate(original_model_features) if feature in features_to_drop]
        X_train = np.delete(X_train, cols_to_drop, axis=1)
        X_test = np.delete(X_test, cols_to_drop, axis=1)

        # drop the col positions from the original model features list
        model_features = np.delete(original_model_features, cols_to_drop)
        
        if X_train.shape[1] == 0:
            print(f"No more features to drop for {metadata_col} iter {i}")
            break
        print(f"Processing metadata column: {metadata_col} iter {i}")

        try:
            model, oh = train_adelie_model(X_train, y_train, n_threads=args.n_threads)
        except Exception as e:
            print(f"Failed to train model for {metadata_col} iter {i}: {e}")
            return
        
        yhat = model.predict(X_test.astype(np.float64))
        if len(np.unique(yhat)) < 2:
            print(f"Skipping {metadata_col} as there are not enough unique predictions...")
            break
        yhat_2d = np.zeros((yhat.size, yhat.max() + 1))
        yhat_2d[np.arange(yhat.size), yhat] = 1
        yhat = yhat_2d

        try:
            y_pred = oh.inverse_transform(yhat).flatten()
            cm = confusion_matrix(y_test, y_pred)
            print(f"Confusion matrix for {metadata_col} iter {i}")
            print(cm)
        except Exception as e:
            print(f"Failed to transform predictions for {metadata_col} iter {i}: {e}")
            return
        
        # extract the nonzero coefficients
        coef = model.coef_
        # get the feature names for the nonzero coefficients
        metadata_categories = oh.categories_[0]
        model_features = [f"{feature}+{category}" for feature in model_features for category in metadata_categories]

        # get the names and values of the nonzero coefficients
        model_features = pd.DataFrame(model_features, columns=["feature"])
        model_features["coefficient"] = coef.toarray().flatten()
        model_features = model_features[model_features["coefficient"] != 0]

        # separate the feature names into the feature and category
        model_features["feature"] = model_features["feature"].str.split("+")
        model_features["category"] = model_features["feature"].str[1]
        model_features["feature"] = model_features["feature"].str[0]

        # gather the coefficients and categorie names into two columns grouping by feature
        model_features = model_features.groupby("feature").agg({"coefficient": list, "category": list}).reset_index()
        model_features["classes"] = model_features["category"].apply(lambda x: "[" + ",".join(map(str, x)) + "]")
        model_features["coefficients"] = model_features["coefficient"].apply(lambda x: "[" + ",".join(map(str, x)) + "]")

        # add column for metadata category
        model_features["metadata_category"] = metadata_col

        # assuming cm can be larger than 2x2
        if cm.shape[0] > 2:
            accuracy = np.trace(cm) / np.sum(cm)
            specificity = None
            sensitivity = None
        else:
            accuracy = (cm[0][0] + cm[1][1]) / sum(sum(row) for row in cm)
            specificity = cm[0][0] / (cm[0][0] + cm[0][1])
            sensitivity = cm[1][1] / (cm[1][0] + cm[1][1])

        if specificity is None:
            print(f"Accuracy: {accuracy:.2f}")
        else:
            print(f"Accuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}")
        
        # add the accuracy, specificity, and sensitivity to the model features
        model_features["accuracy"] = accuracy
        model_features["specificity"] = specificity if specificity is not None else "NA"
        model_features["sensitivity"] = sensitivity if sensitivity is not None else "NA"
        model_features = model_features[["metadata_category", "feature", "accuracy", "sensitivity", "specificity", "classes", "coefficients"]]
        model_features["iteration"] = i

        # join with the larger set of model features
        all_model_features = pd.concat([all_model_features, model_features], axis=0)

    # add a blank line 
    print()

    if all_model_features is None:
        # Write an empty output TSV with just column names
        columns = ["metadata_category", "feature", "accuracy", "sensitivity", "specificity", "classes", "coefficients"]
        pd.DataFrame(columns=columns).to_csv(output_coef, sep="\t", index=False)
    else:
        # output the nonzero coefficients to a tsv file
        all_model_features.to_csv(output_coef, sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
