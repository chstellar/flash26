from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.metrics import confusion_matrix
from sklearn.datasets import (
    load_breast_cancer,
    load_diabetes,
    load_digits,
)
from sklearn.preprocessing import OneHotEncoder
import adelie as ad
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st

import pyarrow.feather as feather

import pandas as pd

import argparse

np.random.seed(42)

def parse_args():
    parser = argparse.ArgumentParser(description="Train a model to predict antibiotic resistance")
    parser.add_argument("--data", type=str, help="Path to the data file")
    parser.add_argument("--metadata", type=str, help="Path to the metadata file")
    parser.add_argument("--output_prefix", type=str, help="Prefix for the output files")
    return parser.parse_args()

def read_feather_data(file_path):
    return feather.read_feather(file_path)

def read_metadata(file_path):
    metadata = pd.read_table(file_path)
    return metadata

def get_metadata_columns(metadata, min_samples=50):
    """
    Returns the columns of the metadata file except for the sample_name column
    Filters the columns so it will only return metadata that have more than two
    discrete values with greater than min_samples per category
    """
    metadata_columns = metadata.columns[metadata.columns != "sample_name"]
    # filter out columns with less than 2 unique values
    metadata_columns = metadata_columns[metadata.apply(lambda x: len(x.unique()) > 2, axis=0)]
    # only grab columns with two or more columns that have more than min_samples
    metadata_columns = metadata_columns[metadata_columns.apply(lambda x: sum(x.value_counts() > min_samples) > 1)]
    return metadata_columns


def merge_and_split_data(data, metadata, metadata_col):
    merged_data = pd.merge(data, metadata, on='sample_name', how='left')
    merged_data = merged_data.dropna(subset=[metadata_col])
    merged_data = merged_data[merged_data[metadata_col] != "I"]
    
    num_to_keep = merged_data[metadata_col].value_counts()['S']
    merged_data = merged_data.groupby(metadata_col).apply(lambda x: x.sample(n=num_to_keep, replace=False))
    
    X = merged_data.drop(["sample_name", metadata_col], axis=1)
    y = merged_data[metadata_col]
    
    X_train, X_test, y_train, y_test = train_test_split(X.to_numpy(), y.to_numpy(), test_size=0.5, stratify=y)
    
    return np.asfortranarray(X_train), np.asfortranarray(X_test), y_train, y_test

def train_adelie_model(X_train, y_train):
    oh = OneHotEncoder(sparse_output=False)
    y_train2 = oh.fit_transform(y_train[:, np.newaxis])
    
    model = ad.GroupElasticNet(solver="cv_grpnet", family="multinomial")
    model.fit(X_train.astype(np.float64), y_train2.astype(np.float64))
    
    return model

def main():
    data = read_feather_data("/scratch/users/dcotter1/metaSPLASH_workflows/eFaecium-CollEtAl/eFaecium-CollEtAl_ohe_features_for_glmnet_filter3_shiftDist-keepTopES_top50000_k54_s54.feather")
    metadata = read_metadata("/oak/stanford/groups/horence/dcotter1/utility_files/metadata/metaSPLASH_metadata/E_faecium_cleaned_resistance_metadata.tsv")
    
    X_train, X_test, y_train, y_test = merge_and_split_data(data, metadata)
    
    model = train_adelie_model(X_train, y_train)
    
    yhat = model.predict(X_test.astype(np.float64))
    
    y_pred = [['R', 'S'][x] for x in yhat]
    cm = confusion_matrix(y_test, y_pred)
    
    print(cm)
    
    accuracy = (cm[0][0] + cm[1][1]) / sum(sum(row) for row in cm)
    print(f"Accuracy: {accuracy}")

if __name__ == "__main__":
    main()