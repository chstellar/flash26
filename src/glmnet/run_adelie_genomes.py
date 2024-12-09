from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import OneHotEncoder

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import adelie as ad

import numpy as np
import scipy.stats as st

import pyarrow.feather as feather
import pandas as pd

import argparse

np.random.seed(42)

def parse_args():
    parser = argparse.ArgumentParser(description="Train a model to predict antibiotic resistance")
    parser.add_argument("--train_features", type=str, help="Path to the training features file", required=True)
    parser.add_argument("--train_metadata", type=str, help="Path to the training metadata file", required=True)
    parser.add_argument("--test_features", type=str, help="Path to the test features file", required=True)
    parser.add_argument("--test_metadata", type=str, help="Path to the test metadata file", required=True)
    parser.add_argument("--output_prefix", type=str, help="Prefix for the output files", required=True)
    parser.add_argument("--min_samples", type=int, default=100, help="Minimum number of samples per category to keep")
    parser.add_argument("--n_threads", type=int, default=1, help="Number of threads to use for training the model")
    return parser.parse_args()

def read_feather_data(file_path):
    return feather.read_feather(file_path)

def read_metadata(file_path):
    metadata = pd.read_table(file_path)
    return metadata

def get_metadata_columns(metadata, min_samples=50):
    filtered_metadata = metadata.loc[:, metadata.columns != "sample_name"]
    filtered_metadata = filtered_metadata.loc[:, filtered_metadata.apply(lambda x: len(x.unique()) > 2, axis=0)]
    filtered_metadata = filtered_metadata.loc[:, filtered_metadata.apply(lambda x: sum(x.value_counts() > min_samples) > 1, axis=0)]
    return filtered_metadata.columns

def merge_data(data, metadata, metadata_col, min_samples=50, even_samples=False):
    metadata = metadata[["sample_name", metadata_col]]
    merged_data = pd.merge(data, metadata, on='sample_name', how='left')
    merged_data = merged_data.dropna(subset=[metadata_col])

    class_counts = merged_data[metadata_col].value_counts()
    class_counts = class_counts[class_counts >= min_samples]
    classes_to_keep = class_counts.index
    if len(classes_to_keep) < 2:
        return None, None, None
    merged_data = merged_data[merged_data[metadata_col].isin(classes_to_keep)]

    if even_samples:
        num_to_keep = class_counts.min()
        indices_to_keep = merged_data.groupby(metadata_col).apply(lambda x: x.sample(n=num_to_keep, replace=False).index, include_groups=False).explode()
        merged_data = merged_data.loc[indices_to_keep]

    X = merged_data.drop(["sample_name", metadata_col], axis=1)
    y = merged_data[metadata_col].to_numpy()
    return np.asfortranarray(X), y, X.columns

def train_adelie_model(X_train, y_train, n_threads=1):
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
    
    train_features = read_feather_data(args.train_features)
    train_metadata = read_metadata(args.train_metadata)
    test_features = read_feather_data(args.test_features)
    test_metadata = read_metadata(args.test_metadata)
    
    train_metadata_columns = get_metadata_columns(train_metadata, min_samples=args.min_samples)
    test_metadata_columns = get_metadata_columns(test_metadata, min_samples=0)
    metadata_columns = [col for col in train_metadata_columns if col in test_metadata_columns]
    
    all_model_features = None

    with PdfPages(output_pdf) as pdf:
        for metadata_col in metadata_columns:
            print(f"Processing metadata column: {metadata_col}")
            print()
            
            X_train, y_train, model_features = merge_data(train_features, train_metadata, metadata_col, min_samples=args.min_samples, even_samples=True)
            X_test, y_test, _ = merge_data(test_features, test_metadata, metadata_col, min_samples=0)
            if X_test is not None:
                test_classes_to_keep = np.isin(y_test, np.unique(y_train))
                X_test = X_test[test_classes_to_keep]
                y_test = y_test[test_classes_to_keep]

            if X_train is None or X_test is None:
                print(f"Skipping {metadata_col} as there are not enough samples after merging and filtering...")
                print()
                continue

            try:
                model, oh = train_adelie_model(X_train, y_train, n_threads=args.n_threads)
            except Exception as e:
                print(f"Failed to train model for {metadata_col}: {e}")
                continue

            yhat = model.predict(X_test.astype(np.float64))
            st.contingency.crosstab(y_test, yhat).count
            if len(np.unique(yhat)) < 2:
                print(f"Skipping {metadata_col} as there are not enough unique predictions...")
                continue
            yhat_2d = np.zeros((yhat.size, yhat.max() + 1))
            yhat_2d[np.arange(yhat.size), yhat] = 1
            yhat = yhat_2d
            
            y_pred = oh.inverse_transform(yhat).flatten()
            cm = confusion_matrix(y_test, y_pred)
            print(f"Confusion matrix for {metadata_col}")
            print(cm)

            coef = model.coef_
            metadata_categories = oh.categories_[0]
            model_features = [f"{feature}+{category}" for feature in model_features for category in metadata_categories]
            
            model_features = pd.DataFrame(model_features, columns=["feature"])
            model_features["coefficient"] = coef.toarray().flatten()
            model_features = model_features[model_features["coefficient"] != 0]

            model_features["feature"] = model_features["feature"].str.split("+")
            model_features["category"] = model_features["feature"].str[1]
            model_features["feature"] = model_features["feature"].str[0]
            
            model_features = model_features.groupby("feature").agg({"coefficient": list, "category": list}).reset_index()
            model_features["classes"] = model_features["category"].apply(lambda x: "[" + ",".join(map(str, x)) + "]")
            model_features["coefficients"] = model_features["coefficient"].apply(lambda x: "[" + ",".join(map(str, x)) + "]")

            model_features["metadata_category"] = metadata_col

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

            model_features["accuracy"] = accuracy
            model_features["specificity"] = specificity if specificity is not None else "NA"
            model_features["sensitivity"] = sensitivity if sensitivity is not None else "NA"
            model_features = model_features[["metadata_category", "feature", "accuracy", "sensitivity", "specificity", "classes", "coefficients"]]

            if all_model_features is None:
                all_model_features = model_features
            else:
                all_model_features = pd.concat([all_model_features, model_features], axis=0)

            plt.figure()
            plt.imshow(cm / cm.sum(axis=1)[:, np.newaxis], cmap='viridis', vmin=0, vmax=1)
            plt.colorbar()
            for i in range(cm.shape[0]):
                for j in range(cm.shape[1]):
                    plt.text(j, i, f"{cm[i, j]}", ha='center', va='center')
            if cm.shape[0] > 2:
                plt.title(f"{metadata_col}\nAccuracy: {accuracy:.2f}")
            else:
                plt.title(f"{metadata_col}\nAccuracy: {accuracy:.2f}, Specificity: {specificity:.2f}, Sensitivity: {sensitivity:.2f}")
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
            plt.text(0.5, 0.5, "No metadata columns were processed", ha='center', va='center')
            plt.axis('off')
            pdf.savefig()
            plt.close()
        
            columns = ["metadata_category", "feature", "accuracy", "sensitivity", "specificity", "classes", "coefficients"]
            pd.DataFrame(columns=columns).to_csv(output_coef, sep="\t", index=False)
        else:
            all_model_features.to_csv(output_coef, sep="\t", index=False, float_format="%.4f")

if __name__ == "__main__":
    main()
