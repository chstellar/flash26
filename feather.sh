#!/bin/bash

ml python/3.12.1 py-pyarrow/18.1.0_py312 py-pandas/2.2.1_py312

python3 - <<'PY'
# import pandas as pd
# import pyarrow.feather as feather

# x = feather.read_feather("results/260621-00-temnothorax-challenge/260621-00-temnothorax-challenge_hyena_top_variance_features_for_glmnet_filter1_shiftDist-levFilter_top20000_target1_k40_s40_normalized.feather")
# m = pd.read_table("/scratch/users/jiamuyu/proj_botryllus/splash2/260620_00_temnothorax_challenge/metadata.tsv")

# print("feature matrix shape:", x.shape)
# print("metadata shape:", m.shape)
# print("feature sample examples:", x["sample_name"].head().tolist())
# print("metadata sample examples:", m["sample_name"].head().tolist())

# overlap = set(x["sample_name"]) & set(m["sample_name"])
# print("overlap:", len(overlap))
# print("feature-only samples:", len(set(x["sample_name"]) - set(m["sample_name"])))
# print("metadata-only samples:", len(set(m["sample_name"]) - set(x["sample_name"])))

# num = x.drop(columns=["sample_name"])
# print("NA cells:", num.isna().sum().sum())
# print("constant columns:", (num.nunique(dropna=False) <= 1).sum())
# print("all-zero columns:", (num.eq(0).all()).sum())

# merged = x[["sample_name"]].merge(m, on="sample_name", how="left")
# for col in m.columns:
#     if col != "sample_name":
#         print("\n", col)
#         print(merged[col].value_counts(dropna=False))

import pandas as pd

ordering = pd.read_csv("results/260621-00-temnothorax-challenge/260621-00-temnothorax-challenge_decomposed_kmers_filter1_shiftDist-levFilter_top20000_target1_k40_s40_kmer_ordering.tsv",
                       sep="\t", header=None,
                       names=["sample_name", "seq", "kmer", "start", "end"])
emb = pd.read_csv("results/260621-00-temnothorax-challenge/260621-00-temnothorax-challenge_hyena-embeddings_filter1_shiftDist-levFilter_top20000_target1_k40_s40.tsv",
                  sep="\t", header=None, usecols=[0], names=["kmer"])

print("ordering kmers:", ordering["kmer"].nunique())
print("embedding kmers:", emb["kmer"].nunique())
print("overlap:", len(set(ordering["kmer"]) & set(emb["kmer"])))
print("ordering examples:", ordering["kmer"].head().tolist())
print("embedding examples:", emb["kmer"].head().tolist())
PY