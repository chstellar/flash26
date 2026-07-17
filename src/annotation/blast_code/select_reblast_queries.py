import argparse
import os
import re

import pandas as pd
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select queried FASTA records that did not receive usable BLAST annotation."
    )
    parser.add_argument("--query_folder", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--output_fasta", required=True)
    parser.add_argument("--mode", choices=["blast", "blastp"], required=True)
    return parser.parse_args()


def read_query_records(query_folder):
    records = []
    for filename in os.listdir(query_folder):
        if filename.endswith(".fasta"):
            records.extend(list(SeqIO.parse(os.path.join(query_folder, filename), "fasta")))
    return records


def has_text(value):
    if pd.isna(value):
        return False
    value = str(value).strip()
    return value not in {"", "NA", "NaN", "None", "none", "[]"}


def feature_has_label(value):
    if not has_text(value):
        return False
    value = str(value)
    if re.search(r"['\"](?:gene|product)['\"]\s*:\s*\[[^\]]*[A-Za-z][^\]]*\]", value):
        return True
    return False


def annotated_queries(annotation_file, mode):
    if not os.path.exists(annotation_file) or os.path.getsize(annotation_file) == 0:
        return set()

    try:
        annotations = pd.read_csv(annotation_file, sep="\t")
    except pd.errors.EmptyDataError:
        return set()

    if "query" not in annotations.columns or annotations.empty:
        return set()

    if mode == "blastp":
        label_cols = [
            col
            for col in ["stitle", "NCBI_protein_accession", "UniProt_accession", "GO"]
            if col in annotations.columns
        ]
        if not label_cols:
            return set()
        label_mask = annotations[label_cols].apply(
            lambda row: any(has_text(value) for value in row), axis=1
        )
    else:
        feature_cols = [
            col for col in annotations.columns if col == "features" or col.startswith("features_")
        ]
        if not feature_cols:
            return set()
        label_mask = annotations[feature_cols].apply(
            lambda row: any(feature_has_label(value) for value in row), axis=1
        )

    return set(annotations.loc[label_mask, "query"].dropna().astype(str))


def main():
    args = parse_args()
    records = read_query_records(args.query_folder)
    annotated = annotated_queries(args.annotations, args.mode)
    selected = [record for record in records if record.id not in annotated]

    os.makedirs(os.path.dirname(os.path.abspath(args.output_fasta)), exist_ok=True)
    SeqIO.write(selected, args.output_fasta, "fasta")
    print(
        f"Selected {len(selected)}/{len(records)} {args.mode} query sequences for unrestricted reblast."
    )


if __name__ == "__main__":
    main()
