from Bio import Entrez, SeqIO
import pandas as pd
import time
import argparse
import sys
import os
from os.path import join, basename
from concurrent.futures import ThreadPoolExecutor
from math import floor

Entrez.email = "dcotter1@stanford.edu"
MAX_RETRIES = 100

def fetch_sequence(seq_id):
    print(f"Fetching sequence {seq_id}")
    for i in range(MAX_RETRIES):
        try:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            break
        except:
            time.sleep(5)
    return record

def find_overlapping_features(record, window_start, window_end):
    overlapping_features = []
    for feature in record.features:
        if feature.type == "source":
            continue
        feature_start = feature.location.start
        feature_end = feature.location.end

        overlap = (feature_start <= window_end) and (feature_end >= window_start)

        if overlap:
            overlapping_features.append({
                "type": feature.type,
                "start": str(feature_start),
                "end": str(feature_end),
                "gene": feature.qualifiers.get("gene"),
                "product": feature.qualifiers.get("product"),
                "protein_seq": feature.qualifiers.get("translation")
            })

    return overlapping_features

def featurize_blast_out(blast_out, window=5000):
    df = pd.read_csv(blast_out, sep="\t", header=None)
    df.columns = ["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",\
                     "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "sgi", \
                        "sacc", "slen", "staxids", "stitle"]
    df["features"] = None
    df[f"features_{window}_window"] = None
    sacc_records = {}
    for sacc in df["sacc"].unique():
        sacc_records[sacc] = fetch_sequence(sacc)

    for index, row in df.iterrows():
        record = sacc_records[row["sacc"]]
        features = find_overlapping_features(record, row["s_start"], row["s_end"])
        df.at[index, "features"] = features
        window_start = max(row["s_start"] - window, 0)
        window_end = row["s_end"] + window
        features = find_overlapping_features(record, window_start, window_end)
        df.at[index, f"features_{window}_window"] = features
    
    return df[["query", "identity", "features", f"features_{window}_window"]]

def process_blast_file(blast_out, blast_feat_out, blast_window):
    if os.path.exists(blast_feat_out) and os.path.getsize(blast_feat_out) > 0:
        print(f"Output file {blast_feat_out} exists and has data. Skipping.")
        return
    df_features = featurize_blast_out(blast_out, blast_window)
    df_features.to_csv(blast_feat_out, index=None, sep="\t")
    print(f"Featurize blast output complete for {blast_out}. Output file: {blast_feat_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Featurize BLAST output")
    parser.add_argument("--blast_folder", required=True, help="Path to the folder containing BLAST output files")
    parser.add_argument("--blast_window", type=int, default=10000, help="Window size for feature extraction")
    parser.add_argument("--output_file", required=True, help="Path to the output file for concatenated results")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum number of workers for parallel processing")
    args = parser.parse_args()

    blast_folder = args.blast_folder
    blast_window = args.blast_window
    output_file = args.output_file
    max_workers = args.max_workers

    if not os.path.exists(blast_folder):
        print(f"Blast folder {blast_folder} does not exist. Exiting.")
        sys.exit(0)

    blast_outs = [join(blast_folder, f) for f in os.listdir(blast_folder) if f.endswith(".blastout.tsv")]
    blast_feat_outs = [join(blast_folder, basename(f).split(".")[0] + ".blastfeatout.tsv") for f in blast_outs]
    print(f"Total blast output files: {len(blast_outs)}")

    with ThreadPoolExecutor(max_workers=max(2,floor(max_workers/2))) as executor:
        futures = [executor.submit(process_blast_file, blast_out, blast_feat_out, blast_window) for blast_out, blast_feat_out in zip(blast_outs, blast_feat_outs)]
        for future in futures:
            future.result()

    # Concatenate all .blastfeatout.tsv files into the output file
    all_feat_outs = [pd.read_csv(join(blast_folder, f), sep="\t") for f in os.listdir(blast_folder) if f.endswith(".blastfeatout.tsv")]
    concatenated_df = pd.concat(all_feat_outs, ignore_index=True)
    concatenated_df.to_csv(output_file, index=None, sep="\t")
    print(f"All .blastfeatout.tsv files have been concatenated into {output_file}")
