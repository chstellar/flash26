import os
import pandas as pd
from Bio import SeqIO
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Lookup Table Tool on a FASTA file."
    )
    parser.add_argument(
        "--cluster_seqs",
        type=str,
        help="Input cluster file.",
    )
    parser.add_argument(
        "--lookup_table",
        type=str,
        help="Path to lookup table file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file.",
    )
    parser.add_argument(
        "--temp_dir",
        type=str,
        help="This folder will be created and moved-to.",
    )
    parser.add_argument(
        "--splash_bin",
        default="/oak/stanford/groups/horence/dcotter1/splash-2.6.1/",
        type=str,
        help="Path to splash bin",
    )
    return parser.parse_args()


def create_temp_fasta(cluster_file, fasta_file):
    df = pd.read_csv(cluster_file, sep='\t')
    seqs = df['seq'].unique()
    
    with open(fasta_file, 'w') as f:
        for seq in seqs:
            f.write(f">{seq}\n{seq}\n")

def run_lookup(anchor_fasta, lookup_file, output_file, splash_bin):
    # run the lookup command
    lookup_table = Path(splash_bin) / "lookup_table"
    lookup_cmd = (
        f"{lookup_table} query "
        f"--kmer_skip 1 --truncate_paths --stats_fmt with_stats "
        f"{Path(anchor_fasta).resolve()} {Path(lookup_file).resolve()} {Path(output_file).resolve()}"
    )
    print(f"Running command: {lookup_cmd}")
    os.system(lookup_cmd)
    return None

def merge_results(cluster_file, lookup_table_file, output_file):
    df_cluster = pd.read_csv(cluster_file, sep='\t')
    df_lookup = pd.read_csv(lookup_table_file, sep='\t')
    # there is no common key but the order is the same
    df_merged = pd.merge(df_cluster, df_lookup, left_index=True, right_index=True)
    df_merged.to_csv(output_file, sep='\t', index=False)

def main():
    # parse arguments
    args = parse_args()
    cluster_file = args.cluster_seqs
    lookup_table = args.lookup_table
    final_output_file = args.output
    temp_fasta_file = os.path.join(args.temp_dir, 'temp.fasta')
    lookup_table_output = os.path.join(args.temp_dir, 'lookup_table_output.txt')
    splash_bin = args.splash_bin
    
    # create the temp dir if it doesn't exist
    if not os.path.exists(args.temp_dir):
        os.mkdir(args.temp_dir)
    
    # create the temp fasta file
    create_temp_fasta(cluster_file, temp_fasta_file)

    # run lookup table on the temp fasta file
    run_lookup(temp_fasta_file, lookup_table, lookup_table_output, splash_bin)

    # merge the original cluster file with the output of the lookup table
    merge_results(cluster_file, lookup_table_output, final_output_file)
    
    os.remove(temp_fasta_file)
    os.remove(lookup_table_output)

if __name__ == "__main__":
    main()
