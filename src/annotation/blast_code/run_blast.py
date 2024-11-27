import os
from os.path import join, basename
import subprocess
from SeqUtils.seq_utils import read_fasta, split_fasta
import shutil
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

SPLIT_THRESH = 100
SPLIT_EACH = 50

def run_blast(splitted_fasta, blast_folder, max_workers):

    fmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseqid sgi sacc slen staxids stitle"

    def run_single_blast(f):
        blast_out = join(blast_folder, basename(f).split(".")[0] + ".blastout.tsv")
        # skip if tsv file already exists and is not empty
        if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
            print(f"Skipping {f} as blast output already exists")
            return
        cmd = f"blastn -outfmt '{fmt}' -query {f} -remote -db core_nt -out {blast_out} -evalue 0.1 -task blastn -dust no -word_size 24 -reward 1 -penalty -3 -max_target_seqs 4"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Blast complete for {f}")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_single_blast, f) for f in splitted_fasta]
        for future in as_completed(futures):
            future.result()  # Raise any exceptions that occurred

def parse_args():
    parser = argparse.ArgumentParser(description="Run BLAST on input fasta file")
    parser.add_argument('--input_file', required=True, help='Path to the input fasta file')
    parser.add_argument('--split_folder', required=True, help='Path to the folder to store split fasta files')
    parser.add_argument('--blast_folder', required=True, help='Path to the folder to store BLAST output')
    parser.add_argument('--max_workers', type=int, default=4, help='Number of concurrent BLAST commands')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if not os.path.exists(args.blast_folder):
        print("Not running blast as the blast folder does not exist")
        sys.exit(0)
    if len(read_fasta(args.input_file)) > SPLIT_THRESH:
        split_fasta(args.input_file, args.split_folder, SPLIT_EACH)
    else:
        shutil.copy(args.input_file, args.split_folder)
    splitted_fasta = [join(args.split_folder, f) for f in os.listdir(args.split_folder) if f.endswith(".fasta")]
    run_blast(splitted_fasta, args.blast_folder, args.max_workers)
