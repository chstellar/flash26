import os
from os.path import join, basename
import subprocess
import shutil
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import pandas as pd

SPLIT_THRESH = 20  # 100
SPLIT_EACH = 10  # 50


def format_taxids(taxid):
    taxid = str(taxid).strip().strip("\"'").strip()
    if taxid in {"", "0", "NA", "NaN", "nan", "None", "none"}:
        return ""

    taxid = taxid.strip("{}[]()")
    taxids = [
        item.strip().strip("\"'").strip()
        for item in taxid.replace(";", ",").split(",")
        if item.strip()
    ]
    if not taxids:
        return ""

    invalid_taxids = [item for item in taxids if not item.isdigit()]
    if invalid_taxids:
        raise ValueError(f"Invalid taxid value(s): {', '.join(invalid_taxids)}")

    return f"-taxids {','.join(taxids)}"


def read_fasta(fasta_file, output_type="dict"):
    """
    Read a fasta file and return a dictionary with the sequence id as key and the sequence as value.
    """
    if output_type == "list":
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
        return sequences
    elif output_type == "dict":
        sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = record.seq
        return sequences
    elif output_type == "pandas":
        sequences = []
        description = []
        ids = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
            description.append(record.description)
            ids.append(record.id)
        return pd.DataFrame(
            {"ID": ids, "Description": description, "Sequence": sequences}
        )


def split_fasta(fasta_file, output_dir, num_seq=1):
    """
    Split a fasta file into multiple files.
    """
    os.makedirs(output_dir, exist_ok=True)
    if num_seq == 1:
        for record in SeqIO.parse(fasta_file, "fasta"):
            output_file = os.path.join(output_dir, record.id + ".fasta")
            with open(output_file, "w") as f:
                f.write(">" + record.description + "\n")
                f.write(str(record.seq) + "\n")
    else:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        num_files = len(records) // num_seq
        for i in range(num_files):
            output_file = os.path.join(output_dir, f"split_{i}.fasta")
            with open(output_file, "w") as f:
                for record in records[i * num_seq : (i + 1) * num_seq]:
                    f.write(">" + record.description + "\n")
                    f.write(str(record.seq) + "\n")


def run_blast(splitted_fasta, blast_folder, max_workers, taxid, local_blast_db=""):
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue qcovs sgi sacc slen staxids stitle"
    taxid = format_taxids(taxid)
    if local_blast_db:
        remote_flag = ""
        local_db_export = f"export BLASTDB={local_blast_db}; "
    else:
        taxid = ""  # cannot use -remote with taxids, so we skip taxid
        remote_flag = "-remote"
        local_db_export = ""

    def run_single_blast(f):
        if taxid:
            print(f"Using taxonomy flag: {taxid}")

        blast_out = join(blast_folder, basename(f).split(".")[0] + ".blastout.tsv")

        # skip if tsv file already exists and is not empty
        if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
            print(f"Skipping {f} as blast output already exists")
            return

        # create the blast command
        params = [
            local_db_export,
            "blastn",
            f"-outfmt '{fmt}'",
            f"-query {f}",
            remote_flag,
            "-db core_nt",
            f"-out {blast_out}",
            "-evalue 0.1",
            "-task blastn",
            "-dust no",
            "-word_size 24",
            "-reward 1",
            "-penalty -3",
            taxid,
            "-max_target_seqs 5",
        ]
        # join only non-empty parts to avoid extra spaces
        cmd = " ".join(p for p in params if p)

        subprocess.run(cmd, shell=True, check=True)
        print(f"Blast complete for {f}")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_single_blast, f) for f in splitted_fasta]
        for future in as_completed(futures):
            future.result()  # Raise any exceptions that occurred


def parse_args():
    parser = argparse.ArgumentParser(description="Run BLAST on input fasta file")
    parser.add_argument(
        "--input_file", required=True, help="Path to the input fasta file"
    )
    parser.add_argument(
        "--split_folder",
        required=True,
        help="Path to the folder to store split fasta files",
    )
    parser.add_argument(
        "--blast_folder", required=True, help="Path to the folder to store BLAST output"
    )
    parser.add_argument(
        "--max_workers", type=int, default=4, help="Number of concurrent BLAST commands"
    )
    parser.add_argument(
        "--taxid",
        type=str,
        default="0",
        help="What tax id to restrict to when searching BLAST",
    )
    parser.add_argument(
        "--local_blast_db",
        type=str,
        default="",
        help="Path to local BLAST database folder (if using local databases)",
    )
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
    splitted_fasta = [
        join(args.split_folder, f)
        for f in os.listdir(args.split_folder)
        if f.endswith(".fasta")
    ]
    run_blast(
        splitted_fasta,
        args.blast_folder,
        args.max_workers,
        args.taxid,
        args.local_blast_db,
    )
