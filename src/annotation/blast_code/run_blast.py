import os
from os.path import join, basename
import subprocess
import shutil
import sys
import argparse
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import pandas as pd

SPLIT_THRESH = 20  # 100
SPLIT_EACH = 10  # 50
NUCLEOTIDE_SUFFIX_RE = re.compile(r"^(?P<cluster>.+)_[ACGTNacgtn]+$")


def format_taxids(taxid):
    taxid = str(taxid).strip().strip("\"'").strip()
    if taxid in {"", "0", "NA", "NaN", "nan", "None", "none"}:
        return ""

    taxid = taxid.strip("{}[]()")
    taxids = [
        item.strip().strip("\"'").strip()
        for item in taxid.replace(";", ",").replace("+", ",").split(",")
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


def get_record_cluster(record):
    match = NUCLEOTIDE_SUFFIX_RE.match(record.id)
    if match:
        return match.group("cluster")
    cluster_match = re.search(r"(cluster_\d+)", record.id)
    if cluster_match:
        return cluster_match.group(1)
    return record.id


def filter_top_n_per_cluster(fasta_file, output_file, top_n):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    total = len(records)
    if top_n is None or top_n <= 0:
        print(f"Sending {total}/{total} sequences to BLAST query (no per-cluster cap).")
        return fasta_file, total, total

    cluster_counts = {}
    selected = []
    for record in records:
        cluster = get_record_cluster(record)
        count = cluster_counts.get(cluster, 0)
        if count < top_n:
            selected.append(record)
        cluster_counts[cluster] = count + 1

    SeqIO.write(selected, output_file, "fasta")
    print(
        f"Sending {len(selected)}/{total} sequences to BLAST query "
        f"(top {top_n} per cluster; {len(cluster_counts)} clusters observed)."
    )
    return output_file, len(selected), total


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
        for i in range(0, len(records), num_seq):
            output_file = os.path.join(output_dir, f"split_{i}.fasta")
            with open(output_file, "w") as f:
                for record in records[i : i + num_seq]:
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
    parser.add_argument(
        "--top_n_sequences_per_cluster",
        type=int,
        default=0,
        help="Only BLAST the first N significant sequences per cluster. 0 means no cap.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.exists(args.blast_folder):
        print("Not running blast as the blast folder does not exist")
        sys.exit(0)

    query_fasta = args.input_file
    if args.top_n_sequences_per_cluster > 0:
        query_fasta = os.path.join(args.split_folder, "top_n_per_cluster_query.fasta")
        os.makedirs(args.split_folder, exist_ok=True)
    query_fasta, _, _ = filter_top_n_per_cluster(
        args.input_file, query_fasta, args.top_n_sequences_per_cluster
    )

    if len(read_fasta(query_fasta)) > SPLIT_THRESH:
        split_fasta(query_fasta, args.split_folder, SPLIT_EACH)
        if query_fasta != args.input_file and os.path.exists(query_fasta):
            os.remove(query_fasta)
    else:
        if os.path.dirname(os.path.abspath(query_fasta)) != os.path.abspath(args.split_folder):
            shutil.copy(query_fasta, args.split_folder)
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
