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


def get_first_coef_abs(coef_string):
    coef_string = str(coef_string).strip().strip("[]")
    if not coef_string:
        return 0.0
    first_coef = coef_string.split(",")[0].strip()
    try:
        return abs(float(first_coef))
    except ValueError:
        return 0.0


def get_feature_cluster(feature):
    cluster_match = re.search(r"(cluster_\d+)", str(feature))
    if cluster_match:
        return cluster_match.group(1)
    kmer_match = re.search(r"(\w+_kmer_\d+)", str(feature))
    if kmer_match:
        return kmer_match.group(1)
    return None


def get_cluster_index(cluster):
    cluster_match = re.search(r"cluster_(\d+)", str(cluster))
    if cluster_match:
        return int(cluster_match.group(1))
    return None


def get_plot_selected_clusters(coefficients_file, num_hits):
    coefficients = pd.read_csv(coefficients_file, sep="\t")
    required_columns = {"metadata_category", "feature", "coefficients"}
    missing_columns = required_columns - set(coefficients.columns)
    if missing_columns:
        raise ValueError(
            "Coefficient file is missing required columns for plot-selected BLAST: "
            + ", ".join(sorted(missing_columns))
        )

    coefficients = coefficients.copy()
    coefficients["cluster"] = coefficients["feature"].apply(get_feature_cluster)
    coefficients = coefficients.dropna(subset=["cluster"])
    coefficients["max_coefficient"] = coefficients["coefficients"].apply(get_first_coef_abs)

    selected_clusters = set()
    for _, category_dt in coefficients.groupby("metadata_category", sort=False):
        top_dt = (
            category_dt[["cluster", "feature", "max_coefficient"]]
            .drop_duplicates()
            .sort_values("max_coefficient", ascending=False)
            .head(num_hits)
        )
        selected_clusters.update(top_dt["cluster"].astype(str))
    return selected_clusters


def get_observed_sequences_by_cluster(sample_sequences_file, selected_clusters, cluster_length):
    if not sample_sequences_file:
        return None

    observed = {cluster: set() for cluster in selected_clusters}
    cluster_indices = {
        cluster: get_cluster_index(cluster)
        for cluster in selected_clusters
    }
    cluster_indices = {
        cluster: index
        for cluster, index in cluster_indices.items()
        if index is not None
    }
    if not cluster_indices:
        return None

    sample_sequences = pd.read_csv(sample_sequences_file, sep="\t")
    if sample_sequences.shape[1] < 2:
        raise ValueError(
            f"Sample sequence file must have at least two columns: {sample_sequences_file}"
        )
    sequence_col = sample_sequences.columns[1]
    concatenated_sequences = sample_sequences[sequence_col].astype(str)
    for cluster, index in cluster_indices.items():
        start = index * cluster_length
        end = start + cluster_length
        observed[cluster].update(
            seq[start:end]
            for seq in concatenated_sequences
            if len(seq) >= end
        )
        observed[cluster].discard("")
    return observed


def filter_plot_selected_clusters(
    fasta_file,
    output_file,
    coefficients_file,
    num_hits,
    sample_sequences_file="",
    cluster_length=0,
):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    total = len(records)
    selected_clusters = get_plot_selected_clusters(coefficients_file, num_hits)
    observed_sequences = None
    if sample_sequences_file and cluster_length > 0:
        observed_sequences = get_observed_sequences_by_cluster(
            sample_sequences_file, selected_clusters, cluster_length
        )

    selected = []
    for record in records:
        cluster = get_record_cluster(record)
        if cluster not in selected_clusters:
            continue
        if observed_sequences is not None:
            cluster_observed = observed_sequences.get(cluster, set())
            if cluster_observed and str(record.seq) not in cluster_observed:
                continue
        selected.append(record)

    SeqIO.write(selected, output_file, "fasta")
    observed_message = ""
    if observed_sequences is not None:
        observed_message = "; restricted to observed sample-sequence variants"
    print(
        f"Sending {len(selected)}/{total} sequences to BLAST query "
        f"(plot-selected clusters: {len(selected_clusters)} clusters from top {num_hits} plotted coefficients per metadata category{observed_message})."
    )
    return output_file, len(selected), total


def filter_plot_selected_and_top_clusters(
    fasta_file,
    output_file,
    coefficients_file,
    num_hits,
    top_n,
    sample_sequences_file="",
    cluster_length=0,
):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    total = len(records)
    selected_clusters = get_plot_selected_clusters(coefficients_file, num_hits)
    observed_sequences = None
    if sample_sequences_file and cluster_length > 0:
        observed_sequences = get_observed_sequences_by_cluster(
            sample_sequences_file, selected_clusters, cluster_length
        )

    selected = []
    nonplot_cluster_counts = {}
    plot_selected_count = 0
    top_n_selected_count = 0
    for record in records:
        cluster = get_record_cluster(record)
        if cluster in selected_clusters:
            if observed_sequences is not None:
                cluster_observed = observed_sequences.get(cluster, set())
                if cluster_observed and str(record.seq) not in cluster_observed:
                    continue
            selected.append(record)
            plot_selected_count += 1
            continue

        if top_n is None or top_n <= 0:
            continue
        count = nonplot_cluster_counts.get(cluster, 0)
        if count < top_n:
            selected.append(record)
            top_n_selected_count += 1
        nonplot_cluster_counts[cluster] = count + 1

    SeqIO.write(selected, output_file, "fasta")
    observed_message = ""
    if observed_sequences is not None:
        observed_message = "; plot-selected clusters restricted to observed sample-sequence variants"
    print(
        f"Sending {len(selected)}/{total} sequences to BLAST query "
        f"(plot-selected-and-top mode: {plot_selected_count} plot-selected sequences from "
        f"{len(selected_clusters)} plotted clusters{observed_message}; "
        f"{top_n_selected_count} fallback sequences from non-plotted clusters, top {top_n} per cluster)."
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
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue qcovs sgi sacc slen staxids sscinames stitle"
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
    parser.add_argument(
        "--blast_selection_mode",
        type=str,
        default="all",
        choices=[
            "all",
            "top_n_per_cluster",
            "plot_selected",
            "plot_selected_and_top",
        ],
        help="Which significant sequences to send to BLAST.",
    )
    parser.add_argument(
        "--coefficients",
        type=str,
        default="",
        help="Nonzero coefficient TSV. Required for plot-selected BLAST modes.",
    )
    parser.add_argument(
        "--num_plot_hits",
        type=int,
        default=10,
        help="Number of top coefficient clusters per metadata category used by the BLAST plot.",
    )
    parser.add_argument(
        "--sample_sequences",
        type=str,
        default="",
        help="Prepared sample sequence TSV used to restrict plot-selected BLAST to observed dot sequences.",
    )
    parser.add_argument(
        "--cluster_length",
        type=int,
        default=0,
        help="Concatenated anchor-target sequence length used with --sample_sequences.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.exists(args.blast_folder):
        print("Not running blast as the blast folder does not exist")
        sys.exit(0)

    query_fasta = args.input_file
    if args.blast_selection_mode == "plot_selected":
        if not args.coefficients:
            raise ValueError("--coefficients is required with --blast_selection_mode plot_selected")
        query_fasta = os.path.join(args.split_folder, "plot_selected_query.fasta")
        os.makedirs(args.split_folder, exist_ok=True)
        query_fasta, _, _ = filter_plot_selected_clusters(
            args.input_file,
            query_fasta,
            args.coefficients,
            args.num_plot_hits,
            args.sample_sequences,
            args.cluster_length,
        )
    elif args.blast_selection_mode == "plot_selected_and_top":
        if not args.coefficients:
            raise ValueError(
                "--coefficients is required with --blast_selection_mode plot_selected_and_top"
            )
        query_fasta = os.path.join(args.split_folder, "plot_selected_and_top_query.fasta")
        os.makedirs(args.split_folder, exist_ok=True)
        query_fasta, _, _ = filter_plot_selected_and_top_clusters(
            args.input_file,
            query_fasta,
            args.coefficients,
            args.num_plot_hits,
            args.top_n_sequences_per_cluster,
            args.sample_sequences,
            args.cluster_length,
        )
    elif (
        args.blast_selection_mode == "top_n_per_cluster"
        or args.top_n_sequences_per_cluster > 0
    ):
        query_fasta = os.path.join(args.split_folder, "top_n_per_cluster_query.fasta")
        os.makedirs(args.split_folder, exist_ok=True)
        query_fasta, _, _ = filter_top_n_per_cluster(
            args.input_file, query_fasta, args.top_n_sequences_per_cluster
        )
    else:
        query_fasta, _, _ = filter_top_n_per_cluster(args.input_file, query_fasta, 0)

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
