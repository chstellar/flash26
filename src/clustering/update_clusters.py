# add new anchors to existing clusters 
# usage: update_anchors.py <cluster_file> <anchor_file> <output_clusters> --threads <num_threads>
# cluster_file: the clusters file to be updated
# anchor_file: the new anchors to be added
# output_clusters: the updated clusters file
# num_threads: the number of threads to use

# import modules
import argparse
from collections import defaultdict
from multiprocessing import Pool
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
import random
from collections import Counter


# define functions
def parse_args():
    parser = argparse.ArgumentParser(description="Add new anchors to existing clusters")
    parser.add_argument("cluster_file", type=str, help="the clusters file to be updated")
    parser.add_argument("anchor_file", type=str, help="the new anchors to be added")
    parser.add_argument("output_clusters", type=str, help="the updated clusters file")
    parser.add_argument("--threads", type=int, default=1, help="the number of threads to use")
    return parser.parse_args()


def read_clusters(cluster_file):
    clusters = defaultdict(list)
    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            clusters[line[0]].append(line[1])
    return clusters


def translate_sequence(sequence, translation_table=1):
    translations = []
    sequence = Seq(sequence)
    for frame in range(3):
        trans = str(Seq(sequence[frame:]).translate(to_stop=True, table=translation_table, cds=False))
        translations.append(trans)
    for frame in range(3):
        trans = str(Seq(sequence.reverse_complement()[frame:]).translate(to_stop=True, table=translation_table, cds=False))
        translations.append(trans)
    translations = [trans for trans in translations if len(trans) >= 8]
    return translations


def read_anchors(anchor_file):
    anchors = []
    if anchor_file.endswith(".fasta") or anchor_file.endswith(".fa"):
        for record in SeqIO.parse(anchor_file, "fasta"):
            anchors.append(str(record.seq))
    else:
        with open(anchor_file, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                anchors.append(line[0])
    return anchors


def create_hashed_cluster_dict(clusters, m=4, N=300):
    """
    for each anchor, create a dictionary with the cluster_id as the value and the anchor as the key. 
    Instead of inserting the anchor insert a the anchor with m characters masked as N. Do this for each 
    anchor N times, so we have a dictionary that is num_anchors * N in size.
    """
    hashed_clusters = defaultdict(list)
    for cluster_id, anchors in clusters.items():
        for anchor in anchors:
            translated_anchors = translate_sequence(anchor)
            for tran_anch in translated_anchors:
                for i in range(N):
                    masked_anchor = list(tran_anch)
                    indices = random.sample(range(len(tran_anch)), m)
                    for index in indices:
                        masked_anchor[index] = "N"
                    masked_anchor = "".join(masked_anchor)
                    if cluster_id not in hashed_clusters[masked_anchor]:
                        hashed_clusters[masked_anchor].append(cluster_id)
    return hashed_clusters

def add_anchor_to_cluster(anchor, cluster_dict, m=4, N=300):
    """
    For each anchor, translate it and mask m characters as Ns. Do this N times. For each masked anchor,
    check if it is in the cluster dictionary. If it is, add the cluster_id to the list of clusters for that anchor.
    At the end return the anchor and the list of clusters it belongs to.
    """
    anchor_translations = translate_sequence(anchor)
    cluster_ids = list()
    for tran_anch in anchor_translations:
        for i in range(N):
            masked_anchor = list(tran_anch)
            indices = random.sample(range(len(tran_anch)), m)
            for index in indices:
                masked_anchor[index] = "N"
            masked_anchor = "".join(masked_anchor)
            cur_clusters = cluster_dict.get(masked_anchor, [])
            cluster_ids.extend(cur_clusters)
    count = Counter(cluster_ids)
    if len(count) == 0:
        return anchor, None
    elif Counter(cluster_ids).most_common(1)[0][1] <= 3:
        return anchor, None
    else:
        most_common_cluster = Counter(cluster_ids).most_common(1)[0][0]
        return anchor, most_common_cluster


# def find_best_match(clusters, anchors, threads, m=4, N=300):
#     """
#     For each anchor in the anchor list, find the most common cluster it belongs to. This is the cluster
#     that we will assign the anchor to for downstream processing.
#     """
#     cluster_lookup = create_hashed_cluster_dict(clusters, m=m, N=N)
#     cluster_assignments = dict()
#     with Pool(threads) as p:
#         results = list(tqdm(p.imap(lambda x: add_anchor_to_cluster(x, cluster_lookup, m=m, N=N), anchors), total=len(anchors)))
#     for anchor, cluster_id in results:
#         if cluster_id is not None:
#             cluster_assignments[anchor] = cluster_id
#     return cluster_assignments


def main():
    m = 4
    N = 300
    args = parse_args()
    clusters = read_clusters(args.cluster_file)
    anchors = read_anchors(args.anchor_file)
    cluster_lookup = create_hashed_cluster_dict(clusters, m=m, N=N)
    cluster_assignments = dict()
    with Pool(args.threads) as p:
        results = list(tqdm(p.imap(lambda x: add_anchor_to_cluster(x, cluster_lookup, m=m, N=N), anchors), total=len(anchors)))
    for anchor, cluster_id in results:
        if cluster_id is not None:
            cluster_assignments[anchor] = cluster_id
    with open(args.output_clusters, "w") as f:
        for anchor, cluster_id in cluster_assignments.items():
            if cluster_id is not None:
                f.write(f"{cluster_id}\t{anchor}\n")
            else:
                f.write(f"NA\t{anchor}\n")


if __name__ == "__main__":
    main()
