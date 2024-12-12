# add new anchors to existing clusters 
# usage: update_anchors.py <cluster_file> <anchor_file> <output_clusters> --threads <num_threads>
# cluster_file: the clusters file to be updated
# anchor_file: the new anchors to be added
# output_clusters: the updated clusters file
# num_threads: the number of threads to use

# import modules
import argparse
from multiprocessing import Pool
import Bio.SeqIO as SeqIO
import random


# define functions
def parse_args():
    """
    Parse the command line arguments and return the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Add new anchors to existing clusters")
    parser.add_argument("anchor_file", type=str, help="the new anchors to be added")
    parser.add_argument("output_clusters", type=str, help="the updated clusters file")
    return parser.parse_args()

def read_anchors(anchor_file):
    """
    Read in the anchors from a fasta or tab-delimited file and return a list of the anchors.
    """
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

def cluster_anchors(anchors, m=4, N=300, j=5):
    """
    Create a dictionary of clusters, then shuffle the input list of anchors and assign them to clusters.
    Fore each anchor, mask m random characters as N and check to see if the masked anchor is in the cluster 
    dictionary. If it is, add the anchor to the cluster. If it is not, assign the anchor to a new cluster.
    Additionally drop the first and last 1:j charaters from each anchor and check that those are in the cluster dictionary.
    This can account for shifts in the anchor sequences.
    m: the number of characters to mask
    N: the number of masked anchors to generate
    j: the number of characters to drop up to from the beginning and end of each anchor
    """
    lookup_dict = dict()
    clusters = dict()
    for id, anchor in enumerate(anchors):
        if id % 1000 == 0:
            print(f"Processing anchor {id}/{len(anchors)}")
        masked_anchors = []
        for i in range(N):
            masked_anchor = list(anchor)
            indices = random.sample(range(len(anchor)), m)
            for idx in indices:
                masked_anchor[idx] = "N"
            masked_anchor = "".join(masked_anchor)
            if masked_anchor in lookup_dict:
                cluster_id = lookup_dict[masked_anchor]
                clusters[cluster_id].append(anchor)
                break
        for i in range(j):
            i = i + 1
            front_trimmed = anchor[i:len(anchor)]
            back_trimmed = anchor[0:len(anchor)-i]
            if front_trimmed in lookup_dict:
                cluster_id = lookup_dict[front_trimmed]
                clusters[cluster_id].append(anchor)
                break
            if back_trimmed in lookup_dict:
                cluster_id = lookup_dict[back_trimmed]
                clusters[cluster_id].append(anchor)
                break
            masked_anchors.append(front_trimmed)
            masked_anchors.append(back_trimmed)
        # if there are no matches, create a new cluster
        # update the lookup dictionary with the masked anchors and a new cluster id
        for masked_anchor in masked_anchors:
            lookup_dict[masked_anchor] = len(clusters)
        clusters[len(clusters)] = [anchor]
    # return the clusters dictionary with the cluster id as the key and the list of anchors as the value
    return clusters


def main():
    m = 4
    N = 1000
    j = 5
    args = parse_args()
    print("Reading anchors...")
    anchors = read_anchors(args.anchor_file)
    random.shuffle(anchors)
    print("Clustering anchors...")
    clusters = cluster_anchors(anchors, m, N, j)
    print("Writing clusters...")
    with open(args.output_clusters, "w") as f:
        for cluster_id, cluster in clusters.items():
            for anchor in cluster:
                f.write(f"{cluster_id}\t{anchor}\n")


if __name__ == "__main__":
    main()
