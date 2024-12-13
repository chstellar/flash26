# add new anchors to existing clusters 
# usage: update_anchors.py <cluster_file> <anchor_file> <output_clusters> --threads <num_threads>
# cluster_file: the clusters file to be updated
# anchor_file: the new anchors to be added
# output_clusters: the updated clusters file
# num_threads: the number of threads to use

# usage: python <fuzzy_lookup_clustering_AA.py> <anchor_file> <output_clusters> --translation_table <translation_table> --protein_db <protein_db>

# import modules
import argparse
from multiprocessing import Pool
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
import random
from fuzzysearch import find_near_matches

# define functions
def parse_args():
    """
    Parse the command line arguments and return the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Add new anchors to existing clusters")
    parser.add_argument("anchor_file", type=str, help="the new anchors to be added")
    parser.add_argument("output_clusters", type=str, help="the updated clusters file")
    parser.add_argument("--translation_table", type=int, default=1, help="the translation table to use")
    parser.add_argument("--protein_db", type=str, help="a protein database to filter the translations", default=None)
    return parser.parse_args()

def read_anchors(anchor_file):
    """
    Read in the anchors from a fasta or tab-delimited file and return a list of the anchors.
    """
    anchors = []
    if anchor_file.endswith(".fasta") or anchor_file.endswith(".fa"):
        # Read fasta file
        for record in SeqIO.parse(anchor_file, "fasta"):
            anchors.append(str(record.seq))
    else:
        # Read tab-delimited file
        with open(anchor_file, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                anchors.append(line[0])
    return anchors

def translate_anchor(anchor, translation_table=1, protein_db=None):
    """
    Translate the anchor sequence in all 6 reading frames and return the translated sequences.
    If a protein database is provided, filter out the translations that are not found in the database.
    """
    translations = []
    anchor = Seq(anchor)
    for frame in range(3):
        # Translate in forward direction
        translated = str(Seq(anchor[frame:]).translate(to_stop=True, table=translation_table, cds=False))
        # Translate in reverse direction
        translated_reverse = str(Seq(anchor.reverse_complement()[frame:]).translate(to_stop=True, table=translation_table, cds=False))
        translations.append(translated)
        translations.append(translated_reverse)
    # Filter out translations that are less than 7 characters
    translations = [translated for translated in translations if len(translated) >= 7]
    # Filter out translations that are not found in a protein database
    if protein_db:
        out_translations = []
        for protein in translations:
            matches = find_near_matches(protein.encode(), protein_db, max_l_dist=1)
            if matches:
                out_translations.append(protein)
    else:
        out_translations = translations
    return out_translations

def cluster_anchors(anchors, m=3, N=300, j=2, translation_table=1, protein_db=None):
    """
    Create a dictionary of clusters, then shuffle the input list of anchors and assign them to clusters.
    For each anchor, mask m random characters as N and check to see if the masked anchor is in the cluster 
    dictionary. If it is, add the anchor to the cluster. If it is not, assign the anchor to a new cluster.
    Additionally drop the first and last 1:j characters from each anchor and check that those are in the cluster dictionary.
    This can account for shifts in the anchor sequences.
    m: the number of characters to mask
    N: the number of masked anchors to generate
    j: the number of characters to drop up to from the beginning and end of each anchor
    """
    lookup_dict = dict()  # Dictionary to store masked anchors and their cluster IDs
    clusters = dict()  # Dictionary to store clusters
    aa_matches = dict()  # List to store anchor and translation matches
    for id, anchor in enumerate(anchors):
        if id % 10 == 0:
            print(f"Processing anchor {id}/{len(anchors)}")
        translations = translate_anchor(anchor, translation_table=translation_table, protein_db=protein_db)
        if not translations:
            #print(f"No translations found for anchor {anchor}")
            continue
        masked_anchors = []
        found_cluster = False
        for tran_anch in translations:
            if not found_cluster:
                for i in range(N):
                    masked_anchor = list(tran_anch)
                    indices = random.sample(range(len(tran_anch)), m)
                    for idx in indices:
                        masked_anchor[idx] = "N"
                    masked_anchor = "".join(masked_anchor)
                    if masked_anchor in lookup_dict:
                        cluster_id = lookup_dict[masked_anchor]
                        clusters[cluster_id].append(anchor)
                        aa_matches[cluster_id].append([anchor, tran_anch])
                        found_cluster = True
                        masked_anchors.append(masked_anchor)
                        # update the lookup dictionary with all the masked anchors so far
                        for mask_anch in masked_anchors:
                            lookup_dict[mask_anch] = cluster_id
                        break
                    masked_anchors.append(masked_anchor)
                if found_cluster:
                    break
                for i in range(j):
                    i = i + 1
                    front_trimmed = tran_anch[i:len(tran_anch)]
                    back_trimmed = tran_anch[0:len(tran_anch)-i]
                    if front_trimmed in lookup_dict:
                        cluster_id = lookup_dict[front_trimmed]
                        clusters[cluster_id].append(anchor)
                        aa_matches[cluster_id].append([anchor, tran_anch])
                        found_cluster = True
                        masked_anchors.append(front_trimmed)
                        # update the lookup dictionary with all the masked anchors so far
                        for mask_anch in masked_anchors:
                            lookup_dict[mask_anch] = cluster_id
                        break
                    if back_trimmed in lookup_dict:
                        cluster_id = lookup_dict[back_trimmed]
                        clusters[cluster_id].append(anchor)
                        aa_matches[cluster_id].append([anchor, tran_anch])
                        found_cluster = True
                        masked_anchors.append(back_trimmed)
                        # update the lookup dictionary with all the masked anchors so far
                        for mask_anch in masked_anchors:
                            lookup_dict[mask_anch] = cluster_id
                        break
                    masked_anchors.append(front_trimmed)
                    masked_anchors.append(back_trimmed)
        if not found_cluster:
            for masked_anchor in masked_anchors:
                lookup_dict[masked_anchor] = len(clusters)
            aa_matches[len(clusters)] = [[anchor, ";".join(translations)]]
            clusters[len(clusters)] = [anchor]

    # Return the clusters dictionary with the cluster id as the key and the list of anchors as the value
    return clusters, aa_matches

def main():
    m = 3
    N = 300
    j = 2
    args = parse_args()
    print("Reading anchors...")
    anchors = read_anchors(args.anchor_file)
    print("Clustering anchors...")
    if args.protein_db:
        with open(args.protein_db, "rb") as f:
            database = f.read()
            # remove newline characters
            database = database.replace(b"\n", b"")
        clusters, aa_matches = cluster_anchors(anchors, m, N, j, translation_table=args.translation_table, protein_db=database)
    else:
        clusters, aa_matches = cluster_anchors(anchors, m, N, j, translation_table=args.translation_table)
    print("Writing clusters...")
    with open(args.output_clusters, "w") as f:
        for cluster_id, cluster in clusters.items():
            for anchor in cluster:
                f.write(f"{cluster_id}\t{anchor}\n")
    # Replace the file extension with _aa_matches.tsv
    aa_matches_file = args.output_clusters.split(".")[0] + "_aa_matches.tsv"
    with open(aa_matches_file, "w") as f:
        for cluster_id, matches in aa_matches.items():
            for match in matches:
                f.write(f"{cluster_id}\t{match[0]}\t{match[1]}\n")

if __name__ == "__main__":
    main()
