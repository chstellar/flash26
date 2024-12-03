# add new anchors to existing clusters 
# usage: update_anchors.py <cluster_file> <anchor_file> <output_clusters> --threads <num_threads>
# cluster_file: the clusters file to be updated
# anchor_file: the new anchors to be added
# output_clusters: the updated clusters file
# num_threads: the number of threads to use

# import modules
import argparse
from collections import defaultdict
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from multiprocessing import Pool
import Bio.SeqIO as SeqIO
from tqdm import tqdm


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
    for frame in range(3):
        trans = str(SeqIO.Seq(sequence[frame:]).translate(to_stop=True, table=translation_table, cds=False))
        translations.append(trans)
    for frame in range(3):
        trans = str(SeqIO.Seq(sequence.reverse_complement()[frame:]).translate(to_stop=True, table=translation_table, cds=False))
        translations.append(trans)
    return translations


def read_anchors(anchor_file):
    anchors = {}
    for record in SeqIO.parse(anchor_file, "fasta"):
        anchors[record.id] = translate_sequence(str(record.seq))
    return anchors


def find_best_match(clusters, anchor_translations, num_threads, chunk_size=100):
    def process_chunk(chunk):
        chunk_best_matches = {}
        for anchor_id, translations in chunk:
            best_score = 0
            best_cluster = None
            for cluster_id, cluster_seqs in clusters.items():
                for cluster_seq in cluster_seqs:
                    cluster_translations = translate_sequence(cluster_seq)
                    for anchor_translation in translations:
                        for cluster_translation in cluster_translations:
                            alignments = pairwise2.align.globalxx(anchor_translation, cluster_translation)
                            score = alignments[0][2] if alignments else 0
                            if score > best_score:
                                best_score = score
                                best_cluster = cluster_id
            chunk_best_matches[anchor_id] = best_cluster
        return chunk_best_matches

    anchor_items = list(anchor_translations.items())
    chunks = [anchor_items[i:i + chunk_size] for i in range(0, len(anchor_items), chunk_size)]
    
    best_matches = {}
    with Pool(processes=num_threads) as pool:
        for result in tqdm(pool.imap(process_chunk, chunks), total=len(chunks), desc="Processing chunks"):
            best_matches.update(result)
    
    return best_matches

def update_clusters(clusters, best_matches, anchor_translations):
    for anchor_id, cluster_id in best_matches.items():
        if cluster_id:
            clusters[cluster_id].append(anchor_id)
        else:
            clusters[anchor_id] = anchor_translations[anchor_id]
    return clusters

def write_clusters(clusters, output_file):
    with open(output_file, "w") as f:
        for cluster_id, seq_ids in clusters.items():
            for seq_id in seq_ids:
                f.write(f"{cluster_id}\t{seq_id}\n")

def main():
    args = parse_args()
    clusters = read_clusters(args.cluster_file)
    anchor_translations = read_anchors(args.anchor_file)
    best_matches = find_best_match(clusters, anchor_translations, args.threads)
    updated_clusters = update_clusters(clusters, best_matches, anchor_translations)
    write_clusters(updated_clusters, args.output_clusters)

if __name__ == "__main__":
    main()