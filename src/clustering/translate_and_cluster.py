import numpy as np
import pandas as pd
from Bio.Seq import Seq
import os
import argparse
import glob

# Parse arguments
def parse_args():
    """
    Daniel Cotter.
    """
    parser = argparse.ArgumentParser(
        description="Translate nucleotide sequences in a FASTA file to amino acid sequences."
    )
    parser.add_argument(
        "-t",
        "--translation_table",
        default=11,
        type=int,
        help="Translation table to use for translating nucleotide sequences to amino acid sequences.",
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input FASTA file containing nucleotide sequences to be translated and clustered.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file reporting input, output, AA sequences and their cluster IDs.",
    )
    parser.add_argument(
        "--strand",
        default="sense",
        type=str,
        help="Strandedness of input sequences; options: [ sense , antisense , both ].",
    )
    parser.add_argument(
        "--temp_dir",
        type=str,
        help="This folder will be created and moved-to.",
    )
    return parser.parse_args()

args = parse_args()

os.system('rm -r '+ args.temp_dir)
os.mkdir(args.temp_dir)
os.chdir(args.temp_dir)

"""
Load FASTA.
"""
input_file = pd.read_csv(args.input,header=None)
input_file = input_file[input_file[0].str[0]!='>']\
    .reset_index(drop=True)\
    .rename(columns={0:'from_nt'})

"""
Assistant.
"""
def reverse_complement(sequence):
    # Define a dictionary for complementing the bases
    complement = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A' }
    # Create the reverse complement using a list comprehension
    return ''.join(complement[base] for base in reversed(sequence))
def translate_sequence(nucleotide_seq, translation_table):
    """
    Daniel Cotter.
    Translates a nucleotide sequence to an amino acid sequence
    using the specified translation table
    """
    seq = Seq(nucleotide_seq)
    return ''.join(seq.translate(table=translation_table))
def detranslate_sequence(sequence):
    # For converting from AA to a stable codon representation. 
    dtd = \
    { 
        "A" :  "AAA"  , "L" : "AAC",
        "R" :  "AAG"  , "K" : "AAT",
        "N" :  "ACA"  , "M" : "ACC",
        "D" :  "ACG"  , "F" : "ACT", 
        "C" :  "ATA"  , "P" : "ATC",
        "E" :  "ATG"  , "S" : "ATT",
        "Q" :  "AGA"  , "T" : "AGC",
        "G" :  "AGG"  , "W" : "AGT",
        "H" :  "TAA"  , "Y" : "TAC",
        "I" :  "TAG"  , "V" : "TAT", "*" : "CAA"
    }
    return ''.join(dtd.get(base,base) for base in sequence)

"""
Perform translations.
"""
if args.strand == "sense": 
    cols_to_operate = ['from_nt']
if args.strand == "antisense":
    cols_to_operate = ['from_nt_rc']
if args.strand == "both": 
    cols_to_operate = ['from_nt','from_nt_rc']
if 'from_nt_rc' in cols_to_operate: 
    input_file['from_nt_rc'] = input_file['from_nt'].apply(lambda x: reverse_complement(x))
record = pd.DataFrame({'from_nt':[], 
                       'frame':[],
                       'aa_seq':[]})
### For each frame: 
for i in [0,1,2]: 
    ### For each strand, as requested:
    for column in cols_to_operate:
        record_i = pd.DataFrame()
        ### Convert to aa_seq using the desired translation table. 
        record_i['aa_seq'] = \
            input_file[column].str[i:]\
                .apply(
                lambda x: translate_sequence(x, args.translation_table)
                )
        ### Report the frame.
        record_i['frame'] = i + 1
        if '_rc' in column:
            record_i['frame'] = record_i['frame'] * -1
        record_i['from_nt'] = input_file['from_nt']
        ### Collect results in a wide-form matrix.
        record = pd.concat([record,record_i])
record = record.reset_index(drop=True)
record['frame'] = record['frame'].astype(int)
del record_i

"""
Perform de-translations.
"""
nonredundant_aa = pd.DataFrame(record['aa_seq'].drop_duplicates())
nonredundant_aa['to_nt'] = nonredundant_aa['aa_seq'].apply(lambda x: detranslate_sequence(x))
record = record.merge(nonredundant_aa,how='left')
del nonredundant_aa
record['aa_seq'] = record['aa_seq'].str.replace('*','X')

#### Write nucleotide sequences to a FASTA. 
os.system('rm detranslated_sequences.fasta')
file = open('detranslated_sequences.fasta','a')
for i in record['to_nt']:
    file.write('>'+i+'\n'+i+'\n')
file.close()

#### Run MMseqs2 easy-clust.
os.system('rm -r tmp')
os.system("/oak/stanford/groups/horence/dcotter1/software/mmseqs/bin/mmseqs easy-cluster detranslated_sequences.fasta nt_cluster nt_tmp")

#### Write amino-acid sequences to a FASTA. 
os.system('rm translated_sequences.fasta')
file = open('translated_sequences.fasta','a')
for i in record['aa_seq']:
    file.write('>'+i+'\n'+i+'\n')
file.close()

#### Run MMseqs2 easy-clust.
os.system('rm -r tmp')
os.system("/oak/stanford/groups/horence/dcotter1/software/mmseqs/bin/mmseqs easy-cluster translated_sequences.fasta aa_cluster aa_tmp")


"""
Read in the clustering result. Add cluster id's. Merge with the existing records.
"""
clust_out = pd.read_csv('nt_cluster_cluster.tsv',sep='\t',header=None)
clust_out = clust_out.rename(columns={0:'nt_centroid',1:'to_nt'})
t = clust_out[['nt_centroid']].drop_duplicates()
t['cluster_id'] = range(len(t))
clust_out = clust_out.merge(t,how='left')
record = record.merge(clust_out)
del t

if 'aa_cluster_cluster.tsv' in glob.glob("*tsv"): 
    clust_out = pd.read_csv('aa_cluster_cluster.tsv',sep='\t',header=None)
    clust_out = clust_out.rename(columns={0:'aa_centroid',1:'aa_seq'})
    t = clust_out[['aa_centroid']].drop_duplicates()
    t['aa_cluster_id'] = range(len(t))
    clust_out = clust_out.merge(t,how='left')
    record = record.merge(clust_out)
    del t  

"""
Write the records (from_nt, frame, aa_seq, to_nt, centroid, cluster_id). 
"""
# First filter the records to only include the from_nt in a single cluster.
record = record.sort_values(by=['from_nt','cluster_id'])
record = record.drop_duplicates(subset=['from_nt'],keep='first')
record = record.reset_index(drop=True)
record = record.sort_values(by=['cluster_id'])
record = record.reset_index(drop=True)

# sort clusters by number of from_nt members.
cluster_sizes = record.groupby('cluster_id').size()
cluster_sizes = cluster_sizes.sort_values(ascending=False)
cluster_sizes = cluster_sizes.reset_index()
cluster_sizes = cluster_sizes.rename(columns={0:'size'})
record = record.merge(cluster_sizes)
record = record.sort_values(by=['size','cluster_id'])
# filter to only include clusters with more than 1 member.
record = record[record['size']>1]
record = record.reset_index(drop=True)

# only write cluster id and from_nt to file
record = record[['cluster_id','from_nt']]

# write to file Drop headers
record.to_csv(args.output,header=False,index=False,sep='\t')
