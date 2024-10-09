# id_kmers_per_cluster.R
# Daniel Cotter
# 2024-09-18

# This script takes in the ordering file and the kmers fasta file and 
# outputs a script in long format with the cluster number, the kmer id, and the kmer sequence.


## import packages --------
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))



## parse arguments --------
# define command line arguments
# define command line arguments
option_list <- list(
  make_option(c("-o", "--ordering"), help="Ordering file", type="character"),
  make_option(c("-k", "--kmers"), help="Kmers fasta file", type="character"),
  make_option(c("-p", "--output"), help="Output file.", type="character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$ordering) | is.null(opt$output) | !file.exists(opt$kmers)) {
  stop("Must provide kmers, ordering, and output file")
}

## print a summary of the arguments
cat("\n####################\n")
cat("Running id_kmers_per_cluster.R with the following arguments:\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Kmers file: ", opt$kmers, "\n")
cat("Output file: ", opt$output, "\n")
cat("####################\n\n")

# load the ordering file
cat("Loading the ordering file...\n")
ordering <- fread(opt$ordering, header=F, sep="\t", 
                  col.names = c("sample_name", "seq", "kmer", "start", "end")) 
ordering <- ordering %>% select(sample_name, kmer)

# labeling clusters in the ordering file (they appear one per row per sample)
ordering <- ordering %>% group_by(sample_name) %>% mutate(cluster=row_number()) %>%
  mutate(cluster=cluster-1) %>% ungroup()

cluster_to_kmer_mapping <- ordering %>% select(cluster, kmer) %>% 
  distinct()
cluster_to_kmer_mapping <- cluster_to_kmer_mapping %>% 
  arrange(cluster, kmer) %>%
  mutate(cluster=paste("cluster", cluster, sep="_"))

# read in the kmers file using Biostrings
dna <- readDNAStringSet(opt$kmers)
dna <- data.frame(kmer=names(dna), seq=as.character(dna), stringsAsFactors=F)

# join the sequences to the cluster to kmer mapping
cluster_to_kmer_mapping <- cluster_to_kmer_mapping %>% 
  left_join(dna, by="kmer")

# write out the cluster to kmer mapping
cluster_to_kmer_mapping_file <- opt$output
cat("Writing cluster to kmer mapping to ", cluster_to_kmer_mapping_file, "\n")
write_tsv(cluster_to_kmer_mapping, cluster_to_kmer_mapping_file)

