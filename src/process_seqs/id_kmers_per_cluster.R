# id_kmers_per_cluster.R
# Daniel Cotter
# 2024-09-18

# This script takes in one embeddings file and the ordering file and outputs a new
# embeddings tsv with samples as rows and the top 10 embeddings by variance per cluster as columns.


## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(resample))
suppressPackageStartupMessages(library(furrr))


## parse arguments --------
# define command line arguments
# define command line arguments
option_list <- list(
  make_option(c("-o", "--ordering"), help="Ordering file", type="character"),
  make_option(c("-p", "--output"), help="Output file.", type="character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$ordering) | is.null(opt$output)) {
  stop("Must provide embeddings, ordering, and output prefix")
}

## print a summary of the arguments
cat("\n####################\n")
cat("Running grab_top_embeddings_by_variance.R with the following arguments:\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Output file: ", opt$output, "\n")
cat("####################\n\n")

# load the ordering file
cat("Loading the ordering file...\n")
ordering <- fread(ordering_temp, header=F, sep="\t", 
                  col.names = c("sample_name", "seq", "kmer", "start", "end")) 
ordering <- ordering %>% select(sample_name, kmer)

# labeling clusters in the ordering file (they appear one per row per sample)
ordering <- ordering %>% group_by(sample_name) %>% mutate(cluster=row_number()) %>%
  mutate(cluster=cluster-1) %>% ungroup()

cluster_to_kmer_mapping <- ordering %>% select(cluster, kmer) %>% distinct()
cluster_to_kmer_mapping <- cluster_to_kmer_mapping %>% 
  arrange(cluster, kmer)

# write out the cluster to kmer mapping
cluster_to_kmer_mapping_file <- opt$output
cat("Writing cluster to kmer mapping to ", cluster_to_kmer_mapping_file, "\n")
write_tsv(cluster_to_kmer_mapping, cluster_to_kmer_mapping_file

