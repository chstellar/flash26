# Load necessary libraries
library(optparse)
library(data.table)
library(tidyverse)
#library(Biostrings)

# Define command line options
option_list <- list(
    make_option(c("-a", "--annotations"), type = "character", default = NULL, 
                            help = "Path to the annotations CSV file", metavar = "character"),
    make_option(c("-c", "--coefficients"), type = "character", default = NULL, 
                            help = "Path to the coefficients CSV file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL, 
                            help = "Path to the output CSV file", metavar = "character")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$annotations) || is.null(opt$coefficients) || is.null(opt$output)) {
    print_help(opt_parser)
    stop("All arguments must be supplied", call. = FALSE)
}

# Read in the data
annotations <- fread(opt$annotations, header = TRUE)

flat_annotations <- annotations %>% 
  mutate(summary = paste0("[[", kmer, ": ", seq, " : ", str_trunc(query,width=200,side="right",), " : ", stats, " ]]")) %>% 
  group_by(cluster) %>% summarise(anno = str_c(summary, sep=";", collapse=";"))
coefficients <- fread(opt$coefficients, header = TRUE)
coefficients <- coefficients %>% mutate(cluster = str_extract(feature, "(cluster_\\d+)_", group = 1))

# Merge the data on the 'cluster' column 
# there may be multiple entries in annotations that match
# in this case, the coefficients will be repeated for each match
merged_data <- left_join(coefficients, flat_annotations, by = "cluster")

# Write the merged data to a new CSV file
write_tsv(merged_data, opt$output, col_names = T, quote="needed")

# Print a message indicating the script has finished
cat("Merging complete. Output saved to", opt$output, "\n")
if (nrow(merged_data %>% filter(!is.na(anno))) == 0) {
  out_file <- gsub("nonzero_coefficients_annotated.tsv", "significant_sequences.fasta", opt$output)
  system(paste("touch", out_file))
} else {
  # get fasta of all nonzero_seqeunces 
  sig_seqs <- merged_data %>% filter(!grepl("rcept", feature)) %>% 
    rowwise() %>% 
    mutate(max_coef=max(abs(as.numeric(
      strsplit(gsub("\\[|\\]", "", coefficients), split=",", fixed=TRUE)[[1]])))) %>% 
    group_by(metadata_category) %>% 
    slice_max(max_coef, n=10) %>% ungroup() %>% select(cluster, anno) %>% 
    distinct(cluster, .keep_all = T) %>% 
    mutate(seqs = sapply(str_extract_all(anno, "[ACTGN]{30,70}"), \(x) toString(unique(x)))) %>% 
    select(-anno) %>% 
    separate_longer_delim(seqs, ", ") %>% 
    filter(!grepl("NNNNNNNNN", seqs)) %>%
    filter(!is.na(cluster))
  
  sig_seqs <- sig_seqs %>% mutate(cluster=paste(cluster, seqs, sep="_"))
  
  library(Biostrings)
  
  dna <- DNAStringSet(sig_seqs$seqs)
  names(dna) <- sig_seqs$cluster
  writeXStringSet(dna, gsub("nonzero_coefficients_annotated.tsv", "significant_sequences.fasta", opt$output))
}
