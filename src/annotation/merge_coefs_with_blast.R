# Load necessary libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
#library(Biostrings)

# Define command line options
option_list <- list(
  make_option(c("-a", "--blast_annotations"), type = "character", default = NULL, 
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
if (is.null(opt$blast_annotations) || is.null(opt$coefficients) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All arguments must be supplied", call. = FALSE)
}

# Read in the data
annotations <- fread(opt$blast_annotations, header = TRUE)

if (str_detect(opt$blast_annotations, "blastp")) {
  annotations <- annotations %>% select(query, evalue, identity, qcovs, qframe, stitle, `NCBI_protein_accession`, UniProt_accession, method, GO) %>% 
    mutate(cluster=str_extract(query, "(cluster_\\d+)_", group = 1))
} else {
  annotations <- annotations %>% select(query, evalue, identity, qcovs, features, contains("window")) %>% 
    mutate(cluster=str_extract(query, "(cluster_\\d+)_", group = 1))
}



coefficients <- fread(opt$coefficients, header = TRUE)
coefficients <- coefficients %>% mutate(cluster = str_extract(feature, "(cluster_\\d+)_", group = 1))

merged_data <- coefficients %>% full_join(annotations, by ="cluster", relationship="many-to-many")

merged_data <- merged_data %>% rowwise() %>% mutate(max_coef=max(abs(as.numeric(
  strsplit(gsub("\\[|\\]", "", coefficients), split=",", fixed=TRUE)[[1]])))) %>% 
  arrange(metadata_category, -abs(max_coef), -identity) %>% select(-max_coef)
# Merge the data on the 'cluster' column 
# there may be multiple entries in annotations that match
# in this case, the coefficients will be repeated for each match

# Write the merged data to a new CSV file
write_tsv(merged_data, opt$output, col_names = T, quote="needed")