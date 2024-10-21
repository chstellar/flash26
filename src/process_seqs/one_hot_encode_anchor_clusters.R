## one_hot_encode_anchor_clusters.R
# Daniel Cotter
# takes in the processed sample sequences and writes out a one hot encoded matrix
# that can be passed into glmnet for training and testing

## Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(mltools))
suppressPackageStartupMessages(library(feather))


## Parse arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to the sample sequences file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to the output feature matrix"),
  make_option(c("-k", "--kmer_width"), type="integer", default=54, 
              help="Length of the kmer to one hot encode (default=3)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Please provide a path to the sample sequences file and output file")
}

# read in the sample sequences file
sample_sequences <- Biostrings::readDNAStringSet(opt$input)
df <- data.frame(sample_name=names(sample_sequences), sequence=sample_sequences, row.names =NULL)

# one hot encode the sample sequences every k bases
kmer <- opt$kmer_width

# determine the widths vector to separate the sequences into
widths = rep(kmer, ceiling(nchar(df[1,]$sequence)/kmer))
cluster_names = paste0("cluster_", 1:ceiling(nchar(df[1,]$sequence)/kmer))
names(widths) = cluster_names

# split the sequences into kmers every 54 positions
df_wide <- df %>% 
  separate_wider_position(sequence, widths = widths) %>%
  mutate(across(starts_with("cluster"), \(x) as.factor(x)))

# one hot encode each column (adding a column for each unique seqeunce in the cluster column)
df_ohe <- mltools::one_hot(as.data.table(df_wide%>%select(-sample_name)))
df_ohe <- cbind(df_wide %>% select(sample_name), df_ohe)  

# write out the one hot encoded matrix as a feather dataframe
feather::write_feather(df_ohe, opt$output)