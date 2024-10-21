# cluster_with_mmseq.R
# Daniel Cotter
# takes in a path to a list of anchor sequences and uses mmseqs to cluster them
# converts the output to a fasta file with the record name the same as the sequence
# then uses the output tsv file to format the cluster file 
# (2 columns, 1 = cluster_id and 2 = sequence)

## Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))


## Parse arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to the anchors file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to the output file"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character", default=NULL),
  make_option(c("--mmseqs"), "Path to the mmseqs executable", 
              type="character", default="/oak/stanford/groups/horence/dcotter1/software/mmseqs/bin/mmseqs")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Please provide a path to the anchors file, output file")
}

# create a temporary directory to store intermediate files
if (!is.null(opt$temp_dir)) {
  temp_dir <- ifelse(grepl("/$", opt$temp_dir),
                     opt$temp_dir, 
                     paste0(opt$temp_dir, "/"))
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output_prefix), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

mmseqs_temp <- file.path(temp_dir, "mmseqs_temp")
system(paste("mkdir -p", mmseqs_temp))

# read in the anchors file and write it out as a temporary fasta file
anchors <- fread(opt$input, header=FALSE, sep="\t", col.names=c("sequence"))
anchors_dna <- DNAStringSet(anchors$sequence)
names(anchors_dna) <- anchors$sequence
writeXStringSet(anchors_dna, file=file.path(temp_dir, "all_anchors.fasta"))

# run mmseqs to cluster the anchors
mmseqs_cmd <- paste(opt$mmseqs, "easy-cluster", 
                    file.path(temp_dir, "all_anchors.fasta"), 
                    paste(mmseqs_temp, "cluster", sep="_"),
                    mmseqs_temp)

system(mmseqs_cmd)

# read in the mmseqs output and write it out as a temporary tsv file
mmseqs_out <- file.path(temp_dir, "mmseqs_temp_cluster_cluster.tsv")
clusters <- fread(mmseqs_out, header=F, col.names = c("rep", "seq"))

# group by the representatives and add a unique group id column for each cluster
clusters <- clusters %>% group_by(rep) %>% 
  mutate(cluster_id=cur_group_id()) %>% 
  mutate(cluster_id = cluster_id - 1) %>%
  arrange(cluster_id) %>% ungroup()

# write out the clusters as a tsv file
clusters %>% select(cluster_id, seq) %>% write_tsv(opt$output, col_names=FALSE, quote="none")
