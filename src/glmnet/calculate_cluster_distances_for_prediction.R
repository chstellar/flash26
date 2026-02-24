# calculate_cluster_distances_for_prediction.R
# Daniel Cotter
# 2026-06-15

# This script takes in one embeddings file and the ordering file and processes
# each cluster one at a time to calculate distances between samples for
# downstream prediction tasks. With the --interactions flag, it will also
# calculate pairwise interaction features between clusters as the angle between
# the pair of embedding vectors per sample.
# The output will match that of grab_top_embeddings_by_variance.R but here
# the embeddings are collapsed to distances rather than kept as raw embeddings.

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
  make_option(c("-e", "--embeddings"), 
              help = "Embeddings file", type = "character"),
  make_option(c("-o", "--ordering"),
              help = "Ordering file", type = "character"),
  make_option(c("-p", "--output"), help = "Output file.", type = "character"),
  make_option(c("--temp_dir"), 
              help = "Temporary directory to store intermediate files", 
              type = "character"),
  make_option(
    c("--num_threads"),
    help = "Number of threads to use for parallel operations",
    type = "integer",
    default = 1
  ),
  make_option(
    c("--normalized"),
    help = "Whether to normalize the distances after calculating distances.",
    action = "store_true",
    default = FALSE
  ),
  make_option(
    c("--interactions"),
    help = "Whether to calculate pairwise interaction features between clusters.",
    action = "store_true",
    default = FALSE
  )
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$embeddings) |
    !file.exists(opt$ordering) | is.null(opt$output)) {
  stop("Must provide embeddings, ordering, and output prefix")
}

# create a temporary directory to store intermediate files
if (!is.null(opt$temp_dir)) {
  temp_dir <- ifelse(grepl("/$", opt$temp_dir),
                     opt$temp_dir,
                     paste0(opt$temp_dir, "/"))
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

# define future plans
setDTthreads(opt$num_threads)
plan(multicore, workers = opt$num_threads)
options(future.globals.maxSize = 8000 * 1024^2)

## print a summary of the arguments
cat("\n####################\n")
cat("Running grab_top_embeddings_by_variance.R with the following arguments:\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Embeddings file: ", opt$embeddings, "\n")
cat("Output file: ", opt$output, "\n")
cat("Temporary directory: ", temp_dir, "\n")
cat("####################\n\n")

## load data --------
# load embeddings
cat("\nLoading embeddings...\n")
# copy the embeddings file to the temp directory to speed up I/O
embeddings_temp <- file.path(temp_dir, "raw_embeddings_temp.tsv")
system(paste("cp", opt$embeddings, embeddings_temp))
embeddings <- fread(embeddings_temp, header = F)
colnames(embeddings) <- c("kmer", paste0("embedding_", 1:(ncol(embeddings) - 1)))

# load the ordering file
cat("Loading the ordering file...\n")
# copy the ordering file to the temp directory to speed up I/O
ordering_temp <- file.path(temp_dir, "ordering_temp.tsv")
system(paste("cp", opt$ordering, ordering_temp))
ordering <- fread(
  ordering_temp,
  header = F,
  sep = "\t",
  col.names = c("sample_name", "seq", "kmer", "start", "end")
)

# process ordering file and assign clusters (since the file is already ordered
# by cluster, we can just assign cluster IDs based on row number per sample)
ordering <- ordering %>% select(sample_name, seq, kmer)

ordering <- ordering %>%
  group_by(sample_name) %>%
  mutate(cluster = row_number()) %>%
  mutate(cluster = cluster - 1) %>%
  ungroup()

# arrange in cluster and kmer order
ordering <- ordering %>% arrange(cluster, kmer)

# split into a list of clusters to post process
clusters <- ordering %>%
  group_by(cluster) %>%
  group_split()

# takes in a cluster data frame and the embeddings data frame and write a new
# file with the embeddings for that cluster joined on 
join_and_write_clusters <- function(cluster_df, filename, all_embeddings) {
  cluster_df %>%
    select(-cluster) %>%
    left_join(all_embeddings, by = "kmer") %>%
    select(-kmer) %>%
    fwrite(file = filename,
           nThread = 1,
           col.names = T)
}

# create a temporary directory to store per-cluster embeddings
temp_embeddings_dir <- file.path(temp_dir, "embeddings_per_cluster/")
system(paste0("rm -r ", temp_embeddings_dir))
system(paste("mkdir -p", temp_embeddings_dir))
cluster_files <- paste0(temp_embeddings_dir,
                        "embeddings_cluster_",
                        0:(length(clusters) - 1),
                        ".csv")

# Cleaning up temp directory
cat("Writing all clusters and their embeddings out to file in: ",
    temp_embeddings_dir)

# this condition will not recalculate if the files already exist, however,
# we delete the temp directory above and so this will always run
if (!sum(file.exists(cluster_files)) == length(cluster_files)) {
  future_walk2(
    clusters,
    cluster_files,
    \(x, y) join_and_write_clusters(x, y, all_embeddings = embeddings)
  )
}

cat("Formatting the embeddings for downstream use...\n")

# For each cluster we have a matrix of samples x embedding dimensions (and a 
# column for sequences and sample name). We want to select which sequence is the most 
# abundant per sample and then calculate the distance between each sample's embedding
# vector and the embedding vector of the most abundant sequence.
calculate_distance_to_most_abundant_embedding <- function(in_file, normalized=FALSE) {
  # grab the cluster number from the file name
  cluster_num <- str_extract(in_file, "cluster_(\\d+).csv", group = 1) %>% as.integer()
  # read in the file and keep the header
  temp_dt <- fread(in_file, header = T, nThread = 1)
  temp_dt <- temp_dt %>% arrange(sample_name)
  most_abundant_seq <- temp_dt %>% select(seq) %>%
    group_by(seq) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    slice(1) %>% pull(seq)
  
  # grab the embeddings for the most abundant sequence
  most_abundant_seq_embeddings <- temp_dt %>% filter(seq== most_abundant_seq) %>%
    slice(1) %>% 
    select(starts_with("embedding")) %>%
    as.matrix() %>% as.numeric()
  
  # grab the matrix of all embeddings for all samples
  all_embeddings_matrix <- temp_dt %>%
    select(starts_with("embedding")) %>%
    as.matrix()
  
  # calculate l2 norm between most_abundant_seq_embeddings and each row of 
  # all_embeddings_matrix. Use matrix operations for speed
  diffs <- all_embeddings_matrix - matrix(
    most_abundant_seq_embeddings,
    nrow = nrow(all_embeddings_matrix),
    ncol = ncol(all_embeddings_matrix),
    byrow = TRUE
  )
  
  distances <- sqrt(rowSums(diffs^2))
  
  # if the user wants normalized distances, scale and center them
  if (normalized) {
    distances <- scale(distances)
  }
  
  # create a data frame with sample names and distances
  distance_dt <- data.frame(
    sample_name = temp_dt$sample_name,
    !!paste0("cluster_", cluster_num, "_relative_distance") := distances
  )
}


calculate_distance_to_average_embedding <- function(in_file, normalized=FALSE) {
  # grab the cluster number from the file name
  cluster_num <- str_extract(in_file, "cluster_(\\d+).csv", group = 1) %>% as.integer()
  # read in the file and keep the header
  temp_dt <- fread(in_file, header = T, nThread = 1)
  temp_dt <- temp_dt %>% arrange(sample_name)
  
  # calculate the average embedding 
}



grab_top_variance_columns <- function(in_file, num_cols, normalized) {
  cluster_num <- str_extract(in_file, "cluster_(\\d+).csv", group = 1) %>% as.integer()
  temp_dt <- fread(in_file, header = T, nThread = 1) %>%
    select(sample_name, starts_with("embedding")) # first filter for only one cluster
  colnames(temp_dt) <- ifelse(
    grepl("embedding", colnames(temp_dt)),
    yes = paste0("cluster_", cluster_num, "_", colnames(temp_dt)),
    no = colnames(temp_dt)
  )
  if (normalized) {
    temp_dt <- temp_dt %>% mutate(across(starts_with("cluster"), scale))
  }
  top_var_cols <- resample::colVars(temp_dt %>% select(starts_with("cluster"))) %>%
    enframe() %>%
    slice_max(n = num_cols, value, with_ties = FALSE) %>%
    pull(name)
  temp_dt <- temp_dt %>%
    select(sample_name, all_of(top_var_cols)) %>%
    arrange(sample_name)
  if (cluster_num != 0) {
    temp_dt <- temp_dt %>% select(-sample_name)
  }
  return(temp_dt)
}

cat("Calculating top variance components per cluster...\n")
top_var_dt <- future_map_dfc(
  cluster_files,
  \(x) grab_top_variance_columns(
    x,
    num_cols = opt$num_to_keep,
    normalized = opt$normalized_embeddings
  )
)
top_var_dt <- top_var_dt %>% relocate(sample_name)

# write out the embeddings matrix to a temp feather file
embeddings_feather <- opt$output
cat("Writing top variance embeddings to ", embeddings_feather, "\n")
feather::write_feather(top_var_dt, embeddings_feather)

# Cleaning up temp directory
system(paste0("rm -r", temp_embeddings_dir))
