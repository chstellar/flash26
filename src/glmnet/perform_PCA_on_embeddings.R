# perform_PCA_on_embeddings.R
# Daniel Cotter
# 2026-06-15

# This script takes in one embedding file and the ordering file and processes
# each cluster one at a time to calculate distances between samples for
# downstream prediction tasks.
# The output will match that of grab_top_embeddings_by_variance.R but here
# the embeddings are collapsed to the top N PCs per cluster rather than
# kept as raw embeddings.

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(resample))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(RhpcBLASctl)) # for controlling BLAS threads with furrr


## parse arguments --------
# define command line arguments
# define command line arguments
option_list <- list(
  make_option(c("-e", "--embeddings"), help = "Embeddings file", type = "character"),
  make_option(c("-o", "--ordering"), help = "Ordering file", type = "character"),
  make_option(c("-p", "--output"), help = "Output file.", type = "character"),
  make_option(c("--temp_dir"), help = "Temporary directory to store intermediate files", type = "character"),
  make_option(
    c("--num_threads"),
    help = "Number of threads to use for parallel operations",
    type = "integer",
    default = 1
  ),
  make_option(
    c("--normalized"),
    help = "Whether to scale and center embeddings before calculating PCs.",
    action = "store_true",
    default = FALSE
  ),
  make_option(
    c("--num_pcs"),
    help = "Number of principal components to retain per cluster.",
    type = "integer",
    default = 5
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
    paste0(opt$temp_dir, "/")
  )
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

# define future plans
# determine number of futures and number of threads per future to use
num_threads <- opt$num_threads
max_futures <- num_threads %/% 4 - 1
if (max_futures < 1) {
  max_futures <- 1
}
threads_per_future <- ceiling(num_threads / max_futures)

# set the future plan
plan(multisession, workers = max_futures) # multisession is better for many small tasks
options(future.globals.maxSize = 8000 * 1024^2)
RhpcBLASctl::blas_set_num_threads(threads_per_future) # to avoid oversubscribing threads

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
    fwrite(
      file = filename,
      nThread = threads_per_future,
      col.names = T
    )
}

# create a temporary directory to store per-cluster embeddings
temp_embeddings_dir <- file.path(temp_dir, "embeddings_per_cluster/")
system(paste0("rm -r ", temp_embeddings_dir))
system(paste("mkdir -p", temp_embeddings_dir))
cluster_files <- paste0(
  temp_embeddings_dir,
  "embeddings_cluster_",
  0:(length(clusters) - 1),
  ".csv"
)

# Cleaning up temp directory
cat(
  "Writing all clusters and their embeddings out to file in: ",
  temp_embeddings_dir
)

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
calculate_PCA_by_cluster <- function(in_file,
                                     normalized = FALSE,
                                     num_pcs = 5) {
  # grab the cluster number from the file name
  cluster_num <- str_extract(in_file, "cluster_(\\d+).csv", group = 1) %>% as.integer()
  # read in the file and keep the header
  temp_dt <- fread(in_file, header = T, nThread = threads_per_future)
  temp_dt <- temp_dt %>% arrange(sample_name)

  # get a vector that will be an indicator variable for missing sequences
  temp_dt <- temp_dt %>%
    mutate(is_missing = grepl(pattern = "NNNN", x = temp_dt$seq))

  # filter the embedding data table for only samples with non-missing data
  temp_dt_filtered <- temp_dt %>% filter(!is_missing)

  embedding_matrix <- temp_dt_filtered %>%
    column_to_rownames("sample_name") %>%
    select(starts_with("embedding")) %>%
    as.matrix()

  # normalize if specified
  if (normalized) {
    embedding_matrix <- scale(embedding_matrix, center = TRUE, scale = TRUE)
  }

  # calculate PCA
  pca_res <- prcomp(embedding_matrix, center = FALSE, scale. = FALSE)

  # grab top N PCs as a data frame with sample_name as a column
  pcs_dt <- as.data.frame(pca_res$x[, 1:num_pcs]) %>%
    rownames_to_column("sample_name")

  # rename the PC columns to indicate cluster
  colnames(pcs_dt) <- ifelse(
    grepl("PC", colnames(pcs_dt)),
    yes = paste0("cluster_", cluster_num, "_PC_", 1:num_pcs),
    no = colnames(pcs_dt)
  )

  # now we need to reinsert rows for missing samples with NA values
  if (any(temp_dt$is_missing)) {
    missing_samples <- temp_dt %>% filter(is_missing) %>%
      # create a column for each PC named the same as above
      select(sample_name, is_missing)
    PC_cols <- paste0("cluster_", cluster_num, "_PC_", 1:num_pcs)
    for (col in PC_cols) {
      missing_samples[[col]] <- 0
    }
    missing_samples <- missing_samples %>%
      select(sample_name, all_of(PC_cols), is_missing)
    pcs_dt <- bind_rows(pcs_dt %>%
      mutate(is_missing = FALSE), missing_samples) %>%
      arrange(sample_name) %>%
      mutate(is_missing = ifelse(is_missing, 1, 0))
  } else {
    pcs_dt <- pcs_dt %>%
      mutate(is_missing = 0)
  }

  pcs_dt <- pcs_dt %>%
    relocate(sample_name, is_missing) %>%
    dplyr::rename(!!paste0("cluster_", cluster_num, "_is_missing") := is_missing)

  # return the pcs data table
  pcs_dt <- pcs_dt %>%
    arrange(sample_name)

  if (cluster_num != 0) {
    pcs_dt <- pcs_dt %>%
      select(-sample_name)
  }

  return(pcs_dt)
}

cat("Calculating top PCs per cluster...\n")
all_pc_dt <- future_map(cluster_files,
  function(x) {
    calculate_PCA_by_cluster(
      x,
      normalized = opt$normalized,
      num_pcs = opt$num_pcs
    )
  },
  .progress = T
)

# bind all the pc data tables together and relocate sample_name to front
all_pc_dt <- bind_cols(all_pc_dt) %>% relocate(sample_name)

# write out the pc matrix to a temp feather file
pc_feather <- opt$output
cat("Writing top variance embeddings to ", pc_feather, "\n")
feather::write_feather(all_pc_dt, pc_feather)

# Cleaning up temp directory
system(paste0("rm -r", temp_embeddings_dir))
