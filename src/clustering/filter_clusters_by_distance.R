# filter_clusters_by_distance.R
# Daniel Cotter
# 2024-09-12

# This script takes in the anchor clusters file and a SPLASH stats file and processes
# the clusters by picking one anchor per cluster based on the mean effect size and number of nonzero samples
# and then using a distance metric to filter out anchors that are too far away from the seed anchor.

## import packages --------
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(furrr))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input_anchor_clusters"), "Input anchor clusters", type="character"),
  make_option(c("-s", "--splash_stats"), "SPLASH stats file", type="character"),
  make_option(c("--effect_size_cutoff"), "Effect size cutoff for filtering anchors from the SPLASH stats file",
              type="double", default=0.5),
  make_option(c("-o", "--output"), "Output file", type="character"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character"),
  make_option(c("--distance_metric"), "Distance metric to use for filtering", 
              type="character", default="lev"),
  make_option(c("--num_cores"), "Number of cores to use for parallel processing"
              , type="integer", default=1),
  make_option(c("--distance_threshold"), "Distance threshold for filtering anchors",
              type="integer", default=5),
  make_option(c("--max_clusters_to_process"), "Maximum number of clusters to process",
              type="integer", default=100000)
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input_anchor_clusters) | !file.exists(opt$splash_stats) | is.null(opt$output)) {
  stop("Must provide input and output files")
}

# set up parallel processing
# multisession works much faster than multicore for this step
plan(multisession, workers = opt$num_cores)

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

# cat to screen the input files, output files, temp dir and paramaters
cat("\n###################################################################\n")
cat("Running reorder_anchor_clusters.R with the following parameters:\n")
cat(paste("Input clusters:", opt$input_anchor_clusters))
cat("\n")
cat(paste("SPLASH stats file:", opt$splash_stats))
cat("\n")
cat(paste("Output file:", opt$output))
cat("\n")
cat(paste("Temporary directory:", temp_dir))
cat("\n")
cat(paste("Distance metric:", opt$distance_metric))
cat("\n")
cat("###################################################################\n\n")

# read in the anchor cluster file
cat("Reading in the anchor clusters file\n")
anchor_clusters <- fread(opt$input_anchor_clusters,
                         header = FALSE,
                         col.names = c("cluster_id", "anchor"))


# read in the splash stats file (grepping for the anchors in the anchor file)
cat("Reading in the splash stats file\n")
## load the data
# read in the headers of the input file to identify the effect_size_bin column
headers <- fread(opt$splash_stats, nrows = 1, header = TRUE)
effect_size_bin_col <- grep("effect_size_bin", names(headers))
nonzero_samples_col <- grep("number_nonzero_samples", names(headers))
max_col_to_read <- max(effect_size_bin_col, nonzero_samples_col) + 1

# load the input file using awk to filter out rows with effect size < 0.7
EFFECT_SIZE_CUTOFF = opt$effect_size_cutoff
if (grepl(".gz$", opt$splash_stats)) {
  load_cmd <- paste0("cut -f2 ", opt$input_anchor_clusters, " | zgrep -Ff - ",
                     opt$splash_stats, " | ",
                     "cut -f1-", max_col_to_read, " ",
                     " | awk '{OFS=\"\t\"}{if ($",
                     effect_size_bin_col, " >= ", EFFECT_SIZE_CUTOFF,
                     ") print $0}'")
} else {
  load_cmd <- paste0("cut -f2 ", opt$input_anchor_clusters, " | grep -Ff - ",
                     opt$splash_stats, " | ",
                     "cut -f1-", max_col_to_read, " ",
                     " | awk '{OFS=\"\t\"}{if ($",
                     effect_size_bin_col, " >= ", EFFECT_SIZE_CUTOFF,
                     ") print $0}'")
}

splash_stats <- fread(cmd = load_cmd, header = FALSE,
                      col.names = names(headers)[1:max_col_to_read])
splash_stats <- splash_stats %>% filter(anchor %in% anchor_clusters$anchor)

# join the splash stats file with the anchor clusters file
cat("Joining the splash stats file with the anchor clusters file\n")
anchor_clusters_with_stats <- anchor_clusters %>%
  left_join(splash_stats, by = "anchor") %>% select(cluster_id, anchor, effect_size_bin, number_nonzero_samples) %>%
  mutate(effect_size_bin = as.numeric(effect_size_bin),
         number_nonzero_samples = as.numeric(number_nonzero_samples))

# reorder the cluster ids based on the mean effect size and number of nonzero samples
cat("Reordering the cluster ids based on the mean effect size and number of nonzero samples\n")
anchor_clusters_with_stats <- anchor_clusters_with_stats %>%
  group_by(cluster_id) %>%
  mutate(mean_effect_size = mean(effect_size_bin),
         mean_nonzero_samples = mean(number_nonzero_samples)) %>%
  mutate(sort_val = mean_effect_size * mean_nonzero_samples) %>%
  arrange(desc(sort_val), cluster_id)  %>%
  ungroup() %>% 
  mutate(new_cluster_id = as.numeric(factor(cluster_id, levels = unique(cluster_id)))) %>% 
  ungroup() %>% group_by(new_cluster_id) %>%
  arrange(new_cluster_id, desc(effect_size_bin)) %>% 
  ungroup() %>% select(new_cluster_id, anchor)

# if we have --distance_metric noCluster, write out the
# anchor clusters file with the new cluster ids and exit
if (opt$distance_metric == "noCluster") {
  cat("Writing out the anchor clusters file with the new cluster ids\n")
  write_tsv(anchor_clusters_with_stats, opt$output, col_names = F, quote="none")
  cat("Done!\n")
  quit(status = 0)
}

# keep only one anchor per cluster (the one with the highest effect size)
one_anchor_per_cluster <- anchor_clusters_with_stats %>%
  group_by(new_cluster_id) %>%
  dplyr::slice(1) %>%
  ungroup()

# now that we have the new clusters and their representative (seed) anchors, 
# we can filter out the anchors in each cluster that are too far away from the seed anchor
representative_anchors <- one_anchor_per_cluster$anchor

# separate the anchor cluster df into a list of dfs, one per cluster
all_anchor_sets <- anchor_clusters_with_stats %>%
  group_split(new_cluster_id, .keep = T)


# Define a function to return only the anchors that are within a certain distance of the representative anchor
# 
# Args:
#   anchor_set: A data frame containing the set of anchors to be filtered.
#   representative_anchor: A string representing the anchor to compare against.
#   distance_metric: A string specifying the distance metric to use ("lev" for Levenshtein distance, "ham" for Hamming distance).
#   distance_threshold: An integer specifying the maximum allowable distance (default is 5).
#
# Returns:
#   A filtered data frame containing only the anchors within the specified distance of the representative anchor.
#
# Raises:
#   Error if the distance metric is not recognized.
#
distance_filter <- function(anchor_set, representative_anchor, distance_metric, distance_threshold=5) {
  # Check if the distance metric is Levenshtein distance
  if (distance_metric == "lev") {
    # Calculate the Levenshtein distance between each anchor and the representative anchor
    lev_dist <- stringdist::stringdistmatrix(anchor_set$anchor, representative_anchor, method = "lv")
    # Create a logical vector indicating which anchors are within the distance threshold
    dist_filter <- as.logical(lev_dist <= distance_threshold)
    # Filter the anchor set based on the distance threshold
    anchor_set <- anchor_set[dist_filter,]
    # Return the filtered anchor set
    return(anchor_set)
  # Check if the distance metric is Hamming distance
  } else if (distance_metric == "ham") {
    # Calculate the Hamming distance between each anchor and the representative anchor
    ham_dist <- stringdist::stringdistmatrix(anchor_set$anchor, representative_anchor, method = "hamming")
    # Create a logical vector indicating which anchors are within the distance threshold
    dist_filter <- as.logical(ham_dist <= distance_threshold)
    # Filter the anchor set based on the distance threshold
    anchor_set <- anchor_set[dist_filter,]
    # Return the filtered anchor set
    return(anchor_set)
  # Raise an error if the distance metric is not recognized
  } else {
    stop("Distance metric not recognized")
  }
}


# apply the distance filter to each cluster
cat("Applying the distance filter to each cluster\n")
cat("Using the ", opt$distance_metric, " distance metric\n")

# only apply the distance filter to the first 100000 clusters
# this is to speed up the process
max_anchor_clusters <- opt$max_clusters_to_process
max_anchor_clusters <- min(max_anchor_clusters, length(all_anchor_sets)) # make sure we don't go over the number of clusters
all_anchor_sets <- all_anchor_sets[1:max_anchor_clusters]
representative_anchors <- representative_anchors[1:max_anchor_clusters]


filtered_anchor_sets <- future_map2(all_anchor_sets, representative_anchors, 
                                    \(x, y) distance_filter(x, y, distance_metric=opt$distance_metric, 
                                                            distance_threshold=opt$distance_threshold))

# combine the filtered anchor sets into a single df
filtered_anchor_clusters <- bind_rows(filtered_anchor_sets) %>% 
  arrange(new_cluster_id)

# write out the new anchor clusters file
cat("Writing out the new anchor clusters file\n")
write_tsv(filtered_anchor_clusters, opt$output, col_names = F, quote="none")

cat("Done!\n")


