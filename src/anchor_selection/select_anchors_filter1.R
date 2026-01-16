# select_anchors_filter1.R
# Daniel Cotter
# 2024-09-12

# This script selects the most important anchors from the output of SPLASH
# it prioritizes anchors that have a high effect size and also have a large
# number of nonzero samples. It also filters out anchors that have lookup
# table hits to artifcats. The output is a list of anchor sequences that can
# be used in downstream analyses.

# Filter 1: effect size >= 0.6
# Filter 1: number of nonzero samples > 10th percentile
# Filter 1: no lookup table hits to artifacts
# Filter 1: select top 1,000,000 anchors by number of nonzero samples

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input"), "Input file", type = "character"),
  make_option(c("-o", "--output"), "Output file", type ="character"),
  make_option(c("-n", "--num_anchors"), "Number of anchors to select",
              type = "integer", default = 1000000),
  make_option(c("-e", "--effect_size"), "Effect size threshold",
              type = "numeric", default = 0.6),
  make_option(c("-l", "--lookup_table"), "Lookup table file",
              type = "character"),
  make_option(c("--splash_bin"), "Path to SPLASH binary folder",
              type = "character", default="splash-2.11.6/"),
  make_option(c("--temp_dir"),
              "Temporary directory to store intermediate files",
              type = "character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input) || is.null(opt$output)) {
  stop("Must provide input and output files")
}
if (!file.exists(opt$lookup_table)) {
  stop("Must provide path to lookup table")
}

# define output files
anchors_only_out <- opt$output

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
cat("Running select_important_anchors.R with the following parameters:\n")
cat(paste("Input file:", opt$input))
cat("\n")
cat(paste("Output anchors:", anchors_only_out))
cat("\n")
cat(paste("Number of anchors to select:", opt$num_anchors))
cat("\n")
cat(paste("Effect size threshold:", opt$effect_size))
cat("\n")
cat(paste("Lookup table file:", opt$lookup_table))
cat("\n")
cat(paste("SPLASH binary folder:", opt$splash_bin))
cat("\n")
cat(paste("Temporary directory:", temp_dir))
cat("\n")
cat("###################################################################\n\n")


## load the data
# read in the headers of the input file to identify the effect_size_bin column
## load the data
# read in the headers of the input file to identify the effect_size_bin column
headers <- fread(opt$input, nrows = 2, header = TRUE)

if (nchar(headers$anchor[1]) < 12) {
  opt$effect_size <- 0.1 # change if using short anchors
}

effect_size_bin_col <- grep("effect_size_bin", names(headers))
nonzero_samples_col <- grep("number_nonzero_samples", names(headers))
max_col_to_read <- max(effect_size_bin_col, nonzero_samples_col) + 1

# load the input file using awk to filter out rows with effect size < 0.7
EFFECT_SIZE_CUTOFF = opt$effect_size

load_cmd <- paste0("cut -f1-", max_col_to_read, " ",
                   opt$input, " | awk '{OFS=\"\t\"}{if ($",
                   effect_size_bin_col, " >= ", EFFECT_SIZE_CUTOFF,
                   ") print $0}'")
if (grepl(".gz$", opt$input)) {
  load_cmd <- paste0("z", load_cmd)
}

# only select the columns up to number_nonzero_samples
cat(paste("Reading in anchors from",
          opt$input, "with effect size >=",
          opt$effect_size))
cat("\n\n")
dt <- fread(cmd = load_cmd, col.names = names(headers)[1:max_col_to_read]) %>%
  mutate(effect_size_bin = as.numeric(effect_size_bin),
         number_nonzero_samples = as.numeric(number_nonzero_samples))

# filter out the bottom 10% of anchors by number nonzero samples
sample_cutoff <- quantile(dt$number_nonzero_samples, 0.1)
dt <- dt %>% filter(number_nonzero_samples >= sample_cutoff)

## run lookup table to filter out artifacts -----------
# write all anchors to a file
anchors_to_keep <- dt %>% select(anchor) %>% distinct() 
anchor_file <- file.path(temp_dir, "all_anchors.txt") %>% gsub("//", "/", .)
anchors_to_keep %>%
  write.table(anchor_file, row.names=FALSE, col.names=FALSE, quote=FALSE)

# write anchor to a fasta file as well (for lookuptable)
anchors_fasta_file <- file.path(temp_dir, "all_anchors.fasta") %>%
  gsub("//", "/", .)
anchors_dna <- DNAStringSet(anchors_to_keep$anchor)
names(anchors_dna) <- anchors_to_keep$anchor
writeXStringSet(anchors_dna, anchors_fasta_file)

# run lookup table
out_lookup_stats <- file.path(temp_dir, "lookup_stats.txt") %>%
  gsub("//", "/", .)
lookup_cmd <- paste0(file.path(opt$splash_bin, "lookup_table"),
                     " query --kmer_skip 1 ",
                     "--truncate_paths --stats_fmt with_stats ",
                     anchors_fasta_file,
                     " ", opt$lookup_table, " ", out_lookup_stats)
system(paste("rm", out_lookup_stats))
if (file.exists(out_lookup_stats)) {
  cat("Lookup stats already in tmp directory. Delete tmp to force run again.\n")
  cat("Filtering anchors for artifacts...\n")
} else {
  cat("Running lookup table...\n")
  system(lookup_cmd)
  cat("Finished Querying lookup table. Filtering anchors for artifacts...\n\n")
}

result <- tryCatch({
  # Read in the lookup table stats
  lookup_stats <- fread(out_lookup_stats,
                        header = FALSE,
                        col.names = c("query", "stats"), sep = "\t")
  lookup_stats <- lookup_stats %>% mutate(anchor = anchors_to_keep$anchor)
  # Filter out anchors that have lookup table hits to artifacts
  artifact_pattern <- paste0("plas|illum|syn|arp|RF|JUNK|",
                             "Ral|purge|P,|Univec|cattle|chicken")
  anchors_to_keep <- lookup_stats %>%
    filter(!grepl(artifact_pattern, query, ignore.case = TRUE)) %>%
    select(anchor)
  # Continue with the rest of your code if no error occurs
  anchors_to_keep

  cat(paste0("Finished filtering. Kept ",
             nrow(anchors_to_keep),
             " anchors out of ",
             nrow(lookup_stats),
             " total anchors.\n"))
  cat(paste0("Keeping the top ",
             opt$num_anchors,
             " by number_nonzero_samples for further analysis...\n\n"))

}, error = function(e) {
  message("An error occurred: ", e$message)
  # Check the length of anchors_to_keep
  if (exists("anchors_to_keep") && nchar(anchors_to_keep[1,]) < 11) {
    message("Continuing execution without modifying anchors_to_keep.")
    # Continue with the rest of your code
    # For example, you can return anchors_to_keep or perform other operations
    return(anchors_to_keep)
  } else {
    message("Stopping execution as anchors_to_keep has 11 or more elements.")
    stop(e)  # Stop execution and raise the error
  }


  cat(paste0("Finished filtering. ",
             "Lookup table failed likely due to short anchors. ",
             "Kept all anchors.\n"))
  cat(paste0("Keeping the top ",
             opt$num_anchors,
             " by number_nonzero_samples for further analysis...\n\n"))

})

## select the most important anchors -----------
anchors_to_keep <- anchors_to_keep %>% 
  left_join(dt %>% select(anchor, number_nonzero_samples), by="anchor") %>%
  arrange(desc(number_nonzero_samples)) %>%
  head(opt$num_anchors)

## write the output -----------
cat(paste("Writing", nrow(anchors_to_keep), "anchors to", anchors_only_out))
cat("\n")
anchors_to_keep$anchor %>%
  write.table(anchors_only_out, row.names = FALSE,
              col.names = FALSE, quote = FALSE)
