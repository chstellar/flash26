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
# These numbers can be changed via command line arguments
# and by changing the varialbes in the config file for the pipeline.

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input"), "Input file", type = "character"),
  make_option(c("-o", "--output"), "Output file", type = "character"),
  make_option(c("-n", "--num_anchors"), "Number of anchors to select",
    type = "integer", default = 1000000
  ),
  make_option(c("-e", "--effect_size"), "Effect size threshold",
    type = "numeric", default = 0.6
  ),
  make_option(c("--target_rank"), "Target rank to use when selecting rank-specific SPLASH score columns",
    type = "integer", default = 1
  ),
  make_option(c("--pre_lookup_multiplier"), "Multiplier for the candidate buffer to send through artifact lookup",
    type = "numeric", default = 3
  ),
  make_option(c("--sort_threads"), "Threads to use for GNU sort",
    type = "integer", default = 1
  ),
  make_option(c("-l", "--lookup_table"), "Lookup table file",
    type = "character"
  ),
  make_option(c("--splash_bin"), "Path to SPLASH binary folder",
    type = "character", default = "splash-2.11.6/"
  ),
  make_option(c("--temp_dir"),
    "Temporary directory to store intermediate files",
    type = "character"
  )
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

cat("select_anchors_filter1.R build: disk-streaming-topn-v1\n")

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
    paste0(opt$temp_dir, "/")
  )
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
} else {
  temp_dir <- file.path(dirname(opt$output), "tmp/")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
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
cat(paste("Target rank:", opt$target_rank))
cat("\n")
cat(paste("Pre-lookup candidate multiplier:", opt$pre_lookup_multiplier))
cat("\n")
cat(paste("Sort threads:", opt$sort_threads))
cat("\n")
cat(paste("Lookup table file:", opt$lookup_table))
cat("\n")
cat(paste("SPLASH binary folder:", opt$splash_bin))
cat("\n")
cat(paste("Temporary directory:", temp_dir))
cat("\n")
cat("###################################################################\n\n")


read_score_lines <- function(path, n = 2) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, open = "rt")
  } else {
    con <- file(path, open = "rt")
  }
  on.exit(close(con))
  readLines(con, n = n, warn = FALSE)
}

score_header_lines <- read_score_lines(opt$input, n = 2)
if (length(score_header_lines) < 2) {
  stop("Input score file must contain a header and at least one data row")
}
header_names <- strsplit(score_header_lines[1], "\t", fixed = TRUE)[[1]]
first_data_fields <- strsplit(score_header_lines[2], "\t", fixed = TRUE)[[1]]
headers <- as.data.frame(as.list(first_data_fields), stringsAsFactors = FALSE)
names(headers) <- header_names

if (nchar(headers$anchor[1]) < 12) {
  opt$effect_size <- 0.1 # change if using short anchors
}

choose_ranked_column <- function(headers, pattern, target_rank, label) {
  candidates <- grep(pattern, names(headers), ignore.case = TRUE)
  if (length(candidates) == 0) {
    stop(paste0("Could not find ", label, " column matching pattern: ", pattern))
  }
  candidate_names <- names(headers)[candidates]
  rank_patterns <- c(
    paste0("target[^[:alnum:]]*", target_rank, "([^[:digit:]]|$)"),
    paste0("rank[^[:alnum:]]*", target_rank, "([^[:digit:]]|$)"),
    paste0("(^|[^[:digit:]])", target_rank, "$")
  )
  rank_match <- rep(FALSE, length(candidate_names))
  for (rank_pattern in rank_patterns) {
    rank_match <- rank_match | grepl(rank_pattern, candidate_names, ignore.case = TRUE)
  }
  if (sum(rank_match) == 1) {
    return(candidates[which(rank_match)])
  }
  if (sum(rank_match) > 1) {
    return(candidates[which(rank_match)[1]])
  }
  if (target_rank <= length(candidates)) {
    return(candidates[target_rank])
  }
  warning(paste0(
    "Requested target_rank=", target_rank,
    " but only found ", length(candidates), " ", label,
    " columns. Falling back to the last matching column: ",
    candidate_names[length(candidate_names)]
  ))
  candidates[length(candidates)]
}

effect_size_bin_col <- choose_ranked_column(
  headers, "effect_size_bin", opt$target_rank, "effect_size_bin"
)
nonzero_samples_col <- choose_ranked_column(
  headers, "number_nonzero_samples", opt$target_rank, "number_nonzero_samples"
)
cols_to_read <- c(1, effect_size_bin_col, nonzero_samples_col)
col_aliases <- setNames(names(headers)[cols_to_read], c("anchor", "effect_size_bin", "number_nonzero_samples"))
cat(paste0(
  "Using SPLASH columns: anchor=", names(headers)[1],
  ", effect_size_bin=", names(headers)[effect_size_bin_col],
  ", number_nonzero_samples=", names(headers)[nonzero_samples_col],
  "\n"
))

## stream-filter anchors without materializing the full SPLASH score table -----
EFFECT_SIZE_CUTOFF <- opt$effect_size
input_reader <- ifelse(grepl("\\.gz$", opt$input), "gzip -dc", "cat")
pre_lookup_n <- max(opt$num_anchors, ceiling(opt$num_anchors * opt$pre_lookup_multiplier))
sort_threads <- max(1, opt$sort_threads)

run_cmd <- function(cmd, label) {
  status <- system(cmd)
  if (!isTRUE(status == 0)) {
    stop(paste0(label, " failed with exit status ", status))
  }
}

stream_prefix <- paste(input_reader, shQuote(opt$input))
nonzero_values_file <- file.path(temp_dir, "number_nonzero_samples_after_effect.txt") %>% gsub("//", "/", .)
nonzero_frequency_file <- file.path(temp_dir, "number_nonzero_samples_frequency.tsv") %>% gsub("//", "/", .)
candidate_ranked_file <- file.path(temp_dir, "anchor_candidates_ranked.tsv") %>% gsub("//", "/", .)
anchor_file <- file.path(temp_dir, "all_anchors.txt") %>% gsub("//", "/", .)
anchors_fasta_file <- file.path(temp_dir, "all_anchors.fasta") %>% gsub("//", "/", .)
artifact_filtered_file <- file.path(temp_dir, "anchor_candidates_artifact_filtered.tsv") %>% gsub("//", "/", .)
out_lookup_stats <- file.path(temp_dir, "lookup_stats.txt") %>% gsub("//", "/", .)

cat(paste(
  "Streaming anchors from",
  opt$input, "with effect size >=",
  opt$effect_size
))
cat("\n\n")

effect_filter_cmd <- paste0(
  stream_prefix,
  " | awk -F'\\t' -v e=", EFFECT_SIZE_CUTOFF,
  " -v ec=", effect_size_bin_col,
  " -v nc=", nonzero_samples_col,
  " 'NR>1 && ($ec+0) >= e {print $nc+0}' > ",
  shQuote(nonzero_values_file)
)
run_cmd(effect_filter_cmd, "Writing nonzero-sample values")

n_effect <- as.numeric(system(paste("wc -l <", shQuote(nonzero_values_file)), intern = TRUE))
if (is.na(n_effect) || n_effect < 1) {
  stop("No anchors passed the effect size cutoff.")
}
quantile_rank <- max(1, ceiling(n_effect * 0.1))
sample_cutoff <- as.numeric(system(paste0(
  "sort -n ", shQuote(nonzero_values_file),
  " | awk -v target=", quantile_rank,
  " 'NR==target{print; found=1; exit} END{if(!found) print 0}'"
), intern = TRUE))
if (is.na(sample_cutoff)) {
  stop("Could not compute number_nonzero_samples cutoff.")
}
cat(paste0(
  "Effect-size filter retained ", n_effect,
  " rows; 10th percentile number_nonzero_samples cutoff is ",
  sample_cutoff, ".\n"
))

frequency_cmd <- paste0(
  "awk -v s=", sample_cutoff,
  " '($1+0) >= s {count[$1+0]++} END{for (v in count) print v \"\\t\" count[v]}' ",
  shQuote(nonzero_values_file),
  " | LC_ALL=C sort -k1,1nr > ",
  shQuote(nonzero_frequency_file)
)
run_cmd(frequency_cmd, "Writing nonzero-sample frequencies")

nonzero_frequency <- fread(
  nonzero_frequency_file,
  header = FALSE,
  col.names = c("number_nonzero_samples", "n")
)
nonzero_frequency[, cumulative_n := cumsum(n)]
candidate_sample_cutoff <- nonzero_frequency[
  cumulative_n >= pre_lookup_n,
  number_nonzero_samples[1]
]
if (is.na(candidate_sample_cutoff)) {
  candidate_sample_cutoff <- sample_cutoff
}
candidate_sample_cutoff <- max(sample_cutoff, candidate_sample_cutoff)
cat(paste0(
  "Top-buffer number_nonzero_samples cutoff is ",
  candidate_sample_cutoff,
  " for pre-lookup cap ", pre_lookup_n,
  ".\n"
))

candidate_cmd <- paste0(
  stream_prefix,
  " | awk -F'\\t' -v OFS='\\t' -v e=", EFFECT_SIZE_CUTOFF,
  " -v s=", candidate_sample_cutoff,
  " -v ec=", effect_size_bin_col,
  " -v nc=", nonzero_samples_col,
  " 'NR>1 && ($ec+0) >= e && ($nc+0) >= s {print $1, $nc+0}' ",
  " | LC_ALL=C sort --parallel=", sort_threads,
  " -S 80% -k2,2nr -T ", shQuote(temp_dir),
  " | head -n ", pre_lookup_n,
  " > ", shQuote(candidate_ranked_file)
)
run_cmd(candidate_cmd, "Writing ranked anchor candidates")

n_candidates <- as.numeric(system(paste("wc -l <", shQuote(candidate_ranked_file)), intern = TRUE))
cat(paste0(
  "Ranked and retained the top ", n_candidates,
  " candidate anchors after nonzero-sample filtering for artifact lookup",
  " (pre-lookup cap: ", pre_lookup_n, ").\n"
))

write_anchor_cmd <- paste0(
  "awk -F'\\t' '{print $1}' ", shQuote(candidate_ranked_file),
  " > ", shQuote(anchor_file)
)
run_cmd(write_anchor_cmd, "Writing anchor list")

write_fasta_cmd <- paste0(
  "awk -F'\\t' '{print \">\"$1\"\\n\"$1}' ",
  shQuote(candidate_ranked_file),
  " > ", shQuote(anchors_fasta_file)
)
run_cmd(write_fasta_cmd, "Writing anchor FASTA")

## run lookup table to filter out artifacts -----------
lookup_cmd <- paste0(
  shQuote(file.path(opt$splash_bin, "lookup_table")),
  " query --kmer_skip 1 ",
  "--truncate_paths --stats_fmt with_stats ",
  shQuote(anchors_fasta_file),
  " ", shQuote(opt$lookup_table), " ", shQuote(out_lookup_stats)
)
if (file.exists(out_lookup_stats)) {
  unlink(out_lookup_stats)
}

artifact_pattern <- paste0(
  "plas|illum|syn|arp|RF|JUNK|",
  "Ral|purge|P,|Univec|cattle|chicken"
)

lookup_ok <- TRUE
if (nchar(headers$anchor[1]) < 11) {
  lookup_ok <- FALSE
  cat("Skipping lookup table because anchors are shorter than 11 bp.\n")
} else {
  cat("Running lookup table...\n")
  lookup_status <- system(lookup_cmd)
  lookup_ok <- isTRUE(lookup_status == 0) && file.exists(out_lookup_stats)
  if (!lookup_ok) {
    cat("Lookup table failed; continuing without artifact filtering.\n")
  } else {
    cat("Finished querying lookup table. Filtering anchors for artifacts...\n\n")
  }
}

if (lookup_ok) {
  filter_artifact_cmd <- paste0(
    "paste ", shQuote(candidate_ranked_file), " ", shQuote(out_lookup_stats),
    " | awk -F'\\t' -v OFS='\\t' -v pat=", shQuote(artifact_pattern),
    " 'BEGIN{IGNORECASE=1} $3 !~ pat {print $1, $2}' ",
    " > ", shQuote(artifact_filtered_file)
  )
  run_cmd(filter_artifact_cmd, "Filtering artifact anchors")
} else {
  if (!file.copy(candidate_ranked_file, artifact_filtered_file, overwrite = TRUE)) {
    stop("Could not copy unfiltered candidate anchors after lookup failure/skip.")
  }
}

n_after_lookup <- as.numeric(system(paste("wc -l <", shQuote(artifact_filtered_file)), intern = TRUE))
cat(paste0(
  "Kept ", n_after_lookup,
  " anchors after artifact filtering. Keeping the top ",
  opt$num_anchors, " by number_nonzero_samples.\n\n"
))
if (!is.na(n_after_lookup) && n_after_lookup < opt$num_anchors) {
  warning(paste0(
    "Only ", n_after_lookup,
    " anchors remained after artifact filtering, fewer than requested num_anchors=",
    opt$num_anchors,
    ". Increase --pre_lookup_multiplier if you need a fuller output."
  ))
}

write_output_cmd <- paste0(
  "head -n ", opt$num_anchors, " ", shQuote(artifact_filtered_file),
  " | awk -F'\\t' '{print $1}' > ", shQuote(anchors_only_out)
)
dir.create(dirname(anchors_only_out), recursive = TRUE, showWarnings = FALSE)
run_cmd(write_output_cmd, "Writing selected anchors")

n_out <- as.numeric(system(paste("wc -l <", shQuote(anchors_only_out)), intern = TRUE))
cat(paste("Writing", n_out, "anchors to", anchors_only_out))
cat("\n")
