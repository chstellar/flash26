suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(msa))


option_list <- list(
  make_option(c("--nonzero_annotations"), type = "character", default = NULL,
              help = "Path to the nonzero annotations tsv file", metavar = "character"),
  make_option(c("--clusters"), type = "character", default = NULL,
              help = "Path to the clusters tsv file", metavar = "character"),
  make_option(c("--feather_file"), type = "character", default = NULL,
              help = "Path to the X matrix feather file", metavar = "character"),
  make_option(c("--sample_seqs"), type = "character", default = NULL,
              help = "Path to the sample sequences file", metavar = "character"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Path to the metadata tsv file", metavar = "character"),
  make_option(c("--reblastp_annotations"), type = "character", default = NULL,
              help = "Optional unrestricted reblastp annotated TSV", metavar = "character"),
  make_option(c("--reblast_annotations"), type = "character", default = NULL,
              help = "Optional unrestricted reblast annotated TSV", metavar = "character"),
  make_option(c("--compactor_summary"), type = "character", default = NULL,
              help = "Optional compactor-filled plot summary TSV used as rescue labels", metavar = "character"),
  make_option(c("--compactor_anchor_len"), type = "integer", default = 31,
              help = "Anchor suffix length used for compactor rescue matching", metavar = "integer"),
  make_option(c("--target_vars"), type = "character", default = "",
              help = "Optional semicolon-delimited residual target settings used by run_adelie", metavar = "character"),
  make_option(c("--confound_vars"), type = "character", default = "",
              help = "Optional semicolon-delimited residual confounder settings used by run_adelie", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Path to set of output plots", metavar = "character"),
  make_option(c("--products"), type= "logical", default=FALSE, action="store_true",
              help = "default to using products for column names instead of genes"),
  make_option(c("--num_hits"), type="numeric", default=10,
              help = "num nonzero coefficients to plot", metavar = "numeric"),
  make_option(c("--cluster_length"), type="integer", default=NULL,
              help = "Length of each concatenated anchor-target sequence. Defaults to k value parsed from input filename, then 54.", metavar = "integer")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$cluster_length)) {
  inferred_cluster_length <- suppressWarnings(as.integer(str_extract(opt$nonzero_annotations, "_k(\\d+)_s", group = 1)))
  opt$cluster_length <- ifelse(is.na(inferred_cluster_length), 54L, inferred_cluster_length)
}

# Check if all required arguments are provided
if (is.null(opt$nonzero_annotations) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All arguments must be supplied", call. = FALSE)
}

# set known_causes to be empty (can be changed for interactive experimentation on specific datasets)
known_causes = "NNNNNNNNNNNNNNN"

# # testing
# setwd("/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline")
# opt$nonzero_annotations = "results/tuberculosis-PZAres//filter1/shiftDist-levFilter/hyena/normalized/tuberculosis-PZAres_hyena_adelie_results_top20000_k54_s54_nonzero_coefficients_blastp_annotated.tsv"
# opt$clusters = "results/tuberculosis-PZAres/filter1/shiftDist-levFilter/tuberculosis-PZAres_sequences_per_cluster_top20000-clusters_k54_s54.tsv"
# opt$feather = "/scratch/users/dcotter1/metaSPLASH_workflows_v2/tuberculosis-PZAres/tuberculosis-PZAres_hyena_top_variance_features_for_glmnet_filter1_shiftDist-levFilter_top20000_k54_s54_normalized.feather"
# opt$sample_seqs = "/scratch/users/dcotter1/metaSPLASH_workflows_v2/tuberculosis-PZAres/tuberculosis-PZAres_prepared_sequences_filter1_shiftDist-levFilter_top20000_sample_sequences.tsv"
# opt$metadata = "/oak/stanford/groups/horence/dcotter1/utility_files/metadata/metaSPLASH_metadata/tb_kim_et_al_cleaned_metadata.tsv"
# opt$output = "/oak/stanford/groups/horence/dcotter1/share/250506/test_eFac_more_blast_hits_out.pdf"

filename = data.frame(path=opt$nonzero_annotations)

filename <- filename %>%
  mutate(num_clusters = str_extract(path, "top(\\d+)", group=1)) %>%
  mutate(path=dirname(gsub("^results/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  dplyr::rename(paramater_set=path) %>%
  mutate(paramater_set=str_replace(paramater_set, "_", "/")) %>%
  separate(paramater_set, into=c("dataset", "paramater_set"), sep="/") %>%
  mutate(model=str_extract(paramater_set,
                           'hyenaHG38_normalized|hyenaHG38_unnormalized|hyenaMarlowe_normalized|hyenaMarlowe_unnormalized|esm_normalized|esm_unnormalized|hyena_normalized|hyena_unnormalized|ohe')) %>%
  mutate(filter = str_extract(paramater_set, "(filter\\d)_",group=1)) %>%
  mutate(cluster_approach = str_extract(paramater_set, "filter\\d_([A-Za-z-2]+)_", group=1))

paramaters <- filename %>% pivot_longer(everything(), names_to="paramater", values_to="value") %>% deframe()


# Define Function
get_max_abs_value <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    max(abs(nums), na.rm = TRUE)
  })
}


get_first_coef <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    nums[1]
  })
}

get_first_class <- function(x) {
  sapply(x, function(str) {
    classes <- strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]]
    classes[1]
  })
}

get_nth_coef <- function(x, n=1) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    nums[n]
  })
}

get_nth_class <- function(x,n=1) {
  sapply(x, function(str) {
    classes <- strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]]
    classes[n]
  })
}

clean_blast_label <- function(x) {
  x <- replace_na(x, "")
  x <- str_replace_all(x, "LOC\\d+[- ]*", "")
  x <- str_replace_all(x, "\\s+isoform\\s+X\\d+\\b", "")
  x <- str_replace_all(x, "\\s+transcript\\s+variant\\s+X?\\d+\\b", "")
  x <- str_replace_all(x, "\\s+variant\\s+X?\\d+\\b", "")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "\\s*[,;]\\s*$", "")
  x <- str_trim(x)
  ifelse(nchar(x) == 0, NA_character_, x)
}

collapse_blast_labels <- function(x) {
  labels <- clean_blast_label(unlist(str_split(replace_na(x, ""), ";|,")))
  labels <- labels[!is.na(labels) & nchar(labels) > 1]
  if (length(labels) == 0) {
    return(NA_character_)
  }
  paste(unique(labels), collapse=";")
}

extract_feature_qualifier <- function(features, qualifier) {
  pattern <- paste0("['\"]", qualifier, "['\"]:\\s*(?:\\[([^\\]]*)\\]|([^,}\\]]+))")
  matches <- str_match_all(replace_na(features, ""), pattern)
  map_chr(matches, function(match) {
    if (nrow(match) == 0) {
      return(NA_character_)
    }
    values <- c(match[, 2], match[, 3])
    values <- values[!is.na(values)]
    values <- unlist(str_extract_all(values, "(?<=['\"])[^'\"]+(?=['\"])"))
    values <- values[!values %in% c("None", "NA", "")]
    values <- clean_blast_label(values)
    values <- values[!is.na(values) & nchar(values) > 1]
    if (length(values) == 0) {
      return(NA_character_)
    }
    paste(unique(values), collapse=";")
  })
}

choose_feature_label <- function(products, genes, prefer_products = FALSE) {
  products <- clean_blast_label(products)
  genes <- clean_blast_label(genes)
  genes_are_loc <- !is.na(genes) & str_detect(genes, "^LOC\\d+$")
  use_products <- prefer_products | genes_are_loc | is.na(genes) | nchar(genes) < 2
  label <- ifelse(use_products & !is.na(products) & nchar(products) > 1, products, genes)
  label <- ifelse((is.na(label) | nchar(label) < 2) & !is.na(products), products, label)
  clean_blast_label(label)
}

combine_blast_labels <- function(blastp_label, blast_label) {
  pmap_chr(list(blastp_label, blast_label), function(x, y) {
    labels <- c(x, y)
    labels <- labels[!is.na(labels)]
    labels <- labels[!labels %in% c("NO MATCH", "NO PROTEIN/GENE HIT", "UNANNOTATED")]
    collapsed <- collapse_blast_labels(paste(labels, collapse=";"))
    if (is.na(collapsed)) {
      return(NA_character_)
    }
    collapsed
  })
}

make_histogram_label <- function(labels) {
  labels <- clean_blast_label(labels)
  labels <- labels[!is.na(labels) & nchar(labels) > 1]
  labels <- str_replace(labels, "^NO PROTEIN/GENE HIT$", "UNANNOTATED")
  labels <- str_replace(labels, "^NO BLAST$", "NO MATCH")
  special_labels <- intersect(c("NO TARGET", "NO MATCH", "UNANNOTATED"), unique(labels))
  real_labels <- labels[!labels %in% c("NO TARGET", "NO MATCH", "UNANNOTATED")]
  if (length(real_labels) > 0) {
    label_keys <- str_to_lower(str_replace_all(real_labels, "[^[:alnum:]]+", " "))
    real_labels <- real_labels[!duplicated(str_squish(label_keys))]
  }
  if (length(real_labels) == 0) {
    if (length(special_labels) == 0) {
      return("NO MATCH")
    }
    return(paste(special_labels, collapse=", "))
  }
  generic <- str_detect(real_labels, "(?i)hypothetical|uncharacterized|predicted protein|unnamed")
  real_labels <- real_labels[order(generic, nchar(real_labels), real_labels)]
  labels_to_show <- c(head(real_labels, 6), head(special_labels, max(0, 6 - length(real_labels))))
  paste(labels_to_show, collapse=", ")
}

format_model_metric <- function(metric_value, classes, train_only = FALSE) {
  is_regression <- any(map_lgl(classes, \(x) length(x) == 1 && x[1] == "residual"))
  metric_value <- suppressWarnings(as.numeric(metric_value))
  metric_value <- metric_value[!is.na(metric_value)]
  if (length(metric_value) == 0) {
    metric_value <- NA_real_
  } else {
    metric_value <- unique(metric_value)[1]
  }
  if (is_regression) {
    label <- ifelse(train_only, "Train R2:", "Test R2:")
    value <- ifelse(is.na(metric_value), "NA", sprintf("%.3f", metric_value))
  } else {
    label <- ifelse(train_only, "Train Accuracy:", "Accuracy:")
    value <- ifelse(is.na(metric_value), "NA", scales::label_percent(accuracy = 0.01)(metric_value))
  }
  paste(label, value)
}

format_residual_adjustment <- function(confounders) {
  if (is.null(confounders)) {
    return("")
  }
  confounders <- unique(na.omit(as.character(confounders)))
  confounders <- confounders[confounders != "" & confounders != "NA"]
  if (length(confounders) == 0) {
    return("")
  }
  str_wrap(paste("Adjusted for:", paste(confounders, collapse="; ")), width=115)
}

infer_residual_source_col <- function(category, metadata_cols) {
  if (category %in% metadata_cols) {
    return(category)
  }
  if (!str_detect(category, "_+residual")) {
    return(category)
  }
  inferred_col <- str_replace(category, "_+residual.*$", "")
  if (inferred_col %in% metadata_cols) {
    return(inferred_col)
  }
  category
}

infer_residual_focus_class <- function(category, metadata_source_col, metadata_values) {
  if (!str_detect(category, "_+residual")) {
    return(NA_character_)
  }
  prefix <- paste0(metadata_source_col, "_residual_")
  if (!startsWith(category, prefix)) {
    return(NA_character_)
  }
  focus <- substr(category, nchar(prefix) + 1, nchar(category))
  if (focus == "" || str_detect(focus, "^adjustment\\d+$")) {
    return(NA_character_)
  }
  metadata_classes <- sort(unique(na.omit(as.character(metadata_values))))
  if (focus %in% metadata_classes) {
    return(focus)
  }
  NA_character_
}

make_plotmath_other_taxa_label <- function(label) {
  label <- replace_na(label, "")
  label <- str_replace_all(label, "\\s*\\(OTHER TAXA\\)\\s*$", "")
  label <- str_replace_all(label, "\n", " ")
  label <- str_squish(label)
  label <- str_replace_all(label, "\\\\", "\\\\\\\\")
  label <- str_replace_all(label, "\"", "\\\\\"")
  paste0("atop(\"", label, "\", italic(\"(OTHER TAXA)\"))")
}

preserve_compactor_suffix <- function(label, width = 80) {
  label <- replace_na(as.character(label), "")
  is_compactor <- str_detect(label, "\\s*\\(COMPACTOR\\)\\s*$")
  base_label <- str_replace(label, "\\s*\\(COMPACTOR\\)\\s*$", "")
  ifelse(is_compactor,
         paste0(str_trunc(base_label, width=max(10, width - 12), side="right"), " (COMPACTOR)"),
         str_trunc(label, width=width, side="right"))
}

make_safe_label <- function(value) {
  label <- as.character(value)
  label <- str_replace_all(label, "[ /\\\\]", "_")
  label <- str_replace_all(label, "[^[:alnum:]._-]", "_")
  label <- paste(unlist(str_split(label, "_"))[unlist(str_split(label, "_")) != ""], collapse="_")
  ifelse(nchar(label) == 0, "value", label)
}

make_unique_name <- function(base_name, existing_names) {
  if (!base_name %in% existing_names) {
    return(base_name)
  }
  suffix <- 2
  candidate <- paste0(base_name, "_", suffix)
  while (candidate %in% existing_names) {
    suffix <- suffix + 1
    candidate <- paste0(base_name, "_", suffix)
  }
  candidate
}

target_residual_names <- function(metadata, target_col) {
  if (!target_col %in% colnames(metadata)) {
    return(character())
  }
  values <- metadata[[target_col]]
  values_chr <- str_trim(as.character(values))
  values_chr[values_chr %in% c("", "nan", "NaN", "NA", "None")] <- NA_character_
  values_num <- suppressWarnings(as.numeric(values_chr))
  nonmissing <- !is.na(values_chr)
  is_numeric <- any(nonmissing) && all(!is.na(values_num[nonmissing]))
  if (is_numeric) {
    return(paste0(target_col, "_residual"))
  }
  categories <- sort(unique(na.omit(values_chr)))
  if (length(categories) <= 2) {
    return(paste0(target_col, "_residual"))
  }
  paste0(target_col, "_residual_", map_chr(categories, make_safe_label))
}

parse_residual_adjustment_map <- function(target_vars, confound_vars, metadata) {
  target_vars <- str_trim(str_remove_all(replace_na(target_vars, ""), "^['\"]|['\"]$"))
  confound_vars <- str_trim(str_remove_all(replace_na(confound_vars, ""), "^['\"]|['\"]$"))
  if (target_vars == "" || confound_vars == "") {
    return(tibble(metadata_category=character(), confounders=character()))
  }
  targets <- str_split(target_vars, ";", simplify=FALSE)[[1]] %>% str_trim()
  confound_groups <- str_split(confound_vars, ";", simplify=FALSE)[[1]] %>% str_trim()
  if (length(confound_groups) < length(targets)) {
    confound_groups <- c(confound_groups, rep("", length(targets) - length(confound_groups)))
  }
  residual_specs <- tibble(target=character(), confounds=list())
  for (i in seq_along(targets)) {
    target <- targets[[i]]
    if (target == "" || str_to_lower(target) == "all") {
      next
    }
    confounds <- str_split(confound_groups[[i]], ",", simplify=FALSE)[[1]] %>% str_trim()
    confounds <- confounds[confounds != ""]
    if (length(confounds) > 0) {
      residual_specs <- bind_rows(residual_specs, tibble(target=target, confounds=list(confounds)))
    }
  }
  if (nrow(residual_specs) == 0) {
    return(tibble(metadata_category=character(), confounders=character()))
  }
  repeated_targets <- residual_specs %>% count(target) %>% filter(n > 1) %>% pull(target)
  existing_names <- colnames(metadata)
  out <- tibble(metadata_category=character(), confounders=character())
  for (i in seq_len(nrow(residual_specs))) {
    target <- residual_specs$target[[i]]
    confounds <- residual_specs$confounds[[i]]
    if (any(str_to_lower(confounds) == "all")) {
      confounds <- setdiff(colnames(metadata), c("sample_name", target))
      confounds <- confounds[!str_detect(str_replace_all(confounds, "__", "_"), "_residual")]
    }
    confound_label <- paste(confounds, collapse=", ")
    base_names <- target_residual_names(metadata, target)
    if (target %in% repeated_targets) {
      base_names <- paste0(base_names, "_adjustment", i)
    }
    for (base_name in base_names) {
      category_name <- make_unique_name(base_name, existing_names)
      existing_names <- c(existing_names, category_name)
      out <- bind_rows(out, tibble(metadata_category=category_name, confounders=confound_label))
    }
  }
  out
}

first_numeric_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  x[1]
}

single_line_text <- function(x) {
  x <- replace_na(as.character(x), "")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  ifelse(nchar(x) == 0, NA_character_, x)
}

metadata_count_vector <- function(metadata_string) {
  entries <- unlist(str_split(replace_na(metadata_string, ""), "/"))
  entries <- entries[nchar(entries) > 0]
  if (length(entries) == 0) {
    return(numeric())
  }
  parsed <- str_match(entries, "^(.+):([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)$")
  parsed <- parsed[!is.na(parsed[, 1]), , drop=FALSE]
  if (nrow(parsed) == 0) {
    return(numeric())
  }
  names <- parsed[, 2]
  values <- suppressWarnings(as.numeric(parsed[, 3]))
  keep <- !is.na(values) &
    !str_detect(names, "^(mean_|median_|sd_|n$|n_)")
  values <- values[keep]
  names(values) <- names[keep]
  values
}

metadata_entropy_stats <- function(metadata_string, total_samples) {
  counts <- metadata_count_vector(metadata_string)
  counts <- counts[counts > 0]
  total <- suppressWarnings(as.numeric(total_samples))
  if (length(total) == 0 || is.na(total)) {
    total <- sum(counts)
  }
  if (length(counts) == 0 || total <= 0) {
    return(tibble(
      metadata_entropy = NA_real_,
      metadata_normalized_entropy = NA_real_,
      metadata_specificity_score = NA_real_,
      dominant_metadata = NA_character_,
      dominant_metadata_count = NA_real_,
      dominant_metadata_fraction = NA_real_
    ))
  }
  probs <- counts / sum(counts)
  entropy <- -sum(probs * log(probs))
  normalized_entropy <- ifelse(length(counts) > 1, entropy / log(length(counts)), 0)
  dominant_idx <- which.max(counts)
  tibble(
    metadata_entropy = entropy,
    metadata_normalized_entropy = normalized_entropy,
    metadata_specificity_score = (1 - normalized_entropy) * total,
    dominant_metadata = names(counts)[dominant_idx],
    dominant_metadata_count = as.numeric(counts[dominant_idx]),
    dominant_metadata_fraction = as.numeric(counts[dominant_idx]) / sum(counts)
  )
}

read_optional_tsv <- function(path) {
  if (is.null(path) || is.na(path) || nchar(path) == 0 || !file.exists(path) || file.info(path)$size == 0) {
    return(data.table())
  }
  fread(path)
}

clean_sequence_candidate <- function(x) {
  x <- replace_na(as.character(x), "")
  x <- str_replace_all(x, "[^ACGTNacgtn]", "")
  toupper(x)
}

sequence_suffix <- function(x, n) {
  x <- clean_sequence_candidate(x)
  ifelse(nchar(x) >= n, str_sub(x, -n), x)
}

make_compactor_summary_label_dt <- function(path) {
  compactor_dt <- read_optional_tsv(path)
  if (nrow(compactor_dt) == 0 || !"compactor_annotation" %in% colnames(compactor_dt)) {
    return(tibble(metadata_category=character(), cluster=character(), feature=character(),
                  sequence=character(), compactor_label=character(),
                  compactor_anchor=character(),
                  compactor_identity=numeric(), compactor_qcovs=numeric()))
  }
  if (!"Blast Label" %in% colnames(compactor_dt)) {
    compactor_dt$`Blast Label` <- NA_character_
  }
  if (!"sequence" %in% colnames(compactor_dt)) {
    compactor_dt$sequence <- NA_character_
  }
  compactor_dt %>%
    mutate(sequence = clean_sequence_candidate(sequence)) %>%
    mutate(compactor_anchor = sequence_suffix(sequence, opt$compactor_anchor_len)) %>%
    mutate(compactor_label = ifelse(has_restricted_label(compactor_annotation),
                                    compactor_annotation, `Blast Label`)) %>%
    filter(has_restricted_label(compactor_label), !is.na(sequence), nchar(sequence) > 0) %>%
    mutate(compactor_identity = first_numeric_or_na(identity),
           compactor_qcovs = first_numeric_or_na(qcovs)) %>%
    group_by(metadata_category, cluster, feature, sequence, compactor_anchor) %>%
    summarise(compactor_label = collapse_blast_labels(paste(unique(na.omit(compactor_label)), collapse=";")),
              compactor_identity = first_numeric_or_na(compactor_identity),
              compactor_qcovs = first_numeric_or_na(compactor_qcovs),
              .groups="drop")
}

has_restricted_label <- function(label) {
  !is.na(label) & nchar(label) > 1 &
    !label %in% c("NO MATCH", "NO TARGET", "UNANNOTATED", "NO PROTEIN/GENE HIT")
}

make_reblastp_label_dt <- function(reblastp_dt, category) {
  if (nrow(reblastp_dt) == 0 || !"query" %in% colnames(reblastp_dt)) {
    return(tibble(sequence=character(), outside_taxid_label=character(),
                  outside_taxid_identity=numeric(), outside_taxid_qcovs=numeric()))
  }
  reblastp_dt %>%
    filter(metadata_category == category) %>%
    mutate(sequence = str_remove(query, "^cluster_\\d+_")) %>%
    mutate(outside_taxid_label = clean_blast_label(str_remove_all(stitle, "\\[.+\\]$|MULTISPECIES:\\s|, partial"))) %>%
    group_by(sequence) %>%
    summarise(outside_taxid_label = collapse_blast_labels(paste(unique(na.omit(outside_taxid_label)), collapse=";")),
              outside_taxid_identity = first_numeric_or_na(identity),
              outside_taxid_qcovs = first_numeric_or_na(qcovs),
              .groups="drop")
}

make_reblast_label_dt <- function(reblast_dt, category) {
  if (nrow(reblast_dt) == 0 || !"query" %in% colnames(reblast_dt)) {
    return(tibble(sequence=character(), outside_taxid_label=character(),
                  outside_taxid_identity=numeric(), outside_taxid_qcovs=numeric()))
  }
  reblast_dt %>%
    filter(metadata_category == category) %>%
    mutate(sequence = str_remove(query, "^cluster_\\d+_")) %>%
    mutate(feature_text = paste(replace_na(as.character(features), ""),
                                replace_na(as.character(features_all), ""),
                                sep=";")) %>%
    mutate(outside_products = extract_feature_qualifier(feature_text, "product")) %>%
    mutate(outside_genes = extract_feature_qualifier(feature_text, "gene")) %>%
    mutate(outside_taxid_label = choose_feature_label(outside_products, outside_genes, opt$products)) %>%
    group_by(sequence) %>%
    summarise(outside_taxid_label = collapse_blast_labels(paste(unique(na.omit(outside_taxid_label)), collapse=";")),
              outside_taxid_identity = first_numeric_or_na(identity),
              outside_taxid_qcovs = first_numeric_or_na(qcovs),
              .groups="drop")
}

# function to read the nth cluster out of the sample sequences file
read_nth_cluster <- function(file_path, n, cluster_length) {
  # Calculate start and end positions for the nth cluster
  start <- n * cluster_length + 1  # n is now 0-indexed
  end <- (n + 1) * cluster_length

  # Construct the awk command
  awk_command <- sprintf("awk 'NR>1 {print $1 \"\t\" substr($2, %d, %d)}' %s", start, cluster_length, file_path)

  # Execute the command and read the result
  result <- read.table(text = system(awk_command, intern = TRUE),
                       col.names = c("sample_name", "sequence"),
                       stringsAsFactors = FALSE)

  return(result)
}

# Function to calculate distance to the most abundant sequence after MSA of unique sequences
calculate_distance_and_align <- function(sequences) {
  # Remove sequences with long strings of N's
  valid_sequences <- sequences[!grepl("N{10,}", sequences)]

  if (length(valid_sequences) == 0) {
    return(list(distances = rep(-1, length(sequences)),
                aligned_seqs = rep("", length(sequences))))
  }

  # Get unique sequences and their counts
  seq_table <- table(valid_sequences)
  unique_seqs <- names(seq_table)

  # Perform MSA on unique sequences
  msa_result <- msa(DNAStringSet(unique_seqs), method = "ClustalOmega", order = "input")

  # Convert MSA result to character vectors
  aligned_seqs <- as.character(msa_result)

  # Function to count leading and trailing dashes
  count_leading_dashes <- function(seq) {
    nchar(str_extract(seq, "^(\\-+)[ACTGN]", group=1))
  }

  count_trailing_dashes <- function(seq) {
    nchar(str_extract(seq, "[ACTGN]+(\\-+)$", group=1))
  }

  # Determine the maximum number of leading and trailing dashes
  max_leading_dashes <- max(c(count_leading_dashes(aligned_seqs),0), na.rm=T)
  max_trailing_dashes <- max(c(count_trailing_dashes(aligned_seqs),0), na.rm=T)

  # Trim the determined number of dashes from each sequence
  trim_dashes <- function(seq) {
    substr(seq, max_leading_dashes + 1, nchar(seq) - max_trailing_dashes)
  }

  trimmed_aligned_seqs <- sapply(aligned_seqs, trim_dashes)

  # Find the most abundant sequence
  most_abundant <- trimmed_aligned_seqs[which.max(seq_table)]

  # Calculate distances for unique sequences
  unique_distances <- stringdist(trimmed_aligned_seqs, most_abundant, method = "lv")

  # Map distances and aligned sequences back to all valid sequences
  all_valid_distances <- unique_distances[match(valid_sequences, unique_seqs)]
  all_valid_aligned_seqs <- aligned_seqs[match(valid_sequences, unique_seqs)]

  # Map distances and aligned sequences back to original sequences, including those with N's
  all_distances <- rep(-1, length(sequences))
  all_aligned_seqs <- rep("", length(sequences))
  valid_indices <- !grepl("N{10,}", sequences)
  all_distances[valid_indices] <- all_valid_distances
  all_aligned_seqs[valid_indices] <- all_valid_aligned_seqs

  return(list(distances = all_distances, aligned_seqs = all_aligned_seqs))
}

# read in input files
dt <- fread(opt$nonzero_annotations)
if (TRUE) {dt2 <- fread(gsub("blastp_annotated", "blast_annotated", opt$nonzero_annotations))}
for (compactor_col in c("compactor_annotation", "compactor_query", "compactor_length",
                        "compactor_exact_support", "compactor_raw_annotation")) {
  if (!compactor_col %in% colnames(dt)) {
    dt[[compactor_col]] <- NA_character_
  }
  if (!compactor_col %in% colnames(dt2)) {
    dt2[[compactor_col]] <- NA_character_
  }
}
if (!"confounders" %in% colnames(dt)) {
  dt$confounders <- NA_character_
}
if (!"confounders" %in% colnames(dt2)) {
  dt2$confounders <- NA_character_
}
dt_reblastp <- read_optional_tsv(opt$reblastp_annotations)
dt_reblast <- read_optional_tsv(opt$reblast_annotations)
compactor_summary_label_dt <- make_compactor_summary_label_dt(opt$compactor_summary)
message("Compactor summary rescue labels loaded: ", nrow(compactor_summary_label_dt))
if (nrow(dt_reblastp) > 0 && !"qcovs" %in% colnames(dt_reblastp)) {
  dt_reblastp$qcovs <- NA
}
if (nrow(dt_reblast) > 0) {
  if (!"qcovs" %in% colnames(dt_reblast)) dt_reblast$qcovs <- NA
  if (!"features" %in% colnames(dt_reblast)) dt_reblast$features <- NA
  if (!"features_all" %in% colnames(dt_reblast)) dt_reblast$features_all <- NA
}
all_clusters <- fread(opt$clusters) %>% select(-kmer)
feather_dt <- feather::read_feather(opt$feather)
all_metadata <- fread(opt$metadata)
residual_adjustment_map <- parse_residual_adjustment_map(opt$target_vars, opt$confound_vars, all_metadata)
if (nrow(residual_adjustment_map) > 0) {
  dt <- dt %>%
    left_join(residual_adjustment_map %>% dplyr::rename(confounders_from_args=confounders),
              by="metadata_category") %>%
    mutate(confounders = coalesce(confounders, confounders_from_args)) %>%
    select(-confounders_from_args)
  dt2 <- dt2 %>%
    left_join(residual_adjustment_map %>% dplyr::rename(confounders_from_args=confounders),
              by="metadata_category") %>%
    mutate(confounders = coalesce(confounders, confounders_from_args)) %>%
    select(-confounders_from_args)
}

categories <- bind_rows(
  dt %>% select(any_of(c("metadata_category", "accuracy"))),
  dt2 %>% select(any_of(c("metadata_category", "accuracy")))
) %>%
  filter(!is.na(metadata_category)) %>%
  mutate(accuracy = suppressWarnings(as.numeric(accuracy))) %>%
  group_by(metadata_category) %>%
  summarise(accuracy = max(accuracy, na.rm=TRUE), .groups="drop") %>%
  mutate(accuracy = ifelse(is.infinite(accuracy), NA_real_, accuracy)) %>%
  arrange(desc(accuracy), metadata_category) %>%
  pull(metadata_category)

if (str_detect(opt$nonzero_annotations, "adelie-train-only")) {
  dt <- dt %>% mutate(accuracy=train_accuracy)
}

#category = "ampicillin_RIS"

pdf(opt$output, width=12, height=8)
all_features_summary <- data.table()
all_blastp_summary <- data.table()
all_blast_summary <- data.table()

# write a title page first
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, paramaters['dataset'])
text(5, 7, paramaters['filter'])
text(5, 6, paramaters['cluster_approach'])
text(5, 5, paramaters['model'])
text(5, 4, paste("At most", paramaters['num_clusters'], "clusters"))
text(5,3, paste(Sys.Date()))

for (category in categories) {
  tryCatch({
    new_dt <- dt %>% filter(is.na(query)) %>% select(-query) %>% left_join(all_clusters %>% mutate(query = paste0(cluster, "_", seq)) %>% select(-seq), by="cluster", relationship="many-to-many") %>%
      filter(!is.na(query))
    summ_dt <- bind_rows(dt, new_dt) %>% filter(!is.na(query)) %>% filter(metadata_category==category) %>%
      mutate(first_coef=get_first_coef(coefficients)) %>% mutate(max_coefficient=abs(first_coef)) %>%
      mutate(second_coef = get_nth_coef(coefficients,2), second_class=get_nth_class(classes, 2)) %>%
      arrange(-max_coefficient) %>% mutate(first_class=get_first_class(classes)) %>%
      mutate(annotation = str_remove_all(stitle, "\\[.+\\]$|MULTISPECIES:\\s|, partial")) %>%
      mutate(annotation = ifelse(has_restricted_label(compactor_annotation),
                                 compactor_annotation, annotation)) %>%
      rowwise() %>%
      mutate(classes=list(str_split_1(gsub("\\[|\\]", "", classes),pattern=","))) %>%
      ungroup() %>%
      select(metadata_category, accuracy, classes, first_class, first_coef, second_coef, second_class, max_coefficient, cluster, feature, query, identity, qcovs, annotation, confounders) %>%
      mutate(query = str_remove(query, "cluster_\\d+_")) %>%
      distinct(cluster,annotation,query,.keep_all = T) %>% group_by(cluster)


    summ_dt <- summ_dt %>% group_by(cluster,query) %>%
      mutate(label=ifelse(!is_empty(unique(na.omit(annotation))), paste0(unique(na.omit(annotation)),collapse=";"), NA)) %>%
      distinct(cluster, query, label, .keep_all=T) %>% ungroup()

    summ_dt_blastp_only <- summ_dt

    if (TRUE) {
      if (!"qcovs" %in% colnames(dt2)) {
        dt2$qcovs <- NA
      }
      if (!"features" %in% colnames(dt2)) {
        dt2$features <- NA
      }
      if (!"features_all" %in% colnames(dt2)) {
        dt2$features_all <- NA
      }
      summ_dt2 <- dt2 %>% filter(metadata_category==category) %>%
        mutate(feature_text = paste(replace_na(as.character(features), ""),
                                    replace_na(as.character(features_all), ""),
                                    sep=";")) %>%
        separate_longer_delim(feature_text, delim = "},") %>%
        mutate(products=extract_feature_qualifier(feature_text, "product")) %>%
        mutate(genes=extract_feature_qualifier(feature_text, "gene")) %>%
        select(-feature_text) %>% mutate(first_coef=get_first_coef(coefficients)) %>% mutate(max_coefficient=abs(first_coef)) %>%
        arrange(-max_coefficient) %>% mutate(first_class=get_first_class(classes)) %>%
        rowwise() %>%
        mutate(classes=list(str_split_1(gsub("\\[|\\]", "", classes),pattern=","))) %>%
        ungroup() %>%
        select(metadata_category, accuracy, classes, first_class, first_coef, max_coefficient, cluster, feature, query, identity, qcovs, products, genes, confounders, compactor_annotation) %>%
        mutate(query = str_remove(query, "^cluster_\\d+_")) %>%
        group_by(cluster) %>%
        ungroup() %>%
        distinct(cluster,products,query,genes,.keep_all = T) %>% group_by(cluster)

      summ_dt2 <- summ_dt2 %>%
        mutate(label = choose_feature_label(products, genes, opt$products)) %>%
        mutate(label = ifelse(has_restricted_label(compactor_annotation),
                              compactor_annotation, label)) %>%
        group_by(cluster,query) %>%
        mutate(label=ifelse(!is_empty(unique(na.omit(label))), paste0(unique(na.omit(label)),collapse=";"), NA)) %>%
        distinct(cluster, query, label, .keep_all=T) %>% ungroup()
      summ_dt2 <- summ_dt2 %>% select(cluster, query, identity, qcovs, label) %>% dplyr::rename(label2=label)
      summ_dt <- summ_dt %>% left_join(summ_dt2 %>%
                                         select(cluster, query, identity, qcovs, label2), by=c("cluster", "query")) %>%
        dplyr::rename(identity=`identity.x`, qcovs=`qcovs.x`) %>%
        mutate(identity = ifelse(is.na(label) & !is.na(label2), `identity.y`, identity)) %>%
        mutate(qcovs = ifelse(is.na(label) & !is.na(label2), `qcovs.y`, qcovs)) %>%
        mutate(label_blastp = label, label_blast = label2) %>%
        mutate(label = combine_blast_labels(label_blastp, label_blast)) %>%
        mutate(label = ifelse(is.na(label) & !is.na(label_blast), label_blast, label)) %>%
        mutate(label = ifelse(is.na(label) & !is.na(label_blastp), label_blastp, label)) %>%
        mutate(label = ifelse(is.na(label) | nchar(label)<2, NA, label))
    }

    summ_dt <- summ_dt %>% group_by(feature) %>% mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)) & (is.na(label)) & (!is.na(identity) | !is.na(identity.y)), "NO PROTEIN/GENE HIT", label))

    summ_dt <- summ_dt %>% group_by(cluster) %>%
      mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)), "NO MATCH", label)) %>%
      mutate(hypothetical=length(unique(label))>1 & sum(str_detect(label, "(?i)hypothetical|uncharacterized"))>0) %>%
      mutate(hypothetical=replace_na(hypothetical, FALSE)) %>%
      rowwise() %>%
      mutate(label = str_replace(label, " ,", ", ") %>% str_replace(" ;", "; ")) %>%
      mutate(label = map2_vec(label, hypothetical, \(x,y) if (y) {str_c(str_trim(unlist(str_split(x, ";"))[str_detect(unlist(str_split(x, ";")), "(?i)hypothetical|uncharact", negate=T)]),sep = ",", collapse=",")} else {x})) %>%
      ungroup()

    reblastp_label_dt <- make_reblastp_label_dt(dt_reblastp, category) %>%
      dplyr::rename(outside_taxid_label_blastp = outside_taxid_label,
                    outside_taxid_identity_blastp = outside_taxid_identity,
                    outside_taxid_qcovs_blastp = outside_taxid_qcovs)
    reblast_label_dt <- make_reblast_label_dt(dt_reblast, category) %>%
      dplyr::rename(outside_taxid_label_blast = outside_taxid_label,
                    outside_taxid_identity_blast = outside_taxid_identity,
                    outside_taxid_qcovs_blast = outside_taxid_qcovs)
    outside_taxid_label_dt <- full_join(reblastp_label_dt, reblast_label_dt, by="sequence") %>%
      mutate(outside_taxid_label = combine_blast_labels(outside_taxid_label_blastp, outside_taxid_label_blast)) %>%
      mutate(outside_taxid_label = ifelse(is.na(outside_taxid_label) & !is.na(outside_taxid_label_blast), outside_taxid_label_blast, outside_taxid_label)) %>%
      mutate(outside_taxid_label = ifelse(is.na(outside_taxid_label) & !is.na(outside_taxid_label_blastp), outside_taxid_label_blastp, outside_taxid_label)) %>%
      mutate(outside_taxid_identity = coalesce(outside_taxid_identity_blastp, outside_taxid_identity_blast),
             outside_taxid_qcovs = coalesce(outside_taxid_qcovs_blastp, outside_taxid_qcovs_blast)) %>%
      select(sequence, outside_taxid_label, outside_taxid_identity, outside_taxid_qcovs)

    category_compactor_label_dt <- compactor_summary_label_dt %>%
      filter(metadata_category == category) %>%
      select(cluster, feature, sequence, compactor_label, compactor_identity, compactor_qcovs)
    category_compactor_anchor_label_dt <- compactor_summary_label_dt %>%
      filter(metadata_category == category) %>%
      select(cluster, feature, compactor_anchor, compactor_anchor_label=compactor_label,
             compactor_anchor_identity=compactor_identity, compactor_anchor_qcovs=compactor_qcovs) %>%
      group_by(cluster, feature, compactor_anchor) %>%
      summarise(compactor_anchor_label = collapse_blast_labels(paste(unique(na.omit(compactor_anchor_label)), collapse=";")),
                compactor_anchor_identity = first_numeric_or_na(compactor_anchor_identity),
                compactor_anchor_qcovs = first_numeric_or_na(compactor_anchor_qcovs),
                .groups="drop")
    if (nrow(category_compactor_label_dt) > 0) {
      message("Compactor rescue labels for ", category, ": ", nrow(category_compactor_label_dt))
    }

    blastp_all_dt <- summ_dt_blastp_only %>%
      group_by(feature) %>%
      mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)) & (is.na(label)) & !is.na(identity), "NO PROTEIN/GENE HIT", label)) %>%
      group_by(cluster) %>%
      mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)), "NO MATCH", label)) %>%
      ungroup() %>%
      mutate(blast_source = "blastp") %>%
      mutate(classes = map_chr(classes, \(x) paste(x, collapse=","))) %>%
      select(any_of(c("metadata_category", "accuracy", "classes", "first_class", "first_coef", "second_coef", "second_class",
                      "max_coefficient", "cluster", "feature", "query", "identity", "qcovs", "label", "blast_source")))

    blast_all_dt <- summ_dt %>%
      mutate(blast_source = "blast") %>%
      mutate(label = combine_blast_labels(label_blastp, label_blast)) %>%
      mutate(label = ifelse(is.na(label) & !is.na(label_blast), label_blast, label)) %>%
      mutate(label = ifelse(is.na(label) & !is.na(label_blastp), label_blastp, label)) %>%
      mutate(label = ifelse(is.na(label) & (!is.na(`identity.y`) | !is.na(`qcovs.y`)),
                            "UNANNOTATED", label)) %>%
      mutate(label = ifelse(is.na(label), "NO MATCH", label)) %>%
      mutate(classes = map_chr(classes, \(x) paste(x, collapse=","))) %>%
      select(any_of(c("metadata_category", "accuracy", "classes", "first_class", "first_coef", "second_coef", "second_class",
                      "max_coefficient", "cluster", "feature", "query", "identity.y", "qcovs.y", "label", "blast_source"))) %>%
      dplyr::rename(identity = `identity.y`, qcovs = `qcovs.y`)

    all_blastp_summary <- bind_rows(all_blastp_summary, blastp_all_dt)
    all_blast_summary <- bind_rows(all_blast_summary, blast_all_dt)

    hist_direct_label_dt <- dt2 %>%
      filter(metadata_category == category) %>%
      mutate(feature_text = paste(replace_na(as.character(features), ""),
                                  replace_na(as.character(features_all), ""),
                                  sep=";")) %>%
      mutate(hist_products = extract_feature_qualifier(feature_text, "product")) %>%
      mutate(hist_genes = extract_feature_qualifier(feature_text, "gene")) %>%
      mutate(hist_label = choose_feature_label(hist_products, hist_genes, opt$products)) %>%
      mutate(hist_label = ifelse(has_restricted_label(compactor_annotation),
                                 compactor_annotation, hist_label)) %>%
      mutate(hist_label = case_when(
        str_detect(query, "NNNNNNNN") ~ "NO TARGET",
        !is.na(hist_label) & nchar(hist_label) > 1 ~ hist_label,
        !is.na(identity) | !is.na(qcovs) ~ "UNANNOTATED",
        TRUE ~ NA_character_
      )) %>%
      mutate(first_coef=get_first_coef(coefficients)) %>%
      mutate(max_coefficient=abs(first_coef)) %>%
      select(cluster, feature, max_coefficient, label=hist_label)

    hist_label_dt <- bind_rows(
      summ_dt %>%
        mutate(label=ifelse(label=="",annotation,label)) %>%
        mutate(sequence = query) %>%
        left_join(outside_taxid_label_dt, by="sequence") %>%
        mutate(label = case_when(
          str_detect(query, "NNNNNNNN") ~ "NO TARGET",
          !is.na(label) & nchar(label) > 1 ~ label,
          !is.na(outside_taxid_label) & nchar(outside_taxid_label) > 1 ~ outside_taxid_label,
          !is.na(identity) | !is.na(qcovs) | !is.na(`identity.y`) | !is.na(`qcovs.y`) ~ "UNANNOTATED",
          TRUE ~ "NO MATCH"
        )) %>%
        select(cluster, feature, max_coefficient, label),
      hist_direct_label_dt
      ,
      category_compactor_label_dt %>%
        select(cluster, feature, label=compactor_label) %>%
        left_join(summ_dt %>% ungroup() %>%
                    select(cluster, feature, max_coefficient) %>%
                    distinct(),
                  by=c("cluster", "feature"))
    )

    plot_dt <- hist_label_dt %>%
      mutate(largest_coef=max(max_coefficient)) %>%
      mutate(coef_mag=max_coefficient/largest_coef) %>%
      group_by(cluster, feature, coef_mag) %>%
      summarise(label=make_histogram_label(label), .groups="drop") %>%
      mutate(label = str_replace(label, " ,", ", ") %>% str_replace(" ;", "; ")) %>%
      arrange(-coef_mag) %>%
      mutate(rank=row_number()) %>%
      mutate(color="no_blast") %>%
      mutate(has_real_label = !str_detect(label, "^(NO TARGET|NO MATCH|UNANNOTATED)(,\\s*(NO TARGET|NO MATCH|UNANNOTATED))*$")) %>%
      mutate(color=ifelse(has_real_label, "blast", color)) %>%
      mutate(color = ifelse(grepl(known_causes, label, ignore.case=T), "known_cause", color)) %>%
      mutate(label = str_wrap(preserve_compactor_suffix(label, width=100), width = 25)) %>%
      mutate(label = replace_na(label, ""))

    metric_subtitle <- format_model_metric(
      summ_dt$accuracy,
      summ_dt$classes,
      str_detect(opt$nonzero_annotations, "adelie-train-only")
    )
    adjustment_subtitle <- format_residual_adjustment(summ_dt$confounders)
    hist_subtitle <- paste(c(metric_subtitle, adjustment_subtitle)[c(metric_subtitle, adjustment_subtitle) != ""],
                           collapse="\n")

    dataset <- str_extract(opt$nonzero_annotations, "results/([A-Za-z\\d-]+)/filter", group=1)
    make_title <- paste(category, "in", dataset)

    p <- plot_dt %>% head(opt$num_hits) %>%
      ggplot(aes(x=rank, y=coef_mag, fill=color, label=label)) +
      geom_col() +
      geom_text(aes(y=coef_mag + 0.05,hjust=0),angle=45,size=3) +
      scale_y_continuous("Magnitude relative to\nlargest nonzero coefficient",
                         limits = c(0,1.5),
                         labels=scales::label_percent(), breaks=seq(0,1,0.25),
                         expand=c(0,0)) +
      xlab("Rank of nonzero coefficient (by magnitude)") +
      scale_fill_manual(breaks=c("known_cause","blast","no_blast"), values=c("forestgreen", "pink", "grey"),
                        labels=c("Known Cause", "Blast hit", "No blast/annotation")) +
      scale_x_continuous(breaks=seq(1,10,1)) +
      ggtitle(make_title,
              subtitle = hist_subtitle) +
      theme_pubr() +
      theme(legend.position="none",
            plot.subtitle = element_text(size=8, lineheight=0.95))

    print(p)

    # select the top N coefficients
    interesting_clusters <- summ_dt %>% distinct(cluster, feature, max_coefficient) %>%
      arrange(-max_coefficient) %>% head(opt$num_hits)

    # filter metadata for only current category. Regression residual names may not
    # exist in the raw metadata sheet, so infer the raw source column when possible.
    metadata_source_col <- infer_residual_source_col(category, colnames(all_metadata))
    if (!metadata_source_col %in% colnames(all_metadata)) {
      message(paste("Skipping detailed plots for", category, "because no matching metadata column was found."))
      next
    }
    my_metadata <- as_tibble(all_metadata) %>% select(sample_name, !!metadata_source_col) %>% dplyr::rename(metadata:=!!metadata_source_col)
    residual_focus_class <- infer_residual_focus_class(category, metadata_source_col, my_metadata$metadata)
    metadata_numeric <- suppressWarnings(as.numeric(my_metadata$metadata))
    category_is_regression <- any(map_lgl(summ_dt$classes, \(x) length(x) == 1 && x[1] == "residual"))
    is_quantitative_target <- category_is_regression &&
      all(is.na(my_metadata$metadata) == is.na(metadata_numeric)) &&
      length(unique(na.omit(metadata_numeric))) >= 2
    if (is_quantitative_target) {
      my_metadata <- my_metadata %>% mutate(metadata = as.numeric(metadata))
    }
    for (j in 1:nrow(interesting_clusters)) {

      my_cluster = interesting_clusters[j,]$cluster
      my_feature = interesting_clusters[j,]$feature


      dt_sub <- summ_dt %>% filter(cluster==my_cluster, feature==my_feature)

      first_beta = unique(dt_sub$first_coef)
      first_class = unique(dt_sub$first_class)
      all_classes = dt_sub[1,]$classes %>% unlist()
      classes_to_plot <- all_classes
      if (length(all_classes) == 1 && all_classes[1] == "residual" && !is_quantitative_target) {
        all_classes <- sort(unique(na.omit(as.character(my_metadata$metadata))))
        if (!is.na(residual_focus_class)) {
          classes_to_plot <- residual_focus_class
        } else if (length(all_classes) > 0) {
          classes_to_plot <- all_classes
        }
        first_class <- all_classes[1]
      }

      seq_sub <- read_nth_cluster(opt$sample_seqs, as.numeric(str_extract(my_cluster, "\\d+")), opt$cluster_length)
      emb_sub <- feather_dt %>% select(sample_name, !!my_feature)

      seq_sub <- seq_sub %>%
        left_join(emb_sub, by="sample_name") %>%
        dplyr::rename(embedding:=!!my_feature) %>%
        mutate(embedding=embedding*first_beta)

      distances <- calculate_distance_and_align(seq_sub$sequence)

      seq_sub$aligned_sequence <- distances$aligned_seq
      seq_sub$lev_dist <- distances$distances

      direct_blast_label_dt <- dt2 %>%
        filter(metadata_category == category, cluster == my_cluster, feature == my_feature) %>%
        mutate(sequence = str_remove(query, "^cluster_\\d+_")) %>%
        mutate(feature_text = paste(replace_na(as.character(features), ""),
                                    replace_na(as.character(features_all), ""),
                                    sep=";")) %>%
        mutate(direct_products = extract_feature_qualifier(feature_text, "product")) %>%
        mutate(direct_genes = extract_feature_qualifier(feature_text, "gene")) %>%
        mutate(direct_blast_label = choose_feature_label(direct_products, direct_genes, opt$products)) %>%
        mutate(direct_blast_label = ifelse(has_restricted_label(compactor_annotation),
                                           compactor_annotation, direct_blast_label)) %>%
        group_by(sequence) %>%
        summarise(direct_blast_label = collapse_blast_labels(paste(unique(na.omit(direct_blast_label)), collapse=";")),
                  direct_identity = first_numeric_or_na(identity),
                  direct_qcovs = first_numeric_or_na(qcovs),
                  .groups="drop")

      label_dt <- dt_sub %>%
        mutate(any_identity = coalesce(`identity.y`, identity),
               any_qcovs = coalesce(`qcovs.y`, qcovs)) %>%
        mutate(combined_label = combine_blast_labels(label_blastp, label_blast)) %>%
        mutate(combined_label = ifelse(is.na(combined_label) & !is.na(label_blast), label_blast, combined_label)) %>%
        mutate(combined_label = ifelse(is.na(combined_label) & !is.na(label_blastp), label_blastp, combined_label)) %>%
        select(query, accuracy, any_identity, any_qcovs, combined_label, label2) %>%
        dplyr::rename(label=combined_label) %>%
        dplyr::rename(identity=any_identity, qcovs=any_qcovs) %>%
        dplyr::rename(sequence=query)

      if (is_quantitative_target) {
        summ_sub_dt <- seq_sub %>% left_join(my_metadata,by="sample_name") %>%
          filter(!is.na(metadata)) %>%
          group_by(sequence,embedding,lev_dist,aligned_sequence) %>%
          summarise(total_samples=n(),
                    mean_metadata=mean(metadata, na.rm=TRUE),
                    median_metadata=median(metadata, na.rm=TRUE),
                    sd_metadata=sd(metadata, na.rm=TRUE),
                    .groups="drop") %>%
          relocate(aligned_sequence, .after="sequence")
      } else {
        summ_sub_dt <- seq_sub %>% left_join(my_metadata,by="sample_name") %>%
          group_by(sequence,embedding,lev_dist,aligned_sequence,metadata) %>%
          summarise(metadata_count=n(), .groups="drop") %>%
          filter(metadata %in% all_classes) %>%
          pivot_wider(id_cols=everything(), names_from=metadata, values_from=metadata_count) %>%
          relocate(aligned_sequence, .after="sequence")
      }

      summ_sub_dt <- summ_sub_dt %>%
        mutate(compactor_anchor = sequence_suffix(sequence, opt$compactor_anchor_len)) %>%
        left_join(label_dt, by="sequence") %>%
        left_join(direct_blast_label_dt, by="sequence") %>%
        left_join(category_compactor_label_dt %>%
                    filter(cluster == my_cluster, feature == my_feature) %>%
                    select(sequence, compactor_label, compactor_identity, compactor_qcovs),
                  by="sequence") %>%
        left_join(category_compactor_anchor_label_dt %>%
                    filter(cluster == my_cluster, feature == my_feature),
                  by=c("cluster", "feature", "compactor_anchor")) %>%
        left_join(outside_taxid_label_dt, by="sequence")

      p_sub <- summ_sub_dt %>%
        mutate(compactor_label = coalesce(compactor_label, compactor_anchor_label),
               compactor_identity = coalesce(compactor_identity, compactor_anchor_identity),
               compactor_qcovs = coalesce(compactor_qcovs, compactor_anchor_qcovs)) %>%
        mutate(label = ifelse(has_restricted_label(compactor_label), compactor_label, label)) %>%
        mutate(identity = coalesce(compactor_identity, identity),
               qcovs = coalesce(compactor_qcovs, qcovs)) %>%
        mutate(label = ifelse((is.na(label) | label == "NO MATCH") &
                                 !is.na(direct_blast_label) & nchar(direct_blast_label) > 1,
                              direct_blast_label, label)) %>%
        mutate(identity = coalesce(identity, direct_identity),
               qcovs = coalesce(qcovs, direct_qcovs)) %>%
        mutate(outside_taxid_only = !has_restricted_label(label) &
                 !is.na(outside_taxid_label) & nchar(outside_taxid_label) > 1) %>%
        mutate(label = ifelse(outside_taxid_only, outside_taxid_label, label)) %>%
        mutate(identity = ifelse(outside_taxid_only & is.na(identity), outside_taxid_identity, identity),
               qcovs = ifelse(outside_taxid_only & is.na(qcovs), outside_taxid_qcovs, qcovs)) %>%
        mutate(has_blast_hit = !is.na(identity) | !is.na(qcovs) |
                 (!is.na(label) & label != "NO MATCH") |
                 (!is.na(label2) & nchar(label2) > 1)) %>%
        mutate(label = ifelse((is.na(label) | nchar(label)<3) & !is.na(label2) & nchar(label2)>3,
                              label2, label)) %>%
        mutate(label = ifelse(str_detect(sequence, "NNNNNNNN"), "NO TARGET", label)) %>%
        mutate(label = ifelse(is.na(label) & has_blast_hit, "UNANNOTATED", label)) %>%
        ungroup() %>%
        mutate(label = map_chr(label, collapse_blast_labels)) %>%
        mutate(label=gsub(",Pbp5","",label)) %>%
        mutate(label = single_line_text(label)) %>%
        mutate(label=ifelse(is.na(label), "NO MATCH", label)) %>%
        dplyr::rename(`Blast Label` = label) %>%
        mutate(label_identity = ifelse(identity == 100 | is.na(identity), "-", paste0(round(identity,2), "%"))) %>%
        mutate(label_coverage = ifelse(qcovs == 100 | is.na(qcovs), "-", paste0(round(qcovs,2), "%"))) %>%
        mutate(label_identity = replace_na(label_identity, "-")) %>%
        mutate(label_coverage = replace_na(label_coverage, "-")) %>%
        rowwise() %>%
        mutate(label_quality = ifelse(label_coverage != "-" | label_identity != "-",
                                      paste0("I:", str_replace(label_identity, "-", "100%"),
                                             "; C:", str_replace(label_coverage, "-", "100%")), "")) %>%
        mutate(point_label = case_when(
          `Blast Label` %in% c("NO MATCH", "NO TARGET") ~ as.character(`Blast Label`),
          `Blast Label` %in% c("NO PROTEIN/GENE HIT", "UNANNOTATED") &
            nchar(label_quality) > 0 ~ paste(`Blast Label`, label_quality, sep="\n"),
          `Blast Label` %in% c("NO PROTEIN/GENE HIT", "UNANNOTATED") ~ as.character(`Blast Label`),
          TRUE ~ as.character(`Blast Label`)
        )) %>%
        mutate(point_label = ifelse(outside_taxid_only,
                                    paste0(str_replace(point_label, "\\s*\\(OTHER TAXA\\)\\s*$", ""),
                                           " (OTHER TAXA)"),
                                    point_label)) %>%
        mutate(point_label = ifelse(outside_taxid_only,
                                    paste0(str_trunc(str_replace(point_label, "\\s*\\(OTHER TAXA\\)\\s*$", ""),
                                                     width=65),
                                           " (OTHER TAXA)"),
                                    preserve_compactor_suffix(point_label, width=80))) %>%
        mutate(point_label = str_wrap(point_label, width=28)) %>%
        mutate(point_label_expr = ifelse(outside_taxid_only,
                                         make_plotmath_other_taxa_label(point_label),
                                         NA_character_)) %>%
        mutate(point_label_color = ifelse(outside_taxid_only, "outside_taxid", "regular")) %>%
        ungroup()

      if (is_quantitative_target) {
        p_sub <- p_sub %>% mutate(color_value = mean_metadata)
        p2 <- p_sub %>%
          ggplot(aes(x=embedding, y=lev_dist, color=color_value, size=total_samples,
                     label=point_label)) +
          geom_vline(xintercept = 0, lty="dashed") +
          geom_point(stroke=1.4) +
          scale_y_continuous(breaks=scales::breaks_width(1)) +
          scale_size_continuous(trans = "log", name = "Total Samples",
                                breaks = c(1, 10, 100, 1000, 10000),
                                limits = c(1, 10000), labels = scales::label_log()) +
          ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color == "regular"),
                                   size=3.2, color="black", max.overlaps=Inf,
                                   min.segment.length=0, segment.color="grey45",
                                   segment.size=0.25, box.padding=0.5,
                                   point.padding=0.55, force=3) +
          ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color == "outside_taxid"),
                                   aes(label=point_label_expr),
                                   size=3.2, color="grey25", max.overlaps=Inf,
                                   min.segment.length=0, segment.color="grey45",
                                   segment.size=0.25, box.padding=0.5,
                                   point.padding=0.55, force=3, parse=TRUE) +
          scale_color_gradient(paste0("Mean\n", metadata_source_col),
                               low = "blue", high = "red") +
          theme_minimal() + xlab(expression("Embedding" ~ "\u00D7" ~ beta)) +
          ylab("Levenshtein Distance\n(to most abundant anchor-target)") +
          ggtitle(paste(category, my_cluster, sep=" | "),
                  subtitle=paste("Color: mean observed", metadata_source_col))
        print(p2)
        summ_out_dt <- p_sub %>% select(-label2) %>%
          mutate(metadata = paste0("mean_", metadata_source_col, ":", round(mean_metadata, 6),
                                   "/median_", metadata_source_col, ":", round(median_metadata, 6),
                                   "/n:", total_samples),
                 metadata_category = category, cluster=my_cluster, feature=my_feature)
      } else {
        missing_class_cols <- setdiff(all_classes, colnames(p_sub))
        for (missing_class_col in missing_class_cols) {
          p_sub[[missing_class_col]] <- 0
        }
        p_sub <- p_sub %>%
          mutate(across(all_of(all_classes), \(x) replace_na(x, 0))) %>%
          mutate(total_samples = rowSums(across(all_of(all_classes))))

        if (length(classes_to_plot) == 0) {
          classes_to_plot <- first_class
        }
        for (class_to_plot in classes_to_plot) {
          if (!class_to_plot %in% colnames(p_sub)) {
            next
          }
          p_class <- p_sub %>%
            mutate(class_proportion = ifelse(total_samples > 0,
                                             !!sym(class_to_plot) / total_samples,
                                             NA_real_))
          p2 <- p_class %>%
            ggplot(aes(x=embedding, y=lev_dist, color=class_proportion,
                       size=total_samples, label=point_label)) +
            geom_vline(xintercept = 0, lty="dashed") +
            geom_point(stroke=1.4) +
            scale_y_continuous(breaks=scales::breaks_width(1)) +
            scale_size_continuous(trans = "log", name = "Total Samples",
                                  breaks = c(1, 10, 100, 1000, 10000),
                                  limits = c(1, 10000), labels = scales::label_log()) +
            ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color == "regular"),
                                     size=3.2, color="black", max.overlaps=Inf,
                                     min.segment.length=0, segment.color="grey45",
                                     segment.size=0.25, box.padding=0.5,
                                     point.padding=0.55, force=3) +
            ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color == "outside_taxid"),
                                     aes(label=point_label_expr),
                                     size=3.2, color="grey25", max.overlaps=Inf,
                                     min.segment.length=0, segment.color="grey45",
                                     segment.size=0.25, box.padding=0.5,
                                     point.padding=0.55, force=3, parse=TRUE) +
            scale_color_gradient(paste0("Proportion\n", class_to_plot),
                                 low = "blue", high = "red", limits = c(0, 1)) +
            theme_minimal() + xlab(expression("Embedding" ~ "\u00D7" ~ beta)) +
            ylab("Levenshtein Distance\n(to most abundant anchor-target)") +
            ggtitle(paste(category, my_cluster, sep=" | "),
                    subtitle=paste("Color: proportion", class_to_plot))
          print(p2)
        }

        summ_out_dt <- p_sub %>% select(-label2) %>% ungroup() %>%
          mutate(across(all_of(all_classes), \(x) paste(cur_column(), replace_na(x, 0),sep=":"))) %>%
          unite(col=metadata, all_of(all_classes), sep="/") %>%
          mutate(metadata_category = category, cluster=my_cluster, feature=my_feature)
      }

      all_features_summary <- bind_rows(all_features_summary, summ_out_dt)
    }

  }, error = function(e) {
    message(paste("Error processing category:", category, "\n", e$message))
    # Optionally log the error or take other actions
  })

}

dev.off()

if (!"metadata_category" %in% colnames(all_features_summary)) {
  stop("No feature plot summary rows were generated; check category-level errors above.", call. = FALSE)
}

all_features_summary <- all_features_summary %>%
  mutate(across(where(is.character), single_line_text))
all_blastp_summary <- all_blastp_summary %>%
  mutate(across(where(is.character), single_line_text))
all_blast_summary <- all_blast_summary %>%
  mutate(across(where(is.character), single_line_text))

unannotated_summary <- all_features_summary %>%
  filter(is.na(`Blast Label`) |
           `Blast Label` %in% c("NO MATCH", "UNANNOTATED", "NO PROTEIN/GENE HIT")) %>%
  mutate(total_samples = suppressWarnings(as.numeric(total_samples))) %>%
  mutate(entropy_stats = map2(metadata, total_samples, metadata_entropy_stats)) %>%
  unnest(entropy_stats) %>%
  arrange(desc(metadata_specificity_score), desc(total_samples), metadata_normalized_entropy) %>%
  relocate(metadata_category, feature, cluster)

all_features_summary %>% relocate(metadata_category, feature, cluster) %>% write_tsv(file = str_replace(opt$output, ".pdf", "_summary.tsv"))
unannotated_summary %>% write_tsv(file = str_replace(opt$output, ".pdf", "_unannotated.tsv"))
all_blastp_summary %>% relocate(metadata_category, feature, cluster) %>% write_tsv(file = str_replace(opt$nonzero_annotations, "blastp_annotated.tsv$", "blastp_all.tsv"))
all_blast_summary %>% relocate(metadata_category, feature, cluster) %>% write_tsv(file = str_replace(gsub("blastp_annotated", "blast_annotated", opt$nonzero_annotations), "blast_annotated.tsv$", "blast_all.tsv"))
