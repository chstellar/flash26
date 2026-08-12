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
  make_option(c("--target_vars"), type = "character", default = "",
              help = "Optional semicolon-delimited residual target settings used by run_adelie", metavar = "character"),
  make_option(c("--confound_vars"), type = "character", default = "",
              help = "Optional semicolon-delimited residual confounder settings used by run_adelie", metavar = "character"),
  make_option(c("--compactor_summary"), type = "character", default = "",
              help = "Optional compactor-filled plot summary TSV used only as an annotation-label lookup", metavar = "character"),
  make_option(c("--compactor_selected"), type = "character", default = "",
              help = "Optional selected compactor TSV used to backfill compactor_sequence by anchor suffix", metavar = "character"),
  make_option(c("--compactor_seed_annotations"), type = "character", default = "",
              help = "Optional seed annotation TSV used to backfill compactor_sequence by exact plotted extendor", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Path to set of output plots", metavar = "character"),
  make_option(c("--products"), type= "logical", default=FALSE, action="store_true",
              help = "default to using products for column names instead of genes"),
  make_option(c("--num_hits"), type="numeric", default=10,
              help = "num nonzero coefficients to plot", metavar = "numeric"),
  make_option(c("--cluster_length"), type="integer", default=NULL,
              help = "Length of each concatenated anchor-target sequence. Defaults to k value parsed from input filename, then 54.", metavar = "integer"),
  make_option(c("--taxid_name_cache"), type="character", default=NULL,
              help = "Optional TSV cache for resolving BLAST staxids to scientific species names.", metavar = "character")
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

message("plot_blast_annotations_each_feature.R build: compactor-seed-seq-v11")

if (is.null(opt$taxid_name_cache) || is.na(opt$taxid_name_cache) || nchar(opt$taxid_name_cache) == 0) {
  opt$taxid_name_cache <- file.path(dirname(opt$output), "blast_taxid_species_cache.tsv")
}

is_missing_path <- function(path) {
  is.null(path) || is.na(path) || nchar(path) == 0 || !file.exists(path) || file.info(path)$size == 0
}

infer_compactor_aux_path <- function(summary_path, suffix_glob) {
  if (is.null(summary_path) || is.na(summary_path) || nchar(summary_path) == 0) {
    return("")
  }
  summary_dir <- dirname(summary_path)
  summary_base <- basename(summary_path)
  run_prefix <- str_replace(summary_base, "_nonzero_coefficients_blast_annotated_plots_summary_compactor\\.tsv$", "")
  if (identical(run_prefix, summary_base)) {
    return("")
  }
  matches <- sort(Sys.glob(file.path(summary_dir, paste0(run_prefix, suffix_glob))))
  matches <- matches[file.exists(matches) & file.info(matches)$size > 0]
  if (length(matches) == 0) {
    return("")
  }
  matches[[1]]
}

if (is_missing_path(opt$compactor_selected)) {
  inferred_selected <- infer_compactor_aux_path(opt$compactor_summary, "_compactor_*_selected.tsv")
  if (nchar(inferred_selected) > 0) {
    opt$compactor_selected <- inferred_selected
  }
}
if (is_missing_path(opt$compactor_seed_annotations)) {
  inferred_seed_annotations <- infer_compactor_aux_path(opt$compactor_summary, "_compactor_*_seed_annotations.tsv")
  if (nchar(inferred_seed_annotations) > 0) {
    opt$compactor_seed_annotations <- inferred_seed_annotations
  }
}
message(paste0("Compactor selected sequence lookup path: ",
               ifelse(is_missing_path(opt$compactor_selected), "<missing>", opt$compactor_selected)))
message(paste0("Compactor seed sequence lookup path: ",
               ifelse(is_missing_path(opt$compactor_seed_annotations), "<missing>", opt$compactor_seed_annotations)))

# set known_causes to be empty (can be changed for interactive experimentation on specific datasets)
known_causes = "NNNNNNNNNNNNNNN"
within_taxid_label_color <- "#E64B35"
outside_taxid_label_color <- "#4DBBD5"
histogram_bar_color <- "#F8766D"
no_taxon_label_color <- "#7A7A7A"

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
coerce_scalar_text <- function(value) {
  if (length(value) == 0 || all(is.na(value))) {
    return(NA_character_)
  }
  if (is.list(value)) {
    value <- unlist(value, recursive = TRUE, use.names = FALSE)
  }
  value <- value[!is.na(value)]
  if (length(value) == 0) {
    return(NA_character_)
  }
  value <- paste(as.character(value), collapse = ",")
  value <- str_trim(value)
  if (nchar(value) == 0 || value %in% c("NA", "NaN", "NULL")) {
    return(NA_character_)
  }
  value
}

parse_coef_values <- function(value) {
  value <- coerce_scalar_text(value)
  if (is.na(value)) {
    return(numeric(0))
  }
  # Coefficients appear in a few serialized shapes across old/new FLASH runs:
  # [-0.1,0.2], c(-0.1, 0.2), "-0.1 0.2", or list-like scalar text.
  # Extract numeric tokens directly so plotting does not depend on one delimiter.
  nums_chr <- str_extract_all(
    value,
    "[-+]?(?:\\d*\\.\\d+|\\d+\\.?\\d*)(?:[eE][-+]?\\d+)?"
  )[[1]]
  nums <- suppressWarnings(as.numeric(nums_chr))
  nums[!is.na(nums)]
}

parse_class_values <- function(value) {
  value <- coerce_scalar_text(value)
  if (is.na(value)) {
    return(character(0))
  }
  value <- gsub("^\\[|\\]$", "", value)
  classes <- str_trim(strsplit(value, ",", fixed = TRUE)[[1]])
  classes[!is.na(classes) & nchar(classes) > 0]
}

get_max_abs_value <- function(x) {
  vapply(seq_along(x), function(i) {
    nums <- parse_coef_values(x[[i]])
    if (length(nums) == 0) {
      return(NA_real_)
    }
    max(abs(nums), na.rm = TRUE)
  }, numeric(1))
}


get_first_coef <- function(x) {
  vapply(seq_along(x), function(i) {
    nums <- parse_coef_values(x[[i]])
    if (length(nums) < 1) {
      return(NA_real_)
    }
    nums[1]
  }, numeric(1))
}

get_first_class <- function(x) {
  vapply(seq_along(x), function(i) {
    classes <- parse_class_values(x[[i]])
    if (length(classes) < 1) {
      return(NA_character_)
    }
    classes[1]
  }, character(1))
}

get_nth_coef <- function(x, n=1) {
  vapply(seq_along(x), function(i) {
    nums <- parse_coef_values(x[[i]])
    if (length(nums) < n) {
      return(NA_real_)
    }
    nums[n]
  }, numeric(1))
}

get_nth_class <- function(x,n=1) {
  vapply(seq_along(x), function(i) {
    classes <- parse_class_values(x[[i]])
    if (length(classes) < n) {
      return(NA_character_)
    }
    classes[n]
  }, character(1))
}

clean_blast_label <- function(x) {
  x <- as.character(x)
  x <- replace_na(x, "")
  x <- str_replace_all(x, "\\s*\\[[^\\]]+\\]\\s*$", "")
  x <- str_replace(x, "^RecName:\\s*Full=([^;]+).*$", "\\1")
  x <- str_replace(x, "^SubName:\\s*Full=([^;]+).*$", "\\1")
  x <- str_replace(x, "^AltName:\\s*Full=([^;]+).*$", "\\1")
  x <- str_replace(x, "(^|;\\s*)Short=([^;]+).*$", "\\2")
  x <- str_replace_all(x, "\\bRecName:\\s*Full=", "")
  x <- str_replace_all(x, "\\bAltName:\\s*Full=", "")
  x <- str_replace_all(x, "\\bSubName:\\s*Full=", "")
  x <- str_replace_all(x, "\\bFlags:\\s*[^;]+;?", "")
  x <- str_replace_all(x, "LOC\\d+[- ]*", "")
  x <- str_replace_all(x, "\\s+isoform\\s+X\\d+\\b", "")
  x <- str_replace_all(x, "\\s+transcript\\s+variant\\s+X?\\d+\\b", "")
  x <- str_replace_all(x, "\\s+variant\\s+X?\\d+\\b", "")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "\\s*[,;]\\s*$", "")
  x <- str_trim(x)
  x <- ifelse(str_detect(x, regex("uncharacteri[sz]ed protein|hypothetical protein|predicted protein|unnamed protein",
                                  ignore_case=TRUE)),
              "UNCHARACTERISED", x)
  ifelse(nchar(x) == 0, NA_character_, x)
}

collapse_blast_labels <- function(x) {
  labels <- clean_blast_label(unlist(str_split(replace_na(as.character(x), ""), ";|,")))
  labels <- labels[!is.na(labels) & nchar(labels) > 1]
  if (length(labels) == 0) {
    return(NA_character_)
  }
  paste(unique(labels), collapse=";")
}

extract_blast_species <- function(x) {
  x <- replace_na(as.character(x), "")
  matches <- str_match_all(x, "\\[([^\\]]+)\\]")
  map_chr(matches, function(match) {
    if (nrow(match) == 0) {
      return(NA_character_)
    }
    clean_blast_label(match[nrow(match), 2])
  })
}

extract_taxid_values <- function(x) {
  x <- coerce_scalar_text(x)
  if (is.na(x)) {
    return(character(0))
  }
  vals <- str_extract_all(x, "\\d+")[[1]]
  unique(vals[!is.na(vals) & nchar(vals) > 0])
}

taxid_cache_env <- new.env(parent = emptyenv())
taxid_cache_env$loaded <- FALSE
taxid_cache_env$map <- setNames(character(0), character(0))

load_taxid_species_cache <- function(cache_path) {
  if (isTRUE(taxid_cache_env$loaded)) {
    return(invisible(NULL))
  }
  taxid_cache_env$loaded <- TRUE
  if (!is.null(cache_path) && file.exists(cache_path) && file.info(cache_path)$size > 0) {
    cache_dt <- tryCatch(fread(cache_path, colClasses = "character"), error = function(e) data.table())
    if (nrow(cache_dt) > 0 && all(c("taxid", "species") %in% colnames(cache_dt))) {
      cache_dt <- cache_dt %>%
        filter(!is.na(taxid), nchar(taxid) > 0, !is.na(species), nchar(species) > 0) %>%
        distinct(taxid, .keep_all = TRUE)
      taxid_cache_env$map <- setNames(cache_dt$species, cache_dt$taxid)
    }
  }
  invisible(NULL)
}

save_taxid_species_cache <- function(cache_path) {
  if (is.null(cache_path) || is.na(cache_path) || nchar(cache_path) == 0) {
    return(invisible(NULL))
  }
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  cache_dt <- tibble(taxid = names(taxid_cache_env$map),
                     species = unname(taxid_cache_env$map)) %>%
    filter(!is.na(taxid), nchar(taxid) > 0, !is.na(species), nchar(species) > 0) %>%
    distinct(taxid, .keep_all = TRUE) %>%
    arrange(suppressWarnings(as.numeric(taxid)))
  fwrite(cache_dt, cache_path, sep = "\t")
  invisible(NULL)
}

fetch_taxid_species_from_ncbi <- function(taxids) {
  taxids <- unique(taxids[!is.na(taxids) & nchar(taxids) > 0])
  if (length(taxids) == 0) {
    return(setNames(character(0), character(0)))
  }
  fetched <- setNames(rep(NA_character_, length(taxids)), taxids)
  batches <- split(taxids, ceiling(seq_along(taxids) / 100))
  for (batch in batches) {
    url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
      "db=taxonomy&retmode=xml&id=",
      URLencode(paste(batch, collapse = ","), reserved = FALSE)
    )
    xml <- tryCatch({
      old_timeout <- getOption("timeout")
      options(timeout = max(20, min(old_timeout, 60)))
      on.exit(options(timeout = old_timeout), add = TRUE)
      paste(readLines(url, warn = FALSE), collapse = "\n")
    }, error = function(e) {
      message("Could not resolve NCBI taxonomy names for taxids ",
              paste(batch, collapse = ","), ": ", e$message)
      NA_character_
    })
    if (is.na(xml) || nchar(xml) == 0) {
      next
    }
    taxa <- str_split(xml, "<Taxon>", simplify = FALSE)[[1]]
    for (taxon_xml in taxa) {
      taxid <- str_match(taxon_xml, "<TaxId>([^<]+)</TaxId>")[,2]
      species <- str_match(taxon_xml, "<ScientificName>([^<]+)</ScientificName>")[,2]
      if (!is.na(taxid) && taxid %in% batch && !is.na(species) && nchar(species) > 0) {
        fetched[[taxid]] <- species
      }
    }
    Sys.sleep(0.34)
  }
  fetched[!is.na(fetched) & nchar(fetched) > 0]
}

taxids_to_species <- function(staxids, cache_path = opt$taxid_name_cache) {
  load_taxid_species_cache(cache_path)
  taxid_list <- lapply(staxids, extract_taxid_values)
  all_taxids <- unique(unlist(taxid_list, use.names = FALSE))
  missing_taxids <- setdiff(all_taxids, names(taxid_cache_env$map))
  if (length(missing_taxids) > 0) {
    fetched <- fetch_taxid_species_from_ncbi(missing_taxids)
    if (length(fetched) > 0) {
      taxid_cache_env$map[names(fetched)] <- fetched
      save_taxid_species_cache(cache_path)
      message("Resolved ", length(fetched), " BLAST taxids to species names via NCBI taxonomy.")
    }
  }
  vapply(taxid_list, function(ids) {
    vals <- unname(taxid_cache_env$map[ids])
    vals <- unique(na.omit(vals[nchar(vals) > 0]))
    if (length(vals) == 0) {
      return(NA_character_)
    }
    paste(vals, collapse = ";")
  }, character(1))
}

species_from_title_or_taxid <- function(stitle, staxids) {
  title_species <- extract_blast_species(stitle)
  taxid_species <- taxids_to_species(staxids)
  coalesce(title_species, taxid_species)
}

species_from_blast_fields <- function(sscinames, stitle, staxids) {
  blast_species <- clean_blast_label(sscinames)
  blast_species <- ifelse(is.na(blast_species) | nchar(blast_species) == 0, NA_character_, blast_species)
  coalesce(blast_species, species_from_title_or_taxid(stitle, staxids))
}

collapse_species_values <- function(x) {
  vals <- clean_blast_label(x)
  vals <- unique(na.omit(vals[nchar(vals) > 0]))
  if (length(vals) == 0) {
    return(NA_character_)
  }
  paste(vals, collapse=";")
}

first_text_or_na <- function(x) {
  vals <- as.character(x)
  vals <- vals[!is.na(vals) & nchar(vals) > 0 & !toupper(vals) %in% c("NA", "NAN", "NONE")]
  if (length(vals) == 0) {
    return(NA_character_)
  }
  vals[1]
}

text_or_na_vec <- function(x) {
  vals <- as.character(x)
  vals[is.na(vals) | nchar(vals) == 0 | toupper(vals) %in% c("NA", "NAN", "NONE")] <- NA_character_
  vals
}

coalesce_text_cols <- function(...) {
  cols <- list(...)
  if (length(cols) == 0) {
    return(character(0))
  }
  dplyr::coalesce(!!!lapply(cols, text_or_na_vec))
}

extract_feature_qualifier <- function(features, qualifier) {
  pattern <- paste0("['\"]", qualifier, "['\"]:\\s*(?:\\[([^\\]]*)\\]|([^,}\\]]+))")
  matches <- str_match_all(replace_na(as.character(features), ""), pattern)
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
    labels <- labels[!labels %in% c("NO MATCH", "NO PROTEIN/GENE HIT", "UNANNOTATED", "UNCHARACTERISED")]
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
  special_labels <- intersect(c("NO TARGET", "NO MATCH", "UNANNOTATED", "UNCHARACTERISED"), unique(labels))
  real_labels <- labels[!labels %in% c("NO TARGET", "NO MATCH", "UNANNOTATED", "UNCHARACTERISED")]
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
  label <- as.character(label)
  label <- replace_na(label, "")
  label <- str_replace_all(label, "\\s*\\(OTHER TAXA\\)", "")
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

ensure_columns <- function(tbl, cols) {
  tbl <- as_tibble(tbl)
  for (col in cols) {
    if (!col %in% colnames(tbl)) {
      tbl[[col]] <- rep(NA_character_, nrow(tbl))
    }
  }
  tbl
}

normalize_reblast_table <- function(tbl, mode=c("blast", "blastp")) {
  mode <- match.arg(mode)
  common_cols <- c("metadata_category", "query", "stitle", "staxids", "sscinames", "identity", "qcovs",
                   "subject", "sacc", "NCBI_protein_accession", "UniProt_accession")
  blast_cols <- c("features", "features_all", "features_10000_window")
  cols <- if (mode == "blast") c(common_cols, blast_cols) else common_cols
  tbl <- ensure_columns(tbl, cols)
  for (col in cols) {
    tbl[[col]] <- vapply(seq_len(nrow(tbl)), function(i) coerce_scalar_text(tbl[[col]][[i]]), character(1))
  }
  tbl$identity <- suppressWarnings(as.numeric(tbl$identity))
  tbl$qcovs <- suppressWarnings(as.numeric(tbl$qcovs))
  tbl
}

ensure_annotation_columns <- function(tbl) {
  fallback_cols <- c("metadata_category", "feature", "cluster", "query", "accuracy",
                     "classes", "coefficients", "identity", "qcovs", "confounders",
                     "stitle", "staxids", "sscinames", "features", "features_10000_window", "features_all",
                     "compactor_species", "compactor_staxids", "compactor_sscinames",
                     "subject", "sacc", "NCBI_protein_accession", "UniProt_accession",
                     "compactor_blast_subject_id", "compactor_blast_accession",
                     "compactor_ncbi_protein_accession", "compactor_uniprot_accession",
                     "first_class", "first_coef", "second_coef", "second_class",
                     "max_coefficient")
  for (col in fallback_cols) {
    if (!col %in% colnames(tbl)) {
      tbl[[col]] <- NA_character_
    }
  }
  for (col in c("metadata_category", "feature", "cluster", "query", "classes",
                "coefficients", "confounders", "stitle", "staxids", "features",
                "sscinames", "features_10000_window", "features_all", "compactor_species", "compactor_staxids",
                "compactor_sscinames",
                "subject", "sacc", "NCBI_protein_accession", "UniProt_accession",
                "compactor_blast_subject_id", "compactor_blast_accession",
                "compactor_ncbi_protein_accession", "compactor_uniprot_accession")) {
    tbl[[col]] <- vapply(seq_len(nrow(tbl)), function(i) coerce_scalar_text(tbl[[col]][[i]]), character(1))
  }
  tbl$accuracy <- suppressWarnings(as.numeric(tbl$accuracy))
  tbl$identity <- suppressWarnings(as.numeric(tbl$identity))
  tbl$qcovs <- suppressWarnings(as.numeric(tbl$qcovs))
  tbl$first_coef <- suppressWarnings(as.numeric(tbl$first_coef))
  tbl$second_coef <- suppressWarnings(as.numeric(tbl$second_coef))
  tbl$max_coefficient <- suppressWarnings(as.numeric(tbl$max_coefficient))
  tbl
}

with_plot_coefficients <- function(tbl) {
  tbl %>%
    mutate(parsed_first_coef = get_first_coef(coefficients),
           parsed_second_coef = get_nth_coef(coefficients, 2),
           parsed_first_class = get_first_class(classes),
           parsed_second_class = get_nth_class(classes, 2)) %>%
    mutate(first_coef = coalesce(suppressWarnings(as.numeric(first_coef)), parsed_first_coef),
           second_coef = coalesce(suppressWarnings(as.numeric(second_coef)), parsed_second_coef),
           first_class = coalesce(first_class, parsed_first_class),
           second_class = coalesce(second_class, parsed_second_class),
           max_coefficient = coalesce(suppressWarnings(as.numeric(max_coefficient)), abs(first_coef))) %>%
    select(-parsed_first_coef, -parsed_second_coef, -parsed_first_class, -parsed_second_class)
}

has_restricted_label <- function(label) {
  label <- str_squish(as.character(label))
  !is.na(label) & nchar(label) > 1 &
    !str_to_upper(label) %in% c("", "NA", "NAN", "NONE", "NO MATCH", "NO TARGET",
                                "UNANNOTATED", "UNCHARACTERISED", "NO PROTEIN/GENE HIT", "BLAST", "BLASTP")
}

compactor_plot_label <- function(label) {
  label <- clean_blast_label(label)
  ifelse(has_restricted_label(label) & !str_detect(label, "\\s*\\(COMPACTOR\\)\\s*$"),
         paste0(label, " (COMPACTOR)"),
         label)
}

make_compactor_summary_label_dt <- function(path) {
  empty_dt <- tibble(metadata_category=character(), feature=character(), cluster=character(),
                     sequence=character(), compactor_summary_label=character(),
                     compactor_summary_identity=numeric(), compactor_summary_qcovs=numeric(),
                     compactor_summary_species=character(), compactor_summary_staxids=character(),
                     compactor_summary_sequence=character(),
                     compactor_summary_subject_id=character(), compactor_summary_accession=character(),
                     compactor_summary_ncbi_protein_accession=character(),
                     compactor_summary_uniprot_accession=character())
  if (is.null(path) || is.na(path) || nchar(path) == 0 || !file.exists(path) || file.info(path)$size == 0) {
    return(empty_dt)
  }
  compactor_dt <- fread(path)
  compactor_cols <- c("metadata_category", "feature", "cluster", "sequence", "compactor_annotation",
                      "Blast Label", "identity", "qcovs", "compactor_species", "compactor_staxids",
                      "compactor_sscinames", "compactor_sequence",
                      "compactor_blast_subject_id", "compactor_blast_accession",
                      "compactor_ncbi_protein_accession", "compactor_uniprot_accession")
  compactor_dt <- ensure_columns(compactor_dt, compactor_cols)
  for (col in compactor_cols) {
    compactor_dt[[col]] <- vapply(seq_len(nrow(compactor_dt)),
                                  function(i) coerce_scalar_text(compactor_dt[[col]][[i]]),
                                  character(1))
  }
  compactor_dt %>%
    mutate(sequence = str_remove_all(str_remove(as.character(sequence), "^cluster_\\d+_"), "-")) %>%
    mutate(compactor_summary_label = ifelse(has_restricted_label(compactor_annotation),
                                            compactor_annotation, `Blast Label`)) %>%
    mutate(compactor_summary_label = compactor_plot_label(compactor_summary_label)) %>%
    filter(has_restricted_label(compactor_summary_label) |
             (!is.na(compactor_sequence) & nchar(compactor_sequence) > 0 & compactor_sequence != "NA")) %>%
    group_by(metadata_category, feature, cluster, sequence) %>%
    summarise(compactor_summary_label = collapse_blast_labels(paste(unique(na.omit(compactor_summary_label)), collapse=";")),
              compactor_summary_identity = first_numeric_or_na(identity),
              compactor_summary_qcovs = first_numeric_or_na(qcovs),
              compactor_summary_species = collapse_species_values(coalesce(compactor_species, compactor_sscinames)),
              compactor_summary_staxids = first_text_or_na(compactor_staxids),
              compactor_summary_sequence = first_text_or_na(compactor_sequence),
              compactor_summary_subject_id = first_text_or_na(compactor_blast_subject_id),
              compactor_summary_accession = first_text_or_na(compactor_blast_accession),
              compactor_summary_ncbi_protein_accession = first_text_or_na(compactor_ncbi_protein_accession),
              compactor_summary_uniprot_accession = first_text_or_na(compactor_uniprot_accession),
              .groups="drop")
}

make_compactor_selected_dt <- function(path) {
  empty_dt <- tibble(anchor=character(), compactor_selected_sequence=character())
  if (is_missing_path(path)) {
    return(empty_dt)
  }
  selected_dt <- fread(path)
  selected_dt <- ensure_columns(selected_dt, c("anchor", "compactor_sequence", "compactor"))
  selected_dt %>%
    mutate(anchor = str_remove_all(as.character(anchor), "-")) %>%
    mutate(compactor_selected_sequence = coalesce_text_cols(compactor_sequence, compactor)) %>%
    filter(!is.na(anchor), nchar(anchor) > 0,
           !is.na(compactor_selected_sequence), nchar(compactor_selected_sequence) > 0) %>%
    group_by(anchor) %>%
    summarise(compactor_selected_sequence = first_text_or_na(compactor_selected_sequence),
              .groups="drop")
}

make_compactor_seed_dt <- function(path) {
  empty_dt <- tibble(sequence=character(), compactor_seed_sequence=character())
  if (is_missing_path(path)) {
    return(empty_dt)
  }
  seed_dt <- fread(path)
  seed_dt <- ensure_columns(seed_dt, c("seed_extendor", "compactor_sequence", "compactor"))
  seed_dt %>%
    mutate(sequence = str_remove_all(as.character(seed_extendor), "-")) %>%
    mutate(compactor_seed_sequence = coalesce_text_cols(compactor_sequence, compactor)) %>%
    filter(!is.na(sequence), nchar(sequence) > 0,
           !is.na(compactor_seed_sequence), nchar(compactor_seed_sequence) > 0) %>%
    group_by(sequence) %>%
    summarise(compactor_seed_sequence = first_text_or_na(compactor_seed_sequence),
              .groups="drop")
}

make_reblastp_label_dt <- function(reblastp_dt, category) {
  if (nrow(reblastp_dt) == 0 || !"query" %in% colnames(reblastp_dt)) {
    return(tibble(sequence=character(), outside_taxid_label=character(),
                  outside_taxid_identity=numeric(), outside_taxid_qcovs=numeric(),
                  outside_taxid_species=character(), outside_taxid_staxids=character(),
                  outside_taxid_subject_id=character(), outside_taxid_accession=character()))
  }
  reblastp_dt <- ensure_columns(reblastp_dt, c("metadata_category", "query", "stitle", "staxids", "sscinames", "identity", "qcovs",
                                               "subject", "sacc", "NCBI_protein_accession", "UniProt_accession"))
  reblastp_dt %>%
    filter(metadata_category == category) %>%
    mutate(sequence = str_remove(query, "^cluster_\\d+_")) %>%
    mutate(outside_taxid_label = clean_blast_label(str_remove_all(stitle, "\\[.+\\]$|MULTISPECIES:\\s|, partial"))) %>%
    mutate(outside_taxid_species = species_from_blast_fields(sscinames, stitle, staxids)) %>%
    mutate(outside_taxid_subject_id = coalesce(subject, NCBI_protein_accession, sacc)) %>%
    mutate(outside_taxid_accession = coalesce(NCBI_protein_accession, UniProt_accession, sacc, subject)) %>%
    group_by(sequence) %>%
    summarise(outside_taxid_label = collapse_blast_labels(paste(unique(na.omit(outside_taxid_label)), collapse=";")),
              outside_taxid_identity = first_numeric_or_na(identity),
              outside_taxid_qcovs = first_numeric_or_na(qcovs),
              outside_taxid_species = collapse_species_values(outside_taxid_species),
              outside_taxid_staxids = first_text_or_na(staxids),
              outside_taxid_subject_id = first_text_or_na(outside_taxid_subject_id),
              outside_taxid_accession = first_text_or_na(outside_taxid_accession),
              .groups="drop")
}

make_reblast_label_dt <- function(reblast_dt, category) {
  if (nrow(reblast_dt) == 0 || !"query" %in% colnames(reblast_dt)) {
    return(tibble(sequence=character(), outside_taxid_label=character(),
                  outside_taxid_identity=numeric(), outside_taxid_qcovs=numeric(),
                  outside_taxid_species=character(), outside_taxid_staxids=character(),
                  outside_taxid_subject_id=character(), outside_taxid_accession=character()))
  }
  reblast_dt <- ensure_columns(reblast_dt, c("metadata_category", "query", "stitle", "staxids", "sscinames", "identity", "qcovs",
                                             "subject", "sacc", "NCBI_protein_accession", "UniProt_accession",
                                             "features", "features_all", "features_10000_window"))
  reblast_dt %>%
    filter(metadata_category == category) %>%
    mutate(sequence = str_remove(query, "^cluster_\\d+_")) %>%
    mutate(feature_text = paste(replace_na(as.character(features), ""),
                                replace_na(as.character(features_all), ""),
                                sep=";")) %>%
    mutate(outside_products = extract_feature_qualifier(feature_text, "product")) %>%
    mutate(outside_genes = extract_feature_qualifier(feature_text, "gene")) %>%
    mutate(outside_taxid_label = choose_feature_label(outside_products, outside_genes, opt$products)) %>%
    mutate(outside_taxid_species = species_from_blast_fields(sscinames, stitle, staxids)) %>%
    mutate(outside_taxid_subject_id = coalesce(subject, sacc, NCBI_protein_accession)) %>%
    mutate(outside_taxid_accession = coalesce(sacc, NCBI_protein_accession, UniProt_accession, subject)) %>%
    group_by(sequence) %>%
    summarise(outside_taxid_label = collapse_blast_labels(paste(unique(na.omit(outside_taxid_label)), collapse=";")),
              outside_taxid_identity = first_numeric_or_na(identity),
              outside_taxid_qcovs = first_numeric_or_na(qcovs),
              outside_taxid_species = collapse_species_values(outside_taxid_species),
              outside_taxid_staxids = first_text_or_na(staxids),
              outside_taxid_subject_id = first_text_or_na(outside_taxid_subject_id),
              outside_taxid_accession = first_text_or_na(outside_taxid_accession),
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
dt <- ensure_annotation_columns(dt)
dt2 <- ensure_annotation_columns(dt2)
for (compactor_col in c("compactor_annotation", "compactor_query", "compactor_sequence", "compactor_length",
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
message(paste0("Compactor-labeled annotation rows seen by plotter: blastp=",
               sum(has_restricted_label(dt$compactor_annotation), na.rm=TRUE),
               ", blast=",
               sum(has_restricted_label(dt2$compactor_annotation), na.rm=TRUE)))
compactor_summary_label_dt <- make_compactor_summary_label_dt(opt$compactor_summary)
if (nrow(compactor_summary_label_dt) > 0) {
  message(paste0("Compactor labels loaded from plot summary lookup: ",
                 nrow(compactor_summary_label_dt)))
}
compactor_selected_dt <- make_compactor_selected_dt(opt$compactor_selected)
compactor_selected_anchor_len <- if (nrow(compactor_selected_dt) > 0) {
  max(nchar(compactor_selected_dt$anchor), na.rm=TRUE)
} else {
  NA_integer_
}
if (nrow(compactor_selected_dt) > 0) {
  message(paste0("Selected compactor sequences loaded for anchor lookup: ",
                 nrow(compactor_selected_dt), " anchors of length ",
                 compactor_selected_anchor_len))
} else {
  message("Selected compactor sequence lookup loaded 0 usable rows.")
}
compactor_seed_dt <- make_compactor_seed_dt(opt$compactor_seed_annotations)
if (nrow(compactor_seed_dt) > 0) {
  message(paste0("Seed compactor sequences loaded for exact extendor lookup: ",
                 nrow(compactor_seed_dt), " plotted extendors"))
} else {
  message("Seed compactor sequence lookup loaded 0 usable rows.")
}
dt_reblastp <- read_optional_tsv(opt$reblastp_annotations)
dt_reblast <- read_optional_tsv(opt$reblast_annotations)
dt_reblastp <- normalize_reblast_table(dt_reblastp, "blastp")
dt_reblast <- normalize_reblast_table(dt_reblast, "blast")
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
      with_plot_coefficients() %>%
      arrange(-max_coefficient) %>%
      mutate(annotation = str_remove_all(stitle, "\\[.+\\]$|MULTISPECIES:\\s|, partial")) %>%
      mutate(annotation = ifelse(has_restricted_label(compactor_annotation),
                                 compactor_plot_label(compactor_annotation), annotation)) %>%
      rowwise() %>%
      mutate(classes=list(str_split_1(gsub("\\[|\\]", "", classes),pattern=","))) %>%
      ungroup() %>%
      select(metadata_category, accuracy, classes, first_class, first_coef, second_coef, second_class, max_coefficient,
             cluster, feature, query, identity, qcovs, stitle, staxids, sscinames,
             subject, sacc, NCBI_protein_accession, UniProt_accession,
             annotation, confounders, compactor_sequence) %>%
      mutate(query = str_remove(query, "cluster_\\d+_")) %>%
      distinct(cluster,annotation,query,.keep_all = T) %>% group_by(cluster)


    summ_dt <- summ_dt %>% group_by(cluster,query) %>%
      mutate(label=ifelse(!is_empty(unique(na.omit(annotation))), paste0(unique(na.omit(annotation)),collapse=";"), NA)) %>%
      mutate(blastp_species=collapse_species_values(species_from_blast_fields(sscinames, stitle, staxids))) %>%
      mutate(blastp_staxids=first_text_or_na(staxids)) %>%
      mutate(blastp_subject_id=first_text_or_na(coalesce(subject, NCBI_protein_accession, sacc))) %>%
      mutate(blastp_accession=first_text_or_na(coalesce(NCBI_protein_accession, UniProt_accession, sacc, subject))) %>%
      mutate(blastp_compactor_sequence=first_text_or_na(compactor_sequence)) %>%
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
        select(-feature_text) %>% with_plot_coefficients() %>%
        arrange(-max_coefficient) %>%
        rowwise() %>%
        mutate(classes=list(str_split_1(gsub("\\[|\\]", "", classes),pattern=","))) %>%
        ungroup() %>%
        select(metadata_category, accuracy, classes, first_class, first_coef, max_coefficient,
               cluster, feature, query, identity, qcovs, stitle, staxids, sscinames,
               subject, sacc, NCBI_protein_accession, UniProt_accession,
               products, genes, confounders, compactor_annotation, compactor_species, compactor_staxids,
               compactor_sscinames, compactor_sequence,
               compactor_blast_subject_id, compactor_blast_accession,
               compactor_ncbi_protein_accession, compactor_uniprot_accession) %>%
        mutate(query = str_remove(query, "^cluster_\\d+_")) %>%
        group_by(cluster) %>%
        ungroup() %>%
        distinct(cluster,products,query,genes,.keep_all = T) %>% group_by(cluster)

      summ_dt2 <- summ_dt2 %>%
        mutate(label = choose_feature_label(products, genes, opt$products)) %>%
        mutate(blast_species = coalesce(species_from_blast_fields(sscinames, stitle, staxids), compactor_species, compactor_sscinames)) %>%
        mutate(blast_staxids = coalesce(staxids, compactor_staxids)) %>%
        mutate(blast_subject_id = coalesce(subject, sacc, NCBI_protein_accession,
                                           compactor_blast_subject_id,
                                           compactor_ncbi_protein_accession,
                                           compactor_blast_accession)) %>%
        mutate(blast_accession = coalesce(sacc, NCBI_protein_accession, UniProt_accession,
                                          compactor_blast_accession,
                                          compactor_ncbi_protein_accession,
                                          compactor_uniprot_accession,
                                          subject)) %>%
        mutate(label = ifelse(has_restricted_label(compactor_annotation),
                              compactor_plot_label(compactor_annotation), label)) %>%
        group_by(cluster,query) %>%
        mutate(label=ifelse(!is_empty(unique(na.omit(label))), paste0(unique(na.omit(label)),collapse=";"), NA)) %>%
        mutate(blast_species=collapse_species_values(blast_species)) %>%
        mutate(blast_staxids=first_text_or_na(blast_staxids)) %>%
        mutate(blast_subject_id=first_text_or_na(blast_subject_id)) %>%
        mutate(blast_accession=first_text_or_na(blast_accession)) %>%
        mutate(compactor_sequence = first_text_or_na(compactor_sequence)) %>%
        distinct(cluster, query, label, .keep_all=T) %>% ungroup()
      summ_dt2 <- summ_dt2 %>% select(cluster, query, identity, qcovs, label, blast_species, blast_staxids,
                                      blast_subject_id, blast_accession, compactor_sequence) %>%
        dplyr::rename(label2=label, blast_species2=blast_species, blast_staxids2=blast_staxids,
                      blast_subject_id2=blast_subject_id, blast_accession2=blast_accession,
                      blastn_compactor_sequence=compactor_sequence)
      summ_dt <- summ_dt %>% left_join(summ_dt2 %>%
                                         select(cluster, query, identity, qcovs, label2, blast_species2, blast_staxids2,
                                                blast_subject_id2, blast_accession2, blastn_compactor_sequence), by=c("cluster", "query")) %>%
        dplyr::rename(identity=`identity.x`, qcovs=`qcovs.x`) %>%
        mutate(identity = ifelse(is.na(label) & !is.na(label2), `identity.y`, identity)) %>%
        mutate(qcovs = ifelse(is.na(label) & !is.na(label2), `qcovs.y`, qcovs)) %>%
        mutate(direct_blast_species = coalesce(blastp_species, blast_species2),
               direct_blast_staxids = coalesce(blastp_staxids, blast_staxids2),
               direct_blast_subject_id = coalesce(blastp_subject_id, blast_subject_id2),
               direct_blast_accession = coalesce(blastp_accession, blast_accession2)) %>%
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
                    outside_taxid_qcovs_blastp = outside_taxid_qcovs,
                    outside_taxid_species_blastp = outside_taxid_species,
                    outside_taxid_staxids_blastp = outside_taxid_staxids,
                    outside_taxid_subject_id_blastp = outside_taxid_subject_id,
                    outside_taxid_accession_blastp = outside_taxid_accession)
    reblast_label_dt <- make_reblast_label_dt(dt_reblast, category) %>%
      dplyr::rename(outside_taxid_label_blast = outside_taxid_label,
                    outside_taxid_identity_blast = outside_taxid_identity,
                    outside_taxid_qcovs_blast = outside_taxid_qcovs,
                    outside_taxid_species_blast = outside_taxid_species,
                    outside_taxid_staxids_blast = outside_taxid_staxids,
                    outside_taxid_subject_id_blast = outside_taxid_subject_id,
                    outside_taxid_accession_blast = outside_taxid_accession)
    outside_taxid_label_dt <- full_join(reblastp_label_dt, reblast_label_dt, by="sequence") %>%
      mutate(outside_taxid_label = combine_blast_labels(outside_taxid_label_blastp, outside_taxid_label_blast)) %>%
      mutate(outside_taxid_label = ifelse(is.na(outside_taxid_label) & !is.na(outside_taxid_label_blast), outside_taxid_label_blast, outside_taxid_label)) %>%
      mutate(outside_taxid_label = ifelse(is.na(outside_taxid_label) & !is.na(outside_taxid_label_blastp), outside_taxid_label_blastp, outside_taxid_label)) %>%
      mutate(outside_taxid_identity = coalesce(outside_taxid_identity_blastp, outside_taxid_identity_blast),
             outside_taxid_qcovs = coalesce(outside_taxid_qcovs_blastp, outside_taxid_qcovs_blast),
             outside_taxid_species = coalesce(outside_taxid_species_blastp, outside_taxid_species_blast),
             outside_taxid_staxids = coalesce(outside_taxid_staxids_blastp, outside_taxid_staxids_blast),
             outside_taxid_subject_id = coalesce(outside_taxid_subject_id_blastp, outside_taxid_subject_id_blast),
             outside_taxid_accession = coalesce(outside_taxid_accession_blastp, outside_taxid_accession_blast)) %>%
      mutate(outside_taxid_species = coalesce(outside_taxid_species,
                                              taxids_to_species(outside_taxid_staxids))) %>%
      select(sequence, outside_taxid_label, outside_taxid_identity, outside_taxid_qcovs,
             outside_taxid_species, outside_taxid_staxids, outside_taxid_subject_id, outside_taxid_accession)
    category_compactor_summary_dt <- compactor_summary_label_dt %>%
      filter(metadata_category == category)

    blastp_all_dt <- summ_dt_blastp_only %>%
      group_by(feature) %>%
      mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)) & (is.na(label)) & !is.na(identity), "NO PROTEIN/GENE HIT", label)) %>%
      group_by(cluster) %>%
      mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)), "NO MATCH", label)) %>%
      ungroup() %>%
      mutate(blast_source = "blastp") %>%
      mutate(classes = map_chr(classes, \(x) paste(x, collapse=","))) %>%
      select(any_of(c("metadata_category", "accuracy", "classes", "first_class", "first_coef", "second_coef", "second_class",
                      "max_coefficient", "cluster", "feature", "query", "identity", "qcovs",
                      "blastp_species", "blastp_staxids", "blastp_subject_id", "blastp_accession",
                      "label", "blast_source")))

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
                      "max_coefficient", "cluster", "feature", "query", "identity.y", "qcovs.y",
                      "direct_blast_species", "direct_blast_staxids", "direct_blast_subject_id",
                      "direct_blast_accession", "label", "blast_source"))) %>%
      dplyr::rename(identity = `identity.y`, qcovs = `qcovs.y`,
                    blast_species = direct_blast_species,
                    blast_staxids = direct_blast_staxids,
                    blast_subject_id = direct_blast_subject_id,
                    blast_accession = direct_blast_accession)

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
                                 compactor_plot_label(compactor_annotation), hist_label)) %>%
      mutate(hist_label = case_when(
        str_detect(query, "NNNNNNNN") ~ "NO TARGET",
        !is.na(hist_label) & nchar(hist_label) > 1 ~ hist_label,
        !is.na(identity) | !is.na(qcovs) ~ "UNANNOTATED",
        TRUE ~ NA_character_
      )) %>%
      with_plot_coefficients() %>%
      select(cluster, feature, max_coefficient, label=hist_label)
    hist_compactor_label_dt <- category_compactor_summary_dt %>%
      left_join(summ_dt %>% distinct(cluster, feature, max_coefficient),
                by=c("cluster", "feature")) %>%
      filter(!is.na(max_coefficient)) %>%
      select(cluster, feature, max_coefficient, label=compactor_summary_label)

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
      hist_direct_label_dt,
      hist_compactor_label_dt
    )

    largest_coef <- suppressWarnings(max(hist_label_dt$max_coefficient, na.rm=TRUE))
    if (!is.finite(largest_coef) || largest_coef <= 0) {
      largest_coef <- 1
    }
    plot_dt <- hist_label_dt %>%
      filter(!is.na(max_coefficient)) %>%
      mutate(coef_mag=max_coefficient/largest_coef) %>%
      group_by(cluster, feature, coef_mag) %>%
      summarise(label=make_histogram_label(label), .groups="drop") %>%
      mutate(label = str_replace(label, " ,", ", ") %>% str_replace(" ;", "; ")) %>%
      arrange(-coef_mag) %>%
      mutate(rank=row_number()) %>%
      mutate(color="no_taxon") %>%
      mutate(has_real_label = !str_detect(label, "^(NO TARGET|NO MATCH|UNANNOTATED|UNCHARACTERISED)(,\\s*(NO TARGET|NO MATCH|UNANNOTATED|UNCHARACTERISED))*$")) %>%
      mutate(color=ifelse(has_real_label, "within_taxid", color)) %>%
      mutate(color=ifelse(str_detect(label, "\\(OTHER TAXA\\)"), "outside_taxid", color)) %>%
      mutate(label_size = case_when(
        str_detect(label, "^(NO TARGET|NO MATCH|UNANNOTATED|UNCHARACTERISED)$") ~ 3.35,
        nchar(label) <= 55 ~ 2.85,
        nchar(label) <= 95 ~ 2.55,
        TRUE ~ 2.25
      )) %>%
      mutate(label = str_wrap(preserve_compactor_suffix(label, width=120), width = 38)) %>%
      mutate(label = replace_na(label, ""))
    message(paste0("Histogram rows for ", category, ": ",
                   nrow(plot_dt), " finite rows from ",
                   nrow(hist_label_dt), " labels."))
    compactor_hist_count <- sum(str_detect(plot_dt$label, "\\(COMPACTOR\\)"), na.rm=TRUE)
    if (compactor_hist_count > 0) {
      message(paste0("Compactor labels in histogram for ", category, ": ", compactor_hist_count))
    }

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

    plot_dt_top <- plot_dt %>% head(opt$num_hits)
    full_cluster_extendors <- map_dfr(seq_len(nrow(plot_dt_top)), function(i) {
      cluster_name <- plot_dt_top$cluster[[i]]
      feature_name <- plot_dt_top$feature[[i]]
      rank_value <- plot_dt_top$rank[[i]]
      cluster_idx <- suppressWarnings(as.numeric(str_extract(cluster_name, "\\d+")))
      if (!is.finite(cluster_idx)) {
        return(tibble(cluster=character(), feature=character(), rank=integer(), sequence=character()))
      }
      seq_dt <- read_nth_cluster(opt$sample_seqs, cluster_idx, opt$cluster_length)
      if (nrow(seq_dt) == 0 || !"sequence" %in% colnames(seq_dt)) {
        return(tibble(cluster=character(), feature=character(), rank=integer(), sequence=character()))
      }
      seq_dt %>%
        transmute(cluster=cluster_name,
                  feature=feature_name,
                  rank=rank_value,
                  sequence=as.character(sequence)) %>%
        filter(!is.na(sequence), nchar(sequence) > 0) %>%
        distinct(cluster, feature, rank, sequence)
    })
    extendor_annotation_dt <- summ_dt %>%
      mutate(sequence = query) %>%
      left_join(outside_taxid_label_dt, by="sequence") %>%
      left_join(category_compactor_summary_dt %>%
                  select(cluster, feature, sequence, compactor_summary_label),
                by=c("cluster", "feature", "sequence")) %>%
      mutate(has_within_taxid_signal =
               has_restricted_label(label) |
               has_restricted_label(label_blastp) |
               has_restricted_label(label_blast) |
               has_restricted_label(annotation) |
               has_restricted_label(compactor_summary_label) |
               !is.na(identity) | !is.na(qcovs) |
               !is.na(`identity.y`) | !is.na(`qcovs.y`),
             has_outside_taxid_signal = has_restricted_label(outside_taxid_label)) %>%
      group_by(cluster, feature, sequence) %>%
      summarise(has_within_taxid_signal=any(has_within_taxid_signal, na.rm=TRUE),
                has_outside_taxid_signal=any(has_outside_taxid_signal, na.rm=TRUE),
                .groups="drop")
    compactor_annotation_dt <- category_compactor_summary_dt %>%
      mutate(has_compactor_annotation=has_restricted_label(compactor_summary_label)) %>%
      group_by(cluster, feature, sequence) %>%
      summarise(has_compactor_annotation=any(has_compactor_annotation, na.rm=TRUE),
                .groups="drop")
    plot_stack_top <- full_cluster_extendors %>%
      left_join(extendor_annotation_dt, by=c("cluster", "feature", "sequence")) %>%
      left_join(compactor_annotation_dt, by=c("cluster", "feature", "sequence")) %>%
      mutate(across(c(has_within_taxid_signal, has_outside_taxid_signal,
                      has_compactor_annotation), \(x) replace_na(x, FALSE))) %>%
      mutate(taxon_source = case_when(
        has_within_taxid_signal | has_compactor_annotation ~ "within_taxid",
        has_outside_taxid_signal ~ "outside_taxid",
        TRUE ~ "no_taxon"
      )) %>%
      count(cluster, feature, rank, taxon_source, name="extendor_n")
    message(paste0("Stacked extendor-count histogram rows for ", category, ": ",
                   sum(plot_stack_top$extendor_n), " unique extendors across ",
                   n_distinct(paste(plot_stack_top$cluster, plot_stack_top$feature)), " plotted features."))
    p <- ggplot(plot_dt_top, aes(x=rank, y=coef_mag)) +
      geom_col(fill=histogram_bar_color) +
      geom_text(data=plot_dt_top,
                aes(x=rank, y=coef_mag + 0.05, label=label, hjust=0, size=label_size),
                angle=45, color="#222222", show.legend=FALSE) +
      scale_size_identity() +
      scale_y_continuous("Magnitude relative to\nlargest nonzero coefficient",
                         limits = c(0,1.6),
                         labels=scales::label_percent(), breaks=seq(0,1,0.25),
                         expand=c(0,0)) +
      xlab("Rank of nonzero coefficient (by magnitude)") +
      scale_x_continuous(breaks=seq(1,10,1)) +
      ggtitle(make_title,
              subtitle = hist_subtitle) +
      theme_pubr() +
      theme(plot.subtitle = element_text(size=8, lineheight=0.95))

    print(p)

    p_stack <- ggplot(plot_stack_top,
                      aes(x=rank, y=extendor_n, fill=taxon_source)) +
      geom_col(position="stack") +
      scale_y_continuous("Number of extendors",
                         breaks=scales::breaks_pretty(n=6),
                         expand=expansion(mult=c(0, 0.08))) +
      xlab("Rank of nonzero coefficient (same order as previous histogram)") +
      scale_x_continuous(breaks=seq(1,10,1)) +
      scale_fill_manual(breaks=c("within_taxid","outside_taxid","no_taxon"),
                        values=c("within_taxid"=histogram_bar_color,
                                 "outside_taxid"=outside_taxid_label_color,
                                 "no_taxon"=no_taxon_label_color),
                        labels=c("Within requested taxids", "Outside requested taxids",
                                 "No hit"),
                        name="Taxon source") +
      ggtitle(make_title,
              subtitle = hist_subtitle) +
      theme_pubr() +
      theme(legend.position="right",
            plot.subtitle = element_text(size=8, lineheight=0.95))

    print(p_stack)

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

      first_beta <- suppressWarnings(as.numeric(unique(dt_sub$first_coef)))
      first_beta <- first_beta[is.finite(first_beta)]
      if (length(first_beta) == 0) {
        message(paste("Skipping detailed plot for", category, my_cluster, my_feature,
                      "because no finite coefficient was available."))
        next
      }
      first_beta <- first_beta[which.max(abs(first_beta))]
      first_class = unique(na.omit(dt_sub$first_class))
      if (length(first_class) == 0) {
        first_class <- NA_character_
      } else {
        first_class <- first_class[1]
      }
      all_classes = dt_sub[1,]$classes %>% unlist()
      classes_to_plot <- all_classes
      if (length(all_classes) == 1 && all_classes[1] == "residual" && !is_quantitative_target) {
        all_classes <- sort(unique(na.omit(as.character(my_metadata$metadata))))
        if (length(all_classes) > 0) {
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
        mutate(direct_blast_species = coalesce(species_from_blast_fields(sscinames, stitle, staxids), compactor_species, compactor_sscinames)) %>%
        mutate(direct_blast_staxids = coalesce(staxids, compactor_staxids)) %>%
        mutate(direct_blast_subject_id = coalesce(subject, sacc, NCBI_protein_accession,
                                                  compactor_blast_subject_id,
                                                  compactor_ncbi_protein_accession,
                                                  compactor_blast_accession)) %>%
        mutate(direct_blast_accession = coalesce(sacc, NCBI_protein_accession, UniProt_accession,
                                                 compactor_blast_accession,
                                                 compactor_ncbi_protein_accession,
                                                 compactor_uniprot_accession,
                                                 subject)) %>%
        mutate(direct_blast_label = ifelse(has_restricted_label(compactor_annotation),
                                           compactor_plot_label(compactor_annotation), direct_blast_label)) %>%
        group_by(sequence) %>%
        summarise(direct_blast_label = collapse_blast_labels(paste(unique(na.omit(direct_blast_label)), collapse=";")),
                  direct_identity = first_numeric_or_na(identity),
                  direct_qcovs = first_numeric_or_na(qcovs),
                  detail_direct_blast_species = collapse_species_values(direct_blast_species),
                  detail_direct_blast_staxids = first_text_or_na(direct_blast_staxids),
                  detail_direct_blast_subject_id = first_text_or_na(direct_blast_subject_id),
                  detail_direct_blast_accession = first_text_or_na(direct_blast_accession),
                  .groups="drop")
      compactor_detail_label_dt <- category_compactor_summary_dt %>%
        filter(cluster == my_cluster, feature == my_feature) %>%
        select(sequence, compactor_summary_label,
               compactor_summary_identity, compactor_summary_qcovs,
               compactor_summary_species, compactor_summary_staxids,
               compactor_summary_sequence,
               compactor_summary_subject_id, compactor_summary_accession,
               compactor_summary_ncbi_protein_accession,
               compactor_summary_uniprot_accession)

      label_dt <- dt_sub %>%
        mutate(any_identity = coalesce(`identity.y`, identity),
               any_qcovs = coalesce(`qcovs.y`, qcovs)) %>%
        mutate(combined_label = combine_blast_labels(label_blastp, label_blast)) %>%
        mutate(combined_label = ifelse(is.na(combined_label) & !is.na(label_blast), label_blast, combined_label)) %>%
        mutate(combined_label = ifelse(is.na(combined_label) & !is.na(label_blastp), label_blastp, combined_label)) %>%
        select(query, accuracy, any_identity, any_qcovs, combined_label, label2,
               direct_blast_species, direct_blast_staxids,
               direct_blast_subject_id, direct_blast_accession,
               any_of(c("blastp_compactor_sequence", "blastn_compactor_sequence"))) %>%
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
        left_join(label_dt, by="sequence") %>%
        left_join(direct_blast_label_dt, by="sequence") %>%
        left_join(compactor_detail_label_dt, by="sequence") %>%
        left_join(outside_taxid_label_dt, by="sequence") %>%
        left_join(compactor_seed_dt, by="sequence")
      if (nrow(compactor_selected_dt) > 0 && !is.na(compactor_selected_anchor_len)) {
        summ_sub_dt <- summ_sub_dt %>%
          mutate(compactor_selected_anchor = str_sub(as.character(sequence), -compactor_selected_anchor_len)) %>%
          left_join(compactor_selected_dt, by=c("compactor_selected_anchor"="anchor")) %>%
          select(-compactor_selected_anchor)
      }
      summ_sub_dt <- ensure_columns(summ_sub_dt,
                                    c("compactor_summary_sequence",
                                      "compactor_seed_sequence",
                                      "blastp_compactor_sequence",
                                      "blastn_compactor_sequence",
                                      "compactor_selected_sequence"))

      p_sub <- summ_sub_dt %>%
        mutate(label = ifelse((is.na(label) | label == "NO MATCH") &
                                 !is.na(direct_blast_label) & nchar(direct_blast_label) > 1,
                              direct_blast_label, label)) %>%
        mutate(identity = coalesce(identity, direct_identity),
               qcovs = coalesce(qcovs, direct_qcovs),
               direct_blast_species = coalesce(direct_blast_species, detail_direct_blast_species),
               direct_blast_staxids = coalesce(direct_blast_staxids, detail_direct_blast_staxids),
               direct_blast_subject_id = coalesce(direct_blast_subject_id, detail_direct_blast_subject_id),
               direct_blast_accession = coalesce(direct_blast_accession, detail_direct_blast_accession)) %>%
        mutate(label = ifelse(!has_restricted_label(label) &
                                has_restricted_label(compactor_summary_label),
                              compactor_summary_label, label)) %>%
        mutate(compactor_sequence = coalesce_text_cols(compactor_summary_sequence,
                                                       compactor_seed_sequence,
                                                       blastp_compactor_sequence,
                                                       blastn_compactor_sequence,
                                                       compactor_selected_sequence)) %>%
        mutate(identity = ifelse(has_restricted_label(compactor_summary_label) & is.na(identity),
                                 compactor_summary_identity, identity),
               qcovs = ifelse(has_restricted_label(compactor_summary_label) & is.na(qcovs),
                              compactor_summary_qcovs, qcovs),
               direct_blast_species = ifelse(has_restricted_label(compactor_summary_label) &
                                               (is.na(direct_blast_species) | nchar(direct_blast_species) == 0),
                                             compactor_summary_species, direct_blast_species),
               direct_blast_staxids = ifelse(has_restricted_label(compactor_summary_label) &
                                               (is.na(direct_blast_staxids) | nchar(direct_blast_staxids) == 0),
                                             compactor_summary_staxids, direct_blast_staxids),
               direct_blast_subject_id = ifelse(has_restricted_label(compactor_summary_label) &
                                                  (is.na(direct_blast_subject_id) | nchar(direct_blast_subject_id) == 0),
                                                coalesce(compactor_summary_subject_id,
                                                         compactor_summary_ncbi_protein_accession,
                                                         compactor_summary_accession),
                                                direct_blast_subject_id),
               direct_blast_accession = ifelse(has_restricted_label(compactor_summary_label) &
                                                 (is.na(direct_blast_accession) | nchar(direct_blast_accession) == 0),
                                               coalesce(compactor_summary_accession,
                                                        compactor_summary_ncbi_protein_accession,
                                                        compactor_summary_uniprot_accession,
                                                        compactor_summary_subject_id),
                                               direct_blast_accession)) %>%
        mutate(outside_taxid_only = !has_restricted_label(label) &
                 !is.na(outside_taxid_label) & nchar(outside_taxid_label) > 1) %>%
        mutate(label = ifelse(outside_taxid_only, outside_taxid_label, label)) %>%
        mutate(identity = ifelse(outside_taxid_only & is.na(identity), outside_taxid_identity, identity),
               qcovs = ifelse(outside_taxid_only & is.na(qcovs), outside_taxid_qcovs, qcovs),
               direct_blast_species = ifelse(outside_taxid_only, outside_taxid_species, direct_blast_species),
               direct_blast_staxids = ifelse(outside_taxid_only, outside_taxid_staxids, direct_blast_staxids),
               direct_blast_subject_id = ifelse(outside_taxid_only, outside_taxid_subject_id, direct_blast_subject_id),
               direct_blast_accession = ifelse(outside_taxid_only, outside_taxid_accession, direct_blast_accession)) %>%
        mutate(blast_species_origin = direct_blast_species,
               blast_staxids_origin = direct_blast_staxids,
               blast_subject_id_origin = direct_blast_subject_id,
               blast_accession_origin = direct_blast_accession,
               blast_annotation_source_id = coalesce(direct_blast_accession,
                                                     direct_blast_subject_id,
                                                     direct_blast_staxids)) %>%
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
        mutate(label_identity = ifelse(is.na(identity), "-", paste0(round(identity,2), "%"))) %>%
        mutate(label_coverage = ifelse(is.na(qcovs), "-", paste0(round(qcovs,2), "%"))) %>%
        mutate(label_identity = replace_na(label_identity, "-")) %>%
        mutate(label_coverage = replace_na(label_coverage, "-")) %>%
        rowwise() %>%
        mutate(label_quality = ifelse(label_coverage != "-" | label_identity != "-",
                                      paste0("I:", str_replace(label_identity, "-", "NA"),
                                             "; C:", str_replace(label_coverage, "-", "NA")), "")) %>%
        mutate(point_label = case_when(
          `Blast Label` %in% c("NO MATCH", "NO TARGET") ~ as.character(`Blast Label`),
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
        mutate(point_label = ifelse(!`Blast Label` %in% c("NO MATCH", "NO TARGET") &
                                      nchar(label_quality) > 0 &
                                      label_quality != "I:100%; C:100%",
                                    paste(point_label, label_quality, sep="\n"),
                                    point_label)) %>%
        mutate(point_label_expr = ifelse(outside_taxid_only,
                                         make_plotmath_other_taxa_label(point_label),
                                         NA_character_)) %>%
        mutate(point_label_color = case_when(
          outside_taxid_only ~ "outside_taxid",
          `Blast Label` %in% c("NO MATCH", "NO TARGET") ~ "no_taxon",
          TRUE ~ "within_taxid"
        )) %>%
        ungroup()
      if (!"total_samples" %in% colnames(p_sub)) {
        missing_class_cols <- setdiff(all_classes, colnames(p_sub))
        for (missing_class_col in missing_class_cols) {
          p_sub[[missing_class_col]] <- 0
        }
        p_sub <- p_sub %>%
          mutate(across(all_of(all_classes), \(x) replace_na(x, 0))) %>%
          mutate(total_samples = rowSums(across(all_of(all_classes))))
      }
      message(paste0("Detail rows for ", category, " ", my_cluster, " ", my_feature,
                     ": rows=", nrow(p_sub),
                     ", finite_embedding=", sum(is.finite(p_sub$embedding), na.rm=TRUE),
                     ", finite_lev_dist=", sum(!is.na(p_sub$lev_dist), na.rm=TRUE),
                     ", finite_total_samples=", sum(is.finite(p_sub$total_samples), na.rm=TRUE)))
      p_sub <- p_sub %>%
        filter(is.finite(embedding), !is.na(lev_dist), is.finite(total_samples), total_samples > 0)
      if (nrow(p_sub) == 0) {
        message(paste("Skipping detailed plot for", category, my_cluster, my_feature,
                      "because no finite plotting rows remained after filtering."))
        next
      }
      compactor_point_count <- sum(str_detect(p_sub$point_label, "\\(COMPACTOR\\)"), na.rm=TRUE)
      if (compactor_point_count > 0) {
        message(paste0("Compactor labels in detail plot for ", category, " ",
                       my_cluster, " ", my_feature, ": ", compactor_point_count))
      }

      if (is_quantitative_target) {
        p_sub <- p_sub %>% mutate(color_value = mean_metadata)
        p2 <- p_sub %>%
          ggplot(aes(x=embedding, y=lev_dist, fill=color_value, size=total_samples,
                     label=point_label)) +
          geom_vline(xintercept = 0, lty="dashed") +
          geom_point(shape=21, color="grey25", stroke=0.35) +
          scale_y_continuous(breaks=scales::breaks_width(1), minor_breaks=NULL) +
          scale_size_continuous(trans = "log", name = "Total Samples",
                                breaks = c(1, 10, 100, 1000, 10000),
                                limits = c(1, 10000), labels = scales::label_log()) +
          ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color != "outside_taxid"),
                                   aes(color=point_label_color),
                                   size=3.2, max.overlaps=Inf,
                                   min.segment.length=0, segment.color="grey45",
                                   segment.size=0.25, box.padding=0.75,
                                   point.padding=0.8, force=6, force_pull=0.08) +
          ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color == "outside_taxid"),
                                   aes(label=point_label_expr, color=point_label_color),
                                   size=3.2, max.overlaps=Inf,
                                   min.segment.length=0, segment.color="grey45",
                                   segment.size=0.25, box.padding=0.75,
                                   point.padding=0.8, force=6, force_pull=0.08, parse=TRUE) +
          scale_fill_gradient(paste0("Mean\n", metadata_source_col),
                              low = "blue", high = "red") +
          scale_color_manual(breaks=c("within_taxid","outside_taxid","no_taxon"),
                             values=c("within_taxid"=histogram_bar_color,
                                      "outside_taxid"=outside_taxid_label_color,
                                      "no_taxon"=no_taxon_label_color),
                             labels=c("Within requested taxids", "Outside requested taxids",
                                      "No hit"),
                             name="Taxon source") +
          guides(fill=guide_colorbar(order=1),
                 size=guide_legend(order=2,
                                   override.aes=list(shape=16, color="black",
                                                     fill="black", alpha=1)),
                 color=guide_legend(order=3,
                                    override.aes=list(label="A", size=4))) +
          theme_minimal() + xlab(expression("Embedding" ~ "\u00D7" ~ beta)) +
          ylab("Levenshtein Distance\n(to most abundant anchor-target)") +
          ggtitle(paste(category, my_cluster, sep=" | "),
                  subtitle=paste("Color: mean observed", metadata_source_col)) +
          theme(panel.grid.minor.y = element_blank())
        print(p2)
        summ_out_dt <- p_sub %>% select(-any_of("label2")) %>%
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
            ggplot(aes(x=embedding, y=lev_dist, fill=class_proportion,
                       size=total_samples, label=point_label)) +
            geom_vline(xintercept = 0, lty="dashed") +
            geom_point(shape=21, color="grey25", stroke=0.35) +
            scale_y_continuous(breaks=scales::breaks_width(1), minor_breaks=NULL) +
            scale_size_continuous(trans = "log", name = "Total Samples",
                                  breaks = c(1, 10, 100, 1000, 10000),
                                  limits = c(1, 10000), labels = scales::label_log()) +
            ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color != "outside_taxid"),
                                     aes(color=point_label_color),
                                     size=3.2, max.overlaps=Inf,
                                     min.segment.length=0, segment.color="grey45",
                                     segment.size=0.25, box.padding=0.75,
                                     point.padding=0.8, force=6, force_pull=0.08) +
            ggrepel::geom_text_repel(data = \(x) filter(x, point_label_color == "outside_taxid"),
                                     aes(label=point_label_expr, color=point_label_color),
                                     size=3.2, max.overlaps=Inf,
                                     min.segment.length=0, segment.color="grey45",
                                     segment.size=0.25, box.padding=0.75,
                                     point.padding=0.8, force=6, force_pull=0.08, parse=TRUE) +
            scale_fill_gradient(paste0("Proportion\n", class_to_plot),
                                low = "blue", high = "red", limits = c(0, 1)) +
            scale_color_manual(breaks=c("within_taxid","outside_taxid","no_taxon"),
                               values=c("within_taxid"=histogram_bar_color,
                                        "outside_taxid"=outside_taxid_label_color,
                                        "no_taxon"=no_taxon_label_color),
                               labels=c("Within requested taxids", "Outside requested taxids",
                                        "No hit"),
                               name="Taxon source") +
            guides(fill=guide_colorbar(order=1),
                   size=guide_legend(order=2,
                                     override.aes=list(shape=16, color="black",
                                                       fill="black", alpha=1)),
                   color=guide_legend(order=3,
                                      override.aes=list(label="A", size=4))) +
            theme_minimal() + xlab(expression("Embedding" ~ "\u00D7" ~ beta)) +
            ylab("Levenshtein Distance\n(to most abundant anchor-target)") +
            ggtitle(paste(category, my_cluster, sep=" | "),
                    subtitle=paste("Color: proportion", class_to_plot)) +
            theme(panel.grid.minor.y = element_blank())
          print(p2)
        }

        summ_out_dt <- p_sub %>% select(-any_of("label2")) %>% ungroup() %>%
          mutate(across(all_of(all_classes), \(x) paste(cur_column(), replace_na(x, 0),sep=":"))) %>%
          unite(col=metadata, all_of(all_classes), sep="/") %>%
          mutate(metadata_category = category, cluster=my_cluster, feature=my_feature)
      }

      summary_internal_cols <- c(
        "direct_blast_species.x", "direct_blast_staxids.x",
        "direct_blast_subject_id.x", "direct_blast_accession.x",
        "direct_blast_species.y", "direct_blast_staxids.y",
        "direct_blast_subject_id.y", "direct_blast_accession.y",
        "detail_direct_blast_species", "detail_direct_blast_staxids",
        "detail_direct_blast_subject_id", "detail_direct_blast_accession",
        "compactor_seed_sequence",
        "blastp_compactor_sequence", "blastn_compactor_sequence",
        "compactor_selected_sequence",
        "compactor_summary_label", "compactor_summary_identity",
        "compactor_summary_qcovs", "compactor_summary_species",
        "compactor_summary_staxids", "compactor_summary_sequence",
        "compactor_summary_subject_id", "compactor_summary_accession",
        "compactor_summary_ncbi_protein_accession",
        "compactor_summary_uniprot_accession",
        "point_label", "point_label_expr", "point_label_color",
        "color_value", "mean_metadata", "median_metadata", "sd_metadata",
        "outside_taxid_only"
      )
      summ_out_dt <- summ_out_dt %>%
        select(-any_of(summary_internal_cols)) %>%
        select(-matches("\\.(x|y)$"))

      all_features_summary <- bind_rows(all_features_summary, summ_out_dt)
    }

  }, error = function(e) {
    parent_msg <- ""
    if (!is.null(e$parent) && !is.null(e$parent$message)) {
      parent_msg <- paste0("\nParent error: ", e$parent$message)
    }
    message(paste0("Error processing category: ", category, "\n", e$message, parent_msg))
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
message(paste0("Compactor labels in generated plot summary: ",
               sum(str_detect(all_features_summary$`Blast Label`, "\\(COMPACTOR\\)"), na.rm=TRUE)))
if ("compactor_sequence" %in% colnames(all_features_summary)) {
  message(paste0("Rows with compactor_sequence in generated plot summary: ",
                 sum(!is.na(text_or_na_vec(all_features_summary$compactor_sequence)), na.rm=TRUE)))
}

unannotated_summary <- all_features_summary %>%
  filter(is.na(`Blast Label`) |
           `Blast Label` %in% c("NO MATCH", "UNANNOTATED", "UNCHARACTERISED", "NO PROTEIN/GENE HIT")) %>%
  mutate(total_samples = suppressWarnings(as.numeric(total_samples))) %>%
  mutate(entropy_stats = map2(metadata, total_samples, metadata_entropy_stats)) %>%
  unnest(entropy_stats) %>%
  arrange(desc(metadata_specificity_score), desc(total_samples), metadata_normalized_entropy) %>%
  relocate(metadata_category, feature, cluster)

default_summary_path <- str_replace(opt$output, ".pdf", "_summary.tsv")
if (nchar(opt$compactor_summary) > 0) {
  all_features_summary %>%
    relocate(metadata_category, feature, cluster) %>%
    write_tsv(file = opt$compactor_summary)
  if (normalizePath(opt$compactor_summary, mustWork = FALSE) !=
      normalizePath(default_summary_path, mustWork = FALSE) &&
      file.exists(default_summary_path)) {
    unlink(default_summary_path)
  }
  message(paste0("Wrote enriched compactor summary with resolved species names to: ",
                 opt$compactor_summary))
} else {
  all_features_summary %>% relocate(metadata_category, feature, cluster) %>% write_tsv(file = default_summary_path)
}
unannotated_summary %>% write_tsv(file = str_replace(opt$output, ".pdf", "_unannotated.tsv"))
all_blastp_summary %>% relocate(metadata_category, feature, cluster) %>% write_tsv(file = str_replace(opt$nonzero_annotations, "blastp_annotated.tsv$", "blastp_all.tsv"))
all_blast_summary %>% relocate(metadata_category, feature, cluster) %>% write_tsv(file = str_replace(gsub("blastp_annotated", "blast_annotated", opt$nonzero_annotations), "blast_annotated.tsv$", "blast_all.tsv"))
