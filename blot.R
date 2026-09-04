suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--project_dir"), type = "character", default = ".",
              help = "FLASH project directory. Default: current directory."),
  make_option(c("--results_dir"), type = "character", default = NULL,
              help = "Run directory containing the annotated nonzero coefficient tables."),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Metadata TSV/CSV passed to the blast plotter."),
  make_option(c("--metadata_column"), type = "character", default = NULL,
              help = "Selected metadata_category/metadata column to replot."),
  make_option(c("--clusters"), type = "character", default = NULL,
              help = "Comma-separated numeric cluster IDs to replot, e.g. 8143,5883. Numeric IDs are matched exactly to cluster_NNN entries."),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output PDF path."),
  make_option(c("--plotter"), type = "character",
              default = "src/annotation/blast_code/plot_blast_annotations_each_feature.R",
              help = "Underlying blast plotter R script."),
  make_option(c("--nonzero_annotations"), type = "character", default = "",
              help = "Optional explicit blastp_annotated TSV."),
  make_option(c("--compactor_summary"), type = "character", default = "",
              help = "Optional explicit blast_annotated_plots_summary_compactor TSV."),
  make_option(c("--clusters_file"), type = "character", default = "",
              help = "Optional explicit sequences_per_cluster TSV."),
  make_option(c("--feather_file"), type = "character", default = "",
              help = "Optional explicit feature matrix feather file."),
  make_option(c("--sample_seqs"), type = "character", default = "",
              help = "Optional explicit prepared sample_sequences TSV."),
  make_option(c("--num_hits"), type = "integer", default = 100,
              help = "Number of selected rows to let the plotter render after filtering. Default: 100."),
  make_option(c("--cluster_length"), type = "integer", default = NA,
              help = "Optional cluster length override passed through to the plotter."),
  make_option(c("--products"), action = "store_true", default = FALSE,
              help = "Prefer product names over gene names in labels.")
)

opt <- parse_args(OptionParser(option_list = option_list))

is_blank <- function(x) {
  is.null(x) || is.na(x) || !nzchar(x)
}

stop_if_blank <- function(value, name) {
  if (is_blank(value)) {
    stop(name, " must be supplied.", call. = FALSE)
  }
}

existing_file <- function(path, label) {
  if (is_blank(path) || !file.exists(path) || file.info(path)$size == 0) {
    stop(label, " does not exist or is empty: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

pick_one <- function(paths, label) {
  paths <- unique(paths[file.exists(paths) & file.info(paths)$size > 0])
  if (length(paths) == 0) {
    stop("Could not find ", label, ".", call. = FALSE)
  }
  if (length(paths) > 1) {
    message("Multiple ", label, " candidates found; using: ", paths[[1]])
  }
  normalizePath(paths[[1]], mustWork = TRUE)
}

list_matching_files <- function(root, pattern) {
  if (is_blank(root) || !dir.exists(root)) {
    return(character())
  }
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
}

normalize_cluster_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^cluster_", "", x, ignore.case = TRUE)
  x <- sub("\\.0+$", "", x)
  paste0("cluster_", x)
}

parse_cluster_list <- function(x) {
  clusters <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  clusters <- clusters[nzchar(clusters)]
  if (length(clusters) == 0) {
    stop("--clusters did not contain any comma-delimited cluster IDs.", call. = FALSE)
  }
  normalize_cluster_id(clusters)
}

detect_metadata_col <- function(dt) {
  candidates <- c("metadata_category", "metadata_column", "metadata")
  found <- candidates[candidates %in% colnames(dt)]
  if (length(found) > 0) {
    return(found[[1]])
  }
  stop("Could not find a metadata column. Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
}

detect_cluster_col <- function(dt, preferred_index = NA_integer_) {
  cluster_like_counts <- vapply(dt, function(col) {
    sum(grepl("^cluster_[0-9]+$", trimws(as.character(col)), ignore.case = TRUE), na.rm = TRUE)
  }, integer(1))
  if (!is.na(preferred_index) && ncol(dt) >= preferred_index &&
      cluster_like_counts[[preferred_index]] > 0) {
    return(colnames(dt)[[preferred_index]])
  }
  if (ncol(dt) >= 20 && cluster_like_counts[[20]] > 0) {
    return(colnames(dt)[[20]])
  }
  if (max(cluster_like_counts) > 0) {
    return(names(which.max(cluster_like_counts)))
  }
  candidates <- c("cluster", "cluster_id", "feature_cluster")
  found <- candidates[candidates %in% colnames(dt)]
  if (length(found) > 0) {
    return(found[[1]])
  }
  if (ncol(dt) >= 20) {
    return(colnames(dt)[[20]])
  }
  stop("Could not find a cluster column by name, column 20, or cluster_NNN pattern.", call. = FALSE)
}

filter_compactor_summary <- function(input_path, output_path, metadata_column, selected_clusters) {
  dt <- fread(input_path)
  metadata_col <- detect_metadata_col(dt)
  if (ncol(dt) < 3) {
    stop("Compactor summary must have cluster in column 3, but it has fewer than 3 columns: ",
         input_path, call. = FALSE)
  }
  cluster_col <- colnames(dt)[[3]]
  dt[, cluster_normalized := normalize_cluster_id(get(cluster_col))]
  filtered <- dt[get(metadata_col) == metadata_column & cluster_normalized %in% selected_clusters]
  filtered[, cluster_normalized := NULL]
  if (nrow(filtered) == 0) {
    available_clusters <- unique(dt[get(metadata_col) == metadata_column, normalize_cluster_id(get(cluster_col))])
    available_clusters <- head(available_clusters[!is.na(available_clusters)], 12)
    stop("No rows matched metadata_column=", metadata_column,
         " and clusters=", paste(selected_clusters, collapse = ","),
         " in ", input_path,
         ". Detected metadata column '", metadata_col,
         "' and cluster column '", cluster_col, "'.",
         ifelse(length(available_clusters) > 0,
                paste0(" Example clusters for this metadata column: ",
                       paste(available_clusters, collapse = ","), "."),
                " No rows were found for this metadata column."),
         call. = FALSE)
  }
  if (cluster_col != "cluster") {
    setnames(filtered, cluster_col, "cluster")
  }
  matched_counts <- dt[get(metadata_col) == metadata_column &
                       cluster_normalized %in% selected_clusters,
                       .N,
                       by = cluster_normalized]
  missing_clusters <- setdiff(selected_clusters, matched_counts$cluster_normalized)
  message("Matched rows by requested cluster in ", basename(input_path), ": ",
          paste(paste0(matched_counts$cluster_normalized, "=", matched_counts$N), collapse = ", "))
  if (length(missing_clusters) > 0) {
    message("Requested cluster(s) with no rows for ", metadata_column, ": ",
            paste(missing_clusters, collapse = ","))
  }
  fwrite(filtered, output_path, sep = "\t")
  filtered
}

find_compactor_summary <- function(results_dir) {
  paths <- list.files(
    results_dir,
    pattern = "_nonzero_coefficients_blast_annotated_plots_summary_compactor\\.tsv$",
    full.names = TRUE
  )
  pick_one(paths, "compactor plot summary")
}

write_blot_plotter <- function(plotter, output_path) {
  code <- readLines(plotter, warn = FALSE)
  hook <- c(
    "",
    "blot_selection_path <- Sys.getenv('BLOT_SELECTION_TSV')",
    "blot_normalize_cluster_id <- function(x) {",
    "  x <- trimws(as.character(x))",
    "  x <- sub('^cluster_', '', x, ignore.case = TRUE)",
    "  x <- sub('\\\\.0+$', '', x)",
    "  paste0('cluster_', x)",
    "}",
    "blot_selection <- data.table()",
    "if (nzchar(blot_selection_path) && file.exists(blot_selection_path)) {",
    "  blot_selection <- fread(blot_selection_path) %>%",
    "    mutate(metadata_category = as.character(metadata_category),",
    "           cluster = blot_normalize_cluster_id(cluster),",
    "           feature = as.character(feature)) %>%",
    "    distinct(metadata_category, cluster, feature)",
    "  message('BLOT selection loaded: ', nrow(blot_selection), ' metadata/cluster/feature row(s).')",
    "}",
    "blot_filter_annotation_table <- function(tbl, label) {",
    "  if (nrow(blot_selection) == 0) {",
    "    return(tbl)",
    "  }",
    "  if (!all(c('metadata_category', 'feature') %in% colnames(tbl))) {",
    "    stop('BLOT cannot filter ', label, ' table because it lacks metadata_category or feature.')",
    "  }",
    "  keep <- tbl %>%",
    "    mutate(metadata_category = as.character(metadata_category),",
    "           feature = as.character(feature)) %>%",
    "    inner_join(blot_selection, by = c('metadata_category', 'feature'), relationship = 'many-to-many')",
    "  if ('cluster.x' %in% colnames(keep)) {",
    "    keep <- keep %>% select(-cluster.x) %>% rename(cluster = cluster.y)",
    "  }",
    "  if ('cluster.y' %in% colnames(keep)) {",
    "    keep <- keep %>% rename(cluster = cluster.y)",
    "  }",
    "  message('BLOT kept ', nrow(keep), ' row(s) from ', label, ' using compactor-summary-selected features.')",
    "  keep",
    "}",
    ""
  )
  parse_idx <- which(code == "opt <- parse_args(opt_parser)")
  if (length(parse_idx) != 1) {
    stop("Could not find parse_args hook point in plotter: ", plotter, call. = FALSE)
  }
  code <- append(code, hook, after = parse_idx)
  code <- sub(
    "^dt <- fread\\(opt\\$nonzero_annotations\\)$",
    "dt <- fread(opt$nonzero_annotations)\ndt <- blot_filter_annotation_table(dt, 'blastp annotated nonzero')",
    code
  )
  code <- sub(
    "^if \\(TRUE\\) \\{dt2 <- fread\\(gsub\\(\"blastp_annotated\", \"blast_annotated\", opt\\$nonzero_annotations\\)\\)\\}$",
    "if (TRUE) {dt2 <- fread(gsub(\"blastp_annotated\", \"blast_annotated\", opt$nonzero_annotations))}\ndt2 <- blot_filter_annotation_table(dt2, 'blast annotated nonzero')",
    code
  )
  if (!any(grepl("dt <- blot_filter_annotation_table\\(dt,", code)) ||
      !any(grepl("dt2 <- blot_filter_annotation_table\\(dt2,", code))) {
    stop("Could not inject BLOT filters into plotter: ", plotter, call. = FALSE)
  }
  writeLines(code, output_path)
  normalizePath(output_path, mustWork = TRUE)
}

derive_run_tokens <- function(nonzero_path) {
  base <- basename(nonzero_path)
  list(
    top = sub(".*_results_top([0-9]+).*", "\\1", base),
    target = sub(".*_target([0-9]+).*", "\\1", base),
    k = sub(".*_k([0-9]+)_s.*", "\\1", base),
    s = sub(".*_s([0-9]+)_trainProp.*", "\\1", base)
  )
}

token_is_valid <- function(x) {
  !is_blank(x) && grepl("^[0-9]+$", x)
}

find_by_context <- function(project_dir, results_dir, nonzero_path, kind) {
  results_dir <- normalizePath(results_dir, mustWork = TRUE)
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  rel_parts <- strsplit(gsub("\\\\", "/", results_dir), "/", fixed = TRUE)[[1]]
  dataset <- rel_parts[length(rel_parts) - 4]
  select_type <- rel_parts[length(rel_parts) - 3]
  cluster_type <- rel_parts[length(rel_parts) - 2]
  model <- rel_parts[length(rel_parts) - 1]
  normalize <- rel_parts[length(rel_parts)]
  tokens <- derive_run_tokens(nonzero_path)
  dataset_root <- file.path(project_dir, "results", dataset)

  if (kind == "clusters" &&
      all(vapply(tokens[c("top", "target", "k", "s")], token_is_valid, logical(1)))) {
    candidate <- file.path(project_dir, "results", dataset, select_type, cluster_type,
                           paste0(dataset, "_sequences_per_cluster_top", tokens$top,
                                  "-clusters_target", tokens$target, "_k", tokens$k, "_s", tokens$s, ".tsv"))
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  if (kind == "feather" &&
      all(vapply(tokens[c("top", "target", "k", "s")], token_is_valid, logical(1)))) {
    candidate <- file.path(project_dir, "results", dataset,
                           paste0(dataset, "_", model,
                                  "_top_variance_features_for_glmnet_", select_type, "_", cluster_type,
                                  "_top", tokens$top, "_target", tokens$target,
                                  "_k", tokens$k, "_s", tokens$s, "_", normalize, ".feather"))
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  if (kind == "sample_seqs" &&
      all(vapply(tokens[c("top", "target", "k", "s")], token_is_valid, logical(1)))) {
    candidate <- file.path(project_dir, "results", dataset,
                           paste0(dataset, "_prepared_sequences_", select_type, "_", cluster_type,
                                  "_top", tokens$top, "_target", tokens$target,
                                  "_k", tokens$k, "_s", tokens$s, "_sample_sequences.tsv"))
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  pattern <- switch(
    kind,
    clusters = "sequences_per_cluster.*clusters.*\\.tsv$",
    feather = "features_for_glmnet.*\\.feather$",
    sample_seqs = "prepared_sequences.*sample_sequences\\.tsv$",
    stop("Unknown kind: ", kind, call. = FALSE)
  )
  pick_one(list_matching_files(dataset_root, pattern), kind)
}

stop_if_blank(opt$results_dir, "--results_dir")
stop_if_blank(opt$metadata, "--metadata")
stop_if_blank(opt$metadata_column, "--metadata_column")
stop_if_blank(opt$clusters, "--clusters")
stop_if_blank(opt$output, "--output")

project_dir <- normalizePath(opt$project_dir, mustWork = TRUE)
results_dir <- normalizePath(opt$results_dir, mustWork = TRUE)
metadata <- existing_file(opt$metadata, "--metadata")
plotter_path <- if (grepl("^(/|[A-Za-z]:[/\\\\])", opt$plotter)) opt$plotter else file.path(project_dir, opt$plotter)
plotter <- existing_file(plotter_path, "--plotter")
selected_clusters <- parse_cluster_list(opt$clusters)

nonzero <- if (!is_blank(opt$nonzero_annotations)) {
  existing_file(opt$nonzero_annotations, "--nonzero_annotations")
} else {
  pick_one(list.files(results_dir, pattern = "_nonzero_coefficients_blastp_annotated\\.tsv$",
                      full.names = TRUE), "blastp annotated nonzero table")
}
blastn <- sub("blastp_annotated", "blast_annotated", nonzero, fixed = TRUE)
blastn <- existing_file(blastn, "matching blast annotated nonzero table")
compactor_summary <- if (!is_blank(opt$compactor_summary)) {
  existing_file(opt$compactor_summary, "--compactor_summary")
} else {
  find_compactor_summary(results_dir)
}

clusters_file <- if (!is_blank(opt$clusters_file)) {
  existing_file(opt$clusters_file, "--clusters_file")
} else {
  find_by_context(project_dir, results_dir, nonzero, "clusters")
}
feather_file <- if (!is_blank(opt$feather_file)) {
  existing_file(opt$feather_file, "--feather_file")
} else {
  find_by_context(project_dir, results_dir, nonzero, "feather")
}
sample_seqs <- if (!is_blank(opt$sample_seqs)) {
  existing_file(opt$sample_seqs, "--sample_seqs")
} else {
  find_by_context(project_dir, results_dir, nonzero, "sample_seqs")
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
work_dir <- tempfile("blot_", tmpdir = dirname(opt$output))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

filtered_blastp <- file.path(work_dir, basename(nonzero))
filtered_blastn <- file.path(work_dir, basename(blastn))
filtered_compactor_summary <- file.path(work_dir, "selected_nonzero_coefficients_blast_annotated_plots_summary_compactor.tsv")
blot_plotter <- file.path(work_dir, "plot_blast_annotations_each_feature_blot.R")
filtered_summary <- filter_compactor_summary(
  compactor_summary,
  filtered_compactor_summary,
  opt$metadata_column,
  selected_clusters
)
file.copy(nonzero, filtered_blastp, overwrite = TRUE)
file.copy(blastn, filtered_blastn, overwrite = TRUE)
blot_plotter <- write_blot_plotter(plotter, blot_plotter)

message("Selected ", nrow(filtered_summary), " compactor summary row(s) for ",
        opt$metadata_column, " / ", paste(selected_clusters, collapse = ","))
message("Using clusters: ", clusters_file)
message("Using feather: ", feather_file)
message("Using sample sequences: ", sample_seqs)
message("Using filtered compactor summary: ", filtered_compactor_summary)
message("Using temporary regular plotter with BLOT selection hook: ", blot_plotter)

cmd_args <- c(
  "--vanilla", blot_plotter,
  "--nonzero_annotations", filtered_blastp,
  "--clusters", clusters_file,
  "--feather_file", feather_file,
  "--sample_seqs", sample_seqs,
  "--metadata", metadata,
  "--output", opt$output,
  "--compactor_summary", filtered_compactor_summary,
  "--num_hits", as.character(opt$num_hits)
)
if (isTRUE(opt$products)) {
  cmd_args <- c(cmd_args, "--products")
}
if (!is.na(opt$cluster_length)) {
  cmd_args <- c(cmd_args, "--cluster_length", as.character(opt$cluster_length))
}

status <- system2("Rscript", cmd_args, env = paste0("BLOT_SELECTION_TSV=", filtered_compactor_summary))
if (!identical(status, 0L)) {
  stop("Underlying blast plotter failed with exit status ", status, call. = FALSE)
}
message("Wrote ", opt$output)
