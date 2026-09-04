suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))

option_list <- list(
  make_option(c("--project_dir"), type = "character", default = ".",
              help = "FLASH project directory. Accepted for wrapper compatibility."),
  make_option(c("--results_dir"), type = "character", default = NULL,
              help = "Run directory containing the compactor summary."),
  make_option(c("--metadata"), type = "character", default = "",
              help = "Metadata TSV/CSV. Accepted for wrapper compatibility."),
  make_option(c("--metadata_column"), type = "character", default = NULL,
              help = "metadata_category to replot."),
  make_option(c("--clusters"), type = "character", default = NULL,
              help = "Comma-separated numeric cluster IDs, e.g. 8143,5883."),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Output PDF path."),
  make_option(c("--compactor_summary"), type = "character", default = "",
              help = "Explicit blast_annotated_plots_summary_compactor TSV."),
  make_option(c("--num_hits"), type = "integer", default = 100,
              help = "Maximum number of selected cluster/feature plots. Default: 100."),
  make_option(c("--width"), type = "numeric", default = 12,
              help = "PDF width. Default: 12."),
  make_option(c("--height"), type = "numeric", default = 8,
              help = "PDF height. Default: 8."),
  make_option(c("--label_size"), type = "numeric", default = 3.2,
              help = "Repelled point-label text size. Default: 3.2."),
  make_option(c("--hist_label_scale"), type = "numeric", default = 1,
              help = "Multiplier for histogram label text sizes. Default: 1."),
  make_option(c("--products"), action = "store_true", default = FALSE,
              help = "Accepted for wrapper compatibility."),
  make_option(c("--plotter"), type = "character", default = "",
              help = "Accepted for wrapper compatibility."),
  make_option(c("--nonzero_annotations"), type = "character", default = "",
              help = "Accepted for wrapper compatibility."),
  make_option(c("--clusters_file"), type = "character", default = "",
              help = "Accepted for wrapper compatibility."),
  make_option(c("--feather_file"), type = "character", default = "",
              help = "Accepted for wrapper compatibility."),
  make_option(c("--sample_seqs"), type = "character", default = "",
              help = "Accepted for wrapper compatibility."),
  make_option(c("--cluster_length"), type = "integer", default = NA,
              help = "Accepted for wrapper compatibility.")
)

opt <- parse_args(OptionParser(option_list = option_list))

within_taxid_label_color <- "#E64B35"
outside_taxid_label_color <- "#4DBBD5"
histogram_bar_color <- "#F8766D"
no_taxon_label_color <- "#7A7A7A"

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

find_compactor_summary <- function(results_dir) {
  pick_one(
    list.files(
      results_dir,
      pattern = "_nonzero_coefficients_blast_annotated_plots_summary_compactor\\.tsv$",
      full.names = TRUE
    ),
    "compactor plot summary"
  )
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
  stop("Could not find a metadata/category column in the compactor summary.", call. = FALSE)
}

coerce_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_text <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x %in% c("NA", "NaN", "NULL")] <- ""
  x
}

first_nonempty <- function(x) {
  x <- safe_text(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

collapse_unique_labels <- function(x, sep = ";") {
  x <- safe_text(x)
  x <- unlist(strsplit(paste(x, collapse = sep), ";|,", perl = TRUE), use.names = FALSE)
  x <- trimws(x)
  x <- unique(x[nzchar(x)])
  if (length(x) == 0) {
    return("NO MATCH")
  }
  paste(x, collapse = sep)
}

is_placeholder_label <- function(x) {
  x <- safe_text(x)
  x == "" |
    str_detect(x, regex("^(NO MATCH|NO TARGET|UNANNOTATED|UNCHARACTERIZED|UNCHARACTERISED|NO PROTEIN/GENE HIT)$",
                        ignore_case = TRUE))
}

format_metric <- function(metric_value) {
  metric <- coerce_num(metric_value)
  metric <- metric[is.finite(metric)]
  if (length(metric) == 0) {
    return("")
  }
  paste0("Accuracy: ", scales::percent(max(metric, na.rm = TRUE), accuracy = 0.01))
}

parse_classes <- function(value) {
  value <- first_nonempty(value)
  if (is.na(value)) {
    return(character())
  }
  value <- gsub("^\\[|\\]$", "", value)
  classes <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
  classes[nzchar(classes)]
}

parse_metadata_counts <- function(value) {
  value <- first_nonempty(value)
  if (is.na(value)) {
    return(setNames(numeric(0), character(0)))
  }
  parts <- unlist(strsplit(value, "/", fixed = TRUE), use.names = FALSE)
  out <- numeric(0)
  for (part in parts) {
    if (!grepl(":", part, fixed = TRUE)) {
      next
    }
    key <- sub(":.*$", "", part)
    val <- sub("^[^:]*:", "", part)
    num <- suppressWarnings(as.numeric(val))
    if (!is.na(num) && nzchar(key)) {
      out[[key]] <- num
    }
  }
  out
}

metadata_count_table <- function(metadata_strings) {
  parsed <- lapply(metadata_strings, parse_metadata_counts)
  classes <- unique(unlist(lapply(parsed, names), use.names = FALSE))
  classes <- classes[nzchar(classes)]
  if (length(classes) == 0) {
    return(tibble(.row_id = seq_along(metadata_strings)))
  }
  rows <- lapply(seq_along(parsed), function(i) {
    vals <- setNames(rep(0, length(classes)), classes)
    vals[names(parsed[[i]])] <- parsed[[i]]
    as_tibble(as.list(vals)) %>% mutate(.row_id = i)
  })
  bind_rows(rows)
}

class_interest_score <- function(label) {
  normalized <- str_to_lower(as.character(label)) %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()
  negative_pattern <- paste(
    c("uninfected", "untreated", "unexposed", "noninfected", "non infected",
      "not infected", "not treated", "not exposed", "control", "mock", "placebo",
      "healthy", "negative", "none", "absent", "baseline", "wild type",
      "wildtype", "susceptible", "not given", "notgiven", "no", "no fungus"),
    collapse = "|"
  )
  positive_pattern <- paste(
    c("infected", "infection", "treated", "treatment", "exposed", "case", "diseased",
      "positive", "present", "resistant", "mutant", "yes", "fungus"),
    collapse = "|"
  )
  if (str_detect(normalized, paste0("(^| )(", negative_pattern, ")($| )"))) {
    return(-1)
  }
  if (str_detect(normalized, paste0("(^| )(", positive_pattern, ")($| )"))) {
    return(1)
  }
  0
}

choose_interesting_binary_class <- function(classes) {
  classes <- unique(classes[!is.na(classes) & nzchar(classes)])
  if (length(classes) != 2) {
    return(classes)
  }
  scores <- vapply(classes, class_interest_score, numeric(1))
  if (sum(scores == max(scores)) == 1) {
    return(classes[[which.max(scores)]])
  }
  sort(classes)[[2]]
}

prepare_point_labels <- function(dt) {
  for (col in c("Blast Label", "label_identity", "label_coverage", "label_quality")) {
    if (!col %in% colnames(dt)) {
      dt[[col]] <- ""
    }
  }
  dt %>%
    mutate(
      blast_label = safe_text(`Blast Label`),
      label_identity = safe_text(label_identity),
      label_coverage = safe_text(label_coverage),
      label_quality = safe_text(label_quality),
      label_quality = ifelse(
        label_quality == "" & (label_identity != "" | label_coverage != ""),
        paste0("I:", ifelse(label_identity == "", "NA", label_identity),
               "; C:", ifelse(label_coverage == "", "NA", label_coverage)),
        label_quality
      ),
      blast_label = ifelse(blast_label == "", "NO MATCH", blast_label),
      point_label = str_wrap(str_trunc(blast_label, width = 80, side = "right"), width = 28),
      point_label = ifelse(
        !is_placeholder_label(blast_label) &
          nzchar(label_quality) &
          label_quality != "I:100%; C:100%",
        paste(point_label, label_quality, sep = "\n"),
        point_label
      ),
      point_label_color = case_when(
        str_detect(blast_label, "\\(OTHER TAXA\\)") ~ "outside_taxid",
        is_placeholder_label(blast_label) ~ "no_taxon",
        TRUE ~ "within_taxid"
      )
    )
}

filter_summary <- function(path, metadata_column, selected_clusters) {
  dt <- fread(path)
  metadata_col <- detect_metadata_col(dt)
  if (ncol(dt) < 3) {
    stop("Compactor summary must have cluster in column 3: ", path, call. = FALSE)
  }
  cluster_col <- colnames(dt)[[3]]
  if (!"feature" %in% colnames(dt)) {
    stop("Compactor summary is missing required column: feature", call. = FALSE)
  }
  dt[, cluster_normalized := normalize_cluster_id(get(cluster_col))]
  filtered <- dt[get(metadata_col) == metadata_column & cluster_normalized %in% selected_clusters]
  if (nrow(filtered) == 0) {
    available_clusters <- unique(dt[get(metadata_col) == metadata_column, cluster_normalized])
    available_clusters <- head(available_clusters[!is.na(available_clusters)], 20)
    stop(
      "No rows matched metadata_column=", metadata_column,
      " and clusters=", paste(selected_clusters, collapse = ","),
      " in ", path,
      ". Example clusters for this metadata column: ",
      paste(available_clusters, collapse = ","),
      call. = FALSE
    )
  }
  matched_counts <- filtered[, .N, by = cluster_normalized]
  missing_clusters <- setdiff(selected_clusters, matched_counts$cluster_normalized)
  message("Matched rows by requested cluster in ", basename(path), ": ",
          paste(paste0(matched_counts$cluster_normalized, "=", matched_counts$N), collapse = ", "))
  if (length(missing_clusters) > 0) {
    message("Requested cluster(s) with no rows for ", metadata_column, ": ",
            paste(missing_clusters, collapse = ","))
  }
  filtered[, cluster_normalized := NULL]
  if (metadata_col != "metadata_category") {
    setnames(filtered, metadata_col, "metadata_category")
  }
  if (cluster_col != "cluster") {
    setnames(filtered, cluster_col, "cluster")
  }
  as_tibble(filtered)
}

make_histogram_label <- function(label) {
  label <- collapse_unique_labels(label, sep = ",")
  label <- str_replace(label, " ,", ", ") %>% str_replace(" ;", "; ")
  str_wrap(str_trunc(label, width = 120, side = "right"), width = 38)
}

make_histogram_data <- function(dt, num_hits) {
  largest_coef <- max(abs(coerce_num(dt$max_coefficient)), na.rm = TRUE)
  if (!is.finite(largest_coef) || largest_coef <= 0) {
    largest_coef <- 1
  }
  dt %>%
    mutate(max_coefficient = abs(coerce_num(max_coefficient))) %>%
    filter(is.finite(max_coefficient)) %>%
    group_by(cluster, feature) %>%
    summarise(
      max_coefficient = max(max_coefficient, na.rm = TRUE),
      label = make_histogram_label(`Blast Label`),
      accuracy = first_nonempty(accuracy),
      classes = first_nonempty(classes),
      .groups = "drop"
    ) %>%
    arrange(desc(max_coefficient), cluster, feature) %>%
    mutate(
      rank = row_number(),
      coef_mag = max_coefficient / largest_coef,
      taxon_source = case_when(
        str_detect(label, "\\(OTHER TAXA\\)") ~ "outside_taxid",
        is_placeholder_label(label) ~ "no_taxon",
        TRUE ~ "within_taxid"
      ),
      label_size = case_when(
        is_placeholder_label(label) ~ 3.35,
        nchar(label) <= 55 ~ 2.85,
        nchar(label) <= 95 ~ 2.55,
        TRUE ~ 2.25
      ) * opt$hist_label_scale
    ) %>%
    head(num_hits)
}

plot_histograms <- function(dt, hist_dt, title, subtitle) {
  if (nrow(hist_dt) == 0) {
    return(invisible(NULL))
  }
  p <- ggplot(hist_dt, aes(x = rank, y = coef_mag)) +
    geom_col(fill = histogram_bar_color) +
    geom_text(aes(y = coef_mag + 0.05, label = label, hjust = 0, size = label_size),
              angle = 45, color = "#222222", show.legend = FALSE) +
    scale_size_identity() +
    scale_y_continuous(
      "Magnitude relative to\nlargest nonzero coefficient",
      limits = c(0, 1.6),
      labels = scales::label_percent(),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    xlab("Rank of nonzero coefficient (by magnitude)") +
    scale_x_continuous(breaks = seq_len(max(10, max(hist_dt$rank)))) +
    ggtitle(title, subtitle = subtitle) +
    theme_pubr() +
    theme(plot.subtitle = element_text(size = 8, lineheight = 0.95))
  print(p)

  stack_dt <- dt %>%
    distinct(cluster, feature, sequence, `Blast Label`) %>%
    inner_join(hist_dt %>% select(cluster, feature, rank), by = c("cluster", "feature")) %>%
    mutate(taxon_source = case_when(
      str_detect(safe_text(`Blast Label`), "\\(OTHER TAXA\\)") ~ "outside_taxid",
      is_placeholder_label(`Blast Label`) ~ "no_taxon",
      TRUE ~ "within_taxid"
    )) %>%
    count(rank, taxon_source, name = "extendor_n")

  p_stack <- ggplot(stack_dt, aes(x = rank, y = extendor_n, fill = taxon_source)) +
    geom_col(position = "stack") +
    scale_y_continuous("Number of extendors",
                       breaks = scales::breaks_pretty(n = 6),
                       expand = expansion(mult = c(0, 0.08))) +
    xlab("Rank of nonzero coefficient (same order as previous histogram)") +
    scale_x_continuous(breaks = seq_len(max(10, max(hist_dt$rank)))) +
    scale_fill_manual(
      breaks = c("within_taxid", "outside_taxid", "no_taxon"),
      values = c("within_taxid" = histogram_bar_color,
                 "outside_taxid" = outside_taxid_label_color,
                 "no_taxon" = no_taxon_label_color),
      labels = c("Within requested taxids", "Outside requested taxids", "No hit"),
      name = "Taxon source"
    ) +
    ggtitle(title, subtitle = subtitle) +
    theme_pubr() +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 8, lineheight = 0.95))
  print(p_stack)
}

prepare_detail_data <- function(dt) {
  if (!"total_samples" %in% colnames(dt)) {
    dt$total_samples <- NA_real_
  }
  dt <- dt %>%
    mutate(
      .row_id = row_number(),
      embedding = coerce_num(embedding),
      lev_dist = coerce_num(lev_dist),
      total_samples = coerce_num(total_samples)
    )
  counts <- metadata_count_table(dt$metadata)
  dt <- left_join(dt, counts, by = ".row_id")
  class_cols <- setdiff(colnames(counts), ".row_id")
  if (length(class_cols) > 0) {
    dt <- dt %>%
      mutate(total_samples = ifelse(
        is.finite(total_samples),
        total_samples,
        rowSums(across(all_of(class_cols)), na.rm = TRUE)
      ))
  }
  list(data = prepare_point_labels(dt), class_cols = class_cols)
}

plot_detail <- function(detail_dt, class_cols, category, cluster_id, feature_id) {
  p_sub <- detail_dt %>%
    filter(cluster == cluster_id, feature == feature_id) %>%
    filter(is.finite(embedding), !is.na(lev_dist), is.finite(total_samples), total_samples > 0)
  if (nrow(p_sub) == 0) {
    message("Skipping detailed plot for ", category, " ", cluster_id, " ", feature_id,
            " because no finite plotting rows remained.")
    return(invisible(NULL))
  }

  classes_from_model <- parse_classes(first_nonempty(p_sub$classes))
  classes_to_plot <- intersect(classes_from_model, class_cols)
  if (length(classes_to_plot) == 0) {
    classes_to_plot <- class_cols
  }
  if (length(classes_to_plot) == 2) {
    focus <- choose_interesting_binary_class(classes_to_plot)
    classes_to_plot <- focus
    message("Binary target ", category, ": plotting only class '", focus, "'.")
  }

  for (class_to_plot in classes_to_plot) {
    if (!class_to_plot %in% colnames(p_sub)) {
      next
    }
    p_class <- p_sub %>%
      mutate(class_proportion = ifelse(total_samples > 0,
                                       .data[[class_to_plot]] / total_samples,
                                       NA_real_))
    p2 <- ggplot(
      p_class,
      aes(x = embedding, y = lev_dist, fill = class_proportion,
          size = total_samples, label = point_label)
    ) +
      geom_vline(xintercept = 0, lty = "dashed") +
      geom_point(shape = 21, color = "grey25", stroke = 0.35) +
      scale_y_continuous(breaks = scales::breaks_width(1), minor_breaks = NULL) +
      scale_size_continuous(
        trans = "log",
        name = "Total Samples",
        breaks = c(1, 10, 100, 1000, 10000),
        limits = c(1, 10000),
        labels = scales::label_log()
      ) +
      geom_text_repel(
        data = function(x) filter(x, point_label_color != "outside_taxid"),
        aes(color = point_label_color),
        size = opt$label_size,
        max.overlaps = Inf,
        min.segment.length = 0,
        segment.color = "grey45",
        segment.size = 0.25,
        box.padding = 0.75,
        point.padding = 0.8,
        force = 6,
        force_pull = 0.08
      ) +
      geom_text_repel(
        data = function(x) filter(x, point_label_color == "outside_taxid"),
        aes(color = point_label_color),
        size = opt$label_size,
        max.overlaps = Inf,
        min.segment.length = 0,
        segment.color = "grey45",
        segment.size = 0.25,
        box.padding = 0.75,
        point.padding = 0.8,
        force = 6,
        force_pull = 0.08
      ) +
      scale_fill_gradient(paste0("Proportion\n", class_to_plot),
                          low = "blue", high = "red", limits = c(0, 1)) +
      scale_color_manual(
        breaks = c("within_taxid", "outside_taxid", "no_taxon"),
        values = c("within_taxid" = histogram_bar_color,
                   "outside_taxid" = outside_taxid_label_color,
                   "no_taxon" = no_taxon_label_color),
        labels = c("Within requested taxids", "Outside requested taxids", "No hit"),
        name = "Taxon source"
      ) +
      guides(
        fill = guide_colorbar(order = 1),
        size = guide_legend(
          order = 2,
          override.aes = list(shape = 16, color = "black", fill = "black", alpha = 1)
        ),
        color = guide_legend(order = 3, override.aes = list(label = "A", size = 4))
      ) +
      theme_minimal() +
      xlab(expression("Embedding" ~ "\u00D7" ~ beta)) +
      ylab("Levenshtein Distance\n(to most abundant anchor-target)") +
      ggtitle(paste(category, cluster_id, sep = " | "),
              subtitle = paste("Color: proportion", class_to_plot)) +
      theme(panel.grid.minor.y = element_blank())
    print(p2)
  }
}

stop_if_blank(opt$results_dir, "--results_dir")
stop_if_blank(opt$metadata_column, "--metadata_column")
stop_if_blank(opt$clusters, "--clusters")
stop_if_blank(opt$output, "--output")

results_dir <- normalizePath(opt$results_dir, mustWork = TRUE)
summary_path <- if (!is_blank(opt$compactor_summary)) {
  existing_file(opt$compactor_summary, "--compactor_summary")
} else {
  find_compactor_summary(results_dir)
}
selected_clusters <- parse_cluster_list(opt$clusters)

plot_dt <- filter_summary(summary_path, opt$metadata_column, selected_clusters)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

hist_dt <- make_histogram_data(plot_dt, opt$num_hits)
detail <- prepare_detail_data(plot_dt)
detail_dt <- detail$data
class_cols <- detail$class_cols
selected_features <- hist_dt %>% select(cluster, feature)

message("Selected ", nrow(plot_dt), " compactor summary row(s) for ",
        opt$metadata_column, " / ", paste(selected_clusters, collapse = ","))
message("Selected ", nrow(selected_features), " cluster/feature plot(s).")
message("Using compactor summary as standalone plotting input: ", summary_path)

pdf(opt$output, width = opt$width, height = opt$height)
plot(0:10, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
text(5, 8, paste("BLOT manuscript subset"))
text(5, 7, paste("Summary:", basename(summary_path)))
text(5, 6, paste("Metadata:", opt$metadata_column))
text(5, 5, paste("Clusters:", paste(selected_clusters, collapse = ", ")))
text(5, 4, paste(Sys.Date()))

dataset <- str_extract(summary_path, "results/([^/]+)/")
dataset <- ifelse(is.na(dataset), basename(dirname(summary_path)), str_remove_all(dataset, "^results/|/$"))
title <- paste(opt$metadata_column, "in", dataset)
subtitle <- format_metric(plot_dt$accuracy)
plot_histograms(plot_dt, hist_dt, title, subtitle)

for (i in seq_len(nrow(selected_features))) {
  plot_detail(
    detail_dt,
    class_cols,
    opt$metadata_column,
    selected_features$cluster[[i]],
    selected_features$feature[[i]]
  )
}
dev.off()

summary_out <- sub("\\.pdf$", "_summary.tsv", opt$output)
fwrite(as.data.table(plot_dt), summary_out, sep = "\t")
message("Wrote ", opt$output)
message("Wrote ", summary_out)
