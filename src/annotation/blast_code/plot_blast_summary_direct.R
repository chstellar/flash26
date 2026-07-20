suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

option_list <- list(
  make_option(c("--summary"), type="character", default=NULL,
              help="Compactor-filled blast plot summary TSV", metavar="character"),
  make_option(c("--output"), type="character", default=NULL,
              help="Output PDF path", metavar="character"),
  make_option(c("--num_hits"), type="integer", default=10,
              help="Number of top features per category to show in histograms", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$summary) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--summary and --output are required", call.=FALSE)
}

first_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  ifelse(length(x) == 0, NA_real_, x[1])
}

clean_label <- function(x) {
  x <- replace_na(as.character(x), "")
  x <- str_replace_all(x, "[\r\n\t]+", " ")
  x <- str_squish(x)
  x <- ifelse(x %in% c("", "NA", "NO PROTEIN/GENE HIT"), NA_character_, x)
  x
}

preserve_suffix <- function(label, width=85) {
  label <- clean_label(label)
  is_compactor <- str_detect(replace_na(label, ""), "\\s*\\(COMPACTOR\\)\\s*$")
  base <- str_replace(replace_na(label, ""), "\\s*\\(COMPACTOR\\)\\s*$", "")
  ifelse(is_compactor,
         paste0(str_trunc(base, width=max(10, width - 12)), " (COMPACTOR)"),
         str_trunc(label, width=width))
}

metadata_count_vector <- function(metadata_string) {
  entries <- unlist(str_split(replace_na(metadata_string, ""), "/"))
  parsed <- str_match(entries, "^(.+):([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)$")
  parsed <- parsed[!is.na(parsed[, 1]), , drop=FALSE]
  if (nrow(parsed) == 0) {
    return(numeric())
  }
  values <- suppressWarnings(as.numeric(parsed[, 3]))
  names(values) <- parsed[, 2]
  values[!is.na(values)]
}

is_quantitative_metadata <- function(metadata_string) {
  str_detect(replace_na(metadata_string, ""), "^(mean_|median_)")
}

make_hist_label <- function(labels) {
  labels <- clean_label(labels)
  labels <- labels[!is.na(labels)]
  labels <- unique(labels)
  if (length(labels) == 0) {
    return("NO MATCH")
  }
  labels <- labels[order(!str_detect(labels, "\\(COMPACTOR\\)$"),
                         str_detect(labels, "(?i)hypothetical|uncharacterized|predicted"),
                         nchar(labels), labels)]
  paste(head(labels, 5), collapse="; ")
}

summary_dt <- fread(opt$summary, fill=TRUE)
if (nrow(summary_dt) == 0) {
  stop("Summary file is empty: ", opt$summary, call.=FALSE)
}

required <- c("metadata_category", "feature", "cluster", "sequence", "embedding", "lev_dist", "metadata")
missing_required <- setdiff(required, colnames(summary_dt))
if (length(missing_required) > 0) {
  stop("Summary file is missing required columns: ", paste(missing_required, collapse=", "), call.=FALSE)
}

if (!"Blast Label" %in% colnames(summary_dt)) {
  summary_dt$`Blast Label` <- NA_character_
}
if (!"compactor_annotation" %in% colnames(summary_dt)) {
  summary_dt$compactor_annotation <- NA_character_
}
if (!"total_samples" %in% colnames(summary_dt)) {
  summary_dt$total_samples <- NA_real_
}
if (!"identity" %in% colnames(summary_dt)) {
  summary_dt$identity <- NA_real_
}
if (!"qcovs" %in% colnames(summary_dt)) {
  summary_dt$qcovs <- NA_real_
}

plot_dt <- summary_dt %>%
  mutate(embedding = suppressWarnings(as.numeric(embedding)),
         lev_dist = suppressWarnings(as.numeric(lev_dist)),
         total_samples = suppressWarnings(as.numeric(total_samples)),
         total_samples = ifelse(is.na(total_samples), 1, total_samples),
         label = ifelse(!is.na(clean_label(compactor_annotation)),
                        compactor_annotation, `Blast Label`),
         label = clean_label(label),
         label = ifelse(is.na(label), "NO MATCH", label),
         point_label = str_wrap(preserve_suffix(label, width=85), width=28),
         max_abs_embedding = abs(embedding)) %>%
  filter(!is.na(metadata_category), !is.na(feature), !is.na(cluster), !is.na(sequence))

message("Direct summary plot rows loaded: ", nrow(plot_dt))
message("Rows with COMPACTOR labels: ", sum(str_detect(plot_dt$label, "\\(COMPACTOR\\)$"), na.rm=TRUE))

pdf(opt$output, width=11, height=8.5)
plot.new()
text(0.5, 0.65, "Compactor-rescued BLAST summary plots", cex=1.5)
text(0.5, 0.50, opt$summary, cex=0.65)
text(0.5, 0.40, paste("Rows with COMPACTOR labels:", sum(str_detect(plot_dt$label, "\\(COMPACTOR\\)$"), na.rm=TRUE)), cex=0.9)

for (category in unique(plot_dt$metadata_category)) {
  category_dt <- plot_dt %>% filter(metadata_category == category)
  feature_dt <- category_dt %>%
    group_by(feature, cluster) %>%
    summarise(max_abs_embedding = max(abs(embedding), na.rm=TRUE),
              label = make_hist_label(label),
              .groups="drop") %>%
    arrange(desc(max_abs_embedding)) %>%
    mutate(rank=row_number(),
           label = str_wrap(preserve_suffix(label, width=100), width=25),
           color = ifelse(str_detect(label, "\\(COMPACTOR\\)"), "compactor",
                          ifelse(str_detect(label, "^(NO MATCH|NO TARGET|UNANNOTATED)"), "no_blast", "blast")))

  if (nrow(feature_dt) > 0) {
    p_hist <- feature_dt %>% head(opt$num_hits) %>%
      ggplot(aes(x=rank, y=max_abs_embedding / max(max_abs_embedding, na.rm=TRUE), fill=color, label=label)) +
      geom_col() +
      geom_text(aes(y=max_abs_embedding / max(max_abs_embedding, na.rm=TRUE) + 0.04),
                hjust=0, angle=45, size=2.8) +
      scale_fill_manual(values=c(compactor="#5B4B8A", blast="#E7A1B0", no_blast="grey70"),
                        breaks=c("compactor", "blast", "no_blast"),
                        labels=c("Compactor rescue", "Blast hit", "No blast/annotation")) +
      scale_x_continuous(breaks=seq_len(min(opt$num_hits, nrow(feature_dt)))) +
      coord_cartesian(ylim=c(0, 1.45), clip="off") +
      labs(title=category,
           subtitle="Top features from compactor-filled plot summary",
           x="Rank of nonzero coefficient/embedding effect",
           y="Magnitude relative to largest plotted feature",
           fill=NULL) +
      theme_bw() +
      theme(plot.margin=margin(10, 35, 10, 10),
            legend.position="bottom")
    print(p_hist)
  }

  for (feature_name in unique(category_dt$feature)) {
    sub_dt <- category_dt %>% filter(feature == feature_name)
    if (nrow(sub_dt) == 0) {
      next
    }
    quantitative <- any(is_quantitative_metadata(sub_dt$metadata))
    if (quantitative) {
      sub_dt <- sub_dt %>%
        mutate(color_value = suppressWarnings(as.numeric(str_match(metadata, "mean_[^:]+:([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)")[,2])))
      p <- ggplot(sub_dt, aes(x=embedding, y=lev_dist, color=color_value,
                              size=total_samples, label=point_label)) +
        geom_vline(xintercept=0, lty="dashed") +
        geom_point(alpha=0.85) +
        ggrepel::geom_text_repel(size=3, max.overlaps=Inf, min.segment.length=0,
                                 segment.color="grey45", box.padding=0.5,
                                 point.padding=0.45, force=3) +
        scale_color_viridis_c(option="C", na.value="grey70") +
        labs(title=paste(category, feature_name),
             subtitle="Color encodes mean metadata value",
             x="embedding * beta",
             y="Edit distance to most abundant anchor-target",
             color="Mean value",
             size="Samples") +
        theme_bw()
    } else {
      class_dt <- sub_dt %>%
        rowwise() %>%
        mutate(counts=list(metadata_count_vector(metadata))) %>%
        ungroup() %>%
        mutate(dominant_class = map_chr(counts, \(x) ifelse(length(x) == 0, NA_character_, names(x)[which.max(x)])),
               dominant_fraction = map_dbl(counts, \(x) ifelse(length(x) == 0 || sum(x) == 0, NA_real_, max(x) / sum(x))))
      p <- ggplot(class_dt, aes(x=embedding, y=lev_dist, color=dominant_class,
                                alpha=dominant_fraction, size=total_samples,
                                label=point_label)) +
        geom_vline(xintercept=0, lty="dashed") +
        geom_point() +
        ggrepel::geom_text_repel(size=3, max.overlaps=Inf, min.segment.length=0,
                                 segment.color="grey45", box.padding=0.5,
                                 point.padding=0.45, force=3) +
        scale_alpha_continuous(range=c(0.45, 1), na.value=0.65) +
        labs(title=paste(category, feature_name),
             subtitle="Color encodes dominant metadata class for the extendor",
             x="embedding * beta",
             y="Edit distance to most abundant anchor-target",
             color="Dominant class",
             alpha="Dominant fraction",
             size="Samples") +
        theme_bw()
    }
    print(p)
  }
}

dev.off()
message("Wrote direct compactor summary PDF to ", opt$output)
