suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggtext))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(fossil))
suppressPackageStartupMessages(library(dbscan))

option_list <- list(
  make_option(
    c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "Path to metadata file",
    metavar = "character"
  ),
  make_option(
    c("-f", "--embeddings"),
    type = "character",
    default = NULL,
    help = "Path to feather file",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Output file for UMAP plot",
    metavar = "character"
  ),
  make_option(
    c("-n", "--num_PCs"),
    type = "integer",
    default = 10,
    help = "Number of principal components for UMAP",
    metavar = "integer"
  ),
  make_option(
    c("-c", "--metadata_column"),
    type = "character",
    default = NA,
    help = "Metadata column to use for coloring",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$metadata) || is.null(opt$embeddings)) {
  stop("Both metadata and feather file paths must be provided.")
}

if (is.null(opt$output)) {
  stop("Output file path must be provided.")
}

metadata_path <- opt$metadata
feather_path <- opt$embeddings
num_pcs <- opt$num_PCs

# Load metadata
metadata <- fread(metadata_path,
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(sample_name = as.character(sample_name)) %>%
  select(sample_name, everything()) # Ensure sample_name is the first column

# read the feather file
full_dt <- feather::read_feather(feather_path)
dataset_name <- tools::file_path_sans_ext(basename(feather_path))
# Extract the dataset name before the first underscore
dataset_name <- str_extract(dataset_name, "^[^_]+")
dt <- full_dt

pca <- dt %>%
  column_to_rownames("sample_name") %>%
  as.matrix() %>%
  t() %>%
  prcomp()

pca_data <- pca$rotation %>%
  as.data.frame() %>%
  rownames_to_column("sample_name")

pca_sub <- pca_data
n <- 10 # num PCs for umap

if (n > ncol(pca_sub) - 1) {
  stop("Number of PCs requested exceeds available PCs in the data.")
}
custom.settings <- umap.defaults
custom.settings$n_epochs <- 500

if (nrow(pca$rotation) > 15) {
  umap_res <- pca_sub %>%
    ungroup() %>%
    select(PC1:!!as.name(paste0("PC", n))) %>%
    umap(config = custom.settings)

  layout <- umap_res$layout %>%
    as.data.frame() %>%
    rename(UMAP1 = V1, UMAP2 = V2)
  # use hdbscan to cluster the UMAP results
  hdbscan_result <- hdbscan(layout[, c("UMAP1", "UMAP2")], minPts = 5) # Adjust minPts as needed
  layout$Cluster <- hdbscan_result$cluster

  umap_df <- pca_sub %>%
    select(-starts_with("PC")) %>%
    cbind(layout)

  if (is.na(opt$metadata_column)) {
    my_metadata <- "HDBSCAN_Cluster"
    umap_df <- umap_df %>%
      mutate(Cluster = as.factor(Cluster)) %>%
      select(sample_name, UMAP1, UMAP2, Cluster) %>%
      rename(!!as.name(my_metadata) := Cluster)
    rand_idx <- NA
  } else {
    my_metadata <- opt$metadata_column
    umap_df <- umap_df %>%
      left_join(metadata, by = "sample_name") %>%
      select(sample_name, UMAP1, UMAP2, Cluster, !!as.name(my_metadata))

    # calculate the random index between the metadata and the clusters
    rand_idx <- round(fossil::rand.index(as.numeric(factor(umap_df[[my_metadata]])), umap_df$Cluster), 3)
  }

  # label the top 8 most abundant classes and change the rest to "other"
  classes <- umap_df %>%
    group_by(!!as.name(my_metadata)) %>%
    summarise(Count = n()) %>%
    top_n(8, Count) %>%
    pull(!!as.name(my_metadata))

  umap_df <- umap_df %>%
    mutate(class = ifelse(
      !!as.name(my_metadata) %in% classes, !!as.name(my_metadata),
      NA
    )) %>%
    mutate(class = factor(class, levels = c(classes, NA)))

  p <- umap_df %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = class)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    geom_point(alpha = 0.8, size = 3) +
    scale_color_brewer(
      paste(my_metadata),
      type = "qual",
      palette = "Set2",
      na.value = "grey",
      labels = function(breaks) {
        breaks[is.na(breaks)] <- "other"
        breaks
      }
    ) +
    ggtitle(
      paste("UMAP of", my_metadata, "in", dataset_name),
      subtitle = paste("rand index:", rand_idx, "\nN_PCs:", n)
    ) +
    theme_minimal()
} else {
  stop("Not enough PCs to run UMAP. Please increase the number of PCs or check your data.")
}

# Save the UMAP plot
ggsave(
  plot = p,
  filename = opt$output,
  width = 6,
  height = 5
)
