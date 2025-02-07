# summary_plots.R
# Daniel Cotter 
# takes in a path to a dataset folder and produces a summary plot of the accuracies of 
# all of the models run on the dataset 

## Load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

## Parse arguments
option_list = list(
  make_option(c("-o", "--output_folder"), type="character", default=NULL, 
              help="Path to the output folder")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$output_folder)) {
  stop("Please provide an output location")
}

## Load data
# strip tailing slash from dataset path if it exists
library(ggpubr)

## plot all data in one plot
load_cmd1 <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients.tsv")
data1 <- fread(cmd=load_cmd1, header=FALSE, sep="\t", 
               col.names = c("path_metadata", "feature", "accuracy",
                             "sensitivity", "specificity", "classes", "coefficients"))

load_cmd2 <- paste0("grep cluster ", "results/*/filter*/*/*/*_nonzero_coefficients.tsv")
data2 <- fread(cmd=load_cmd2, header=FALSE, sep="\t", 
               col.names = c("path_metadata", "feature", "accuracy",
                             "sensitivity", "specificity", "classes", "coefficients"))

data <- rbind(data1,data2)
data <- data %>% 
  distinct(path_metadata, .keep_all = T) %>%
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity)

# extract the dataset name from the path
data <- data %>% mutate(num_clusters = str_extract(path, "top(\\d+)", group=1))
data <- data %>% mutate(path=dirname(gsub("^results/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  dplyr::rename(paramater_set=path)
data <- data %>% mutate(paramater_set=str_replace(paramater_set, "_", "/")) %>% 
  separate(paramater_set, into=c("dataset", "paramater_set"), sep="/")
data <- data %>% mutate(model=str_extract(paramater_set, 
                                          'hyenaHG38_normalized|hyenaHG38_unnormalized|hyenaMarlowe_normalized|hyenaMarlowe_unnormalized|esm_normalized|esm_unnormalized|hyena_normalized|hyena_unnormalized|ohe')) %>% arrange(desc(model))

data <- data %>% mutate(filter = str_extract(paramater_set, "(filter\\d)_",group=1)) %>% 
  mutate(cluster_approach = str_extract(paramater_set, "filter\\d_([A-Za-z-2]+)_", group=1))

data %>% filter(filter =="filter1" ) %>% filter(cluster_approach=="shiftDist-keepTopES") %>% filter(num_clusters=="10000") %>% ggplot(aes(x=dataset, y=specificity)) + geom_boxplot(aes(fill=model)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_brewer(type="qual") + scale_fill_brewer(type="qual")

data %>% distinct(dataset) %>% arrange(dataset) %>% filter(dataset!="bartlau-phage-infection") %>% pull(dataset) -> unique_datasets

for (i in 1:length(unique_datasets)) {
  DATASET=unique_datasets[i]
  accuracy_data <- data %>% 
    filter(dataset==DATASET) %>% 
    filter(num_clusters=="10000") %>%
    select(filter, cluster_approach, metadata, model, accuracy, num_clusters) %>% 
    distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
    mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) %>% 
    pivot_wider(names_from = model, values_from = accuracy)
  specificity_data <- data %>% 
    filter(dataset==DATASET) %>% 
    filter(num_clusters=="10000") %>%
    select(filter, cluster_approach, metadata, model, specificity, num_clusters) %>% 
    distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
    mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) %>% 
    pivot_wider(names_from = model, values_from = specificity)
  sensitivity_data <- data %>% 
    filter(dataset==DATASET) %>% 
    filter(num_clusters=="10000") %>%
    select(filter, cluster_approach, metadata, model, sensitivity, num_clusters) %>% 
    distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
    mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) %>% 
    pivot_wider(names_from = model, values_from = sensitivity)
  all_cols <- colnames(accuracy_data %>% select(-c(filter, cluster_approach, metadata, `ohe (10000 clusters)`)))
  pdf(file.path(opt$output_folder, paste0(DATASET, "_model_comparisons_", format(Sys.time(), "%Y%m%d"), ".pdf")))
  for (i in all_cols) {
    y_col = sym(i)
    p <- accuracy_data %>% 
      ggplot(aes(x=`ohe (10000 clusters)`,y=!!y_col)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
      geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
      ggrepel::geom_text_repel(aes(label=metadata)) +
      ggtitle("Accuracy")
    print(p)
    p2 <- specificity_data %>% 
      ggplot(aes(x=`ohe (10000 clusters)`,y=!!y_col)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
      geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
      ggrepel::geom_text_repel(aes(label=metadata)) +
      ggtitle("Specificity")
    print(p2)
    p3 <- sensitivity_data %>% 
      ggplot(aes(x=`ohe (10000 clusters)`,y=!!y_col)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
      geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
      ggrepel::geom_text_repel(aes(label=metadata)) +
      ggtitle("Sensitivity")
    print(p3)
  }
  dev.off()
}
