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

system(paste("mkdir -p", opt$output_folder))

## Load data
# strip tailing slash from dataset path if it exists
library(ggpubr)


## plot all data in one plot
load_cmd1 <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients_blast_annotated.tsv")
data1 <- fread(cmd=load_cmd1, header=FALSE, sep="\t", 
               col.names = c("path_metadata", "feature", "accuracy",
                             "sensitivity", "specificity", "classes", "coefficients", 
                             "cluster", "query", "identity", "features", "features_10000_window"))

load_cmd2 <- paste0("grep cluster ", "results/*/filter*/*/*/*_nonzero_coefficients_blast_annotated.tsv")
data2 <- fread(cmd=load_cmd2, header=FALSE, sep="\t", 
               col.names = c("path_metadata", "feature", "accuracy",
                             "sensitivity", "specificity", "classes", "coefficients", 
                             "cluster", "query", "identity", "features", "features_10000_window"))


data <- rbind(data1,data2) %>% select(-features_10000_window)

data <- data %>% separate_longer_delim(features, delim = "},") %>% 
  mutate(products=str_extract(features, "'product': \\['([\\w\\s-]+)'\\]", group=1)) %>% 
  mutate(genes=str_extract(features, "'gene': \\['([\\w\\s-]+)'\\]", group=1)) %>% 
  select(-features) %>%
  group_by(path_metadata, accuracy, sensitivity, specificity) %>% 
  summarise(products=list(unique(products)),
            genes=list(unique(genes))) %>% ungroup()

data <- data %>% 
  distinct(path_metadata, .keep_all = T) %>%
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity, products, genes) %>%
  filter(metadata != "metadata_category")

# extract the dataset name from the path
data <- data %>% mutate(num_clusters = str_extract(path, "top(\\d+)", group=1))
data <- data %>% mutate(path=dirname(gsub("^results/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  dplyr::rename(paramater_set=path)
data <- data %>% mutate(paramater_set=str_replace(paramater_set, "_", "/")) %>% 
  separate(paramater_set, into=c("dataset", "paramater_set"), sep="/")
data <- data %>% mutate(model=str_extract(paramater_set, 
                                          'hyenaHG38_normalized|hyenaHG38_unnormalized|hyenaMarlowe_normalized|hyenaMarlowe_unnormalized|esm_normalized|esm_unnormalized|hyena_normalized|hyena_unnormalized|ohe')) %>% arrange(desc(model))

data <- data %>% mutate(accuracy=as.numeric(accuracy), sensitivity=as.numeric(sensitivity), specificity=as.numeric(specificity))

data <- data %>% mutate(filter = str_extract(paramater_set, "(filter\\d)_",group=1)) %>% 
  mutate(cluster_approach = str_extract(paramater_set, "filter\\d_([A-Za-z-2]+)_", group=1))

accuracy_dfs <- list()

system_time <- format(Sys.time(), "%Y%m%d")

pdf(file.path(opt$output_folder, paste0("all_datsets", "_model_comparisons_", system_time, ".pdf")))

merged_accuracy_data <- data %>% 
  filter(num_clusters=="20000") %>%
  filter(!(dataset %in% c("prism-metabolomics-SRP129027", "vibrio-cholerae-PRJNA723557"))) %>%
  filter(cluster_approach %in% c("shiftDist-levFilter", "masked-aa-clustered")) %>%
  filter(model %in% c("ohe", "hyena_normalized", "hyena_unnormalized", "esm_normalized", "esm_unnormalized")) %>%
  select(dataset, filter, cluster_approach, metadata, model, accuracy, num_clusters, products, genes) %>% 
  distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
  mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) 

merged_accuracy_data %>% ungroup() %>% select(-products, -genes) %>% 
  ggplot(aes(x=dataset,y=accuracy,fill=model)) + 
  geom_boxplot() + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_brewer(type="qual") + scale_fill_brewer(type="qual") +
  theme(legend.position="bottom")

merged_accuracy_data %>% ungroup() %>% select(-products, -genes) %>% 
  ggplot(aes(x=dataset,y=accuracy,fill=cluster_approach)) + 
  geom_boxplot() + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_brewer(type="qual") + scale_fill_brewer(type="qual") +
  theme(legend.position="bottom")

merged_accuracy_data %>% ungroup() %>% select(-products, -genes) %>% 
  pivot_wider(id_cols=dataset:model, names_from = cluster_approach, values_from=accuracy) %>%
  ggplot(aes(x=`masked-aa-clustered`,y=`shiftDist-levFilter`,color=dataset,shape=model)) + 
  geom_abline(intercept=0,slope=1,lty="dashed") +
  geom_point() + theme_minimal() + coord_cartesian(xlim=c(0,1),ylim=c(0,1))

unique_datasets <- merged_accuracy_data %>% pull(dataset) %>% unique()

for (i in 1:length(unique_datasets)) {
  DATASET=unique_datasets[i]
  accuracy_data <- data %>% 
    filter(dataset==DATASET) %>% 
    filter(num_clusters=="20000") %>%
    select(filter, cluster_approach, metadata, model, accuracy, num_clusters) %>% 
    distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
    mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) %>% 
    pivot_wider(names_from = model, values_from = accuracy)
  specificity_data <- data %>% 
    filter(dataset==DATASET) %>% 
    filter(num_clusters=="20000") %>%
    select(filter, cluster_approach, metadata, model, specificity, num_clusters) %>% 
    distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
    mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) %>% 
    pivot_wider(names_from = model, values_from = specificity)
  sensitivity_data <- data %>% 
    filter(dataset==DATASET) %>% 
    filter(num_clusters=="20000") %>%
    select(filter, cluster_approach, metadata, model, sensitivity, num_clusters) %>% 
    distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
    mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) %>% 
    pivot_wider(names_from = model, values_from = sensitivity)
  all_cols <- colnames(accuracy_data %>% select(-c(filter, cluster_approach, metadata, `ohe (20000 clusters)`)))
  #pdf(file.path(opt$output_folder, paste0(DATASET, "_model_comparisons_", format(Sys.time(), "%Y%m%d"), ".pdf")))
  for (j in all_cols) {
    y_col = sym(j)
    p <- accuracy_data %>% 
      ggplot(aes(x=`ohe (20000 clusters)`,y=!!y_col)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
      geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
      ggrepel::geom_text_repel(aes(label=metadata)) +
      ggtitle(paste(DATASET, "- Accuracy"))
    print(p)
    if (sum(!is.na(specificity_data %>% pull(y_col)))>0) {
      p2 <- specificity_data %>% 
        ggplot(aes(x=`ohe (20000 clusters)`,y=!!y_col)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
        geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
        ggrepel::geom_text_repel(aes(label=metadata)) +
        ggtitle(paste(DATASET, "-Specificity"))
      print(p2)
    }
    if (sum(!is.na(sensitivity_data %>% pull(y_col)))>0) {
      p3 <- sensitivity_data %>% 
        ggplot(aes(x=`ohe (20000 clusters)`,y=!!y_col)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
        geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
        ggrepel::geom_text_repel(aes(label=metadata)) +
        ggtitle(paste(DATASET, "-Sensitivity"))
      print(p3)
    }
  }
  #dev.off()
  
  summ_accuracy <- data %>% 
    filter(dataset==DATASET) %>% 
    filter(num_clusters=="20000") %>%
    filter(cluster_approach %in% c("shiftDist-levFilter", "test-aa-clustered-v2")) %>%
    filter(model %in% c("ohe", "hyena_normalized", "hyena_unnormalized", "esm_normalized", "esm_unnormalized")) %>%
    select(dataset, filter, cluster_approach, metadata, model, accuracy, num_clusters) %>% 
    distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
    mutate(model=paste0(model, "_accuracy")) %>%
    pivot_wider(names_from = model, values_from = accuracy)
  
  accuracy_dfs[[i]] <- summ_accuracy
}
dev.off()

full_df <- bind_rows(accuracy_dfs)

full_df <- full_df %>% filter(!(metadata %in% c("date","geo_loc_name", "LibrarySelection", "reference", "studyaccession")))

write_tsv(full_df, file.path(opt$output_folder, paste0("all_datsets", "_model_comparisons_table_", system_time, ".tsv")), col_names = T, quote="needed")



full_df %>% group_by(dataset) %>% summarise(across(-c(filter:num_clusters), \(x) mean(x,na.rm=T))) %>%
  write_tsv(file.path(opt$output_folder, paste0("all_datsets_SUMMARY", "_model_comparisons_table_", system_time, ".tsv")), col_names = T, quote="needed")


full_df %>% group_by(dataset) %>% slice_max(ohe_accuracy+hyena_normalized_accuracy+esm_normalized_accuracy+esm_unnormalized_accuracy+hyena_unnormalized_accuracy,n=2) %>% 
  summarise(across(-c(filter:num_clusters), \(x) mean(x,na.rm=T))) %>% 
  write_tsv(file.path(opt$output_folder, paste0("all_datsets_SUMMARY_top2_per_dataset", "_model_comparisons_table_", system_time, ".tsv")), col_names = T, quote="needed")

full_df %>% group_by(dataset) %>% slice_min(ohe_accuracy+hyena_normalized_accuracy+esm_normalized_accuracy+esm_unnormalized_accuracy+hyena_unnormalized_accuracy,n=2) %>% 
  summarise(across(-c(filter:num_clusters), \(x) mean(x,na.rm=T))) %>% 
  write_tsv(file.path(opt$output_folder, paste0("all_datsets_SUMMARY_bottom2_per_dataset", "_model_comparisons_table_", system_time, ".tsv")), col_names = T, quote="needed")
