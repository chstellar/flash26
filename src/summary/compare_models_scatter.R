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

if (is.null(opt$dataset_path) || is.null(opt$output) || is.null(opt$name)) {
  stop("Please provide a dataset path, output path, and dataset name")
}

## Load data
# strip tailing slash from dataset path if it exists
library(ggpubr)

## plot all data in one plot
load_cmd1 <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients.tsv")
data1 <- fread(cmd=load_cmd1a, header=FALSE, sep="\t", 
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
                                          'hyenaMarlowe_normalized|hyenaMarlowe_unnormalized|esm_normalized|esm_unnormalized|hyena_normalized|hyena_unnormalized|ohe')) %>% arrange(desc(model))

data <- data %>% mutate(filter = str_extract(paramater_set, "(filter\\d)_",group=1)) %>% 
  mutate(cluster_approach = str_extract(paramater_set, "filter\\d_([A-Za-z-]+)_", group=1))

data %>% filter(filter =="filter1" ) %>% filter(cluster_approach=="shiftDist-keepTopES") %>% filter(num_clusters=="10000") %>% ggplot(aes(x=dataset, y=specificity)) + geom_boxplot(aes(fill=model)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_brewer(type="qual") + scale_fill_brewer(type="qual")

data %>% distinct(dataset) %>% arrange(dataset) -> unique_datasets

DATASET="eFaecium-CollEtAl"
DATASET="eColi-arcadia-amr"
plot_data <- data %>% 
  filter(dataset==DATASET) %>% 
  select(filter, cluster_approach, metadata, model, accuracy, num_clusters) %>% 
  distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
  mutate(model = paste0(model, " (", num_clusters, " clusters)")) %>% select(-num_clusters) %>% 
  pivot_wider(names_from = model, values_from = accuracy)
all_cols <- colnames(plot_data %>% select(-c(filter, cluster_approach, metadata, `ohe (50000 clusters)`)))
pdf("results/summary_plots/eColi_accuracy_comparison_using_adelie_results_241121.pdf")
for (i in all_cols) {
  y_col = sym(i)
  p <- plot_data %>% 
    ggplot(aes(x=`ohe (50000 clusters)`,y=!!y_col)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
    geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0.4,1), xlim=c(0.4,1)) +
    ggrepel::geom_text_repel(aes(label=metadata)) +
    ggtitle("Accuracy")
  print(p)
}
dev.off()

data %>% 
  filter(dataset==DATASET) %>% 
  select(filter, cluster_approach, metadata, model, accuracy, sensitivity) %>% 
  distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% 
  pivot_wider(names_from = model, values_from = c(accuracy, sensitivity), names_glue = "{model}_{.value}") %>% write_tsv("/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline/results/summary_plots/eFaecium_comparison_data_241119.tsv", col_names = T)

# plot genome data 
load_cmd="grep cluster /oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline/results/*/filter*/*/esm/genomes/normalized/*_glmnet_genomes_results_*_nonzero_coefficients.tsv"
genome_data <- fread(cmd=load_cmd, header=FALSE, sep="\t", 
                     col.names = c("path_metadata", "feature", "accuracy",
                                   "sensitivity", "specificity", "classes", "coefficients"))

genome_data <- genome_data %>% 
  distinct(path_metadata, .keep_all = T) %>%
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity)

genome_data %>% ggplot(aes(x=metadata, y=accuracy)) + geom_point(position="jitter") + 
  geom_boxplot(aes(fill=metadata)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


# read in annotations
load_cmd="grep cluster_ /oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline/results/*/filter*/*/esm/genomes/normalized/*_glmnet_genomes_results_*_nonzero_coefficients_annotated.tsv"
genome_data <- fread(cmd=load_cmd, header=FALSE, sep="\t", 
                     col.names = c("path_metadata", "feature", "accuracy",
                                   "sensitivity", "specificity", "classes", "coefficients", "cluster", "annotations"))
genome_data <- genome_data %>% 
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity, annotations)

# grab annotation data
setwd("/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline")
load_cmd <- paste0("grep cluster_ ", "results/*/filter*/*/*/*malized/*_nonzero_coefficients_annotated.tsv")
data <- fread(cmd=load_cmd, header=FALSE, sep="\t", 
              col.names = c("path_metadata", "feature", "accuracy",
                            "sensitivity", "specificity", "classes", "coefficients", "cluster", "anno"))

data <- data %>%  separate(path_metadata, into=c("path", "metadata"), sep=":") %>% mutate(path=dirname(gsub("^results/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  dplyr::rename(paramater_set=path)
data <- data %>% mutate(paramater_set=str_replace(paramater_set, "_", "/")) %>% 
  separate(paramater_set, into=c("dataset", "paramater_set"), sep="/")
data <- data %>% mutate(model=ifelse(grepl("esm_normalized", paramater_set), yes="esm_normalized", no=
                                       ifelse(grepl("esm_unnormalized", paramater_set), yes="esm_unnormalized",
                                              ifelse(grepl("hyena", paramater_set), yes="hyena", no="ohe"))))

data %>% filter(accuracy>0.75) %>% filter(dataset=="eFaecium-CollEtAl") %>% 
  filter(grepl("anti|IS|mge|ice", anno, ignore.case=T)) %>% 
  distinct(cluster, .keep_all = T) %>% 
  mutate(seqs = sapply(str_extract_all(anno, "[ACTGN]{54}"), toString)) %>% 
  mutate(anno = sapply(str_extract_all(anno, "antibiotic|IS|ICE|ice|mge|MGE|Antibiotic"), \(x) toString(unique(x)))) %>% 
  select(seqs, anno) %>% separate_longer_delim(seqs, delim=", ") %>% write_tsv("logs/summary_annos.tsv", quote="needed", col_names=T)

#######
load_cmd1a <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients.tsv")
data1a <- fread(cmd=load_cmd1a, header=FALSE, sep="\t", 
                col.names = c("path_metadata", "feature", "accuracy",
                              "sensitivity", "specificity", "classes", "coefficients"))
load_cmd1b <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*glmnet*_nonzero_coefficients.tsv")
data1b <- fread(cmd=load_cmd1b, header=FALSE, sep="\t", 
                col.names = c("path_metadata", "feature", "accuracy",
                              "sensitivity", "specificity", "classes", "coefficients"))

data1a$measure <- "adelie"
data1b$measure <- "glmnet"
data<- rbind(data1a,data1b)
data <- data %>% 
  distinct(path_metadata, .keep_all = T) %>%
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity, measure)

data <- data %>% mutate(num_clusters = str_extract(path, "top(\\d+)", group=1))
data <- data %>% mutate(path=dirname(gsub("^results/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  dplyr::rename(paramater_set=path)
data <- data %>% mutate(paramater_set=str_replace(paramater_set, "_", "/")) %>% 
  separate(paramater_set, into=c("dataset", "paramater_set"), sep="/")
data <- data %>% mutate(model=str_extract(paramater_set, 
                                          'hyenaMarlowe_normalized|hyenaMarlowe_unnormalized|esm_normalized|esm_unnormalized|hyena_normalized|hyena_unnormalized|ohe')) %>% arrange(desc(model))

data <- data %>% mutate(filter = str_extract(paramater_set, "(filter\\d)_",group=1)) %>% 
  mutate(cluster_approach = str_extract(paramater_set, "filter\\d_([A-Za-z-]+)_", group=1))

data <- data %>% filter(as.integer(num_clusters)>=10000)
DATASET="eFaecium-CollEtAl"
#DATASET="eColi-arcadia-amr"
plot_data <- data %>% 
  filter(dataset==DATASET) %>% 
  select(filter, cluster_approach, metadata, model, specificity, num_clusters, measure) %>% 
  pivot_wider(names_from = measure, values_from = specificity)

models <- plot_data %>% filter(!is.na(adelie)) %>% filter(!is.na(glmnet)) %>% pull(model) %>% unique()

pdf("results/summary_plots/eFaecium_specificity_adelie_vs_glmnet_comparison_241120.pdf")
for (i in models) {
  p <- plot_data %>% filter(model==i) %>%
    ggplot(aes(x=glmnet,y=adelie)) + geom_abline(slope=1, intercept = 0, color="black", lty="dashed") + 
    geom_point(aes(shape=filter, color=cluster_approach), size=2) + theme_minimal() + coord_cartesian(ylim=c(0.4,1), xlim=c(0.4,1)) +
    ggrepel::geom_text_repel(aes(label=metadata)) +
    ggtitle(paste("Specificity (10k clusters):", i))
  print(p)
}
dev.off()


## compare e faecium results to the lancet 
lancet_data <- fread("/oak/stanford/groups/horence/dcotter1/utility_files/lancet_RIS_predictions.csv")
lancet_data <- lancet_data %>% mutate(metadata_category=paste0(str_to_lower(metadata_category), "_RIS")) %>% rename(metadata=metadata_category)
comp_data <-data %>% filter(dataset =="eFaecium-CollEtAl") %>% left_join(lancet_data)

comp_data %>% ggplot(aes(x=lancet_accuracy, y=accuracy)) + geom_abline(intercept=0, slope=1, color="black", lty="dashed") +
  geom_point(aes(color=filter, shape=model)) + theme_minimal() +
  coord_cartesian(xlim=c(0.6,1), ylim=c(0.6,1)) + ggtitle("Accuracy") + geom_text(y=0.95,aes(x=lancet_accuracy, label=metadata), angle=90, size=3)
comp_data %>% ggplot(aes(x=lancet_sensitivity, y=sensitivity)) + geom_abline(intercept=0, slope=1, color="black", lty="dashed") +
  geom_point(aes(color=filter, shape=model)) + theme_minimal() +
  coord_cartesian(xlim=c(0.6,1), ylim=c(0.6,1)) + ggtitle("Sensitivity") + geom_text(y=0.95,aes(x=lancet_sensitivity, label=metadata), angle=90, size=3)
comp_data %>% ggplot(aes(x=lancet_specificity, y=specificity)) + geom_abline(intercept=0, slope=1, color="black", lty="dashed") +
  geom_point(aes(color=filter, shape=model)) + theme_minimal() +
  coord_cartesian(xlim=c(0.6,1), ylim=c(0.6,1)) + ggtitle("Specificity") + geom_text(y=0.95,aes(x=lancet_specificity, label=metadata), angle=90, size=3)

