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
  make_option(c("-d", "--dataset_path"), type="character", default=NULL, 
              help="Path to the dataset folder"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to the output folder"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="Name of the dataset")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$dataset_path) || is.null(opt$output) || is.null(opt$name)) {
  stop("Please provide a dataset path, output path, and dataset name")
}

## Load data
# strip tailing slash from dataset path if it exists
library(ggpubr)
dataset_path = gsub("/$", "", opt$dataset_path)

# generate the glob pattern for the dataset
load_cmd <- paste0("grep cluster ", opt$dataset_path, "/filter*/*/*/*malized/*_nonzero_coefficients.tsv")
data <- fread(cmd=load_cmd, header=FALSE, sep="\t", 
              col.names = c("path_metadata", "feature", "accuracy",
                            "sensitivity", "specificity", "classes", "coefficients"))

data <- data %>% 
  distinct(path_metadata, .keep_all = T) %>%
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity)

# extract the dataset name from the path
data <- data %>% mutate(path=dirname(gsub("^results/[A-Za-z0-9-]+/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  rename(paramater_set=path)

# create a plot for each metadata category where each point refers to the parameter set
# facet by metadata category all in one column
p <- data %>% ggplot(aes(x=paramater_set, y=accuracy)) +
  # add a line to each plot at the maximum accuracy for that metadata category
  geom_hline(data=data %>% group_by(metadata) %>% summarise(max_accuracy=max(accuracy)), 
             aes(yintercept=max_accuracy), linetype="dashed", color="red") +
  # add a line to each plot at the mean accuracy for that metadata category
  geom_hline(data=data %>% group_by(metadata) %>% summarise(mean_accuracy=mean(accuracy)), 
             aes(yintercept=mean_accuracy), linetype="dotted", color="blue") +
  geom_point() +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  facet_wrap(~metadata, ncol=1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title=paste("Accuracy of models on", opt$name, "dataset"),
       x="Parameter set", y="Accuracy") +
  theme(plot.margin=margin(t=1,b=1,l=3,r=1, unit="cm"))

# save the plot (with height scaled to the number of metadata categories)
num_metadata <- data %>% distinct(metadata) %>% nrow()
ggsave(opt$output, height=ifelse(num_metadata<10, 10, 2*num_metadata), width=10, plot=p)



## plot all data in one plot
load_cmd1a <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients.tsv")
data1a <- fread(cmd=load_cmd1a, header=FALSE, sep="\t", 
              col.names = c("path_metadata", "feature", "accuracy",
                            "sensitivity", "specificity", "classes", "coefficients"))
load_cmd1b <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*glmnet*_nonzero_coefficients.tsv")
data1b <- fread(cmd=load_cmd1b, header=FALSE, sep="\t", 
                col.names = c("path_metadata", "feature", "accuracy",
                              "sensitivity", "specificity", "classes", "coefficients"))

load_cmd2 <- paste0("grep cluster ", "results/*/filter*/*/*/*_nonzero_coefficients.tsv")
data2 <- fread(cmd=load_cmd2, header=FALSE, sep="\t", 
              col.names = c("path_metadata", "feature", "accuracy",
                            "sensitivity", "specificity", "classes", "coefficients"))

load_cmd3 <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*_important_features.tsv")
data3 <- fread(cmd=load_cmd3, header=FALSE, sep="\t", 
               col.names = c("path_metadata", "feature", "accuracy",
                             "sensitivity", "specificity", "Gini"))

data1 <- data1a
data <- rbind(data1,data2)
data <- data %>% 
  distinct(path_metadata, .keep_all = T) %>%
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity)

data_rf <- data3 %>% distinct(path_metadata, .keep_all = T) %>%
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity)

data$measure <- "glmnet"
data_rf$measure <- "randomForests"

data <- rbind(data, data_rf)

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

data %>% ggplot(aes(x=dataset, y=specificity)) + geom_boxplot(aes(fill=filter)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_brewer(type="qual") + scale_fill_brewer(type="qual")

data %>% distinct(dataset, metadata) -> unique_sets

pdf("all_summary_plots_with_random_forests.pdf")
for (i in 1:nrow(unique_sets)) {
  temp_data <- data %>% filter(dataset==unique_sets[i,]$dataset, metadata==unique_sets[i,]$metadata)
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 10, paste(unique_sets[i,]$dataset), cex=1.5, font=2)
  text(5,5, paste(unique_sets[i,]$metadata), cex=1.2, font=2)
  p1 <- temp_data %>% ggplot(aes(x=paramater_set, y=accuracy, shape=filter, color=model)) +
    geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + coord_cartesian(ylim=c(0,1)) + ggtitle(unique_sets[i,]$metadata)
  p2 <- temp_data %>% ggplot(aes(x=paramater_set, y=accuracy, shape=filter, color=cluster_approach)) + 
    geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0,1)) + ggtitle(unique_sets[i,]$metadata)
  # p3 <- temp_data %>% ggplot(aes(x=paramater_set, y=specificity, shape=filter, color=model)) + 
  #   geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + coord_cartesian(ylim=c(0,1)) + ggtitle(unique_sets[i,]$metadata)
  # p4 <- temp_data %>% ggplot(aes(x=paramater_set, y=specificity, shape=filter, color=cluster_approach)) + 
  #   geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0,1)) + ggtitle(unique_sets[i,]$metadata)
  # p5 <- temp_data %>% ggplot(aes(x=paramater_set, y=sensitivity, shape=filter, color=model)) + 
  #   geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + coord_cartesian(ylim=c(0,1)) + ggtitle(unique_sets[i,]$metadata)
  # p6 <- temp_data %>% ggplot(aes(x=paramater_set, y=sensitivity, shape=filter, color=cluster_approach)) + 
  #   geom_point(size=3) + theme_minimal() + theme(axis.text.x = element_blank()) + coord_cartesian(ylim=c(0,1)) + ggtitle(unique_sets[i,]$metadata)
  p7 <- temp_data %>% ggplot(aes(x=filter, y=accuracy, fill=measure)) +
    geom_boxplot(position="dodge", color="black") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) + coord_cartesian(ylim=c(0,1)) + ggtitle(unique_sets[i,]$metadata)
  print(p1)
  print(p2)
  # print(p3)
  # print(p4)
  # print(p5)
  # print(p6)
  print(p7)
}
dev.off()
write_tsv(data, "all_summary_data_with_randForests.tsv", col_names = T, quote="needed")

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

