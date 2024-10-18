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
