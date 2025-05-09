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
load_cmd1 <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients_AMR-annotated.tsv")
data1 <- fread(cmd=load_cmd1, header=FALSE, sep="\t", 
               col.names = c("path_metadata", "feature", "accuracy",
                             "sensitivity", "specificity", "classes", "coefficients", 
                             "cluster", "anno"))

load_cmd2 <- paste0("grep cluster ", "results/*/filter*/*/*/*_nonzero_coefficients_AMR-annotated.tsv")
data2 <- fread(cmd=load_cmd2, header=FALSE, sep="\t", 
               col.names = c("path_metadata", "feature", "accuracy",
                             "sensitivity", "specificity", "classes", "coefficients", 
                             "cluster", "anno"))


get_max_abs_value <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    max(abs(nums), na.rm = TRUE)
  })
}

data <- rbind(data1,data2) %>% filter(feature!="feature") %>% 
  mutate(coefficients=get_max_abs_value(coefficients)) %>% 
  rename(max_coefficient=coefficients)


data <- data %>% distinct(path_metadata,feature,.keep_all=T) %>% mutate(AMR=grepl("megares|MEG_", anno)) %>% 
  separate(path_metadata, into=c("path", "metadata"), sep=":") %>%
  select(path, metadata, accuracy, sensitivity, specificity, feature, max_coefficient, AMR, anno)

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

data <- data %>% relocate(num_clusters:cluster_approach, .after=metadata) %>% select(-paramater_set)


data <- data %>%
  mutate(meg_entries = str_extract_all(
    anno,
    "MEG_\\d+\\|[\\w-]+\\|[\\w-]+\\|[\\w-]+\\|[\\w-]+"
  ))

data <- data %>% select(-anno)

bacterial_data <- data %>% filter(dataset %in% c("eColi-arcadia-amr", "eFaecium-CollEtAl", "pneumo-ERP001505", "klebsiella-AMR-PRJEB42462"))

non_amr_metadata_categories <- c("reference", "studyaccession", "species", "virulence_score", "resistance_score", "num_resistance_classes",
                                 "num_resistance_genes", "date", "main_st", "major_wards", "neonate", "GPSC", "GPSC type", "Continent", "Predicted_PEN_MIC", 
                                 "Predicted_PEN_MIC_CLSI", "folA", "folP", "Predicted_Cotrimoxazole_susceptibility", "cat1", "Predicted_Chloramphenicol_susceptibility",
                                 "ermB1", "mef1", "Predicted_Erythromycin_susceptibility", "tet", "Predicted_Tetracycline_susceptibility",
                                 "No_of_nonsusceptible", "MDR")

bacterial_data <- bacterial_data %>% filter(!(metadata %in% non_amr_metadata_categories))

bacterial_data <- bacterial_data %>% filter(model=="hyena_unnormalized", 
                          cluster_approach=="masked-aa-clustered", 
                          filter=="filter1") %>% 
  group_by(dataset, metadata) %>% arrange(desc(max_coefficient)) %>% 
  mutate(num_amr_features=sum(AMR),
         prop_amr_features=sum(AMR)/n())

bacterial_data %>% distinct(dataset,metadata,.keep_all = T) %>%
  ggplot(aes(x=dataset,y=num_amr_features)) +
  geom_point(position=position_jitter(width=0.05, height = 0)) +
  geom_boxplot(color="red", fill=NA, outlier.alpha = 0) + 
  theme_pubclean() + ylab("Number of AMR features in Prediction") + 
  xlab("Dataset")

bacterial_data %>% distinct(dataset,metadata,.keep_all = T) %>% 
  filter(num_amr_features > 0) %>% 
  ggplot(aes(x=num_amr_features, y=accuracy)) +
  geom_point() + geom_smooth(method="lm") + 
  theme_pubr()

mechanism_data <- bacterial_data %>% rowwise() %>% mutate(meg_entries = paste(meg_entries, collapse="+")) %>%
  separate_longer_delim(meg_entries, "+") %>%
  separate(meg_entries, into=c("meg_id", "type", "class", "mechanism", "group"), sep="\\|", fill="right") %>% 
  distinct(dataset, metadata, meg_id, .keep_all=T) %>% 
  filter(nchar(meg_id)>2) %>% 
  arrange(dataset,metadata) %>% 
  group_by(dataset,mechanism) %>%
  summarise(count=n())

mechanism_order <- mechanism_data %>%group_by(mechanism) %>% summarise(n=sum(count)) %>% filter(n>1) %>% arrange(desc(n))

mechanism_data %>% filter(mechanism %in% mechanism_order$mechanism) %>% 
  mutate(mechanism=factor(mechanism, levels=mechanism_order$mechanism)) %>%
  ggplot(aes(x=mechanism,y=count,fill=dataset)) + geom_col() + theme_pubr() +
  theme(axis.text.x = element_text(angle=45,hjust=1),
        plot.margin = margin(t=0.1,r=0.1,b=0.1,l=2,unit="cm"),
        legend.position = c(0.8,0.7)) +
  ylab("Distinct occurrences\nper metadata category") + 
  xlab("MEGARes Mechanism")


class_data <- bacterial_data %>% rowwise() %>% mutate(meg_entries = paste(meg_entries, collapse="+")) %>%
  separate_longer_delim(meg_entries, "+") %>%
  separate(meg_entries, into=c("meg_id", "type", "class", "mechanism", "group"), sep="\\|", fill="right") %>% 
  distinct(dataset, metadata, meg_id, .keep_all=T) %>% 
  filter(accuracy>=0.6) %>%
  filter(nchar(meg_id)>2) %>% 
  arrange(dataset,metadata) %>% 
  group_by(dataset,class) %>%
  summarise(count=n())

class_order <- class_data %>%group_by(class) %>% summarise(n=sum(count)) %>% filter(n>1) %>% arrange(desc(n))

class_data %>% filter(class %in% class_order$class) %>% 
  mutate(class=factor(class, levels=class_order$class)) %>%
  mutate(dataset=str_extract(dataset, "\\w+")) %>%
  ggplot(aes(x=class,y=count,fill=dataset)) + geom_col() + theme_pubr() +
  theme(axis.text.x = element_text(angle=45,hjust=1),
        plot.margin = margin(t=0.1,r=0.1,b=0.1,l=2,unit="cm"),
        legend.position = c(0.8,0.7)) +
  ylab("Distinct occurrences\nper metadata category") + 
  xlab("MEGARes class") +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=12))

ggsave("results/summary/all_amr_annotations_geq60percent.pdf", width=6, height=4.5)

  
