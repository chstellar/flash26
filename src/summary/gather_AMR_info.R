# summary_plots.R
# Daniel Cotter 
# takes in a path to a dataset folder and produces a summary plot of the accuracies of 
# all of the models run on the dataset 

## Load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
library(ComplexHeatmap)

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
all_files <- fread(cmd="ls results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients_AMR-annotated.tsv",
                   header=FALSE) %>% pull(V1)

all_files <- c(all_files, fread(cmd="ls results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients_AMR-annotated.tsv",
                                header=FALSE) %>% pull(V1))

all_files <- all_files[str_detect(all_files, "eColi|eFaec|pneum|kleb")]


dt <- map(all_files, \(x,y) fread(x, sep = "\t") %>% mutate(path=x)) %>% bind_rows()


get_max_abs_value <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    max(abs(nums), na.rm = TRUE)
  })
}

data <- dt %>% filter(feature!="feature") %>% 
  mutate(coefficients=get_max_abs_value(coefficients)) %>% 
  dplyr::rename(max_coefficient=coefficients)


data <- data %>% distinct(path,metadata_category,feature,.keep_all=T) %>% mutate(AMR=grepl("MEG_|.fasta", anno)) %>% 
  select(path, metadata_category, accuracy, sensitivity, specificity, feature, max_coefficient, AMR, anno)

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

data <- data %>% relocate(num_clusters:cluster_approach, .after=metadata_category) %>% select(-paramater_set)


data <- data %>% 
  mutate(meg_entries = str_extract_all(
    anno,
    "(?<=\\\"\\\"\\\"\\\")[\\(\\)\\d\\w\\.\\-_]+(?=\\\"\\\"\\\"\\\")"
  ))

data <- data %>% select(-anno)

bacterial_data <- data %>% filter(dataset %in% c("eColi-arcadia-amr", "eFaecium-CollEtAl", "pneumo-ERP001505", "klebsiella-AMR-PRJEB42462"))

non_amr_metadata_categories <- c("reference", "studyaccession", "species", "virulence_score", "resistance_score", "num_resistance_classes",
                                 "num_resistance_genes", "date", "main_st", "major_wards", "neonate", "GPSC", "GPSC type", "Continent", "Predicted_PEN_MIC", 
                                 "Predicted_PEN_MIC_CLSI", "folA", "folP", "Predicted_Cotrimoxazole_susceptibility", "cat1", "Predicted_Chloramphenicol_susceptibility",
                                 "ermB1", "mef1", "Predicted_Erythromycin_susceptibility", "tet", "Predicted_Tetracycline_susceptibility",
                                 "No_of_nonsusceptible", "MDR", "PEN-SXT-CHL-ERY-TET")

bacterial_data <- bacterial_data %>% filter(!(metadata_category %in% non_amr_metadata_categories))

bacterial_data_non_pneumo <- bacterial_data %>% filter(model=="hyena_normalized", 
                                                       cluster_approach=="shiftDist-levFilter", 
                                                       filter=="filter1",
                                                       dataset != "pneumo-ERP001505") %>% 
  group_by(dataset, metadata_category) %>% arrange(desc(max_coefficient)) %>% 
  mutate(num_amr_features=sum(AMR),
         prop_amr_features=sum(AMR)/n())

bacterial_data_pneumo <- bacterial_data %>% filter(model=="hyena_unnormalized", 
                          cluster_approach=="masked-aa-clustered", 
                          filter=="filter1",
                          dataset=="pneumo-ERP001505") %>% 
  group_by(dataset, metadata_category) %>% arrange(desc(max_coefficient)) %>% 
  mutate(num_amr_features=sum(AMR),
         prop_amr_features=sum(AMR)/n())

bacterial_data <- bind_rows(bacterial_data_non_pneumo, bacterial_data_pneumo)

bacterial_data %>% distinct(dataset,metadata_category,.keep_all = T) %>%
  ggplot(aes(x=dataset,y=num_amr_features)) +
  geom_point(position=position_jitter(width=0.05, height = 0)) +
  geom_boxplot(color="red", fill=NA, outlier.alpha = 0) + 
  theme_pubclean() + ylab("Number of AMR features in Prediction") + 
  xlab("Dataset")

# bacterial_data %>% distinct(dataset,metadata_category,.keep_all = T) %>% 
#   filter(num_amr_features > 0) %>% 
#   ggplot(aes(x=num_amr_features, y=accuracy)) +
#   geom_point() + geom_smooth(method="lm") + 
#   theme_pubr()

# mechanism_data <- bacterial_data %>% rowwise() %>% mutate(meg_entries = paste(meg_entries, collapse="+")) %>%
#   separate_longer_delim(meg_entries, "+") %>%
#   separate(meg_entries, into=c("meg_id", "type", "class", "mechanism", "group"), sep="\\|", fill="right") %>% 
#   distinct(dataset, metadata_category, meg_id, .keep_all=T) %>% 
#   filter(nchar(meg_id)>2) %>% 
#   arrange(dataset,metadata_category) %>% 
#   group_by(dataset,mechanism) %>%
#   summarise(count=n())
# 
# mechanism_order <- mechanism_data %>%group_by(mechanism) %>% summarise(n=sum(count)) %>% filter(n>1) %>% arrange(desc(n))
# 
# mechanism_data %>% filter(mechanism %in% mechanism_order$mechanism) %>% 
#   mutate(mechanism=factor(mechanism, levels=mechanism_order$mechanism)) %>%
#   ggplot(aes(x=mechanism,y=count,fill=dataset)) + geom_col() + theme_pubr() +
#   theme(axis.text.x = element_text(angle=45,hjust=1),
#         plot.margin = margin(t=0.1,r=0.1,b=0.1,l=2,unit="cm"),
#         legend.position = c(0.8,0.7)) +
#   ylab("Distinct occurrences\nper metadata category") + 
#   xlab("MEGARes Mechanism")
# 
# 
# class_data <- bacterial_data %>% rowwise() %>% mutate(meg_entries = paste(meg_entries, collapse="+")) %>%
#   separate_longer_delim(meg_entries, "+") %>%
#   separate(meg_entries, into=c("meg_id", "type", "class", "mechanism", "group"), sep="\\|", fill="right") %>% 
#   distinct(dataset, metadata_category, meg_id, .keep_all=T) %>% 
#   filter(accuracy>=0.6) %>%
#   filter(nchar(meg_id)>2) %>% 
#   arrange(dataset,metadata_category) %>% 
#   group_by(dataset,class) %>%
#   summarise(count=n())
# 
# class_order <- class_data %>%group_by(class) %>% summarise(n=sum(count)) %>% filter(n>1) %>% arrange(desc(n))
# 
# class_data %>% filter(class %in% class_order$class) %>% 
#   mutate(class=factor(class, levels=class_order$class)) %>%
#   mutate(dataset=str_extract(dataset, "\\w+")) %>%
#   ggplot(aes(x=class,y=count,fill=dataset)) + geom_col() + theme_pubr() +
#   theme(axis.text.x = element_text(angle=45,hjust=1),
#         plot.margin = margin(t=0.1,r=0.1,b=0.1,l=2,unit="cm"),
#         legend.position = c(0.8,0.7)) +
#   ylab("Distinct occurrences\nper metadata category") + 
#   xlab("MEGARes class") +
#   theme(axis.text = element_text(size=8), axis.title = element_text(size=12))
# 
# ggsave("results/summary/all_amr_annotations_geq60percent.pdf", width=6, height=4.5)

col = "mechanism"

bacterial_data <- bacterial_data %>% rowwise() %>% mutate(meg_entries = paste(meg_entries, collapse="+") %>% unlist()) %>%
  separate_longer_delim(meg_entries, "+") %>%
  mutate(meg_entries=str_replace_all(meg_entries,"\\(|\\)|MEG_\\d+\\.|\\.fasta|\\\"\\\"\\\"\\\"", "")) %>%
  separate(meg_entries, into=c("type", "class", "mechanism", "group", "info"), sep="\\.", fill="right") %>% 
  filter(AMR)

bact_summ <- bacterial_data %>%
  distinct(dataset, metadata_category, feature, !!as.name(col), .keep_all=T) %>% 
  arrange(dataset,metadata_category) %>%
  group_by(dataset,metadata_category,!!as.name(col)) %>%
  summarise(count=n())

bact_summ <- bact_summ %>%  mutate(dataset = c("eColi-arcadia-amr"="E. coli",
                                               "eFaecium-CollEtAl"="E. faecium", 
                                               "klebsiella-AMR-PRJEB42462"="Klebsiella pneumoniae", 
                                               "pneumo-ERP001505"= "S. pneumoniae")[dataset])

dataset_mapping <- bact_summ %>% distinct(dataset,metadata_category) %>% select(metadata_category, dataset) %>% deframe()

bact_mat <- bact_summ %>% ungroup() %>% select(-dataset) %>% 
  pivot_wider(id_cols=!!col,names_from="metadata_category", values_from="count") %>% 
  mutate(across(everything(), ~replace_na(.,-1))) %>%
  column_to_rownames(col) %>% as.matrix()

heatmap_matrix <- bact_mat

column_names = colnames(heatmap_matrix)
column_names = tolower(column_names)
column_names = str_remove(column_names, "_ris|phenotypic_|_\\d+.+|_\\d+")
column_names = str_replace(column_names, "_|/", ".")
names(column_names) = colnames(heatmap_matrix)
column_names = str_to_sentence(column_names)

# Generate a color palette dynamically
unique_datasets <- unique(dataset_mapping)

colors <- RColorBrewer::brewer.pal(length(unique_datasets), "Set1") # You can choose other palettes as needed

names(colors) <- unique_datasets # Name the colors according to the datasets

#heatmap_matrix <- heatmap_matrix[rowSums(heatmap_matrix)>1,]

row_mapping <- data.frame(mechanism=rownames(heatmap_matrix)) %>% left_join(bacterial_data %>% distinct(class, mechanism))

row_order <- row_mapping %>% arrange(class) %>% pull(mechanism)

heatmap_matrix <- heatmap_matrix[row_order,]

# Create an annotation object

annotation_df <- data.frame(Dataset = dataset_mapping[colnames(heatmap_matrix)])

ha <- HeatmapAnnotation(name = "Species", Dataset = annotation_df$Dataset, 
                        
                        col = list(Dataset = colors))

unique(row_mapping$class)

p <- Heatmap(heatmap_matrix, 
        
        name = "Number of\nFeature Clusters", # Name of the heatmap legend
        
        col = circlize::colorRamp2(c(-1, 0, 5), c("white", "red", "darkred")),# Color scale
        
        show_row_names = TRUE, # Show gene names
        
        show_column_names = TRUE, # Show annotation labels
        row_names_side = "left",
        row_dend_side = "right",
        rect_gp = gpar(col = "black", lwd = 0.5),
        column_names_side = "top",
        column_labels = column_names,
        column_dend_side = "bottom",
        cluster_rows = TRUE, # Cluster rows (genes)
        cluster_columns = TRUE, # Cluster columns (annotations)
        heatmap_legend_param = list(at = c(-1, 5, 10), # Specify the positions for the labels
                                    labels = c("0", "5", ">10")),
        top_annotation = ha) # Legend title

png("/oak/stanford/groups/horence/dcotter1/FLASH_paper_figures/heatmap/amr_metadata_megares_lookup_table.png", height=15, width=15, res=400,units = "in")
draw(p, padding=unit(c(10,150,20,10), "pt"))
dev.off()


