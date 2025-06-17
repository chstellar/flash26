library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(furrr)

plan(multisession, workers=8)

# Define Function 
get_max_abs_value <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    max(abs(nums), na.rm = TRUE)
  })
}

all_files <- fread(cmd="ls results/*/filter1/shiftDist-levFilter/hyena/normalized/*_hyena_adelie_results_top20000_k54_s54_nonzero_coefficients_blastp_annotated.tsv",
                   header=FALSE) %>% pull(V1)

all_files <- all_files[!str_detect(all_files, "pneumo")]

all_files2 <- fread(cmd="ls results/pneumo*/filter1/masked-aa-clustered/hyena/unnormalized/*_hyena_adelie_results_top20000_k54_s54_nonzero_coefficients_blastp_annotated.tsv",
                    header=FALSE) %>% pull(V1)

all_files <- c(all_files, all_files2)

datasets <- str_extract(all_files, "results/(.+)/filter", group=1)

dt <- map2(all_files, datasets, \(x,y) fread(x, sep = "\t") %>% mutate(dataset=y)) %>% bind_rows()
#dt <- fread("results/eFaecium-CollEtAl/filter1/shiftDist-levFilter/hyena/normalized/eFaecium-CollEtAl_hyena_adelie_results_top20000_k54_s54_nonzero_coefficients_blastp_annotated.tsv")

dt_mod <- dt %>% group_by(dataset, metadata_category) %>% mutate(max_coefficient=get_max_abs_value(coefficients)) %>% 
  arrange(dataset, metadata_category, -max_coefficient) %>% 
  mutate(annotation = str_remove_all(stitle, "\\[.+\\]$|MULTISPECIES:\\s|, partial"))

non_amr_metadata_categories <- c("reference", "studyaccession", "species", "virulence_score", "resistance_score", "num_resistance_classes",
                                 "num_resistance_genes", "date", "main_st", "major_wards", "neonate", "GPSC", "GPSC type", "Continent", "Predicted_PEN_MIC", 
                                 "Predicted_PEN_MIC_CLSI", "folA", "folP", "Predicted_Cotrimoxazole_susceptibility", "cat1", "Predicted_Chloramphenicol_susceptibility",
                                 "ermB1", "mef1", "Predicted_Erythromycin_susceptibility", "tet", "Predicted_Tetracycline_susceptibility",
                                 "No_of_nonsusceptible", "MDR", "PEN-SXT-CHL-ERY-TET")

dt_amr <- dt_mod %>% filter(dataset %in% c("eFaecium-CollEtAl", "pneumo-ERP001505", "klebsiella-AMR-PRJEB42462", "eColi-arcadia-amr")) %>%
  filter(!(metadata_category %in% non_amr_metadata_categories))

dt_amr <- dt_amr %>% group_by(dataset, metadata_category, feature) %>% distinct(annotation, .keep_all=T) %>% ungroup() %>%
  mutate(annotation = tolower(annotation)) %>%
  mutate(annotation = str_replace(annotation, "n\\(6\\)", "n-6")) %>%
  mutate(annotation = str_remove_all(annotation, "\\(|\\)")) %>%
  mutate(annotation_label = str_extract(annotation, "penicillin[ -]binding|IS\\w+|\\w+[ -]transporter|ion channel|transposase|transcriptional regulator|^.+domain[ -]containing[ -]protein|[\\w-]+[ -]\\w+ase|[\\w-]+[ -]family[ -]protein|.+-like protein|\\w+[ -]\\w+[ -]protein|\\w+[ -]protein")) %>%
  mutate(annotation_label = str_remove(annotation_label, " domain-containing protein| family protein")) %>%
  mutate(annotation_label = str_remove(annotation_label, "^family")) %>%
  mutate(annotation_label = str_trim(annotation_label)) %>% 
  mutate(annotation_label = ifelse(str_detect(annotation_label, "n-6-methyltransferase"), "n-6-methyltransferase", annotation_label)) %>%
  filter(!annotation_label == "f")

dt_amr_counts <- dt_amr %>% mutate(annotation_label=ifelse(is.na(annotation_label), annotation, annotation_label)) %>% 
  filter(!is.na(annotation_label)) %>% distinct(dataset, metadata_category, feature, annotation_label, .keep_all=T) %>% 
  select(dataset, metadata_category, feature, max_coefficient, annotation_label) %>% group_by(metadata_category, annotation_label) %>% 
  summarize(sum=sum(max_coefficient), n_features=n_distinct(feature)) %>% ungroup()

sum_mat <- dt_amr_counts %>% select(-n_features) %>% group_by(annotation_label) %>%
  filter(n() > 4) %>% ungroup() %>% 
  pivot_wider(names_from="annotation_label", values_from="sum") %>% 
  mutate(across(everything(), \(x) replace_na(x, 0))) %>%
  column_to_rownames("metadata_category")

count_mat <- dt_amr_counts %>% select(-sum) %>% group_by(annotation_label) %>%
  filter(n() > 4) %>% ungroup() %>% 
  pivot_wider(names_from="annotation_label", values_from="n_features") %>% 
  mutate(across(everything(), \(x) replace_na(x, 0))) %>%
  column_to_rownames("metadata_category")

heatmap_matrix <- as.matrix(t(count_mat))


# Create the heatmap

# Create an annotation object

dataset_mapping <- dt_amr %>% filter(metadata_category %in% unique(dt_amr_counts$metadata_category)) %>% distinct(dataset, metadata_category) %>% select(metadata_category, dataset) %>% deframe()



# Generate a color palette dynamically
unique_datasets <- unique(dataset_mapping)

colors <- RColorBrewer::brewer.pal(length(unique_datasets), "Set1") # You can choose other palettes as needed

names(colors) <- unique_datasets # Name the colors according to the datasets

ordered_columns <- order(colnames(heatmap_matrix))

heatmap_matrix <- heatmap_matrix[, ordered_columns]

column_names = colnames(heatmap_matrix)

column_names = tolower(column_names)
column_names = str_remove(column_names, "_ris|phenotypic_|_\\d+.+|_\\d+")
column_names = str_replace(column_names, "_|/", ".")
names(column_names) = colnames(heatmap_matrix)

ordered_columns <- order(column_names)

heatmap_matrix = heatmap_matrix[, ordered_columns]

# Create an annotation object

annotation_df <- data.frame(Dataset = dataset_mapping[colnames(heatmap_matrix)])

ha <- HeatmapAnnotation(Dataset = annotation_df$Dataset, 
                        
                        col = list(Dataset = colors))


p <- Heatmap(heatmap_matrix, 
        
        name = "Gene Expression", # Name of the heatmap legend
        
        col = circlize::colorRamp2(c(0, 4), c("white", "red")),# Color scale
        
        show_row_names = TRUE, # Show gene names
        
        show_column_names = TRUE, # Show annotation labels
        column_labels = column_names,
        row_names_side = "left",
        row_dend_side = "right",
        rect_gp = gpar(col = "black", lwd = 0.5),
        column_names_side = "top",
        column_dend_side = "bottom",
        cluster_rows = TRUE, # Cluster rows (genes)
        cluster_columns = FALSE, # Cluster columns (annotations)
        heatmap_legend_param = list(title = "Count Different \nNonzero Coefficients",
                                    at = c(0, 4), # Specify the positions for the labels
                                    labels = c("0", ">4")),
        top_annotation = ha) # Legend title

pdf("/oak/stanford/groups/horence/dcotter1/FLASH_paper_figures/heatmap/all_amr_heatmap.pdf", height=18, width=15)
draw(p, padding=unit(c(10,40,10,10), "pt"))
dev.off()
