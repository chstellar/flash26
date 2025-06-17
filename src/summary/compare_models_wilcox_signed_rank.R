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


data_files <- list.files("results", recursive = T, "nonzero_coefficients.tsv", full.names = T)

all_data <- map(data_files, \(x) fread(x) %>%mutate(path=x), .progress=T)
discard(all_data, \(x) nrow(x) == 0) -> all_data

# ## plot all data in one plot
# load_cmd1 <- paste0("grep cluster ", "results/*/filter*/*/*/*malized/*adelie*_nonzero_coefficients.tsv")
# data1 <- fread(cmd=load_cmd1, header=FALSE, sep="\t", 
#                col.names = c("path_metadata", "feature", "accuracy",
#                              "sensitivity", "specificity", "classes", "coefficients"))
# 
# load_cmd2 <- paste0("grep cluster ", "results/*/filter*/*/*/*_nonzero_coefficients.tsv")
# data2 <- fread(cmd=load_cmd2, header=FALSE, sep="\t", 
#                col.names = c("path_metadata", "feature", "accuracy",
#                              "sensitivity", "specificity", "classes", "coefficients"))


# data <- rbind(data1,data2)

data <- bind_rows(all_data)
data <- data %>% 
  distinct(path, metadata_category, .keep_all = T) %>%
  dplyr::rename(metadata=metadata_category) %>%
  select(path, metadata, classes, accuracy, sensitivity, specificity) %>%
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

data <- data %>% rowwise() %>% 
  mutate(classes = length(str_split(gsub("\\[|\\]", "", classes), ",") %>% unlist()))

merged_accuracy_data <- data %>% 
  filter(classes ==2) %>%
  filter(num_clusters=="20000") %>%
  filter(filter=="filter1") %>%
  filter(!(dataset %in% c("prism-metabolomics-SRP129027", "vibrio-cholerae-PRJNA723557", "sepsis-PRJNA507824"))) %>%
  filter(cluster_approach %in% c("shiftDist-levFilter", "masked-aa-clustered", "masked-nucleotide-clustered", "shiftDist-keepTopES")) %>%
  filter(model %in% c("ohe", "hyena_normalized", "hyena_unnormalized")) %>%
  select(dataset, filter, cluster_approach, metadata, model, accuracy, sensitivity, specificity, num_clusters) %>% 
  distinct(filter, cluster_approach, metadata, model, .keep_all=T) %>% select(-num_clusters) 

# filter down for simpler figure
if (FALSE){
  merged_accuracy_data <- merged_accuracy_data %>% filter(model %in% c("hyena_normalized", "esm_normalized", "ohe")) %>%
    filter(cluster_approach %in% c("shiftDist-levFilter", "masked-aa-clustered"))
}

merged_accuracy_data <- merged_accuracy_data %>% filter(!(dataset %in% c("y1000-genomes-withinOrder")))

merged_accuracy_data <- merged_accuracy_data %>% mutate(dataset=ifelse(str_detect(dataset,"y1000"), "y1000", dataset))
datasets <- unique(merged_accuracy_data$dataset)
datasets <- c(datasets, "ALL")

# Create empty dataframes to store results
cluster_approach_results <- data.frame()
model_results <- data.frame()

# For each dataset
for(ds in datasets) {
  cat("\n\nAnalyzing dataset:", ds, "\n")
  if (ds=="ALL") {
    dataset_data <- merged_accuracy_data
  } else {
    dataset_data <- merged_accuracy_data %>% filter(dataset == ds)
  }
  
  # 1. Compare clustering approaches
  cat("\n--- Clustering Approach Comparisons ---\n")
  # Get all pairs of cluster approaches
  cluster_approaches <- unique(dataset_data$cluster_approach)
  if (length(cluster_approaches) >=2) {
    pairs <- combn(cluster_approaches, 2, simplify = FALSE)
    
    for(pair in pairs) {
      ca1 <- pair[1]
      ca2 <- pair[2]
      
      # Get paired data (same filter, metadata, model)
      data1 <- dataset_data %>% filter(cluster_approach == ca1) %>% 
        select(filter, metadata, model, accuracy)
      data2 <- dataset_data %>% filter(cluster_approach == ca2) %>% 
        select(filter, metadata, model, accuracy)
      
      # Join to get pairs
      paired_data <- inner_join(data1, data2, by = c("filter", "metadata", "model"), 
                                suffix = c("_ca1", "_ca2"))
      
      if(nrow(paired_data) > 0) {
        # Perform Wilcoxon signed-rank test
        test_result <- wilcox.test(paired_data$accuracy_ca1, paired_data$accuracy_ca2, 
                                   paired = TRUE, alternative = "two.sided")
        
        # Calculate median difference and effect size
        median_diff <- median(paired_data$accuracy_ca1 - paired_data$accuracy_ca2, na.rm=T)
        
        # Store results
        result_row <- data.frame(
          dataset = ds,
          comparison_type = "cluster_approach",
          group1 = ca1,
          group2 = ca2,
          p_value = test_result$p.value,
          median_diff = median_diff,
          n_pairs = nrow(paired_data)
        )
        
        cluster_approach_results <- rbind(cluster_approach_results, result_row)
        
        cat(sprintf("%s vs %s: p-value = %.4f, median diff = %.4f (n=%d)\n", 
                    ca1, ca2, test_result$p.value, median_diff, nrow(paired_data)))
      }
    }
  }
  
  # 2. Compare models
  cat("\n--- Model Comparisons ---\n")
  # Get all pairs of models
  models <- unique(dataset_data$model)
  if (length(models) >=2) {
    pairs <- combn(models, 2, simplify = FALSE)
    
    for(pair in pairs) {
      mdl1 <- pair[1]
      mdl2 <- pair[2]
      
      # Get paired data (same filter, metadata, cluster_approach)
      data1 <- dataset_data %>% filter(model == mdl1) %>% 
        select(filter, metadata, cluster_approach, accuracy)
      data2 <- dataset_data %>% filter(model == mdl2) %>% 
        select(filter, metadata, cluster_approach, accuracy)
      
      # Join to get pairs
      paired_data <- inner_join(data1, data2, by = c("filter", "metadata", "cluster_approach"), 
                                suffix = c("_mdl1", "_mdl2"))
      
      if(nrow(paired_data) > 0) {
        # Perform Wilcoxon signed-rank test
        test_result <- wilcox.test(paired_data$accuracy_mdl1, paired_data$accuracy_mdl2, 
                                   paired = TRUE, alternative = "two.sided")
        
        # Calculate median difference
        median_diff <- median(paired_data$accuracy_mdl1 - paired_data$accuracy_mdl2, na.rm=T)
        
        # Store results
        result_row <- data.frame(
          dataset = ds,
          comparison_type = "model",
          group1 = mdl1,
          group2 = mdl2,
          p_value = test_result$p.value,
          median_diff = median_diff,
          n_pairs = nrow(paired_data)
        )
        
        model_results <- rbind(model_results, result_row)
        
        cat(sprintf("%s vs %s: p-value = %.4f, median diff = %.4f (n=%d)\n", 
                    mdl1, mdl2, test_result$p.value, median_diff, nrow(paired_data)))
      }
    }
  }
}

# Apply multiple testing correction
cluster_approach_results$p_adj <- p.adjust(cluster_approach_results$p_value, method = "BH")
model_results$p_adj <- p.adjust(model_results$p_value, method = "BH")

# Create summary tables
cluster_summary <- cluster_approach_results %>%
  arrange(dataset, p_adj) %>%
  select(dataset, group1, group2, p_value, p_adj, median_diff, n_pairs)
write_tsv(cluster_summary, paste0(opt$output_folder, "/cluster_approaches_wilcox_signed_rank_comparisons.tsv"), col_names = T, quote="needed")

model_summary <- model_results %>%
  arrange(dataset, p_adj) %>%
  select(dataset, group1, group2, p_value, p_adj, median_diff, n_pairs)
write_tsv(model_summary, paste0(opt$output_folder, "/models_wilcox_signed_rank_comparisons.tsv"), col_names = T, quote="needed")

# # Display summaries
# cat("\n\n=== CLUSTER APPROACH COMPARISON SUMMARY ===\n")
# print(cluster_summary)
# 
# cat("\n\n=== MODEL COMPARISON SUMMARY ===\n")
# print(model_summary)


create_forest_plot <- function(results_df, title) {
  # Add significance stars
  results_df <- results_df %>%
    mutate(significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  # Create comparison labels and determine which group performed better
  results_df <- results_df %>%
    mutate(
      better_group = ifelse(median_diff > 0, group1, group2),
      comparison = paste0(group1, " vs ", group2),
      # Create a cleaner y-axis label without "vs"
      plot_label = comparison
    )
  
  # Plot
  p <- ggplot(results_df, aes(x = median_diff, y = reorder(plot_label, median_diff), color=significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = significance), size = 3) +
    geom_errorbarh(aes(xmin = median_diff - 1.96*median_diff/sqrt(n_pairs),
                       xmax = median_diff + 1.96*median_diff/sqrt(n_pairs)),
                   height = 0.3) +
    # Add direction labels on both sides
    # annotate("text", x = min(results_df$median_diff, -0.01) * 1.2, 
    #          y = Inf, label = "← Group 2 better", 
    #          hjust = 0, vjust = 1.5, color = "blue4", fontface = "bold") +
    # annotate("text", x = max(results_df$median_diff, 0.01) * 1.2, 
    #          y = Inf, label = "Group 1 better →", 
    #          hjust = 1, vjust = 1.5, color = "red4", fontface = "bold") +
    # Add group names as text on opposite sides
    geom_text(aes(label = group2), 
              x = -0.2,
              hjust = 0, size = 3, color = "blue4") +
    geom_text(aes(label = group1), 
              x=0.2,
              hjust = 1, size = 3, color = "red4") +
    facet_wrap(~ dataset, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("ns" = "gray60", "*" = "#0072B2", 
                                  "**" = "#D55E00", "***" = "#CC0000")) +
    coord_cartesian(xlim=c(-0.2,0.2)) +
    labs(title = title,
         x = "Median accuracy difference",
         y = "",
         color = "Significance") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "gray90"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.subtitle = element_text(face = "italic"),
      panel.spacing = unit(1, "lines"),
      axis.text.y=element_blank()
    )
  
  # Determine appropriate x-axis limits with some padding
  max_abs_diff <- max(abs(results_df$median_diff), na.rm = TRUE)
  padding <- max_abs_diff * 0.6
  p <- p + scale_x_continuous(limits = c(-max_abs_diff - padding, max_abs_diff + padding))
  
  return(p)
}

# Create and save the improved forest plots
forest_cluster <- create_forest_plot(cluster_approach_results, 
                                     "Effect Size of Clustering Approach Comparisons")
forest_model <- create_forest_plot(model_results, 
                                   "Effect Size of Model Comparisons")
# Save the plots
ggsave(paste0(opt$output_folder, "/forest_cluster_approaches_comparison.png"), 
       forest_cluster, width = 8, height = length(datasets)+3)
ggsave(paste0(opt$output_folder, "/forest_models_comparison_slimmed.png"), 
       forest_model, width = 8, height = length(datasets)+3)
