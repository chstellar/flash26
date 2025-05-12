suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(msa))


option_list <- list(
  make_option(c("--nonzero_annotations"), type = "character", default = NULL, 
              help = "Path to the nonzero annotations tsv file", metavar = "character"),
  make_option(c("--clusters"), type = "character", default = NULL, 
              help = "Path to the clusters tsv file", metavar = "character"),
  make_option(c("--feather_file"), type = "character", default = NULL, 
              help = "Path to the X matrix feather file", metavar = "character"),
  make_option(c("--sample_seqs"), type = "character", default = NULL, 
              help = "Path to the sample sequences file", metavar = "character"),
  make_option(c("--metadata"), type = "character", default = NULL, 
              help = "Path to the metadata tsv file", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL, 
              help = "Path to set of output plots", metavar = "character"),
  make_option(c("--products"), type= "logical", default=FALSE, action="store_true",
              help = "default to using products for column names instead of genes"),
  make_option(c("--num_hits"), type="numeric", default=10,
              help = "num nonzero coefficients to plot", metavar = "numeric")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$nonzero_annotations) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All arguments must be supplied", call. = FALSE)
}

# set known_causes to be empty (can be changed for interactive experimentation on specific datasets)
known_causes = "NNNNNNNNNNNNNNN"

# # testing
# setwd("/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline")
# opt$nonzero_annotations = "results/eFaecium-CollEtAl/filter1/shiftDist-levFilter/hyena/normalized/eFaecium-CollEtAl_hyena_adelie_results_top20000_k54_s54_nonzero_coefficients_blast_annotated.tsv"
# opt$clusters = "results/eFaecium-CollEtAl/filter1/shiftDist-levFilter/eFaecium-CollEtAl_sequences_per_cluster_top20000-clusters_k54_s54.tsv"
# opt$feather = "/scratch/users/dcotter1/metaSPLASH_workflows_v2/eFaecium-CollEtAl/eFaecium-CollEtAl_hyena_top_variance_features_for_glmnet_filter1_shiftDist-levFilter_top20000_k54_s54_normalized.feather"
# opt$sample_seqs = "/scratch/users/dcotter1/metaSPLASH_workflows_v2/eFaecium-CollEtAl/eFaecium-CollEtAl_prepared_sequences_filter1_shiftDist-levFilter_top20000_sample_sequences.tsv"
# opt$metadata = "/oak/stanford/groups/horence/dcotter1/utility_files/metadata/metaSPLASH_metadata/E_faecium_cleaned_resistance_metadata.tsv"
# opt$output = "/oak/stanford/groups/horence/dcotter1/share/250506/test_eFac_more_blast_hits_out.pdf"

filename = data.frame(path=opt$nonzero_annotations)

filename <- filename %>% 
  mutate(num_clusters = str_extract(path, "top(\\d+)", group=1)) %>%
  mutate(path=dirname(gsub("^results/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  dplyr::rename(paramater_set=path) %>%
  mutate(paramater_set=str_replace(paramater_set, "_", "/")) %>% 
  separate(paramater_set, into=c("dataset", "paramater_set"), sep="/") %>%
  mutate(model=str_extract(paramater_set,
                           'hyenaHG38_normalized|hyenaHG38_unnormalized|hyenaMarlowe_normalized|hyenaMarlowe_unnormalized|esm_normalized|esm_unnormalized|hyena_normalized|hyena_unnormalized|ohe')) %>%
  mutate(filter = str_extract(paramater_set, "(filter\\d)_",group=1)) %>% 
  mutate(cluster_approach = str_extract(paramater_set, "filter\\d_([A-Za-z-2]+)_", group=1)) 

paramaters <- filename %>% pivot_longer(everything(), names_to="paramater", values_to="value") %>% deframe()


# Define Function 
get_max_abs_value <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    max(abs(nums), na.rm = TRUE)
  })
}


get_first_coef <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    nums[1]
  })
}

get_first_class <- function(x) {
  sapply(x, function(str) {
    classes <- strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]]
    classes[1]
  })
}

# function to read the nth cluster out of the sample sequences file 
read_nth_cluster <- function(file_path, n) {
  # Calculate start and end positions for the nth cluster
  cluster_length <- 54
  start <- n * cluster_length + 1  # n is now 0-indexed
  end <- (n + 1) * cluster_length
  
  # Construct the awk command
  awk_command <- sprintf("awk 'NR>1 {print $1 \"\t\" substr($2, %d, %d)}' %s", start, cluster_length, file_path)
  
  # Execute the command and read the result
  result <- read.table(text = system(awk_command, intern = TRUE), 
                       col.names = c("sample_name", "sequence"), 
                       stringsAsFactors = FALSE)
  
  return(result)
}

# Function to calculate distance to the most abundant sequence after MSA of unique sequences
calculate_distance_and_align <- function(sequences) {
  # Remove sequences with long strings of N's
  valid_sequences <- sequences[!grepl("N{10,}", sequences)]
  
  if (length(valid_sequences) == 0) {
    return(list(distances = rep(-1, length(sequences)),
                aligned_seqs = rep("", length(sequences))))
  }
  
  # Get unique sequences and their counts
  seq_table <- table(valid_sequences)
  unique_seqs <- names(seq_table)
  
  # Perform MSA on unique sequences
  msa_result <- msa(DNAStringSet(unique_seqs), method = "ClustalOmega")
  
  # Convert MSA result to character vectors
  aligned_seqs <- as.character(msa_result)
  
  # Find the most abundant sequence
  most_abundant <- aligned_seqs[which.max(seq_table)]
  
  # Calculate distances for unique sequences
  unique_distances <- stringdist(aligned_seqs, most_abundant, method = "lv")
  
  # Map distances and aligned sequences back to all valid sequences
  all_valid_distances <- unique_distances[match(valid_sequences, unique_seqs)]
  all_valid_aligned_seqs <- aligned_seqs[match(valid_sequences, unique_seqs)]
  
  # Map distances and aligned sequences back to original sequences, including those with N's
  all_distances <- rep(-1, length(sequences))
  all_aligned_seqs <- rep("", length(sequences))
  valid_indices <- !grepl("N{10,}", sequences)
  all_distances[valid_indices] <- all_valid_distances
  all_aligned_seqs[valid_indices] <- all_valid_aligned_seqs
  
  return(list(distances = all_distances, aligned_seqs = all_aligned_seqs))
}

# read in input files
dt <- fread(opt$nonzero_annotations)
all_clusters <- fread(opt$clusters) %>% select(-kmer)
feather_dt <- feather::read_feather(opt$feather)
all_metadata <- fread(opt$metadata)

categories <- dt %>% select(metadata_category, accuracy) %>% distinct() %>% arrange(-accuracy) %>% pull(metadata_category)

#category = "ampicillin_RIS"

pdf(opt$output, width=12, height=8)

# write a title page first
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, paramaters['dataset'])
text(5, 7, paramaters['filter'])
text(5, 6, paramaters['cluster_approach'])
text(5, 5, paramaters['model'])
text(5, 4, paste("At most", paramaters['num_clusters'], "clusters"))

for (category in categories) {
  tryCatch({
    summ_dt <- dt %>% filter(metadata_category==category) %>%
      separate_longer_delim(features, delim = "},") %>% 
      mutate(products=str_extract(features, "'product': \\['([\\w\\s-]+)'\\]", group=1)) %>% 
      mutate(genes=str_extract(features, "'gene': \\['([\\w\\s-]+)'\\]", group=1)) %>% 
      select(-features) %>% mutate(first_coef=get_first_coef(coefficients)) %>% mutate(max_coefficient=abs(first_coef)) %>% 
      arrange(-max_coefficient) %>% mutate(first_class=get_first_class(classes)) %>%
      rowwise() %>%
      mutate(classes=list(str_split_1(gsub("\\[|\\]", "", classes),pattern=","))) %>%
      ungroup() %>%
      select(metadata_category, accuracy, classes, first_class, first_coef, max_coefficient, cluster, feature, query, identity, products, genes) %>%
      mutate(query = str_remove(query, "cluster_\\d+_")) %>%
      group_by(cluster) %>%
      ungroup() %>%
      distinct(cluster,products,query,genes,.keep_all = T) %>% group_by(cluster) 
    
    if (opt$products) {
      summ_dt <- summ_dt %>% group_by(cluster,query) %>%
        mutate(label=ifelse(!is_empty(unique(na.omit(products))), paste(unique(na.omit(products)),collapse=","), paste(unique(na.omit(genes)), collapse=","))) %>% 
        distinct(cluster, query, label, .keep_all=T) %>% ungroup()
    } else {
      summ_dt <- summ_dt %>% group_by(cluster,query) %>%
        mutate(label=ifelse(!is_empty(unique(na.omit(genes))), paste(unique(na.omit(genes)), collapse=","), paste(unique(na.omit(products)),collapse=","))) %>% 
        distinct(cluster, query, label, .keep_all=T) %>% ungroup()
    }
    
    plot_dt <- summ_dt %>% 
      mutate(largest_coef=max(max_coefficient)) %>%
      mutate(coef_mag=max_coefficient/largest_coef) %>% 
      group_by(cluster, coef_mag) %>% 
      summarise(label=paste(unique(label), collapse=",")) %>%
      ungroup() %>%
      arrange(-coef_mag) %>%
      mutate(rank=row_number()) %>%
      mutate(color=NA) %>%
      mutate(color=ifelse(nchar(label)>1, "blast", color)) %>%
      mutate(color = ifelse(grepl(known_causes, label, ignore.case=T), "known_cause", color)) %>%
      mutate(label = str_wrap(str_trunc(label, width =100, side="right"), width = 30)) %>% 
      mutate(label = replace_na(label, ""))
    
    accuracy <- summ_dt$accuracy %>% unique()
    
    dataset <- str_extract(opt$nonzero_annotations, "results/([A-Za-z\\d-]+)/filter", group=1)
    make_title <- paste(category, "in", dataset)
    
    p <- plot_dt %>% head(opt$num_hits) %>%
      ggplot(aes(x=rank, y=coef_mag, fill=color, label=label)) +
      geom_col() + 
      geom_text(aes(y=coef_mag + 0.05,hjust=0),angle=45) +
      scale_y_continuous("Magnitude relative to\nlargest nonzero coefficient", 
                         limits = c(0,1.5), 
                         labels=scales::label_percent(), breaks=seq(0,1,0.25),
                         expand=c(0,0)) +
      xlab("Rank of nonzero coefficient (by magnitude)") +
      scale_fill_manual(breaks=c("known_cause","blast",NA), values=c("forestgreen", "pink", "grey"),
                        labels=c("Known Cause", "Blast hit", NA)) +
      scale_x_continuous(breaks=seq(1,10,1)) +
      ggtitle(make_title,
              subtitle = paste("Accuracy:", scales::label_percent(accuracy = 0.01)(accuracy))) +
      theme_pubr() + 
      theme(legend.position="none")
    
    print(p)
    
    # select the top N coefficients
    interesting_clusters <- summ_dt %>% distinct(cluster, feature, max_coefficient) %>%
      arrange(-max_coefficient) %>% head(opt$num_hits)
    
    # filter metadata for only current category
    my_metadata <- all_metadata %>% select(sample_name, !!category) %>% rename(metadata:=!!category)
    
    for (j in 1:nrow(interesting_clusters)) {
      
      my_cluster = interesting_clusters[j,]$cluster
      my_feature = interesting_clusters[j,]$feature
      
      
      dt_sub <- summ_dt %>% filter(cluster==my_cluster)
      
      first_beta = unique(dt_sub$first_coef)
      first_class = unique(dt_sub$first_class)
      all_classes = dt_sub[1,]$classes %>% unlist()
      
      seq_sub <- read_nth_cluster(opt$sample_seqs, as.numeric(str_extract(my_cluster, "\\d+")))
      emb_sub <- feather_dt %>% select(sample_name, !!my_feature)
      
      seq_sub <- seq_sub %>% 
        left_join(emb_sub, by="sample_name") %>% 
        dplyr::rename(embedding:=!!my_feature) %>% 
        mutate(embedding=embedding*first_beta)
      
      distances <- calculate_distance_and_align(seq_sub$sequence)
      
      seq_sub$aligned_sequence <- distances$aligned_seq
      seq_sub$lev_dist <- distances$distances
      
      summ_sub_dt <- seq_sub %>% left_join(my_metadata,by="sample_name") %>% 
        group_by(sequence,embedding,lev_dist,aligned_sequence,metadata) %>% 
        summarise(metadata_count=n()) %>%
        filter(metadata %in% all_classes) %>%
        pivot_wider(id_cols=everything(), names_from=metadata, values_from=metadata_count) %>% 
        relocate(aligned_sequence, .after="sequence")
      
      summ_sub_dt <- summ_sub_dt %>% left_join(dt_sub %>% select(query,accuracy,identity,label) %>% 
                                                 rename(sequence=query), by="sequence")
      
      p_sub <- summ_sub_dt %>%
        mutate(across(all_of(all_classes), \(x) replace_na(x, 0))) %>%  # Replace NA with 0 for all specified columns
        mutate(prop_first_class = !!sym(first_class) / rowSums(across(all_of(all_classes)))) %>%  # Calculate proportion
        mutate(total_samples = rowSums(across(all_of(all_classes)))) %>% 
        mutate(label = ifelse(is.na(label), "NO BLAST HIT", label)) %>%
        mutate(label = ifelse(nchar(label) == 0, "MISSING", label)) %>%
        ungroup() %>%
        mutate(label = str_wrap(label, width = 40)) %>%
        mutate(label=gsub(",Pbp5","",label)) %>%
        rename(`Blast Label` = label) %>%
        mutate(label_identity = ifelse(identity == 100 | is.na(identity), "", paste0(round(identity,2), "%"))) %>% 
        mutate(label_identity = replace_na(label_identity, ""))
      
      if (length(unique(p_sub$`Blast Label`)) <= 6) {
        p2 <- p_sub %>% 
          ggplot(aes(x=embedding, y=lev_dist, color=prop_first_class, 
                     shape=`Blast Label`, size=total_samples,
                     label=label_identity)) +
          geom_vline(xintercept = 0, lty="dashed") +
          geom_point(stroke=1.4) + 
          scale_y_continuous(breaks=scales::breaks_width(1)) +
          scale_size_continuous(
            trans = "log", 
            name = "Total Samples",
            breaks = c(1, 10, 100, 1000, 10000),  # Specify breaks for the legend
            limits = c(1, 10000), # Set limits for the size scale
            labels = scales::label_log()
          ) +
          ggrepel::geom_text_repel(aes(label = label_identity),
                                   size = 3,        # Adjust the size of the text
                                   hjust = 0,       # Horizontal justification (0 = left, 0.5 = center, 1 = right)
                                   vjust = 0,        # Vertical justification (0 = bottom, 0.5 = center, 1 = top)
                                   color="black"
          ) +       
          scale_color_gradient(paste0("Proportion\n", first_class), 
                               low = "blue", high = "red", limits = c(0, 1)) +
          theme_minimal() + xlab(expression("Embedding" ~ "\u00D7" ~ beta)) +
          ylab("Levenshtein Distance\n(to most abundant anchor-target)") +
          ggtitle(my_cluster)
      } else {
        p2 <- p_sub %>% 
          mutate(`Blast Label`=ifelse(nchar(label_identity) >1, 
                                      paste0(`Blast Label`, " (", label_identity, ")"),
                                      `Blast Label`)) %>%
          ggplot(aes(x=embedding, y=lev_dist, color=prop_first_class, 
                     label=`Blast Label`, size=total_samples)) +
          geom_vline(xintercept = 0, lty="dashed") +
          geom_point(stroke=1.4) +
          scale_y_continuous(breaks=scales::breaks_width(1)) +
          scale_size_continuous(
            trans = "log", 
            name = "Total Samples",
            breaks = c(1, 10, 100, 1000, 10000),  # Specify breaks for the legend
            limits = c(1, 10000),  # Set limits for the size scale
            labels = scales::label_log()
          ) +
          ggrepel::geom_text_repel(size=4) +
          scale_color_gradient(paste0("Proportion\n", first_class), 
                               low = "blue", high = "red", limits = c(0, 1)) +
          theme_minimal() + xlab(expression("Embedding" ~ "\u00D7" ~ beta)) +
          ylab("Levenshtein Distance\n(to most abundant anchor-target)") +
          ggtitle(my_cluster)
      }
      
      print(p2)
    }
    
  }, error = function(e) {
    message(paste("Error processing category:", category, "\n", e$message))
    # Optionally log the error or take other actions
  })
  
}

dev.off()
# ggsave(filename="/oak/stanford/groups/horence/dcotter1/share/250501/test_plot_eFac_ampicillin_blast_hits.pdf")
