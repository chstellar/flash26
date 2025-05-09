suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggpubr))

option_list <- list(
  make_option(c("-a", "--nonzero_annotations"), type = "character", default = NULL, 
              help = "Path to the nonzero annotations tsv file", metavar = "character"),
  make_option(c("-c", "--output"), type = "character", default = NULL, 
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

# testing
#opt$nonzero_annotations = "results/eFaecium-CollEtAl/filter1/shiftDist-levFilter/hyena/normalized/eFaecium-CollEtAl_hyena_adelie_results_top20000_k54_s54_nonzero_coefficients_blast_annotated.tsv"
#opt$output = "test.pdf"

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


# read in input file
dt <- fread(opt$nonzero_annotations)

categories <- dt %>% select(metadata_category, accuracy) %>% distinct() %>% arrange(-accuracy) %>% pull(metadata_category)

get_max_abs_value <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    max(abs(nums), na.rm = TRUE)
  })
}


pdf(opt$output, width=8, height=8)

# write a title page first
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, paramaters['dataset'])
text(5, 7, paramaters['filter'])
text(5, 6, paramaters['cluster_approach'])
text(5, 5, paramaters['model'])
text(5, 4, paste("At most", paramaters['num_clusters'], "clusters"))

for (category in categories) {
  
summ_dt <- dt %>% filter(metadata_category==category) %>%
  separate_longer_delim(features, delim = "},") %>% 
  mutate(products=str_extract(features, "'product': \\['([\\w\\s-]+)'\\]", group=1)) %>% 
  mutate(genes=str_extract(features, "'gene': \\['([\\w\\s-]+)'\\]", group=1)) %>% 
  select(-features) %>% mutate(max_coefficient=get_max_abs_value(coefficients)) %>% arrange(-max_coefficient) %>%
  select(metadata_category, accuracy, max_coefficient, cluster, feature, query, products, genes) %>%
  mutate(query = str_remove(query, "cluster_\\d+_")) %>%
  group_by(cluster) %>% mutate(query=paste(unique(na.omit(query)), collapse = ",")) %>%
  ungroup() %>%
  distinct(cluster,products,genes,.keep_all = T) %>% group_by(cluster) 

if (opt$products) {
  summ_dt <- summ_dt %>% 
    mutate(label=ifelse(!is_empty(unique(na.omit(products))), paste(unique(na.omit(products)),collapse=","), paste(unique(na.omit(genes)), collapse=","))) %>% 
    distinct(cluster, label, .keep_all=T) %>% select(-products) %>% ungroup()
} else {
  summ_dt <- summ_dt %>% 
    mutate(label=ifelse(!is_empty(unique(na.omit(genes))), paste(unique(na.omit(genes)), collapse=","), paste(unique(na.omit(products)),collapse=","))) %>% 
    distinct(cluster, label, .keep_all=T) %>% select(-products) %>% ungroup()
}

summ_dt <- summ_dt %>% 
  mutate(rank=row_number()) %>%
  mutate(largest_coef=max(max_coefficient)) %>%
  mutate(coef_mag=max_coefficient/largest_coef) %>% 
  mutate(color=NA) %>%
  mutate(color=ifelse(!is.na(label), "blast", color)) %>%
  mutate(color = ifelse(grepl(known_causes, label, ignore.case=T), "known_cause", color)) %>%
  mutate(label = str_wrap(str_trunc(label, width =100, side="right"), width = 30)) %>% 
  mutate(label = replace_na(label, ""))

accuracy <- summ_dt$accuracy %>% unique()

dataset <- str_extract(opt$nonzero_annotations, "results/([A-Za-z\\d-]+)/filter", group=1)
make_title <- paste(category, "in", dataset)

p <- summ_dt %>% head(opt$num_hits) %>%
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
}

dev.off()
