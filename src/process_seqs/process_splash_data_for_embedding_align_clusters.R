# process_splash_data_for_embedding.R
# Daniel Cotter
# 2024-09-13

# This script takes in a path to SPLASH SATC files as well as a list of anchors
# and a list of anchor clusters. It then dumps the anchors from the SATC files 
# and reformats a sequence for each sample that can be used in downstream
# analyses. The output is a fasta file and a tsv file with the sample id and 
# the sequence.

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(stringdist))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-a", "--anchor_file"), "File with list of anchors", type="character"),
  make_option(c("-c", "--cluster_file"), "File with list of anchor clusters", type="character"),
  make_option(c("i", "--id_mapping"), "File with sample id mapping", type="character"),
  make_option(c("-s", "--satc_files"), "Path to SPLASH SATC file directory", type="character"),
  make_option(c("-o", "--output_prefix"), "Output prefix.", type="character"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character"),
  make_option(c("--satc_util_bin"), "Path to SATC Util binary folder",
              type="character", default="/oak/stanford/groups/horence/dcotter1/satc_utils/"),
  make_option(c("--num_cores"), "Number of cores to use", type="integer", default = 8)
)


# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set up parallel processing
plan(multicore, workers = opt$num_cores)

# check that user specified all files 
if (!file.exists(opt$anchor_file) | !file.exists(opt$cluster_file) | is.null(opt$satc_files) | is.null(opt$output_prefix) | !file.exists(opt$id_mapping)) {
  stop("Must provide anchor file, cluster file, satc files, id mapping file, and output prefix")
}

# create a temporary directory to store intermediate files
if (!is.null(opt$temp_dir)) {
  temp_dir <- ifelse(grepl("/$", opt$temp_dir),
                     opt$temp_dir, 
                     paste0(opt$temp_dir, "/"))
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output_prefix), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

# read in the anchor cluster file
anchor_clusters <- fread(opt$cluster_file, 
                    header = F, col.names = c("cluster_id", "anchor")) %>%
  group_by(cluster_id) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# list all of the .satc files in the result_satc folder
satc_files <- list.files(opt$satc_files, pattern = ".satc", full.names = T)
satc_files <- data.frame(satc_file=satc_files) %>% 
  mutate(satc_dump = gsub(".satc", ".satc.dump", 
                          file.path(opt$temp_dir, 'dumped', basename(satc_file))))
system(paste("mkdir -p", file.path(opt$temp_dir, "dumped")))

# declare a satc file for the output of all the dump files
all_satc_file <- file.path(opt$temp_dir, "all_satc_merged.txt")

if (!file.exists(all_satc_file)) {
  # dump the satc files
  future_walk2(satc_files$satc_file, satc_files$satc_dump, \(x,y) system(
    paste0(file.path(opt$satc_util_bin, "satc_dump"), " --anchor_list ", opt$anchor_file,
           " --sample_names ", opt$id_mapping, " ", x, " ", y)))
  
  # merge them into one file
  walk(satc_files$satc_dump, \(x) system(paste("cat", x, ">>", all_satc_file)))

  # remove any lines that start or end in [ACTG]
  system(paste("grep -v '^[ACTG]' ", all_satc_file, " | grep -v '[ACTG]$' > ", 
               file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt")))
  system(paste("mv", file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt"), all_satc_file))
}

# undump the satc files
all_satc_undumped <- file.path(opt$temp_dir, "all_satc_merged.undumped.txt")
all_satc_temp_mapping <- file.path(opt$temp_dir, "all_satc_merged.temp_mapping.txt")
if (!file.exists(all_satc_undumped) | !file.exists(all_satc_temp_mapping)) {
  system(paste(file.path(opt$satc_util_bin, "satc_undump"), "-i", all_satc_file, 
               "-o", all_satc_undumped, "-m", all_satc_temp_mapping))
}

# filter the satc files 
all_satc_filtered <- file.path(opt$temp_dir, "all_satc.filtered.txt")
if (!file.exists(all_satc_filtered)) {
  system(paste(file.path(opt$satc_util_bin, "satc_filter"), 
               "-i", all_satc_undumped, "-o", all_satc_filtered,
               "-d", opt$anchor_file, "-n", 1))
}

# redump the satc file 
all_satc_filtered_dump <- file.path(opt$temp_dir, "all_satc.filtered.dump")
if (!file.exists(all_satc_filtered_dump)) {
  system(paste(file.path(opt$satc_util_bin, "satc_dump"),
               "--sample_names", opt$id_mapping,
               all_satc_filtered, all_satc_filtered_dump))
}

# read in the dumped satc file
satc_dt <- fread(all_satc_filtered_dump, header=F,
                 col.names=c("sample", "anchor", "target", "count"))

# grab the top anchor per cluster as a representative anchor
representative_anchors <- anchor_clusters %>% group_by(cluster_id) %>% 
  distinct(cluster_id, .keep_all=T) %>% ungroup() %>% pull(anchor)

# read in the satc and pivot it wider 
wide_satc <- merge(satc_dt, anchor_clusters, by="anchor", all.x=TRUE)
wide_satc <- as.data.table(wide_satc)
wide_satc <- wide_satc[order(cluster_id, rank)]
wide_satc <- unique(wide_satc, by=c("sample", "cluster_id"))

wide_satc <- wide_satc %>% mutate(seq=str_c(anchor, target, sep="")) %>% 
  select(sample, cluster_id, seq)

wide_satc <- dcast(wide_satc, sample ~ cluster_id, value.var="seq")
wide_satc <- as.data.frame(wide_satc)

# write a function that operates on every cluster column of the wide satc
# and aligns the nonrepresentative anchors to the representative anchor
align_to_representative <- function(x, colname, representative_anchor) {
  # x is a one column data frame
  x <- data.frame(sequence=x)
  # grab the unique sequences in the column
  unique_seqs <- x %>% group_by(sequence) %>% summarise(n=n()) %>% ungroup()
  
  # if there is only one anchor return the original column
  n_anchors <- unique_seqs %>% mutate(sequence=substr(sequence, 1, 27)) %>% 
    filter(!is.na(sequence)) %>% 
    pull(sequence) %>% unique() %>% length()
  if (n_anchors == 1) {
    colnames(x) <- c(colname)
    return(x)
  }
  
  # grab the sequence to align to (the most abundant one containing the anchor)
  representative_seq <- unique_seqs %>% filter(str_detect(sequence, representative_anchor)) %>% 
    filter(n == max(n)) %>% head(1) %>% pull(sequence)
  
  # format the unique sequences as a character vector
  unique_seqs <- unique_seqs %>% pull(sequence)
  
  # function to perform alignment on a unique set of sequences
  perform_alignment <- function(seq, ref_seq) {
    # Calculate the Levenshtein edit moves
    lev_moves <- attr(adist(seq, ref_seq, counts = TRUE), "trafos")[[1]]
    
    # Initialize empty aligned sequence
    aligned_seq <- ""
    
    # Process each move backwards
    # split lev moves into a vector one letter each
    lev_moves <- strsplit(lev_moves, "")[[1]]
    # perform the moves on the seq to get the aligned seq (skip S)
    # if there is an M keep the corresponding letter from the seq
    # if there is a D, delete the corresponding letter from the seq
    # if there is a I, insert an N in the aligned seq
    # if there is an S, do nothing
    for (move in lev_moves) {
      if (move=="M") {
        aligned_seq <- paste0(aligned_seq, substr(seq, 1, 1))
        seq <- substr(seq, 2, nchar(seq))
      } else if (move=="D") {
        seq <- substr(seq, 2, nchar(seq))
      } else if (move=="I") {
        aligned_seq <- paste0(aligned_seq, "N")
      }
    }
    
    # first trim the leading Ns from the aligned seq
    aligned_seq <- substr(aligned_seq, 4, nchar(aligned_seq))
    
    # trim or pad with Ns the aligned sequence to be the same length as the reference
    if (nchar(aligned_seq) < nchar(ref_seq)) {
      aligned_seq <- str_pad(aligned_seq, nchar(ref_seq), pad="N", side="right")
    } else if (nchar(aligned_seq) > nchar(ref_seq)) {
      aligned_seq <- substr(aligned_seq, 1, nchar(ref_seq))
    }
    
    # Return the final aligned sequence
    return(aligned_seq)
  }
  
  # Function to align all sequences against the reference using the above alignment logic
  align_sequences <- function(seqs, ref_seq) {
    aligned_seqs <- sapply(seqs, function(seq) perform_alignment(seq, ref_seq))
    return(aligned_seqs)
  }
  
  # align the unique sequences to the representative anchor
  aligned_seqs <- align_sequences(unique_seqs, representative_seq)
  
  # mutate the column in x mapping each unique sequence to its aligned sequence
  x <- x %>% mutate(sequence=map_chr(sequence, \(x) ifelse(is.na(x), NA, aligned_seqs[x])))
  colnames(x) <- c(colname)
  
  return(x)
}

# apply the alignment function to each cluster column
wide_satc <- cbind(wide_satc[1], 
                   future_pmap_dfc(list(wide_satc[,2:3],
                                        colnames(wide_satc[,2:3]),
                                        representative_anchors[1:2]),
                                   \(x,y,z) align_to_representative(x, y, z)))

# add the representative anchors with Ns to the wide satc where there are NAs
wide_satc <- cbind(wide_satc[1],
                   map2_df(wide_satc[,2:ncol(wide_satc)], 
                           1:length(representative_anchors), 
                           \(x,y) ifelse(is.na(x), str_c(representative_anchors[y], strrep("N", 27), sep = ""), x)))
wide_satc <- wide_satc %>% ungroup()

# join the columns together into one sequence and write to a tsv
seqs <- wide_satc %>% unite(seq, -sample, sep="")
seqs %>% write_tsv(paste0(opt$output_prefix, "_sample_sequences.tsv"))

# write the data to a fasta file
dna <- Biostrings::DNAStringSet(seqs$seq)
names(dna) <- seqs$sample
writeXStringSet(dna, paste0(opt$output_prefix, "_sample_sequences.fasta"))

