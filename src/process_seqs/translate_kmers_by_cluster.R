# translate_kmers_by_cluster.R
# Daniel Cotter
# 2025-05-12

# This script takes in the kmers per cluster file and translates each cluster
# ensuring the kmers that are translated are in the same frame in each cluster


## import packages --------
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(msa))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(furrr))



## parse arguments --------
# define command line arguments
# define command line arguments
option_list <- list(
  make_option(c("-k", "--kmers"), help="Kmers fasta file", type="character"),
  make_option(c("-t", "--translation_table"), help="translation table to use",
              default=1, type="integer"),
  make_option(c("--num_cores"), "Number of cores to use", type="integer", default = 8),
  make_option(c("-p", "--output"), help="Output fasta file.", type="character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# testing
# opt$kmers <- "/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline/results/eFaecium-CollEtAl/filter1/masked-aa-clustered/eFaecium-CollEtAl_sequences_per_cluster_top20000-clusters_k54_s54.tsv"
# opt$translation_table <- 11
# opt$output <- "/scratch/users/dcotter1/test_translated_output.fasta"

# check that user specified all files
if (is.null(opt$output) | !file.exists(opt$kmers)) {
  stop("Must provide kmers and output file")
}

# set up parallel processing
plan(multisession, workers = opt$num_cores)

## print a summary of the arguments
cat("\n####################\n")
cat("Running transalte_kmers_by_cluster.R with the following arguments:\n")
cat("Kmers file: ", opt$kmers, "\n")
cat("Translation table: ", opt$translation_table, "\n")
cat("Output file: ", opt$output, "\n")
cat("####################\n\n")

## functions 

# Function to determine the shift needed for each sequence in a cluster
align_sequences <- function(cluster_sequences) {
  # Extract the first 25 nucleotides of each sequence
  first_25 <- sapply(cluster_sequences, function(seq) substr(seq, 1, 25))
  
  # Determine the consensus sequence (most common starting sequence)
  consensus <- names(sort(table(first_25), decreasing = TRUE))[1]
  
  # Calculate the shift needed for each sequence to match the consensus
  shifts <- sapply(first_25, function(seq) {
    # Find the best alignment by shifting left or right
    max_shift <- 2
    best_shift <- 0
    best_score <- -Inf
    
    for (shift in -max_shift:max_shift) {
      # Shift the sequence
      shifted_seq <- if (shift < 0) {
        substr(seq, abs(shift) + 1, nchar(seq))
      } else {
        paste0(rep("-", shift), seq)
      }

      # Calculate alignment score (number of matching characters)
      score <- sum(sapply(1:nchar(consensus), function(i) {
        substr(shifted_seq, i, i) == substr(consensus, i, i)
      }))
      
      # Update best shift if score is improved
      if (score > best_score) {
        best_score <- score
        best_shift <- shift
      }
    }
    
    return(best_shift)
  })
  
  return(shifts)
}

# Function to translate a sequence in a specific frame with padding and genetic code
translate_in_frame <- function(seq, frame, genetic_code) {
  # Adjust sequence based on frame and pad the end
  adjusted_seq <- substr(seq, frame, nchar(seq))
  if (frame == 2) {
    adjusted_seq <- substr(adjusted_seq, 1, nchar(adjusted_seq)-2)
  } else if (frame == 3) {
    adjusted_seq <- substr(adjusted_seq, 1, nchar(adjusted_seq)-1)
  }
  
  # Convert to DNAString and translate using the specified genetic code
  dna_seq <- DNAString(adjusted_seq)
  translated_seq <- translate(dna_seq, genetic.code = genetic_code)
  
  return(as.character(translated_seq))
}

# Function to find the best translation frame for a set of sequences
find_best_translation_frame <- function(cluster_sequences, shifts, genetic_code) {
  cluster_sequences <- str_remove(cluster_sequences, "N+")
  # Evaluate all frames for each sequence
  frame_scores <- sapply(1:length(cluster_sequences), function(i) {
    seq <- cluster_sequences[i]
    shift <- shifts[i]
    
    # Adjust sequence based on shift
    adjusted_seq <- if (shift < 0) {
      substr(seq, abs(shift) + 1, nchar(seq) - (3-(abs(shift) %% 3)))
    } else if (shift > 0) {
      substr(seq, 3-(abs(shift) %% 3) + 1, nchar(seq) - abs(shift))
    } else {
      seq
    }
  
    # Translate in all three frames
    translations <- sapply(1:3, function(j) translate_in_frame(adjusted_seq, j, genetic_code))
    
    # Score translations (e.g., by length of sequence without stop codons)
    scores <- sapply(translations, function(trans) {
      # Example scoring: count amino acids before the first stop codon
      stop_pos <- regexpr("\\*", trans)
      if (stop_pos == -1) {
        return(nchar(trans))
      } else {
        return(stop_pos - 1)
      }
    })
    
    return(scores)
  })
  
  # Sum scores across sequences to find the best overall frame
  total_scores <- rowSums(frame_scores)
  best_frame <- which.max(total_scores)
  
  # Translate all sequences using the best frame
  best_translations <- sapply(1:length(cluster_sequences), function(i) {
    seq <- cluster_sequences[i]
    shift <- shifts[i]
    
    # Adjust sequence based on shift
    adjusted_seq <- if (shift < 0) {
      substr(seq, abs(shift) + 1, nchar(seq) - (3-(abs(shift) %% 3)))
    } else if (shift > 0) {
      substr(seq, 3-(abs(shift) %% 3) + 1, nchar(seq) - abs(shift))
    } else {
      seq
    }
    
    # Translate in the best frame
    best_translation <- translate_in_frame(adjusted_seq, best_frame, genetic_code)
    return(best_translation)
  })
  
  best_translations %>% str_replace("\\*.+", "") %>% str_replace("\\*$", "")
  return(best_translations %>% str_replace("\\*.+", "") %>% str_replace("\\*$", ""))
}

translate_in_df <- function(sub_df, genetic_code) {
  my_seqs <- sub_df %>% pull(seq)
  seqs_shifted <- align_sequences(my_seqs)
  translated <- find_best_translation_frame(names(seqs_shifted), 
                                            seqs_shifted, genetic_code = genetic_code)
  sub_df$translation <- translated
  return(sub_df)
}


## read in the clusters file 
cat("Loading data...\n")
dt <- fread(opt$kmers)

dt_groups <- dt %>% group_split(cluster)

my_genetic_code = getGeneticCode("11")

cat("translating all clusters...\n")
all_translations <- future_map(dt_groups, \(x) translate_in_df(x, genetic_code = my_genetic_code)) %>% bind_rows()

cat("Finished translating...\n")
out_AA <- AAStringSet(all_translations %>% select(kmer,translation) %>% deframe())

cat("Writing to file...")
writeXStringSet(out_AA, opt$output)
write_tsv(all_translations, gsub(".tsv", "_translated.tsv", opt$kmers), quote="none", col_names = T)
