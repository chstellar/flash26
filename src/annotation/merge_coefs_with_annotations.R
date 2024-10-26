# Load necessary libraries
library(optparse)
library(data.table)

# Define command line options
option_list <- list(
    make_option(c("-a", "--annotations"), type = "character", default = NULL, 
                            help = "Path to the annotations CSV file", metavar = "character"),
    make_option(c("-c", "--coefficients"), type = "character", default = NULL, 
                            help = "Path to the coefficients CSV file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL, 
                            help = "Path to the output CSV file", metavar = "character")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$annotations) || is.null(opt$coefficients) || is.null(opt$output)) {
    print_help(opt_parser)
    stop("All arguments must be supplied", call. = FALSE)
}

# Read in the data
annotations <- fread(opt$annotations, header = TRUE)
coefficients <- fread(opt$coefficients, header = TRUE)
coefficients <- coefficients %>% mutate(cluster = str_extract(feature, "(cluster_\\d+)_", group = 1))

# Merge the data on the 'cluster' column 
# there may be multiple entries in annotations that match
# in this case, the coefficients will be repeated for each match
merged_data <- merge(annotations, coefficients, by = "cluster")

# Write the merged data to a new CSV file
fwrite(merged_data, opt$output, row.names = FALSE)

# Print a message indicating the script has finished
cat("Merging complete. Output saved to", opt$output, "\n")