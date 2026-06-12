# ============================================================================
# blastx_add_taxid_Version2.r
# ----------------------------------------------------------------------------
# Reads Diamond BLASTX output (outfmt 6 with staxids), expands multi-taxid
# rows, and joins with NCBI scientific names.
#
# Usage: Rscript blastx_add_taxid_Version2.r --file input.txt --output output.txt
# ============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description = "Process BLAST output and join with scientific names")
parser$add_argument("--file", type = "character", help = "The BLAST file to process")
parser$add_argument("--output", type = "character", help = "The output file to save")
args <- parser$parse_args()

# Read the names.txt (located in same directory as this script)
names_scientific <- read_delim("/home/576/mw9045/Jobs/UKBombus_BLAST/names.txt",
                               delim = "\t",
                               col_names = c("taxid", "name"),
                               trim_ws = TRUE,
                               col_types = "cc") |>
  mutate(taxid = as.character(taxid))

# Read the BLAST output file
sample_blast <- read_delim(args$file, delim = "\t", col_names = FALSE)

# Assign column names matching Diamond outfmt 6
colnames(sample_blast) <- c("query", "hit_id", "identity", "alignment_length", "evalue", "score", "taxid")

# Expand rows to handle multiple taxIDs (semicolon-separated)
expanded_blast <- sample_blast |>
  separate_rows(taxid, sep = ";\\s*") |>
  mutate(taxid = as.character(taxid))

# Join with scientific names
blast_with_names <- expanded_blast |>
  left_join(names_scientific, by = "taxid")

# Save output
write_delim(blast_with_names, args$output, delim = "\t")
