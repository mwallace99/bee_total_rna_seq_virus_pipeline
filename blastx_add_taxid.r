# ============================================================================
# blastx_add_taxid.r
# ----------------------------------------------------------------------------
# Reads DIAMOND BLASTX output (outfmt 6 with staxids), expands multi-taxid
# rows, and joins with NCBI scientific names.
#
# Usage:
#   Rscript blastx_add_taxid.r --file input.tsv --names names.txt --output out.tsv
#
# The taxid->name table (--names) is tab separated: <taxid> <TAB> <name>.
# (Path provided on the command line so the script is not tied to one project;
# default falls back to $NAMES_TXT in the environment if --names is omitted.)
# ============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(argparse)

parser <- ArgumentParser(description = "Process BLAST output and join with scientific names")
parser$add_argument("--file", type = "character", help = "The BLAST file to process")
parser$add_argument("--names", type = "character",
                    default = Sys.getenv("NAMES_TXT", ""),
                    help = "Tab-separated taxid<TAB>name table")
parser$add_argument("--output", type = "character", help = "The output file to save")
args <- parser$parse_args()

if (is.null(args$names) || args$names == "") {
  stop("No names table provided: pass --names <file> or set NAMES_TXT in the environment")
}

# Read the taxid -> scientific name table
names_scientific <- read_delim(args$names,
                               delim = "\t",
                               col_names = c("taxid", "name"),
                               trim_ws = TRUE,
                               col_types = "cc") |>
  mutate(taxid = as.character(taxid))

# Read the BLAST output file
sample_blast <- read_delim(args$file, delim = "\t", col_names = FALSE)

# Assign column names matching the DIAMOND outfmt 6 used in stage 14
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
