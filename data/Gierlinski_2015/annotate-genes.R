# ==============================================================================
# Meta
# ==============================================================================
# annotate-genes
# ------------------------------------------------
# Author: James Hawley
# Description: annotate yeast genes for easy plotting


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

r64_1_1 <- fread(
	file.path("..", "Kallisto_s-cerevisiae_Index", "Saccharomyces_cerevisiae.R64-1-1.96.gtf"),
	skip = 5,
	sep = "\t",
	header = FALSE,
	col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
)

# ==============================================================================
# Analysis
# ==============================================================================

# 1. Extract gene information
# ------------------------------------------------
loginfo("Extracting gene information")

# convert ;-separated values into their own columns for certain attributes
genes <- r64_1_1[feature == "gene"]
genes_split <- genes[, tstrsplit(attribute, "; ")]

# get gene IDs
genes$gene_id <- genes_split[, tstrsplit(V1, "\"")][, V2]

# some genes don't have names. if this is the case, adjust columns
genes_with_names <- genes_split[, grepl("gene_name", V2)]
genes_split[!genes_with_names, V4 := V3]
genes_split[!genes_with_names, V3 := V2]

# for genes without a name, use the ID in place of the name
genes_split[!genes_with_names, V2 := V1]
genes$gene_name <- genes_split[, tstrsplit(V2, "\"")][, V2]

# get gene biotype
genes$gene_biotype <- genes_split[, tstrsplit(V4, "\"")][, V2]

# drop the following columns
genes[, `:=`(
	source = NULL,
	feature = NULL,
	score = NULL,
	frame = NULL,
	attribute = NULL
)]

# 2. Extract transcripts information
# ------------------------------------------------
loginfo("Extracting transcript information")

# convert ;-separated values into their own columns for certain attributes
tx <- r64_1_1[feature == "transcript"]
tx_split <- tx[, tstrsplit(attribute, "; ")]

# get gene and transcript IDs
tx$gene_id <- tx_split[, tstrsplit(V1, "\"")][, V2]
tx$transcript_id <- tx_split[, tstrsplit(V2, "\"")][, V2]

# some transcripts don't have gene names. if this is the case, adjust columns
tx_split[grepl("transcript_biotype", V6), V7 := V6]
tx_split[grepl("transcript_source", V5), V6 := V5]
tx_split[grepl("gene_biotype", V4), V5 := V4]
tx_split[grepl("gene_source", V3), V4 := V3]
# for genes without a name, use the ID in place of the name
tx_split[grepl("gene_source", V3), V3 := V1]

# get relevant information
tx$gene_name <- tx_split[, tstrsplit(V3, "\"")][, V2]
tx$gene_biotype <- tx_split[, tstrsplit(V5, "\"")][, V2]
tx$transcript_biotype <- tx_split[, tstrsplit(V7, "\"")][, V2]

# drop the following columns
tx[, `:=`(
	source = NULL,
	feature = NULL,
	score = NULL,
	frame = NULL,
	attribute = NULL
)]

# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")

fwrite(
    genes,
    "S-cerevisiae.genes.tsv",
    sep = "\t",
    col.names = TRUE
)

fwrite(
    tx,
    "S-cerevisiae.transcripts.tsv",
    sep = "\t",
    col.names = TRUE
)
