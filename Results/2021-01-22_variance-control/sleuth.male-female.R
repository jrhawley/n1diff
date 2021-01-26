# ==============================================================================
# Meta
# ==============================================================================
# sleuth
# ------------------------------------------------
# Author: James Hawley
# Description: Differential expression analysis with all samples


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))


# ==============================================================================
# Functions
# ==============================================================================


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Preparing sleuth objects")
so <- sleuth_prep(
    design,
    extra_bootstrap_summary = TRUE,
    num_cores = 8
)

logging("Fitting model and hypothesis testing")
so <- sleuth_fit(so, ~condition, "full")

# perform differential analysis
so <- sleuth_wt(so, "conditionmale")

# extract results
so_transcripts <- as.data.table(sleuth_results(
    so,
    "conditionmale",
    "wt",
    show_all = FALSE,
    pval_aggregate = FALSE
))
so_genes <- as.data.table(sleuth_results(
    so,
    "conditionmale",
    "wt",
    show_all = FALSE,
    pval_aggregate = TRUE
))


# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")
fwrite(
    so_transcripts,
    "sleuth.transcripts.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    so_genes,
    "sleuth.genes.tsv",
    sep = "\t",
    col.names = TRUE
)

