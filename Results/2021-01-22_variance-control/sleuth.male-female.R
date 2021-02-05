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
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# load transcript ID to gene ID mapping
t2g_map <- fread(
    "../../Data/Kallisto_Index/transcripts_to_genes.txt",
    sep = "\t",
    header = FALSE,
    col.names = c("target_id", "ens_gene", "ext_gene")
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Preparing sleuth objects")
so <- sleuth_prep(
    meta,
    extra_bootstrap_summary = TRUE,
    num_cores = 8,
    target_mapping = t2g_map,
    aggregation_column = "ens_gene",
    # filter out transcripts that have < 10 reads in > 60% of the samples
    filter_fun = function(x) length(which(x > 10)) / length(x) > 0.6
)

loginfo("Fitting model and hypothesis testing")
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
saveRDS(so, "sleuth.object.rds")

