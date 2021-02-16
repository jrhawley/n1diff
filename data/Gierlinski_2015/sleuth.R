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
meta <- fread("experimental-design.tsv")

# set the "wt" as "control" so that it is the baseline for comparison
meta[condition == "wt", condition := "ctrl"]


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Preparing sleuth objects")
so <- sleuth_prep(
    meta,
    extra_bootstrap_summary = TRUE,
    num_cores = 8
)

loginfo("Fitting model and hypothesis testing")
so <- sleuth_fit(so, ~condition, "full")

# perform differential analysis
so <- sleuth_wt(so, "conditionmu")

# save data
saveRDS(so, file.path("Sleuth", "sleuth-object.rds"))

# extract results
so_genes <- as.data.table(sleuth_results(
    so,
    "conditionmu",
    "wt",
    show_all = FALSE
))

# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")
fwrite(
    so_genes,
    file.path("Sleuth", "genes.tsv"),
    sep = "\t",
    col.names = TRUE
)

