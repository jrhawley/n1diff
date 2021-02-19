# ==============================================================================
# Meta
# ==============================================================================
# compare-swi-snf
# ------------------------------------------------
# Author: James Hawley
# Description: Compare differential analysis results of SWI/SNF complex subunits


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))

RESULT_DIR <- file.path("..", "..", "results", "2021-02-10_fewer-transcripts")
N_ITERS <- 30


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# load SWI/SNF complex subunit annotations
swi_snf <- fread("swi-snf-complex-subunits.tsv")

# load fully balanced simulations
bal <- rbindlist(lapply(
    1:N_ITERS,
    function(i) {
        dt <- fread(
            file.path(RESULT_DIR, "Iterations", "SWI-SNF", "balanced", paste0(i, ".genes.tsv")),
            sep = "\t",
            header = TRUE
        )
        # only keep SWI/SNF subunits
        dt <- dt[target_id %in% swi_snf$transcript_id]
        dt[, Iteration := i]
        dt[, Method := "Balanced"]
        return(dt)
    }
))
# load unbalanced OLS simulations
unbal_ols <- rbindlist(lapply(
    1:N_ITERS,
    function(i) {
        dt <- fread(
            file.path(RESULT_DIR, "Iterations", "SWI-SNF", "unbalanced", paste0(i, ".ols.genes.tsv")),
            sep = "\t",
            header = TRUE
        )
        # only keep SWI/SNF subunits
        dt <- dt[target_id %in% swi_snf$transcript_id]
        dt[, Iteration := i]
        dt[, Method := "Unbalanced OLS"]
        return(dt)
    }
))
# load unbalanced JS simulations
unbal_js <- rbindlist(lapply(
    1:N_ITERS,
    function(i) {
        dt <- fread(
            file.path(RESULT_DIR, "Iterations", "SWI-SNF", "unbalanced", paste0(i, ".jse.genes.tsv")),
            sep = "\t",
            header = TRUE
        )
        # only keep SWI/SNF subunits
        dt <- dt[target_id %in% swi_snf$transcript_id]
        dt[, Iteration := i]
        dt[, Method := "Unbalanced JS"]
        return(dt)
    }
))

sims <- rbindlist(list(
    bal[, .SD, .SDcols = c("target_id", "pval", "qval", "b", "Iteration", "Method")],
    unbal_ols[, .SD, .SDcols = c("target_id", "pval", "qval", "b", "Iteration", "Method")],
    unbal_js[, .SD, .SDcols = c("target_id", "pval", "qval", "b", "Iteration", "Method")]
))

# load full set of differential analysis data
full <- fread(
    file.path("..", "..", "data", "Gierlinski_2015", "Sleuth", "genes.tsv"),
    sep = "\t",
    header = TRUE
)
full <- full[target_id %in% swi_snf$transcript_id]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Comparing simulations to full dataset")
sims_comp <- merge(
    x = sims,
    y = full[, .SD, .SDcols = c("target_id", "pval", "qval", "b")],
    by = "target_id",
    suffixes = c("_small", "_full")
)

# calculate square error
sims_comp[, Sq_Err := (b_small - b_full) ^ 2]

# calculate MSE
err <- sims_comp[,
    .(
        Mean_SE = mean(Sq_Err),
        SD_SE = sd(Sq_Err)
    ),
    by = c("target_id", "Method")
]


# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")

# ensure that the directory exists
dest_dir <- file.path(RESULT_DIR, "Comparison", "SWI-SNF")
if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
}

fwrite(
    err,
    file.path(dest_dir, "err.tsv"),
    sep = "\t",
    col.names = TRUE
)
