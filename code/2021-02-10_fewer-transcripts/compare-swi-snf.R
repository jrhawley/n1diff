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
            file.path(RESULT_DIR, "Iterations", "SWI-SNF", paste0("iter_", i, ".balanced.genes.tsv")),
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
            file.path(RESULT_DIR, "Iterations", "SWI-SNF", paste0("iter_", i, ".unbalanced-ols.genes.tsv")),
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
            file.path(RESULT_DIR, "Iterations", "SWI-SNF", paste0("iter_", i, ".unbalanced-jse.genes.tsv")),
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

# merge full results with gene annotations
full <- merge(
    x = full,
    y = swi_snf,
    by.x = "target_id",
    by.y = "transcript_id"
)

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

# calculate error
sims_comp[, `:=`(
    Err = b_small - b_full,
    Abs_Err = abs(b_small - b_full),
    Sq_Err = (b_small - b_full) ^ 2,
    Frac_Err = (b_small - b_full) / b_full
)]

# recast to wide format for easier calculations
sims_comp_wide <- dcast(
    sims_comp[, .(
        Iteration,
        target_id,
        Method = gsub(" ", "_", Method),
        Sq_Err
    )],
    Iteration + target_id ~ Method,
    value.var = "Sq_Err"
)

# calculate percent difference in error between OLS and JSE
sims_comp_wide[, Delta_Err := (Unbalanced_JS - Unbalanced_OLS)]
sims_comp_wide[, Frac_Delta := (Unbalanced_JS - Unbalanced_OLS) / Unbalanced_OLS]

# calculate MSE
err <- sims_comp_wide[,
    .(
        Median_Delta_Err = median(Delta_Err, na.rm = TRUE),
        Mean_Delta_Err = mean(Delta_Err, na.rm = TRUE),
        SD_Delta_Err = sd(Delta_Err, na.rm = TRUE),
        Median_Frac_Delta = median(Frac_Delta, na.rm = TRUE),
        Mean_Frac_Delta = mean(Frac_Delta, na.rm = TRUE),
        SD_Frac_Delta = sd(Frac_Delta, na.rm = TRUE)
    )
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
    sims_comp_wide,
    file.path(dest_dir, "simulations.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    err,
    file.path(dest_dir, "err.tsv"),
    sep = "\t",
    col.names = TRUE
)
