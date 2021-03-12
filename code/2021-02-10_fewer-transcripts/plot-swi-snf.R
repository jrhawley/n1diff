# ==============================================================================
# Meta
# ==============================================================================
# plot-swi-snf
# ------------------------------------------------
# Author: James Hawley
# Description: Plot statistics and results from SWI/SNF differential expression analyses


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))

RESULT_DIR <- file.path("..", "..", "results", "2021-02-10_fewer-transcripts")
PLOT_DIR <- file.path(RESULT_DIR, "Plots")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load comparison data
err <- fread(file.path(RESULT_DIR, "Comparison", "SWI-SNF", "err.tsv"))
comp <- fread(file.path(RESULT_DIR, "Comparison", "SWI-SNF", "simulations.tsv"))

# load SWI/SNF complex information
swi_snf <- fread("swi-snf-complex-subunits.tsv")

# ==============================================================================
# Analysis
# ==============================================================================
# add human-friendly gene names
comp <- merge(
    x = comp,
    y = swi_snf[, .SD, .SDcols = c("transcript_id", "gene_id", "gene_name")],
    by.x = "target_id",
    by.y = "transcript_id"
)


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
# ensure that the directory exists
if (!dir.exists(PLOT_DIR)) {
    dir.create(PLOT_DIR, recursive = TRUE)
}

gg <- (
    ggplot(data = comp)
    + geom_point(aes(x = Unbalanced_OLS, y = Unbalanced_JS))
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "swi-snf.unbalanced-comparison.png"),
    gg,
    width = 12,
    height = 8,
    units = "cm"
)

# gg <- (
#     ggplot(data = comp)
#     + geom_col(aes(x = Method, y = Mean_SE, fill = Method))
#     + scale_x_discrete(
#         name = NULL
#     )
#     + scale_y_continuous(
#         name = "Fold-Change MSE"
#     )
#     + facet_wrap(~ gene_name, scales = "free_y")
#     + guides(fill = FALSE)
#     + theme_minimal()
#     + theme(
#         axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1, vjust = 0.5)
#     )
# )
# ggsave(
#     file.path(PLOT_DIR, "swi-snf.mse.png"),
#     gg,
#     width = 16,
#     height = 10,
#     units = "cm"
# )

gg <- (
    ggplot(data = comp)
    + geom_col(aes(x = gene_name, y = 100 * (Unbalanced_JS - Unbalanced_OLS) / Unbalanced_OLS))
    + scale_x_discrete(
        name = NULL
    )
    + scale_y_continuous(
        name = "MSE % Difference",
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank()
    )
)
ggsave(
    file.path(PLOT_DIR, "swi-snf.mse.percent-change.png"),
    gg,
    width = 16,
    height = 10,
    units = "cm"
)
