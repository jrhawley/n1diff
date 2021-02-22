# ==============================================================================
# Meta
# ==============================================================================
# plot-random
# ------------------------------------------------
# Author: James Hawley
# Description: Plot statistics and results from differential expression analyses with randomly selected transcripts


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
comp <- fread(file.path(RESULT_DIR, "Comparison", "random", "err.tsv"))

# depending on which genes were filtered, not all 3 methods may have the same transcripts
comp <- comp[complete.cases(comp)]


# ==============================================================================
# Analysis
# ==============================================================================
# convert to wide format for plotting
comp_wide_mean <- dcast(
    # replace spaces in Method column with underscores on the fly for easier column access after casting
    comp[, .(
        target_id,
        Total,
        Method = gsub(" ", "_", Method),
        Mean_SE
    )],
    target_id + Total ~ Method,
    value.var = "Mean_SE"
)
comp_wide_sd <- dcast(
    # replace spaces in Method column with underscores on the fly for easier column access after casting
    comp[, .(
        target_id,
        Total,
        Method = gsub(" ", "_", Method),
        SD_SE
    )],
    target_id + Total ~ Method,
    value.var = "SD_SE"
)

# merge Mean and SD calculations into the same table
comp_wide <- merge(
    x = comp_wide_mean,
    y = comp_wide_sd,
    suffixes = c("_mean", "_sd")
)
comp_wide <- comp_wide[complete.cases(comp_wide)]

# calculate relevant statistics
err <- 

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
# ensure that the directory exists
if (!dir.exists(PLOT_DIR)) {
    dir.create(PLOT_DIR, recursive = TRUE)
}

gg <- (
    ggplot(data = comp_wide_mean)
    + geom_point(aes(x = Unbalanced_OLS, y = Unbalanced_JS))
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    + facet_wrap(~ Total)
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "random.unbalanced-comparison.png"),
    gg,
    width = 12,
    height = 8,
    units = "cm"
)

gg <- (
    ggplot(data = comp)
    + geom_col(aes(x = Method, y = Mean_SE, fill = Method))
    + scale_x_discrete(
        name = NULL
    )
    + scale_y_continuous(
        name = "Fold-Change MSE"
    )
    + facet_wrap(~ Total, scales = "free_y")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1, vjust = 0.5)
    )
)
ggsave(
    file.path(PLOT_DIR, "random.mse.png"),
    gg,
    width = 16,
    height = 10,
    units = "cm"
)

gg <- (
    ggplot(data = comp_wide_mean)
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
    file.path(PLOT_DIR, "mse.percent-change.png"),
    gg,
    width = 16,
    height = 10,
    units = "cm"
)
