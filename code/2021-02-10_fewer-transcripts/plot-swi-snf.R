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

mse <- comp[,
    .(
        Balanced = mean(Balanced, na.rm = TRUE),
        Unbalanced_JS = mean(Unbalanced_JS, na.rm = TRUE),
        Unbalanced_OLS = mean(Unbalanced_OLS, na.rm = TRUE)
    ),
    by = "Iteration"
]

# convert to long format
mse_long <- melt(
	mse,
	id.vars = "Iteration",
	variable.name = "Method",
	value.name = "MSE"
)

# calculate MSE relative to the mean balanced MSE
mean_balanced_mse <- mse_long[
	Method == "Balanced",
	.(Mean_Balanced = mean(MSE)),
	by = "Total"
]

mse_long <- merge(
	x = mse_long,
	y = mean_balanced_mse,
	by = "Total"
)
mse_long[, Rel_MSE := MSE / Mean_Balanced]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
# ensure that the directory exists
if (!dir.exists(PLOT_DIR)) {
    dir.create(PLOT_DIR, recursive = TRUE)
}

gg <- (
	ggplot(data = mse)
	+ geom_point(aes(x = Unbalanced_OLS, y = Unbalanced_JS))
	+ geom_abline(slope = 1, intercept = 0, linetype = "dashed")
	+ scale_x_continuous(
		name = "Mean Square Error (OLS)"
	)
	+ scale_y_continuous(
		name = "Mean Square Error (JS)"
	)
	# + facet_wrap(~ Total, scales = "free", labeller = as_labeller(appender))
	+ theme_minimal()
)
ggsave(
	file.path(PLOT_DIR, "swi-snf.unbalanced-comparison.png"),
	gg,
	width = 12,
	height = 8,
	units = "cm"
)
