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
comp <- fread(file.path(RESULT_DIR, "Comparison", "random", "simulations.tsv"))

# depending on which genes were filtered, not all 3 methods may have the same transcripts
comp <- comp[complete.cases(comp)]

# ==============================================================================
# Analysis
# ==============================================================================
# calculate MSE for each method
mse <- comp[,
	.(
		Balanced = mean(Balanced, na.rm = TRUE),
		Unbalanced_OLS = mean(Unbalanced_OLS, na.rm = TRUE),
		Unbalanced_JS = mean(Unbalanced_JS, na.rm = TRUE)
	),
	by = c("Iteration", "Total")
]

# convert to long format
mse_long <- melt(
	mse,
	id.vars = c("Iteration", "Total"),
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

# perform ANOVA test to estimate coefficients for each factor
htest <- aov(Rel_MSE ~ Total + Method + Iteration, data = mse_long[Method != "Balanced"])

# save coefficients
htest_summ <- as.data.table(
	unclass(summary(htest))[[1]],
	keep.rownames = TRUE
)
colnames(htest_summ) <- c(
	"Statistic",
	"df",
	"sum_sq",
	"mean_sq",
	"F",
	"pval"
)
fwrite(
	htest_summ,
	file.path(RESULT_DIR, "Comparison", "random", "anova.tsv"),
	sep = "\t",
	col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
# ensure that the directory exists
if (!dir.exists(PLOT_DIR)) {
	dir.create(PLOT_DIR, recursive = TRUE)
}

appender <- function(s) {
	paste(s, "Transcripts")
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
	+ facet_wrap(~ Total, scales = "free", labeller = as_labeller(appender))
	+ theme_minimal()
)
ggsave(
	file.path(PLOT_DIR, "random.unbalanced-comparison.png"),
	gg,
	width = 20,
	height = 12,
	units = "cm"
)

gg <- (
	ggplot(
		data = mse_long[Method != "Balanced"],
		aes(
			x = as.factor(Total),
			y = Rel_MSE,
			colour = Method,
			group = paste(Method, Total)
		)
	)
	+ geom_point(
		alpha = 0.4,
		position = position_jitterdodge(
			jitter.width = 0.2,
			jitter.height = 0,
			dodge.width = 0.75
		)
	)
	+ geom_boxplot(
		fill = NA,
		outlier.shape = NA
	)
	+ scale_x_discrete(
		name = "Number of Transcripts"
	)
	+ scale_y_continuous(
		name = "Mean Square Error\n(Relative to Balanced)",
	)
	+ theme_minimal()
	+ theme(
		legend.position = "top",
		axis.text.x = element_text(
			colour = "#000000",
			angle = 90,
			hjust = 1,
			vjust = 0.5
		),
		panel.grid.major.x = element_blank()
	)
)
ggsave(
	file.path(PLOT_DIR, "random.mse.all.png"),
	gg,
	width = 20,
	height = 12,
	units = "cm"
)
