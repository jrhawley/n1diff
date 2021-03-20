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
suppressWarnings(library("scales"))

RESULT_DIR <- file.path("..", "..", "results", "fewer-transcripts")
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
mse[, Frac_Delta := (Unbalanced_JS - Unbalanced_OLS) / Unbalanced_OLS]

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

htest_summ <- rbindlist(lapply(
	c(3, 10, 25, 50, 100, 250),
	function(i) {
		htest <- t.test(
			x = mse_long[(Total == i) & (Method == "Unbalanced_OLS"), Rel_MSE],
			y = mse_long[(Total == i) & (Method == "Unbalanced_JS"), Rel_MSE],
			alternative = "greater",
			paired = TRUE
		)
		return(data.table(
			"Total" = i,
			"t" = htest$statistic,
			"p" = htest$p.value,
			"shift" = htest$estimate
		))
	}
))
htest_summ[, q := p.adjust(p, method = "fdr")]

fwrite(
	htest_summ,
	file.path(RESULT_DIR, "Comparison", "random", "pairwise-tests.tsv"),
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

gg_comp <- (
	ggplot(
		data = mse[Total < 500],
		mapping = aes(
			x = Unbalanced_OLS,
			y = Unbalanced_JS,
			colour = as.factor(Total)
		)
	)
	+ geom_point()
	+ geom_abline(slope = 1, intercept = 0, linetype = "dashed")
	+ scale_x_continuous(
		name = "Mean Square Error (OLS)",
		limits = function(x) c(0, max(0.1, x))
	)
	+ scale_y_continuous(
		name = "Mean Square Error (JS)",
		limits = function(x) c(0, max(0.1, x))
	)
	+ scale_colour_viridis_d()
	+ guides(colour = FALSE)
	+ facet_wrap(
		~ Total,
		scales = "free",
		labeller = as_labeller(appender),
		ncol = 2
	)
	+ theme_minimal()
)
ggsave(
	file.path(PLOT_DIR, "random.unbalanced-comparison.png"),
	gg_comp,
	width = 8,
	height = 12,
	units = "cm"
)
ggsave(
	file.path(PLOT_DIR, "random.unbalanced-comparison.pdf"),
	gg_comp,
	width = 8,
	height = 12,
	units = "cm"
)

gg_diff <- (
	ggplot(
		data = mse[Total < 500],
		mapping = aes(
			x = as.factor(Total),
			y = 100 * Frac_Delta,
			colour = as.factor(Total)
		)
	)
	+ geom_point(
		alpha = 0.4,
		position = position_jitter(
			width = 0.2,
			height = 0
		)
	)
	+ geom_boxplot(
		fill = NA,
		outlier.shape = NA
	)
	+ geom_text(
		data = htest_summ,
		mapping = aes(
			x = as.factor(Total),
			y = 100,
			label = paste0(
				format(round(100 * shift, 2), nsmall = 2), "%",
				"\nq = ", scientific(p)
			)
		),
		vjust = -0.5
	)
	+ scale_x_discrete(
		name = "Number of Transcripts"
	)
	+ scale_y_continuous(
		name = "Mean Square Error\n(% Difference)",
		breaks = c(-100, -50, 0, 50, 100, 150),
		labels = paste(c(-100, -50, 0, 50, 100, 150), "%")
	)
	+ scale_colour_viridis_d()
	+ guides(colour = FALSE)
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
	file.path(PLOT_DIR, "random.mse.frac-delta.png"),
	gg_diff,
	width = 20,
	height = 12,
	units = "cm"
)
ggsave(
	file.path(PLOT_DIR, "random.mse.frac-delta.pdf"),
	gg_diff,
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
ggsave(
	file.path(PLOT_DIR, "random.mse.all.pdf"),
	gg,
	width = 20,
	height = 12,
	units = "cm"
)
