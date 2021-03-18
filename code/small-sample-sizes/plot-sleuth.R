# ==============================================================================
# Meta
# ==============================================================================
# plot-sleuth
# ------------------------------------------------
# Author: James Hawley
# Description: Plot comparisons between full scale dataset and smaller sample sizes

# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

DATA_DIR <- file.path("..", "..", "data", "Gierlinski_2015", "Sleuth")
RESULT_DIR <- file.path("..", "..", "results", "small-sample-sizes")
PLOT_DIR <- file.path(RESULT_DIR, "Plots")
# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load data
confusion <- fread(
	file.path(RESULT_DIR, "Comparison", "confusion.tsv"),
	sep = "\t",
	header = TRUE
)
rates <- fread(
	file.path(RESULT_DIR, "Comparison", "rates.tsv"),
	sep = "\t",
	header = TRUE
)

mse_all <- fread(
	file.path(RESULT_DIR, "Comparison", "mse.all.tsv"),
	sep = "\t",
	header = TRUE
)

roc <- fread(
	file.path(RESULT_DIR, "Comparison", "roc.tsv"),
	sep = "\t",
	header = TRUE
)

prc <- fread(
	file.path(RESULT_DIR, "Comparison", "prc.tsv"),
	sep = "\t",
	header = TRUE
)

# load full dataset
full <- fread(
	file.path(DATA_DIR, "genes.tsv"),
	sep = "\t",
	header = TRUE
)

# calculate PPV of an uninformative binary classifier
full_qval_thresh <- 0.01
uninform_ppv <- full[qval <= full_qval_thresh, .N] / full[, .N]

# data table for making plots more intuitive
stat_better_direction <- data.table(
	Stat = c(
		"ACC",
		"BA",
		"MCC",
		"NPV",
		"PPV",
		"TNR",
		"TPR"
	),
	Direction = c(
		"higher",
		"higher",
		"higher",
		"higher",
		"higher",
		"higher",
		"higher"
	)
)


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
for (s in stat_better_direction$Stat) {
	loginfo(s)

	better_dir <- stat_better_direction[Stat == s, Direction]

	gg <- (
		ggplot(
			data = rates,
			aes(x = Total, y = get(s), colour = Test_Condition)
		)
		+ geom_boxplot(
			mapping = aes(
				fill = Test_Condition,
				group = paste(Total, Test_Condition)
			),
			colour = "#000000",
			alpha = 0.5,
			outlier.shape = 21,
			position = position_dodge(width = 1.5)
		)
		+ scale_x_continuous(
			name = "Total Number of Samples",
			breaks = seq(4, 24, 2),
			labels = seq(4, 24, 2)
		)
		+ scale_y_continuous(
			name = paste0(s," (", better_dir, " is better)")
			# limits = c(bounds[1], bounds[2])
		)
		+ scale_colour_manual(
			name = "Design",
			breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
			labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
			values = c("#4B4B4B", "#DA8FD6", "#FA8072")
		)
		+ scale_fill_manual(
			name = "Design",
			breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
			labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
			values = c("#4B4B4B", "#DA8FD6", "#FA8072")
		)
		+ theme_minimal()
		+ theme(
			panel.grid.minor = element_blank(),
			legend.position = "bottom",
			axis.text.x = element_text(
				colour = "#000000",
				angle = 90
			),
			axis.text.y = element_text(
				colour = "#000000"
			)
		)
	)
	ggsave(
		file.path(PLOT_DIR, paste0(s, ".png")),
		width = 12,
		height = 8,
		units = "cm"
	)
}

loginfo("Confusion matrix")
gg <- (
	ggplot(
		data = confusion,
		aes(x = Total, y = N, colour = Test_Condition)
	)
	+ geom_boxplot(
		mapping = aes(
			fill = Test_Condition,
			group = paste(Total, Test_Condition)
		),
		colour = "#000000",
		alpha = 0.5,
		outlier.shape = 21,
		position = position_dodge(width = 1.5)
	)
	+ scale_x_continuous(
		name = "Total Number of Samples",
		breaks = seq(4, 24, 2),
		labels = seq(4, 24, 2)
	)
	+ scale_y_continuous(
		name = "Frequency"
		# limits = c(bounds[1], bounds[2])
	)
	+ scale_colour_manual(
		name = "Design",
		breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
		labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
		values = c("#4B4B4B", "#DA8FD6", "#FA8072")
	)
	+ scale_fill_manual(
		name = "Design",
		breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
		labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
		values = c("#4B4B4B", "#DA8FD6", "#FA8072")
	)
	+ facet_wrap(~ Result, scales = "free_y")
	+ theme_minimal()
	+ theme(
		panel.grid.minor = element_blank(),
		legend.position = "bottom",
		axis.text.x = element_text(
			colour = "#000000",
			angle = 90
		),
		axis.text.y = element_text(
			colour = "#000000"
		)
	)
)
ggsave(
	file.path(PLOT_DIR, "confusion.png"),
	width = 12,
	height = 8,
	units = "cm"
)

loginfo("MSE")
gg_mse_all <- (
	ggplot(
		data = mse_all,
		mapping = aes(x = Test_Condition, y = MSE, fill = Test_Condition)
		)
	+ geom_col()
	+ scale_x_discrete(
		name = NULL,
		breaks = c("Balanced", "Unbalanced JS", "Unbalanced OLS"),
		labels = c("Even split", "1-vs-N James-Stein", "1-vs-N")
	)
	+ scale_y_continuous(
		name = "Mean Square Error"
	)
	+ scale_fill_manual(
		name = "Design",
		breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
		labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
		values = c("#4B4B4B", "#DA8FD6", "#FA8072")
	)
	+ guides(fill = FALSE)
	+ theme_minimal()
	+ theme(
		axis.text.x = element_text(
			colour = "#000000",
			angle = 0
		),
		axis.text.y = element_text(
			colour = "#000000"
		)
	)
)
ggsave(
	file.path(PLOT_DIR, "mse.all.png"),
	width = 12,
	height = 8,
	units = "cm"
)

loginfo("ROC")
gg <- (
	ggplot(data = roc)
	+ geom_path(
		aes(
			x = FPR,
			y = TPR,
			colour = Test_Condition
		)
	)
	+ geom_abline(
		slope = 1,
		intercept = 0,
		linetype = "dashed",
		colour = "#000000"
	)
	+ facet_wrap(~ Total)
	+ scale_x_continuous(
		name = "False positive rate",
		limits = c(0, 1)
	)
	+ scale_y_continuous(
		name = "True positive rate",
		limits = c(0, 1)
	)
	+ theme_minimal()
	+ theme(
		legend.position = "bottom",
		axis.text.x = element_text(
			colour = "#000000",
			angle = 90
		),
		axis.text.y = element_text(
			colour = "#000000"
		)
	)
)
ggsave(
	file.path(PLOT_DIR, "auroc.png"),
	width = 12,
	height = 12,
	units = "cm"
)

loginfo("PRC")
gg <- (
	ggplot(data = prc)
	+ geom_path(
		aes(
			x = TPR,
			y = PPV,
			colour = Test_Condition
		)
	)
	+ geom_abline(
		slope = 0,
		intercept = uninform_ppv,
		linetype = "dashed",
		colour = "#000000"
	)
	+ facet_wrap(~ Total)
	+ scale_x_continuous(
		name = "True positive rate",
		limits = c(0, 1)
	)
	+ scale_y_continuous(
		name = "Positive predictive value",
		limits = c(0, 1)
	)
	+ theme_minimal()
	+ theme(
		legend.position = "bottom",
		axis.text.x = element_text(
			colour = "#000000",
			angle = 90
		),
		axis.text.y = element_text(
			colour = "#000000"
		)
	)
)
ggsave(
	file.path(PLOT_DIR, "auprc.png"),
	width = 12,
	height = 12,
	units = "cm"
)

