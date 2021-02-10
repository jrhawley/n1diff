# ==============================================================================
# Meta
# ==============================================================================
# sleuth
# ------------------------------------------------
# Author: James Hawley
# Description: Differential expression analysis with some samples


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))

# ==============================================================================
# Functions
# ==============================================================================
#' Calculate TP/FP/TN/FN, and derivative metrics, for a provided q-value threshold
#' dt: data.table test comparisons
#' qval_thresh: numeric q-value threshold
assess_jse <- function(dt, qval_thresh = 0.01) {
	# q-value threshold for ideal full results
	full_qval_thresh <- 0.01
	# classify points based on the threshold
	dt[(qval_full <= full_qval_thresh) & (qval_small <= qval_thresh), Result := "TP"]
	dt[(qval_full <= full_qval_thresh) & (qval_small > qval_thresh), Result := "FN"]
	dt[(qval_full > full_qval_thresh) & (qval_small <= qval_thresh), Result := "FP"]
	dt[(qval_full > full_qval_thresh) & (qval_small > qval_thresh), Result := "TN"]
	dt[, Result := factor(Result, levels = c("TP", "FN", "FP", "TN"))]

	# calculate confusion matrices
	confusion <- dt[, .N, by = c("Result", "Total", "Test_Condition", "Iteration")]
	# ensure that all 4 values (TP/FN/FP/TN) exist in the table for each combo of iteration, total, and balance
	setkey(confusion, Total, Test_Condition, Iteration, Result)
	## uses a cross-join to calculate all combinations of the four columns
	confusion <- confusion[
		CJ(
			unique(Total),
			levels(Test_Condition),
			unique(Iteration),
			levels(Result)
		),
		allow.cartesian = TRUE
	]
	confusion[is.na(N), N := 0L]

	# calculate various rates from the confusion matrices
	rates <- dcast(
		confusion,
		Total + Test_Condition + Iteration ~ Result,
		value.var = "N"
	)

	rates[,
		`:=`(
			TPR = TP / (TP + FN),
			TNR = TN / (TN + FP),
			PPV = TP / (TP + FP),
			NPV = TN / (TN + FN),
			ACC = (TP + TN) / (TP + TN + FP + FN),
			BA = ( TP / (TP + FN) + TN / (TN + FP) ) / 2,
			F1 = 2 * TP / (2 * TP + FP + FN),
			MCC = (TP * TN - FP * FN) / (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))
		)
	]

	# compare methods with hypothesis testing
	unbal_rates <- rates[Test_Condition != "Balanced"]
	unbal_rates[, Total := factor(Total)]
	unbal_rates[, Test_Condition := factor(Test_Condition)]

	comp_unbal <- list(
		"TPR" = aov(
			formula = TPR ~ Test_Condition + Total,
			data = unbal_rates,
			contrasts = list("Test_Condition" = contr.treatment(c("OLS", "JS")))
		),
		"TNR" = aov(
			formula = TNR ~ Test_Condition + Total,
			data = unbal_rates,
			contrasts = list("Test_Condition" = contr.treatment(c("OLS", "JS")))
		),
		"PPV" = aov(
			formula = PPV ~ Test_Condition + Total,
			data = unbal_rates,
			contrasts = list("Test_Condition" = contr.treatment(c("OLS", "JS")))
		),
		"NPV" = aov(
			formula = NPV ~ Test_Condition + Total,
			data = unbal_rates,
			contrasts = list("Test_Condition" = contr.treatment(c("OLS", "JS")))
		),
		"ACC" = aov(
			formula = ACC ~ Test_Condition + Total,
			data = unbal_rates,
			contrasts = list("Test_Condition" = contr.treatment(c("OLS", "JS")))
		),
		"BA" = aov(
			formula = BA ~ Test_Condition + Total,
			data = unbal_rates,
			contrasts = list("Test_Condition" = contr.treatment(c("OLS", "JS")))
		),
		"MCC" = aov(
			formula = MCC ~ Test_Condition + Total,
			data = unbal_rates,
			contrasts = list("Test_Condition" = contr.treatment(c("OLS", "JS")))
		)
	)

	comp_unbal_stats <- rbindlist(lapply(
		names(comp_unbal),
		function(s) {
			data.table(
				"Statistic" = s,
				# extracting the coefficient
				# the alphabetical ordering has the contrast set to when Test_ConditionJS == 0
				# we have a James-Stein estimate. When Test_ConditionJS == 1, it's the OLS.
				# To get the effect of the James-Stime estimate, we flip the sign of the coefficient
				# to have the OLS method as the baseline
				"Value" = -comp_unbal[[s]]$coefficients["Test_ConditionJS"],
				# extracting the p-value from the summary
				# it's a little convoluted, but Test_Condition1 is the first row of the summary table
				# this corresponds to the Test_Condition factor (aka the method)
				"pval" = summary(comp_unbal[[s]])[[1]][["Pr(>F)"]][1]
			)
		}
	))
	comp_unbal_stats[, qval := p.adjust(pval, method = "fdr")]

	# 2. Calculate MSE between fold-change estimates themselves
	# ------------------------------------------------
	mse <- dt[,
		.(
			MSE = mean(Sq_Err),
			SE_SD = sd(Sq_Err)
		),
		by = c("target_id", "Test_Condition")
	]
	mse_all <- dt[,
		.(
			MSE = mean(Sq_Err),
			SE_SD = sd(Sq_Err)
		),
		by = "Test_Condition"
	]

	return(list(
		"qval_thresh" = qval_thresh,
		"results" = dt,
		"rates" = rates,
		"confusion" = confusion,
		"jse_comp" = comp_unbal,
		"coefficients" = comp_unbal_stats,
		"mse" = mse,
		"mse_all" = mse_all
	))
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# replace "wt" with "ctrl" to ensure that its the baseline control, alphabetically
meta[condition == "wt", condition := "ctrl"]

# load full scale calculations
full <- fread(
	file.path("..", "..", "data", "Gierlinski_2015", "Sleuth", "genes.tsv"),
	sep = "\t",
	header = TRUE
)

# load sample calculations
total_reps <- 30

tests <- rbindlist(lapply(
	seq(4, 24, 2),
	function(total) {
		bal <- rbindlist(lapply(
			1:total_reps,
			function(i) {
				dt2 <- fread(
					file.path("Iterations", paste0("Total_", total), "balanced", paste0(i, ".genes.tsv")),
					sep = "\t",
					header = TRUE
				)
				dt2[, Iteration := i]
				return(dt2)
			}
		))
		bal[, Test_Condition := "Balanced"]
		unbal_naive <- rbindlist(lapply(
			1:total_reps,
			function(i) {
				dt2 <- fread(
					file.path("Iterations", paste0("Total_", total), "unbalanced", paste0(i, ".genes.tsv")),
					sep = "\t",
					header = TRUE
				)
				dt2[, Iteration := i]
				return(dt2)
			}
		))
		unbal_naive[, Test_Condition := "Unbalanced OLS"]
		unbal_jse <- rbindlist(lapply(
			1:total_reps,
			function(i) {
				dt2 <- fread(
					file.path("Iterations", paste0("Total_", total), "unbalanced", paste0(i, ".genes.jse.tsv")),
					sep = "\t",
					header = TRUE
				)
				dt2[, Iteration := i]
				return(dt2)
			}
		))
		unbal_jse[, Test_Condition := "Unbalanced JS"]
		combined_dt <- rbindlist(
			list(
				bal,
				unbal_naive,
				unbal_jse
			),
			fill = TRUE
		)
		combined_dt[, Total := total]
		return(combined_dt)
	}
))


# ==============================================================================
# Analysis
# ==============================================================================

# 1. Calculate accuracy and other measures of TP/FP/TN/FN between methods at a given q-value threshold
# ------------------------------------------------
loginfo("Comparing methods")

# compare to the full dataset
test_comparisons <- merge(
	x = tests[, .SD, .SDcols = c("target_id", "b", "pval", "qval", "Total", "Iteration", "Test_Condition")],
	y = full[, .SD, .SDcols = c("target_id", "b", "pval", "qval")],
	by = "target_id",
	suffixes = c("_small", "_full")
)
# force to be a factor for easier calculations later
test_comparisons[, Test_Condition := factor(
	Test_Condition,
	levels = c("Balanced", "Unbalanced OLS", "Unbalanced JS")
)]
# calculate square-error of estimates from the smaller samples
test_comparisons[, Sq_Err := (b_full - b_small) ^ 2]

# label the result as "TP", "FN", "TN", "FP" depending on the q-value of the full comparison
# this is a "methodologically true/false", not a "biologically true/false"
iterations <- lapply(
	seq(0, 0.99, 0.01),
	function(q) {
		assess_jse(test_comparisons, q)
	}
)

# calculate receiver-operater characteristic
roc <- rbindlist(lapply(
	iterations,
	function(it) {
		dt <- it$rates[,
			.(
				TPR = mean(TPR),
				FPR = 1 - mean(TNR)
			),
			by = c("Total", "Test_Condition")
		]
		dt[, q := it$qval_thresh]
		return(dt)
	}
))

# calculate precision-recall curve
prc <- rbindlist(lapply(
	iterations,
	function(it) {
		dt <- it$rates[,
			.(
				PPV = mean(PPV),
				TPR = mean(TPR)
			),
			by = c("Total", "Test_Condition")
		]
		dt[, q := it$qval_thresh]
		return(dt)
	}
))

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

# write results with the q-value threshold of 0.01
fwrite(
	iterations[[2]]$results,
	file.path("Comparison", "tests.tsv"),
	sep = "\t",
	col.names = TRUE
)
fwrite(
	iterations[[2]]$confusion,
	file.path("Comparison", "confusion.tsv"),
	sep = "\t",
	col.names = TRUE
)
fwrite(
	iterations[[2]]$rates,
	file.path("Comparison", "rates.tsv"),
	sep = "\t",
	col.names = TRUE
)
fwrite(
	iterations[[2]]$coefficients,
	file.path("Comparison", "jse-comp-coeffs.tsv"),
	sep = "\t",
	col.names = TRUE
)

fwrite(
	iterations[[2]]$mse,
	file.path("Comparison", "mse.by-transcript.tsv"),
	sep = "\t",
	col.names = TRUE
)
fwrite(
	iterations[[2]]$mse_all,
	file.path("Comparison", "mse.all.tsv"),
	sep = "\t",
	col.names = TRUE
)

fwrite(
	roc,
	file.path("Comparison", "roc.tsv"),
	sep = "\t",
	col.names = TRUE
)

fwrite(
	prc,
	file.path("Comparison", "prc.tsv"),
	sep = "\t",
	col.names = TRUE
)

