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
	file.path("..", "..", "Data", "Gierlinski_2015", "Sleuth", "genes.tsv"),
	sep = "\t",
	header = TRUE
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Compare to full scale analyses")

# load sample calculations
tests <- rbindlist(lapply(
	seq(4, 24, 2),
	function(total) {
		dt <- rbindlist(lapply(
			c("unbalanced", "balanced"),
			function(bal) {
				dt2 <- rbindlist(lapply(
					1:10,
					function(i) {
						dt3 <- fread(
							file.path("Iterations", paste0("Total_", total), bal, paste0(i, ".genes.tsv")),
							sep = "\t",
							header = TRUE
						)
						dt3[, Iteration := i]
						return(dt3)
					}
				))
				dt2[, Balanced := bal]
				return(dt2)
			}
		))
		dt[, Total := total]
		return(dt)
	}
))

# compare to the full dataset
test_comparisons <- merge(
	x = tests[, .SD, .SDcols = c("target_id", "b", "pval", "qval", "Total", "Iteration", "Balanced")],
	y = full[, .SD, .SDcols = c("target_id", "b", "pval", "qval")],
	by = "target_id",
	suffixes = c("_small", "_full")
)

# label the result as "TP", "FN", "TN", "FP" depending on the q-value of the full comparison
# this is a "methodologically true/false", not a "biologically true/false"
qval_thresh <- 0.01
test_comparisons[(qval_full <= qval_thresh) & (qval_small <= qval_thresh), Result := "TP"]
test_comparisons[(qval_full <= qval_thresh) & (qval_small > qval_thresh), Result := "FN"]
test_comparisons[(qval_full > qval_thresh) & (qval_small <= qval_thresh), Result := "FP"]
test_comparisons[(qval_full > qval_thresh) & (qval_small > qval_thresh), Result := "TN"]
test_comparisons[, Balanced := factor(Balanced, levels = c("balanced", "unbalanced"))]
test_comparisons[, Result := factor(Result, levels = c("TP", "FN", "FP", "TN"))]

# calculate confusion matrices
confusion <- test_comparisons[, .N, by = c("Result", "Total", "Balanced", "Iteration")]
# ensure that all 4 values (TP/FN/FP/TN) exist in the table for each combo of iteration, total, and balance
setkey(confusion, Total, Balanced, Iteration, Result)
## uses a cross-join to calculate all combinations of the four columns
confusion <- confusion[CJ(unique(Total), levels(Balanced), unique(Iteration), levels(Result)), allow.cartesian = TRUE]
confusion[is.na(N), N := 0L]

# calculate various rates from the confusion matrices
rates <- dcast(
	confusion,
	Total + Balanced + Iteration ~ Result,
	value.var = "N"
)

rates[,
	`:=`(
		TPR = TP / (TP + FN),
		TNR = TN / (TN + FP),
		PPV = TP / (TP + FP),
		NPV = TN / (TN + FN),
		FNR = FN / (FN + TP),
		FDR = FP / (FP + TP),
		FOR = FN / (FN + TN),
		ACC = (TP + TN) / (TP + TN + FP + FN),
		F1 = 2 * TP / (2 * TP + FP + FN),
		MCC = (TP * TN - FP * FN) / (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))
	)
]


# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
	rates,
	file.path("Comparison", "rates.tsv"),
	sep = "\t",
	col.names = TRUE
)

