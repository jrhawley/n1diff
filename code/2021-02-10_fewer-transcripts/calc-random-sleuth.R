# ==============================================================================
# Meta
# ==============================================================================
# generic-small-sleuth
# ------------------------------------------------
# Author: James Hawley
# Description: Differential expression analysis for randomly selected transcripts


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("argparse"))

if (!interactive()) {
	parser <- ArgumentParser(description = "Differential expression analysis for randomly selected transcripts")
	parser$add_argument(
		"n",
		type = "integer",
		help = "Number of transcripts to sample"
	)
	cli_args <- parser$parse_args()
} else {
	cli_args <- list(
		n = 3
	)
}

suppressWarnings(library("logging"))

loginfo("Loading packages")
suppressWarnings(library("data.table"))
suppressWarnings(library("sleuth"))
source(file.path("..", "jse-shrinkage", "jse.R"))

RESULT_DIR <- file.path("..", "..", "results", "2021-02-10_fewer-transcripts")

# randomly select groups this many times before performing comparison
TOTAL_REPS <- 30


# ==============================================================================
# Functions
# ==============================================================================
# Helper function to randomly select a total number of samples and return the design matrices
select_samples <- function(metadata, total) {
	if (total %% 2 != 0) {
		stop("`total` argument must be even.")
	}
	half_total <- total / 2
	# pick balanced samples
	balanced <- lapply(
		1:2,
		function(i) {
			mu_ids <- metadata[condition == "mu", sample(sample, half_total)]
			ctrl_ids <- metadata[condition == "ctrl", sample(sample, half_total)]
			return(rbindlist(list(
				metadata[sample %in% mu_ids],
				metadata[sample %in% ctrl_ids]
			)))
		}
	)
	# pick unbalanced samples 1 mu vs total - 1 ctrls
	usm_mu_ids <- metadata[condition == "mu", sample(sample, 1)]
	usm_ctrl_ids <- metadata[condition == "ctrl", sample(sample, total - 1)]
	unbalanced_single_mu <- rbindlist(list(
		metadata[sample %in% usm_mu_ids],
		metadata[sample %in% usm_ctrl_ids]
	))
	# pick unbalanced samples 1 ctrl vs total - 1 mus
	usf_mu_ids <- metadata[condition == "mu", sample(sample, total - 1)]
	usf_ctrl_ids <- metadata[condition == "ctrl", sample(sample, 1)]
	unbalanced_single_ctrl <- rbindlist(list(
		metadata[sample %in% usf_mu_ids],
		metadata[sample %in% usf_ctrl_ids]
	))

	# aggregate and return entire list
	return(list(
		list(
			"balanced" = balanced[[1]],
			"unbalanced" = unbalanced_single_ctrl
		),
		list(
			"balanced" = balanced[[2]],
			"unbalanced" = unbalanced_single_mu
		)
	))
}

# helper function for verbose printing
vprint <- function(s, prefix = NULL) {
	loginfo(paste(prefix, s))
}

# Helper function for determining which is the single sample in an unbalanced design
which_single <- function(small_meta) {
	single_cdn <- small_meta[, .N, by = "condition"][N == 1, condition]
	single_sample <- small_meta[condition == single_cdn, sample]
	return(single_sample)
}

# Helper function for selecting random set of transcripts to assess
sample_transcripts <- function(obj, n_tx = 3, limit_to = NULL) {
	# only consider transcripts that have not been filtered out due to low read counts or other QC filters
	possible_tx <- names(which(obj$filter_bool))

	# specify if the total set of transcripts to choose from needs to be limited
	if (!is.null(limit_to)) {
		possible_tx <- intersect(limit_to, possible_tx)
	}

	# return the selected set of transcripts
	return(sample(possible_tx, n_tx, replace = FALSE))
}

# Helper function to perform differential analysis
diff_ge <- function(small_meta, iter_idx = "") {
	iter_prefix <- paste0("(", iter_idx, " / ", TOTAL_REPS, ")")
	so <- sleuth_prep(
		small_meta,
		extra_bootstrap_summary = TRUE,
		num_cores = 8,
		# filter out transcripts with < 10 reads in any sample
		filter_fun = function(row, min_reads = 10, min_prop = 0.83) { mean(row >= min_reads) >= min_prop }
	)
	vprint("Fitting model and hypothesis testing", iter_prefix)
	so <- sleuth_fit(so, ~condition, "full")

	# perform differential analysis
	so <- sleuth_wt(so, "conditionmu")

	# extract results
	so_genes <- as.data.table(sleuth_results(
		so,
		"conditionmu",
		"wt",
		show_all = FALSE
	))

	# return this set of items for future saving
	return(list(
		"object" = so,
		"results" = so_genes,
		"design" = small_meta
	))
}

# Helper function for performing James-Stein shrinkage
perform_jse <- function(so, small_meta, subset_tx = NULL, iter_idx = "") {
	iter_prefix <- paste0("(", iter_idx, " / ", TOTAL_REPS, ")")
	vprint("Calculating James-Stein estimates", iter_prefix)
	# specifiy the sample with the single condition
	single_sample <- which_single(small_meta)
	so_jse <- jse_shrinkage(
		so,
		which_sample = single_sample,
		which_beta = "conditionmu",
		# only select transcripts in the desired subset
		s_targets = subset_tx
	)
	so_jse_results <- jse_wald_test(so_jse)

	# return this set of items for future saving
	return(list(
		"object" = so_jse,
		"results" = so_jse_results,
		"design" = small_meta
	))
}

# Helper function for saving differential analysis data
save_dge_data <- function(obj, res, design, prefix) {
	# ensure that the directory exists
	dest_dir <- dirname(prefix)
	if (!dir.exists(dest_dir)) {
		dir.create(dest_dir, recursive = TRUE)
	}

	# save sleuth object
	saveRDS(
		obj,
		paste0(prefix, ".sleuth-object.rds")
	)
	# save differential analysis results
	fwrite(
		res,
		paste0(prefix, ".genes.tsv"),
		sep = "\t",
		col.names = TRUE
	)
	# save experimental design for troubleshooting
	fwrite(
		design,
		paste0(prefix, ".design.tsv"),
		sep = "\t",
		col.names = TRUE
	)
}

# Helper function for determining if everything in the James-Stein shrinkage worked properly
is_jse_error <- function(jse_obj) {
	return(
		is.na(jse_obj$object$shrinkage_coef)
		|| (any(is.na(jse_obj$results$b)))
		|| (any(is.na(jse_obj$results$se_b)))
	)
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# replace "wt" with "ctrl" to ensure that its the baseline control, alphabetically
meta[condition == "wt", condition := "ctrl"]

# load full dataset to determine the total set of transcripts to select from
full <- fread(
	file.path("..", "..", "data", "Gierlinski_2015", "Sleuth", "genes.tsv"),
	sep = "\t",
	header = TRUE
)
allowed_tx <- full$target_id

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# set random seed for reproducible results
set.seed(10000 + cli_args$n)

for (i in 1:(TOTAL_REPS / 2)) {
	loginfo("\tPreparing sleuth objects")
	# balanced samples
	sampled_metadata <- select_samples(meta, 6)
	for (j in 1:2) {
		iter_idx <- (i - 1) * 2 + j
		# perform differential analysis of balanced experimental design
		vprint("Balanced", iter_idx)
		comp_bal <- diff_ge(
			small_meta = sampled_metadata[[j]][["balanced"]],
			iter_idx = iter_idx
		)

		# set default values to initiate while loop
		subset_tx <- NA
		comp_unbal_jse <- list(
			"object" = list(
				"shrinkage_coef" = NA
			),
			"results" = list(
				"b" = NA,
				"se_b" = NA
			)
		)

		# check that
		# 1. the balanced and unbalanced designs share >= cli_args$n transcripts that aren't filtered out
		# 2. the trace-lambda condition for James-Stein shrinkage is met
		while_counter <- 0
		while (is.na(subset_tx) || is_jse_error(comp_unbal_jse)) {
			vprint(while_counter)
			# perform differential analysis of unbalanced experimental design (only OLS estimates)
			vprint("Unbalanced OLS", iter_idx)
			comp_unbal <- diff_ge(
				small_meta = sampled_metadata[[j]][["unbalanced"]],
				iter_idx = iter_idx
			)

			# limit analyses to randomly select from transcripts that weren't filtered out in both the balanced and unbalanced design
			# this ensures that all iterations of the same index will use the same transcripts
			# (i.e. paired experimental comparison)
			allowed_filtered_tx <- intersect(
				intersect(
					names(which(comp_bal$object$filter_bool)),
					names(which(comp_unbal$object$filter_bool))
				),
				allowed_tx
			)
			if (length(allowed_filtered_tx) >= cli_args$n) {
				subset_tx <- sample(allowed_filtered_tx, cli_args$n, replace = FALSE)
			}

			vprint("Unbalanced JS", iter_idx)
			comp_unbal_jse <- perform_jse(
				so = comp_unbal$object,
				small_meta = sampled_metadata[[j]][["unbalanced"]],
				subset_tx = subset_tx,
				iter_idx = iter_idx
			)
			# need to re-select samples if
			# 1. trace-lambda condition isn't met
			# 2. any transcript has an NA as an folc change estimate
			# 3. any transcript has an NA as a variance estimate
			if (is_jse_error(comp_unbal_jse)) {
				vprint("Failed shrinkage, trying with new samples", iter_idx)
				# if not, re-select samples to work with for unbalanced experimental design
				sampled_metadata[[j]][["unbalanced"]] <- select_samples(meta, 6)[[j]][["unbalanced"]]
			}
			while_counter <- while_counter + 1
		}

		vprint("Saving data", iter_idx)
		# save balanced design
		save_dge_data(
			comp_bal$object,
			comp_bal$results[target_id %in% subset_tx, .SD, keyby = "target_id"],
			comp_bal$design,
			file.path(
				RESULT_DIR,
				"Iterations",
				"random",
				paste0("total_", cli_args$n),
				paste0("iter_", iter_idx, ".balanced")
			)
		)
		# save unbalanced design with OLS results
		save_dge_data(
			comp_unbal$object,
			comp_unbal$results[target_id %in% subset_tx, .SD, keyby = "target_id"],
			comp_unbal$design,
			file.path(
				RESULT_DIR,
				"Iterations",
				"random",
				paste0("total_", cli_args$n),
				paste0("iter_", iter_idx, ".unbalanced-ols")
			)
		)
		# save unbalanced design with James-Stein results
		save_dge_data(
			comp_unbal_jse$object,
			comp_unbal_jse$results[, .SD, keyby = "target_id"],
			comp_unbal$design,
			file.path(
				RESULT_DIR,
				"Iterations",
				"random",
				paste0("total_", cli_args$n),
				paste0("iter_", iter_idx, ".unbalanced-jse")
			)
		)
	}
}

