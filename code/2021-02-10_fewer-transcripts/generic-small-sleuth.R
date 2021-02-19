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
	parser$add_arg(
		"n",
		type = "integer",
		description = "Number of transcripts to sample"
	)
	cli_args <- parser$parse_args()
} else {
	cli_args <- list(
		n = 3
	)
)

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
sample_transcripts <- function(obj, n_tx = 3) {
	# only consider transcripts that have not been filtered out due to low read counts or other QC filters
	possible_tx <- names(which(obj$filter_bool))
	return(sample(possible_tx, n_tx, replace = FALSE))
}

# Helper function to perform differential analysis
diff_ge <- function(small_meta, unbalanced = FALSE, iter_idx = "", n_tx = 3) {
	iter_prefix <- paste0("(", iter_idx, " / ", TOTAL_REPS, ")")
	so <- sleuth_prep(
		small_meta,
		extra_bootstrap_summary = TRUE,
		num_cores = 8
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
	# pick a random subset of transcripts
	subset_tx <- sample_transcripts(so, n_tx)
	so_genes <- so_genes[transcript_id %in% subset_tx]

	# perform James-Stein shrinkage if desired
	if (unbalanced) {
		vprint("Calculating James-Stein estimates", iter_prefix)
		# specifiy the sample with the single condition
		single_sample <- which_single(small_meta)
		so_jse <- jse_shrinkage(
			so,
			which_sample = single_sample,
			which_beta = "conditionmu",
			# only select SWI/SNF complex subunits
			s_targets = subset_tx
		)
		so_jse_results <- jse_wald_test(so_jse)
	} else {
		so_jse <- NA
		so_jse_results <- NA
	}

	# return this set of items for future saving
	return(list(
		"object" = so,
		"results" = so_genes,
		"design" = small_meta,
		"shrunk_object" = so_jse,
		"shrunk_results" = so_jse_results
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

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# replace "wt" with "ctrl" to ensure that its the baseline control, alphabetically
meta[condition == "wt", condition := "ctrl"]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# set random seed for reproducible results
set.seed(10002)

for (i in 1:(TOTAL_REPS / 2)) {
	loginfo("\tPreparing sleuth objects")
	# balanced samples
	sampled_metadata <- select_samples(meta, 6)
	for (j in 1:2) {
		iter_idx <- (i - 1) * 2 + j
		# perform differential analysis of balanced experimental design
		comp_bal <- diff_ge(
			sampled_metadata[[j]][["balanced"]],
			FALSE,
			iter_idx,
			cli_args$n
		)

		# perform differential analysis of unbalanced experimental design
		# both OLS and James-Stein estimates
		comp_unbal <- diff_ge(
			sampled_metadata[[j]][["unbalanced"]],
			TRUE,
			iter_idx,
			cli_args$n
		)

		# check that the trace-lambda condition for James-Stein shrinkage is met
		while (is.na(comp_unbal$shrunk_object$shrinkage_coef)) {
			loginfo("Failed shrinkage, trying with new samples")
			# if not, re-select samples to work with
			sampled_metadata[[j]][["unbalanced"]] <- select_samples(meta, 6)[[j]][["unbalanced"]]
			comp_unbal <- diff_ge(
				sampled_metadata[[j]][["unbalanced"]],
				TRUE,
				iter_idx,
				cli_args$n
			)
		}

		vprint("Saving data", iter_idx)
		# save balanced design
		save_dge_data(
			comp_bal$object,
			comp_bal$results,
			comp_bal$design,
			file.path(RESULT_DIR, "Iterations", "random", paste0("total_", cli_args$n), "balanced", iter_idx)
		)
		# save unbalanced design with OLS results
		save_dge_data(
			comp_unbal$object,
			comp_unbal$results,
			comp_unbal$design,
			file.path(RESULT_DIR, "Iterations", "random", paste0("total_", cli_args$n), "unbalanced", paste0(iter_idx, ".ols"))
		)
		# save unbalanced design with James-Stein results
		save_dge_data(
			comp_unbal$shrunk_object,
			comp_unbal$shrunk_results,
			comp_unbal$design,
			file.path(RESULT_DIR, "Iterations", "random", paste0("total_", cli_args$n), "unbalanced", paste0(iter_idx, ".jse"))
		)
	}
}
