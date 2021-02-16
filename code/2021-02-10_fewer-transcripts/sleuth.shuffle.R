# ==============================================================================
# Meta
# ==============================================================================
# swi-snf-sleuth
# ------------------------------------------------
# Author: James Hawley
# Description: Differential expression analysis for SWI/SNF complex


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))
source(file.path("..", "jse-shrinkage", "jse.R"))

RESULT_DIR <- file.path("..", "..", "results", "2021-02-19_fewer-transcripts")

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


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# replace "wt" with "ctrl" to ensure that its the baseline control, alphabetically
meta[condition == "wt", condition := "ctrl"]

# SWI/SNF complex subunits
swi_snf <- fread("swi-snf-complex-subunits.tsv", sep = "\t", header = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

for (i in 1:(TOTAL_REPS / 2)) {
	loginfo("\tPreparing sleuth objects")
	# balanced samples
	sampled_metadata <- select_samples(meta, 6)
	for (j in 1:2) {
		iter_idx <- (i - 1) * 2 + j
		for (k in c("balanced", "unbalanced")) {
			iter_prefix <- paste0("(", iter_idx, " / ", TOTAL_REPS, ")")
			so <- sleuth_prep(
				sampled_metadata[[j]][[k]],
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

			if (k == "unbalanced") {
				vprint("Calculating James-Stein estimates", iter_prefix)
				# specifiy the sample with the single condition
				single_cdn <- sampled_metadata[[j]][[k]][, .N, by = "condition"][N == 1, condition]
				single_sample <- sampled_metadata[[j]][[k]][condition == single_cdn, sample]
				so_jse <- jse_shrinkage(
					so,
					which_sample = single_sample,
					which_beta = "conditionmu",
					# only select SWI/SNF complex subunits
					s_targets = swi_snf$transcript_id
				)
				saveRDS(
					so_jse,
					file.path(
						RESULT_DIR,
						"Iterations",
						paste0("Total_", cli_args$n),
						k,
						paste0(iter_idx, ".jse-object.rds")
					)
				)
				so_jse_results <- jse_wald_test(so_jse)
				fwrite(
					so_jse_results,
					file.path(
						RESULT_DIR,
						"Iterations",
						paste0("Total_", cli_args$n),
						k,
						paste0(iter_idx, ".genes.jse.tsv")
					),
					sep = "\t",
					col.names = TRUE
				)
			}

			vprint("Saving data", iter_prefix)
			fwrite(
				so_genes,
				file.path(
					RESULT_DIR,
					"Iterations",
					paste0("Total_", cli_args$n),
					k,
					paste0(iter_idx, ".genes.tsv")
				),
				sep = "\t",
				col.names = TRUE
			)
			saveRDS(
				so,
				file.path(
					RESULT_DIR,
					"Iterations",
					paste0("Total_", cli_args$n),
					k,
					paste0(iter_idx, ".sleuth-object.rds")
				)
			)
		}
	}
}
