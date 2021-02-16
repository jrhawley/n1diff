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
suppressMessages(library("argparse"))

if (!interactive()) {
	parser <- ArgumentParser(description = "Differential expression analysis with some samples")
	parser$add_argument(
		"n",
		type = "integer",
		help = "Total number of samples to include in the comparison. Must be even"
	)
	cli_args <- parser$parse_args()
} else {
	cli_args <- list(
		"n" = 4
	)
}

if (cli_args$n %% 2 != 0) {
	stop("`n` argument must be even.")
}

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))

RESULT_DIR <- file.path("..", "..", "results", "small-sample-sizes")
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

jse_shrinkage <- function(obj, s_targets = NULL, which_sample, which_beta, which_fit = "full", shrinkage_coef = NULL, which_var = "obs_counts") {
	# shorthand for Not-In function
	'%ni%' <- Negate('%in%')

	# check input arguments

	# if no specific list of targets is given, take all targets, by default
	if (is.null(s_targets)) {
		s_targets <- names(obj$fits[[which_fit]]$beta_covars)
	}

	# only keep targets within s_targets
	fit <- obj$fits[[which_fit]]
	mes <- fit$models[names(fit$models) %in% s_targets]
	summary <- fit$summary[fit$summary$target_id %in% s_targets, ]
	covariances <- summary$smooth_sigma_sq_pmax + summary$sigma_q_sq
	names(covariances) <- summary$target_id
	rm(summary)

	# get OLS estimate (if any value is NULL, will return a list, otherwise a numeric vector)
	ols <- sapply(
		s_targets,
		function(tid) {
			mes[[tid]]$ols_fit$coefficients[[which_beta]]
		}
	)

	# check if any coefficients are NULL, and remove them if so
	null_counter <- sapply(ols, is.null)
	tx_to_remove <- c()
	if (any(null_counter)) {
		tx_to_remove <- names(ols)[null_counter]
		s_targets <- setdiff(s_targets, tx_to_remove)
		mes <- mes[names(mes) %ni% tx_to_remove]
		covariances <- covariances[names(covariances) %ni% tx_to_remove]
		ols <- as.vector(ols[names(ols) %ni% tx_to_remove], mode = "numeric")
	}

	# get covariance matrix
	sigma <- diag(covariances)
	sigma_inv <- solve(sigma)
	trace_sigma <- sum(covariances)
	lambda_L <- max(covariances)

	# check that conditions for JSE are met (i.e. Trace(Sigma) <= 2 * \lambda_L )
	if (trace_sigma <= 2 * lambda_L) {
		stop('Conditions for James-Stein criteria are not met: Tr(Sigma) > 2 * lambda_L')
	}

	# matrix multiplication creates a 1x1 matrix, not a scalar
	denom <- as.numeric(t(ols) %*% sigma_inv %*% ols)

	# maximize shrinkage coefficient if not given
	if (is.null(shrinkage_coef)) {
		shrinkage_numerator <- max(0, 2 * (trace_sigma / lambda_L - 2))
	} else if (shrinkage_numerator > 2 * (trace_sigma / lambda_L - 2)) {
		warning('`shrinkage_numerator` > 2 * (Tr(Sigma) / lambda_L - 2). Capping at max value')
		shrinkage_numerator <- 2 * (trace_sigma / lambda_L - 2)
	}
	shrinkage_coef <- (1 - shrinkage_numerator / denom)
	b1 <- shrinkage_coef * ols
	b1_var <- (
		shrinkage_coef ^ 2 * sigma
		+ ( shrinkage_coef ^ 2 + 2 * shrinkage_numerator / (denom ^ 2) - 1 ) * (b1 %*% t(b1))
	)

	return(list(
		targets = s_targets,
		na_targets = tx_to_remove,
		which_var = which_var,
		sigma = sigma,
		trace_sigma = trace_sigma,
		lambda_L = lambda_L,
		shrinkage_numerator = shrinkage_numerator,
		shrinkage_denom = denom,
		shrinkage_coef = shrinkage_coef,
		ols = ols,
		jse = b1,
		jse_var = b1_var
	))
}

jse_wald_test <- function(ols_list) {
	# shorthand for estimators
	targets <- ols_list$targets
	b1_ols <- ols_list$ols
	b1_jse <- ols_list$jse
	v <- ols_list$jse_var

	# calculate the Wald test statistic for each target
	w <- sapply(
		1:length(b1_jse),
		function(i) {
			b1_jse[i] ^ 2 / v[i, i]
		}
	)

	# calculate p-values
	pvals <- pchisq(
		q = w,
		df = 1,
		ncp = 0,
		# small p-value when w >> 0
		lower.tail = FALSE
	)
	# calculate multiple test correction
	qvals <- p.adjust(pvals, method = "fdr")

	# return data.table with the results
	return(data.table(
		target_id = targets,
		W = w,
		pval = pvals,
		qval = qvals,
		b = b1_jse,
		se_b = sapply(1:length(b1_jse), function(i) sqrt(v[i, i]))
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

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")
# randomly select groups this many times before performing comparison
total_reps <- 30

for (i in 1:(total_reps / 2)) {
	loginfo("\tPreparing sleuth objects")
	# balanced samples
	sampled_metadata <- select_samples(meta, cli_args$n)
	for (j in 1:2) {
		iter_idx <- (i - 1) * 2 + j
		for (k in c("balanced", "unbalanced")) {
			iter_prefix <- paste0("(", iter_idx, " / ", total_reps, ")")
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
					which_beta = "conditionmu"
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
