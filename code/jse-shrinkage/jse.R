# ==============================================================================
# Meta
# ==============================================================================
# jse
# ------------------------------------------------
# Author: James Hawley
# Description: Functions for calculating James-Stein fold-change estimates and hypothesis testing


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Functions
# ==============================================================================
# Perform James-Stein shrinkage
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

# Perform James-Stein hypothesis testing with the Wald test
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