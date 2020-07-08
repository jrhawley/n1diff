# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Perform differential expression analysis with only 1 sample in the treatment group"
    )
    PARSER$add_argument(
        "sampleID",
        type = "character",
        help = "Single SampleID to compare against other group"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        sampleID = "PCa13266"
    )
}

KALLISTO_DIR <- file.path("..", "..", "Data", "CPC-GENE")

# ==============================================================================
# Functions
# ==============================================================================
vprint <- function(msg) {
    underline <- paste0(rep("-", nchar(msg)), collapse = "")
    cat(msg, "\n", underline, "\n", sep = "")
}

remove_terms <- function(form, term) {
  fterms <- terms(form)
  fac <- attr(fterms, "factors")
  lab <- attr(fterms, "term.labels")
  idx <- which(as.logical(fac[term, ]))
  if (length(lab[-idx]) == 0) {
    return(~1)
  }
  new_fterms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
  return(formula(new_fterms))
}

me_model_by_row <- function(obj, design, bs_summary, which_var = 'obs_counts', ignore_col = c()) {
  which_var <- match.arg(which_var, c('obs_counts', 'obs_tpm'))
  if (which_var == "obs_counts")
    sigma_var <- "sigma_q_sq"
  else
    sigma_var <- "sigma_q_sq_tpm"

  stopifnot( all.equal(names(bs_summary[[sigma_var]]), rownames(bs_summary[[which_var]])) )
  stopifnot( length(bs_summary[[sigma_var]]) == nrow(bs_summary[[which_var]]))

  models <- lapply(1:nrow(bs_summary[[which_var]]),
    function(i) {
      sleuth:::me_model(design, bs_summary[[which_var]][i, -ignore_col], bs_summary[[sigma_var]][i])
    })
  names(models) <- rownames(bs_summary[[which_var]])

  models
}

jse_fit <- function(obj, formula = NULL, fit_name = 'jse', s_targets, which_beta, shrinkage_coef = NULL,...) {
  stopifnot( is(obj, 'sleuth') )
  stopifnot( sleuth:::check_norm_status(obj) )

  extra_opts <- list(...)
  if ('which_var' %in% names(extra_opts)) {
    which_var <- extra_opts$which_var
  } else {
    which_var <- 'obs_counts'
  }
  if ('n_bins' %in% names(extra_opts)) {
    n_bins <- extra_opts$n_bins
  } else {
    n_bins <- 100
  }
  if ('lwr' %in% names(extra_opts)) {
    lwr <- extra_opts$lwr
  } else {
    lwr <- 0.25
  }
  if ('upr' %in% names(extra_opts)) {
    upr <- extra_opts$lwr
  } else {
    upr <- 0.75
  }

  which_var <- match.arg(which_var, c('obs_counts', 'obs_tpm'))

  if (is.null(obj$bs_summary[[which_var]])) {
    if (which_var == "obs_tpm") {
      stop(which_var, " does not exist. Make sure sleuth_prep was used with 'read_bootstrap_tpm'",
           " set to TRUE")
    } else {
      stop(which_var, " does not exist. Make sure sleuth_prep was used with 'extra_bootstrap_summary'",
           " set to TRUE")
    }
  }

  if ( is.null(formula) ) {
    formula <- obj$full_formula
    if (is.null(formula)) {
      stop("'formula' was not specified and the 'full' model was not specified in `sleuth_prep`.",
        " Please specify a formula and a label.")
    }
  } else if ( !is(formula, 'formula') && !is(formula, 'matrix') ) {
    stop("'", substitute(formula), "' is not a valid 'formula' or 'matrix'")
  }

  if ( is.null(fit_name) ) {
    fit_name <- 'jse'
  } else if ( !is(fit_name, 'character') ) {
    stop("'", substitute(fit_name), "' is not a valid 'character'")
  }

  if ( length(fit_name) > 1 ) {
    stop("'", substitute(fit_name), "' is of length greater than one.",
      " Please only supply one string.")
  }

  # stop if dimension of estimate is too small
  if ( length(s_targets) <= 2 ) {
      stop('Selected 2 or fewer targets, for which the James-Stein estimator provides no benefit over the standard model')
  }

  # stop if shrinkage_coef is negative
  if (!is.null(shrinkage_coef) && shrinkage_coef < 0) {
    stop('`shrinkage_coef` must be non-negative.')
  }

  # TODO: check if model matrix is full rank
  X <- NULL
  if ( is(formula, 'formula') ) {
    X <- model.matrix(formula, obj$sample_to_covariates)
  } else {
    if ( is.null(colnames(formula)) ) {
      stop("If matrix is supplied, column names must also be supplied.")
    }
    X <- formula
  }
  rownames(X) <- obj$sample_to_covariates$sample
  A <- solve(t(X) %*% X)

  # remove covariate(s) related to single sample
  if ( is(formula, 'formula') ) {
    # identify single sample to be calculated separately
    which_label <- which(attr(terms(formula), "term.labels") == which_beta)
    design_table <- table(X[, which_label + 1]) # +1 because of the (Intercept) in the first column
    d_val <- as.integer(names(design_table[design_table == 1]))
    single_sample <- names(which(X[, which_label + 1] == d_val))
    if (length(single_sample) > 1) {
      stop("`which_beta` is not a factor that contains a single sample in one group. Select a different covariate, or rerun `sleuth_prep` with a new design matrix.")
    }

    formula_null <- remove_terms(formula, which_beta)
    other_sample_to_covariates <- obj$sample_to_covariates[obj$sample_to_covariates$sample != single_sample, ]
    X_null <- model.matrix(formula_null, other_sample_to_covariates)
    rownames(X_null) <- other_sample_to_covariates$sample
  } else {
    if ( is.null(colnames(formula)) ) {
      stop("If matrix is supplied, column names must also be supplied.")
    }
    X_null <- X[, -which_beta]
  }
  A_null <- solve(t(X_null) %*% X_null)

  sleuth:::msg("fitting measurement error models")
  which_col <- which(colnames(obj$bs_summary[[which_var]]) == single_sample)
  mes <- me_model_by_row(obj, X_null, obj$bs_summary, which_var, ignore_col = which_col)
  # mes <- sleuth:::me_model_by_row(obj, X_null, obj$bs_summary[[which_var]][, -which_col], which_var)
  tid <- names(mes)

  mes_df <- dplyr::bind_rows(lapply(mes,
    function(x) {
      data.frame(rss = x$rss, sigma_sq = x$sigma_sq, sigma_q_sq = x$sigma_q_sq,
        mean_obs = x$mean_obs, var_obs = x$var_obs)
    }))

  mes_df$target_id <- tid
  rm(tid)

  mes_df <- dplyr::mutate(mes_df, sigma_sq_pmax = pmax(sigma_sq, 0))

  # FIXME: sometimes when sigma is negative the shrinkage estimation becomes NA
  # this is for the few set of transcripts, but should be able to just do some
  # simple fix
  sleuth:::msg('shrinkage estimation')
  swg <- sleuth:::sliding_window_grouping(mes_df, 'mean_obs', 'sigma_sq_pmax',
    n_bins = n_bins, lwr = lwr, upr = upr, ignore_zeroes = TRUE)
  l_smooth <- sleuth:::shrink_df(swg, sqrt(sqrt(sigma_sq_pmax)) ~ mean_obs, 'iqr')
  l_smooth <- dplyr::select(
    dplyr::mutate(l_smooth, smooth_sigma_sq = shrink ^ 4),
    -shrink)

  l_smooth <- dplyr::mutate(l_smooth,
    smooth_sigma_sq_pmax = pmax(smooth_sigma_sq, sigma_sq))

  sleuth:::msg('computing variance of betas')
  beta_covars <- lapply(1:nrow(l_smooth),
    function(i) {
      row <- l_smooth[i, ]
      with(
        row,
        sleuth:::covar_beta(smooth_sigma_sq_pmax + sigma_q_sq, X_null, A_null)
      )
    })
  names(beta_covars) <- l_smooth$target_id
  
  # only keep targets within s_targets
  mes <- mes[names(mes) %in% s_targets]
  l_smooth <- l_smooth[l_smooth$target_id %in% s_targets, ]
  beta_covars <- beta_covars[names(beta_covars) %in% s_targets]

  # calculate James-Stein estimator

  # get intercept coefficients
  beta_0 <- sapply(
    names(mes),
    function(tid) {
      mes[[tid]]$ols_fit$coefficients["(Intercept)"] # you need to generalize this step a bit more, include the entire null model not just the intercept
    }
  )
  names(beta_0) <- names(mes)

  # get observed counts
  delta <- obj$bs_summary[[which_var]][rownames(obj$bs_summary[[which_var]]) %in% s_targets, which_col]
  # get covariance matrix
  sigma <- diag(beta_covars)
  sigma_inv <- solve(sigma)
  trace_sigma <- sum(unlist(beta_covars))
  lambda_L <- max(unlist(beta_covars))
  
  # check that conditions for JSE are met (i.e. Trace(Sigma) <= 2 * \lambda_L )
  if (trace_sigma <= 2 * lambda_L) {
    stop('Conditions for James-Stein criteria are not met: Tr(Sigma) > 2 * lambda_L')
  }

  # calculate normalized, observed value \nu
  nu <- as.vector(sqrt(sigma_inv) %*% (delta - beta_0))
  names(nu) <- names(mes)

  # maximize shrinkage coefficient
  if (is.null(shrinkage_coef)) {
    shrinkage_coef <- max(0, 2 * (trace_sigma / lambda_L - 2))
  } else if (shrinkage_coef > 2 * (trace_sigma / lambda_L - 2)) {
    warning('`shrinkage_coef` > 2 * (Tr(Sigma) / lambda_L - 2). Capping at max value')
    shrinkage_coef <- 2 * (trace_sigma / lambda_L - 2)
  }

  b1 <- (1 - shrinkage_coef / sum(nu^2)) * nu

  # swap signs on the coefficient if d_val is 0 (this keeps the "mutated" samples with coefficient 1, non-mutated with 0)
  if (d_val == 0) {
    b1 <- -b1
  }

  if ( is.null(obj$fits) ) {
    obj$fits <- list()
  }

  obj$fits[[fit_name]] <- list(
    models = mes,
    summary = l_smooth,
    beta_covars = beta_covars,
    formula = formula,
    design_matrix = X,
    transform_synced = TRUE,
    which_var = which_var,
    delta = delta,
    sigma = sigma,
    trace_sigma = trace_sigma,
    lambda_L = lambda_L,
    shrinkage_coef = shrinkage_coef,
    b0 = beta_0,
    b1 = b1,
    nu = nu
  )

  class(obj$fits[[fit_name]]) <- 'sleuth_model'

  obj
}

jse_shrinkage <- function(obj, s_targets, which_sample, which_beta, which_fit = "full", shrinkage_coef = NULL, which_var = "obs_counts") {
  # only keep targets within s_targets
  fit <- so$fits[[which_fit]]
  mes <- fit$models[names(fit$models) %in% s_targets]
  l_smooth <- fit$summary[fit$summary$target_id %in% s_targets, ]
  summary <- fit$summary[fit$summary$target_id %in% s_targets, ]
  covariances <- summary$smooth_sigma_sq_pmax + summary$sigma_q_sq
  names(covariances) <- summary$target_id

  # get null design that does not include this covariate
  X <- fit$design_matrix
  X_null <- as.matrix(X[which(rownames(X) == which_sample), -which(colnames(X) == which_beta)])

  # get intercept coefficients
  mu_0 <- sapply(
    s_targets,
    function(tid) {
      other_covariate_coefs <- mes[[tid]]$ols_fit$coefficients[-which(names(mes[[tid]]$ols_fit$coefficients) == which_beta)]
      return(as.numeric(X_null %*% other_covariate_coefs))
    }
  )
  names(mu_0) <- names(mes)
  
  naive <- sapply(
    s_targets,
    function(tid) {
      mes[[tid]]$ols_fit$coefficients[names(mes[[tid]]$ols_fit$coefficients) == which_beta]
    }
  )
  names(naive) <- names(mes)

  # get observed counts
  which_col <- which(colnames(obj$bs_summary[[which_var]]) == which_sample)
  delta <- obj$bs_summary[[which_var]][rownames(obj$bs_summary[[which_var]]) %in% s_targets, which_col]
  # get covariance matrix
  sigma <- diag(covariances)
  sigma_inv <- solve(sigma)
  trace_sigma <- sum(covariances)
  lambda_L <- max(covariances)
  
  # check that conditions for JSE are met (i.e. Trace(Sigma) <= 2 * \lambda_L )
  if (trace_sigma <= 2 * lambda_L) {
    stop('Conditions for James-Stein criteria are not met: Tr(Sigma) > 2 * lambda_L')
  }

  # calculate normalized, observed value \nu
  nu <- as.vector(sqrt(sigma_inv) %*% (delta - mu_0))
  names(nu) <- names(mes)

  # maximize shrinkage coefficient
  if (is.null(shrinkage_coef)) {
    shrinkage_coef <- max(0, 2 * (trace_sigma / lambda_L - 2))
  } else if (shrinkage_coef > 2 * (trace_sigma / lambda_L - 2)) {
    warning('`shrinkage_coef` > 2 * (Tr(Sigma) / lambda_L - 2). Capping at max value')
    shrinkage_coef <- 2 * (trace_sigma / lambda_L - 2)
  }

  b1 <- (1 - shrinkage_coef / sum(nu^2)) * (delta - mu_0)

  # swap signs on the coefficient if coefficient in design is 0
  # (this keeps the "mutated" samples with coefficient 1, non-mutated with 0)
  if (X[rownames(X) == which_sample, which_beta] == 0) {
    b1 <- -b1
  }
  return(list(
    which_var = which_var,
    sigma = sigma,
    trace_sigma = trace_sigma,
    lambda_L = lambda_L,
    shrinkage_coef = shrinkage_coef,
    delta = delta,
    b0 = mu_0,
    nu = nu,
    b1 = b1,
    naive = naive
  ))
}

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    file.path(KALLISTO_DIR, "config.tsv"),
    sep = "\t",
    header = TRUE
)
SAMPLES <- metadata$SampleID

# add kallisto paths
metadata[, Path := file.path(KALLISTO_DIR, SampleID)]

design <- metadata[, .(sample = SampleID, condition = T2E_Status, path = Path)]

# ==============================================================================
# Analysis
# ==============================================================================
# load naive sleuth object
so <- readRDS(
  file.path("..", "2020-07-04_naive-analysis", "sleuth", paste0(ARGS$sampleID, ".sleuth-object.rds"))
)

# calculate James-Stein estimator for single mutated sample
s_targets <- unique(head(so$obs_norm_filt, 210)$target_id)
jse <- jse_shrinkage(so, s_targets = s_targets, which_sample = ARGS$sampleID, which_beta = "conditionYes")

# ==============================================================================
# Save results
# ==============================================================================
# save sleuth object
saveRDS(jse, paste0("sleuth/", ARGS$sampleID, ".jse-object.rds"))
