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
