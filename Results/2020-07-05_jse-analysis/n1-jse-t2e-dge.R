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

sleuth_jse <- function(obj, s_targets, which_beta, which_model = 'reduced') {
  stopifnot( is(obj, 'sleuth') )

  if ( length(s_targets) <= 2 ) {
      stop('Selected 2 or fewer targets, for which the James-Stein estimator provides no benefit over the standard model')
  }

  if ( !sleuth:::model_exists(obj, which_model) ) {
    stop("'", sleuth:::which_model, "' is not a valid model. Please see models(",
      substitute(obj), ") for a list of fitted models")
  }

  if(!obj$fits[[which_model]]$transform_synced) {
    stop("Model '", sleuth:::which_model, "' was not computed using the sleuth object's",
         " current transform function. Please rerun sleuth_fit for this model.")
  }

  d_matrix <- obj$fits[[which_model]]$design_matrix

  # identify single sample to be calculate separately
  design_table <- table(d_matrix[, which_beta])
  d_val <- as.integer(names(design_table[design_table == 1]))
  single_sample <- names(which(d_matrix[, which_beta] == d_val))

  # get the beta index
  beta_i <- which(colnames(d_matrix) == "(Intercept)")

  if ( length(beta_i) == 0 ) {
    stop(paste0("'", which_beta,
        "' doesn't appear in your design. Try one of the following:\n",
        paste(colnames(d_matrix), collapse = ' ')))
  } else if ( length(beta_i) > 1 ) {
    stop(paste0("Sorry. '", which_beta, "' is ambiguous for columns: ",
        paste(colnames(d_matrix[beta_i]), collapse = ' ')))
  }

  b <- sapply(
    obj$fits[[ which_model ]]$models,
    function(x) {
      x$ols_fit$coefficients[ beta_i ]
    })
  names(b) <- names(obj$fits[[ which_model ]]$models)

  res <- obj$fits[[ which_model ]]$summary
  res$target_id <- as.character(res$target_id)
  res <- res[match(names(b), res$target_id), ]

  stopifnot( all.equal(res$target_id, names(b)) )

  se <- sapply(obj$fits[[ which_model ]]$beta_covars,
    function(x) {
      x[beta_i, beta_i]
    })
  se <- sqrt( se )
  se <- se[ names(b) ]

  stopifnot( all.equal(names(b), names(se)) )

  res <- dplyr::mutate(
    res,
    b_0 = b,
    se_b_0 = se
  )

  # only keep values in s_targets
  res <- res[res[, "target_id"] %in% s_targets, ]
  
  # get count estimates for single sample
  these_counts <- obj$obs_norm_filt[obj$obs_norm_filt[, "sample"] == single_sample, ]
  these_counts <- these_counts[these_counts[, "target_id"] %in% s_targets, ]
  these_counts
  
#   obj <- add_test(obj, res, which_beta, 'jse', which_model)

  #obj
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
# only keep specified sample and all samples of the other T2E_Status group
status_of_specified <- design[sample == ARGS$sampleID, condition]
design <- design[(sample == ARGS$sampleID) | (condition != status_of_specified), .SD]

# create sleuth object
so <- sleuth_prep(
    design,
    extra_bootstrap_summary = TRUE,
    num_cores = 2
)

# fit model to estimate \beta_{0,t} and variances
so <- sleuth_fit(so, ~condition, "full")
so <- sleuth_fit(so, ~1, "reduced")

# calculate James-Stein estimator for single mutated sample
s_targets <- unique(head(so$obs_norm_filt, 210)$target_id)
jse <- sleuth_jse(so, s_targets, "conditionYes")
# perform differential analysis
so <- sleuth_wt(so, "conditionYes")

# extract results
so_transcripts <- as.data.table(sleuth_results(so, "conditionYes", "wt", show_all = FALSE, pval_aggregate = FALSE))

# ==============================================================================
# Save results
# ==============================================================================
# save sleuth object
saveRDS(so, paste0("sleuth/", ARGS$sampleID, ".sleuth-object.rds"))

# save annotated table
fwrite(
    so_transcripts,
    paste0("sleuth/", ARGS$sampleID, ".results.tsv"),
    sep = "\t",
    col.names = TRUE
)
