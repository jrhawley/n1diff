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
    PARSER$add_argument(
      "-t", "--target",
      type = "character",
      help = "Transcript ENSEMBL ID to test",
      nargs = "+"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        sampleID = "PCa13848",
        target = c(
          "ENST00000409312.5",
          "ENST00000409948.1",
          "ENST00000433917.5",
          "ENST00000457259.5",
          "ENST00000477296.5",
          "ENST00000509587.1",
          "ENST00000526421.5",
          "ENST00000543215.5",
          "ENST00000613154.4",
          "ENST00000644286.1"
        )
    )
}

KALLISTO_DIR <- file.path("..", "..", "Data", "CPC-GENE")

# ==============================================================================
# Functions
# ==============================================================================
'%ni%' <- Negate('%in%')

jse_shrinkage <- function(obj, s_targets, which_sample, which_beta, which_fit = "full", shrinkage_coef = NULL, which_var = "obs_counts") {
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
        shrinkage_coef <- (1 - shrinkage_numerator / denom)
    } else if (shrinkage_numerator > 2 * (trace_sigma / lambda_L - 2)) {
        warning('`shrinkage_numerator` > 2 * (Tr(Sigma) / lambda_L - 2). Capping at max value')
        shrinkage_numerator <- 2 * (trace_sigma / lambda_L - 2)
    }

    b1 <- shrinkage_coef * ols

    return(list(
        targets = s_targets,
        na_targets = tx_to_remove,
        which_var = which_var,
        sigma = sigma,
        trace_sigma = trace_sigma,
        lambda_L = lambda_L,
        shrinkage_numerator = shrinkage_numerator,
        shrinkage_coef = shrinkage_coef,
        ols = ols,
        b1 = b1
    ))
}

mse <- function(jse, full_so, targets) {
    # get estimates from full experiment
    full_targets <- sapply(
        targets,
        function(target_id) as.numeric(full_so$fits$full$models[[target_id]]$b1)
    )
    return(list(
        "OLS" = sum((jse$ols - full_targets) ^ 2),
        "JSE" = sum((jse$b1 - full_targets) ^ 2)
    ))
}

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread("../../Data/CPC-GENE/config.tsv", sep = "\t")
T2E_SAMPLES <- metadata[T2E_Status == "Yes", SampleID]
NONT2E_SAMPLES <- metadata[T2E_Status == "No", SampleID]

# load naive sleuth object
cat("Loading single objects\n")
sobjs <- lapply(
    metadata[, SampleID],
    function(s) {
        readRDS(
            file.path("..", "2020-07-04_naive-analysis", "sleuth", paste0(ARGS$sampleID, ".sleuth-object.rds"))
        )
    }
)
names(sobjs) <- metadata[, SampleID]

# load full sleuth object
cat("Loading full object\n")
so_full <- readRDS("../2020-07-03_differential-expression/sleuth-object.rds")

target_mapping <- rownames(so_full$bs_quants[[1]]$est_counts)

# ==============================================================================
# Analysis
# ==============================================================================
cat("Simulations\n")
# number of simulations
M <- 1e2
successful_tests <- 0
unsuccessful_tests <- 0
# range of number of transcripts to collect
n_transcripts <- 20
results <- list()
# random seed
set.seed(100)

# removing T2E+ sample
while (successful_tests < M) {
    single_t2e_sample <- sample(T2E_SAMPLES, 1)
    targets <- sample(target_mapping, n_transcripts)
    jse <- tryCatch(
        {
            jse_shrinkage(
                sobjs[[single_t2e_sample]],
                s_targets = targets,
                which_sample = single_t2e_sample,
                which_beta = "conditionYes"
            )
        },
        error = function(cond) {
            return(NA)
        }
    )
    if (is.na(jse)) {
        unsuccessful_tests <- unsuccessful_tests + 1
        if (unsuccessful_tests > 10 * M) break
        next
    }
    err <- mse(jse, so_full, targets)
    successful_tests <- successful_tests + 1
    results[[successful_tests]] <- data.table(Index = successful_tests, JSE = err$JSE, OLS = err$OLS)
    cat(successful_tests, "\n")
}

res_dt <- rbindlist(results)
res_dt[, Diff := OLS - JSE]

fwrite(res_dt, paste0("simulation.removing-T2E.", n_transcripts, ".tsv"), sep = "\t")
