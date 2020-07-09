# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

KALLISTO_DIR <- file.path("..", "..", "Data", "CPC-GENE")

# ==============================================================================
# Data
# ==============================================================================
metadata <- fread(
    file.path(KALLISTO_DIR, "config.tsv"),
    sep = "\t",
    header = TRUE
)
SAMPLES <- metadata$SampleID

jse_dt <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        jse <- readRDS(file.path("sleuth", paste0(s, ".jse-object.rds")))
        n_valid <- length(jse$targets)
        n_invalid <- length(jse$na_targets)
        dt <- data.table(
            SampleID = s,
            target_id = rep(c(jse$targets, jse$na_targets), 2),
            beta_1 = c(
                jse$b1,
                rep(NA, n_invalid),
                jse$ols,
                rep(NA, n_invalid)
            ),
            method = rep(c("JSE", "OLS"), each = n_valid + n_invalid)
        )
        return(dt)
    }
))

shrinkage_dt <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        jse <- readRDS(file.path("sleuth", paste0(s, ".jse-object.rds")))
        dt <- data.table(
            SampleID = s,
            shrinkage = jse$shrinkage_coef
        )
        return(dt)
    }
))

s_targets <- jse_dt[, unique(target_id)]

full_obj <- readRDS(file.path("..", "2020-07-03_differential-expression", "sleuth-object.rds"))
full_dt <- data.table(
    SampleID = "All",
    target_id = s_targets,
    beta_1 = sapply(
        s_targets,
        function(tid) {
            v <- full_obj$fits$full$models[[tid]]$ols_fit$coefficients["conditionYes"]
            # if the transcript has been filtered out, this will return NULL, which needs to be handled properly
            return(ifelse(is.null(v), NA, v))
        }
    ),
    method = "Full"
)

# ==============================================================================
# Analysis
# ==============================================================================
all_dt <- rbindlist(list(jse_dt, full_dt))
# replace any NULLs with NAs
all_dt[is.null(beta_1), "beta_1"] <- NA

# ==============================================================================
# Save data
# ==============================================================================
fwrite(all_dt, "combined-results.tsv", sep = "\t", col.names = TRUE)
fwrite(shrinkage_dt, "shrinkage-coefficients.tsv", sep = "\t", col.names = TRUE)
