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

# ==============================================================================
# Data
# ==============================================================================
vprint("Loading metadata")
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
vprint("Normalizing counts")
# only keep specified sample and all samples of the other T2E_Status group
status_of_specified <- design[sample == ARGS$sampleID, condition]
design <- design[(sample == ARGS$sampleID) | (condition != status_of_specified), .SD]

# create sleuth object, aggregating transcripts to genes
so <- sleuth_prep(
    design,
    extra_bootstrap_summary = TRUE,
    num_cores = 2
)

vprint("Fitting model and hypothesis testing")
# fit full model
so <- sleuth_fit(so, ~condition, "full")

# perform differential analysis
so <- sleuth_wt(so, "conditionYes")

# extract results
so_transcripts <- as.data.table(sleuth_results(so, "conditionYes", "wt", show_all = FALSE, pval_aggregate = FALSE))

# ==============================================================================
# Save results
# ==============================================================================
vprint("Saving results")
# save sleuth object
saveRDS(so, paste0("sleuth/", ARGS$sampleID, ".sleuth-object.rds"))

# save annotated table
fwrite(
    so_transcripts,
    paste0("sleuth/", ARGS$sampleID, ".results.tsv"),
    sep = "\t",
    col.names = TRUE
)
