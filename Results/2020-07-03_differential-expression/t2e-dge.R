# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))

KALLISTO_DIR <- file.path("..", "..", "Data", "CPC-GENE")

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
# create sleuth object, aggregating transcripts to genes
so <- sleuth_prep(
    design,
    extra_bootstrap_summary = TRUE,
    num_cores = 2
)

# fit full model
so <- sleuth_fit(so, ~condition, "full")

# perform differential analysis
so <- sleuth_wt(so, "conditionYes")

# extract results
so_transcripts <- as.data.table(sleuth_results(so, "conditionYes", "wt", show_all = FALSE, pval_aggregate = FALSE))

# ==============================================================================
# Save results
# ==============================================================================
# save sleuth object
saveRDS(so, "sleuth-object.rds")

# write full kallisto count table
full_table <- as.data.table(kallisto_table(so))
fwrite(
    full_table,
    "abundance.all.tsv",
    sep = "\t",
    col.names = TRUE
)
# save annotated table
fwrite(so_transcripts, "results.transcripts.tsv", sep = "\t", col.names = TRUE)
