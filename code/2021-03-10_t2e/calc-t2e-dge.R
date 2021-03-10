# ==============================================================================
# Meta
# ==============================================================================
# calc-t2e-dge
# ------------------------------------------------
# Author: James Hawley
# Description: Differential expression analysis for genes between TMPRSS2 and ERG genes


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))

loginfo("Loading packages")
suppressWarnings(library("data.table"))
suppressWarnings(library("sleuth"))
source(file.path("..", "jse-shrinkage", "jse.R"))

RESULT_DIR <- file.path("..", "..", "results", "2021-03-10_t2e")
KALLISTO_DIR <- file.path("..", "..", "data", "CPC-GENE")


# ==============================================================================
# Functions
# ==============================================================================
# helper function for verbose printing
vprint <- function(s, prefix = NULL) {
	loginfo(paste(prefix, s))
}

# Helper function for determining which is the single sample in an unbalanced design
which_single <- function(small_meta) {
	single_cdn <- small_meta[, .N, by = "condition"][N == 1, condition]
	single_sample <- small_meta[condition == single_cdn, sample]
	return(single_sample)
}

# Helper function to perform differential analysis
diff_ge <- function(small_meta) {
	so <- sleuth_prep(
		small_meta,
		extra_bootstrap_summary = TRUE,
		num_cores = 4
	)
	vprint("Fitting model and hypothesis testing", iter_prefix)
	so <- sleuth_fit(so, ~condition, "full")

	# perform differential analysis
	so <- sleuth_wt(so, "conditionYes")

	# extract results
	so_genes <- as.data.table(sleuth_results(
		so,
		"conditionYes",
		"wt",
		show_all = FALSE
	))

	# return this set of items
	return(list(
		"object" = so,
		"results" = so_genes,
		"design" = small_meta
	))
}

# Helper function for performing the James-Stein shrinkage
perform_jse <- function(so, small_meta, single_sample, subset_tx = NULL) {
	so_jse <- jse_shrinkage(
		so,
		which_sample = single_sample,
		which_beta = "conditionYes",
		# only select the desired transcripts
		s_targets = subset_tx
	)
	so_jse_results <- jse_wald_test(so_jse)

	# return this set of items
	return(list(
		"object" = so_jse,
		"results" = so_jse_results,
		"design" = small_meta
	))
}

# Helper function for saving differential analysis data
save_dge_data <- function(obj, res, design, prefix) {
	# ensure that the directory exists
	dest_dir <- dirname(prefix)
	if (!dir.exists(dest_dir)) {
		dir.create(dest_dir, recursive = TRUE)
	}

	# save sleuth object
	saveRDS(
		obj,
		paste0(prefix, ".sleuth-object.rds")
	)
	# save differential analysis results
	fwrite(
		res,
		paste0(prefix, ".genes.tsv"),
		sep = "\t",
		col.names = TRUE
	)
	# save experimental design for troubleshooting
	fwrite(
		design,
		paste0(prefix, ".design.tsv"),
		sep = "\t",
		col.names = TRUE
	)
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# convert into a design matrix
design <- meta[, .(
	sample = Label,
	condition = T2E_Status, # Yes = T2E+
	path = file.path(KALLISTO_DIR, Label)
)]

# gene annotation
hg38 <- fread(
	file.path("..", "..", "data", "Kallisto_GRCh38_Index", "Homo_sapiens.GRCh38.96.gtf"),
	skip = 5,
	sep = "\t",
	header = FALSE,
	col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
)

# convert ;-separated values into their own columns for certain attributes
hg38_genes <- hg38[feature == "gene"]
hg38_genes_split <- hg38_genes[, tstrsplit(attribute, "; ")]
hg38_genes$gene_id <- hg38_genes_split[, tstrsplit(V1, "\"")][, V2]
hg38_genes$gene_name <- hg38_genes_split[, tstrsplit(V3, "\"")][, V2]
hg38_genes$gene_biotype <- hg38_genes_split[, tstrsplit(V5, "\"")][, V2]

# drop the following columns
hg38_genes[, `:=`(
	source = NULL,
	feature = NULL,
	score = NULL,
	frame = NULL,
	attribute = NULL
)]

hg38_tx <- hg38[feature == "transcript"]
hg38_tx_split <- hg38_tx[, tstrsplit(attribute, "; ")]
hg38_tx$gene_id <- hg38_tx_split[, tstrsplit(V1, "\"")][, V2]
hg38_tx$transcript_id <- hg38_tx_split[, tstrsplit(V3, "\"")][, V2]
hg38_tx$gene_name <- hg38_tx_split[, tstrsplit(V5, "\"")][, V2]
hg38_tx$gene_biotype <- hg38_tx_split[, tstrsplit(V7, "\"")][, V2]
hg38_tx$transcript_name <- hg38_tx_split[, tstrsplit(V8, "\"")][, V2]
hg38_tx$gene_biotype <- hg38_tx_split[, tstrsplit(V10, "\"")][, V2]

# drop the following columns
hg38_tx[, `:=`(
	source = NULL,
	feature = NULL,
	score = NULL,
	frame = NULL,
	attribute = NULL
)]

# only collect genes between ERG and TMPRSS2
# these two genes are also on TAD boundaries, so we don't need to extend past them
t2e_genes <- hg38_genes[(chr == 21) & (start >= 38380027) & (end <= 41531116)]
t2e_tx <- hg38_tx[(chr == 21) & (start >= 38380027) & (end <= 41531116)]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations full calculation")

full_dge <- diff_ge(design)
agg_res <- list(
	full_dge$results[
		target_id %in% t2e_tx$transcript_id,
		.(target_id, b, se_b, pval, qval, Set = "Full", Single_Sample = NA)
	]
)

# redo the same analyses with 1 T2E+ sample and 6 T2E- samples
for (i in 1:6) {
	vprint(paste0("Single T2E+ (iteration ", i, " / 6)"))
	t2e_idx <- design[, which(condition == "Yes")]
	single_sample_idx <- t2e_idx[i]
	single_sample <- design[single_sample_idx, Label]
	small_meta <- rbindlist(list(
		design[single_sample_idx],
		design[-t2e_idx]
	))
	small_dge <- diff_ge(small_meta)
	small_jse <- perform_jse(small_dge$object, small_meta, single_sample, t2e_tx)
	agg_res[[length(agg_res) + 1]] <- small_dge$results[
		target_id %in% t2e_tx$transcript_id,
		.(target_id, b, se_b, pval, qval, Set = "Single_T2E", Single_Sample = single_sample)
	]
}

# redo the same analyses with 6 T2E+ sample and 1 T2E- samples
for (i in 1:6) {
	vprint(paste0("Single T2E- (iteration ", i, " / 6)"))
	t2e_idx <- design[, which(condition == "No")]
	single_sample_idx <- t2e_idx[i]
	single_sample <- design[single_sample_idx, Label]
	small_meta <- rbindlist(list(
		design[single_sample_idx],
		design[-t2e_idx]
	))
	small_dge <- diff_ge(small_meta)
	agg_res[[length(agg_res) + 1]] <- small_dge$results[
		target_id %in% t2e_tx$transcript_id,
		.(target_id, b, se_b, pval, qval, Set = "Single_nonT2E", Single_Sample = single_sample)
	]
}

if (!dir.exists(RESULT_DIR)) {
	dir.create(RESULT_DIR, recursive = TRUE)
}
saveRDS(agg_res, file.path(RESULT_DIR, "agg-res.rds"))
