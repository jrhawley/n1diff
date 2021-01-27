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
            male_ids <- metadata[condition == "male", sample(sample, half_total)]
            female_ids <- metadata[condition == "female", sample(sample, half_total)]
            return(rbindlist(list(
                metadata[sample %in% male_ids],
                metadata[sample %in% female_ids]
            )))
        }
    )
    # pick unbalanced samples 1 male vs total - 1 females
    usm_male_ids <- metadata[condition == "male", sample(sample, 1)]
    usm_female_ids <- metadata[condition == "female", sample(sample, total - 1)]
    unbalanced_single_male <- rbindlist(list(
        metadata[sample %in% usm_male_ids],
        metadata[sample %in% usm_female_ids]
    ))
    # pick unbalanced samples 1 female vs total - 1 males
    usf_male_ids <- metadata[condition == "male", sample(sample, 1)]
    usf_female_ids <- metadata[condition == "female", sample(sample, total - 1)]
    unbalanced_single_female <- rbindlist(list(
        metadata[sample %in% usf_male_ids],
        metadata[sample %in% usf_female_ids]
    ))

    # aggregate and return entire list
    return(list(
        list(
            "balanced" = balanced[[1]],
            "unbalanced" = unbalanced_single_female
        ),
        list(
            "balanced" = balanced[[2]],
            "unbalanced" = unbalanced_single_male
        )
    ))
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# load transcript ID to gene ID mapping
t2g_map <- fread(
    "../../Data/Kallisto_Index/transcripts_to_genes.txt",
    sep = "\t",
    header = FALSE,
    col.names = c("target_id", "ens_gene", "ext_gene")
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# randomly select groups 10 times before performing comparison
for (i in 1:5) {
    loginfo("\tPreparing sleuth objects")
    # balanced samples
    sampled_metadata <- select_samples(meta, cli_args$n)
    for (j in 1:2) {
        iter_idx <- (i - 1) * 2 + j
        for (k in c("balanced", "unbalanced")) {
            so <- sleuth_prep(
                sampled_metadata[[j]][[k]],
                extra_bootstrap_summary = TRUE,
                num_cores = 8,
                target_mapping = t2g_map,
                aggregation_column = "ens_gene"
            )
            loginfo("\tFitting model and hypothesis testing")
            so <- sleuth_fit(so, ~condition, "full")

            # perform differential analysis
            so <- sleuth_wt(so, "conditionmale")

            # extract results
            so_transcripts <- as.data.table(sleuth_results(
                so,
                "conditionmale",
                "wt",
                show_all = FALSE,
                pval_aggregate = FALSE
            ))
            so_genes <- as.data.table(sleuth_results(
                so,
                "conditionmale",
                "wt",
                show_all = FALSE,
                pval_aggregate = TRUE
            ))

            loginfo("\tSaving data")
            fwrite(
                so_transcripts,
                file.path(
                    "Iterations",
                    paste0("Total_", cli_args$n),
                    k,
                    paste0(iter_idx, ".transcripts.tsv")
                ),
                sep = "\t",
                col.names = TRUE
            )
            fwrite(
                so_genes,
                file.path(
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
                    "Iterations",
                    paste0("Total_", cli_args$n),
                    k,
                    paste0(iter_idx, ".sleuth-object.rds")
                )
            )
        }
    }
}
