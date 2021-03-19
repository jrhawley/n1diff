# ==============================================================================
# Meta
# ==============================================================================
# plot-volcanoes
# ------------------------------------------------
# Author: James Hawley
# Description: Plot differential gene expression results


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))

RESULT_DIR <- file.path("..", "..", "results", "small-sample-sizes")
PLOT_DIR <- file.path(RESULT_DIR, "Plots")

if (!dir.exists(PLOT_DIR)) {
    dir.create(PLOT_DIR, recursive = TRUE)
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load gene annotation
tx_annot <- fread(
    file.path("..", "..", "data", "Gierlinski_2015", "S-cerevisiae.transcripts.tsv"),
    sep = "\t",
    header = TRUE
)
named_tx <- tx_annot[gene_id != gene_name]

so_res <- list(
    "full" = fread(
        file.path("..", "..", "data", "Gierlinski_2015", "Sleuth", "genes.tsv"),
        sep = "\t",
        header = TRUE
    ),
    "even" = fread(
        file.path(RESULT_DIR, "Iterations", "Total_4", "balanced", "1.genes.tsv"),
        sep = "\t",
        header = TRUE
    ),
    "ols" = fread(
        file.path(RESULT_DIR, "Iterations", "Total_4", "unbalanced", "1.genes.tsv"),
        sep = "\t",
        header = TRUE
    ),
    "js" = fread(
        file.path(RESULT_DIR, "Iterations", "Total_4", "unbalanced", "1.genes.jse.tsv"),
        sep = "\t",
        header = TRUE
    )
)

# set FDR threshold
qval_thresh <- 0.01
log2fc_thresh <- 1


# ==============================================================================
# Analysis
# ==============================================================================
# combine into a single table for comparisons
so_tx <- rbindlist(
    lapply(
        names(so_res),
        function(comparison) {
            dt <- so_res[[comparison]]
            dt[, `:=`(
                # add Comparison column
                Comparison = comparison,
                # calculate log2 FC
                log2b = b / log(2)
            )]
            return(dt)
        }
    ),
    fill = TRUE
)

# convert to wide format for easier method comparison
so_tx_wide <- merge(
    x = dcast(
        so_tx,
        target_id ~ Comparison,
        value.var = "log2b"
    ),
    y = dcast(
        so_tx,
        target_id ~ Comparison,
        value.var = "qval"
    ),
    by = "target_id",
    suffixes = c("_log2b", "_qval")
)

# classify the transcript based on how the fold change and significance change between methods
so_tx_wide[, `:=`(
    All_Sig = (
        (even_qval < qval_thresh)
        & (full_qval < qval_thresh)
        & (ols_qval < qval_thresh)
        & (js_qval < qval_thresh)
    ),
    Unbal_Lost = (
        (even_qval < qval_thresh)
        & (full_qval < qval_thresh)
        & !(ols_qval < qval_thresh)
        & !(js_qval < qval_thresh)
    ),
    JS_Better = (
        (even_qval < qval_thresh)
        & (full_qval < qval_thresh)
        & !(ols_qval < qval_thresh)
        & (js_qval < qval_thresh)
    ),
    JS_Worse = (
        (even_qval < qval_thresh)
        & (full_qval < qval_thresh)
        & (ols_qval < qval_thresh)
        & !(js_qval < qval_thresh)
    )
)]

model_tx_ids <- so_tx_wide[
    (
        (All_Sig == TRUE)
        | (Unbal_Lost == TRUE)
        | (JS_Better == TRUE)
        | (JS_Worse == TRUE)
    ) & (target_id %in% named_tx$transcript_id),
    target_id
]
model_tx <- so_tx[target_id %in% model_tx_ids]

# merge annotation information
model_tx <- merge(
    x = model_tx,
    y = tx_annot[, .SD, .SDcols = c("transcript_id", "gene_name")],
    by.x = "target_id",
    by.y = "transcript_id"
)

so_tx <- merge(
    x = so_tx,
    y = tx_annot[, .SD, .SDcols = c("transcript_id", "gene_name")],
    by.x = "target_id",
    by.y = "transcript_id"
)

# set point colours based on fold change
so_tx[, Colour := ifelse(
    (qval > qval_thresh) | (abs(log2b) < log2fc_thresh),
    "N.S. or |log2(FC)| < 1",
    ifelse(
        b > 0,
        "Up in Mut",
        "Up in WT"
    )
)]

# make Comparison an ordered factor for easier plotting
so_tx[, Comparison := factor(
    Comparison,
    levels = c("full", "even", "ols", "js"),
    ordered = TRUE
)]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

# genes to be labelled in the volcano plots
labelled_genes <- c("SNF2", "PHO12", "TYE7", "TIS11")

# volcano plot of transcript fold changes
gg_vol <- (
    ggplot(
        data = so_tx,
        mapping = aes(
            x = log2b,
            y = -log10(pval),
            colour = Colour,
            label = gene_name
        )
    )
    + geom_point(
        alpha = 0.2
    )
    + geom_label_repel(
        data = so_tx[gene_name %in% labelled_genes],
        min.segment.length = 0
    )
    + scale_x_continuous(
        name = bquote(log[2] * "(" * Delta * "Snf2 / WT" * ")")
    )
    + scale_colour_manual(
        name = "Fold Change",
        breaks = c(
            "N.S. or |log2(FC)| < 1",
            "Up in Mut",
            "Up in WT"
        ),
        labels = c(
            "N.S. or |log2(FC)| < 1",
            bquote("Up in " * Delta * "Snf2"),
            "Up in WT"
        ),
        values = c(
            "#4B4B4B",
            "#F4A460",
            "#6495ED"
        )
    )
    + facet_wrap(
        ~ Comparison,
        scales = "free"
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
ggsave(
    file.path(PLOT_DIR, "volcano.png"),
    gg_vol,
    height = 12,
    width = 16,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "volcano.pdf"),
    gg_vol,
    height = 12,
    width = 16,
    units = "cm"
)
