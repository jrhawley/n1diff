# ==============================================================================
# Meta
# ==============================================================================
# plot-sleuth
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

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
so_genes <- fread(
    file.path("Sleuth", "genes.tsv"),
    sep = "\t",
    header = TRUE
)
so_genes[, log2b := b / log(2)]

# set FDR threshold
qval_thresh <- 0.01
equiv_pval_thresh <- so_genes[qval > qval_thresh, min(pval)]
log2fc_thresh <- 1

# load gene annotation
tx <- fread(
    "S-cerevisiae.transcripts.tsv",
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Analysis
# ==============================================================================
# set point colours based on fold change
so_genes[, Colour := ifelse(
    (qval > qval_thresh) | (abs(log2b) < log2fc_thresh),
    "N.S. or |log2(FC)| < 1",
    ifelse(
        b > 0,
        "Up in Mut",
        "Up in WT"
    )
)]

# merge annotation information
so_genes <- merge(
    x = so_genes,
    y = tx,
    by.x = "target_id",
    by.y = "transcript_id"
)


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
# p-value histograms
gg_pval <- (
    ggplot(data = so_genes)
    + geom_histogram(
        aes(x = pval),
        bins = 100
    )
    + scale_x_continuous(
        name = "p-value"
    )
    + scale_y_continuous(
        name = "Frequency"
    )
    + theme_minimal()
)
ggsave(
    "Plots/p-values.genes.png",
    gg_pval,
    height = 8,
    width = 12,
    units = "cm"
)

# volcano plot of transcript fold changes
labelled_genes <- c("SNF2", "PHO12", "TYE7", "TIS11")
gg_vol <- (
    ggplot(
        data = so_genes,
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
        data = so_genes[gene_name %in% labelled_genes],
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
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/volcano.png",
    gg_vol,
    height = 12,
    width = 16,
    units = "cm"
)
ggsave(
    "Plots/volcano.pdf",
    gg_vol,
    height = 12,
    width = 16,
    units = "cm"
)
