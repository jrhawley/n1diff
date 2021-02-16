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
ggsave("Plots/p-values.genes.png", gg_pval, height = 8, width = 12, units = "cm")

# volcano plot of transcript fold changes
gg_vol <- (
    ggplot(data = so_genes)
    + geom_point(
        aes(x = log2b, y = -log10(pval), colour = Colour),
        alpha = 0.2
    )
    + scale_x_continuous(
        name = bquote(log[2] * "(" * Delta * "snf2 / WT" * ")")
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
            bquote("Up in " * Delta * "snf2"),
            "Up in WT"
        ),
        values = c(
            "#4B4B4B",
            "#F4A460",
            "#6495ED"
        )
    )
    + theme_minimal()
)
ggsave("Plots/volcano.png", gg_vol, height = 8, width = 12, units = "cm")
