# ==============================================================================
# Meta
# ==============================================================================
# plot-sleuth
# ------------------------------------------------
# Author: James Hawley
# Description: Plot differential analyses results


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# load fully balanced 12-vs-12 results
full <- fread("sleuth.transcripts.tsv")
full_genes <- fread("sleuth.genes.tsv")
full_so <- readRDS("sleuth.object.rds")
full_mat <- sleuth_to_matrix(full_so, "obs_norm", "est_counts")

# load smaller sample tests
tests <- rbindlist(lapply(
    c("balanced", "unbalanced"),
    function(bal) {
        dt1 <- rbindlist(lapply(
            c(4, 6, 8, 10, 12),
            function(total) {
                dt2 <- rbindlist(lapply(
                    1:10,
                    function(i) {
                        dt3 <- fread(
                            file.path(
                                "Iterations",
                                paste0("Total_", total),
                                bal,
                                paste0(i, ".transcripts.tsv")
                            ),
                            sep = "\t",
                            header = TRUE
                        )
                        dt3[, Iteration := i]
                        return(dt3)
                    }
                ))
                dt2[, Total := total]
                return(dt2)
            }
        ))
        dt1[, Balanced := (bal == "balanced")]
        return(dt1)
    }
))


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating regularized logarithms")

# convert wide estimated counts matrix to long form
full_mat_long <- melt(
    as.data.table(full_mat),
    measure.vars = colnames(full_mat),
    variable.name = "Sample_ID",
    value.name = "est_counts"
)

# merge metadata
full_mat_long <- merge(
    x = full_mat_long,
    y = meta[, .SD, .SDcols = c("sample", "condition")],
    by.x = "Sample_ID",
    by.y = "sample"
)

rle_mat = data.table()
for (i in 1:dim(full_mat)[2]) {
    v = full_mat[, i]
    med_v = median(v, na.rm = TRUE)
    rle_mat[, V1 := log2((v + 1) / (med_v + 1))]
    colnames(rle_mat)[i] = colnames(full_mat)[i]
}

# convert wide estimated counts matrix to long form
rle_mat_long <- melt(
    rle_mat,
    measure.vars = colnames(rle_mat),
    variable.name = "Sample_ID",
    value.name = "est_counts"
)

# merge metadata
rle_mat_long <- merge(
    x = rle_mat_long,
    y = meta[, .SD, .SDcols = c("sample", "condition")],
    by.x = "Sample_ID",
    by.y = "sample"
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

gg_est_counts <- (
    ggplot(data = rle_mat_long)
    + geom_boxplot(
        aes(x = Sample_ID, y = est_counts, fill = condition),
        outlier.shape = NA
    )
    + scale_x_discrete(
        name = "Sample"
    )
    + scale_y_continuous(
        name = "Estimated Counts"
        # limits = c(0, 10)
    )
    + scale_fill_manual(
        name = "Sex",
        breaks = c("male", "female"),
        labels = c("Male", "Female"),
        values = c("#1e90ff", "#dda0dd")
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/full.counts.png",
    gg_est_counts,
    height = 8,
    width = 12,
    units = "cm"
)

# plot p-values for full comparison
gg_p_full <- (
    ggplot(data = full_genes)
    + geom_histogram(aes(x = pval))
    + scale_x_continuous(
        name = "p-value"
    )
    + scale_y_continuous(
        name = "Frequency"
    )
)
ggsave("Plots/full.p-values.png", gg_p_full, height = 8, width = 12, units = "cm")

# plot p-values for each test

