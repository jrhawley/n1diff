# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

KALLISTO_DIR <- file.path("..", "..", "data", "CPC-GENE")
PLOT_DIR <- "Plots"

# ==============================================================================
# Function
# ==============================================================================
#' Save figures in multiple formats
#'
#' @param gg ggplot object
#' @param prefix Prefix for output file
#' @param ext Output extensions
#' @param dpi DPI resolution
savefig = function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
    for (e in ext) {
        ggsave(
            paste(prefix, e, sep = "."),
            gg,
            height = height,
            width = width,
            units = "cm",
            dpi = dpi
        )
    }
}


# ==============================================================================
# Data
# ==============================================================================
results <- fread("combined-results.tsv", sep = "\t", header = TRUE)
results[, method := factor(method, ordered = TRUE, levels = c("OLS", "JSE", "Full"))]

metadata <- fread(
    file.path(KALLISTO_DIR, "config.tsv"),
    sep = "\t",
    header = TRUE
)
SAMPLES <- metadata$SampleID

sim_results <- rbindlist(lapply(
    seq(3, 15),
    function(i) {
        dt <- fread(paste0("simulation.removing-T2E.", i, ".tsv"))
        dt[, N_transcripts := i]
        # remove outliers
        quartiles <- quantile(dt$Diff, c(0.25, 0.75))
        iqr <- quartiles[2] - quartiles[1]
        bounds <- c(median(dt$Diff) - 1.5 * iqr, median(dt$Diff) + 1.5 * iqr)
        dt[, Outlier := TRUE]
        dt[Diff >= bounds[1] & Diff <= bounds[2], Outlier := FALSE]
        return(dt)
    }
))
sim_results_long <- melt(
    sim_results,
    id.vars = c("Index", "N_transcripts", "Outlier"),
    variable.name = "Fit",
    value.name = "Error"
)

# ==============================================================================
# Plots
# ==============================================================================
for (s in SAMPLES) {
    print(s)
    gg <- (
        ggplot(data = results[(SampleID == s) | (SampleID == "All"), .SD, keyby = "method"])
        + geom_point(aes(x = beta_1, y = method, colour = method))
        + geom_path(aes(x = beta_1, y = method, group = target_id))
        + labs(x = expression(beta[1]), y = NULL, title = s)
        + theme_minimal()
    )
    savefig(gg, file.path(PLOT_DIR, s))
}

gg_sim <- (
    ggplot(data = sim_results[Outlier == FALSE & N_transcripts %in% c(3, 5, 10, 15)])
    + geom_histogram(aes(x = Diff), alpha = 0.8)
    + labs(x = expression(MSE[OLS] - MSE[JSE]), y = "Frequency")
    + facet_wrap(~ N_transcripts)
    + theme_minimal()
)
savefig(gg_sim, "Plots/simulation-diffs")
