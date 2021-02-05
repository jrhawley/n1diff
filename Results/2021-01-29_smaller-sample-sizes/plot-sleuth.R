# ==============================================================================
# Meta
# ==============================================================================
# plot-sleuth
# ------------------------------------------------
# Author: James Hawley
# Description: Plot comparisons between full scale dataset and smaller sample sizes

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

# load data
confusion <- fread(
    file.path("Comparison", "confusion.tsv"),
    sep = "\t",
    header = TRUE
)
rates <- fread(
    file.path("Comparison", "rates.tsv"),
    sep = "\t",
    header = TRUE
)

mse_all <- fread(
    file.path("Comparison", "mse.all.tsv"),
    sep = "\t",
    header = TRUE
)

# data table for making plots more intuitive
stat_better_direction <- data.table(
    Stat = c(
        "ACC",
        "MCC",
        "NPV",
        "PPV",
        "TNR",
        "TPR"
    ),
    Direction = c(
        "higher",
        "higher",
        "lower",
        "lower",
        "lower",
        "higher",
        "higher",
        "higher",
        "higher",
        "higher"
    )
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
for (s in colnames(rates)[8:17]) {

    better_dir <- stat_better_direction[Stat == s, Direction]

    gg <- (
        ggplot(
            data = rates,
            aes(x = Total, y = get(s), colour = Test_Condition)
        )
        + geom_boxplot(
            mapping = aes(
                fill = Test_Condition,
                group = paste(Total, Test_Condition)
            ),
            colour = "#000000",
            alpha = 0.5,
            outlier.shape = NA,
            position = position_dodge(width = 1.5)
        )
        # + geom_point(
        #     position = position_jitterdodge(
        #         jitter.width = 0.6,
        #         dodge.width = 1.5
        #     ),
        #     alpha = 0.5
        # )
        + scale_x_continuous(
            name = "Total Number of Samples",
            breaks = seq(4, 24, 2),
            labels = seq(4, 24, 2)
        )
        + scale_y_continuous(
            name = paste0(s," (", better_dir, " is better)")
            # limits = c(bounds[1], bounds[2])
        )
        + scale_colour_manual(
            name = "Design",
            breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
            labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
            values = c("#4B4B4B", "#DA8FD6", "#FA8072")
        )
        + scale_fill_manual(
            name = "Design",
            breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
            labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
            values = c("#4B4B4B", "#DA8FD6", "#FA8072")
        )
        + theme_minimal()
        + theme(
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        )
    )
    ggsave(
        file.path("Plots", paste0(s, ".png")),
        width = 12,
        height = 8,
        units = "cm"
    )
}

gg <- (
    ggplot(
        data = confusion,
        aes(x = Total, y = N, colour = Test_Condition)
    )
    + geom_boxplot(
        mapping = aes(
            fill = Test_Condition,
            group = paste(Total, Test_Condition)
        ),
        colour = "#000000",
        alpha = 0.5,
        outlier.shape = NA,
        position = position_dodge(width = 1.5)
    )
    # + geom_point(
    #     position = position_jitterdodge(
    #         jitter.width = 0.6,
    #         dodge.width = 1.5
    #     ),
    #     alpha = 0.5
    # )
    + scale_x_continuous(
        name = "Total Number of Samples",
        breaks = seq(4, 24, 2),
        labels = seq(4, 24, 2)
    )
    + scale_y_continuous(
        name = "Frequency"
        # limits = c(bounds[1], bounds[2])
    )
    + scale_colour_manual(
        name = "Design",
        breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
        labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
        values = c("#4B4B4B", "#DA8FD6", "#FA8072")
    )
    + scale_fill_manual(
        name = "Design",
        breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
        labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
        values = c("#4B4B4B", "#DA8FD6", "#FA8072")
    )
    + facet_wrap(~ Result, scales = "free_y")
    + theme_minimal()
    + theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    file.path("Plots", "confusion.png"),
    width = 12,
    height = 8,
    units = "cm"
)

gg_mse <- (
    ggplot(data = mse_all, mapping = aes(x = Test_Condition, y = MSE, fill = Test_Condition))
    + geom_col()
    + scale_x_discrete(
        name = "Method"
    )
    + scale_y_continuous(
        name = "Mean Square Error"
    )
    + scale_fill_manual(
        name = "Design",
        breaks = c("Balanced", "Unbalanced OLS", "Unbalanced JS"),
        labels = c("Even split", "1-vs-N", "1-vs-N James-Stein"),
        values = c("#4B4B4B", "#DA8FD6", "#FA8072")
    )
    + guides(fill = FALSE)
    + theme_minimal()
)
ggsave(
    file.path("Plots", "mse.png"),
    width = 12,
    height = 8,
    units = "cm"
)
