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


# ==============================================================================
# Functions
# ==============================================================================


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")

# load fully balanced 12-vs-12 results
full <- fread("sleuth.transcripts.tsv")

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
loginfo("Comparing test results to full comparison")

# calculate confusion matrices for tests with different totals
conf_mtx <- lapply(
    c(4, 6, 8, 10, 12),
    function(total) {
        # subset by the total number of samples
        dt <- tests[Total == total]
        
    }
)


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")


# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")

