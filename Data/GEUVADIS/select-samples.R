# ==============================================================================
# Meta
# ==============================================================================
# select samples
# ------------------------------------------------
# Author: James Hawley
# Description: Figure out which samples to pick from the entire dataset


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load all sample metadata
meta <- fread("metadata.tsv")

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating")

# count unique individuals by the variables that matter
dt <- unique(meta[,
    .(
        Sex = get("Characteristics[sex]"),
        Lab = get("Characteristics[laboratory]"),
        Ancestry = get("Factor Value[ancestry category]"),
        Performer
    ),
    by = "Factor Value[individual]"
])
dt <- dt[, .N, by = c("Sex", "Lab", "Ancestry", "Performer")]

# extract counts by sex
dt <- dcast(dt, Lab + Ancestry + Performer ~ Sex, value.var = "N")

# identify samples with >= 12 male and >= 12 female samples
dt_chars <- dt[(male >= 12) & (female >= 12)]

# only pick these samples
# Ancestry: Finland
# Performer: UNIGE
# Lab: 1
selected_samples <- meta[
    (get("Factor Value[ancestry category]") == "Finland")
    & (Performer == "UNIGE")
    & (get("Characteristics[laboratory]") == 1),
    .SD,
    by = "Characteristics[sex]"
]

# select the first 12 of each
subset_ids <- c(
    selected_samples[
        get("Characteristics[sex]") == "male",
        unique(get("Factor Value[individual]"))
    ][1:12],
    selected_samples[
        get("Characteristics[sex]") == "female",
        unique(get("Factor Value[individual]"))
    ][1:12]
)

subset_samples <- meta[get("Factor Value[individual]") %in% subset_ids]

# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")
fwrite(subset_samples, "config.tsv", sep = "\t", col.names = TRUE)
