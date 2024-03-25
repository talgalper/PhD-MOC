library(biomaRt)

MOC_CNV_data <- read.csv("../../../../../Desktop/kylie_cnv_pipe/Source data_All CN segments_with IDs.csv")
sample_info <- read.csv("../../RNA-Seq_pipeline/latest/MOC/data/All survival_CN_Aug18.csv")
sample_info <- subset(sample_info, select = c("GAMUT_ID", "Grade", "Stage"))

MOC_CNV_data <- merge(MOC_CNV_data, sample_info, by.x = "GAMuT.ID", by.y = "GAMUT_ID")

CNV_gain <- MOC_CNV_data[MOC_CNV_data$Event == c("CN Gain", "High Copy Gain"), ]

# Function to consolidate stages
consolidate_stages <- function(stages) {
  # Define a function to map stages to their consolidated forms
  map_stage <- function(stage) {
    if (grepl("^I[A-C]$", stage)) {
      return("I")
    } else if (grepl("^II[A-C]$", stage)) {
      return("II")
    } else if (grepl("^III[A-C]$", stage)) {
      return("III")
    } else if (grepl("^III[c]$", stage)) {
      return("III")
    } else {
      return(stage)
    }
  }
  
  # Apply the mapping function to each stage
  consolidated_stages <- sapply(stages, map_stage)
  
  return(consolidated_stages)
}

# add a column for new stage IDs
CNV_gain$stage <- consolidate_stages(CNV_gain$Stage)

# Replace empty strings with NA in the Stage column
CNV_gain$stage[CNV_gain$stage == ""] <- NA

# remove any rows with no gene terms
CNV_gain <- CNV_gain[CNV_gain$Gene.Symbols != "", ]

# Split merged gene symbols
split_genes <- strsplit(CNV_gain$Gene.Symbols, ", ")

# Create a list to store expanded data
expanded_data <- list()

# Expand the data by duplicating rows for each gene symbol
for (i in seq_along(split_genes)) {
  cat("Row:", i, "/", length(split_genes), "\n")
  expanded_data <- c(expanded_data, lapply(split_genes[[i]], function(x) {
    data.frame(
      Gene.Symbols = x,
      Stage = CNV_gain$stage[i]
    )
  }))
}

expanded_df <- do.call(rbind, expanded_data)

expanded_df$Stage[is.na(expanded_df$Stage)] <- "NA"

result <- table(expanded_df$Gene.Symbols, expanded_df$Stage)

# Calculate total counts for each gene
gene_total_counts <- rowSums(result)

# Calculate proportions for each stage
proportions <- prop.table(result, margin = 1) * 100

# Combine total counts and proportions into a dataframe
CN_gain_result_df <- data.frame(
  gene = rownames(proportions),
  counts = gene_total_counts,
  perc_stage_I = proportions[, "I"],
  perc_stage_II = proportions[, "II"],
  perc_stage_III = proportions[, "III"],
  perc_stage_IV = proportions[, "IV"],
  perc_NAs = proportions[, "NA"]
)

result_df_subset <- CN_gain_result_df[CN_gain_result_df$counts > 10, ]


