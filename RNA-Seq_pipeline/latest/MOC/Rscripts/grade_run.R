library(PCSF)
library(plyr)
library(biomaRt)
library(dplyr)
library(tidyr)
library(grid)
library(edgeR)
library(ggplot2)
library(reshape2)
library(progress)
library(RobustRankAggreg)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rm(list = ls()[!ls() %in% c("ensembl")])



sample_info <- read.csv("data/All survival_CN_Aug18.csv")
sample_info <- subset(sample_info, select = c("GAMUT_ID", "Grade", "Stage"))

kylie_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)

# Extract column names
column_names <- colnames(kylie_counts)

# Remove "GAMuT_" part from column names
gamut_ids <- sub("^GAMuT_", "", column_names)


# create list of stage ids 
grade_1 <- c("1")
grade_2 <- c("2")
grade_3 <- c("3")
BDL <- c("BDL")
BEN <- c("BEN")
EOM <- c("EOM")
unlabelled <- c("")


grade_ids <- list(grade_1, grade_2, grade_3, BDL, BEN, EOM, unlabelled)

# common GAMuT IDs between sample_info and kylie_counts
matching_ids <- intersect(gamut_ids, sample_info$GAMUT_ID)

grade_df <- list()

for (grade in grade_ids) {
  # sanity check
  print(grade)
  # retrieve GAMuT ID for the stage
  grade_id <- sample_info[sample_info$Grade %in% grade, 1]
  # common GAMuT IDs between stage and matching_ids
  matching_grade_id <- intersect(matching_ids, grade_id)
  # subset the kylie_counts df based on the IDs from matching_stage_id
  subset_grade_counts <- kylie_counts[, paste0("GAMuT_", matching_grade_id)]
  # add it to the df list
  grade_df <- c(grade_df, list(subset_grade_counts))
}

# retrieve each df from the list and save a separate variable
grade_1 <- grade_df[[1]]
grade_2 <- grade_df[[2]]
grade_3 <- grade_df[[3]]
BDL <- grade_df[[4]]
BEN <- grade_df[[5]]
EOM <- grade_df[[6]]

unlabelled <- as.data.frame(grade_df[[7]])
colnames(unlabelled) <- paste0("GAMuT_", matching_grade_id)
rownames(unlabelled) <- rownames(kylie_counts)



moc_data <- cbind(grade_1, grade_2, grade_3, BDL, BEN, unlabelled)
