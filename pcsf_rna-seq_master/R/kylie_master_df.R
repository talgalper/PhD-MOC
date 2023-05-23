

sample_info <- read.csv("kylie/data/All survival_CN_Aug18.csv")
sample_info <- subset(sample_info, select = c("GAMUT_ID", "Stage"))

kylie_counts <- read.csv("kylie/data/analysis_set_raw_counts.csv", row.names = 1)

# Extract column names
column_names <- colnames(kylie_counts)

# Remove "GAMuT_" part from column names
gamut_ids <- sub("^GAMuT_", "", column_names)


# create list of stage ids 
stage_I <- c("I","IA","IB","IC")
stage_II <- c("II", "IIA", "IIB", "IIC")
stage_III <- c("III", "IIIA", "IIIB", "IIIC", "IIIc")
stage_IV <- c("IV")
stage_ids <- list(stage_I, stage_II, stage_III, stage_IV)

# common GAMuT IDs between sample_info and kylie_counts
matching_ids <- intersect(gamut_ids, sample_info$GAMUT_ID)

stage_df <- list()

for (stage in stage_ids) {
  # sanity check
  print(stage)
  # retrieve GAMuT ID for the stage
  stage_id <- sample_info[sample_info$Stage %in% stage, 1]
  # common GAMuT IDs between stage and matching_ids
  matching_stage_id <- intersect(matching_ids, stage_id)
  # subset the kylie_counts df based on the IDs from matching_stage_id
  subset_stage_counts <- kylie_counts[, paste0("GAMuT_", matching_stage_id)]
  # add it to the df list
  stage_df <- c(stage_df, list(subset_stage_counts))
}

# retrieve each df from the list and save a separate variable
stage_I <- stage_df[[1]]
stage_II <- stage_df[[2]]
stage_III <- stage_df[[3]]
stage_IV <- as.data.frame(stage_df[[4]])
colnames(stage_IV) <- paste0("GAMuT_", matching_stage_id)
rownames(stage_IV) <- rownames(kylie_counts)


write.csv(stage_I, "kylie/stage_I_master_df.csv")
write.csv(stage_II, "kylie/stage_II_master_df.csv")
write.csv(stage_III, "kylie/stage_III_master_df.csv")
write.csv(stage_IV, "kylie/stage_IV_master_df.csv")




