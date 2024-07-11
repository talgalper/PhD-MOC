library(tidyverse)

# targets from API query (profile)
OpenTargets <- read_tsv("OpenTargets_data/breast_carcinoma_known_drugs.tsv")
#OpenTargets <- subset(OpenTargets, select = c("Target Approved Symbol", "Disease Name", "Drug Name", "Action Type", "Mechanism of Action","Drug Type"))

OpenTargets_filtered <- unique(OpenTargets) # drop duplicate rows
#OpenTargets_filtered <- OpenTargets_filtered[OpenTargets_filtered$`Drug Type` == "Small molecule", ] # subset small molecule drugs
#OpenTargets_filtered <- OpenTargets_filtered[OpenTargets_filtered$Status %in% c("Active, not recruiting", "Completed", "Recruiting"), ]

# either filter by unique drug name
OpenTargets_unique_drug <- OpenTargets_filtered[!duplicated(OpenTargets_filtered$`Drug Name`), ]

OpenTargets_unique_drug$NCT_ID <- sub(".*(NCT\\d+).*", "\\1", OpenTargets_unique_drug$URL) # isolate NCT IDs
OpenTargets_unique_drug$NCT_ID <- ifelse(grepl("^NCT", OpenTargets_unique_drug$NCT_ID), OpenTargets_unique_drug$NCT_ID, NA)

write.csv(OpenTargets_unique_drug, "../../../../Desktop/OpenTargets_unique_drug.csv")


# read in data from NCT
NCT_summaries <- read.csv("OpenTargets_data/NCT_summaries_expanded.csv")
colnames(NCT_summaries) <- c("NCT_ID", "trial_summary")

NCT_OpenTargets <- merge(OpenTargets_unique_drug, NCT_summaries, by = "NCT_ID")
NCT_OpenTargets <- NCT_OpenTargets[!duplicated(NCT_OpenTargets$`Drug Name`), ]
NCT_OpenTargets <- subset(NCT_OpenTargets, select = c("NCT_ID", "Phase", "Status", "Disease Name", "Drug Name",
                                                      "Mechanism of Action", "Action Type", "Target ID", "Target Approved Name", 
                                                      "Target Approved Symbol", "Drug Type", "trial_summary"))


write.csv(NCT_OpenTargets, "OpenTargets_data/OpenTargets_unique_drug.csv")










