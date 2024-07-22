
OpenTargets <- read_tsv("OpenTargets_data/breast_carcinoma_known_drugs.tsv")



OpenTargets_unique$anti_cancer <- ifelse(grepl(paste(c("pain", "diarrhoea", "symptoms", "diarrhea", "dermatitis", 
                                                       "vomiting", "manage"), collapse = "|"), 
                                               OpenTargets_unique$trial_summary, ignore.case = T) | 
                                           !grepl(paste(c("breast", "TNBC"), collapse = "|"), OpenTargets_unique$trial_summary, ignore.case = TRUE), "no", "yes")

OpenTargets_unique <- OpenTargets_unique[OpenTargets_unique$anti_cancer == "yes", ]
OpenTargets_unique <- OpenTargets_unique[!duplicated(OpenTargets_unique$Target.Approved.Symbol), ]


OpenTargets_SM <- OpenTargets_unique[OpenTargets_unique$Drug.Type == "Small molecule", ]

write.csv(OpenTargets_SM, "OpenTargets_data/OpenTargets_updated.csv")