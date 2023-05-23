
# read in data
af_drugability <- read.csv("data/druggability/fpocket_druggability.csv")
weighted_pcsf <- read.csv("results/weighted_pcsf_results.csv")
unweighted_pcsf <- read.csv("results/unweighted_pcsf_results.csv")

# merge AF data with PCSF data by uniprot ID
weighted_pcsf_master <- merge(weighted_pcsf, af_drugability, by.x = "uniprot_id")
weighted_pcsf_master <- subset(weighted_pcsf_master, select = c("ensembl_peptide_id", "uniprot_id", "cluster", "betweenness", 
                                                                "degree_centrality", "pocket", "druggability", "num_drug_pockets"))
weighted_pcsf_master <- weighted_pcsf_master[order(-weighted_pcsf_master$druggability), ]

write.csv(weighted_pcsf_master, "results/weighted_PCSF_druggability.csv")

# repeat for unweighted
unweighted_pcsf_master <- merge(unweighted_pcsf, af_drugability, by.x = "uniprot_id")
unweighted_pcsf_master <- subset(unweighted_pcsf_master, select = c("ensembl_peptide_id", "uniprot_id", "cluster", "betweenness", 
                                                                "degree_centrality", "pocket", "druggability", "num_drug_pockets"))
unweighted_pcsf_master <- unweighted_pcsf_master[order(-unweighted_pcsf_master$druggability), ]

write.csv(unweighted_pcsf_master, "results/unweighted_PCSF_druggability.csv")
