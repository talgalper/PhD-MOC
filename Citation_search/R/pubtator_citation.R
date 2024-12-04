library(data.table)
library(fst)

# first time read in and read/write optimisation of pubator data
pubtator3 <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3", sep = "\t")
colnames(pubtator3) <- c("PMID", "Type", "Concept_ID", "Mentions", "Resource")
write_fst(pubtator3, "~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3.fst")



pubtator3 <- read_fst("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3.fst")
pmids <- fread("data/pmids.txt", sep = "\t")
synon <- fread("data/human_syno.csv")

pubtator3_genes <- pubtator3[!duplicated(Mentions)]





KRAS_synon <- synon[ref_term == "KRAS"]
KRAS_subset <- pubtator3[Mentions %in% c(KRAS_synon, "KRAS")]

TP53_synon <- synon[ref_term == "TP53"]
TP53_subset <- pubtator3[Mentions %in% c(TP53_synon, "TP53")]

