library(data.table)
library(fst)

# first time read in and read/write optimisation of pubator data
pubtator3 <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3", sep = "\t")
colnames(pubtator3) <- c("PMID", "Type", "Concept_ID", "Mentions", "Resource")
write_fst(pubtator3, "~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3.fst")


# read in data
pubtator3 <- read_fst("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3.fst") # mac
pubtator3 <- read_fst("/home/ubuntu/Downloads/gene2pubtator3.fst") # ubuntu
pubtator3 <- as.data.table(pubtator3)

pmids <- fread("data/pmids.txt", sep = "\t")
pmids <- pmids$V1
synon <- fread("data/human_syno.csv")

# isolate cancer related studies 
pubtator3_neoplasm <- pubtator3[PMID %in% pmids]



target_set <- c('BRCA1', 'TP53', 'ESR1', 'ERBB2', 'MYC', 'KIT', 'KRAS', 'AR', 'CD4', 'PIK3CA',
                'H4C4', 'GABBR2', 'F8', 'ALDH2', 'COL1A1', 'MYZAP', 'CENPK', 'KIF26B', 'USP25', 'CLOCK')

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

target_set <- getBM(attributes = c( "external_gene_name", "description"), 
                       filters = "external_gene_name", 
                       values = target_set, 
                       mart = ensembl)
target_set$description <- gsub("\\s*\\[.*?\\]", "", target_set$description)



# get PMIDS for 





ESR1_synon <- synon[ref_term == "ESR1"]
ESR1_subset <- pubtator3[Mentions %in% c(ESR1_synon, "ESR1")]

KRAS_synon <- synon[ref_term == "KRAS"]
KRAS_subset <- pubtator3[Mentions %in% c(KRAS_synon, "KRAS")]

TP53_synon <- synon[ref_term == "TP53"]
TP53_subset <- pubtator3[Mentions %in% c(TP53_synon, "TP53")]

