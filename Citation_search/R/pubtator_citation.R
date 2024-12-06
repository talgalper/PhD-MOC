library(data.table)
library(fst)
library(progress)

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

target_set <- getBM(attributes = c( "external_gene_name", "description", "entrezgene_id"), 
                       filters = "external_gene_name", 
                       values = target_set, 
                       mart = ensembl)
target_set$description <- gsub("\\s*\\[.*?\\]", "", target_set$description)



citation_counts <- as.data.table(target_set)
citation_counts <- citation_counts[, count := NA]

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", 
                       total = nrow(target_set))
# Iterate over each row using an index
for (i in 1:nrow(citation_counts)) {
  gene_id <- citation_counts$entrezgene_id[i]
  gene_PMIDS <- pubtator3_neoplasm[Concept_ID == gene_id]
  counts <- uniqueN(gene_PMIDS$PMID)
  citation_counts$count[i] <- counts

  rm(counts, gene_id, gene_PMIDS, i)
  pb$tick()
}
rm(pb)

citation_counts <- citation_counts[order(external_gene_name)]
fwrite(citation_counts, "~/Desktop/temp.csv")



# repeat but for every gene in the set
unique_gene_ids <- unique(pubtator3_neoplasm$Concept_ID)
citation_counts <- data.table(entrezgene_id = unique_gene_ids,
                              counts = integer(length(unique_gene_ids)))

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
                       total = length(unique_gene_ids))

for (i in seq_along(unique_gene_ids)) {
  gene_id <- unique_gene_ids[i]
  gene_PMIDS <- pubtator3_neoplasm[Concept_ID == gene_id, PMID]
  counts <- uniqueN(gene_PMIDS)
  citation_counts$counts[i] <- counts
  
  rm(gene_PMIDS, gene_id, counts, i)
  pb$tick()
}

save(citation_counts, file = "citation_counts_pubtator.RData")

entrezgene_ids <- getBM(attributes = c( "entrezgene_id", "external_gene_name", "description", "gene_biotype"), 
                    filters = "entrezgene_id", 
                    values = citation_counts$entrezgene_id, 
                    mart = ensembl)
entrezgene_ids$description <- gsub("\\s*\\[.*?\\]", "", entrezgene_ids$description)





# create list of gfene name+description+synons
# dont need any more using entrezgene_id
target_set_list <- list()
for (gene in target_set$external_gene_name) {
  synons <- synon$synonym[synon$ref_term == gene]
  description <- target_set$description[target_set$external_gene_name == gene]
  
  target_set_list[[gene]] <- c(gene, description, synons)
  rm(synons, description, gene)
}







ESR1_synon <- synon[ref_term == "ESR1"]
ESR1_subset <- pubtator3[Mentions %in% c(ESR1_synon, "ESR1")]
ESR1_subset2 <- pubtator3[Concept_ID == "2099"]


KRAS_synon <- synon[ref_term == "KRAS"]
KRAS_subset <- pubtator3[Mentions %in% c(KRAS_synon, "KRAS")]

TP53_synon <- synon[ref_term == "TP53"]
TP53_subset <- pubtator3[Mentions %in% c(TP53_synon, "TP53")]

