library(data.table)
library(fst)
library(progress)

# first time read in and read/write optimisation of pubator data
gene2pubtator3 <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3", sep = "\t")
colnames(gene2pubtator3) <- c("PMID", "Type", "Concept_ID", "Mentions", "Resource")
write_fst(gene2pubtator3, "~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3.fst")

species2pubtator3 <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/species2pubtator3", sep = "\t")
colnames(species2pubtator3) <- c("PMID", "Type", "Concept_ID", "Mentions", "Resource")
write_fst(species2pubtator3, "~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/species2pubtator3.fst")

# read in data
gene2pubtator3 <- read_fst("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/gene2pubtator3.fst") # mac
gene2pubtator3 <- read_fst("/home/ubuntu/Downloads/gene2pubtator3.fst") # ubuntu
gene2pubtator3 <- as.data.table(gene2pubtator3)

species2pubtator3 <- read_fst("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/species2pubtator3.fst") # mac
species2pubtator3 <- read_fst("/home/ubuntu/Downloads/species2pubtator3.fst") # ubuntu
species2pubtator3 <- as.data.table(species2pubtator3)

# subset PMIDs for humans
species2pubtator3_human <- species2pubtator3[Concept_ID == "9606"]
species2pubtator3_human_PMIDs <- species2pubtator3_human$PMID

# subset gene mentions for humans
gene2pubtator3_human <- gene2pubtator3[PMID %in% species2pubtator3_human_PMIDs]

# read in PMIDs for "Human[MH] AND Noeplasms[MH]" search
pmids <- fread("data/pmids.txt", sep = "\t")
pmids <- pmids$V1
#synon <- fread("data/human_syno.csv")

# subset Neoplasm related literature
gene2pubtator3_humanNeoplasm <- gene2pubtator3_human[PMID %in% pmids]

# Count number of PMIDs for each entrez_id
unique_gene_ids <- unique(gene2pubtator3_humanNeoplasm$Concept_ID)
citation_counts <- data.table(entrezgene_id = unique_gene_ids,
                              counts = integer(length(unique_gene_ids)))

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
                       total = length(unique_gene_ids))

for (i in seq_along(unique_gene_ids)) {
  gene_id <- unique_gene_ids[i]
  gene_PMIDS <- gene2pubtator3_humanNeoplasm[Concept_ID == gene_id, PMID]
  counts <- uniqueN(gene_PMIDS)
  citation_counts$counts[i] <- counts
  
  rm(gene_PMIDS, gene_id, counts, i)
  pb$tick()
}

save(citation_counts, file = "results/citation_counts_pubtator.RData")
load("results/citation_counts_pubtator.RData")



# separate multi mention-entries
multi_gene <- citation_counts[grepl(";", entrezgene_id)]
citation_counts <- citation_counts[!entrezgene_id %in% multi_gene$entrezgene_id]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 113)

# break up query into batches
entrezgene_ids <- data.frame()
batch_size <- 10000
num_batches <- ceiling(length(citation_counts$entrezgene_id) / batch_size)
pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
                       total = num_batches)
for (i in seq_len(num_batches)) {
  # Define the start and end indices for the current batch
  start_index <- (i - 1) * batch_size + 1
  end_index <- min(i * batch_size, length(citation_counts$entrezgene_id))
  
  # Extract the current batch of gene IDs
  current_batch <- citation_counts$entrezgene_id[start_index:end_index]
  
  batch_result <- getBM(
    attributes = c("entrezgene_id", "entrezgene_accession", "entrezgene_description"),
    filters = "entrezgene_id",
    values = current_batch,
    mart = ensembl
    )
  
  entrezgene_ids <- rbind(entrezgene_ids, batch_result)
  
  # Delay for 5 seconds so they dont freak out
  Sys.sleep(5)
  
  pb$tick()
  
  rm(start_index, end_index, current_batch, batch_result)
}


citation_counts <- merge(entrezgene_ids, citation_counts, by = "entrezgene_id", all.y = T)
citation_counts <- as.data.table(citation_counts)
citation_counts <- citation_counts[order(-counts)]


batch_result <- getBM(
  attributes = c("entrezgene_id", "entrezgene_accession", "entrezgene_description"),
  filters = "entrezgene_id",
  values = current_batch[1:10],
  mart = ensembl
)






# small example from data set
target_set <- c('BRCA1', 'TP53', 'ESR1', 'ERBB2', 'MYC', 'KIT', 'KRAS', 'AR', 'CD4', 'PIK3CA',
                'H4C4', 'GABBR2', 'F8', 'ALDH2', 'COL1A1', 'MYZAP', 'CENPK', 'KIF26B', 'USP25', 'CLOCK')

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

target_set <- getBM(attributes = c( "external_gene_name", "description", "entrezgene_id"), 
                       filters = "external_gene_name", 
                       values = target_set, 
                       mart = ensembl)
target_set$description <- gsub("\\s*\\[.*?\\]", "", target_set$description)
target_set$entrezgene_id <- as.character(target_set$entrezgene_id)


target_set_counts <- merge.data.table(target_set, citation_counts, by = "entrezgene_id", all.x = T)











ESR1_synon <- synon[ref_term == "ESR1"]
ESR1_subset <- pubtator3[Mentions %in% c(ESR1_synon, "ESR1")]
ESR1_subset2 <- pubtator3[Concept_ID == "2099"]


KRAS_synon <- synon[ref_term == "KRAS"]
KRAS_subset <- pubtator3[Mentions %in% c(KRAS_synon, "KRAS")]

TP53_synon <- synon[ref_term == "TP53"]
TP53_subset <- pubtator3[Mentions %in% c(TP53_synon, "TP53")]

