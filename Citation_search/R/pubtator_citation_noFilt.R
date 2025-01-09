library(data.table)
library(progress)
library(fst)

# counts from python script
citation_counts <- fread("results/gene_counts_noFilt.txt", sep = "\t")

# separate multi mention-entries
multi_gene <- citation_counts[grepl(";", entrezgene_id)]
citation_counts <- citation_counts[!entrezgene_id %in% multi_gene$entrezgene_id]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

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
  
  batch_result <- suppressMessages(getBM(
    attributes = c("entrezgene_id", "entrezgene_accession", "entrezgene_description"),
    filters = "entrezgene_id",
    values = current_batch,
    mart = ensembl
  ))
  
  entrezgene_ids <- rbind(entrezgene_ids, batch_result)
  
  # Delay for 3 seconds so they dont freak out and ban me
  Sys.sleep(3)
  
  pb$tick()
  
  rm(start_index, end_index, current_batch, batch_result, i)
}
rm(pb)

save(entrezgene_ids, file = "data/entrezgene_id_converted.RData")

entrezgene_ids <- as.data.table(entrezgene_ids)
entrezgene_ids$entrezgene_id <- as.character(entrezgene_ids$entrezgene_id)

citation_counts <- merge.data.table(entrezgene_ids, citation_counts, by = "entrezgene_id", all.y = T)
citation_counts <- citation_counts[order(-count)]

citation_counts_noFilt <- citation_counts
save(citation_counts_noFilt, file = "results/citation_counts_pubtator_noFilt.RData")


# 2nd pass for mouse genes
# ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
# 
# mouse_genes <- data.table()
# NA_entrezgene_ids <- citation_counts_noFilt$entrezgene_id[is.na(citation_counts_noFilt$entrezgene_accession)]
# batch_size <- 10000
# num_batches <- ceiling(length(NA_entrezgene_ids) / batch_size)
# pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
#                        total = num_batches)
# for (i in seq_len(num_batches)) {
#   # Define the start and end indices for the current batch
#   start_index <- (i - 1) * batch_size + 1
#   end_index <- min(i * batch_size, length(NA_entrezgene_ids))
#   
#   # Extract the current batch of gene IDs
#   current_batch <- NA_entrezgene_ids[start_index:end_index]
#   
#   batch_result <- suppressMessages(getBM(
#     attributes = c("entrezgene_id", "entrezgene_accession", "entrezgene_description"),
#     filters = "entrezgene_id",
#     values = current_batch,
#     mart = ensembl
#   ))
#   
#   mouse_genes <- rbind(mouse_genes, batch_result)
#   
#   # Delay for 3 seconds so they dont freak out and ban me
#   Sys.sleep(3)
#   
#   pb$tick()
#   
# }
# rm(start_index, end_index, current_batch, batch_result, pb, i)
# 
# save(mouse_genes, file = "data/entrezgene_id_converted(mouse).RData")
# 
# mouse_genes <- as.data.table(mouse_genes)
# mouse_genes$entrezgene_id <- as.character(mouse_genes$entrezgene_id)
# 
# # merge with count data
# citation_counts_mouse <- merge.data.table(mouse_genes, citation_counts, by = "entrezgene_id", all.y = T)
# citation_counts_mouse <- citation_counts_mouse[order(-count)]
# 
# save(citation_counts_mouse, file = "results/citation_counts_pubtator_mouse.RData")


# instead of using biomaRt to convert everything, downloaded gene annotation data from entrez
entrezgene_data <- fread("/home/ubuntu/Downloads/All_Data.gene_info", sep = "\t")
temp <- entrezgene_data[1:10,1:16] # preview structure

entrezgene_data_subset <- entrezgene_data[, c(1,2,3,9,10)]
entrezgene_data_subset[] <- lapply(entrezgene_data_subset, as.character)

load("results/citation_counts_pubtator_noFilt.RData")
citaiton_counts_ognsmAnnot <- citation_counts_noFilt
citaiton_counts_ognsmAnnot <- merge.data.table(entrezgene_data_subset, citaiton_counts_ognsmAnnot, by.x = "GeneID", by.y = "entrezgene_id", all.y = T)
citaiton_counts_ognsmAnnot <- citaiton_counts_ognsmAnnot[order(-count)]
citaiton_counts_ognsmAnnot <- citaiton_counts_ognsmAnnot[, c(1:5,8)]
colnames(citaiton_counts_ognsmAnnot) <- c("entrezgene_id", "tax_id", "symbol", "description", "type", "count")
fwrite(citaiton_counts_ognsmAnnot, "results/citaiton_counts_ognsmAnnot.csv")

# subset of missing entries
# seem to have a withdraw notice on entrez page
temp <- citaiton_counts_ognsmAnnot[is.na(tax_id)]


# annotate the tax_ids with the organism scientific names
# NEED TO ADD ERROR HANDLING FOR TAX_IDs THAT DO NOT RETURN DATA
library(rentrez)
tax_ids <- citaiton_counts_ognsmAnnot[,unique(tax_id)]
tax_ids <- na.omit(tax_ids)

get_organism_name <- function(entrez_id) {
  summary <- entrez_summary(db = "taxonomy", id = entrez_id)
  return(summary$scientificname)
}

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta",
                       total = length(tax_ids))
organisms <- data.table(tax_id = character(), organism = character())
for (id in tax_ids) {
  organism <- get_organism_name(id)
  organisms <- rbind(organisms, data.table(tax_id = id, organism = organism))
  
  pb$tick()
}
rm(pb, organism, id)







# annotate tax_id with organism name
library(data.table)
library(progress)

# spelling error in filename (fix later)
citation_counts_ognsmAnnot <- fread("/home/ubuntu/Desktop/pubtator3/citaiton_counts_ognsmAnnot.csv") # ubuntu
citation_counts_ognsmAnnot <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/citaiton_counts_ognsmAnnot.csv") # mac

library(rentrez)
tax_ids <- citation_counts_ognsmAnnot[,unique(tax_id)]
tax_ids <- na.omit(tax_ids)

get_organism_name <- function(entrez_id) {
  tryCatch({
    summary <- entrez_summary(db = "taxonomy", id = entrez_id)
    summary$scientificname
  }, error = function(e) {
    NA  # Return NA if there's an error
  })
}

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta",
                       total = length(tax_ids))

organisms <- data.table(tax_id = character(), organism = character())
for (id in tax_ids) {
  organism <- get_organism_name(id)
  organisms <- rbind(organisms, data.table(tax_id = id, organism = organism))
  
  pb$tick()
}
rm(pb, organism, id)

tmp <- get_organism_name("3248881")

citation_counts_ognsmAnnot$tax_id <- as.character(citation_counts_ognsmAnnot$tax_id)
citation_counts_ognsmAnnot$entrezgene_id <- as.character(citation_counts_ognsmAnnot$entrezgene_id)
citation_counts_ognsmAnnot <- merge.data.table(organisms, citation_counts_ognsmAnnot, by = "tax_id", all.x = T)
citation_counts_ognsmAnnot <- citation_counts_ognsmAnnot[order(-citation_counts_ognsmAnnot$count), ]

temp <- copy(citation_counts_ognsmAnnot)
temp[, combined := ifelse(symbol == "", entrezgene_id, symbol)]

# apparently much faster than a loop. Can confirm it fukn defs is
citation_counts <- temp[
  , .(counts = sum(count)),
  by = combined
]
citation_counts <- citation_counts[order(-citation_counts$counts), ]

rm(temp)
gc()


plot_data <- citation_counts[citation_counts$counts > 50, ]
summary(plot_data$counts)

hist(log10(plot_data$counts))
summary(log10(citation_counts$counts))


