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

save(entrezgene_ids, file = "entrezgene_id_converted.RData")

entrezgene_ids <- as.data.table(entrezgene_ids)
entrezgene_ids$entrezgene_id <- as.character(entrezgene_ids$entrezgene_id)

citation_counts <- merge.data.table(entrezgene_ids, citation_counts, by = "entrezgene_id", all.y = T)
citation_counts <- citation_counts[order(-count)]

citation_counts_noFilt <- citation_counts
save(citation_counts_noFilt, file = "results/citation_counts_pubtator_noFilt.RData")


# 2nd pass for mouse genes
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

mouse_genes <- data.table()
NA_entrezgene_ids <- citation_counts_noFilt$entrezgene_id[is.na(citation_counts_noFilt$entrezgene_accession)]
num_batches <- ceiling(length(NA_entrezgene_ids) / batch_size)
pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
                       total = num_batches)
for (i in seq_len(num_batches)) {
  # Define the start and end indices for the current batch
  start_index <- (i - 1) * batch_size + 1
  end_index <- min(i * batch_size, length(NA_entrezgene_ids))
  
  # Extract the current batch of gene IDs
  current_batch <- NA_entrezgene_ids[start_index:end_index]
  
  batch_result <- suppressMessages(getBM(
    attributes = c("entrezgene_id", "entrezgene_accession", "entrezgene_description"),
    filters = "entrezgene_id",
    values = current_batch,
    mart = ensembl
  ))
  
  mouse_genes <- rbind(mouse_genes, batch_result)
  
  # Delay for 3 seconds so they dont freak out and ban me
  Sys.sleep(3)
  
  pb$tick()
  
  rm(start_index, end_index, current_batch, batch_result, i)
}
rm(pb)












# Resume at batch 275:
for (i in 430:num_batches) {
  start_index <- (i - 1) * batch_size + 1
  end_index <- min(i * batch_size, length(NA_entrezgene_ids))
  
  current_batch <- NA_entrezgene_ids[start_index:end_index]
  
  batch_result <- suppressMessages(getBM(
    attributes = c("entrezgene_id", "entrezgene_accession", "entrezgene_description"),
    filters = "entrezgene_id",
    values = current_batch,
    mart = ensembl
  ))
  
  mouse_genes <- rbind(mouse_genes, batch_result)
  
  Sys.sleep(3)
  
  pb$tick()
  
  # Instead of removing i at the end of each iteration, you can remove it after the loop if needed.
  rm(start_index, end_index, current_batch, batch_result)
}

rm(pb)

save(mouse_genes, file = "entrezgene_id_mouse_429it.RData")
