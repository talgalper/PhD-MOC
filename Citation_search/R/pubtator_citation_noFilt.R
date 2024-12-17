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
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

mouse_genes <- data.table()
NA_entrezgene_ids <- citation_counts_noFilt$entrezgene_id[is.na(citation_counts_noFilt$entrezgene_accession)]
batch_size <- 10000
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
  
}
rm(start_index, end_index, current_batch, batch_result, pb, i)

save(mouse_genes, file = "data/entrezgene_id_converted(mouse).RData")

mouse_genes <- as.data.table(mouse_genes)
mouse_genes$entrezgene_id <- as.character(mouse_genes$entrezgene_id)

# merge with count data
citation_counts_mouse <- merge.data.table(mouse_genes, citation_counts, by = "entrezgene_id", all.y = T)
citation_counts_mouse <- citation_counts_mouse[order(-count)]

save(citation_counts_mouse, file = "results/citation_counts_pubtator_mouse.RData")


# using new library to annotate entrez IDs with organism
library(rentrez)
library(future.apply) # need to parallelise
library(progressr)


# Define the function to fetch organism information for a single Entrez ID
get_organism <- function(entrez_id) {
  tryCatch({
    summary <- entrez_summary(db = "gene", id = entrez_id)
    return(summary$organism$scientificname)
  }, error = function(e) {
    return(NA)  # Return NA if there's an error
  })
}

# Set up parallel processing and progress handler
plan(multisession)  # Use multisession for Mac/RStudio stability
handlers(global = TRUE)  # Enable global progress handlers

# Enable progress bar
with_progress({
  p <- progressor(along = citation_counts_orgnsmAnnot$entrezgene_id)
  
  # Parallelised organism retrieval
  results <- future_lapply(citation_counts_orgnsmAnnot$entrezgene_id, function(entrez_id) {
    p(message = paste("Processing Entrez ID:", entrez_id))
    organism <- get_organism(entrez_id)
    return(data.table(entrezgene_id = entrez_id, organism = organism))
  })
})

# Combine results into a single data.table
organisms <- rbindlist(results)

# Print final results
print(organisms)

# Clean up
plan(sequential)











