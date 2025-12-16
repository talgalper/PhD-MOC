

# load in STRING data outside function because we can just copy it to each dir
STRING_edge <- read.csv("STRING_data/STRING_physical_ENSG.csv") 
STRING_edge <- STRING_edge[, -3]
STRING_edge <- STRING_edge[!duplicated(t(apply(STRING_edge, 1, sort))), ]

preProcessing <- function(subtype) {
  
  # create dirs if they dont exist
  dir.create(paste0("BRCA/subtype_analysis/", subtype, "/data"), recursive = TRUE, showWarnings = FALSE)
  
  # create index and indexed edge list
  gene_index <- as.data.frame(unique(c(STRING_edge$protein1_ENSG, STRING_edge$protein2_ENSG)))
  gene_index$row_num <- seq.int(nrow(gene_index))
  colnames(gene_index) <- c("ensembl_id", "index")
  gene_index <- gene_index[c("index", "ensembl_id")]
  
  edge_list_index <- data.frame(from = match(STRING_edge$protein1_ENSG, gene_index$ensembl_id),
                                to = match(STRING_edge$protein2_ENSG, gene_index$ensembl_id))
  
  # create score file
  DE_data <- read.csv(paste0("../BRCA_pipe/latest_run/subtype_analysis/", subtype, "_difExp.csv"))
  DE_data_abs <- subset(DE_data, select = c("gene_id", "logFC_abs"))
  DE_data_signed <-  subset(DE_data, select = c("gene_id", "logFC"))
  
  write.table(gene_index, paste0("BRCA/subtype_analysis/", subtype, "/data/gene_index.tsv"), col.names = F, row.names = F, sep = "\t")
  write.table(edge_list_index, paste0("BRCA/subtype_analysis/", subtype, "/data/edge_list_index.tsv"), col.names = F, row.names = F, sep = "\t")
  write.table(DE_data_abs, paste0("BRCA/subtype_analysis/", subtype, "/data/logFC_scores_abs.tsv"), col.names = F, row.names = F, sep = "\t")
  write.table(DE_data_signed, paste0("BRCA/subtype_analysis/", subtype, "/data/logFC_scores_signed.tsv"), col.names = F, row.names = F, sep = "\t")
}

subtypes <- c("lumA", "lumB", "Her2", "basal")

for (subtype in subtypes) {
  preProcessing(subtype)
  
  # create intermediate filder for matrix
  dir.create(paste0("BRCA/subtype_analysis/", subtype, "/intermediate/STRING"), recursive = TRUE, showWarnings = FALSE)
}
