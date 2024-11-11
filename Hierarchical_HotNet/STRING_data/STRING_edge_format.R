library(data.table)
# STRING source data saved to onedrive

STRING_edge <- fread("~/Downloads/9606.protein.physical.links.full.v12.0.txt")
STRING_edge <- STRING_edge[, c(1,2,10)]
STRING_info <- fread("~/Downloads/9606.protein.info.v12.0.txt")
STRING_info <- STRING_info[, c(1,2)]
colnames(STRING_info) <- c("ENSP", "symbol")

STRING_edge_symbol <- merge.data.table(STRING_edge, STRING_info, by.x = "protein1", by.y = "ENSP", all.x = T)
colnames(STRING_edge_symbol)[4] <- "protein1_symbol"
STRING_edge_symbol <- merge.data.table(STRING_edge_symbol, STRING_info, by.x = "protein2", by.y = "ENSP", all.x = T, sort = F)
colnames(STRING_edge_symbol)[5] <- "protein2_symbol"
STRING_edge_symbol <- subset(STRING_edge_symbol, select = c("protein1_symbol", "protein2_symbol", "combined_score"))
STRING_edge_symbol <- STRING_edge_symbol[STRING_edge_symbol$combined_score >= 400, ]
STRING_edge_symbol$combined_score <- STRING_edge_symbol$combined_score / 1000

STRING_aliases <- fread("~/Downloads/9606.protein.aliases.v12.0.txt")
STRING_aliases <- STRING_aliases[STRING_aliases$source == "Ensembl_gene", ]
STRING_aliases <- STRING_aliases[, c(1,2)]

STRING_edge_ensemblGene <- merge.data.table(STRING_edge, STRING_aliases, by.x = "protein1", by.y = "#string_protein_id", all.x = T)
colnames(STRING_edge_ensemblGene)[4] <- "protein1_ENSG"
STRING_edge_ensemblGene <- merge.data.table(STRING_edge_ensemblGene, STRING_aliases, by.x = "protein2", by.y = "#string_protein_id", all.x = T, sort = F)
colnames(STRING_edge_ensemblGene)[5] <- "protein2_ENSG"
STRING_edge_ensemblGene <- subset(STRING_edge_ensemblGene, select = c("protein1_ENSG", "protein2_ENSG", "combined_score"))
STRING_edge_ensemblGene <- STRING_edge_ensemblGene[STRING_edge_ensemblGene$combined_score >= 400, ]
STRING_edge_ensemblGene$combined_score <- STRING_edge_ensemblGene$combined_score / 1000

fwrite(STRING_edge_symbol, "STRING_data/STRING_physical_geneSymbol.csv")
fwrite(STRING_edge_ensemblGene, "STRING_data/STRING_physical_ENSG.csv")
