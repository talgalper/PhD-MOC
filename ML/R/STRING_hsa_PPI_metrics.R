library(data.table)
# STRING source data saved to onedrive
setDTthreads(8)

STRING_edge <- fread("~/Downloads/9606.protein.physical.links.full.v12.0.onlyAB.txt")

# STRING_edge <- STRING_edge[, c(1,2,10)]
# STRING_info <- fread("~/Downloads/9606.protein.info.v12.0.txt")
# STRING_info <- STRING_info[, c(1,2)]
# colnames(STRING_info) <- c("ENSP", "uniprot")

STRING_aliases <- fread("~/Downloads/9606.protein.aliases.v12.0.txt.gz")
STRING_aliases <- na.omit(STRING_aliases)
# STRING_aliases <- STRING_aliases[alias %in% feature_matrix$Protein]
# STRING_aliases <- STRING_aliases[!duplicated(STRING_aliases$`#string_protein_id`), ]


wide_table <- dcast(
  STRING_aliases,
  `#string_protein_id` ~ source,
  value.var = "alias",
  fun.aggregate = function(x) paste(unique(x), collapse = ";")
)


# this code just reads in original feature matrix to see which alias type to use
# gene_ids <- feature_matrix$Protein
# source_cols <- setdiff(names(wide_table), "#string_protein_id")
# # Count matches for each column
# matching_counts <- sapply(source_cols, function(col_name) {
#   # Count how many entries in the column are in your gene ID set
#   sum(wide_table[[col_name]] %in% gene_ids, na.rm = TRUE)
# })
# 
# # Find the column with the highest number of matches
# names(matching_counts)[which.max(matching_counts)]
# matching_counts[matching_counts != 0]



STRING_uniprot <- STRING_aliases[, .(
  alias = {
    # Identify if any row has source == "Ensembl_UniProt"
    idx <- which(source == "Ensembl_HGNC_uniprot_ids")
    if (length(idx) > 0) {
      alias[idx[1L]]  # Take the first match if multiple
    } else {
      NA_character_    # Return NA if no match
    }
  }
), by = "#string_protein_id"]


# see if id is in a different column
uni_aliases <- STRING_aliases[
  source == "UniProt_AC", 
  .(alias_UniProt_AC = alias[1L]), 
  by = "#string_protein_id"
]
result <- merge(STRING_uniprot, uni_aliases, by = "#string_protein_id", all.x = TRUE)
result[is.na(alias), alias := alias_UniProt_AC]
result[, alias_UniProt_AC := NULL]

uni_aliases <- STRING_aliases[
  source == "Ensembl_UniProt", 
  .(alias_Ensembl_HGNC_uniprot_ids = alias[1L]), 
  by = "#string_protein_id"
]
result <- merge(result, uni_aliases, by = "#string_protein_id", all.x = TRUE)
result[is.na(alias), alias := alias_Ensembl_HGNC_uniprot_ids]
result[, alias_Ensembl_HGNC_uniprot_ids := NULL]

table(unique(c(STRING_edge$protein1, STRING_edge$protein2)) %in% result$`#string_protein_id`)
# table(feature_matrix$Protein %in% result$alias)





STRING_edge_uniprot <- merge.data.table(STRING_edge, result, by.x = "protein1", by.y = "#string_protein_id", all.x = T)
colnames(STRING_edge_uniprot)[11] <- "protein1_uniprot"
STRING_edge_uniprot <- merge.data.table(STRING_edge_uniprot, result, by.x = "protein2", by.y = "#string_protein_id", all.x = T, sort = F)
colnames(STRING_edge_uniprot)[12] <- "protein2_uniprot"
STRING_edge_uniprot <- STRING_edge_uniprot[STRING_edge_uniprot$combined_score >= 400, ]
STRING_edge_uniprot$combined_score <- STRING_edge_uniprot$combined_score / 1000

temp <- STRING_edge_uniprot[!complete.cases(STRING_edge_uniprot), ]
STRING_edge_uniprot <- na.omit(STRING_edge_uniprot)

STRING_edge_uniprot <- subset(STRING_edge_uniprot, select = c("protein1_uniprot", "protein2_uniprot", "combined_score"))


table(unique(c(STRING_edge$protein1, STRING_edge$protein2)) %in% STRING_aliases$`#string_protein_id`)
table(is.na(STRING_aliases$alias))
temp <- STRING_aliases[is.na(STRING_aliases$alias), ]
table(feature_matrix$Protein %in% unique(c(STRING_edge_uniprot$protein1_uniprot, STRING_edge_uniprot$protein2_uniprot)))



library(igraph)
STRING_net <- graph_from_data_frame(STRING_edge_uniprot, directed = F)

string_df <- data.frame(ENSG = V(STRING_net)$name,
                        degree = degree(STRING_net),
                        betweenness = betweenness(STRING_net),
                        closeness = closeness(STRING_net),
                        eigen_centrality = eigen_centrality(STRING_net)$vector,
                        page_rank = page_rank(STRING_net)$vector)

save(string_df, file = "RData/human_string_PPI_metrics.RData")


STRING_top10pct <- STRING_edge_uniprot[order(-combined_score)]
STRING_top10pct <- STRING_top10pct[1:(nrow(STRING_top10pct)*0.1), ]
STRING_top10pct_net <- graph_from_data_frame(STRING_top10pct, directed = F)

STRING_top10pct_df <- data.frame(ENSG = V(STRING_top10pct_net)$name,
                                 degree = degree(STRING_top10pct_net),
                                 betweenness = betweenness(STRING_top10pct_net),
                                 closeness = closeness(STRING_top10pct_net),
                                 eigen_centrality = eigen_centrality(STRING_top10pct_net)$vector,
                                 page_rank = page_rank(STRING_top10pct_net)$vector)

save(STRING_top10pct_df, file = "RData/human_STRING_top10pct_PPI_metrics.RData")
