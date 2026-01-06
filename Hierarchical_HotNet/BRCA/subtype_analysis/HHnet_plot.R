library(data.table)

STRING_net <- fread("STRING_data/STRING_physical_ENSG.csv")
STRING_net <- STRING_net[!duplicated(t(apply(STRING_net, 1, sort))), ] # get rid of doubled up edges i.e. (a,b) (b,a) same edge weight
STRING_net <- graph_from_data_frame(STRING_net, directed = F)


get_clusters <- function(hh_results, STRING_net, 
                         topx_clusters = 10, 
                         all_viable_clusters = FALSE,   
                         palette = "Set3",     # RColorBrewer palette name or a vector of colours
                         bg_col = "white",
                         plot = "subnet",
                         out_dir) {
  library(igraph)
  
  hh_results <- stringr::str_split(hh_results, pattern = "\t")
  cat("Cluster details:", "\n")
  print(table(lengths(hh_results)))
  
  # option to include all clusters with at least 2 members
  if (isTRUE(all_viable_clusters)) {
    selected <- hh_results[lengths(hh_results) > 1]
  } else {
    selected <- hh_results[seq_len(min(topx_clusters, length(hh_results)))]
  }
  
  # function to create colour palette
  make_cols <- function(n, palette) {
    if (n <= 0) return(character(0))
    
    # user supplied a vector of colours
    if (length(palette) > 1) {
      cols <- palette
      if (length(cols) < n) cols <- grDevices::colorRampPalette(cols)(n)
      return(cols[seq_len(n)])
    }
    
    # palette is a single string -> assume RColorBrewer palette name
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("Package 'RColorBrewer' is required for palette = \"", palette,
           "\". Install it or pass a vector of colours.")
    }
    
    if (!(palette %in% rownames(RColorBrewer::brewer.pal.info))) {
      stop("Unknown RColorBrewer palette: '", palette,
           "'. Try e.g. 'Set3', 'Dark2', 'Paired'.")
    }
    
    max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    base  <- RColorBrewer::brewer.pal(min(max_n, max(3, n)), palette)
    
    if (n <= length(base)) return(base[seq_len(n)])
    grDevices::colorRampPalette(base)(n)
  }
  
  
  cluster_cols <- make_cols(length(selected), palette)
  # Init colours + optional membership attribute
  V(STRING_net)$color <- bg_col
  V(STRING_net)$cluster <- NA_integer_
  
  vnames <- V(STRING_net)$name
  
  for (i in seq_along(selected)) {
    members <- selected[[i]]
    
    idx <- which(vnames %in% members)
    if (length(idx)) {
      V(STRING_net)$color[idx] <- cluster_cols[i]
      V(STRING_net)$cluster[idx] <- i
    }
  }
  
  # subnetwork of cluster genes
  selected_genes <- unique(unlist(selected))
  subnet <- induced_subgraph(STRING_net, vids = V(STRING_net)[name %in% selected_genes])
  
  # add direct neighbors of cluster nodes
  neighbs <- neighborhood(
    STRING_net,
    order = 1,
    nodes = V(STRING_net)[name %in% selected_genes]
  )
  
  neighbs <- purrr::flatten(neighbs)
  neighbs <- names(neighbs)
  neighbs <- unique(neighbs)
  
  subnet_neighs <- induced_subgraph(STRING_net, vids = neighbs)
  #subnet_neighs <- induced_subgraph(STRING_net, V(STRING_net)[V(STRING_net)$name %in% neighbs])
  V(subnet_neighs)$size <- scales::rescale(degree(subnet_neighs), c(1, 7))
  
  if (plot == "subnet") {
    set.seed(1234)
    plot.igraph(subnet, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)
  } else if (plot == "subnet_neighs") {
    set.seed(1234)
    plot.igraph(subnet_neighs, asp = 0, vertex.label = NA, edge.arrow.size = 0.3, edge.curved = F)
  } else if (plot == "both") {
    set.seed(1234)
    plot.igraph(subnet, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)
    plot.igraph(subnet_neighs, asp = 0, vertex.label = NA, edge.arrow.size = 0.3, edge.curved = F)
  } else if (plot == "none") {
    
  } else {
    cat("Invalid plot parameter.","\n", 
        "Availiable options are: 'subnet', 'subnet_neighs', 'both', 'none'")
  } 
  
  # format network metrics for both subnetworks
  df_subnet <- data.frame(ENSG = V(subnet)$name,
                          degree = degree(subnet),
                          betweenness = betweenness(subnet),
                          closeness = closeness(subnet),
                          eigen_centrality = eigen_centrality(subnet)$vector,
                          page_rank = page_rank(subnet)$vector,
                          cluster = V(subnet)$color)
  df_subnet$source <- ifelse(df_subnet$cluster != bg_col, "subnet", "STRING")
  
  
  df_subnetNeighs <- data.frame(ENSG = V(subnet_neighs)$name,
                                degree = degree(subnet_neighs),
                                betweenness = betweenness(subnet_neighs),
                                closeness = closeness(subnet_neighs),
                                eigen_centrality = eigen_centrality(subnet_neighs)$vector,
                                page_rank = page_rank(subnet_neighs)$vector,
                                cluster = V(subnet_neighs)$color)
  df_subnetNeighs$source <- ifelse(df_subnetNeighs$cluster != bg_col, "subnet", "STRING")
  
  # get id converter fucntion and run for both subnetworks
  source("../MOC_pipe/R/functions.R")
  cat("\n", "Annotating with gene symbols...", "\n")
  
  if (exists("ensembl")) {
    df_subnet <- id_annot(ensembl, data = df_subnet, 
                          input_type = "ensembl_gene_id", 
                          convert_to = c("external_gene_name", "description", "gene_biotype"))
    
    df_subnetNeighs <- id_annot(ensembl, data = df_subnetNeighs, 
                                input_type = "ensembl_gene_id", 
                                convert_to = c("external_gene_name", "description", "gene_biotype"))
  } else {
    df_subnet <- id_annot(ensembl, data = df_subnet, 
                          input_type = "ensembl_gene_id", 
                          convert_to = c("external_gene_name", "description", "gene_biotype"))
    
    df_subnetNeighs <- id_annot(data = df_subnetNeighs, 
                                input_type = "ensembl_gene_id", 
                                convert_to = c("external_gene_name", "description", "gene_biotype"))
  }
  

  write_graph(subnet_neighs, file = file.path(out_dir, "hhnet_subnetNeighs.graphml"), format = "graphml")
  write_graph(subnet, file = file.path(out_dir, "hhnet_subnet.graphml"), format = "graphml")
  
  write.csv(df_subnet, file = file.path(out_dir, "df_subnet.csv"), row.names = FALSE)
  write.csv(df_subnetNeighs, file = file.path(out_dir, "df_subnet_neighs.csv"), row.names = FALSE)
  
  #return items
  return(list(
    clusters = selected,
    colours = cluster_cols,
    subnet = subnet,
    subnet_neighs = subnet_neighs,
    df_subnet = df_subnet,
    df_subnetNeighs = df_subnetNeighs
  ))
}


# repeat for each subtype
hh_results <- readr::read_lines("BRCA/subtype_analysis/lumA/results/clusters_STRING_lumA_logFC.tsv", skip = 7)
hh_results <- gsub('[\\\\"]', '', hh_results) # clean up escaped quotes/backslashes. 
lumA_results <- get_clusters(hh_results, STRING_net, topx_clusters = 10,
                    palette = c("tomato", "springgreen", "royalblue", "maroon1", "gold",
                                "orchid", "cyan", "yellowgreen", "mediumseagreen", "saddlebrown"),
                    out_dir = "BRCA/subtype_analysis/lumA/results")

hh_results <- readr::read_lines("BRCA/subtype_analysis/lumB/results/clusters_STRING_lumB_logFC.tsv", skip = 7)
hh_results <- gsub('[\\\\"]', '', hh_results)
lumB_results <- get_clusters(hh_results, STRING_net, topx_clusters = 10,
                    palette = c("tomato", "springgreen", "royalblue", "maroon1", "gold",
                                "orchid", "cyan", "yellowgreen", "mediumseagreen", "saddlebrown"),
                    out_dir = "BRCA/subtype_analysis/lumB/results")

hh_results <- readr::read_lines("BRCA/subtype_analysis/Her2/results/clusters_STRING_Her2_logFC.tsv", skip = 7)
hh_results <- gsub('[\\\\"]', '', hh_results)
Her2_results <- get_clusters(hh_results, STRING_net, topx_clusters = 10,
                    palette = c("tomato", "springgreen", "royalblue", "maroon1", "gold",
                                "orchid", "cyan", "yellowgreen", "mediumseagreen", "saddlebrown"),
                    out_dir = "BRCA/subtype_analysis/Her2/results")

hh_results <- readr::read_lines("BRCA/subtype_analysis/basal/results/clusters_STRING_basal_logFC.tsv", skip = 7)
hh_results <- gsub('[\\\\"]', '', hh_results)
basal_results <- get_clusters(hh_results, STRING_net, topx_clusters = 10,
                    palette = c("tomato", "springgreen", "royalblue", "maroon1", "gold",
                                "orchid", "cyan", "yellowgreen", "mediumseagreen", "saddlebrown"),
                    out_dir = "BRCA/subtype_analysis/basal/results")


save(lumA_results, lumB_results, Her2_results, basal_results, file = "BRCA/subtype_analysis/subtype_results.RData")








df_subnet <- read.csv("MOC/results/df_subnet.csv")
df_subnetNeighs <- read.csv("MOC/results/df_subnetNeighs.csv")


library(clusterProfiler)
library(progress)

# MOCvsBEN <- read.csv("../MOC_pipe/results/HHnet_RS_ML_overlap.csv")
# genes_ENS <- id_annot(ensembl, data = MOCvsBEN$external_gene_name, input_type = "external_gene_name", convert_to = "ensembl_gene_id")
# 
# GO <- enrichGO(MOCvsBEN$external_gene_name, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP", universe = rownames(filtered_MOC_data))
# result <- GO@result

load("../MOC_pipe/DE/MOC_vs_BEN/DE_results.RData")
filtered_MOC_data <- DE_results$input_data

pb <- progress_bar$new(
  format = "  Performing GO Analysis [:bar] :percent eta: :eta",
  total = 10
)

cluster_GO <- list()
clust_no <- 1
for (cluster in hh_results[1:10]) {
  GO <- enrichGO(cluster, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP", universe = rownames(filtered_MOC_data))
  cluster_GO[[clust_no]] <- GO
  clust_no <- clust_no + 1
  pb$tick()
}
rm(clust_no, pb, cluster, GO)

GO_formatted <- data.frame()
clust_no <- 1
for (cluster in cluster_GO) {
  result <- cluster@result
  result_top <- head(result, 10)
  result_top$cluster <- rep(clust_no, nrow(result_top))
  
  GO_formatted <- rbind(GO_formatted, result_top)
  clust_no <- clust_no + 1
}
rm(clust_no, result_top, result)

GO_formatted$colour <- ifelse(GO_formatted$cluster == 1, "tomato", "")
GO_formatted$colour <- ifelse(GO_formatted$cluster == 2, "springgreen", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 3, "royalblue", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 4, "maroon1", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 5, "gold", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 6, "orchid", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 7, "cyan", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 8, "yellowgreen", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 9, "mediumseagreen", GO_formatted$colour)
GO_formatted$colour <- ifelse(GO_formatted$cluster == 10, "saddlebrown", GO_formatted$colour)

save(cluster_GO, GO_formatted, file = "~/OneDrive - RMIT University/PhD/large_git_files/MOC/MOCvsBEN_HHnet_cluster_GO.RData")
fwrite(GO_formatted, "MOC/results/subnet_cluster_GO.csv")

load("~/OneDrive - RMIT University/PhD/large_git_files/HHnet/HHnet_cluster_GO.RData")


# remove redundant GO terms based on semantic similarity
library(rrvgo)

cluster_GO_reduced <- list()
clust_no <- 1
for (cluster in cluster_GO) {
  cat("Processing cluster ", clust_no, "...", sep = "")
  
  if (length(cluster@result$ID) > 100) {
    GO_ids <- cluster@result$ID[1:100]
    qvals <- cluster@result$qvalue[1:100]
  } else {
    GO_ids <- cluster@result$ID
    qvals <- cluster@result$qvalue
  }
  
  simMatrix <- calculateSimMatrix(GO_ids,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(qvals), GO_ids)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  results <- list(reducedTerms = reducedTerms, simMatrix = simMatrix)
  
  cluster_GO_reduced[[clust_no]] <- results
  clust_no <- clust_no + 1
}

rm(simMatrix, scores, reducedTerms, cluster, clust_no, results, GO_ids, qvals)
gc()


GO_formatted_reduced <- data.frame()
clust_no <- 1
for (cluster in cluster_GO_reduced) {
  result <- cluster$reducedTerms
  result_top <- head(result, 10)
  result_top$cluster <- rep(clust_no, nrow(result_top))
  
  GO_formatted_reduced <- rbind(GO_formatted_reduced, result_top)
  clust_no <- clust_no + 1
}
rm(clust_no, result_top, result)
gc()

GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 1, "tomato", "")
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 2, "springgreen", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 3, "royalblue", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 4, "maroon1", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 5, "gold", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 6, "orchid", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 7, "cyan", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 8, "yellowgreen", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 9, "mediumseagreen", GO_formatted$colour)
GO_formatted_reduced$colour <- ifelse(GO_formatted$cluster == 10, "saddlebrown", GO_formatted$colour)


# plot results
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

scatterPlot(simMatrix, reducedTerms)

wordcloudPlot(reducedTerms, min.freq=1, colors="black")



