

# fix up target gene ID mapping
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


genes_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                         filters = "ensembl_gene_id", 
                         values = module_gene_mapping$gene_id, 
                         mart = ensembl)

genes_converted <- merge(module_gene_mapping, genes_converted, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = T)

OpenTargets_mapping <- merge(genes_converted, OpenTargets_raw, by.x = "gene_id", by.y = "Target ID", all.x = T)
OpenTargets_mapping <- subset(OpenTargets_mapping, select = c("gene_id", "Target Approved Symbol", "external_gene_name"))

OpenTargets_mapping$combined <- ifelse(OpenTargets_mapping$external_gene_name == "" | is.na(OpenTargets_mapping$external_gene_name), 
                                       OpenTargets_mapping$`Target Approved Symbol`, OpenTargets_mapping$external_gene_name)
OpenTargets_mapping$combined <- ifelse(is.na(OpenTargets_mapping$combined), OpenTargets_mapping$gene_id, OpenTargets_mapping$combined)
OpenTargets_mapping <- OpenTargets_mapping[!duplicated(OpenTargets_mapping$gene_id) | !duplicated(OpenTargets_mapping$combined), ]


approved_openTargets_EGN <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                                  filters = "ensembl_gene_id", 
                                  values = unique(approved_openTargets$`Target ID`), 
                                  mart = ensembl)
approved_openTargets_EGN <- merge(approved_openTargets, approved_openTargets_EGN, by.x = "Target ID", by.y = "ensembl_gene_id")
approved_openTargets_EGN <- approved_openTargets_EGN[!duplicated(approved_openTargets_EGN$`Target ID`) & !duplicated(approved_openTargets_EGN$external_gene_name), ]



# add gene ensembl to drugbank targets
DrugBank_targets <- read.csv("../BRCA_pipe/DrugBank_targets.csv")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_converted <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                         filters = "external_gene_name", 
                         values = DrugBank_targets$drugBank_target, 
                         mart = ensembl)

genes_converted <- merge(DrugBank_targets, genes_converted, by.x = "drugBank_target", by.y = "external_gene_name", all.x = T)
write.csv(genes_converted, "../BRCA_pipe/DrugBank_targets_ENS.csv", row.names = F)

# check where drugbank targets scored in DE analysis
hits <- DE_results$toptags$table
x <- merge(DrugBank_targets_unique, hits, by.x = "ensembl_gene_id", by.y = "row.names", all.x = T)
table(x$FDR > 0.1) # high FDR genes

# check drug taregts in expanded STRING DE dataset
string_edge_data <- read.table("../../../../Desktop/STRING network (physical) default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("../../../../Desktop/STRING network (physical) default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "display.name"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("display.name.x", "display.name.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(string_edge_data, string_node_data, final_df, file = "BRCA/RData/DE_subset/STRING_DE_subset(+100).RData")

temp <- DrugBank_targets_unique[DrugBank_targets_unique$drugBank_target %in% unique(c(final_df$node_1, final_df$node_2)), ]
DrugBank_targets_unique$drugBank_target[!unique(DrugBank_targets_unique$drugBank_target) %in% unique(temp$drugBank_target)]




# module preservation using seperate expr but combined colours
load("BRCA/RData/all_together/bwnet.RData")
load("../../../../Desktop/WGCNA_BRCA_large_files/data_norm_filt_GTEx.RData") # for ubuntu

tumour_colours <- bwnet$colors[names(bwnet$colors) %in% colnames(tumour_data)]
GTEx_colours <- bwnet$colors[names(bwnet$colors) %in% colnames(control_data)]

x <- tumour_data[, colnames(tumour_data) %in% names(bwnet$colors)]
y <- control_data[, colnames(control_data) %in% names(bwnet$colors)]

multidata <- multiData(Tumour = x,
                       Control = y)
multicolour <- list(Tumour = tumour_colours,
                    Control = GTEx_colours)

library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
preserved_modules <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        dataIsExpr = T,
                                        networkType = "unsigned",
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 100,
                                        referenceNetworks = 1,
                                        maxModuleSize = max(table(bwnet$colors)),
                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time


# plot preserved modules.
plot_preserved_modules <- function(modulePreservation_data) {
  plot_data <- data.frame(
    cluster = rownames(modulePreservation_data$preservation$Z$ref.Tumour$inColumnsAlsoPresentIn.Control),
    moduleSize = modulePreservation_data$preservation$observed$ref.Tumour$inColumnsAlsoPresentIn.Control$moduleSize,
    medianRank.pres = modulePreservation_data$preservation$observed$ref.Tumour$inColumnsAlsoPresentIn.Control$medianRank.pres,
    Zsummary.pres = modulePreservation_data$preservation$Z$ref.Tumour$inColumnsAlsoPresentIn.Control$Zsummary.pres
  )
  
  modColors <- unique(plot_data$cluster) 
  plotData <-  plot_data[, c(2:ncol(plot_data), 1)]
  #plotMods <-  !(modColors %in% c("grey", "gold")) # excludes gold and grey i think
  
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  plot1 <- ggplot(plot_data, aes(x = moduleSize, y = medianRank.pres)) +
    geom_point(aes(fill = factor(cluster)), shape = 21, size = 2.4, colour = modColors) +
    scale_x_log10() +
    labs(x = "Module size", y = "Median Rank") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(aes(label = cluster), position = position_nudge(x = 0.1, y = 0.1), color = "black") +
    geom_hline(yintercept = 8, linetype = "dashed") +
    annotate("text", x = 1.5, y = 8.5, label = "Above", size = 3) +
    scale_fill_manual(values = modColors) 
  
  
  plot2 <- ggplot(plot_data, aes(x = moduleSize, y = Zsummary.pres)) +
    geom_point(aes(fill = factor(cluster)), shape = 21, size = 2.4, colour = modColors) +
    scale_x_log10() +
    labs(x = "Module size", y = "Z Summary") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(aes(label = cluster), position = position_nudge(x = 0.1, y = 0.1), color = "black") +
    geom_hline(yintercept = 10, linetype = "dashed") +
    annotate("text", x = 1.5, y = 9, label = "Below", size = 3) +
    scale_fill_manual(values = modColors)
  
  # Display both plots side by side
  grid.arrange(plot1, plot2, ncol = 2)
  
  return(list(plot_data = list(plot_data = plot_data, modColors = modColors),
              meadianRank_plt = plot1,
              Zsummary_plt = plot2))
}

modulePreservation_plt <- plot_preserved_modules(preserved_modules)


# non-preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres > 8 & plot_data$Zsummary.pres < 10, ]
nonPreservedGenes <- as.data.frame(bwnet$colors)
nonPreservedGenes <- rownames_to_column(nonPreservedGenes)
nonPreservedGenes <- nonPreservedGenes[nonPreservedGenes$`bwnet$colors` %in% non_preserved_modules$cluster, ]

save(preserved_modules, file = "BRCA/RData/all_together/modulePreservation_sepEXPR_combined_colour(n=100).RData")

# get list of genes in each module
gene_modules <- list()
for (colour in non_preserved_modules$cluster) {
  module <- module_gene_mapping[module_gene_mapping$module == colour, ]
  gene_modules[[colour]] <- module
  
  rm(module, colour)
}


# count the number of targets in each module
counts <- c()
for (module in gene_modules) {
  genes <- module$gene_id
  targets_in_genes <- DrugBank_targets_unique[DrugBank_targets_unique$ensembl_gene_id %in% genes, ]
  unique_targets <- unique(targets_in_genes$drugBank_target)
  counts <- append(counts, length(unique_targets))
  
  rm(genes, targets_in_genes, unique_targets)
}

# does not add up to 23 because some gene ensembles = more than one gene symbol
targets_in_modules <- data.frame(module = names(gene_modules),
                                 counts = counts)
rm(counts)



## module preservation using seperate expr but combined colours and on the DE subset data
load("BRCA/RData/all_together/bwnet.RData")
load("../../../../Desktop/WGCNA_BRCA_large_files/DE_subset_data.RData") # for ubuntu
load("../../../../Desktop/WGCNA_BRCA_large_files/data_norm_filt_GTEx.RData") # for ubuntu

tumour_colours <- bwnet$colors[names(bwnet$colors) %in% colnames(tumour_data)]
GTEx_colours <- bwnet$colors[names(bwnet$colors) %in% colnames(control_data)]

tumour_colours <- tumour_colours[names(tumour_colours) %in% DE_genes]
GTEx_colours <- GTEx_colours[names(GTEx_colours) %in% DE_genes]

x <- tumour_data[, colnames(tumour_data) %in% names(bwnet$colors)]
y <- control_data[, colnames(control_data) %in% names(bwnet$colors)]

x <- x[, colnames(x) %in% DE_genes]
y <- y[, colnames(y) %in% DE_genes]

multidata <- multiData(Tumour = x,
                       Control = y)
multicolour <- list(Tumour = tumour_colours,
                    Control = GTEx_colours)

library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
preserved_modules <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        dataIsExpr = T,
                                        networkType = "unsigned",
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 100,
                                        referenceNetworks = 1,
                                        maxModuleSize = max(table(bwnet$colors)),
                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

modulePreservation_plt <- plot_preserved_modules(preserved_modules)

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/all_together/modulePreservation_sepEXPR_combined_colour_DE_subset(n=100).RData")

# non-preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres > 8 & plot_data$Zsummary.pres < 10, ]

nonPreservedGenes <- names(tumour_colours)[tumour_colours %in% non_preserved_modules$cluster]









