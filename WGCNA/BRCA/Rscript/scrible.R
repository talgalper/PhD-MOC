

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



### module preservation using seperate expr but combined colours and on the DE subset data ###
load("../../../../Desktop/WGCNA_BRCA_large_files/DE_subset_data.RData") # for ubuntu

load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_normal.RData")

# load data
GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})

# subtype sample info
control_info <- data.frame(sample = colnames(GTEx_ENS),
                           group = rep("control", ncol(GTEx_ENS)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("basal", ncol(Basal_unstranded)))
tumour_info <- rbind(lumA_info, lumB_info, her2_info, basal_info)
sample_info <- rbind(control_info, tumour_info)
rm(lumA_info, lumB_info, her2_info, basal_info, control_info, tumour_info)


# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded)


# need to re-filter on only tumour set. Make sure to read in WGCNA_functions.R again
wgcna_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                     control_matrix = GTEx_ENS,
                                     sep = F)
# VST normalisation
wgcna_data_norm <- vst_norm(counts_df = wgcna_counts_filt)

# subset DE genes from WGCNA dataset
DE_genes <- DE_results$dif_exp$gene_id
DE_subset <- wgcna_data_norm[, colnames(wgcna_data_norm) %in% DE_genes]

# check for targets using new drugbank targets
DrugBank_targets <- read.csv("../BRCA_pipe/DrugBank_targets_ENS.csv")
DrugBank_targets_unique <- na.omit(DrugBank_targets)
DrugBank_targets_unique <- DrugBank_targets_unique[!duplicated(DrugBank_targets_unique$ensembl_gene_id), ]

temp <- DrugBank_targets_unique[DrugBank_targets_unique$ensembl_gene_id %in% colnames(DE_subset), ]
table(unique(DrugBank_targets_unique$drugBank_target) %in% unique(temp$drugBank_target))

# PCA
PCA_all <- plot_PCA(expr_data = DE_subset,
                    sample_info = sample_info,
                    plot_tree = F,
                    output_plot_data = T)

# choose soft thresholding power
sft_data <- pick_power(WGCNA_data = DE_subset,
                       network_type = "unsigned")

sft_data <- pick_power(WGCNA_data = DE_subset,
                       network_type = "signed")


# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

DE_all_bwnet <- network_modules(WGCNA_data = DE_subset,
                                Power = 6)

save(DE_all_bwnet, file = "BRCA/RData/DE_subset/DE_all_bwnet.RData")

# load in individually filtered and normalised datasets
load("../../../../Desktop/WGCNA_BRCA_large_files/data_norm_filt_GTEx.RData") # for ubuntu

tumour_colours <- DE_all_bwnet$colors[names(DE_all_bwnet$colors) %in% colnames(tumour_data)]
GTEx_colours <- DE_all_bwnet$colors[names(DE_all_bwnet$colors) %in% colnames(control_data)]

#tumour_colours <- tumour_colours[names(tumour_colours) %in% DE_genes]
#GTEx_colours <- GTEx_colours[names(GTEx_colours) %in% DE_genes]

tumour_data <- tumour_data[, colnames(tumour_data) %in% names(DE_all_bwnet$colors)]
control_data <- control_data[, colnames(control_data) %in% names(DE_all_bwnet$colors)]

#x <- x[, colnames(x) %in% DE_genes]
#y <- y[, colnames(y) %in% DE_genes]

multidata <- multiData(Control = control_data,
                       Tumour = tumour_data)
multicolour <- list(Control = GTEx_colours,
                    Tumour = tumour_colours)

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
                                        maxModuleSize = max(table(DE_all_bwnet$colors)),
                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

modulePreservation_plt <- plot_preserved_modules(preserved_modules)

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/DE_subset/modulePreservation_DE_subset_combined(n=100).RData")

# non-preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres > 8 & plot_data$Zsummary.pres < 10, ]

nonPreservedGenes <- names(tumour_colours)[tumour_colours %in% non_preserved_modules$cluster]

temp <- DrugBank_targets_unique[DrugBank_targets_unique$ensembl_gene_id %in% nonPreservedGenes, ]
table(unique(DrugBank_targets_unique$drugBank_target) %in% unique(temp$drugBank_target))


###########################################################


# XENA DE pipeline data GTEx vs TCGA
xena_DE <- read.csv("../../../../Downloads/DEG_results_GTEX vs. TCGA.csv")
xena_DE_targets <- merge(DrugBank_targets_unique, xena_DE, by.x = "drugBank_target", by.y = "X", all.x = T)
xena_DE_targets <- xena_DE_targets[!duplicated(xena_DE_targets$drugBank_target), ]


xena_DE <- read.csv("../../../../Downloads/DEG_results_GTEX vs. TCGA.csv")
xena_DE <- xena_DE[xena_DE$logFC >= 1 | xena_DE$logFC <= -1, ]
xena_DE <- xena_DE[xena_DE$adj.P.Val <= 0.05, ]

xena_DE <- subset(xena_DE, select = c("X", "logFC"))
colnames(xena_DE) <- c("xena_gene", "xena_logFC")

DE_hits <- DE_results$hits
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                         filters = "ensembl_gene_id", 
                         values = DE_hits$gene_id, 
                         mart = ensembl)

genes_converted <- merge(DE_hits, genes_converted, by = "gene_id", all.x = T)

table(xena_DE$xena_gene %in% unique(genes_converted$external_gene_name))





# subtype sample info
control_info <- data.frame(sample = colnames(GTEx_ENS),
                           group = rep("control", ncol(GTEx_ENS)))
normal_info <- data.frame(sample = colnames(normal_unstranded),
                           group = rep("normal", ncol(normal_unstranded)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("basal", ncol(Basal_unstranded)))
sample_info <- rbind(control_info, normal_info, lumA_info, lumB_info, her2_info, basal_info)
rm(lumA_info, lumB_info, her2_info, basal_info, control_info, normal_info)


# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
control <- merge(GTEx_ENS, normal_unstranded, by = "row.names")
control <- column_to_rownames(control, "Row.names")

rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded)


# need to re-filter on only tumour set. Make sure to read in WGCNA_functions.R again
wgcna_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                     control_matrix = control,
                                     sep = F)
# VST normalisation
wgcna_data_norm <- vst_norm(counts_df = wgcna_counts_filt)

all_data_PCA <- plot_PCA(wgcna_data_norm, sample_info = sample_info, plot_tree = F, output_plot_data = T)








# module preservation analysis with only reference colour module. Ends up the same.
multidata <- multiData(Control = control_data, 
                       Tumour = tumour_data)
multicolour <- list(Control = control_bwnet$colors)

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
preserved_modules_2 <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        dataIsExpr = T,
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 100,
                                        maxModuleSize = max(max(table(tumour_bwnet$colors)), 
                                                            max(table(control_bwnet$colors))),                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

# plot results
modulePreservation_plt_2 <- plot_preserved_modules(preserved_modules_2)

plot_data_2 <- modulePreservation_plt_2$plot_data$plot_data
non_preserved_modules_2 <- plot_data_2[plot_data_2$medianRank.pres > 8 & plot_data_2$Zsummary.pres < 10, ]



# non preserved modules with shared colour
multidata <- multiData(Control = control_data, 
                       Tumour = tumour_data)
multicolour <- list(Control = bwnet$colors)

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
preserved_modules_3 <- modulePreservation(multiData = multidata,
                                          multiColor = multicolour,
                                          dataIsExpr = T,
                                          quickCor = 1,
                                          randomSeed = 1234,
                                          verbose = 3,
                                          nPermutations = 100,
                                          maxModuleSize = max(table(bwnet$colors)), 
                                          calculateClusterCoeff = F,
                                          parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

# plot results
modulePreservation_plt_3 <- plot_preserved_modules(preserved_modules_3)

# non preserved modules
plot_data_3 <- modulePreservation_plt_3$plot_data$plot_data
non_preserved_modules_3 <- plot_data_3[plot_data_3$medianRank.pres > 8 & plot_data_3$Zsummary.pres < 10, ]




# run GO enrichment for all ontologies BP, CC & MF
all_GO <- list()
for (module in unique(bwnet$colors)) {
  genes <- names(bwnet$colors)[bwnet$colors %in% module]
  GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")
  result <- GO@result
  split_result <- split(result, result$ONTOLOGY)
  result_first_5 <- lapply(split_result, function(df) head(df, 5))
  result_first_5 <- do.call(rbind, result_first_5)
  
  all_GO[[module]] <- result_first_5
  
  rm(genes, GO, split_result, result_first_5, module)
  pb$tick()
}

save(all_GO, file = "BRCA/RData/all_default/all_WGCNA_GO.RData")


temp <- all_GO$salmon
temp$`log(p.adjust)` <- -log(temp$p.adjust)

ggplot(temp, aes(x = reorder(Description, ONTOLOGY, FUN = identity), y = `log(p.adjust)`, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  ggtitle("Module: black") +
  labs(x = "", y = "p.adjust", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
        plot.margin = margin(l = 50, r = 10, t = 10, b = 10), 
        panel.grid = element_blank())







#### WGCNA default run using DE subset
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(doParallel)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()


# load in data
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_normal.RData")

# load data
GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})

# sample info
control_info <- data.frame(sample = colnames(GTEx_ENS),
                           group = rep("control", ncol(GTEx_ENS)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("basal", ncol(Basal_unstranded)))
sample_info <- rbind(control_info, lumA_info, lumB_info, her2_info, basal_info)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(control_info, lumA_info, lumB_info, her2_info, basal_info)
collectGarbage()

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = GTEx_ENS)

# normalisation (transposes matrix)
all_wgcna_data <- vst_norm(all_subtype_counts_filt)

# plot PCA
PCA_results <- plot_PCA(expr_data = all_wgcna_data,
                        sample_info = sample_info,
                        plot_tree = F,
                        output_plot_data = T)

# plot PCA of DE subset
load("BRCA/RData/DE_subset/dif_exp.RData")
wgcna_data_DE_subset <- all_wgcna_data[, colnames(all_wgcna_data) %in% dif_exp$gene_id]

PCA_results_DE_subset <- plot_PCA(expr_data = wgcna_data_DE_subset,
                        sample_info = sample_info,
                        plot_tree = F,
                        output_plot_data = T)

sft_DE_subset <- pickSoftThreshold(wgcna_data_DE_subset,
                         blockSize = 45000,
                         verbose = 2)
sft_DE_subset <- sft_DE_subset$fitIndices

library(gridExtra)
library(grid)
a1 <- ggplot(sft_DE_subset, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit') +
  theme_classic()

a2 <- ggplot(sft_DE_subset, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
rm(a1, a2)

save(sft, file = "BRCA/RData/all_default/combined_sft.RData")


# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
bwnet_DE_subset <- blockwiseModules(wgcna_data_DE_subset,
                                    maxBlockSize = 45000,
                                    power = 6,
                                    mergeCutHeight = 0.25,
                                    numericLabels = FALSE,
                                    randomSeed = 1234,
                                    verbose = 3,
                                    saveTOMs = FALSE)
elapsed_time <- Sys.time() - start_time
print(paste0("Elapsed time: ", elapsed_time))


# Plot the dendrogram
plotDendroAndColors(bwnet_DE_subset$dendrograms[[1]], cbind(bwnet_DE_subset$unmergedColors, bwnet_DE_subset$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


save(bwnet, sample_info, file = "BRCA/RData/all_default/all_bwnet.RData")
load("BRCA/RData/all_default/all_bwnet.RData")






