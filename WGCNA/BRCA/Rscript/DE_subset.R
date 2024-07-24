library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(doParallel)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

### subset DE genes from network
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
rm(lumA_info, lumB_info, her2_info, basal_info)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded)

# read in functions from "../BRCA_pipe/Rscript/DE_functions.R"
DE_counts_filt <- filter_low_expr(disease_data = all_subtypes, 
                               control_data = GTEx_ENS)

DE_results <- DE_analysis(counts_matrix = DE_counts_filt$counts_filt,
                          sample_info = DE_counts_filt$sample_info)

# save to onedrive (macOS only)
save(DE_results, DE_counts_filt, DE_genes, file = "../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/DE_subset_data.RData")

# check to see if drugs are in DE dataset
OpenTargets <- read.csv("../BRCA_pipe/OpenTargets_data/OpenTargets_unique_drug.csv", row.names = 1)
OpenTargets_raw <- read_tsv("../BRCA_pipe/OpenTargets_data/breast_carcinoma_known_drugs.tsv")
OpenTargets_updated <- read.csv("../BRCA_pipe/OpenTargets_data/OpenTargets_updated.csv")

NIH_targets <- read.table("../BRCA_pipe/NIH_BRCA_approved_drugs.txt", sep = "\t")
colnames(NIH_targets)[1] <- "approved_drugs"
NIH_targets$approved_drugs <- toupper(NIH_targets$approved_drugs)

approved_openTargets <- merge(NIH_targets, OpenTargets_raw, by.x = "approved_drugs", by.y = "Drug Name")
approved_openTargets <- approved_openTargets[!duplicated(approved_openTargets$approved_drugs) | !duplicated(approved_openTargets$`Target ID`), ]

# number of targets in DE dataset
#temp <- merge(approved_openTargets, DE_results$dif_exp, by.x = "Target ID", by.y = "gene_id", all.x = T)
#temp <- temp[!duplicated(temp$`Target Approved Symbol`), ]
#table(is.na(temp$logFC))

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% DE_results$dif_exp$gene_id, ]
table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(temp$`Target Approved Symbol`))

# need to re-filter on only tumour set. Make sure to read in WGCNA_functions.R again
wgcna_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                     control_matrix = GTEx_ENS,
                                     sep = T)

# VST normalisation
tumour_data_norm <- vst_norm(counts_df = wgcna_counts_filt$tumour)

# subset DE genes from WGCNA dataset
DE_genes <- DE_results$dif_exp$gene_id
tumour_DE_subset <- tumour_data_norm[, colnames(tumour_data_norm) %in% DE_genes]

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% colnames(tumour_DE_subset), ]
table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(temp$`Target Approved Symbol`))

# PCA
PCA_tumour <- plot_PCA(expr_data = tumour_DE_subset,
                       sample_info = tumour_info,
                       plot_tree = T,
                       output_plot_data = T)

# choose soft thresholding power
tumour_sft_data <- pick_power(WGCNA_data = tumour_DE_subset,
                              network_type = "unsigned")

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

tumour_bwnet <- network_modules(WGCNA_data = tumour_DE_subset,
                                Power = 6)

save(tumour_bwnet, file = "BRCA/RData/DE_subset/tumour_bwnet.RData")

# intramodular connectivity
colours <- labels2colors(tumour_bwnet$colors)
tumour_kWithin <- intramodularConnectivity.fromExpr(tumour_DE_subset, colours, power = 6)
rownames(tumour_kWithin) <- colnames(tumour_DE_subset)
tumour_kWithin <- tumour_kWithin[order(-tumour_kWithin$kWithin), ]

save(tumour_kWithin, file = "BRCA/RData/DE_subset/tumour_kWithin.RData")

# get top 10 genes for connectivity for each non-preserved module
tumour_topGenes = list()
for (module in unique(tumour_bwnet$colors)) {
  moduleGenes = names(tumour_bwnet$colors)[tumour_bwnet$colors == module]
  moduleKWithin = tumour_kWithin[moduleGenes, ]
  topModuleGenes = head(order(moduleKWithin$kWithin, decreasing = TRUE), 3)
  tumour_topGenes[[module]] = moduleGenes[topModuleGenes]
  
  rm(moduleGenes, module, topModuleGenes)
}

tumour_topGenes <- melt(tumour_topGenes)
colnames(tumour_topGenes) <- c("ensembl_id", "module")

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% tumour_topGenes$ensembl_id, ]

temp <- merge(tumour_kWithin, approved_openTargets, by.x = "row.names", by.y = "Target ID")
temp <- temp[!duplicated(temp$Row.names), ]



topHubGenes <- chooseTopHubInEachModule(tumour_DE_subset, colorh = tumour_bwnet$colors, power = 6, type = "unsigned", omitColors = "grey")








### analyse each cluster individually ###
module_gene_mapping <- as.data.frame(tumour_bwnet$colors)
module_gene_mapping <- rownames_to_column(module_gene_mapping)
colnames(module_gene_mapping) <- c("gene_id", "module")

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% module_gene_mapping$gene_id, ]
table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(temp$`Target Approved Symbol`))

# try converting to external_gene_name first to try and avoid conflict between Ensembl IDs with OpenTargets
# OpenTargets bloke was smoking meth when mapping ensembl IDs to gene symbols
OpenTargets <- read.csv("../BRCA_pipe/OpenTargets_data/OpenTargets_unique_drug.csv", row.names = 1)
OpenTargets_raw <- read_tsv("../BRCA_pipe/OpenTargets_data/breast_carcinoma_known_drugs.tsv")
OpenTargets_updated <- read.csv("../BRCA_pipe/OpenTargets_data/OpenTargets_updated.csv")

NIH_targets <- read.table("../BRCA_pipe/NIH_BRCA_approved_drugs.txt", sep = "\t")
colnames(NIH_targets)[1] <- "approved_drugs"
NIH_targets$approved_drugs <- toupper(NIH_targets$approved_drugs)

approved_openTargets <- merge(NIH_targets, OpenTargets_raw, by.x = "approved_drugs", by.y = "Drug Name")
approved_openTargets <- approved_openTargets[!duplicated(approved_openTargets$approved_drugs) | !duplicated(approved_openTargets$`Target ID`), ]


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



# get list of genes in each module
gene_modules <- list()
for (colour in unique(tumour_bwnet$colors)) {
  module <- module_gene_mapping[module_gene_mapping$module == colour, ]
  gene_modules[[colour]] <- module
  
  rm(module, colour)
}

# map to gene symbols
gene_modules_converted <- list()
for (module in gene_modules) {
  module <- merge(module, OpenTargets_mapping, all.x = T)
  module_genes <- module$combined
  
  colour <- unique(module$module)
  gene_modules_converted[[colour]] <- module_genes
  
  rm(colour, module_genes)
}

## NEED TO FIX OPENTARGETS DATASET FIRST ##
# count the number of targets in each module
counts <- c()
for (module in names(gene_modules_converted)) {
  genes <- gene_modules_converted[[colour]]
  
  targets_in_genes <- approved_openTargets[approved_openTargets$`Target Approved Symbol` %in% genes, ]
  targets_in_ensemblGenes <- approved_openTargets[approved_openTargets$`Target ID` %in% genes, ]
  targets_in_genes <- rbind(targets_in_genes, targets_in_ensemblGenes)
  
  unique_targets <- unique(targets_in_genes$`Target Approved Symbol`)
  counts <- append(counts, length(unique_targets))
  
  rm(genes, unique_targets, targets_in_genes, colour)
}

# does not add up to 23 because some gene ensembles = more than one gene symbol
targets_in_modules <- data.frame(module = names(gene_modules),
                                 counts = counts)
rm(counts)




### preserved modules ###
# module preservation using expr data
# need to filter and normalise GTEx data
GTEx_norm <- vst_norm(wgcna_counts_filt$control)

multidata <- multiData(Tumour = tumour_DE_subset,
                       Control = GTEx_norm)
multicolour <- list(Tumour = tumour_bwnet$colors)

# RESTART R AND LOAD WGCNA ONLY
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
                                        maxModuleSize = max(table(tumour_bwnet$colors)),
                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time


# plot preserved modules
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

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/DE_subset/modulePreservation(n=100).RData")

# non-preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres > 8 & plot_data$Zsummary.pres < 10, ]

nonPreservedGenes <- common_genes[tumour_common_colours %in% non_preserved_modules$cluster]






#### PCSF with DE subset and weighted with druggability scores (using new druggability score formulae) ####

