### This script combines all the methods used in the different WGCNA attempts to hopefully be the best one ###
## firstly, combines TCGA and GTEx data to perform trait correlation ###
### secondly, analyses the disease groups separately for module preservation ###

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

# add sample type to sample_info
load("../BRCA_pipe/RData/TCGA_query.RData")
common <- common[, c(1,3)]
sample_info <- merge(sample_info, common, by.x = "sample", by.y = "cases", all.x = T)
sample_info$sample_type <- ifelse(is.na(sample_info$sample_type), "Healthy", sample_info$sample_type)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(control_info, lumA_info, lumB_info, her2_info, basal_info)
rm(clinical, common, query_TCGA, subtypes)
collectGarbage()

STN_samples <- sample_info$sample[sample_info$sample_type == "Solid Tissue Normal"]
all_subtypes <- all_subtypes[, !colnames(all_subtypes) %in% STN_samples]

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = GTEx_ENS)

# normalisation (transposes matrix)
all_wgcna_data <- vst_norm(all_subtype_counts_filt)

# plot PCA
PCA_results <- plot_PCA(expr_data = all_wgcna_data,
                        sample_info = sample_info,
                        plot_tree = T,
                        output_plot_data = T)


sft <- pickSoftThreshold(all_wgcna_data,
                         blockSize = 45000,
                         verbose = 2,
                         networkType = "signed")

sft <- sft$fitIndices

library(gridExtra)
library(grid)
a1 <- ggplot(sft, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit') +
  theme_classic()

a2 <- ggplot(sft, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
rm(a1, a2)

save(sft, file = "BRCA/RData/STN_filt/combined_sft.RData")
load("BRCA/RData/STN_filt/combined_sft.RData")

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
bwnet <- blockwiseModules(all_wgcna_data,
                          maxBlockSize = 45000,
                          power = 12,
                          networkType = "signed",
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3,
                          saveTOMs = FALSE)
elapsed_time <- Sys.time() - start_time
print(paste0("Elapsed time: ", elapsed_time))


# Plot the dendrogram
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


save(bwnet, sample_info, file = "BRCA/RData/STN_filt/all_bwnet.RData")
load("BRCA/RData/STN_filt/all_bwnet.RData")


## trait correlation
library(tidyverse)
library(CorLevelPlot)

sample_info <- data.frame(row.names = sample_info$sample[!sample_info$sample %in% STN_samples],
                          status = c(rep("tumour", ncol(all_subtypes)), 
                                     rep("control", ncol(GTEx_ENS))),
                          group = c(rep("lumA", ncol(LumA_unstranded)-12),
                                    rep("lumB", ncol(LumB_unstranded)),
                                    rep("Her2", ncol(Her2_unstranded)),
                                    rep("basal", ncol(Basal_unstranded)),
                                    rep("GTEx", ncol(GTEx_ENS))))

rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
collectGarbage()

traits.state <- binarizeCategoricalColumns.forPlots(sample_info$status)
traits.subtype <- binarizeCategoricalColumns.forPlots(sample_info$group)
traits <- cbind(traits.state, traits.subtype)
rownames(traits) <- c(colnames(all_subtypes), colnames(GTEx_ENS))
traits <- subset(traits, select = c("data.control", "data.tumour", "data.lumA", "data.lumB", "data.Her2","data.basal", "data.GTEx")) # reorder columns

moduleTrait_cor <- cor(bwnet$MEs, traits, use = "p")
moduleTrait_cor_pvals <- corPvalueStudent(moduleTrait_cor, nSamples = nrow(all_wgcna_data))

heatmap_data <- merge(bwnet$MEs, traits, by = "row.names")

heatmap_data <- column_to_rownames(heatmap_data, "Row.names")
colnames(heatmap_data) <- gsub("data.", "", colnames(heatmap_data))

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[13:18],
             y = names(heatmap_data)[1:12],
             col = c("blue1", "skyblue", "white", "pink", "red"))

save(heatmap_data, file = "BRCA/RData/STN_filt/heatmap_data.RData")

# intramodular connectivity
colours <- labels2colors(bwnet$colors)
kWithin <- intramodularConnectivity.fromExpr(all_wgcna_data, colours, power = 6)
rownames(kWithin) <- colnames(all_wgcna_data)
kWithin <- kWithin[order(-kWithin$kWithin), ]
save(kWithin, file = "BRCA/RData/STN_filt/all_kwithin.RData")
load("BRCA/RData/all_default/signed/all_kwithin.RData")

# Does the same thing as above (intra modular connectivity)
module.membership.measure <- cor(bwnet$MEs, all_wgcna_data, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nrow(all_wgcna_data))
module.membership.measure.pvals <- as.data.frame(t(module.membership.measure.pvals))

# gene significance
gene.signf.corr <- cor(all_wgcna_data, traits$data.tumour, use = 'p')
gene.signf.corr.pvals <- as.data.frame(corPvalueStudent(gene.signf.corr, nrow(all_wgcna_data)))
gene.signf.corr.pvals <- rownames_to_column(gene.signf.corr.pvals)
gene.signf.corr.pvals <- gene.signf.corr.pvals[order(gene.signf.corr.pvals$V1), ]
rownames(gene.signf.corr.pvals) <- NULL
colnames(gene.signf.corr.pvals) <- c("gene_id", "pvalue")

# get top 10% genes for connectivity for each non-preserved module
top_connectivity_genes = list()
for (module in unique(bwnet$colors)) {
  moduleGenes = names(bwnet$colors)[bwnet$colors == module]
  moduleKWithin = kWithin[moduleGenes, ]
  topNumGenes <- ceiling(length(moduleGenes) * 0.10)
  topModuleGenes = head(order(moduleKWithin$kWithin, decreasing = TRUE), topNumGenes)
  top_connectivity_genes[[module]] = moduleGenes[topModuleGenes]
  
  rm(moduleGenes, module, topModuleGenes, topNumGenes, moduleKWithin)
}

library(reshape2)
top_connectivity_genes <- melt(top_connectivity_genes)
colnames(top_connectivity_genes) <- c("ensembl_id", "module")
#top_connectivity_genes$DE <- top_connectivity_genes$ensembl_id %in% dif_exp$gene_id


# perform GO and pathway analysis
load("BRCA/RData/all_default/signed/all_bwnet.RData")

library(clusterProfiler)
library(org.Hs.eg.db)
library(progress)
pb <- progress_bar$new(
  format = "  Performing GO Analysis [:bar] :percent eta: :eta",
  total = length(unique(bwnet$colors)), clear = FALSE)

# run GO enrichment
# turn this into a function at some point
all_GO <- list()
for (module in unique(bwnet$colors)) {
  genes <- names(bwnet$colors)[bwnet$colors %in% module]
  GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
  
  all_GO[[module]] <- GO
  
  rm(genes, GO, module)
  pb$tick()
}

save(all_GO, file = "../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/all_WGCNA_GO_BP_signed.RData")


GO_formatted <- data.frame()
for (i in seq_along(all_GO)) {
  module <- all_GO[[i]]
  module_name <- names(all_GO)[i]
  
  result <- module@result
  result_top <- head(result, 5)
  result_top$module <- rep(module_name, nrow(result_top))
  result_top$`-log(p.adjust)` <- -log(result_top$p.adjust)
  
  GO_formatted <- rbind(GO_formatted, result_top)
  
  rm(module, i, result, result_top, module_name)
}

# Function to convert GeneRatio to numeric
convert_gene_ratio <- function(gene_ratio) {
  # Split the string by "/"
  parts <- strsplit(gene_ratio, "/")[[1]]
  # Convert the parts to numeric and calculate the ratio
  ratio <- as.numeric(parts[1]) / as.numeric(parts[2])
  return(ratio)
}

GO_formatted$GeneRatio.num <- sapply(GO_formatted$GeneRatio, convert_gene_ratio)
save(GO_formatted, file = "BRCA/RData/all_default/signed/all_WGCNA_GO_BP.RData")

ggplot(data = GO_formatted, aes(x = module, y = Description, 
                                color = `p.adjust`, size = GeneRatio.num)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  labs(size = "Gene Ratio") +
  ggtitle("GO enrichment analysis (BP)")



## KEGG pwathway
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pb <- progress_bar$new(total = length(unique(bwnet$colors)))

all_KEGG <- list()
for (module in unique(bwnet$colors)) {
  genes <- names(bwnet$colors)[bwnet$colors %in% module]
  genes_converted <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                           filters = "ensembl_gene_id", 
                           values = genes, 
                           mart = ensembl,
                           verbose = F)
  
  KEGG <- enrichKEGG(genes_converted$entrezgene_id, organism = "hsa", keyType = "ncbi-geneid")
  
  all_KEGG[[module]] <- KEGG
  
  rm(genes, genes_converted, KEGG, module)
  pb$tick()
}

KEGG_formatted <- data.frame()
for (i in seq_along(all_KEGG)) {
  module <- all_KEGG[[i]]
  module_name <- names(all_KEGG)[i]
  
  result <- module@result
  result_top <- head(result, 5)
  result_top$module <- rep(module_name, nrow(result_top))
  result_top$`-log(p.adjust)` <- -log(result_top$p.adjust)
  
  KEGG_formatted <- rbind(KEGG_formatted, result_top)
  
  rm(module, i, result, result_top)
}

# uses convert_gene_ratio function
KEGG_formatted$GeneRatio.num <- sapply(KEGG_formatted$GeneRatio, convert_gene_ratio)
save(KEGG_formatted, file = "BRCA/RData/all_default/signed/all_WGCNA_KEGG.RData")

ggplot(data = KEGG_formatted, aes(x = module, y = Description, 
                                  color = `p.adjust`, size = GeneRatio.num)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  labs(size = "Gene Ratio") +
  ggtitle("KEGG enrichment analysis")


## cross section of tumour modules and DE genes
load("../BRCA_pipe/latest_run/RData/STN_filt/dif_exp.RData")

DE_genes_bwnet <- bwnet$colors[names(bwnet$colors) %in% dif_exp$gene_id]
DE_genes_bwnet <- as.data.frame(table(DE_genes_bwnet))
colnames(DE_genes_bwnet) <- c("module", "DE_genes")

temp <- as.data.frame(table(bwnet$colors))
colnames(temp) <- c("module", "total_size")
DE_genes_bwnet <- merge(temp, DE_genes_bwnet, by = "module", all = T)
DE_genes_bwnet[is.na(DE_genes_bwnet)] <- 0
DE_genes_bwnet$`proportion(%)` <- DE_genes_bwnet$DE_genes/DE_genes_bwnet$total_size * 100
DE_genes_bwnet <- DE_genes_bwnet[order(-DE_genes_bwnet$`proportion(%)`), ]


## venn diagram for cross 
DE_genes <- dif_exp$gene_id
tumour_associated <- names(bwnet$colors)[bwnet$colors %in% c("blue", "magenta","green", "brown", "yellow", "pink", "greenyellow")]
top_kwithin <- top_connectivity_genes$ensembl_id
top_gene_membership <- gene.signf.corr.pvals$gene_id[1:(length(gene.signf.corr.pvals$gene_id) * 0.1)]
save(DE_genes, tumour_associated, top_kwithin, top_gene_membership, file = "BRCA/RData/STN_filt/venn_data.RData")

kWithin[rownames(kWithin) %in% "ENSG00000141510", ] # TP53

load("BRCA/RData/all_default/signed/venn_data.RData")

library(ggVennDiagram)
ggVennDiagram(x = list(`DE genes` = DE_genes, 
                       `Tumour Associated` = tumour_associated, 
                       `Top10% kWithin` = top_kwithin,
                       `Top10% MM` = top_gene_membership))



#library(VennDiagram)
#venn.diagram(
#  x = list(`DE genes` = DE_genes, 
#           `Tumour Associated` = tumour_associated, 
#           `Top kWithin` = top_kwithin,
#           `Top Membership Pvalue` = top_gene_membership),
#  category.names = c("DE genes", "Tumour associated", "Top10% Kwithin", "Top10% MM"),
#  col = "transparent",  # set the color of the intersections to transparent
#  fill = c("dodgerblue", "goldenrod1", "lightcoral", "mediumseagreen"),  # set colors for each category
#  alpha = 0.5,  # set the transparency level of the circles
#  cat.col = c("dodgerblue", "goldenrod1", "lightcoral", "mediumseagreen"),  # set colors for category labels
#  cat.fontfamily = "Arial",  # set the font family for category labels
#  cat.fontface = "bold",  # set the font face for category labels
#  cat.fontsize = 10,  # set the font size for category labels
#  cex = 1.5,  # increase the size of the circles
#  margin = 0.1,  # set the margin size (proportion of the plot)
#  filename = "BRCA/RData/all_default/signed/consensus_genes.png",
#  disable.logging = TRUE
#)


common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                         filters = "ensembl_gene_id", 
                         values = common_genes, 
                         mart = ensembl,
                         verbose = F)

GO <- enrichGO(common_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO <- GO@result


targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2:4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

temp <- targets[targets$ensembl_gene_id %in% common_genes, ]

targets$drugBank_target[targets$ensembl_gene_id %in% tumour_associated]

###############################################################################
# Will now perform WGCNA on TCGA and GTEx groups separately for module     
# preservation analysis                                                    
###############################################################################
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(doParallel)
library(reshape2)
library(igraph)

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

# add sample type to sample_info
load("../BRCA_pipe/RData/TCGA_query.RData")
common <- common[, c(1,3)]
sample_info <- merge(sample_info, common, by.x = "sample", by.y = "cases", all.x = T)
sample_info$sample_type <- ifelse(is.na(sample_info$sample_type), "Healthy", sample_info$sample_type)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(control_info, lumA_info, lumB_info, her2_info, basal_info)
rm(clinical, common, query_TCGA, subtypes)
collectGarbage()

STN_samples <- sample_info$sample[sample_info$sample_type == "Solid Tissue Normal"]
all_subtypes <- all_subtypes[, !colnames(all_subtypes) %in% STN_samples]

# QC + combines tumour and control samples
counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                               control_matrix = GTEx_ENS,
                               sep = T)

# normalisation (transposes matrix)
tumour_data <- vst_norm(counts_df = counts_filt$tumour)
control_data <- vst_norm(counts_df = counts_filt$control)



# choose soft thresholding power
pick_power <- function(WGCNA_data) {
  sft <- pickSoftThreshold(WGCNA_data,
                           blockSize = 45000,
                           verbose = 2,
                           networkType = "signed")
  sft <- sft$fitIndices
  
  library(gridExtra)
  library(grid)
  a1 <- ggplot(sft, aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    geom_hline(yintercept = 0.8, color = 'red') +
    labs(x = 'Power', y = 'Scale free topology model fit') +
    theme_classic()
  
  a2 <- ggplot(sft, aes(Power, mean.k., label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    labs(x = 'Power', y = 'Mean Connectivity') +
    theme_classic()
  
  grid.arrange(a1, a2, nrow = 2, top = textGrob(deparse(substitute(WGCNA_data))))
  
  return(sft)
}

tumour_sft <- pick_power(tumour_data)
control_sft <- pick_power(control_data)

save(tumour_sft, file = "BRCA/RData/STN_filt/tumour_sft.RData")
save(control_sft, file = "BRCA/RData/STN_filt/control_sft.RData")


# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 16
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

network_modules <- function(WGCNA_data, Power) {
  start_time <- Sys.time()
  bwnet <- blockwiseModules(WGCNA_data,
                            maxBlockSize = 45000,
                            power = Power,
                            networkType = "signed",
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            randomSeed = 1234,
                            verbose = 3,
                            saveTOMs = FALSE)
  elapsed_time <- Sys.time() - start_time
  cat("Elapsed time: ")
  print(elapsed_time)
  
  # Plot the dendrogram
  plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                      c("unmerged", "merged"),
                      dendroLabels = FALSE,
                      addGuide = TRUE,
                      hang= 0.03,
                      guideHang = 0.05)
  
  return(bwnet)
}

tumour_bwnet <- network_modules(tumour_data, Power = 12)
control_bwnet <- network_modules(control_data, Power = 12)

save(tumour_bwnet, file = "BRCA/RData/STN_filt/tumour_bwnet.RData")
save(control_bwnet, file = "BRCA/RData/STN_filt/control_bwnet.RData")

load("BRCA/RData/STN_filt/tumour_bwnet.RData")
load("BRCA/RData/STN_filt/control_bwnet.RData")


# module preservation analysis
multidata <- multiData(Control = control_data, 
                       Tumour = tumour_data)
multicolour <- list(Control = control_bwnet$colors,
                    Tumour = tumour_bwnet$colors)

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 16
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
preserved_modules <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        dataIsExpr = T,
                                        networkType = "signed",
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 100,
                                        maxModuleSize = max(max(table(tumour_bwnet$colors)), 
                                                            max(table(control_bwnet$colors))),                                        
                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

# plot results
modulePreservation_plt <- plot_preserved_modules(preserved_modules)

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/STN_filt/modulePreservation(n=100).RData")
load("BRCA/RData/STN_filt/modulePreservation(n=100).RData")

library(reshape2)
# non preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres > 8 & plot_data$Zsummary.pres < 10, ]


# non preserved genes common with rest of data: kwithin, MM, DE, tumour associated.
non_preserved_genes <- names(control_bwnet$colors)[control_bwnet$colors %in% non_preserved_modules$cluster]
#non_preserved_common <- non_preserved_genes[non_preserved_genes %in% common_genes]

#library(biomaRt)
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#genes_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
#                         filters = "ensembl_gene_id", 
#                         values = non_preserved_common, 
#                         mart = ensembl,
#                         verbose = F)
#genes_converted$description <- gsub("\\[.*?\\", "", genes_converted$description)

# plot cross-tabulation
cross_tab_counts <- preserved_modules$accuracy$observedCounts$ref.Control$inColumnsAlsoPresentIn.Tumour
cross_tab_pvalues <- preserved_modules$accuracy$observedFisherPvalues$ref.Control$inColumnsAlsoPresentIn.Tumour

data <- melt(cross_tab_counts)
colnames(data) <- c("Test", "Reference", "Count")

pvalues_melted <- melt(cross_tab_pvalues)
colnames(pvalues_melted) <- c("Test", "Reference", "PValue")
# Combine the counts and p-values
data$PValue <- pvalues_melted$PValue
data$LogPValue <- -log10(data$PValue)

# Add the total size of each module
control_sizes <- as.data.frame(table(control_bwnet$colors))
colnames(control_sizes) <- c("Module", "ControlSize")

tumour_sizes <- as.data.frame(table(tumour_bwnet$colors))
colnames(tumour_sizes) <- c("Module", "TumourSize")

# Merge the sizes with the data
data <- merge(data, tumour_sizes, by.x = "Test", by.y = "Module")
data <- merge(data, control_sizes, by.x = "Reference", by.y = "Module")

# Create labels with module size
data$Reference <- paste(data$Reference, "(", data$ControlSize, ")", sep = "")
data$Test <- paste(data$Test, "(", data$TumourSize, ")", sep = "")

ggplot(data, aes(x = Test, y = Reference, fill = LogPValue)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = Count), size = 3) +  # Display only counts
  theme_minimal() +
  labs(title = "Control modules (rows) vs. Tumour modules (columns)",
       x = "Test Modules", y = "Reference Modules") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## plot the colour mapping between control and tumour bwnet
control_colors <- control_bwnet$colors
tumour_colors <- tumour_bwnet$colors

common_genes <- intersect(names(control_colors), names(tumour_colors))

# Identify genes unique to control and tumour datasets
unique_control_genes <- setdiff(names(control_colors), common_genes)
unique_tumour_genes <- setdiff(names(tumour_colors), common_genes)

# Create a data frame mapping control and tumour module colors for common genes
color_mapping <- data.frame(
  gene = common_genes,
  control_color = control_colors[common_genes],
  tumour_color = tumour_colors[common_genes]
)

# Add unique control genes with tumour_colour as "white"
unique_control_df <- data.frame(
  gene = unique_control_genes,
  control_color = control_colors[unique_control_genes],
  tumour_color = "white"
)

# Add unique tumour genes with control_color as "white"
unique_tumour_df <- data.frame(
  gene = unique_tumour_genes,
  control_color = "white",
  tumour_color = tumour_colors[unique_tumour_genes]
)

# Combine all data frames
color_mapping <- bind_rows(color_mapping, unique_control_df, unique_tumour_df)

# Count the number of genes in each control module that correspond to each tumour module color
color_counts <- color_mapping %>%
  group_by(control_color, tumour_color) %>%
  summarise(count = n()) %>%
  ungroup()

# Reshape the data for plotting
color_counts_wide <- color_counts %>%
  spread(key = tumour_color, value = count, fill = 0)

# Gather data back into long format for ggplot
color_counts_long <- color_counts_wide %>%
  gather(key = "tumour_color", value = "count", -control_color)

# Define the color palette
color_palette <- unique(c(color_counts_long$control_color, color_counts_long$tumour_color))
color_palette <- setNames(color_palette, color_palette)
color_palette["white"] <- "white"

# Plot the results without major and minor grid lines and without legend
ggplot(color_counts_long, aes(x = control_color, y = count, fill = tumour_color)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Control Module Color", y = "Number of Genes") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")



## perform enrichment on modules
library(clusterProfiler)
library(org.Hs.eg.db)
library(progress)

# function to enrich with GO terms and format results
GO_analysis <- function(bwnet) {
  pb <- progress_bar$new(total = length(unique(bwnet$colors)))
  all_GO <- list()
  for (module in unique(bwnet$colors)) {
    genes <- names(bwnet$colors)[bwnet$colors %in% module]
    GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
    
    all_GO[[module]] <- GO
    
    rm(genes, GO, module)
    pb$tick()
  }
  
  GO_formatted <- data.frame()
  for (i in seq_along(all_GO)) {
    module <- all_GO[[i]]
    module_name <- names(all_GO)[i]
    
    result <- module@result
    result_top <- head(result, 5)
    result_top$module <- rep(module_name, nrow(result_top))
    result_top$`-log(p.adjust)` <- -log(result_top$p.adjust)
    
    GO_formatted <- rbind(GO_formatted, result_top)
  }
  
  # Function to convert GeneRatio to numeric
  convert_gene_ratio <- function(gene_ratio) {
    # Split the string by "/"
    parts <- strsplit(gene_ratio, "/")[[1]]
    # Convert the parts to numeric and calculate the ratio
    ratio <- as.numeric(parts[1]) / as.numeric(parts[2])
    return(ratio)
  }
  
  GO_formatted$GeneRatio.num <- sapply(GO_formatted$GeneRatio, convert_gene_ratio)
  
  return(GO_formatted)
}


control_GO <- GO_analysis(control_bwnet)
tumour_GO <- GO_analysis(tumour_bwnet)

save(control_GO, tumour_GO, file = "BRCA/RData/all_default/signed/split_GO_data.RData")






temp <- control_GO[control_GO$module %in% "green", ]
temp2 <- tumour_GO[tumour_GO$module %in% "green", ]

common <- intersect(names(control_bwnet$colors)[control_bwnet$colors %in% "green"],  names(tumour_bwnet$colors)[tumour_bwnet$colors %in% "green"])
GO <- enrichGO(common, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

GO <- GO@result


temp <- temp[temp$ONTOLOGY == "BP", ]
temp2 <- temp2[temp2$ONTOLOGY == "BP", ]
temp3 <- temp3[temp3$ONTOLOGY == "BP", ]

temp <- rbind(temp, temp2)
temp <- rbind(temp, temp3)

common_modules <- temp[temp$ID %in% temp2$ID, ]


# master table of all preservation statistics
z.pres <- preserved_modules$preservation$Z$ref.Control$inColumnsAlsoPresentIn.Tumour
z.qual <- preserved_modules$quality$Z$ref.Control$inColumnsAlsoPresentIn.Tumour
m.pres <- preserved_modules$preservation$observed$ref.Control$inColumnsAlsoPresentIn.Tumour
m.qual <- preserved_modules$quality$observed$ref.Control$inColumnsAlsoPresentIn.Tumour

temp <- cbind(z.pres, z.qual, m.pres, m.qual)
temp <- temp[, -c(15,22,36)]


# cross section of tumour modules and DE genes
load("BRCA/RData/DE_subset/dif_exp.RData")

DE_genes_bwnet <- bwnet$colors[names(bwnet$colors) %in% dif_exp$gene_id]
DE_genes_bwnet <- as.data.frame(table(DE_genes_bwnet))
colnames(DE_genes_bwnet) <- c("module", "DE_genes")

temp <- as.data.frame(table(bwnet$colors))
colnames(temp) <- c("module", "total_size")
DE_genes_bwnet <- merge(temp, DE_genes_bwnet, by = "module", all = T)
DE_genes_bwnet[is.na(DE_genes_bwnet)] <- 0
DE_genes_bwnet$`proportion(%)` <- DE_genes_bwnet$DE_genes/DE_genes_bwnet$total_size * 100
DE_genes_bwnet <- DE_genes_bwnet[order(-DE_genes_bwnet$`proportion(%)`), ]



targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2:4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

temp <- as.data.frame(control_bwnet$colors[names(control_bwnet$colors) %in% targets$ensembl_gene_id])
temp2 <- as.data.frame(tumour_bwnet$colors[names(tumour_bwnet$colors) %in% targets$ensembl_gene_id])

temp <- merge(temp, temp2, by = "row.names")
colnames(temp) <- c("gene", "control_bwnet", "tumour_bwnet")
temp <- merge(targets, temp, by.x = "ensembl_gene_id", by.y = "gene")

# matches gene ensembl to GO term
#match_genes_return_id <- function(geneID, temp_genes) {
#  # Split the geneID string by '/'
#  gene_list <- unlist(strsplit(geneID, "/"))
#  # Find the matched gene(s)
#  matched_genes <- gene_list[gene_list %in% temp_genes]
#  # If any match is found, return the matched gene(s), otherwise return NA
#  if (length(matched_genes) > 0) {
#    return(paste(matched_genes, collapse = "/"))  # Collapse multiple matches into a string
#  } else {
#    return(NA)
#  }
#}
#
## apply function to control GO
#matched_genes <- sapply(control_GO_formatted$geneID, match_genes_return_id, temp$ensembl_gene_id)
#temp3 <- control_GO_formatted[!is.na(matched_genes), ]
#temp3$matched_gene <- matched_genes[!is.na(matched_genes)]
#
## apply function to tumour GO
#matched_genes <- sapply(tumour_GO$geneID, match_genes_return_id, temp$ensembl_gene_id)
#temp4 <- tumour_GO[!is.na(matched_genes), ]
#temp4$matched_gene <- matched_genes[!is.na(matched_genes)]




# redo GO to get FULL results for tumour and healthy
pb <- progress_bar$new(
  format = "  Performing GO Analysis [:bar] :percent eta: :eta",
  total = length(unique(control_bwnet$colors)), clear = FALSE)
control_GO <- list()
for (module in unique(control_bwnet$colors)) {
  genes <- names(control_bwnet$colors)[control_bwnet$colors %in% module]
  GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
  
  control_GO[[module]] <- GO
  
  rm(genes, GO, module)
  pb$tick()
}

pb <- progress_bar$new(
  format = "  Performing GO Analysis [:bar] :percent eta: :eta",
  total = length(unique(tumour_bwnet$colors)), clear = FALSE)
tumour_GO <- list()
for (module in unique(tumour_bwnet$colors)) {
  genes <- names(tumour_bwnet$colors)[tumour_bwnet$colors %in% module]
  GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
  
  tumour_GO[[module]] <- GO
  
  rm(genes, GO, module)
  pb$tick()
}



control_GO_formatted <- data.frame()
for (i in seq_along(control_GO)) {
  module <- control_GO[[i]]
  module_name <- names(control_GO)[i]
  
  result <- module@result
  result$module <- rep(module_name, nrow(result))
  result$`-log(p.adjust)` <- -log(result$p.adjust)
  
  control_GO_formatted <- rbind(control_GO_formatted, result)
  
  rm(module, i, result, module_name)
}

tumour_GO_formatted <- data.frame()
for (i in seq_along(tumour_GO)) {
  module <- tumour_GO[[i]]
  module_name <- names(tumour_GO)[i]
  
  result <- module@result
  result$module <- rep(module_name, nrow(result))
  result$`-log(p.adjust)` <- -log(result$p.adjust)
  
  tumour_GO_formatted <- rbind(tumour_GO_formatted, result)
  
  rm(module, i, result, module_name)
}


