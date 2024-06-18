library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(gridExtra)
library(doParallel)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
allowWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

load("../BRCA_pipe/RData/paired/paired_subtypes.RData")

subtype_info <- as.data.frame(lumA_paired)
subtype_info <- subtype_info[!duplicated(subtype_info[1]), ] # check for duplicates
colnames(subtype_info) <- sub(".*\\.", "", colnames(subtype_info))

# load in data
load("../BRCA_pipe/RData/TCGA_normal.RData")
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_query.RData")

# merge normal and disease samples
data <- merge(LumA_unstranded, normal_unstranded, by = "row.names")
data <- column_to_rownames(data, var = "Row.names")

data <- data[, colnames(data) %in% c(subtype_info$normal, subtype_info$lumA)]

# remove outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# filter out bad genes. There were no bad samples
data <- data[gsg$goodGenes == TRUE, ]


sample_info <- data.frame(samples = c(subtype_info$normal, subtype_info[[3]]),
                          pat_id = c(subtype_info$bcr_patient_barcode, subtype_info$bcr_patient_barcode),
                          group = c(rep("control", length(subtype_info$normal)), rep("disease", length(subtype_info[[3]]))))
sample_info$group <- factor(sample_info$group, levels = c("control", "disease"))
sample_info$pat_id <- factor(sample_info$pat_id)

counts_filt <- filterByExpr(data, group = sample_info$group)

# passed genes
table(counts_filt)
counts_filt <- data[counts_filt, ]

# normalisation
wgcna_data <- as.matrix(counts_filt)
wgcna_data <- varianceStabilizingTransformation(wgcna_data)
wgcna_data <- as.data.frame(wgcna_data)
wgcna_data <- t(wgcna_data)

# Choose a set of soft-threshold powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(wgcna_data,
                         powerVector = power,
                         networkType = "unsigned",
                         verbose = 5)

sft_data <- sft$fitIndices


# Visualise to pick power
a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

set.seed(1234)
start_time <- Sys.time()
# identify modules. this includes both benign and disease groups.
bwnet <- blockwiseModules(wgcna_data,
                          maxBlockSize = 45000,
                          TOMType = "unsigned",
                          power = 6,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3,
                          saveTOMs = FALSE)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)


# Plot the dendrogram
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

MEs <- bwnet$MEs
module_membership <- cor(bwnet$MEs, wgcna_data, use = "p")
module_membership_pvals <- corPvalueStudent(module_membership, nrow(wgcna_data))

hub_genes <- chooseTopHubInEachModule(wgcna_data, bwnet$colors)


module_gene_mapping <- as.data.frame(bwnet$colors)
module_gene_mapping <- rownames_to_column(module_gene_mapping)
colnames(module_gene_mapping) <- c("gene_id", "module")


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

module_gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = module_gene_mapping$gene_id,
                        mart = ensembl)

module_gene_mapping <- merge(module_gene_id, module_gene_mapping, by.x = "ensembl_gene_id", by.y = "gene_id", all = T)
module_gene_mapping[module_gene_mapping == ""] <- NA


modules <- unique(module_gene_mapping$module)
module_clusters <- list()
for (i in seq_along(modules)) {
  module <- modules[i]
  x <- module_gene_mapping[module_gene_mapping$module == module, ]
  module_clusters[[module]] <- x
}


save(bwnet, sft, wgcna_data, module_clusters, file = "BRCA/RData/lumA_WGCNA.RData")


## module enrichment 
library(clusterProfiler)
library(org.Hs.eg.db)

enrichModule <- function(module_genes) {
  enrichResult <- enrichGO(gene = module_genes, 
                           OrgDb = org.Hs.eg.db, 
                           keyType = "ENSEMBL", 
                           ont = "ALL", 
                           pAdjustMethod = "BH")
  result <- enrichResult@result
  
  return(result)
}

module_enrichment <- list()
for (i in seq_along(module_clusters)) {
  gene_set <- module_clusters[[i]]
  colour <- names(module_clusters)[i]
  print(paste0("Enriching module ", colour, ": ", i, " of ", length(names(module_clusters))))
  
  enrichedModule <- enrichModule(gene_set$ensembl_gene_id)
  module_enrichment[[colour]] <- enrichedModule
  
  rm(gene_set, colour, enrichedModule)
}


save(module_enrichment, file = "BRCA/RData/lumA_module_enrichment.RData")



## cross reference modules with OpenTargets
load("BRCA/RData/OpenTargets_NCT_filtered.RData")
OpenTargets_NCT_unique <- OpenTargets_NCT_filtered[!duplicated(OpenTargets_NCT_filtered$Target.Approved.Symbol), ]

OpenTargets_data <- subset(OpenTargets_NCT_unique, select = c("Target.Approved.Symbol", "Drug.Name"))

module_OpenTargets <- list()
drug_count <- data.frame(module = character(), drug_count = integer(), module_size = integer(), stringsAsFactors = FALSE)
for (i in seq_along(module_clusters)) {
  module <- module_clusters[[i]]
  colour <- names(module_clusters)[i]
  OpenTargets_module <- merge(module, OpenTargets_data, by.x = "external_gene_name", by.y = "Target.Approved.Symbol", all.x = T)
  
  module_OpenTargets[[colour]] <- OpenTargets_module
  
  # count number of drugs in each module
  num_drugs <- nrow(unique(OpenTargets_module[!is.na(OpenTargets_module$Drug.Name), "Drug.Name", drop = FALSE]))
  drug_count <- rbind(drug_count, data.frame(module = colour, drug_count = num_drugs, module_size = nrow(module), stringsAsFactors = FALSE))

  # tidy up
  rm(module, OpenTargets_module, num_drugs)
}

save(module_OpenTargets, drug_count, file = "BRCA/RData/LumA/OpenTargets_modules.RData")



