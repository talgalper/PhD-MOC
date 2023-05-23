library(tidyverse)
library(edgeR)
library(AnnotationDbi)


stage_II <- read.table("attemp_3_data/stage_II_master_df.tsv", header = T, sep = "\t", row.names = "gene_id")
colnames(stage_II) <- c(paste("stage_II_fpkm_unstranded", 1:23))

stage_III <- read.table("attemp_3_data/stage_III_master_df.tsv", header = T, sep = "\t", row.names = "gene_id")
colnames(stage_III) <- c(paste("stage_III_fpkm_unstranded", 1:23))

tcga_data <- cbind(stage_II, stage_III)


group <- factor(c(rep(1, 23), rep(2, 23)))

data <- DGEList(counts = tcga_data, group = group)

design <- model.matrix(~group)

common <- estimateGLMCommonDisp(data, design, verbose = T)

trend <- estimateGLMTrendedDisp(common, design)

tagwise <- estimateGLMTagwiseDisp(trend, design)

plotBCV(tagwise)

fit <- glmFit(tagwise, design)

lrt <- glmLRT(fit, coef = 2)

toptags <- topTags(lrt, n = Inf)

summary(dif_exp <- decideTestsDGE(lrt, p = 0.1, adjust = "fdr", lfc = 1))

dif_exp_genes <- rownames(tagwise)[as.logical(dif_exp)]

plotSmear(lrt, de.tags = dif_exp_genes)

hits <- toptags$table[toptags$table$FDR < 0.1, ]
colnames <- colnames(hits)
hits$gene_id <- rownames(hits)
hits <- hits[,c("gene_id", colnames)]


write_tsv(hits, "attemp_3_data/tcga_stage_II_vs_III_edgeR_dif_exp(p=0.1).tsv")


# volcano plot
data <- read.table("attemp_3_data/tcga_stage_II_vs_III_edgeR_dif_exp(p=0.1).tsv", header = T, sep = "\t")
ggplot(data, aes(x=logFC, y=PValue)) + geom_point()






