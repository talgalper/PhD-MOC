library(edgeR)
library(tidyverse)

tumour_ec <- read.table("RNAseqDB/data/brca-rsem-count-tcga-t.txt", sep = "\t", header = T)
tumour_ec <- column_to_rownames(tumour_ec, "Hugo_Symbol")
tumour_ec <- tumour_ec[,-1]
normal_ec <- read.table("RNAseqDB/data/brca-rsem-count-tcga.txt", sep = "\t", header = T)
normal_ec <- column_to_rownames(normal_ec, "Hugo_Symbol")
normal_ec <- normal_ec[,-1]
GTEx_ec <- read.table("RNAseqDB/data/breast-rsem-count-gtex.txt", sep = "\t", header = T)
GTEx_ec <- column_to_rownames(GTEx_ec, "Hugo_Symbol")
GTEx_ec <- GTEx_ec[,-1]

data <- cbind(GTEx_ec, tumour_ec)

n1 <-  length(names(GTEx_ec))
n2 <-  length(names(tumour_ec))

design <- cbind(Grp1=1, Grp2vs1=c(rep(1,n1), rep(2,n2)) )

v_max <-  max(data)
v_min <-  min(data)

if(v_min >= 0 & v_max >= 1000){  # If the expression data was not normalized before
  # filter out low-count genes
  isexpr <- rowSums(cpm(data) > 5) >= 10
  
  DGE <- DGEList(data[isexpr,])
  DGE <- calcNormFactors(DGE,method =c("TMM"))
  
  # perform voom normalization
  data_norm <- voom(DGE, design, plot=TRUE)
}





