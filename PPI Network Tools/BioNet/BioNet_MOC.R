library(BioNet)
library(DLBCL)
library(PCSF)

df <- read.table("data/wgcna-10th-power-MOC-edgelist-cost-filterDot01.txt", header = T)
interactome <- construct_interactome(df)
interactomne <- as_graphnel(interactome)





pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order = 2, plot = FALSE)

subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet

fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.001)

module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label

plotModule(module, scores = scores, diff.expr = logFC)
