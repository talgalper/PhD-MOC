### R code from vignette source 'Tutorial.Rnw'

###################################################
### code chunk number 1: Tutorial.Rnw:55-57
###################################################
options(width=60)
ps.options(family="sans")


###################################################
### code chunk number 2: Tutorial.Rnw:60-64
###################################################
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)


###################################################
### code chunk number 3: Tutorial.Rnw:68-71
###################################################
pvals <- cbind(t=dataLym$t.pval, s=dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order=2, plot=FALSE)


###################################################
### code chunk number 4: Tutorial.Rnw:75-78
###################################################
subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet


###################################################
### code chunk number 5: Tutorial.Rnw:84-86
###################################################
fb <- fitBumModel(pval, plot=FALSE)
scores <- scoreNodes(subnet, fb, fdr=0.001)


###################################################
### code chunk number 6: Tutorial.Rnw:93-96
###################################################
module <- runFastHeinz(subnet, scores)
logFC <- dataLym$diff
names(logFC) <- dataLym$label


###################################################
### code chunk number 7: Tutorial.Rnw:103-104
###################################################
plotModule(module, scores=scores, diff.expr=logFC)


###################################################
### code chunk number 8: Tutorial.Rnw:159-163
###################################################
library(BioNet)
library(DLBCL)
data(exprLym)
data(interactome)


###################################################
### code chunk number 9: Tutorial.Rnw:171-172
###################################################
exprLym


###################################################
### code chunk number 10: Tutorial.Rnw:178-179
###################################################
interactome


###################################################
### code chunk number 11: Tutorial.Rnw:186-188
###################################################
network <- subNetwork(featureNames(exprLym), interactome)
network


###################################################
### code chunk number 12: Tutorial.Rnw:193-195
###################################################
network <- largestComp(network)
network


###################################################
### code chunk number 13: Tutorial.Rnw:205-209
###################################################
library(genefilter)
library(impute)
expressions <- impute.knn(exprs(exprLym))$data
t.test <- rowttests(expressions, fac=exprLym$Subgroup)


###################################################
### code chunk number 14: Tutorial.Rnw:212-213
###################################################
t.test[1:10, ]


###################################################
### code chunk number 15: Tutorial.Rnw:220-223
###################################################
library(xtable)
top.table <- xtable(t.test[1:10,], display=c("s", "f", "f", "f"))
print(top.table, floating=FALSE)                                    


###################################################
### code chunk number 16: Tutorial.Rnw:233-238
###################################################
data(dataLym)
ttest.pval <- t.test[, "p.value"]
surv.pval <- dataLym$s.pval
names(surv.pval) <- dataLym$label
pvals <- cbind(ttest.pval, surv.pval)


###################################################
### code chunk number 17: Tutorial.Rnw:247-248
###################################################
pval <- aggrPvals(pvals, order=2, plot=FALSE)


###################################################
### code chunk number 18: Tutorial.Rnw:255-257
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb   


###################################################
### code chunk number 19: Tutorial.Rnw:260-265
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)
dev.off()


###################################################
### code chunk number 20: Tutorial.Rnw:280-281
###################################################
plotLLSurface(pval, fb)


###################################################
### code chunk number 21: Tutorial.Rnw:290-291
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=0.001)


###################################################
### code chunk number 22: Tutorial.Rnw:298-301
###################################################
network <- rmSelfLoops(network)                                                                                                     
writeHeinzEdges(network=network, file="lymphoma_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="lymphoma_nodes_001", node.scores = scores)


###################################################
### code chunk number 23: Tutorial.Rnw:318-320
###################################################
datadir <- file.path(path.package("BioNet"), "extdata") 
dir(datadir)


###################################################
### code chunk number 24: Tutorial.Rnw:326-329
###################################################
module <- readHeinzGraph(node.file=file.path(datadir, "lymphoma_nodes_001.txt.0.hnz"), network=network)
diff <- t.test[, "dm"]
names(diff) <- rownames(t.test) 


###################################################
### code chunk number 25: Tutorial.Rnw:332-333
###################################################
plotModule(module, diff.expr=diff, scores=scores)


###################################################
### code chunk number 26: Tutorial.Rnw:355-358
###################################################
sum(scores[nodes(module)])
sum(scores[nodes(module)]>0)
sum(scores[nodes(module)]<0)


###################################################
### code chunk number 27: Tutorial.Rnw:380-385
###################################################
library(BioNet)
library(DLBCL)
library(ALL)
data(ALL)
data(interactome)


###################################################
### code chunk number 28: Tutorial.Rnw:392-393
###################################################
ALL


###################################################
### code chunk number 29: Tutorial.Rnw:401-402
###################################################
interactome


###################################################
### code chunk number 30: Tutorial.Rnw:412-414
###################################################
mapped.eset <- mapByVar(ALL, network=interactome, attr="geneID")
mapped.eset[1:5,1:5]


###################################################
### code chunk number 31: Tutorial.Rnw:419-420
###################################################
length(intersect(rownames(mapped.eset), nodes(interactome)))


###################################################
### code chunk number 32: Tutorial.Rnw:429-434
###################################################
network <- subNetwork(rownames(mapped.eset), interactome)
network
network <- largestComp(network)
network <- rmSelfLoops(network)
network


###################################################
### code chunk number 33: Tutorial.Rnw:446-454
###################################################
library(limma)
design <- model.matrix(~ -1+ factor(c(substr(unlist(ALL$BT), 0, 1))))
colnames(design)<- c("B", "T")
contrast.matrix <- makeContrasts(B-T, levels=design)
contrast.matrix 
fit <- lmFit(mapped.eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


###################################################
### code chunk number 34: Tutorial.Rnw:459-460
###################################################
pval <- fit2$p.value[,1]                                


###################################################
### code chunk number 35: Tutorial.Rnw:471-473
###################################################
fb <- fitBumModel(pval, plot=FALSE)
fb


###################################################
### code chunk number 36: Tutorial.Rnw:476-480
###################################################
dev.new(width=13, height=7)
par(mfrow=c(1,2))
hist(fb)
plot(fb)


###################################################
### code chunk number 37: Tutorial.Rnw:492-493
###################################################
scores <- scoreNodes(network=network, fb=fb, fdr=1e-14)


###################################################
### code chunk number 38: Tutorial.Rnw:499-501
###################################################
writeHeinzEdges(network=network, file="ALL_edges_001", use.score=FALSE)
writeHeinzNodes(network=network, file="ALL_nodes_001", node.scores = scores)


###################################################
### code chunk number 39: Tutorial.Rnw:522-524
###################################################
datadir <- file.path(path.package("BioNet"), "extdata") 
module <- readHeinzGraph(node.file=file.path(datadir, "ALL_nodes_001.txt.0.hnz"), network=network)


###################################################
### code chunk number 40: Tutorial.Rnw:529-534
###################################################
nodeDataDefaults(module, attr="diff") <- ""
nodeData(module, n=nodes(module), attr="diff") <- fit2$coefficients[nodes(module),1]
nodeDataDefaults(module, attr="score") <- ""
nodeData(module, n=nodes(module), attr="score") <- scores[nodes(module)]
nodeData(module)[1]


###################################################
### code chunk number 41: Tutorial.Rnw:540-541
###################################################
saveNetwork(module, file="ALL_module", type="XGMML")


###################################################
### code chunk number 42: Tutorial.Rnw:573-581 (eval = FALSE)
###################################################
## j.repl <- 100
## resampling.pvals <- list()
## for(i in 1:j.repl)
## {
##   resampling.result <- resamplingPvalues(exprMat=mapped.eset, groups=factor(c(substr(unlist(ALL$BT), 0, 1))), resampleMat=FALSE, alternative="two.sided")
##   resampling.pvals[[i]] <- resampling.result$p.values
##   print(i)
## }


###################################################
### code chunk number 43: Tutorial.Rnw:590-596 (eval = FALSE)
###################################################
## fb <- lapply(resampling.pvals, fitBumModel, plot=FALSE, starts=1)
## resampling.scores <- c()
## for(i in 1:j.repl)
## {
##   resampling.scores[[i]] <- scoreNodes(network=network, fb=fb[[i]], fdr=1e-14)
## }


###################################################
### code chunk number 44: Tutorial.Rnw:601-603 (eval = FALSE)
###################################################
## score.mat <- as.data.frame(resampling.scores)
## colnames(score.mat) <- paste("resample", (1:j.repl), sep="")


###################################################
### code chunk number 45: Tutorial.Rnw:608-610 (eval = FALSE)
###################################################
## writeHeinzEdges(network=network, file="ALL_e_resample", use.score=FALSE)
## writeHeinzNodes(network=network, file="ALL_n_resample", node.scores = score.mat)


###################################################
### code chunk number 46: Tutorial.Rnw:623-625
###################################################
datadir <- file.path(path.package("BioNet"), "extdata") 
modules <- readHeinzGraph(node.file=file.path(datadir, "ALL_n_resample.txt.0.hnz"), network=network)


###################################################
### code chunk number 47: Tutorial.Rnw:634-636
###################################################
cons.scores <- consensusScores(modules, network)
writeHeinz(network=network, file="ALL_cons", node.scores=cons.scores$N.scores, edge.scores=cons.scores$E.scores)


###################################################
### code chunk number 48: Tutorial.Rnw:643-648
###################################################
datadir <- file.path(path.package("BioNet"), "extdata") 
cons.module <- readHeinzGraph(node.file=file.path(datadir, "ALL_cons_n.txt.0.hnz"), network=network)
cons.edges <- sortedEdgeList(cons.module)
E.width <- 1+cons.scores$E.freq[cons.edges]*10
N.size <- 1+cons.scores$N.freq[nodes(cons.module)]*10


###################################################
### code chunk number 49: Tutorial.Rnw:652-653
###################################################
plotModule(cons.module, edge.width=E.width, vertex.size=N.size, edge.label=cons.scores$E.freq[cons.edges]*100, edge.label.cex=0.6)

