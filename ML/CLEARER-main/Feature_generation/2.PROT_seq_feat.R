# attach required packages
library(seqinr)
library(stringr)
library(protr)
# read in protein seqeunces
Seq <- read.fasta("/path/to/Sc_prot.fa",forceDNAtolower = F,seqtype = "AA",as.string = T)
names <- names(Seq)
Seq <- getSequence(Seq,as.string = T)
Seq <- as.character(unlist(Seq))
names(Seq) <- names

# remove sequences < 60 amino acids if there are any
#len <- c()
#for (i in 1:length(Seq)){
#  len[i] <- length(s2c(as.character(Seq[i])))
#}
#Seq <- Seq[-c(which(len < 60))]
#names <- names[-c(which(len < 60))]

# seqinr features
names <- c()
Tiny <- c()
Small <- c()
Aliphatic <- c()
Aromatic <- c()
Non.polar <- c()
Polar <- c()
Charged <- c()
Basic <- c()
Acidic <- c()
Pi <- c()

for (i in 1:length(Seq)) {
  names[i] <- as.character(names[i])
  Tiny[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Tiny
  Small[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Small
  Aliphatic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Aliphatic
  Aromatic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Aromatic
  Non.polar[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Non.polar
  Polar[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Polar
  Charged[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Charged
  Basic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Basic
  Acidic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Acidic
  Pi[i] <- AAstat(s2c(as.character(Seq[1])),plot = F)$Pi
}

datalist = list()
for (i in 1:length(Seq)) {
  dat <- AAstat(s2c(as.character(Seq[i])),plot =F)$Compo
  dat$i <- i  
  datalist[[i]] <- dat 
}
big_data = as.table(do.call(rbind, datalist))
dim(big_data)
write.csv(big_data,"AA_numbers_of_genes.csv",row.names = F)
big_data <- read.csv("AA_numbers_of_genes.csv")
sapply(big_data, class)

comb <- as.data.frame(cbind(names(Seq),big_data,Tiny,
                            Small,Aliphatic,Non.polar,Polar,Charged,Basic,Acidic,Pi ))
head(comb)

# protr features
seqs <- Seq
head(seqs)
## remove * character at the end from each sequence ##
seqs <- str_sub(seqs, 1, str_length(seqs)-1)
names(seqs) <- names
seqs <- seqs[(sapply(seqs, protcheck))]

AAC <- t(sapply(seqs, extractAAC))
head(AAC)
dim(AAC)
DC <- t(sapply(seqs, extractDC))
dim(DC)
TC <- t(sapply(seqs, extractTC))
dim(TC)
MoreauBroto <- t(sapply(seqs, extractMoreauBroto))
dim(MoreauBroto)
Moran <- t(sapply(seqs, extractMoran))
dim(Moran)
Geary <- t(sapply(seqs, extractGeary))
dim(Geary)
CTDC <- t(sapply(seqs, extractCTDC))
dim(CTDC)
CTDT <- t(sapply(seqs, extractCTDT))
dim(CTDT)
CTDD <- t(sapply(seqs, extractCTDD))
dim(CTDD)
CTriad <- t(sapply(seqs, extractCTriad))
dim(CTriad)
SOCN <- t(sapply(seqs, extractSOCN))
dim(SOCN)
QSO <- t(sapply(seqs, extractQSO))
dim(QSO)
PAAC <- t(sapply(seqs, extractPAAC))
dim(PAAC)
APAAC <- t(sapply(seqs, extractAPAAC))
dim(APAAC)

# combine and save
colnames(comb)[1] <- "Gene"
comb2 <- as.data.frame(cbind(AAC,DC,TC,MoreauBroto,Moran,Geary,CTDC,CTDT,CTDD,CTriad,SOCN,QSO,PAAC,APAAC))
comb2$Gene <- comb$Gene

results <- merge.data.frame(comb,comb2,"Gene")
dim(results)

write.table(results,"protein_physiochemical_features.csv", row.names = F, sep = "\t")
