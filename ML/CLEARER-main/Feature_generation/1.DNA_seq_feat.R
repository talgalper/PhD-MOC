# load required packages
library(rDNAse)
library(seqinr)

# load DNA sequences
fasta <- read.fasta("/path/to/Sc_DNA.fa",seqtype = "DNA", forceDNAtolower = F,as.string = T)
names <- names(fasta)
head(fasta)
seqs <- unlist(fasta)
names(seqs) <- names
seqs <- seqs[(sapply(seqs, dnacheck))]

# generate features
DAC <- t(sapply(seqs, function(x) extrDAC(x,allprop = T)))
DACC <- t(sapply(seqs, function(x) extrDACC(x,allprop = T)))
DCC <- t(sapply(seqs, function(x) extrDCC(x,allprop = T)))
PseDNC <- t(sapply(seqs, function(x) extrPseDNC(x)))
extrPseKNC <- t(sapply(seqs, function(x) extrPseKNC(x)))
TAC <- t(sapply(seqs, function(x) extrTAC(x,allprop = T)))
TACC <- t(sapply(seqs, function(x) extrTACC(x,allprop = T)))
TCC <- t(sapply(seqs, function(x) extrTCC(x)))
kmer2 <- t(sapply(seqs, function(x) kmer(x, k = 2, upto = FALSE, normalize = FALSE, reverse = FALSE)))
kmer3 <- t(sapply(seqs, function(x) kmer(x, k = 3, upto = FALSE, normalize = FALSE, reverse = FALSE)))
kmer4 <- t(sapply(seqs, function(x) kmer(x, k = 4, upto = FALSE, normalize = FALSE, reverse = FALSE)))
kmer5 <- t(sapply(seqs, function(x) kmer(x, k = 5, upto = FALSE, normalize = FALSE, reverse = FALSE)))
kmer6 <- t(sapply(seqs, function(x) kmer(x, k = 6, upto = FALSE, normalize = FALSE, reverse = FALSE)))
kmer7 <- t(sapply(seqs, function(x) kmer(x, k = 7, upto = FALSE, normalize = FALSE, reverse = FALSE)))

# combine
results <- as.data.frame(cbind(
  DAC,
  DACC,
  DCC,
  PseDNC,
  extrPseKNC,
  TAC,
  TACC,
  TCC,
  kmer2,
  kmer3,
  kmer4,
  kmer5,
  kmer6,
  kmer7))

results$Gene <- rownames(results)

# save
write.table(results,"DNA_Seq_feat.csv", sep = "\t", row.names = F) 
