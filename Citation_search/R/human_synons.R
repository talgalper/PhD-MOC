library(geneSynonym)
library(progress)
library(reshape2)

syno9606 <- as.data.frame(syno9606)


pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", 
                       total = nrow(syno9606))

syno_list <- list()
# Loop through each row of the data frame
for (i in 1:nrow(syno9606)) {
  # Split the string at the '|' character
  terms <- unlist(strsplit(syno9606$syno9606[i], split = "\\|"))
  # The first term is the reference
  reference <- terms[1]
  # The rest are synonyms
  synonyms <- terms[-1]
  # Add to the list
  if (length(synonyms) == 0) {
    syno_list[[reference]] <- reference
  }
  else {
    syno_list[[reference]] <- synonyms
  }
  
  pb$tick()
  rm(i, reference, synonyms, terms)
}
rm(pb)

syno_unlist <- melt(syno_list)
colnames(syno_unlist) <- c("synonym", "ref_term")

write.csv(syno_unlist, "data/human_syno.csv", row.names = F)




