library(tidyverse)

benign <- read.table("attempt_2_data/benign_master_df.tsv", header = T, sep = "\t", row.names = "gene_id")
colnames(benign) <- c(paste("benign_fpkm_unstranded", 1:48))

disease <- read.table("attempt_2_data/disease_master_df.tsv", header = T, sep = "\t", row.names = "gene_id")
colnames(disease) <- c(paste("disease_fpkm_unstranded", 1:254))


# read in list of file names
ben_files <- list.files("attempt_2_data/tcga_files/benign/.","tsv")
dis_files <- list.files("attempt_2_data/tcga_files/disease/.","tsv")

# make same no. of rows. empty = NA
length(ben_files) <- length(dis_files)

# match up rows
compare <- data.frame(dis_files, ben_files)
compare$matched_ben <- compare$ben_files[match(compare$dis_files, compare$ben_files)]
compare <- subset(compare, select = c("dis_files", "matched_ben"))

# delete unmatched rows
compare <- na.omit(compare)

# create list of matched files and copy to new directory
files <- compare$matched_ben

for (file in files) {
  file.copy(from = file.path("attempt_2_data/tcga_files/disease/", file), 
            to = "attempt_2_data/tcga_files/matched/disease/")
}

for (file in files) {
  file.copy(from = file.path("attempt_2_data/tcga_files/benign/", file), 
            to = "attempt_2_data/tcga_files/matched/benign/")
}
