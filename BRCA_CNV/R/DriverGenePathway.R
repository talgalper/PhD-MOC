library(devtools)
library(DriverGenePathway)
library(data.table)

load("RData/BRCA_SNV_data.RData")

# if these files not availiable run the DriverGene function with muation data only to install them.
C <- as.data.frame(data.table::fread(file = "exome_full192.coverage.txt"))
dict <- as.data.frame(data.table::fread(file = "mutation_type_dictionary_file.txt"))
V <- as.data.frame(data.table::fread(file = "gene.covariates.txt"))

driver_genes <- DriverGene(Mutation = SNV_data, 
                           Coverage = C, 
                           Covariate = V,
                           MutationDict = dict,
                           chr_files_directory = "chr_files_hg19/")



