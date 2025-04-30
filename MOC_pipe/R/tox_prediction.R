### toxicity prediction ###

library(data.table)
# GenomeAd data - healthy cell toxicity 
GenomAD <- fread("~/Downloads/gnomad.v4.1.constraint_metrics.tsv", sep = "\t") #

GenomAD_filt <- GenomAD[mane_select == TRUE | canonical == TRUE, ]
GenomAD_filt <- GenomAD_filt[, c("gene", "gene_id", "lof_hc_lc.pLI", "lof.pLI", "lof.oe_ci.upper", 
                                 "lof.oe", "lof.oe_ci.upper_rank", "lof.oe_ci.upper_bin_decile", 
                                 "mis.z_score", "mis.oe", "constraint_flags")]
GenomAD_filt <- GenomAD_filt[order(gene)]
GenomAD_filt <- GenomAD_filt[!duplicated(gene, fromLast = TRUE), ]

# lof.oe_ci.upper: LOEUF: upper bound of 90% confidence interval for o/e ratio for high confidence pLoF variants 
#                  (lower values indicate more constrained)       
# lof.pLI: Probability of loss-of-function intolerance; probability that transcript falls into distribution of 
#          haploinsufficient genes (~21% o/e pLoF ratio;  computed from high confidence pLoF gnomAD data)
GenomAD_filt$essential <- ifelse(GenomAD_filt$lof.pLI >= 0.9 & GenomAD_filt$lof.oe_ci.upper <= 0.35, TRUE, FALSE)


# depmap data - cancer selectivity
effects <- fread("~/Downloads/CRISPRGeneEffect.csv")
models <- fread("~/Downloads/Model.csv")

# subset only cancer mnodels
cancer_ids <- models[OncotreePrimaryDisease != "Non-Cancerous", ]
cancer_ids <- cancer_ids[SampleCollectionSite == "ovary", ]
cancer_ids <- cancer_ids[OncotreeSubtype == "Mucinous Ovarian Cancer", ]
cancer_ids <- cancer_ids$ModelID

eff_can <- as.data.frame(effects[V1 %in% cancer_ids]) # convert to data.frame, data.table being weird
gene_mat <- eff_can[, c(2:ncol(eff_can))]

long_df <- data.frame(
  ModelID = rep(eff_can$V1, each = ncol(gene_mat)),
  gene = rep(colnames(gene_mat), times = nrow(gene_mat)),
  chronos = as.vector(as.matrix(gene_mat)),
  stringsAsFactors = FALSE
)

med_effect <- tapply(long_df$chronos,
                     long_df$gene,
                     median,
                     na.rm=TRUE)

mean_effect <- tapply(long_df$chronos,
                     long_df$gene,
                     mean,
                     na.rm=TRUE)

gene_summary <- data.frame(
  gene = names(med_effect),
  median_effect = as.numeric(med_effect),
  mean_effect = as.numeric(mean_effect),
  row.names = NULL,
  stringsAsFactors = FALSE
)


essentials <- read.csv("~/Downloads/CRISPRInferredCommonEssentials.csv")
gene_summary$comEss <- ifelse(gene_summary$gene %in% essentials$Essentials, TRUE, FALSE)



# merge data together
gene_summary$gene <- sub("\\s*\\(.*", "", gene_summary$gene)

tox_summary <- merge(GenomAD_filt[,c(1,4,5)], gene_summary, by = "gene", all = TRUE)
# remove rows with all NA values
tox_summary <- tox_summary[!which(apply(tox_summary, 1, function(r) all(is.na(r)))), ]

tox_summary$category <- with(tox_summary,
                    ifelse(
                      # Cancer-selective hits
                      mean_effect <= -0.5 &
                        (lof.pLI <  0.9 & lof.oe_ci.upper > 0.35),
                      "cancer-selective",
                      
                      ifelse(
                        # Pan-essential
                        mean_effect <= -0.5 &
                          (lof.pLI >= 0.9 | lof.oe_ci.upper <= 0.35),
                        "pan-essential",
                        
                        ifelse(
                          # Non-essential
                          mean_effect >  -0.2 &
                            (lof.pLI <  0.9 & lof.oe_ci.upper > 0.35),
                          "non-essential",
                          
                          # Everything else
                          "uncertain"
                        ))))


# CHRONOS = 0: no dependency in cancer cell line
# CHRONOS = -1: essential in cancer cell line

# lof.oe = 0: intolerant
# lof.oe = 1: non-essential

# lof.pLI = 0: non-essential
# lof.pLI = 1 intolerance

# cancer selective: mostly harms cancer cells
# pan-essential: harms both cancer and normal cells
# non-essential: no effect on cancer or normal cells
# uncertain: doesnt meet any of the above criteria


# read in results from rank sensitivity analysis
MOCvsBEN_RS <- read.csv("results/HHnet_RS_ML_overlap.csv")
MOCvsGTEx_RS <- read.csv("results/MOC_vs_GTEx/HHnet_RS_ML_overlap.csv")

# combine tox with RS 
MOCvsBEN_RS$relative_rank <- rownames(MOCvsBEN_RS)
MOCvsBEN_RS <- MOCvsBEN_RS[, c(1:4,10,5:9)]
RS_tox <- merge(MOCvsBEN_RS, tox_summary, by.x = "external_gene_name", by.y = "gene", all.x = TRUE)
RS_tox <- RS_tox[order(RS_tox$avg_rank), ]
rownames(RS_tox) <- NULL

# annotate known drugs
load("~/Desktop/TTD_master.RData")
TTD_master <- TTD_master[!is.na(TTD_master$TTDDrugID), ]
pattern <- paste(c("cancer", "carcinoma", "leukemia", "leukaemia", "neoplasm", "metastases", "tumour", "adenoma", "sarcoma"), 
                 collapse = "|")
TTD_master <- TTD_master[grepl(pattern, TTD_master$INDICATION, ignore.case = TRUE), ]


RS_tox_drug <- merge(RS_tox, TTD_master[, c(1,8,10)], by.x = "external_gene_name", by.y = "all_target_genes", all.x = T)
RS_tox_drug <- RS_tox_drug[order(RS_tox_drug$avg_rank), ]
RS_tox_drug <- RS_tox_drug[!is.na(RS_tox_drug$INDICATION) | !duplicated(RS_tox_drug$external_gene_name), ]
rownames(RS_tox_drug) <- NULL

# simplified version
RS_tox$TTD_onco_drug <- ifelse(RS_tox$external_gene_name %in% TTD_master$all_target_genes, TRUE, FALSE)

write.csv(RS_tox, "~/Desktop/temp.csv", row.names = FALSE)





