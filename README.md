# PhD-MOC

This repository encompasses all coding work in Tal Galper's PhD project

## Thesis context
This repository aggregates the computational work supporting my PhD research into the molecular underpinnings and therapeutic vulnerabilities of mucinous ovarian carcinoma (MOC). The code base spans RNA-seq processing, network propagation, druggability profiling, and downstream integration analyses that collectively prioritise candidate targets for experimental follow-up. Note that some large files are not found here and are either on the university OneDrive and/or in my posession. Can be provided on request.

## Pipeline overview
The final pipeline combines multiple analytic layers:
- **Transcriptomic preprocessing (`MOC_pipe/`)** – generates differential expression contrasts (MOC vs. benign ovarian, GTEx, and TCGA cohorts) and quality-control visualisations used to seed downstream analyses.
- **Network diffusion (`Hierarchical_HotNet/`)** – applies Hierarchical HotNet to propagate differential signals through STRING-derived interaction networks and returns subnetworks enriched for MOC-specific perturbations.
- **Druggability scoring (`Druggability_analysis/`)** – integrates curated pharmacology datasets (DrugBank, DGIdb, OpenTargets) with structure-based predictors (Fpocket, PocketMiner) to annotate network hits with targetability metrics.
- **Machine learning prioritisation (`ML/`)** – trains bagged random forest ensembles on aggregated feature matrices to rank genes for small-molecule oncology target potential and evaluate model generalisation to clinical target sets.


## Reproducing MOC Results
1. **MOC expression processing:** launch the R project in [`MOC_pipe/`](MOC_pipe/) and run the differential expression scripts under `R/` to regenerate contrast tables in `results/`.
2. **Subnetwork Identification:** feed the ranked gene lists into [`Hierarchical_HotNet/`](Hierarchical_HotNet/) using `HHnet_run.sh` to derive significant subnet hierarchies saved under `MOC/results/`. Follow the script to enrich the 
likely small subnetwork with nearest neighbours for greater bioloigcal context.
3. **Integration:** combine network metrics from Hierarchical HotNet in [`MOC_pipe/R/MOC_RS.R`](MOC_pipe/R/MOC_RS.R) following the script workflow which pulls citaiton and druggability data from respective directories for post hoc analysis.  

### Optional (fix this section)
If you would like to re-run any of scripts to generate the mmachine learning or druggability data you can do so as follows:

- **Druggability annotation:** open [`Druggability_analysis/`](Druggability_analysis/) and run relevant structural database downlaod scripts if you wish to obtain the latest version of the datasets (optionally visit the database webpages). Run the drug binding pocket prediction tools [`Fpocket/`](Druggability_analysis/Fpocket/) and [`PocketMiner`](Druggability_analysis/PocketMiner/) on the updated structure set(s).

- **Machine learning ranking:** use [`ML/`](ML/) to train or refresh the bagged random forest models (`all_ML_func.R`, `all_ML_run.R`), writing ranked predictions into `results/`.


## Repository layout
- [`MOC_pipe/`](MOC_pipe/) – RNA-seq ingestion, differential expression, PCA visualisations, and survival-linked metadata underpinning the MOC-specific signal set.
- [`BRCA_pipe/`](BRCA_pipe/) – early breast cancer expression and drug repurposing workflow used to benchmark the  pipeline against TCGA BRCA cohorts using known breast cancer targets as validation. Retain for reproducing comparative figures and prototype code used in validation thesis chapters.
- [`Hierarchical_HotNet/`](Hierarchical_HotNet/) – Python/Fortran implementation of the Hierarchical HotNet algorithm together with STRING networks and shell wrappers for batch processing.
- [`Druggability_analysis/`](Druggability_analysis/) – data folders for DGIdb, DrugBank, OpenTargets, structural prediction outputs, and R scripts that synthesise targetability evidence.
- [`ML/`](ML/) – machine learning experiments covering feature matrix construction, target/negative set curation, and ensemble evaluation utilities, with intermediate artefacts in `RData/` and `results/`.
- [`Citation_search/`](Citation_search/) – scripts for PubTator-powered literature mining and synonym harmonisation used to weight targets by publication evidence.
- [`WGCNA/`](WGCNA/) – weighted gene co-expression network analyses (WGCNA), and integration notebooks tying network communities to druggability and literature signals. Details all experimental code used to benchmark and test this tool for thesis purposes.
- `Citation_search/results/`, `WGCNA/`, and `MOC_pipe/results/` – downstream outputs referenced by integration notebooks and figures.

## Archived and legacy workflows
The following directories preserve historical pipelines that informed validation experiments and comparative analyses:
- [`BRCA_CNV/`](BRCA_CNV/) – exploratory copy-number driver modelling that was leveraged to validate CNV calls and prioritisation strategies mirrored in the MOC cohort.
- [`RNA-Seq_pipeline/`](RNA-Seq_pipeline/) – initial generic RNA-seq processing scripts (now superseded by `MOC_pipe/`) kept for provenance and to reproduce baseline processing used in validation appendices.
- [`multiWGCNA-BRCA/`](multiWGCNA-BRCA/) – prototype multi-dataset WGCNA analyses on BRCA cohorts that established module preservation benchmarks for the current MOC network integration results.

These legacy paths should remain untouched unless you need to recreate validation experiments or trace historical decisions referenced in the thesis narrative.
