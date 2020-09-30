# HMO_GeneReaction_pred
A flux prediction and network-based approach to resolving gene-reaction relations through transcriptomic and metabolomic data integration

### Corresponding Publication
Elucidating Human Milk Oligosaccharide biosynthetic genes through network-based multiomics integration
Benjamin P. Kellman, Anne Richelle, Jeong-Yeh Yang, Digantkumar Chapla, Austin W. T. Chiang, Julia Najera, Bokan Bao, Natalia Koga, Mahmoud A. Mohammad, Anders Bech Bruntse, Morey W. Haymond, Kelley W. Moremen, Lars Bode, Nathan E. Lewis
bioRxiv 2020.09.02.278663; doi: https://doi.org/10.1101/2020.09.02.278663

## Dependencies
- Hardware
  - Flux Analysis: macbook pro 16G RAM, 8 cores
  - Flux-Expression Analysis: HP Z620, 96G RAM, 32 cores, Intel Xeon CPU E5-2670 @ 2.6GHz
- Software
  - Flux Analysis: Matlab 2016b, Cobratoolbox 3.0
  - Flux-Expression Analysis: see code/R.env 

## Data
- Source data (data/data_raw)
  - Formatted data (data/data_raw/HMO_expss.xlsx): Expression data: Cohort 1 (GSE36936), Cohort 2 (GSE12669), HMO Concentration data: Cohort 1 (DOI:10.1101/693507), Cohort 2 (DOI:10.1101/2020.09.02.278663)
  - MIRAGE Annotation (data/data_raw/MIRAGE-1601485787699.xlsx): MIRAGE complient description of HPLC methods
  - HPLC measured HMO abundance has also been shared on Glycopost (**TBD**)
  - The exact raw data used is shared through zenodo ( doi: **TBD**)
- Exact raw/intermediate data files used ( zendodo doi: **TBD**)
  - Raw/Input data (raw.zip); data_HMO.xlsx is the central dataset including HMO concentrations and gene expression data. 
Other files are .csv exports from data_HMO.xlsx to circumvent a temporary bug in the xlsx reading package. They should not differ from the original file and are provided to facilitate code readability.
  - (models_split.zip) HMO biosynthesis models generated through network generation ( <github>/code/1.flux/A) and enumeration (<github>/code/1.flux/B)
  - (scores_spllit.zip) Correlation (R) and correlation p-values (P) between estimated flux and observed gene expression for corresponding genes (<github>/code/1.flux/C)
  - (big_mat.<bin/desc>.zip) Intermediate matrixes containing gene-linkage scores, models score, model classifications and other information generated/used by the Flux-Expression Comparision code (<github>/code/3.flux_expression)
- Output data (data/data_out)
  - Summary statistics including preference, and enrichemnt for genes across fundamental linkage reactions and aggregate model scores. These data are parsed and presented in table S2 of the corresponding publication (https://doi.org/10.1101/2020.09.02.278663)

## Code
The codebase is split into three sections: flux analysis, expression analysis and Flux-Expression Comparison. The code was written and run over several years and uses multiple packages. We have prioritized completeness in sharing this codebase to ensure our methods are as fully specified as possible. However, due to the completeness of the codebase and the length of runtimes, we have not been able to dockerize and reproduction will be non-trivial. To minimize the challenge, we have provided all intermediate calculations through zenodo. We are happy to provide support if necessary.

To run code, please first update paths. To make the codebase easier to understand, we have moved several files from their original position in the file structure.

### 1) Flux analysis (runtime: at least 1 month)
- code/1.flux/A_\*/MAINCODE_\* builds the complete HMO biosynthesis network using the fundamental reactions then trims the complete network using FVA to remove irrelevant reactions producing the reduced network (code/1.flux/B_*/Network_Red_ext.mat)
- code/1.flux/B_\*/parallel_enumerationSubNetworks.m runs a Mixed Integer Linear Program to enumerate every unique & minimal network within the reduced HMO biosynthetic network. FBA was run on each unique network to estimate flux.
- code/1.flux/C_\*/parallel_correlations.m computes the correlation between FBA estimted flux for each model (normalized to the flux of immidiate upstream reactions) and the gene expression of corresponding candidate genes. <R/P>_data*_C.mat are example files from a subset of models of the output correlation and p-values provided through zenodo for all models*

### 2) Expression Analysis and Gene Filtering (runtime: <1 min after loading GTEx )
- code/2.expression/GeneInclusion.ipynb code demonstrating the criteria and analysis used to filter gene candidates prior to the flux-expression comparison. (GTEx expression data should be downloaded to run this code, see code for the specific version)

### 3) Flux-Expression Comparison: Model and Gene Selection (runtime: 7-30 days)
- code/3.flux_expression/A_permuted_background/calc_null_distr.r calculates the random background used to estimate the signifiance of the support score
- code/3.flux_expression/B_largeModel_code_par.r uses the prior calculations to determine the consistency between estimate flux and candidate gene expression.
  - run_GLS: calculates the Gene-Linkage Score (GLS) for each reaction across all models
  - run_GLSbin: The gene linkage score is used to calculate the best candidate gene for each linkage-reaction in each model. "Best" is indicated by 1 or 0
  - run_LS: GLS is reduced to a linkage score (LS), the GLS for the best candidate gene for each linkage
  - run_MS: LS is reduced to model score (MS), the average linkage score within a model.
  - run_groups: Identifies the mean and std. dev of model scores for each cohort then designates the top 5% (relative to a normal distribution) as "high performing" models.
  - run_inclusion: (provided through zenodo, do not run) generates a binary matrix indicating the reactions included (1) or excluded (0) from the reducted network in each model
  - run_SC: calculates the proportion of ambiguous structures in all models and high-performing models. Enrichment is performed for over-representation of structures in the high-performing set relative to all models
  - run_Struc_CBX: this is not a relevant analysis but needs to be if the zenodo data is not used to initialize and intermediate data structure
  - run_GC: calculates the proportion of candidate genes selected for each model and high-performing models. Enrichment is performed for over-representation of selected (highest GLS) genes in the high-performing set relative to all models
  - (irrelevant) run_check, vis, cooc, colink, anyl_CBX, cluster, cluster_vis, cluster_vis_movie, CBX_extend, new_sels, gene_v_str, RAND, runPCA, inclusion_rank, quantile_norm: NOT used in this analysis, do not run
- code/3.flux_expression/C_GCagg.ipynb Aggregates Proportion, GLS and Model Contribution Score (MCS) to determine the Support score of each gene for each corresponding reaction. The support score is compared to the null distribution to compute a p-value.
