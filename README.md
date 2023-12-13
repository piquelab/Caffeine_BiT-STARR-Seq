# Characterization of caffeine response regulatory variants in vascular endothelial cells
Contains scripts used in "Characterization of caffeine response regulatory variants in vascular endothelial cells" (doi: https://doi.org/10.1101/2022.11.22.517533)

## Data Analysis Pipeline
The script used for alignment can be found at https://github.com/cakalita/BiT-STARR-seq/blob/main/alignandprep.sh  
`DESeq_script.R` is the script for differential activity analysis.  
To complete ASE and cASE analyses, we first run `bitstarr_ASE_analysis.R` to obtain individual replicate-level results, and then run `manual_calc_meta_analysis_using_quasar_se.R` to complete the meta-analysis. Additionally, we filtered our results to remove variants with less data using `filter_meta_analysis.R`.  
To do the enrichment analyses for open chromatin and artery eQTLs, we use the following scripts respectively: `open_chromatin_OL.R` and `chromatin_OL.sh`, `GTEx_Overlap.R`.  
For the motif analysis, `PWM_Scan.sh` is the script used to scan for motifs and `motif_proptest.R` is the script used to complete the statistical test.  
`intact_OL.R` is the script used to view the overlap between our results and the fine-mapped artery eQTLs, and `cASE_statistics_working_back_from_INTACT.R` is the script used to link our INTACT results with our cASE results.  
