module load bedtools
bedtools intersect -wa -a ASE_for_chromatin_OL.bed -b /wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/open_chromatin_regions_for_OL.bed> ASE_open_chromatin_ol.bed
bedtools intersect -wa -a cASE_for_chromatin_OL.bed -b /wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/open_chromatin_regions_for_OL.bed> cASE_open_chromatin_ol.bed
