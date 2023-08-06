##candidate gene MYB5
intersectBed -a merge_tetraploid_q60_GT_dosage.vcf -b MYB5_10k_region.bed -header > MYB5_chr7_10k_region_freebayes.vcf
cut -f1,2,10- MYB5_chr7_10k_region_freebayes.vcf > freebayes_R_use_format_chr7_MYB5.vcf
##please refer MYB5_analysis.R
