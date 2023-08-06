# population-genetic-analysis-of-alfalfa
This is the script used for analysis the population data of alfalfa

### step1, snp, indel and SV calling pipeline, including tetraploid variants calling

you can see the detail information from here: 
```
snp_and_indel_calling.sh
freebayes_variants_calling.sh
SV_calling.sh
```
### step2, population genetic analysis

This part includes all related content of population genetic-related analysis, including phylogenetic tree, admixture, PCA, IBD, introgression
```
phylogenetic tree: phylogenetic_tree.sh, phylogenetic_tree_plot.R, vcf2other.py
PCA and IBD: PCA_and_IBD.sh, PCA.R, IBD.R, lecture06_07_add_id.pl
admixture: admixture.sh, admixture.R
introgression: introgression.sh, F4-ratio.R, introgression_ABBA_R_plot.R

```
### step3, population dynamics analysis

This part includes all related content of population dynamics analysis, including fastsimcoal, 
```
fastsimcoal: fastsimcoal.sh, 2psfs.py


```

### step4, genotype and environmental factors association analysis


### step5, genetic offset analysis

