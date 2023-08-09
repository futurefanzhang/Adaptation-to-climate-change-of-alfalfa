# Adaptation-to-climate-change-of-alfalfa
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

This part includes all related content of population dynamics analysis, including fastsimcoal, candidate gene analysis, MDS, deleterious variance and nucleotide diversity
```
fastsimcoal: fastsimcoal.sh, 2psfs.py
fst candidate gene analysis: fst_candidate_gene.R, candidate_gene.sh, candidate_gene_filter.R
MYB5 gene analysis: MYB5_gene.sh, MYB5_analysis.R, local_manhattan.sh
MDS analysis: MDS_analysis.R
Deleterious variance: deleterious_variance.sh, deleterious_variance_figure.R, change_vcf_format_genotype.pl
nucleotide diversity analysis: nucleotide_diversity.sh
```

### step4, genotype and environmental factors association analysis

This part includes all related content of population dynamics analysis, including fastsimcoal, candidate gene analysis, MDS, deleterious variance and nucleotide diversity
```
LFMM (latent factor mixed model ) analysis: LFMM_analysis.R, LFMM_results_figure.R, match_indel_position.pl
PBS (population branch statistic ) analysis: PBS_analysis.R

```

### step5, genetic offset analysis

This part includes offset estimation command and plot offset in map
```
step5.1 prepare_maf_matrix.R
step5.2 evaluate_representative_biofactor.R
step5.3 offset_estimation.R
step5.4 future_forward_offset_by_migration.R, local_forward_reverse_offset.R
species distribution region based on bioclimate variables: species_distribution.R
```

### if you have any questions, please send email to me: zhangfan06@caas.cn or zfan887@gmail.com
