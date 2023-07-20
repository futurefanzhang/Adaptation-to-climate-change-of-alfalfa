###iqtree, generate population phylogenetic tree
plink --allow-extra-chr --threads 2 --vcf all_population_SNP_filter_missing0.2_meanDP2-50_allele2.vcf  --recode vcf --out all_712_iqtree_plink_vcf  ##use plink generate simple format vcf file
python vcf2other.py -i all_712_iqtree_plink_vcf.vcf -o SRR12686087_SRR12686087 -f #change vcf file to fasta format, -o set out group, -f  Write a FASTA matrix
iqtree -s all_712_iqtree_plink_vcf.min4.phy -nt 96 -mem 180G -m MF -pre test_sample_tree.nwk -st DNA #find best model, find which model has minimum AIC value
iqtree -s all_712_iqtree_plink_vcf.min4.phy -nt 96 -mem 180G -m GTR+F+G4 -bb 1000 -pre merge_712_sativa.nwk -st DNA #built phylogenetic tree, bootstrap value of 1000, my model is GTR+F+G4

##iqtree results can be plotted in R, please see phylogenetic_tree_plot.R
