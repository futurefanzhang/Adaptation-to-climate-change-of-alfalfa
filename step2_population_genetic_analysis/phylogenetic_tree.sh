###step1, generate population phylogenetic tree
plink --allow-extra-chr --threads 2 --vcf all_population_SNP_filter_missing0.2_meanDP2-50_allele2.vcf  --recode vcf --out all_712_iqtree_plink_vcf  ##use plink generate simple format vcf file
python vcf2other.py -i all_712_iqtree_plink_vcf.vcf -o SRR12686087_SRR12686087 -f #-o set out group, -f  Write a FASTA matrix

