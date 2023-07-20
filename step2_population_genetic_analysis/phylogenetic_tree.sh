###step1, generate population phylogenetic tree
plink --allow-extra-chr --threads 2 --vcf 243_migrate_accession_reorder.vcf --chr 1  --recode vcf --out 243_migrate_accession_reorder_chr1
python vcf2other.py -i all_712_iqtree_plink_vcf.vcf -o SRR12686087_SRR12686087 -f #-o是设置外类群
