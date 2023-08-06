##step 1, filter SNP by LD
perl /public/agis/zhouyongfeng_group/zhangfan02/vcf_file/lecture06_07_add_id.pl merge_712_allele2_meanDP3-62_miss0.2.vcf.recode.vcf merge_712_allele2_meanDP3-62_miss0.2_addid.vcf.recode.vcf
plink --vcf merge_712_allele2_meanDP3-62_miss0.2_addid.vcf.recode.vcf --indep-pairwise 100 50 0.2 --out m.714.snp.ld --allow-extra-chr --make-bed
plink --bfile m.712.snp.ld --extract m.712.snp.ld.prune.in --out all_712_iqtree_plink_vcf --recode vcf --allow-extra-chr
plink --vcf all_712_iqtree_plink_vcf.vcf --keep 702_accessions.txt --allow-extra-chr --threads 1 --recode vcf --out all_702_iqtree_plink_vcf

##step 2, PCA analysis
plink --allow-extra-chr --threads 2 --vcf all_702_iqtree_plink_vcf.vcf --pca 20 --out 702_accession_pca
#please refer PCA.R for plot 

##step 3, IBD analysis
plink --vcf all_702_iqtree_plink_vcf.vcf --const-fid --allow-extra-chr --genome --out IBD_702_result --threads 2
#please refer IBD.R for plot 
