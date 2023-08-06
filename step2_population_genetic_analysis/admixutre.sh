perl /public/agis/zhouyongfeng_group/zhangfan02/vcf_file/lecture06_07_add_id.pl merge_712_allele2_meanDP3-62_miss0.2.vcf.recode.vcf merge_712_allele2_meanDP3-62_miss0.2_addid.vcf.recode.vcf
plink --vcf merge_712_allele2_meanDP3-62_miss0.2_addid.vcf.recode.vcf --indep-pairwise 100 50 0.2 --out m.712.snp.ld --allow-extra-chr --make-bed
plink --bfile m.712.snp.ld --extract m.712.snp.ld.prune.in --out all_712_admixture --recode 12 --allow-extra-chr
admixture --cv prunData.ped 2 >>log.txt ##k=2
admixture --cv prunData.ped 3 >>log.txt ##k=3, same for other k from 4 to 10
grep "CV error" log.txt ##chick CV error for different k value
