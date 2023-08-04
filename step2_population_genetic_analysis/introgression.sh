##the overall introgression percent, Dsuite, F4-ratio (D-statistic)
##step1, filter individuals
python /public/home/wangxu02/software/genomics_general/VCF_processing/parseVCF.py -i merge_712_allele2_meanDP3-62_miss0.2_addid.plink.vcf.vcf -o genomic_general_all_vcf.gz
plink --vcf merge_712_allele2_meanDP3-62_miss0.2_addid.plink.vcf --keep 296_accession.txt --allow-extra-chr --recode vcf --out merge_296_allele2_meanDP3-62_miss0.2_addid.plink --threads 1  ##only extract samples used for F4-ratio test
Dsuite Dtrios  ./merge_296_allele2_meanDP3-62_miss0.2_addid.vcf.recode.vcf.vcf ./296_group.txt -o Dtrios_sativa_falcata_extract_notree ##296_group.txt include two columns, individual name and group info
##The results will be in *.BBAA.txt file, please refer F4-ratio.R for plot
