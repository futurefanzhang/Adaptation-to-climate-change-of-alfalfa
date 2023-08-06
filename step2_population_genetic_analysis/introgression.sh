##One: the overall introgression percent, Dsuite, F4-ratio (D-statistic) (https://github.com/millanek/Dsuite)
##step1, filter individuals
python /public/home/wangxu02/software/genomics_general/VCF_processing/parseVCF.py -i merge_712_allele2_meanDP3-62_miss0.2_addid.plink.vcf.vcf -o genomic_general_all_vcf.gz
plink --vcf merge_712_allele2_meanDP3-62_miss0.2_addid.plink.vcf --keep 296_accession.txt --allow-extra-chr --recode vcf --out merge_296_allele2_meanDP3-62_miss0.2_addid.plink --threads 1  ##only extract samples used for F4-ratio test
##step2, run Dsuite for F4-ratio test
Dsuite Dtrios  ./merge_296_allele2_meanDP3-62_miss0.2_addid.vcf.recode.vcf.vcf ./296_group.txt -o Dtrios_sativa_falcata_296_notree ##296_group.txt include two columns, individual name and group info
##The results will be in *.BBAA.txt file, please refer F4-ratio.R for plot

##Two: the introgression region were tested using ABBA test (https://github.com/simonhmartin/genomics_general)
python /public/home/wangxu02/software/genomics_general/VCF_processing/parseVCF.py -i merge_712_allele2_meanDP3-62_miss0.2_addid.plink.vcf.vcf -o genomic_general_all_vcf.gz 
python /public/home/wangxu02/software/genomics_general/ABBABABAwindows.py -g genomic_general_all_vcf.gz \
-f phased -o m_type0_type1_falcata_ABBA.csv -w 10000 -m 50 -s 10000 -P1 popA -P2 popB -P3 popC -O popD -T 5 --popsFile pops_m_type0_type1_falcata.txt --writeFailedWindows ##windows size 10kb, snp step size 50, step size 10kb
##the results will be in m_type0_type1_falcata_ABBA.csv, please refer introgression_ABBA_R_plot.R for plot
