##this script include the structure variance (SV) calling pipeline
##ref:https://github.com/shangshanzhizhe/Work_flow_of_population_genetics/blob/master/Work_flows/structure_variation.md (in Chinese)
##step1, call sv for every individual
  delly call --genome /public/agis/zhouyongfeng_group/caoshuo/data/zhangfan/index/zm4-all.fasta -o sample1.geno.bcf sample1.sorted.rmdup.bam ##use the bam file from snp_indel_calling.sh pipeline
##step2, merge bcf file of all individual
  delly merge -o all.sites.bcf sample1.geno.bcf sample2.geno.bcf sample3.geno.bcf
##step3, call genotype for every individual. It is called by the merged SV of all individual
  delly call --genome /public/agis/zhouyongfeng_group/caoshuo/data/zhangfan/index/zm4-all.fasta -v all.sites.bcf -o sample1.geno.bcf sample1.bam
##step4, merge all *.geno.bcf
  bcftools merge -m id -O b -o merge.final.bcf sample1.geno.bcf sample2.geno.bcf sample3.geno.bcf
##step5, generate vcf file
  bcftools index merge.final.bcf #build index for bcf file
  delly filter -f germline -o merge.final_filter.bcf merge.final.bcf #filter bcf file
  bcftools view merge.final_filter.bcf > merge.final.bcf.vcf #change bcf to vcf file
  python3 /public/agis/zhouyongfeng_group/zhangfan02/delly_results/filter.py merge.final.bcf.vcf #filter out IMPRECISE and lowQual SV
  vcftools --vcf merge.final.bcf.filter.vcf --max-missing 0.8 --recode --recode-INFO-all --out final_704ind_sv_miss0.2.vcf #filter missing rate of 20%

