##this is the pipeline of tetraploid genotype calling
##bam_file.txt is the absolute path of bam file
##you can check the detail info of command from here:https://github.com/freebayes/freebayes
freebayes -f /public/agis/zhouyongfeng_group/zhangfan02/index_bwa/zm4-all.fasta -p 4 -g 100000  -L bam_file.txt --min-mapping-quality 10 --min-alternate-count 4 --min-coverage 4  --use-best-n-alleles 2 --min-base-quality 10 --region chr1:10000000-11000000 > all_variants_tetraploid_chr1_10-11Mb.vcf ##only run one region at here. If you want to quickly get the results, you can run all regions at the same time.
