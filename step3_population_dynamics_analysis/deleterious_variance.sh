##step 1, make database. in order to improve analysis speed, I sperate genome by each chromosome
seqkit grep -f chr1.txt zm4-all.fasta >  zm4_chr1.fa # only show chr1
gzip zm4_chr1.fa 
grep "chr1" Msa.H..gff3 >zm4_chr1.gff3
gffread zm4_chr1.gff3 -T -o zm4_chr1.gene.gtf
gzip zm4_chr1.gene.gtf #generate gtf file
mkdir zm4_chr1
cd zm4_chr1
mkdir gene-annotation-src
mkdir chr-src
mkdir dbSNP
mv zm4_chr1.fa.gz zm4_chr1/chr-src #move fa file to corresponding folder
mv zm4_chr1.gene.gtf.gz zm4_chr1/gene-annotation-src #move gtf to corresponding folder
perl make-SIFT-db-all.pl -config ./zm4_chr1/alfalfa_config_chr1.txt #generate alfalfa database, my personal server location:/data/home/zhangfan/basic_file/software/sift4g/scripts_to_build_SIFT_db/zm4_chr1

##step 2, move all chr*.gz, chr*regions, and chr*_SIFTDB_stats.txt to one folder
/public/agis/zhouyongfeng_group/caoshuo/tools/plink/plink --vcf merge_712_allele2_meanDP3-62_miss0.2_addid.plink.vcf --keep 267_accession.txt --allow-extra-chr --recode vcf --out merge_267_allele2_meanDP3-62_miss0.2_addid.plink --threads 1 ##挑选需要分析的个体
python /public/home/zhangtianhao/soft/superSFS/superSFS.py 1 10truncatula_outgroup.txt 11 merge_267_allele2_meanDP3-62_miss0.2_addid.plink.vcf merge_267_allele2_meanDP3-62_miss0.2_addid.plink_reverse.vcf  ###功能1，根据外类群反转vcf；10truncatula_outgroup.txt是外类群样本名，11代表10个外类群（20个单倍型）中超过50%单倍型不一致就翻转SNP，
java -Djava.awt.headless=true -jar /public/agis/zhouyongfeng_group/zhangfan02/arm_software/SIFT4G_Annotator.jar -c -i merge_267_allele2_meanDP3-62_miss0.2_addid.plink_reverse.vcf -d zm4_database/ -r ./results/ -t #输入vcf文件，-d就是SIFT的文件夹， -r是结果输出路径， -c代表命令行运行，-t表示从多个转录本提取注释信息

##step 3, extract deleterious SNP and change to R format
grep "DELETERIOUS" merge_267_allele2_meanDP3-62_miss0.2_addid.plink_reverse_SIFTpredictions.vcf > deleterious_vcf_file.vcf ##提取对应的vcf位点信息
head -n 100 merge_267_allele2_meanDP3-62_miss0.2_addid.plink_reverse_SIFTpredictions.vcf> head100.txt #把vcf的表头加上,（这块我用了三步，有更便捷方法请告诉我）
grep '^#' head100.txt >head_info.txt
cat head_info.txt deleterious_vcf_file.vcf > deleterious_vcf_file_all.vcf ##完整的有害变异vcf文件
sed -e "s/\r//g" deleterious_vcf_file_all.vcf > deleterious_vcf_file_all2.vcf ##remove ^M sign of vcf file
perl /public/agis/zhouyongfeng_group/zhangfan02/vcf_file/SIFT4G_zm4_database/change_vcf_format_genotype.pl deleterious_vcf_file_all2.vcf > deleterious_vcf_file_all_dosage.vcf ##把0/0,0/1,1/1转换成0,1,2格式
cut -f3,10- deleterious_vcf_file_all_dosage.vcf > R_use_format.vcf #去掉不需要的列
grep -v '^#' R_use_format.vcf >R_use_format1.vcf #去掉注释信息，生成模R需要用的文件

##step 4 calculate deleterious variance number, generate figure in R
##please refer deleterious_variance_figure.R
