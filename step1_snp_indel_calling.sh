#####this script is the pipeline of SNP and Indel calling, the annotation were added at the end of command

##step1 filter raw reads by trimmomatic,remove adapter and low quality reads.
  trimmomatic PE -threads 30 *.R1.fq.gz *.R2.fq.gz trimmed.sample1_R1.clean.fastq.gz trimmed.sample1_unpaired_R1.clean.fastq.gz trimmed.sample1_R2.clean.fastq.gz trimmed.sample1_unpaired_R2.clean.fastq.gz ILLUMINACLIP:/data/home/zhangfan/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:35:4:12:true LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:50 

##step2 mapping reads to reference genome using bwa-mem, and filter mapping results
  bwa mem -M -t 40 /data/home/zhangfan/basic_file/reference/zm4-all.fasta trimmed.sample1_R1.clean.fastq.gz trimmed.sample1_R2.clean.fastq.gz -o sample1.sam
  samtools view -@ 40 -bS -F 0x100 -q 10 sample1.sam| samtools sort -@ 40 -o sample1.sorted.bam ##filter reads with to filter multiple alignment (-F 0x100) and quality less than 10 (-q 10) and sort bam file

##step3 mark PCR duplicates
  java -jar /data/home/zhangfan/basic_file/software/picard.jar MarkDuplicates INPUT= sample1.sorted.bam OUTPUT= sample1.sorted.rmdup.bam METRICS_FILE= duplicateMatrix REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=tmp/ ##remove duplicated reads

##step4 detected SNP and Indel by GATK
  ~/conda/gatk-4.2.0.0/gatk --java-options "-XX:ConcGCThreads=8 -XX:+UseSerialGC -Xmx35g" AddOrReplaceReadGroups -I sample1.sorted.rmdup.bam -O sample1.addhead.sorted.bam -ID sample1 -LB sample1 -PL illumina -PU samplePU -SM sample1 ##add file head
  samtools index sample1.addhead.sorted.bam  ##add index
  ~/conda/gatk-4.2.0.0/gatk --java-options "-XX:ConcGCThreads=8 -XX:+UseSerialGC -Xmx35g" HaplotypeCaller --native-pair-hmm-threads 8 -R /data/home/zhangfan/basic_file/reference/zm4-all.fasta -I sample1.addhead.sorted.bam -O sample1.gvcf --emit-ref-confidence GVCF -stand-call-conf 30  ##generate gvcf file
##step4.1
  repeat step 1-4 for each sample
  
##step5 combine all samples gvcf file
  mv sample*.gvcf* ./all_gvcf/ #move all gvcf file to one folder
  ls *.gvcf > sample_map.txt  ##generate sample map file, and you need to add sample name info in the first column (I preferred Excel)
  mkdir temp
  gatk --java-options "-XX:ConcGCThreads=8 -XX:+UseSerialGC -Xmx100g -Xms100g" GenomicsDBImport --tmp-dir ./temp/ -L Chr1 --batch-size 50 --reader-threads 8 --sample-name-map sample_map.txt --genomicsdb-workspace-path ./database1 ##generate database for chromosome 1, you need to do this for each chromosome
  gatk  GenotypeGVCFs -L Chr1 -V gendb://database1 -O Chr1.vcf -R /data/home/zhangfan/basic_file/reference/zm4-all.fasta ##generate raw vcf file for chromosome 1
	grep '^#' Chr1.vcf >merge.vcf #extract annotation info from Chr1.vcf 
	grep -h -v '^#' Chr*.vcf >>merge.vcf  #merge all chromosome results
	gatk --java-options "-Xmx60g -Xms60g" SelectVariants -select-type SNP -V merge.vcf -O merge.snp.vcf && gatk VariantFiltration -V merge.snp.vcf -O SNP_filter.vcf --filter-expression "QD<2.0 || FS >60.0 || MQRankSum <- 12.5 || ReadPosRankSum <- 8.0 || SOR >3.0 || MQ<40.0" --filter-name "Fail" && grep -v "Fail" SNP_filter.vcf >SNP_filter_out.vcf
	gatk --java-options "-Xmx60g -Xms60g" SelectVariants -select-type INDEL -V merge.vcf -O merge.indel.vcf && gatk VariantFiltration -V merge.indel.vcf -O Indel_filter.vcf --filter-expression "QD<2.0 || FS >200.0 || MQRankSum <- 12.5 || ReadPosRankSum <- 8.0 || SOR >10.0" --filter-name "Fail" && grep -v "Fail" Indel_filter.vcf >Indel_filter_out.vcf
##step6 filter SNP and indel
	vcftools --vcf SNP_filter_out.vcf --max-missing 0.8 --min-meanDP 3 --max-meanDP 62 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out all_population_SNP_filter_missing0.2_meanDP2-50_allele2.vcf
  vcftools --vcf Indel_filter_out.vcf --max-missing 0.8 --min-meanDP 3 --max-meanDP 62 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out all_population_Indel_filter_missing0.2_meanDP2-50_allele2.vcf
