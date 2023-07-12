##step1 filter raw reads by trimmomatic,remove adapter and low quality reads.
trimmomatic PE -threads 30 *.R1.fq.gz *.R2.fq.gz trimmed.sample1_R1.clean.fastq.gz trimmed.sample1_unpaired_R1.clean.fastq.gz trimmed.sample1_R2.clean.fastq.gz trimmed.sample1_unpaired_R2.clean.fastq.gz ILLUMINACLIP:/data/home/zhangfan/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:35:4:12:true LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:50 

##step2 mapping reads to reference genome using bwa-mem, and filter mapping results
bwa mem -M -t 40 /data/home/zhangfan/basic_file/reference/zm4-all.fasta trimmed.sample1_R1.clean.fastq.gz trimmed.sample1_R2.clean.fastq.gz -o sample1.sam
samtools view -@ 40 -bS -F 0x100 -q 10 sample1.sam| samtools sort -@ 40 -o sample1.sorted.bam ##filter reads with to filter multiple alignment (-F 0x100) and quality less than 10 (-q 10) and sort bam file

##step3 mark PCR duplicates


