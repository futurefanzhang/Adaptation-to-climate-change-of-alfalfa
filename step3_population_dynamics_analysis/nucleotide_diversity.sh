##this command will calculate nucleotide diversity using non-variant sites.
gatk  GenotypeGVCFs -L Chr1 -V gendb://database1 -O Chr1.vcf -R /data/home/zhangfan/basic_file/reference/zm4-all.fasta -all-sites
vcftools --vcf Chr1.vcf --remove-indels --max-missing 0.8 --min-meanDP 3 --max-meanDP 62 --keep 712_sample.txt --recode --stdout | bgzip > my_filtered_chr1_vcf.vcf.gz
tabix my_filtered_chr1_vcf.vcf.gz
pixy --stats pi fst dxy --vcf my_filtered_chr1_vcf.vcf.gz --populations group_712_pixy.txt --window_size 10000 --n_cores 1 --chromosomes 'chr1' --output_prefix chr1_pixy ##run for each chromosome

###these command run in R
##one species
pi=read.table("chr1_pixy_pi.txt",sep="\t",header=T)
pi=na.omit(pi)
pi_species=pi[pi$pop=="M.truncatula",]
pi_species=pi_species[pi_species$count_diffs>0,]
mean_pi=sum(pi_species$count_diffs)/sum(pi_species$count_comparisons)
mean_pi
##all species
all_info=data.frame(matrix(nrow = 0, ncol = 0))
for (i in 1:8){
pi=read.table(paste("chr",i,"_pixy_all_alfalfa_pi.txt",sep = ""),sep="\t",header=T)
pi=na.omit(pi)
pi=pi[pi$count_diffs>0,]
all_info=rbind(all_info,pi)
}
all_info_use=all_info[all_info$pop!="Outgroup",]
library(ggplot2)
p1<-ggplot(all_info_use, aes(x=pop, y=avg_pi, fill=pop))+labs(x="Species", y="Nucleotide diversity (pi)") +
  geom_boxplot(width=0.8)+#ylim(0.7,1.2)+
  theme_classic()+theme(
    legend.position="none",legend.title = element_text(color = "black", size = 10), legend.text = element_text(color = "black", size = 10),axis.text = element_text(color="black",size=10),axis.title.x=element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size=10),axis.title.y=element_text(size=10))+
scale_fill_manual( values = c("#FABB2E","#19B700","#ee0000","#BF3EFF","#42d4f4","#A3A500","#8794FF","#A8422D","#00868B"),
                   breaks=c("M.sativa", "M.caerulea","M.varia","M.falcata.tetraploid","M.falcata.diploid", "M.ruthenica","M.archiducis-nicolai","M.truncatula","M.group"),
                   labels=c("M.sativa", "M.caerulea","M.varia","M.falcata.tetraploid","M.falcata.diploid", "M.ruthenica","M.archiducis-nicolai","M.truncatula","M.group"))

ggsave(plot=p1,"M.sativa_species.genetic_diversity.pdf", width = 15, height = 6)
