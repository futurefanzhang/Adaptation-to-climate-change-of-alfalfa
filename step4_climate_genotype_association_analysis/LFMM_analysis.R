##ref:https://github.com/jingwanglab/Populus_genomic_prediction_climate_vulnerability

plink --vcf final_702ind_sv_miss0.2_maf0.05.vcf.recode.vcf --keep 399_alfalfa_accession.txt --maf 0.1 --allow-extra-chr --recode vcf --out all_399_sv_plink_vcf --threads 2 ##extract vcf for SNP, indel and SV
java -Xmx10g -jar /public/agis/zhouyongfeng_group/zhangfan02/vcf_file/beagle.22Jul22.46e.jar nthreads=2 gt=all_399_sv_plink_vcf.vcf out=all_399_sv_plink.beagle.vcf
gzip -d all_399_sv_plink.beagle.vcf.vcf.gz
plink --vcf all_399_sv_plink.beagle.vcf.vcf --allow-extra-chr --recode  --out all_399_sv_plink --threads 2

##run in R
library(LEA)
pedfile = ped2lfmm("all_397_sv_plink.ped") #inport ped file
BIO_file="LFMM_BIO1.env" ##input bioclim data info (phenotype), no header
project=lfmm(input.file=pedfile,environment.file=BIO_file,CPU=6,K = 3,repetitions = 5, project = "new") ##if you want to speed up, add more CPU
p = lfmm.pvalues(project, K = 3,  d = 1, run = 1)
pvalues1 = p$pvalues
p = lfmm.pvalues(project, K = 3,  d = 1, run = 2)
pvalues2 = p$pvalues
p = lfmm.pvalues(project, K = 3,  d = 1, run = 3)
pvalues3 = p$pvalues
p = lfmm.pvalues(project, K = 3,  d = 1, run = 4)
pvalues4 = p$pvalues
p = lfmm.pvalues(project, K = 3,  d = 1, run = 5)
pvalues5 = p$pvalues
mean_p=cbind(pvalues1,pvalues2,pvalues3,pvalues4,pvalues5)
mean_p1=as.data.frame(apply(mean_p,1,mean))
name <- read.table("all_397_sv_plink.map",header = F,sep = "")
mean_p1$name=name$V2
mean_p1=mean_p1[order(mean_p1$`apply(mean_p, 1, mean)`),]
mean_p1$FDR=p.adjust(mean_p1$`apply(mean_p, 1, mean)`, method = "fdr")
colnames(mean_p1)=c("mean_P_value","SNP","FDR")
write.table(mean_p1,paste("result_lfmm/lfmm_pvalue_maf0.1",BIO_file,"txt",sep="."),row.names = F,quote=F,sep="\t") #export all results
mean_p2=mean_p1[mean_p1$FDR<0.05,]
write.table(mean_p2,paste("result_lfmm/lfmm_0.05fdr_maf0.1",BIO_file,"txt",sep="."),row.names = F,quote=F,sep="\t") #export fdr more than 0.05
