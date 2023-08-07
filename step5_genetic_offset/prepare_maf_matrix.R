##step 1, get maf matrix
library(ENMwizard)#version 0.4.2
obs_data397 <- read.table(file =  "397sample_info.txt",header = T,sep = "\t") ##sampling points
occs=obs_data397[,c("taxa","Longitude","Latitude")]
colnames(occs)=c("species","lon","lat")
occs$species=c("alfalfa")
occs <- occs[!duplicated(occs),] ##remove duplicated location
spp.occ.list <- list(alfalfa = occs) ##change table to list,name as alfalfa
thinned.dataset.batch <- thin_b(loc.data.lst = spp.occ.list) ##get represent sample group,from 397 to 116

##using plink calculate maf info
plink --vcf merge.vcf --make-bed --allow-extra-chr --out merge_order --threads 1 ##merge.vcf include all significant associated variants (SNP, Indel and SV)
plink --bfile merge_order --freq --allow-extra-chr --within 397group.txt --out 116_position_maf_freq.txt ##get maf info of 116 group

##generate maf matrix in R
mydata=read.table("116_position_maf_freq.txt.frq.strat",header = T,sep = "")
df=mydata[,c("SNP","CLST","MAF")]
library(reshape2)
df <- dcast(df,SNP~CLST,fun.aggregate = list,value.var='MAF')
df <- t(df)
write.table(df,"116_position_MAF_freq_R_format.txt",sep = "\t", row.names =TRUE, col.names =FALSE, quote =FALSE)
