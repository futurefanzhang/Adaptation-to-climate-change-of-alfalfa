##step 1, generate fst results: vcftools --vcf merge_712_allele2_meanDP3-62_miss0.2_addid.vcf.recode.vcf --weir-fst-pop East_Asia.txt --weir-fst-pop Africa.txt --out East_asia_vs_africa_Fst --fst-window-size 10000

##step 2, using R do MDS analysis
sativa_falcata<-read.table("East_asia_vs_africa_Fst.windowed.weir.fst",header=TRUE,sep="\t")
all_pop_use=sativa_falcata[sativa_falcata$N_VARIANTS>=10 ,] #only use SNP number more than 10
all_pop_use$MEAN_FST[all_pop_use$MEAN_FST<0] = 0 #change less than zero as zero
#write.csv(all_pop_use,"fst_10bk_use_sativa_falcata.csv")
mean(all_pop_use$MEAN_FST) ##calculate mean Fst for each pairwise group
##step3, do MDS analysis
Fst_data=read.table("Fst_matrix.txt",header = T,row.names = 1,sep = "\t") ##pairwise Fst matrix
MDS=cmdscale(Fst_data,k=2,eig=TRUE)
pdf("MDS_by_Fst.pdf",paper = "a4")
par(omi=c(0.1,0.1,0.1,0.1))
plot(x=MDS$points[,1], y=MDS$points[,2], xlab="MDS1", ylab="MDS2",
     main="Multidimensional Scaling Results", pch=16,cex=1.5,col=c("#ee0000","#BF3EFF","#19B700","#8794FF","#FABB2E","#42d4f4"))
#add row names of data frame as labels
text(x=MDS$points[,1], y=MDS$points[,2], labels=row.names(MDS$points),pos=4)
dev.off()
