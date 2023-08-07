##step 1, add header, these command in linux
cut -f1-3  all_399_indel_plink_addid.vcf > all_399_indel_position.txt ##extract pos info
perl /public/agis/zhouyongfeng_group/zhangfan02/vcf_file/match_indel_position.pl all_399_indel_position.txt ./result_lfmm/lfmm_pvalue_maf0.1_new.LFMM_BIO2.env.txt ./result_lfmm/add_pos_indel_maf0.1_bio2.txt

##step 2 figure in R
library(dplyr)
library(ggplot2)
library(ggrepel)
library(qqman)

##raw figure
BIO="1"
my_indel=read.table(paste("add_pos_indel_maf0.1_bio",BIO,".txt",sep=""),header = F,sep = "\t")
colnames(my_indel)=c("CHR","BP","SNP","P","FDR")
min=min(my_indel[my_indel$P>0,][,4]) ##find minimum p value
my_indel$P[my_indel$P==0] <- min ##change zero p value to minimum p value
mydata1=my_indel[,c(3,1,2,4)]
Bonferroni_sig=-log10(0.05/nrow(mydata1)) ##set bonferroni significant level
my_indel2=my_indel[my_indel$FDR<=0.05,]
fdr_sig=-log10(max(my_indel2$P)) ##set fdr threshold
png(paste("LFMM_manhattan",BIO,"png",sep="."),width=15,height=4,units='in',res=500)
par(mar=c(5,5,2,2))
manhattan(mydata1,main=paste("BIO",BIO,sep=" "),suggestiveline =fdr_sig,genomewideline =Bonferroni_sig,col=c("grey60","#4197d8"),cex.axis = 1.2,cex.lab=1.4)
dev.off()

##filter significant variance
mydata=read.table("add_pos_snp_bio1.txt",header = F,sep = "\t")
colnames(mydata)=c("CHR","BP","SNP","P","FDR")
              ########raw figure, add highlight point############
              library(CMplot)
              mydata1<-mydata[,c(3,1,2,4)]
              colnames(mydata1) <- c("ID","CHR", "BP","P")
              
              my_snp=mydata[mydata$FDR<=0.05,]
              fdr_sig=-log10(max(my_snp$P)) ##set fdr threshold
              mydata1$P=-log10(mydata1$P)
              SNPs <- mydata1$ID[mydata1$CHR=="8" & 
                                   mydata1$BP<=(28347599)&mydata1$BP>=(28347599)]
              CMplot(mydata1,plot.type = "m",cex=c(0.5,0.8,1), col=c("#B2B0E1","#AEE3AD"),highlight=SNPs,highlight.col="red",
                     width=14,height=5,signal.cex=0.8,LOG10=F,threshold = fdr_sig,chr.den.col = NULL, file="jpg",memo="test manhattan", 
                     ylab="Fst", main="Cross vs Self", file.output = TRUE,verbose = TRUE)
              ##finish##########
use_data=mydata[mydata$FDR<=0.05,] #filter FDR<=0.05
use_data$BIO="BIO1"
for (i in 2:19){
  BIO=i
  my_snp=read.table(paste("add_pos_snp_bio",BIO,".txt",sep=""),header = F,sep = "\t")
  colnames(my_snp)=c("CHR","BP","SNP","P","FDR")
  my_snp=my_snp[my_snp$FDR<=0.05,]
  if(nrow(my_snp)!=0){
  my_snp$BIO=paste("BIO",BIO,sep="")
  use_data=rbind(use_data,my_snp)}
}
write.table(use_data,"significant_snp_results_bio1-19.txt",sep="\t",row.names = F,quote = F)

##step 3, generate combine manhattan plot####
##manhattan plot, combine figure SNP, Indel and SV in R
for (i in 2:19){
BIO=i
my_snp=read.table(paste("add_pos_snp_bio",BIO,".txt",sep=""),header = F,sep = "\t")
colnames(my_snp)=c("CHR","BP","SNP","P","FDR")
mydata1=my_snp[,c(3,1,2,4)]
Bonferroni_sig=-log10(0.05/nrow(mydata1)) ##set bonferroni significant level
my_snp=my_snp[my_snp$FDR<=0.05,]
fdr_sig=-log10(max(my_snp$P)) ##set fdr threshold
mydata1$Type="SNP"

my_indel=read.table(paste("add_pos_indel_bio",BIO,".txt",sep=""),header = F,sep = "\t")
colnames(my_indel)=c("CHR","BP","SNP","P","FDR")
mydata2=my_indel[,c(3,1,2,4)]
mydata2$Type="Indel"

my_sv=read.table(paste("add_pos_sv_bio",BIO,".txt",sep=""),header = F,sep = "\t")
colnames(my_sv)=c("CHR","BP","SNP","P","FDR")
mydata3=my_sv[,c(3,1,2,4)]
mydata3$Type="SV"

gwas_data=rbind(mydata1,mydata2,mydata3)
gwas_data=gwas_data[order(gwas_data$CHR,gwas_data$BP),]

data_cum <- gwas_data %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHR, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)
xmin=min(gwas_data$bp_cum)
xmax=max(gwas_data$bp_cum)

#####ADD point##
add_gene=read.table("add_gene_BIO15.txt",header = F,sep = "\t")  ##set by gene annotation info
add_gene$SNP<-paste(add_gene$V1,add_gene$V2, sep="__")
highlight_point=gwas_data[(gwas_data$SNP %in% add_gene$SNP),]
order=match(highlight_point$SNP, add_gene$SNP)
add_gene=add_gene[order,] #change order
highlight_point$gene=add_gene$V4
####

gwas_data$Type =  factor(gwas_data$Type, levels =c("SNP","Indel","SV")) #set order of Type,SNP, point, Indel, triangel SV box
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P), group=Type,
                                  color = as.factor(CHR)),size=0.8) +
  geom_hline(yintercept = c(fdr_sig, Bonferroni_sig), color = c('blue', 'red'), linetype = c("dashed", "dashed")) + 
  geom_point(aes(shape=Type)) +
  scale_shape_manual(values=c(1,2,15))+
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
  #scale_size_manual(values=c(2,3,4))+
  ###############add point#########
  geom_point(data=highlight_point, 
             aes(x = bp_cum, y =-log10(P), color = as.factor(CHR)),color='red')+
  geom_text_repel(data=highlight_point,
                  aes(label = gene),#box.padding = unit(0.35, "lines"),
                  size = 3, #基因字体大小
                  #nudge_y      = 0.1, #箭头的长度
                  #direction    = "x",
                  #angle        = 0,
                  #hjust        = 0.15,
                  #segment.size = 0.3,
                  #max.iter = 1e4, max.time = 1,
                  #arrow = arrow(length = unit(0.02, "npc"),type = "closed"),
                  color="red"
  )+
  ###########point end#########
  scale_x_continuous(label = axis_set$CHR,limits = c(xmin,xmax),expand = expansion(0), breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(0), limits = c(0, ylim),breaks = seq(0,ylim,3)) +
  #scale_color_manual(values = rep(c("grey60","#4197d8"), unique(length(axis_set$CHR)))) +
  scale_color_manual(values = rep(c("#c4e8a8","#72B33D"), unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome", title = paste("BIO",BIO,sep=""),
       y = expression(-log[10](italic(P)))) + 
  theme_classic() +
  theme( 
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 15, vjust = 0.5),
    axis.text.y = element_text(size = 15)
  )

ggsave(plot=manhplot,paste("LFMM_manhattan_merge_SNP_Indel_SV_BIO_test",BIO,".png",sep="") ,width=7.5,height=4,units='in',dpi=500)
}
