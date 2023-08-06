data3=read.table("freebayes_R_use_format_chr7_MYB5.vcf",header =T,sep = "\t")
data4=data3[,-c(1,2)]
group=read.table("group_sample_info.txt",header = T,sep="\t")
gene_info=read.table("MYB5_gene_structure.txt",header = T,sep="\t")
##extract group info
m.falcata=group[group$group=="M.falcata.tetraploid",]
snp_falcata=data4[,colnames(data4) %in% m.falcata$taxa]
m.sativa=group[group$group=="M.sativa",]
snp_sativa=data4[,colnames(data4) %in% m.sativa$taxa]

###m.falcata
data4=snp_falcata
allele0=rowSums(data4 == "0",na.rm = T)
allele1=rowSums(data4 == "1",na.rm = T)
allele2=rowSums(data4 == "2",na.rm = T)
allele3=rowSums(data4 == "3",na.rm = T)
allele4=rowSums(data4 == "4",na.rm = T)
data_f=as.data.frame(cbind(data3$POS, allele0,allele1,allele2,allele3,allele4))
colnames(data_f)=c("POS","Allele0","Allele1","Allele2","Allele3","Allele4")
data_f$sum=data_f$Allele0+data_f$Allele1+data_f$Allele2+data_f$Allele3+data_f$Allele4
data_f$Allele0_ratio=data_f$Allele0/data_f$sum
data_f$Allele1_ratio=data_f$Allele1/data_f$sum
data_f$Allele2_ratio=data_f$Allele2/data_f$sum
data_f$Allele3_ratio=data_f$Allele3/data_f$sum
data_f$Allele4_ratio=data_f$Allele4/data_f$sum

##m.sativa
data4=snp_sativa
allele0=rowSums(data4 == "0/0/0/0",na.rm = T)
allele1=rowSums(data4 == "0/0/0/1",na.rm = T)
allele2=rowSums(data4 == "0/0/1/1",na.rm = T)
allele3=rowSums(data4 == "0/1/1/1",na.rm = T)
allele4=rowSums(data4 == "1/1/1/1",na.rm = T)
data_s=as.data.frame(cbind(data3$POS, allele0,allele1,allele2,allele3,allele4))
colnames(data_s)=c("POS","Allele0","Allele1","Allele2","Allele3","Allele4")
data_s$sum=data_s$Allele0+data_s$Allele1+data_s$Allele2+data_s$Allele3+data_s$Allele4
data_s$Allele0_ratio=data_s$Allele0/data_s$sum
data_s$Allele1_ratio=data_s$Allele1/data_s$sum
data_s$Allele2_ratio=data_s$Allele2/data_s$sum
data_s$Allele3_ratio=data_s$Allele3/data_s$sum
data_s$Allele4_ratio=data_s$Allele4/data_s$sum

#barplot,sativa
data1=t(data_s[,c(8:12)])
data_percentage <- apply(data1, 2, function(x){x*100})
library(RColorBrewer)
coul <- brewer.pal(5, "Pastel2")
pdf(file="sativa_Allele_frequency_information_order.pdf",width = 16,height=4)
barplot(data_percentage, col=coul ,width = 0.1,cex.names=0.6,space = 0,names.arg=gene_info$order, #xaxt='n',
        border="white", xlab="",ylab="M.sativa allele Percent",legend = TRUE,args.legend = list(x = "topright",inset = c(-0.02, -0.18)),ylim = c(0,100))
dev.off()

##falcata
data1=t(data_f[,c(8:12)])
data_percentage <- apply(data1, 2, function(x){x*100})
colnames(data_percentage)=gene_info$order
pdf(file="falcata_Allele_frequency_information.pdf",width = 16,height=4)
##x-axis,names.arg=data3$POS/1000000,
barplot(data_percentage,col=coul ,width = 0.3,cex.names=0.4,space = 0, names.arg=gene_info$order,border="white", xaxt='n',
        xlab="SNP",ylab="M.falcata allele Percent",legend = TRUE,args.legend = list(x = "topright",inset = c(-0.02, -0.18)),ylim = c(0,100))
dev.off()

#association test，关联分析，可分别对SNP，Indel和SV分析
geno=read.table("702_association_geno_snp.txt",header = T,row.names = 1,sep = "\t")
pheno=read.table("457_association_pheno.txt",header = T,sep = "\t")
geno=geno[rownames(geno) %in% pheno$ID,]
order=match(rownames(geno),pheno$ID)
pheno=pheno[order,]
x=data.frame(matrix(nrow = ncol(geno), ncol = 2))
for (i in 1:ncol(geno)){
x[i,1]=colnames(geno)[i]
x[i,2]=cor.test(geno[,i],pheno$flower_color,na.action = "na.exclude")$p.value
}
colnames(x)=c('POS',"pvalue")
write.table(x,"457_association_snp_result.txt",row.names = F,quote = F,sep = "\t") ##can be used for local manhattan, local_manhattan.sh

##plot
data=read.table("457_association_snp_indel_sv_results.txt",header = T,sep = "\t")
data$Type =  factor(data$Type, levels =c("SNP","Indel","SV"))
data=data[data$POS>=101700479 & data$POS<=101705734,]
manhplot <- ggplot(data, aes(x = POS/1000000, y = -log10(pvalue), group=Type,
                                  color =Type),size=0.8) +
  geom_hline(yintercept = c(19.4), color = c('blue'), linetype = c("dashed")) + ##Bonferroni test 0.05/13460992=3.714436e-09
  geom_point(aes(shape=Type)) +
  #scale_x_continuous(limits = c(101701299, 101705347),breaks = seq(101701299, 101705347,1000)) +  #可设置基因位置信息
  scale_color_manual(values = c("#A3A500","#FABB2E","#8794FF")) +
  labs(x = "Chromosome (Mb)", title = "Regional association test",
       y = expression(-log[10](italic(P)))) +
  geom_segment(aes(x = 101702479/1000000, y = 0, xend = 101703734/1000000, yend = 0),color="red",size=2)+
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black",size = 15, vjust = 0.5),
    axis.text.y = element_text(color="black",size = 15),
    axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
  )
ggsave(plot=manhplot,"457_MYB5_manhattan_2kb_merge_SNP_Indel_SV.pdf" ,width=5,height=4)

##candidate SV haplotype
mydata=read.table("genotype101701543_SV.txt",header = T,sep = "\t")
mydata1=mydata[mydata$group=="M.varia",]
group_mean <- aggregate(mydata1$ID, by=list(mydata1$X101701543), FUN=length)
group_mean$x=group_mean$x/sum(group_mean$x)
group_mean ##used for next step

##pie plot，饼图
pie<- ggplot(data1, aes(x="", y=as.numeric(M.group), fill=type))+#set basic figure，y信息每个物种都要做
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ #change to pie plot
  scale_fill_manual(values= c("#FFC000","#A9D18E","#D7A1F9"),name = "Genotype code",
                    labels = c("Reference", "Heterozygotye","Alternative"))+ #change different color for each piece,another color:Zissou1
  labs(x="",y="M.group")+
  theme_void(base_size = 10, #set size of text
             base_family = "",
             base_line_size = base_size/20,
             base_rect_size = base_size/20)+
  geom_text(aes(x=1.3,size=20,label = paste(format(round(as.numeric(M.group)*100,2),nsmall=2),"%",sep="")),color = c("black"),position = position_stack(vjust = 0.5))+ #add text number to pie plot，更换物种记得修改M.group
  theme(plot.margin = margin(t=0, r=0, b=1, l=0, "cm"),
        axis.title = element_text(color="black",size=15),
        axis.text = element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
        legend.title = element_text(color = "black", size = 15), #设置图例标题大小为15，黑色
        legend.text = element_text(color = "black", size = 15),
        legend.position = "none")
ggsave(plot=pie,"Pie_M.group_SV_101701543.pdf" ,width=5,height=4)
