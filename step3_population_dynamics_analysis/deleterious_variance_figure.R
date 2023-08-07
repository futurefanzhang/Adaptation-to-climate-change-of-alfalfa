##this script calculate deleterious variance number by homozygosity or heterozygosity, or total variance number
mydata=read.table("R_use_format1.vcf",header = T,row.names = 1,sep = "\t")
group=read.table("ind_group.txt",header = T,sep = "\t") ##two columns, taxa and group info
het_del=colSums(mydata == "1",na.rm = T)
hom_del=colSums(mydata == "2",na.rm = T)
het_del=as.data.frame(het_del)
hom_del=as.data.frame(hom_del)
combine=cbind(het_del,hom_del)
order=match(rownames(combine),group$taxa)
combine2=cbind(combine,group[order,])
combine2$all_del=combine2$het_del+2*combine2$hom_del
aggregate(combine2$het_del, list(combine2$group), FUN=mean)
aggregate(combine2$hom_del, list(combine2$group), FUN=mean)
aggregate(combine2$all_del, list(combine2$group), FUN=mean)
library(ggpubr)
#t-test
stat.test <- compare_means(
  all_del ~ group, data = combine2,  method = "t.test"
)
stat.test
library(ggplot2)
##all info ,outgroup, introgression species
p1<-ggplot(combine2, aes(x=group, y=all_del, fill=group))+labs(x="Species", y="No. of deleterious \n alleles/accession") +
  geom_boxplot(width=0.8)+#ylim(0.7,1.2)+
  theme_classic()+theme(
    legend.position="none",legend.title = element_text(color = "black", size = 10), legend.text = element_text(color = "black", size = 10),axis.text = element_text(color="black",size=10),axis.title.x=element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size=10),axis.title.y=element_text(size=10))
ggsave(plot=p1,"M.sativa_deleterious.pdf", width = 310, height = 160,units = "mm")

##use different geographical regions
combine3=combine2[combine2$group!="falcata",]
combine3=combine3[combine3$group!="Outgroup",]
p1<-ggplot(combine3, aes(x=group, y=all_del, fill=group))+labs(x="Region", y="No. of deleterious \n alleles/accession") +
  geom_boxplot(width=0.6)+
  scale_fill_manual(
    values = c("#FABB2E","#19B700","#42d4f4","#BF3EFF","#8794FF","#ee0000"),
    breaks=c("Africa", "Central_Asia","East_Asia", "Europe_north_america","South_America","West_Asia"),
    labels=c("Africa", "Central Asia","East Asia", "Europe and North America","South America","West Asia"))+
  theme_classic()+theme(
    legend.position="none",legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    axis.text.y = element_text(color="black",size=10),
    axis.text.x = element_text(color = "black",size = 10, angle = 30,vjust = 0.5),
    axis.title.x=element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size=10),
    axis.title.y=element_text(size=10))
ggsave(plot=p1,"M.sativa_6_region_deleterious.pdf", width=6,height=4)
