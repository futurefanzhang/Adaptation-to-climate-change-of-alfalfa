library(ggplot2)
setwd("D:\\Population_643\\all_761_info\\PCA")
mydata<-read.table("PCA_result1.txt",header=TRUE,sep="\t") 
mydata$group=factor(mydata$group,levels = c("M.sativa", "M.caerulea", "M.varia","M.falcata.tetraploid","M.falcata.diploid","M.ruthenica","M.archiducis-nicolai","M.truncatula","M.group"),ordered = TRUE) ##排序一下分组信息
p1<-ggplot(data=mydata)+
  geom_point(aes(x=PC1,y=PC2,colour=group),data=mydata,size=2,shape=1)+
  labs(x="PC1(10.20%)",y="PC2(7.21%)")+#use 702_accession_pca.eigenval calculate percent
  scale_color_manual(
    values = c("#FABB2E","#19B700","#ee0000","#42d4f4","#BF3EFF","#A3A500","#8794FF","#A8422D","#00868B"),
    breaks=c("M.sativa", "M.caerulea","M.varia","M.falcata.diploid", "M.falcata.tetraploid","M.ruthenica","M.archiducis-nicolai","M.truncatula","M.group"),
    labels=c("M.sativa", "M.caerulea","M.varia","M.falcata.diploid", "M.falcata.tetraploid","M.ruthenica","M.archiducis-nicolai","M.truncatula","M.group"))+
  labs(colour = "Species")+
  theme_classic()+theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
    axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
        axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
        legend.title = element_text(color = "black", size = 15), 
        legend.text = element_text(color = "black", size = 15),
        legend.position = "right")
ggsave(plot=p1,"PC_1_2.pdf", width=7.5,height=4)
