
mydata=read.table("candidate_gene.sprot.blast",header = F,sep = "\t")
colnames(mydata)=c("qseqid", "sseqid",    "pident",    "length",    "mismatch",    "gapopen",    "qstart",    "qend",    "sstart",    "send",    "evalue",    "bitscore")
#filter p-value<1e-10
mydata=mydata[mydata$evalue<1e-10,]
mydata=mydata[order(mydata$pident,decreasing = TRUE),]
mydata=mydata[!duplicated(mydata$qseqid),] #remove duplicated gene info
write.table(mydata,"deg.sprot.blast_use.txt",sep = "\t",quote = F,row.names = F,col.names = F) #export use blast info
write.table(mydata$sseqid,"gene_function_and_Uniprot_info.txt",sep = "\t",quote = F,row.names = F) #export UniProt accessions info, use DAVID web check function,https://david.ncifcrf.gov/list.jsp

##Go enrichment plot
library(ggplot2)
setwd("D:\\Population_643\\all_761_info\\population_analysis\\sativa_falcata_tetraploid")

data1 <- read.table("Go_CC_DIRECT_fst_medicago_falcata_tetraploid.txt", header = T, sep = "\t") #Cellular Component,CC
#data1=data1[data1$FDR<0.05,] #根据自己情况更改
data1$GeneRatio=data1$Count/data1$List.Total
data1=data1[,c("Term","Count","GeneRatio","PValue")]
data1$Type="Cellular Component"

data2 <- read.table("Go_BP_DIRECT_medicago_falcata_tetraploid.txt", header = T, sep = "\t") #Biological Process,BP
#data2=data2[data2$FDR<0.05,] #根据自己情况更改
data2$GeneRatio=data2$Count/data2$List.Total
data2=data2[,c("Term","Count","GeneRatio","PValue")]
data2$Type="Biological Process"

data3 <- read.table("Go_MF_DIRECT_fst_medicago_falcata_tetraploid.txt", header = T, sep = "\t") #Molecular Function,MF
#data3=data3[data3$FDR<0.05,] #根据自己情况更改
data3$GeneRatio=data3$Count/data3$List.Total
data3=data3[,c("Term","Count","GeneRatio","PValue")]
data3$Type="Molecular Function"

data=rbind(data1,data2,data3)
df <- data[order(data$Type,data$GeneRatio),]
library(stringr)
df$Term=str_split_fixed(df$Term, '~', 2)[,2] ##仅使用~之后的信息
df$Term =  factor(df$Term, levels = unique(df$Term))

p1<-ggplot(df,aes(x = GeneRatio, y = Term)) +        # x 轴用GeneRatio, y轴用GO或KEGG注释
  geom_point(aes(color = -log10(PValue), size = Count), pch = 19) +            # 颜色用p.adjust，size设为用count
  scale_color_steps(low = "#023e7d", high = "#ff8000") +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  theme_bw() +
  theme(panel.grid=element_blank(), #去掉网格线
        panel.background= element_rect(size = 0.2, colour = 1),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  facet_wrap(~ data$Type , scales= "free", nrow=3)

ggsave(plot=p1,"GO_enrichment_fst_sativa_falcata.pdf" ,width=6.5,height=5)
