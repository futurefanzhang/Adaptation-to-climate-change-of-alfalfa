setwd("D:\\Population_643\\all_761_info")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtreeExtra")
##R4.2.3
library(ggtreeExtra)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggstar)
library(ggnewscale)
library(dplyr)

tr <- read.tree("merge_712_all_medicago.nwk.treefile") ##the output of iqtree results
data1=read.table("iqtree_group712.txt",header = T,sep = "\t") ##at lest two columns, sample name and group information

dat4 <- data1 %>% select(c("number1", "group"))
colnames(dat4)=c("ID","CLADE")
dat4 <- aggregate(.~CLADE, dat4, FUN=paste, collapse=",")
clades <- lapply(dat4$ID, function(x){unlist(strsplit(x,split=","))})
names(clades) <- dat4$CLADE
tr <- groupOTU(tr, clades, "Clade")
Clade <- NULL
p <- ggtree(tr, aes(color=Clade),branch.length = "none",layout="circular", size=0.2)+ geom_tiplab(size=0.5)+
  theme(legend.key.width = unit(1, 'cm'),legend.position = c(0.2,0.7),legend.text = element_text(size=10),legend.title = element_text(size=10))+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_nodelab(size=1) +
  scale_color_manual(
    values = c("#FABB2E","#19B700","#ee0000","#BF3EFF","#42d4f4","#A3A500","#8794FF","#A8422D","#00868B","gray","gray"),
    breaks=c("M.sativa", "M.caerulea","M.varia","M.falcata.tetraploid","M.falcata.diploid", "M.ruthenica","M.archiducis-nicolai","M.truncatula","M.group","Outgroup","0"),
                      labels=c("M.sativa", "M.caerulea","M.varia","M.falcata.tetraploid","M.falcata.diploid", "M.ruthenica","M.archiducis-nicolai","M.truncatula","M.group","Outgroup","Outgroup"))+
  labs(colour = "Species")
ggsave(plot=p,"tree_clade_iqtree_712_addbootstrap.pdf",  width=15,height=4) 

##add outside circle
cir=as.data.frame(data1[,c(5)])##use life history info
order=match(tr$tip.label,data1$number1)
cir=as.data.frame(cir[order,])
rownames(cir)=tr$tip.label
colnames(cir)="habit"
p1 <- gheatmap(p, cir, offset=-0.1, width=.05,colnames = FALSE,low = "gray", high = "red",
               colnames_angle=95, colnames_offset_y = .05) +
  scale_fill_manual("Growth Habit",values = c("#619CFF","#F8766D","gray"),
                    breaks=c("Perennial", "Annual","Outgroup"),
                    labels=c("Perennial", "Annual","Outgroup"))
ggsave(plot=p1,"tree_clade_iqtree_circular712_add_circle.pdf",  width=15,height=4) 
