##ggplot manhattan refer:https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##step1 use vcftools calculate Fst value 
##vcftools --vcf all712_indel_filter_missing0.2_meanDP3-62_allele2.vcf.recode.vcf --weir-fst-pop m.sativa_399_accession.txt --weir-fst-pop m.falcata.tetra.accession.txt --out m.sativa_vs_m.falcata.tetra_indel_10kb_Fst --fst-window-size 10000
sativa_falcata<-read.table("m.sativa_vs_falcata_tetra_10kb_Fst.windowed.weir.fst",header=TRUE,sep="\t")
all_pop_use=sativa_falcata[sativa_falcata$N_VARIANTS>=10 ,] #only use SNP number more than 10
all_pop_use$MEAN_FST[all_pop_use$MEAN_FST<0] = 0 #change less than zero as zero
write.table(all_pop_use,"use_fst_sativa_vs_falcata_tetra_10kb.txt",quote = F,row.names = F)
mean(all_pop_use$MEAN_FST)
thres0.01<-quantile(all_pop_use$MEAN_FST,0.99,na.rm=T) #top 1% threshold

###calculate Fst with point figure
mydata <-read.table("use_fst_sativa_vs_falcata_tetra_10kb.txt",header=TRUE,sep="\t")
mydata$ID<-paste(mydata$CHROM,mydata$BIN_START, sep="_")
mydata1<-mydata[,c(7,1,2,6)]
colnames(mydata1) <- c("ID","CHR", "BP","Fst")
add_gene=read.table("candidate_gene_region_new_A17.txt",header = F,sep = "\t")  ##please refer candidate_gene.R
add_gene$SNP<-paste(add_gene$V1,add_gene$V2, sep="_")
##ggplot
library(dplyr)
library(ggplot2)
library(ggrepel)
gwas_data=mydata1[order(mydata1$CHR,mydata1$BP),]
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
xmin=min(gwas_data$bp_cum)
xmax=max(gwas_data$bp_cum)

highlight_point=gwas_data[(gwas_data$ID %in% add_gene$SNP),]
order=match(highlight_point$ID, add_gene$SNP)
add_gene=add_gene[order,] #change order
highlight_point$gene=add_gene$V4
thres0.01<-quantile(mydata1$Fst,0.99,na.rm=T)
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = Fst, 
                                  color = as.factor(CHR)),size=1.5) +
  geom_point() +
  geom_hline(yintercept =thres0.01 , color = 'red', linetype = "solid") + 
  geom_point(data=highlight_point, 
             aes(x = bp_cum, y = Fst, color = as.factor(CHR)),color='black')+
  geom_text_repel(data=highlight_point,
    aes(label = highlight_point$gene),#box.padding = unit(0.35, "lines"),
    size = 3,
    nudge_y      = 0.16,
    direction    = "x",
    angle        = 0,
    hjust        = 0.15,
    segment.size = 0.3,
    max.iter = 1e4, max.time = 1,
    arrow = arrow(length = unit(0.02, "npc"),type = "closed"),color="black",
  )+
  scale_x_continuous(label = axis_set$CHR,limits = c(xmin,xmax),expand = expansion(0), breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.5),breaks = seq(0,0.5,0.1)) +
  scale_color_manual(values = rep(c("#FF939E","#FCD966"), unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome", title = "M.sativa vs M.falcata(tetraploid)",
       y = "Fst") + 
  theme_classic() +
  theme( 
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 15, vjust = 0.5),
    axis.text.y = element_text(size = 15)
  )

ggsave(plot=manhplot,"Fst between sativa and falcata(tetraploid).png" ,width=15,height=4,units='in',dpi=500)
ggsave(plot=manhplot,"Fst between sativa and falcata(tetraploid)_color.pdf" ,width=15,height=4)
