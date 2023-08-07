##step 1 run PBS analysis in linux
plink --vcf merge_409_iqtree_plink_addhead.vcf --keep east_west_asia_10truncala_95sample.txt --allow-extra-chr --recode vcf --out east_west_asia_outgroup_pbs_95_sample --threads 1
java -Xmx10g -jar beagle.22Jul22.46e.jar nthreads=1 gt=east_west_asia_outgroup_pbs_95_sample.vcf out=east_west_asia_outgroup_pbs_95_sample_beagle.vcf
gzip -d east_west_asia_outgroup_pbs_95_sample_beagle.vcf.vcf.gz
pbscan -vcf east_west_asia_outgroup_pbs_95_sample_beagle.vcf.vcf -pop1 east_asia48sample.txt -pop2 west_asia37sample.txt -pop3 outgroup10sample.txt -out east_west_raw -win 100 -step 101 -div 2 -min 2 -mc 100 #-win窗口大小100个SNP，-step为SNP的步长，-div代表分化测定方式，2代表Dxy，-min 2 代表每个群体最少2个个体，可选参数：增加-mc参数可以导出带P值的分析结果，代表测试分析循环数

##step 2, draw figure in R
mydata=read.table("europe_north_america.pbs",header = T,sep = "")

mydata$ID<-paste(mydata$Chromo,mydata$Middle, sep="_")
mydata$logp=-log10(mydata$P1)

mydata1<-mydata[,c(12,1,3,13)]
colnames(mydata1) <- c("ID","CHR", "BP","logp")
mydata1$logp[!is.finite(mydata1$logp)]<-4.2 ##change zero P value to -log10(P)=4

##ggplot manhattan

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

thres0.05<-(-log10(0.01))
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = logp, 
                                  color = as.factor(CHR)),size=1.5) +
  geom_point() +
  geom_hline(yintercept =thres0.05 , color = 'red', linetype = "dashed") + 
  scale_x_continuous(label = axis_set$CHR,limits = c(xmin,xmax),expand = expansion(0), breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 5),breaks = seq(0,5,1)) +
  scale_color_manual(values = rep(c("#fde8ad","#FCCC00"), unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome", title = "PBS",
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

ggsave(plot=manhplot,"europe and north america west asia P value.png" ,width=15,height=3,units='in',dpi=500)

ggsave(plot=manhplot,"europe and north america west asia P value.pdf" ,width=7.5,height=3)
