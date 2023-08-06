ABBA<-read.table("m_type0_type1_falcata_ABBA.csv",header=TRUE,sep=",")
ABBA=na.omit(ABBA)
ABBA$fd[ABBA$fd<0] = NA #change less than zero as zero
ABBA$fd[ABBA$fd>1] = NA #change less than zero as zero
ABBA=na.omit(ABBA)
mean(ABBA$fd)
write.csv(ABBA,"m_type0_type1_falcata_ABBA_use.csv",row.names = F)
##density plot
p1=ggplot(ABBA, aes(x = fd)) + geom_histogram(aes(y = ..density..),colour = 1, fill = "white") +
  geom_density(lwd = 1, colour = 4,fill = 4, alpha = 0.25)+xlab("fd")+ylab("Density")+
  theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
        axis.text = element_text(color="black",size=10),
        axis.title.y = element_text(color="black",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(color="black",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
        legend.title = element_text(color = "black", size = 10), #设置图例标题大小为20，黑色
        legend.text = element_text(color = "black", size = 10),
        legend.position = "none")
ggsave(plot=p1,"m_africa21_south_america_25_europe_north_america_69_falcata_ABBA_fd.pdf" ,width = 210, height = 70,units = "mm")
##manhattan
mydata <- ABBA
mydata$ID<-paste(mydata$scaffold,mydata$start, sep="_")
mydata1<-mydata[,c(12,1,2,10)]
colnames(mydata1) <- c("ID","CHR", "BP","Fd")
thread=0.3 #set threshold as 0.3
library(dplyr)
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
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = Fd,
                                  color = as.factor(CHR)),size=1.5) +
  geom_point() +
  geom_hline(yintercept =thread , color = 'red', linetype = "solid") +
  scale_x_continuous(label = axis_set$CHR,limits = c(xmin,xmax),expand = expansion(0), breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 1),breaks = seq(0,1,0.2)) +
  scale_color_manual(values = rep(c("grey60","#4197d8"), unique(length(axis_set$CHR)))) +
  labs(x = "Chromosome", title = "",
       y = "Fd") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 15, vjust = 0.5),
    axis.text.y = element_text(size = 15)
  )
ggsave(plot=manhplot,"Fd africa south america and ENA from M.falcata.png" ,width=15,height=4,units='in',dpi=500)
