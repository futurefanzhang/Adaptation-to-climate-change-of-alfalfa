mydata=read.table("Dtrios_sativa_falcata_296_notree_BBAA.txt",header = T)
##type1,P1:west asia,P3:falcata
mydata2=mydata[mydata$P3=="falcata",]
mydata3=mydata2[mydata2$P1=="West_Asia",]
p<-ggplot(data=mydata3, aes(x=P2, y=f4.ratio,fill=P2)) +
  geom_bar(stat="identity",width = 0.6)+ylim(0,0.15)+
  labs(x = "P2", title = "",
       y = "F4-ratio") +
  scale_fill_manual("",values = c("#19B700","#42d4f4","#BF3EFF"),
                     breaks=c("Central_Asia", "East_Asia", "Europe_north_america"),
                     labels=c("Central Asia", "East Asia", "Europe and North America"))+
  theme_classic() + theme( 
    legend.position = "none",
    axis.title.x = element_text(color = "black",size = 10),
    axis.title.y = element_text(color = "black",size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color = "black",size = 10, angle = 30,vjust = 0.5),
    axis.text.y = element_text(color = "black",size = 10)
  )
ggsave(plot=p,"p1_west_asia_p3_falcata.pdf", width=3,height=5) 
