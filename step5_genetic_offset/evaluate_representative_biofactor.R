##ref:https://github.com/jingwanglab/Populus_genomic_prediction_climate_vulnerability/tree/main/7-Genetic_offset
library(gradientForest)
library(raster)
genotype=read.table("116_position_MAF_freq_R_format.txt",header = T,sep = "\t") ##maf matrix
phe=read.table("116_group_lat_long.txt",header = T,sep = "\t") #group info of 397 sample
order=match(genotype$SNP,phe$group1)
phe=phe[order,]
all_SNPs=genotype[,-c(1)]
present <-list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_2.5m_bio\\", pattern = "*.tif$", full.names=TRUE) ##download from here：https://www.worldclim.org/data/worldclim21.html
presClim <- stack(present)
pred <- data.frame(pop=phe$group1,long=phe$lon, lat=phe$lat, raster::extract(presClim, phe[,c("lon","lat")]),stringsAsFactors=FALSE)
envGF=pred[,-c(1:3)]
preds <- colnames(envGF)
specs <- colnames(all_SNPs)
nSites <- dim(envGF)[1]
nSpecs <- dim(all_SNPs)[2]
maxLevel <- log2(0.368*nrow(envGF)/2)
colnames(envGF)=c("BIO1","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9")
all_gfmod <- gradientForest(cbind(envGF, all_SNPs), predictor.vars=colnames(envGF), response.vars=colnames(all_SNPs), ntree=500, compact=T, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)
require(rlist)
list.save(all_gfmod, 'list.rdata') ##save model results
all_gfmod=list.load("list.rdata")
pdf(file="picture/all_predictoroverallimportance.pdf") ##biofactoer importance order
plot(all_gfmod,plot.type="O")
dev.off()
write.table(importance(all_gfmod),"result/weighted importance of bio19.txt",quote = F,sep = "\t") ##save importance results

##evaluate results
most_important <- names(importance(all_gfmod))[1:19]
#splits density plots
pdf(file="picture/all_splitsdensityplots.pdf")
plot(all_gfmod, plot.type="S",imp.vars = most_important, leg.posn="topright", cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))
dev.off()
#speciescumulativeplot #the legend identifies the top 5 most responsive SNPs for each predictor
pdf(file="picture/all_speciescumulativeplot.pdf")
plot(all_gfmod, plot.type="Cumulative.Importance", imp.vars=most_important, show.overall=T, legend=T,common.scale=T,leg.posn="topleft", leg.nspecies=5, cex.lab=0.7, cex.legend=0.4, cex.axis=0.6, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1),omi=c(0,0.3,0,0)))
dev.off()
#predictorcumulative #show cumulative change in overall composition of the community, where changes occur on the gradient
pdf(file="picture/all_predictorcumulative.pdf")
plot(all_gfmod, plot.type="C", imp.vars=most_important, show.species=F, common.scale=T, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(2.5,1.0,0.1,0.5), omi=c(0,0.3,0,0)))
dev.off()
#R2
pdf(file="picture/all_R2.pdf")
plot(all_gfmod, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)
dev.off()

###19 biofactor correlation
present <-list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_10m_bio\\", pattern = "*.tif$", full.names=TRUE)
###extrat first BIO
presClim=raster::stack(present[1])
presClim=as.data.frame(presClim,xy=TRUE)
##extract rest BIO
for (i in 2:19){
presClim1 <- raster::stack(present[i])
presClim[,i+2]=as.data.frame(presClim1,xy=TRUE)[,3]
}
names=list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_10m_bio\\", pattern = "*.tif$", full.names=F)
colnames(presClim)=c("X","Y",names)
presClim=na.omit(presClim)
present=presClim
present=present[present$X<=180 & present$X>=(-170),]
present=present[present$Y<=90 & present$Y>=(-60),]
##########test correlation between 19 bioclim
colnames(present)=c("X","Y","BIO1","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9")
present_cor=cor(present[,most_important])
#write.table(present_cor,"present_BIO19_all_location_cor.txt",quote = "\t",sep = "\t")
library(corrplot) #based on R 4.0.5
present_cor=as.matrix(present_cor)
pdf("picture/present_BIO19_all_location_cor.pdf",paper = "a4")
par(omi=c(0.1,0.1,0.1,0.1)) ##set cor 0.7 as threshold
corrplot(present_cor,method = 'color',type = 'lower', number.cex = 0.7,tl.srt = 25,addCoef.col = 'black', order = 'original',tl.col = 'black',col.lim = c(-1,1),col = colorRampPalette(c("blue","white","red"))(100),cl.cex=1)
dev.off()

##generate figure for representative biofactor
ata=read.table("picture/weighted importance of bio19_used_for_figure.txt",header = T,sep = "\t")
data$bio <- factor(data$bio, levels=data$bio)
p<-ggplot(data, aes(x=bio, y=importance, fill=type)) +
  geom_bar(stat="identity")+ylim(0,0.01)+
  labs(x = "", title = "",y = "R2 weighted importance") +
  scale_fill_manual("",values = c("#999999", "#E69F00"),
                    breaks=c("not_use","use"),
                    labels=c("not_use","use"))+
  theme_classic() + theme(
    legend.position = "none",
    axis.title.x = element_text(color = "black",size = 10),
    axis.title.y = element_text(color = "black",size = 10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color = "black",size = 10, angle = 30,vjust = 0.5),
    axis.text.y = element_text(color = "black",size = 10))
ggsave(plot=p,"19bio_importance.pdf", width=6,height=4)

##further calculate offset please refer offset_estimation.R
