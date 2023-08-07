##continue for evaluate_representative_biofactor.R
##this part used for offset evaluate by Euclidean distance
most_important=c("BIO19","BIO8","BIO9","BIO15","BIO2","BIO17","BIO18") #use 7 most important bioclim
all_tgrid=cbind(present[,c("X","Y")], predict(all_gfmod,present[,most_important])) #"Longitude","Latitude"
all_PCs <- prcomp(all_tgrid[,most_important])
summary(all_PCs)

# set up a colour palette for the mapping, add color for different geographical regions
a1 <- all_PCs$x[,1]
a2 <- all_PCs$x[,2]
a3 <- all_PCs$x[,3]
r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255
grid=present[,c("X","Y")]
grid$R=r
grid$G=g
grid$B=b
grid$cols=rgb(r,g,b,max=255)
write.table(grid,"BIO7_present_RGB_all_location.txt",row.names = F,quote = F,sep = "\t")
nvs <- dim(all_PCs$rotation)[1] # number of variables
vec <- most_important #only use the top five most important
lv <- length(vec)
vind <- rownames(all_PCs$rotation) %in% vec
scal <- 80
xrng <- range(all_PCs$x[,1], all_PCs$rotation[,1]/scal)*1.1
yrng <- range(all_PCs$x[,2], all_PCs$rotation[,2]/scal)*1.1
pdf(file="picture/all_PCplot_bio7.pdf",width = 7.5,height = 6)
plot((all_PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".",xlab = "PC1 (47.6%)",ylab = "PC2 (32.9%)", cex=7, col=rgb(r,g,b, max = 255), asp=1)
arrows(rep(0,lv), rep(0,lv), all_PCs$rotation[,1]/scal, all_PCs$rotation[,2]/scal, length = 0.04)
jit <- 0.0007
text(all_PCs$rotation[,1]/scal+jit*sign(all_PCs$rotation[,1]), all_PCs$rotation[,2]/scal+jit*sign(all_PCs$rotation[,2]), labels = vec)
dev.off()

pdf("picture/Map_bio7.pdf",width = 7.5,height = 5)
plot(all_tgrid[, c("X","Y")],pch=".", cex=3, xlab = "Longitude",ylab = "Latitude",asp=1, col=rgb(r,g,b, max = 255))
dev.off()

###evaluate offset 
all_tgrid=cbind(present[,c("X","Y")], predict(all_gfmod,present[,most_important])) #X represent "Longitude",Y represent "Latitude"
fut_cli=raster::stack("wc2.1_10m_bioc_ACCESS-CM2_ssp245_2061-2080.tif") ##your can change future model from here
fut_cli=data.frame(X=present$X, Y=present$Y, raster::extract(fut_cli, present[,c("X","Y")]),stringsAsFactors=FALSE)
most_important1=c("X","Y","bio19",'bio08','bio09','bio15','bio02','bio17','bio18')
fut_cli=fut_cli[,most_important1]
colnames(fut_cli)=c("X","Y","BIO19","BIO8","BIO9","BIO15","BIO2","BIO17","BIO18")
fut_cli=na.omit(fut_cli)
future_all=cbind(fut_cli[,c("X","Y")], predict(all_gfmod,fut_cli[,most_important]))
library(dplyr)
df2 <- future_all %>% left_join( all_tgrid,
                                 by=c('X'='X', 'Y'='Y')) #combine two dataframe

genOffsetAll<-sqrt((df2[,3]-df2[,3+7])^2+(df2[,4]-df2[,4+7])^2+(df2[,5]-df2[,5+7])^2+(df2[,6]-df2[,6+7])^2+(df2[,7]-df2[,7+7])^2+(df2[,8]-df2[,8+7])^2+(df2[,9]-df2[,9+7])^2) ##欧式距离Euclidean distance
Offset=cbind(df2[,c("X","Y")],genOffsetAll)
colnames(Offset)[3]<-"offset"
saveRDS(Offset,file="bio7_ssp245_2061_2080_Offset.rds") ##save as R data

##R figure
library(dplyr)
library(rasterVis)
library(RColorBrewer)
mydata397=read.table("present_bio_info_19BIO.txt",header = T,sep="\t") ##use data point
use_location <- sp::SpatialPoints(mydata397[,c("Longitude","Latitude")])
Offset=readRDS("bio7_ssp245_2081_2100_Offset.rds")
mask=raster::raster("D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_10m_bio\\wc2.1_10m_bio_1.tif") #I use present tif file generate raster for plot, it is corresponding with offset x and y columns
mask1=raster::as.data.frame(mask, xy=TRUE)
colnames(Offset)=c("x","y","offset")
res<-quantile(Offset$offset, probs = c(0,0.25,0.5,0.75,0.98,1))
res ##set level strandand
Offset$offset[Offset$offset>0.004]=0.004 ##set max Genetic offset to 0.004, based on 98% data info
mask2 <- mask1 %>% left_join(Offset, by=c('x'='x', 'y'='y')) ##add offset info to world map
mask$wc2.1_10m_bio_1[]=mask2$offset ##change bio1 info to offset
names(mask)="offset"

#pdf("Genetic offset ssp585_2081_2100_bio7_new.pdf",width = 8,height = 4)
png("Genetic offset ssp245_2081_2100_bio7_use.png",width = 8,height = 4,units = "in",res = 500)
rasterVis::levelplot(mask, main = "SSP245 (2081-2100)", margin = FALSE,at=seq(0,0.004, length.out=100),maxpixels = 2e6,colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                     par.settings=rasterVis::rasterTheme(c("#084594","#2171b5","#4292c6","#6baed6","#9ecae1","#c6dbef","#dee6fc","#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"))) +
  latticeExtra::layer(sp.points(use_location, col="black",pch=2,cex=0.5)) ##add sampling point            
dev.off()
