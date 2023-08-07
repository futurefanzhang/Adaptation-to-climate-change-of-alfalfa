setwd("D:\\Population_643\\all_761_info\\local_adaption\\future_climate")
require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
library(sf)
library(rlist)
##############
#Read in population locations & climate data
##############
pops <- read.table("116_group_lat_long.txt",header = T,sep = "\t")
pops=pops[,c(1,2,3)]
colnames(pops)=c('code','long','lat')
#choose predictors
#predNames <- c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')
predNames <- c('wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18')
present <-list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_2.5m_bio\\", pattern = "*.tif$", full.names=TRUE)
presClim <- stack(present)
presClim <- presClim[[predNames]]
#Creates pred data for gf (cols = population name, long, lat, climate data)
pred <- data.frame(pop=pops$code,long=pops$long, lat=pops$lat, raster::extract(presClim, pops[,c("long","lat")]),stringsAsFactors=FALSE)

######################
#GF model
######################
#read in maf data
genotype=read.table("116_position_MAF_freq_R_format.txt",header = T,sep = "\t")
snpsAll=genotype[,-1]
#get snp names
snps <- cbind(pops[,c(1,2,3)],snpsAll)
#merge climate data and maf
snps <- merge(pred, snps, by.x="pop", by.y="code", all.x=TRUE)
#run gradient forest model,save it, next time can directly use it. 
gfMod <- gradientForest(data=snps, predictor.vars=predNames,response.vars=colnames(snpsAll),
                        ntree=500, 
                        maxLevel=log2(0.368*nrow(snps)/2), trace=T, 
                        corr.threshold=0.50)

list.save(gfMod, 'gfMod_forward.rdata')
gfMod=list.load("gfMod_forward.rdata")

#load future climate data
futClims <- stack("wc2.1_10m_bioc_ACCESS-CM2_ssp245_2081-2100.tif") #stack future climate layers
predNames1=c('bio19','bio08','bio09','bio15','bio02',"bio17","bio18")
futClims <- futClims[[predNames1]]
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)
colnames(futClimDat)=c("x","y",'wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18')
#transform future climate data with gradient forest model
futClimDatGF <- data.frame(futClimDat[,c("x","y")], predict(gfMod,futClimDat[,predNames])) 
saveRDS(futClimDatGF,file="futClimDatGF_forward_future_ssp245_2081_2100.rds")

#transform current climate data with gradient forest model
present <-list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_10m_bio\\", pattern = "*.tif$", full.names=TRUE)
presClim <- stack(present)
##extract 7 present bio data by the location of future tif file
predNames2=c("long","lat",'wc2.1_10m_bio_19','wc2.1_10m_bio_8','wc2.1_10m_bio_9','wc2.1_10m_bio_15','wc2.1_10m_bio_2','wc2.1_10m_bio_17','wc2.1_10m_bio_18')
pre_clim=present_all_bio19[,predNames2]
colnames(pre_clim)=c("x","y",'wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18')
popDatGF <- data.frame(pre_clim[,c("x","y")], predict(gfMod, pre_clim[,predNames]))
saveRDS(popDatGF,file="popDatGF_forward_future_present.rds")

#Forward offset calculation
##############
#########This step I use server rather than my personal laptop, more CPU can be used
#setwd('/public/agis/zhouyongfeng_group/zhangfan02/genetic_offset')
require(fields)
require(geosphere)
require(parallel)
#require(gdm)
require(doParallel)
popDatGF=readRDS("popDatGF_forward_future_present.rds")
futClimDatGF=readRDS("futClimDatGF_forward_future_ssp245_2081_2100.rds")
predNames <- c('wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18')

cl <- makeCluster(40)
registerDoParallel(cl)
forwardOffsetGF <- foreach(i = 1:nrow(popDatGF), .packages=c("fields","geosphere")) %dopar%{
  #get the focal population
  onePopGF <- popDatGF[i,]
  #get destination populations and add gf distance
  combinedDatGF <- futClimDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], futClimDatGF[,predNames]))
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  #choose the pixels with the minimum gfOffse
  #############
  combinedDatGF['dists']=distGeo(p1=coordGF, p2=combinedDatGF[,1:2])
  #limit the migrant distance 
  combinedDatGF<-combinedDatGF[combinedDatGF['dists']<1000000,] ##set distance修改迁移距离
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  #calculate the distance to the sites with minimum gfOffse, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dist == min(minCoordsGF$dists)),]
  #if multiple sites have the same gfOffset, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  #get local offset
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y ==coordGF$y),"gfOffset"]
  #get the minimum predicted fst - forward offset in this case
  minValGF <- minCoordsGF$gfOffset
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF, forwardOffset=minValGF, predDist=toGoGF, bearing=bearGF,x2=minPtGF[[1]],y2=minPtGF[[2]])
}
stopCluster(cl)
forwardOffsetGF <- do.call(rbind, forwardOffsetGF)
saveRDS(forwardOffsetGF,file="future_SSP245_2081_1000km_ForwardOffsetGF_7bio.rds")


##plot
library(dplyr)
library(rasterVis)
library(RColorBrewer)

mydata397=read.table("present_bio_info_19BIO.txt",header = T,sep="\t") ##use data point
use_location <- sp::SpatialPoints(mydata397[,c("Longitude","Latitude")])
Offset=readRDS("future_SSP245_2081_unlimit_ForwardOffsetGF_7bio.rds")
mask=raster::raster("wc2.1_10m_bioc_ACCESS-CM2_ssp585_2081-2100.tif") #I use present tif file generate raster for plot, it is corresponding with offset x and y columns
mask1=raster::as.data.frame(mask, xy=TRUE)
Offset=as.data.frame(Offset[,c("x1","y1","forwardOffset")])
colnames(Offset)=c("x","y","offset")
Offset=Offset[Offset$x<=180 & Offset$x>=(-170),]
Offset=Offset[Offset$y<=90 & Offset$y>=(-60),]
res<-quantile(Offset$offset, probs = c(0,0.25,0.5,0.75,0.98,1))
res ##set level strandand
Offset$offset[Offset$offset>0.004]=0.004 ##set max Genetic offset to 0.006, based on 98% data info
mask2 <- mask1 %>% left_join(Offset, by=c('x'='x', 'y'='y')) ##add offset info to world map

mask$bio01[]=mask2$offset ##change bio1 info to offset
names(mask)="offset"

pdf("Genetic offset ssp245_2081_2100_unlimit_bio7.pdf",width = 8,height = 4)
png("Genetic offset ssp245_2081_2100_unlimit_bio7.png",width = 8,height = 4, units="in",res = 300)
rasterVis::levelplot(mask, main = "Genetic offset", margin = FALSE,at=seq(0,0.004, length.out=100),maxpixels = 2e6,colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
                     par.settings=rasterVis::rasterTheme(c("#084594","#2171b5","#4292c6","#6baed6","#9ecae1","#c6dbef","#dee6fc","#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"))) +
                      latticeExtra::layer(sp.points(use_location, col="black",pch=2,cex=0.5))
dev.off()
