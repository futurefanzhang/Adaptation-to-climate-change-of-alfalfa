###R4.2.3
setwd("D:\\Population_643\\all_761_info\\local_adaption\\future_climate") # set working directory
require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
library(sf)
require(rlist)
library(dplyr)
pops <- read.table("116_group_lat_long.txt",header = T,sep = "\t")
pops=pops[,c(1,2,3)]
colnames(pops)=c('code','long','lat')
#choose predictors
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
#run gradient forest model
gfMod <- gradientForest(data=snps, predictor.vars=predNames,response.vars=colnames(snpsAll),
                        ntree=500, 
                        maxLevel=log2(0.368*nrow(snps)/2), trace=T, 
                        corr.threshold=0.50)
library(rlist)
list.save(gfMod, 'gfMod_forward_reverse_offset.rdata')
gfMod=list.load("gfMod_forward_reverse_offset.rdata")

#load future climate data
futClims <- stack("wc2.1_10m_bioc_ACCESS-CM2_ssp585_2081-2100.tif") #stack future climate layers
predNames1=c('bio19','bio08','bio09','bio15','bio02',"bio17","bio18")
futClims <- futClims[[predNames1]]
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)
colnames(futClimDat)=c("x","y",'wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18')
#transform future climate data with gradient forest model
futClimDatGF <- data.frame(futClimDat[,c("x","y")], predict(gfMod,futClimDat[,predNames])) 
saveRDS(futClimDatGF,file="futClimDatGF_forward_reverse_ssp585_2081_2100.rds")

#transform current climate data with gradient forest model
present <-list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_10m_bio\\", pattern = "*.tif$", full.names=TRUE)
presClim <- stack(present)
##extract all present bio19 data by the location of future tif file
present_all_bio19 <- data.frame(long=futClimDatGF$x, lat=futClimDatGF$y,raster::extract(presClim, futClimDatGF[,c("x","y")]),stringsAsFactors=FALSE)

predNames2=c("long","lat",'wc2.1_10m_bio_19','wc2.1_10m_bio_8','wc2.1_10m_bio_9','wc2.1_10m_bio_15','wc2.1_10m_bio_2','wc2.1_10m_bio_17','wc2.1_10m_bio_18')
pre_clim=present_all_bio19[,predNames2]
colnames(pre_clim)=c("x","y",'wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18')
popDatGF <- data.frame(pre_clim[,c("x","y")], predict(gfMod, pre_clim[,predNames]))
saveRDS(popDatGF,file="popDatGF_forward_reverse_present.rds")


###I calculate forward and reverse offset in server
####read previous saved results
setwd('/public/agis/zhouyongfeng_group/zhangfan02/genetic_offset')
require(fields)
require(geosphere)
require(parallel)
require(doParallel)
popDatGF=readRDS("popDatGF_forward_reverse_present.rds") ##导入第一步的当前全球数据结果
futClimDatGF=readRDS("futClimDatGF_forward_reverse_ssp245_2081_2100.rds") ##导入第一步的未来全球数据点信息
predNames <- c('wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18')
cl <- makeCluster(40) ###multiple cores is quicker
registerDoParallel(cl)

###############
#Forward offset calculation
##############
forwardOffsetGF <- foreach(i = 1:nrow(popDatGF), .packages=c("fields","geosphere")) %dopar%{
  onePopGF <- popDatGF[i,]
  #get destination populations and add gf distance
  combinedDatGF <- futClimDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], futClimDatGF[,predNames]))
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  #choose the pixels with the minimum gfOffse
  #############
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
#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#forwardOffset: forward offset
#predDist: distance to site of forward offset
#bearing: bearing to site of forward offset
#x2/y2: coordinate of site of forward offset
saveRDS(forwardOffsetGF,file="future_SSP245_2081_2100_ForwardOffsetGF_7bio.rds")  ##保存RDS格式

###############
#Reverse offset calculation
##############
cl <- makeCluster(40) ###multiple cores is quicker
registerDoParallel(cl)
reverseOffsetGF <- foreach(i = 1:nrow(futClimDatGF), .packages=c("fields","geosphere")) %dopar%{
  #get the focal population in future climate
  onePopGF <- futClimDatGF[i,]
  #make prediction between focal population and current climate
  combinedDatGF <- popDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], popDatGF[,predNames]))
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  #choose the pixels with the minimum offset
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  #get local offset
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y == coordGF$y),"gfOffset"]
  #get the minimum predicted offset - reverse offset in this case
  minValGF <- minCoordsGF$gfOffset
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]],local=offsetGF, reverseOffset=minValGF, predDist=toGoGF, bearing=bearGF, x2=minPtGF[[1]],y2=minPtGF[[2]])
}
stopCluster(cl)
#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#reverseOffset: reverse offset
#predDist: distance to site of reverse offset
#bearing: bearing to site of reverse offset
#x2/y2: coordinate of site of reverse offset
reverseOffsetGF <- do.call(rbind, reverseOffsetGF)
saveRDS(reverseOffsetGF,file="future_SSP245_2081_2100_reverseOffsetGF_7bio.rds")  ##保存RDS格式


############figure##############
#####################input the forward and reverse results from server running results, generate figure
forwardOffsetGF=readRDS("future_SSP245_2081_2100_ForwardOffsetGF_7bio.rds")
forwardOffsetGF=as.data.frame(forwardOffsetGF)
reverseOffsetGF=readRDS("future_SSP245_2081_2100_reverseOffsetGF_7bio.rds")
reverseOffsetGF=as.data.frame(reverseOffsetGF)
#data1 <- forwardOffsetGF %>% inner_join(reverseOffsetGF, by=c('x1'='x1', 'y1'='y1'))
#data1$mean_local=(data1$local.x+data1$local.y)/2

##color
library(sf)
library(ggplot2)
library(dplyr)

creategroup <- function(tiff){
  colnames(tiff)=c("x","y","bio_value")
  tiff$level=
    ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.95),'255',
           ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.9),'220',
                  ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.85),'190',
                         ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.8),'180',
                                ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.75),'140',
                                       ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.7),'135',
                                              ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.65),'130',
                                                     ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.6),'125',
                                                            ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.55),'120',
                                                                   ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.5),'115',
                                                                          ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.45),'110',
                                                                                 ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.4),'105',
                                                                                        ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.35),'100',
                                                                                               ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.3),'95',
                                                                                                      ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.25),'90',
                                                                                                             ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.2),'85',
                                                                                                                    ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.15),'40',
                                                                                                                           ifelse(tiff$bio_value>=quantile(tiff$bio_value, probs = 0.1),'20','0'
                                                                                                                           ))))))))))))))))))
  return(tiff)
}


use_forward=forwardOffsetGF[,c("x1",'y1','local','forwardOffset')]
use_reverse=reverseOffsetGF[,c("x1",'y1','reverseOffset')]
data1 <- use_forward %>% inner_join(use_reverse, by=c('x1'='x1', 'y1'='y1'))
###get r,g,b value from local, forward, reverse offset

data=data1
data=data[data$x>=-170 & data$x<=180,]
data=data[data$y>=-60 & data$y<=90,]

####add RGB color
local=creategroup(data[c("x1",'y1','local')])
colnames(local)=c('x','y','local','local_red')
forward=creategroup(data[c("x1",'y1','forwardOffset')])
colnames(forward)=c('x','y','forward','forward_green')
reverse=creategroup(data[c("x1",'y1','reverseOffset')])
colnames(reverse)=c('x','y','reverse','reverse_blue')
temp=local %>% inner_join(forward, by=c('x'='x', 'y'='y'))
data=temp %>% inner_join(reverse, by=c('x'='x', 'y'='y'))

data$color=NA
for (i in 1:nrow(data)){
  color=rgb(as.numeric(data$local_red[i]),as.numeric(data$forward_green[i]),as.numeric(data$reverse_blue[i]),maxColorValue = 255 ) 
  data$color[i]=color
}
saveRDS(data,file="future_SSP245_2081_2100_color_bio7.rds")

##point plot, correlation among local, forward and reverse
data=readRDS(file="future_SSP585_2081_2100_color_bio7.rds")
library(ggplot2)
color=levels(factor(data$color))
p1<-ggplot(data=data)+
  geom_point(aes(x=local,y=forward,color=as.factor(color)),size=1,shape=19)+
  labs(x="Local offset",y="Forward offset")+
  scale_color_manual(
    values=color
    )+
  geom_abline(slope=1, intercept=0)+
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),
                        axis.text = element_text(color="black",size=15),#axis.line = element_line(colour = "black", size = 1),
                        #axis.title.y = element_text(color="black",size=15,margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        #axis.title.x = element_text(color="black",size=15,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
                        legend.title = element_text(color = "black", size = 15), 
                        legend.text = element_text(color = "black", size = 15),
                        legend.position = "none",
        panel.background = element_rect(fill = "gray43",
                                        colour = "gray43",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray61"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "gray61"))

ggsave(plot=p1,"local_forward_genetic_offset_SSP245.png", width=4,height=4,dpi = 300)

##world map
mydata397=read.table("present_bio_info_19BIO.txt",header = T,sep="\t")
library(ggplot2)
library(ggmap)
library(sp)
library(RColorBrewer)
library(maptools)
library(maps) 
mp<-NULL 
mapworld<-borders("world",colour = "white",fill="white")
mp<-ggplot()+mapworld+ylim(-60,90)+xlim(-170,180)
color=levels(factor(data$color))
p=mp+geom_point(aes(x=x, y=y,color=as.factor(color)),shape=19,size=0.02,data=data)+
  scale_color_manual(values=color)+
  theme(axis.text = element_blank(), 
        legend.position = "none",legend.key.size = unit(1.5,"line"),
        axis.title=element_blank(),axis.ticks = element_blank(),panel.background = element_blank(),
        legend.title = element_text(color = "black", size = 15), 
        legend.text = element_text(color = "black", size = 15))+
  geom_point(aes(x=Longitude, y=Latitude),color="white",shape=17,size=2,data=mydata397)
ggsave(plot=p,"world_map_geography_local_forward_reverse_genetic_offset_SSP585_addpoint.png", width = 15, height = 6)


###partical figure color, I use these command for generate origin center offset info
data=readRDS("future_SSP585_2081_2100_color_bio7.rds")
data=data[data$x>=30 & data$x<=60,]
data=data[data$y>=32 & data$y <=50,]
color=levels(factor(data$color))
country_shp=sf::st_read("World_Countries/World_Countries.shp") ##down shp file from http://tapiquen-sig.jimdo.com.
p1=ggplot()+
  geom_raster(aes(x=x,y=y,fill=as.factor(color)),data=data,stat="identity")+
  geom_sf(fill="transparent",data=country_shp,size=0.5)+
  coord_sf(xlim=c(30, 60),ylim=c(32, 50))+
  scale_fill_manual(values=color)+
  scale_x_continuous(limits = c(30, 60))+
  scale_y_continuous(limits = c(32, 50))+
  labs(x="Longitude",y="Latitude") +
  theme_bw()+
  geom_point(aes(x=Longitude, y=Latitude),color="white",shape=17,size=4,data=mydata397)+
  theme(text=element_text(family="serif"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key.size=unit(0.1,'inch'),
        legend.title = element_text(size=10.5, color="black",vjust=0.5, hjust=0.5),
        legend.position = "none",
        legend.background=element_rect(colour= "grey" ,fill= "white" ,size=0.6),
        legend.text= element_text(size=7.3, color="black",vjust=0.5, hjust=0.5),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(colour = "white",size=0.1,linetype = 4),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  guides(fill=F)
ggsave(p1,file="alfalfa_forward_reverse_local_SSP585_RGB_region.pdf",height=8,width=8)
