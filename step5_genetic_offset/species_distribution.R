##R4.05, 4.2.3
#ref:https://jcoliver.github.io/learn-r/011-species-distribution-models.html
##ref:https://rspatial.org/raster/sdm/2_sdm_occdata.html


library(dismo)
library(raster)
library(geodata)


bioclim_data <-list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_2.5m_bio\\", pattern = "*.tif$", full.names=TRUE)
bioclim_data <- stack(bioclim_data)
predNames <- c('wc2.1_2.5m_bio_19','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18') ##only use some bioclim data
bioclim_data=bioclim_data[[predNames]]


obs_data <- read.table(file = "use_medicago_sativa_USDA_ARS.txt",header = T,sep = "\t")  ###我用的是比较保守的GRIN-ARS网站信息，主要使用经纬度信息
obs_data=obs_data[,c(1,2)] #only use latitude and longitude
obs_data397 <- read.table(file =  "397sample_info.txt",header = T,sep = "\t")  ##采样点信息
# Determine geographic extent of our data (whole world)
max_lat = 90
min_lat = -60
max_lon = 180
min_lon = -170
geographic_extent <- extent(x = c(min_lon, max_lon, min_lat, max_lat))
# Crop the bioclim data to geographic extent of my data
bioclim_data <- raster::crop(x = bioclim_data, y = geographic_extent)


tif_files <- list.files(path = "D:\\Population_643\\all_761_info\\local_adaption\\wc2.1_2.5m_bio\\",
                        pattern = "*.tif$",
                        full.names = TRUE)


# We only need one file, so use the first one in the list of .tif files
mask <- raster(tif_files[1])
set.seed(1)
# Randomly sample points, not overlap with our present data
background <- randomPoints(mask = mask,     # Provides resolution of sampling points
                           n = nrow(obs_data),      # Number of random points
                           ext = geographic_extent, # Spatially restricts sampling
                           extf = 1.25)             # Expands sampling a little bit
# Download data to use for our base map
world_map <- world(resolution = 3,
                   path = "data/")


# Crop the map to our area of interest
my_map <- raster::crop(x = world_map, y = geographic_extent)


##test background point location
pdf(file = "background_point_falcata_distribution.pdf",   width = 7.5, height = 4)
##plot whole world
plot(my_map,col='lightgray', border='lightgray',xlab="Longitude",ylab="Latitude" )


points(background, col = "grey30", pch = 1, cex = 0.6)
# And add those observations
points(x = obs_data$Longitude, y = obs_data$Latitude, col = "red",pch = 2, cex = 0.2)
dev.off()


# Arbitrarily assign group 1 as the testing data group
testing_group <- 1


# Create vector of group memberships,test using random point, real data need use representative data
group_presence <- kfold(x = obs_data, k = 5) # kfold is in dismo package


# Separate observations into training and testing groups
presence_train <- obs_data[group_presence != testing_group, ]
presence_test <- obs_data[group_presence == testing_group, ]


# Repeat the process for pseudo-absence points
group_background <- kfold(x = background, k = 5)


background_train <- background[group_background != testing_group, ]
background_test <- background[group_background == testing_group, ]

##bioclim模型
# Build a model using training data, 
bc_model <- bioclim(x = bioclim_data, p = presence_train)
# Predict presence from model
predict_presence <- dismo::predict(object = bc_model,
                                   x = bioclim_data,
                                   ext = geographic_extent)
# Use testing data for model evaluation
bc_eval <- evaluate(p = presence_test,   # The presence testing data
                    a = background_test, # The absence testing data
                    model = bc_model,    # The model we are evaluating
                    x = bioclim_data)    # Climatic variables for use by model
# Determine minimum threshold for "presence"
bc_threshold <- threshold(x = bc_eval, stat = "spec_sens")

##present distribution
pdf(file = "bioclim_present_distribution2.pdf",   width = 7.5, height = 4)
##plot whole world
plot(my_map, border='gray',col='gray',xlab="Longitude",ylab="Latitude" )
##add distribution region
plot(predict_presence > bc_threshold, legend = FALSE, add = TRUE,col = c(NA, "olivedrab"))
# And add those observations
points(x = obs_data$Longitude, y = obs_data$Latitude, col = "red",pch = 20, cex = 0.2)
dev.off()

##future
forecast_data <- brick(x = "wc2.1_10m_bioc_ACCESS-CM2_ssp245_2061-2080.tif")
names(forecast_data) <- c('wc2.1_2.5m_bio_1','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_3','wc2.1_2.5m_bio_4','wc2.1_2.5m_bio_5','wc2.1_2.5m_bio_6','wc2.1_2.5m_bio_7','wc2.1_2.5m_bio_8','wc2.1_2.5m_bio_9','wc2.1_2.5m_bio_10','wc2.1_2.5m_bio_11','wc2.1_2.5m_bio_12','wc2.1_2.5m_bio_13','wc2.1_2.5m_bio_14','wc2.1_2.5m_bio_15','wc2.1_2.5m_bio_16','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18','wc2.1_2.5m_bio_19')
forecast_presence <- dismo::predict(object = bc_model,
                                    x = forecast_data,
                                    ext = geographic_extent)

pdf(file = "bioclim_falcata_ssp245_2061-2080_distribution.pdf",   width = 7.5, height = 4)
##plot whole world
plot(my_map, border='gray',col='gray',xlab="Longitude",ylab="Latitude" )
# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast_presence > bc_threshold,
     add = TRUE,
     legend = FALSE,
     col = c(NA, "olivedrab"))
dev.off()
##part 2,
###linear model prepare data
colnames(background_train)=c("Longitude","Latitude")
train <- rbind(presence_train, background_train)
pb_train <- c(rep(1, nrow(presence_train)), rep(0, nrow(background_train)))
envtrain <- raster::extract(bioclim_data, train[,c("Longitude","Latitude")])
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
envtrain=na.omit(envtrain) ##remove na
#colnames(envtrain)=c("pa","bio1",'bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9')
testpres <- data.frame( extract(bioclim_data, presence_test) )
testbackg <- data.frame( extract(bioclim_data, background_test) )
testpres=na.omit(testpres)
testbackg=na.omit(testbackg)

##2.1.1 glm model
gm2 <- glm(pa ~ .,
           family = gaussian(link = "identity"), data=envtrain)
ge2 <- evaluate(testpres, testbackg, gm2)
pg <- predict(bioclim_data,gm2, ext=geographic_extent)
tr <- threshold(ge2, 'spec_sens')
pdf(file = "glm_point_distribution.pdf",   width = 7.5, height = 4)
##plot whole world
plot(pg > tr,main='GLM model',xlab="Longitude",ylab="Latitude",legend = FALSE)
# And add those observations
points(x = obs_data$Longitude, y = obs_data$Latitude, col = "red",pch = 20, cex = 0.6)
dev.off()

##2.1.2random forest
library(randomForest)
rf1 <- randomForest(pa~wc2.1_2.5m_bio_19+wc2.1_2.5m_bio_8+wc2.1_2.5m_bio_9+wc2.1_2.5m_bio_15+wc2.1_2.5m_bio_2+wc2.1_2.5m_bio_17+wc2.1_2.5m_bio_18, data=envtrain)
erf <- evaluate(testpres, testbackg, rf1)
erf
pr <- dismo::predict(bioclim_data, rf1, ext=geographic_extent) ##need some time
tr <- threshold(erf, 'spec_sens')
pdf(file = "RF_model_present_distribution_new_use.pdf",   width = 12, height = 6)
plot(pr > (tr*5.5), main='',xlab="Longitude",ylab="Latitude",legend = FALSE,col=c("gray95","lightsteelblue1"),ylim=c(-60,90),xlim=c(-180,180)) ##这个分布阈值可以根据自己的需要更改，找到合理的阈值即可（个人感慨：这个阈值也是根据不同人的偏好设置，所以单纯根据气候模型预测分布区域存在局限性，结合基因型数据更靠谱）
dev.off()

##2.1.3 svm Support Vector Machines model
library(kernlab)
svm <- ksvm(pa ~ ., data=envtrain)
esv <- evaluate(testpres, testbackg, svm)
esv
ps <- predict(bioclim_data, svm, ext=geographic_extent) ##may need long time
tr <- threshold(esv, 'spec_sens')
pdf(file = "SVM_model_present_distribution.pdf",   width = 7.5, height = 4)
plot(ps > tr, main='SVM model',xlab="Longitude",ylab="Latitude",legend = FALSE)
# And add those observations
points(x = obs_data397$Longitude, y = obs_data397$Latitude, col = "red",pch = 20, cex = 0.6)
dev.off()

