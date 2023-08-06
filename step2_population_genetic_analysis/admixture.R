library(ggplot2)
mydata1 <- read.table("CV_error.txt", head = T,sep="") ##CV error info
mydata1$K_value<-as.factor(mydata1$K_value)
p_k<-ggplot(data=mydata1, aes(x=K_value, y=CV_error, group=1)) +
  geom_line(color="orange")+theme_gray(base_size = 6)+
  geom_point(color="darkgreen")+labs(x = 'K value',y= "CV error",
                                     title = '')+
  theme(plot.margin = margin(t=1, r=1, b=1, l=1, "cm"),
        axis.text = element_text(color="black",size=10),axis.line = element_line(colour = "black", size = 1),
        axis.title.y = element_text(color="black",size=10,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(color="black",size=10,margin = margin(t = 10, r = 0, b = 0, l = 0)), #top(t),right(r),bottom(b),left(l)
        legend.title = element_text(color = "black", size = 10), 
        legend.text = element_text(color = "black", size = 10),
        legend.position = "none")
ggsave(plot=p_k,"K_value_CV_error.pdf", width = 15, height = 6)

##admixture plot
library(pophelper)
file_list_admix <- list.files("./admixture2_10/", pattern = ".Q", full.names = T)
#file_list_admix <- list.files("./admixture2_10/", pattern = ".9.Q", full.names = T) ##only use k=9, corresponding ##single K part
info <- read.table("sample_order.txt", header = T, stringsAsFactors = F) #less prunData.ped |awk '{print $1,$2}' >sample_order.txt , you can use these command extract sample name
label_set<-read.table("sample_label_group.txt",sep="\t", header=T,stringsAsFactors=F) #input geography info
one_label<-label_set[,2,drop=FALSE]
colnames(one_label)<-" " #remove the header of geography info
qlist_admix <- readQ(file_list_admix)
for(i in 1:length(qlist_admix)){
  row.names(qlist_admix[[i]]) <- info$sample1
}
order<-label_set[,1,drop=FALSE] #get individual order
index=match(info$sample1,order$taxa) #match initial order with use order
for(i in 1:length(qlist_admix)){
  qlist_admix[[i]] <- qlist_admix[[i]][order(index),]
}
k_order <- vector()
mink <- 1
maxk <- 10 
if(maxk < 10){
  k_order <- 1:length(qlist_admix)
} else if (maxk < 20) {
  end1 <- maxk - 10 + 1
  start2 <- end1 + 1
  k_order <- c(start2:length(qlist_admix), 1:end1)
}
klab <- vector()
if(mink == 1){
  klab <- 2:maxk
} else {
  klab <- mink:maxk
}

##all k from 2 to 10
prefix <- "admix_Q2_10_group702"
height <- 1.2
width <- 15
plotQ(qlist_admix[k_order],imgoutput="join",showindlab=F,showlegend=F, sortind = "label",panelratio=c(1.9,2.1),
      #indlabsize=2,indlabheight=0,indlabspacer=0.05,barbordersize=NA,clustercol=c("#EE3435","#FABB2E","#F79331","#22B45F","#8AC666","#A7432C","#83CABA","#A983B6","#69459C","#4288AD"),
      indlabsize=2,indlabheight=0,indlabspacer=0.05,barbordersize=NA,clustercol=c("#72B33D","#FCCC00","#F05C80","#00827C","#EB051F","#B0D882","#FBA117","#5E4886","#FCB887","#CD5530"),
      outputfilename=prefix,imgtype="pdf", sharedindlab = F,splabcol = "black",grplabsize = 2.5,grplabpos=0.7,grplabheight=1,grplabjust=1,grplab=one_label,grplabangle=30,grplabcol = "black",
      useindlab = F, showyaxis = F, basesize = 10, sppos = "right", showticks = F,
      splab = paste0("K = ",klab), divsize = 0.5,splabsize = 10,exportpath = getwd(),
      width = width, height = height, panelspacer = 0.02, dpi = 600, barbordercolour = NA)

##single K
prefix <- "admix_Q9_group702"
height <- 1.5
width <- 15
plotQ(qlist_admix[k_order],imgoutput="sep",showindlab=F,showlegend=F, sortind = "label",panelratio=c(1.9,2.1),
      #indlabsize=2,indlabheight=0,indlabspacer=0.05,barbordersize=NA,clustercol=c("#EE3435","#FABB2E","#F79331","#22B45F","#8AC666","#A7432C","#83CABA","#A983B6","#69459C","#4288AD"),
      indlabsize=2,indlabheight=0,indlabspacer=0.05,barbordersize=NA,clustercol=c("#72B33D","#FCCC00","#F05C80","#00827C","#EB051F","#B0D882","#FBA117","#5E4886","#FCB887","#CD5530"),
      outputfilename=prefix,imgtype="pdf", sharedindlab = F,splabcol = "black",grplabsize = 2,grplabpos=0.7,grplabheight=1,grplabjust=1,grplab=one_label,grplabangle=30,grplabcol = "black",
      useindlab = F, showyaxis = F, basesize = 8, sppos = "right", showticks = F,
      splab = paste0("K = 9"), divsize = 0.5,splabsize = 8,exportpath = getwd(),
      width = width, height = height, panelspacer = 0.02, dpi = 600, barbordercolour = NA)
