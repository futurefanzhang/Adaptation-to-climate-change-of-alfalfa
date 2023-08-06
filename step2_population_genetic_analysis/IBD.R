IBD<-read.table("IBD_702_result.genome",header=TRUE,sep="")
group_info=read.table("group.txt",header=TRUE,sep="")
colnames(group_info)=c("IID1","group")
library(dplyr)
combine1=inner_join(IBD, group_info, by="IID1")
colnames(group_info)=c("IID2","group2")
combine2=inner_join(combine1, group_info, by="IID2")
result=aggregate(x=combine2$PI_HAT,by=list(combine2$group,combine2$group2),FUN=mean)
colnames(result)=c("group1","group2","mean_PI_HAT")
write.csv(result,"IBD_result.csv",row.names = F)      ##used for cytoscape 
##I use cytoscape for plot, it is a window interface software. please refer:https://cytoscape.org/manual/Cytoscape3_7_0Manual.pdf.
