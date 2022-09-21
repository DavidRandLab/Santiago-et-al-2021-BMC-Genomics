###############################################################################
#### Heatmap of sig. DEGs in OXPHOS KEGG categories
###############################################################################
library(KEGGREST)
library(gplots)

##ID conversion table
convert=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/FBgnConversionTable.csv",row.names=1)

## generated in Figure 2 A-B code
clusterdata=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/MBcluster.seq%20Output.csv",row.names=1)

#### Figure 3A
kegg=keggLink("pathway","dme")
temp=names(kegg)
temp2=convert[convert[,5]!="",]
i=1
while(i<=nrow(temp2)){
  if(length(temp[temp==temp2[i,5]])){
    temp[temp==temp2[i,5]]=as.character(temp2[i,1])
  }
  i=i+1
}
names(kegg)=temp
OXPHOS=names(kegg[grep("00190",kegg)])
OXPHOS=intersect((OXPHOS),row.names(clusterdata))
OXPHOS=setdiff(OXPHOS,convert[convert$mt==1,1])
hmcolors=clusterdata[OXPHOS,"Cluster"]
hmdata=clusterdata[OXPHOS,1:16]
row.names(hmdata)=as.vector(convert[OXPHOS,"Symbol"])
hmorder=order(hmcolors, decreasing=T)
hmdata=hmdata[hmorder,]
hmcolors=hmcolors[hmorder]
cluster=paste("cluster",hmcolors)
basecolors=c("red","green3","black","yellow","blue")
hmcolors=basecolors[(hmcolors)]

col.pan <- colorpanel(100, "black","white","red")
dev.new()
heatmapped=heatmap.2(as.matrix(hmdata),Rowv=F,Colv=F,RowSideColors = hmcolors,trace="none",col=col.pan,scale="row",key=TRUE,cexRow = .3,cexCol = 1,main="OXPHOS DEGs")

#### Figure 3B
##Cluster Line Graph
##gets the mean cluster data for line plots
rawclusterdata=as.matrix(clusterdata[OXPHOS,])
rawclusterdata=rawclusterdata[rawclusterdata[,"Cluster"]==1,]
temp=convert$Symbol
names(temp)=convert$original.FBgn
temp=temp[row.names(rawclusterdata)]
rawclusterdata=hmdata[temp,]
OOC=c(mean(rawclusterdata[,1]),mean(rawclusterdata[,2]),mean(rawclusterdata[,3]),mean(rawclusterdata[,4]))
OOR=c(mean(rawclusterdata[,5]),mean(rawclusterdata[,6]),mean(rawclusterdata[,7]),mean(rawclusterdata[,8]))
SOC=c(mean(rawclusterdata[,9]),mean(rawclusterdata[,10]),mean(rawclusterdata[,11]),mean(rawclusterdata[,12]))
SOR=c(mean(rawclusterdata[,13]),mean(rawclusterdata[,14]),mean(rawclusterdata[,15]),mean(rawclusterdata[,16]))
timelabs=c("0H", "1H", "2H", "4H")

dev.new()
i=1
plot(x=c(1:4),rawclusterdata[i,c(1:4)],type="l",col="grey", ylim=c(min(rawclusterdata[,c(1:14)]),max(rawclusterdata[,c(1:14)])),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,16),xlab="")
title(main="Cluster 1 Data")
title(xlab="Samples Over Time",line=4) 
i=i+1
while(i<=NROW(rawclusterdata)){
  lines(x=c(1:4),rawclusterdata[i,c(1:4)],type="l",col="grey")
  lines(x=c(5:8),rawclusterdata[i,c(5:8)],type="l",col="grey")
  lines(x=c(9:12),rawclusterdata[i,c(9:12)],type="l",col="grey")
  lines(x=c(13:16),rawclusterdata[i,c(13:16)],type="l",col="grey")
  i=i+1
}
line1=OOC
line2=OOR
line3=SOC
line4=SOR
lines(x=c(1:4),line1,type="l",col="blue",lwd=2)
lines(x=c(5:8),line2,type="l",col="red",lwd=2)
lines(x=c(9:12),line3,type="l",col="blue",lwd=2)
lines(x=c(13:16),line4,type="l",col="red",lwd=2)
abline(h=0, col="black", lty=2, lwd=1)
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(5:8),labels=timelabs,col.axis="Black")
axis(1,at=c(9:12),labels=timelabs,col.axis="Black")
axis(1,at=c(13:16),labels=timelabs,col.axis="Black")
axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="blue")
axis(1,at=c(6.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="red")
axis(1,at=c(10.5),labels=c("Control"),padj=2,tick=F,col.axis="blue")
axis(1,at=c(14.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="red")
axis(1,at=c(4.55),labels=c("OreR;OreR"),padj=4,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR"),padj=4,tick=F,col.axis="black")

#### Figure 3C
##Cluster Line Graph
##gets the mean cluster data for line plots
rawclusterdata=as.matrix(clusterdata[OXPHOS,])
rawclusterdata=rawclusterdata[rawclusterdata[,"Cluster"]==5,]
temp=convert$Symbol
names(temp)=convert$original.FBgn
temp=temp[row.names(rawclusterdata)]
rawclusterdata=hmdata[temp,]
OOC=c(mean(rawclusterdata[,1]),mean(rawclusterdata[,2]),mean(rawclusterdata[,3]),mean(rawclusterdata[,4]))
OOR=c(mean(rawclusterdata[,5]),mean(rawclusterdata[,6]),mean(rawclusterdata[,7]),mean(rawclusterdata[,8]))
SOC=c(mean(rawclusterdata[,9]),mean(rawclusterdata[,10]),mean(rawclusterdata[,11]),mean(rawclusterdata[,12]))
SOR=c(mean(rawclusterdata[,13]),mean(rawclusterdata[,14]),mean(rawclusterdata[,15]),mean(rawclusterdata[,16]))
timelabs=c("0H", "1H", "2H", "4H")

dev.new()
i=1
plot(x=c(1:4),rawclusterdata[i,c(1:4)],type="l",col="grey", ylim=c(min(rawclusterdata[,c(1:14)]),max(rawclusterdata[,c(1:14)])),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,16),xlab="")
title(main="Cluster 1 Data")
title(xlab="Samples Over Time",line=4) 
i=i+1
while(i<=NROW(rawclusterdata)){
  lines(x=c(1:4),rawclusterdata[i,c(1:4)],type="l",col="grey")
  lines(x=c(5:8),rawclusterdata[i,c(5:8)],type="l",col="grey")
  lines(x=c(9:12),rawclusterdata[i,c(9:12)],type="l",col="grey")
  lines(x=c(13:16),rawclusterdata[i,c(13:16)],type="l",col="grey")
  i=i+1
}
line1=OOC
line2=OOR
line3=SOC
line4=SOR
lines(x=c(1:4),line1,type="l",col="blue",lwd=2)
lines(x=c(5:8),line2,type="l",col="red",lwd=2)
lines(x=c(9:12),line3,type="l",col="blue",lwd=2)
lines(x=c(13:16),line4,type="l",col="red",lwd=2)
abline(h=0, col="black", lty=2, lwd=1)
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(5:8),labels=timelabs,col.axis="Black")
axis(1,at=c(9:12),labels=timelabs,col.axis="Black")
axis(1,at=c(13:16),labels=timelabs,col.axis="Black")
axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="blue")
axis(1,at=c(6.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="red")
axis(1,at=c(10.5),labels=c("Control"),padj=2,tick=F,col.axis="blue")
axis(1,at=c(14.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="red")
axis(1,at=c(4.55),labels=c("OreR;OreR"),padj=4,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR"),padj=4,tick=F,col.axis="black")