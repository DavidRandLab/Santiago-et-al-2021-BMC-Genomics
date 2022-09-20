library(MBCluster.Seq)
library(edgeR)
source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/Output.r")
source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/KmeansPlus.r")
source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/LGLK.r")
source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/Cluster.r")
source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/Math.r")
##source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/HybridTree-KM.r")
source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/Tree.r")
source("https://raw.githubusercontent.com/cran/MBCluster.Seq/master/R/Distance.r")

## Count table and groups with outlier libraries removed
counts=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Trimmed.Counts.csv",row.names=1)
groups=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Trimmed.Groups.csv",row.names=1)

##ImpulseDE2 output tables
oocc=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3E%20OOC-OOR.csv",row.names=1)
oorapa=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3B%20OOR.csv",row.names=1)
oocontrol=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3A%20OOC.csv",row.names=1)
socc=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3F%20SOC-SOR.csv",row.names=1)
sorapa=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3D%20SOR.csv",row.names=1)
socontrol=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3C%20SOC.csv",row.names=1)
oocsoc=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3G%20OOC-SOC.csv",row.names=1)
oorsor=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3H%20OOR-SOR.csv",row.names=1)

##ImpulseDE2 output tables with significant genes only
OCSC=na.omit(oocsoc[oocsoc$padj<.05,])
ORSR=na.omit(oorsor[oorsor$padj<.05,])
OO=na.omit(oocc[oocc$padj<.05,])
SO=na.omit(socc[socc$padj<.05,])
SOR=na.omit(sorapa[sorapa$padj<.05,])
SOC=na.omit(socontrol[socontrol$padj<.05,])
OOR=na.omit(oorapa[oorapa$padj<.05,])
OOC=na.omit(oocontrol[oocontrol$padj<.05,])

#### edgeR DEG analyses
##normalizing data with edgeR 
x<-DGEList(counts)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 
cpm=cpm(z)

###############################################################################
#### run mbclusterseq on all ImpulseDE1 DEGs using CPM normalized expression data 
###############################################################################
#### Figure 2A
##create heatmaps and files
allDEGs=unique(unlist(lapply(list(OCSC,ORSR,OO,SO,SOR,SOC,OOR,OOC),row.names)))
geneset=allDEGs

mbcounts=z$counts[geneset,]
mbgroups=z$samples
order=c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)

mbcounts=mbcounts[,order]
mbgroups=mbgroups[order,]
mbgroups$Time=c(rep(c(0,0,1,1,2,2,4,4),4))
mbgroups$Treat[c(9,10,25,26)]="Rapa"
mbgroups$Group=paste(substr(mbgroups$Sample,1,2),substr(mbgroups$Time,1,1),"H",substr(mbgroups$Treat,1,1),sep="")
GeneID=row.names(mbcounts)  
Normalizer=mbgroups$norm.factors
Treatment=mbgroups$Group
Treatment=paste(substr(Treatment,1,2),substr(Treatment,5,5),substr(Treatment,3,4),sep="")

mydata=RNASeq.Data(mbcounts,Normalizer,Treatment,GeneID) 

c0=KmeansPlus.RNASeq(mydata,nK=5)$centers

cls=Cluster.RNASeq(data=mydata,model="nbinom",centers=c0,method="DA")$cluster

tr=Hybrid.Tree(data=mydata,cluster0 =cls,model="nbinom") 

image=plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)

clusterdata=cbind(mydata[[6]],cls) 
colnames(clusterdata)=c(levels(as.factor(Treatment)),"Cluster")

dev.new(height=7.7, width=6.6, noRStudioGD = TRUE, units = "inch")
plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)




#### Figure 2B
clusterdata=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/MBcluster.seq%20Output.csv",row.names=1)

cluster1=clusterdata[clusterdata$Cluster==1,]
line1=c(rep(0,16))
i=1
while(i<17){
  line1[i]=mean(cluster1[,i])
  i=i+1
}

cluster2=clusterdata[clusterdata$Cluster==2,]
line2=c(rep(0,16))
i=1
while(i<17){
  line2[i]=mean(cluster2[,i])
  i=i+1
}
allclusters=rbind(line1,line2)

cluster3=clusterdata[clusterdata$Cluster==3,]
line3=c(rep(0,16))
i=1
while(i<17){
  line3[i]=mean(cluster3[,i])
  i=i+1
}
allclusters=rbind(allclusters,line3)

cluster4=clusterdata[clusterdata$Cluster==4,]
line4=c(rep(0,16))
i=1
while(i<17){
  line4[i]=mean(cluster4[,i])
  i=i+1
}
allclusters=rbind(allclusters,line4)

cluster5=clusterdata[clusterdata$Cluster==5,]
line5=c(rep(0,16))
i=1
while(i<17){
  line5[i]=mean(cluster5[,i])
  i=i+1
}
allclusters=rbind(allclusters,line5)
limits=c(min(allclusters),max(allclusters))
allclusters=allclusters[,c(1:4,NA,5:8,NA,9:12,NA,13:16)]
colnames(allclusters)=c("OOC0H","OOC1H","OOC2H","OOC4H",NA,"OOR0H","OOR1H","OOR2H","OOR4H",NA,"SOC0H","SOC1H","SOC2H","SOC4H",NA,"SOR0H","SOR1H","SOR2H","SOR4H")


##Actually plots it. plots the 
dev.new(height=7.7, width=10, noRStudioGD = TRUE, units = "inch")
par(mar=c(5,5,1,1),bg="transparent")
plot(allclusters[1,],col="red", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=limits,lwd=2)
title(main="")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(11:14),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(16:19),labels=c(0,1,2,4),col.axis="Black")

axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="Blue")
axis(1,at=c(7.5),labels=c("Rapa"),padj=2,tick=F,col.axis="Red")
axis(1,at=c(12.5),labels=c("Control"),padj=2,tick=F,col.axis="Blue")
axis(1,at=c(17.5),labels=c("Rapa"),padj=2,tick=F,col.axis="Red")

axis(1,at=c(5),labels=c("OreR;OreR"),padj=3,tick=F,col.axis="black")
axis(1,at=c(15),labels=c("sm21;OreR"),padj=3,tick=F,col.axis="black")

lines(allclusters[2,],type="l",col="green",lwd=2)
lines(allclusters[3,],type="l",col="black",lwd=2)
lines(allclusters[4,],type="l",col="yellow",lwd=2)
lines(allclusters[5,],type="l",col="blue",lwd=2)

abline(h=0, col="black", lty=2, lwd=1)