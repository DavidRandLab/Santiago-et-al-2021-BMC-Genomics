library(visNetwork)
library(goseq)

##C1.KEGG.DEGs
clusterdata=read.csv("/Users/johncsantiago/Google Drive/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv",row.names=1)

convert=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/FBgnConversionTable.csv",row.names=1)

kegg=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/FBgn2KEGG.csv")

cpm=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Trimmed.CPM.data.csv",row.names=1)

##ImpulseDE2 output tables
oocsoc=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3G%20OOC-SOC.csv",row.names=1)
oorsor=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/ImpulseDE%20Outputs/Table%20S3H%20OOR-SOR.csv",row.names=1)

##ImpulseDE2 output tables with significant genes only
OCSC=na.omit(oocsoc[oocsoc$padj<.05,])
ORSR=na.omit(oorsor[oorsor$padj<.05,])

tempkegg=as.vector(kegg[,2])
names(tempkegg)=convert[as.vector(kegg[,1]),"Symbol"]
kegg.sym=as.list(na.omit(tempkegg))

i=1
kegg.names=unique(substr(unlist(kegg),6,13))
names(kegg.names)=as.vector(kegg.names)
while(i<=length(kegg.names)){
  kegg.names[i] <- strsplit(keggGet(kegg.names[i])[[1]]$NAME," - ")[[1]][1]
  i=i+1
}

#### Figure 6B
Cluster1.Genes=row.names(convert[convert$cluster==1&convert$mt==0,])
Cluster1.Genes=intersect(Cluster1.Genes,c(row.names(ORSR),row.names(OCSC)))

temp.genes=cpm
temp.genes[,1]=0
temp.genes[Cluster1.Genes,1]=1
genes=temp.genes[,1]
names(genes)=convert[row.names(temp.genes),"Symbol"]
genes=(genes[-grep("jus",names(genes))[1]])
genes=genes[na.omit(names(genes))]

table(genes)
pwf=nullp(genes,"dm3","geneSymbol")
head(pwf)

##KEGG term enrichment is KEGG
KEGG=goseq(pwf,gene2cat=kegg.sym)
KEGG$adjp=p.adjust(KEGG$over_represented_pvalue,method="BH")
row.names(KEGG)=substr(as.character(KEGG[,1]),6,13)
KEGG$Name=kegg.names[row.names(KEGG)]
##Remove Global and Overview Maps
KEGG=KEGG[substr(row.names(KEGG),1,6)!="dme011",]
Cluster1.KEGG.Table=KEGG[substr(row.names(KEGG),1,6)!="dme012",]
sigCluster1.KEGG.Table=Cluster1.KEGG.Table[Cluster1.KEGG.Table$adjp<.05,]


#### Figure 6A
inKegg=(kegg[,2])
names(inKegg)=kegg[,1]
kegg.genes=rep(0,length(Cluster1.Genes))
names(kegg.genes)=Cluster1.Genes
i=1
while(i<=nrow(sigCluster1.KEGG.Table)){
  temp=grep(sigCluster1.KEGG.Table[i,"category"],inKegg)
  kegg.genes[names(inKegg[temp])]=1
  i=i+1
}

kegg.genes=names(kegg.genes[kegg.genes==1])
kegg.genes=intersect(kegg.genes,Cluster1.Genes)

rawclusterdata=as.matrix(clusterdata[kegg.genes,])
OOC=c(mean(rawclusterdata[,1]),mean(rawclusterdata[,2]),mean(rawclusterdata[,3]),mean(rawclusterdata[,4]))
OOR=c(mean(rawclusterdata[,5]),mean(rawclusterdata[,6]),mean(rawclusterdata[,7]),mean(rawclusterdata[,8]))
SOC=c(mean(rawclusterdata[,9]),mean(rawclusterdata[,10]),mean(rawclusterdata[,11]),mean(rawclusterdata[,12]))
SOR=c(mean(rawclusterdata[,13]),mean(rawclusterdata[,14]),mean(rawclusterdata[,15]),mean(rawclusterdata[,16]))
timelabs=c("0H", "1H", "2H", "4H")

sdev1=c(sd(rawclusterdata[,1]),sd(rawclusterdata[,2]),sd(rawclusterdata[,3]),sd(rawclusterdata[,4]))
sdev2=c(sd(rawclusterdata[,5]),sd(rawclusterdata[,6]),sd(rawclusterdata[,7]),sd(rawclusterdata[,8]))
sdev3=c(sd(rawclusterdata[,9]),sd(rawclusterdata[,10]),sd(rawclusterdata[,11]),sd(rawclusterdata[,12]))
sdev4=c(sd(rawclusterdata[,13]),sd(rawclusterdata[,14]),sd(rawclusterdata[,15]),sd(rawclusterdata[,16]))

dev.new()
plot(OOC,col="blue",ylim=c((min(c(OOR,OOC,SOC,SOR))-max(c(sdev1,sdev2,sdev3,sdev4))),(max(c(OOR,OOC,SOC,SOR))+max(c(sdev1,sdev2,sdev3,sdev4)))),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,9),xlab="",type="l",lwd=2)
lines(x=c(6:9),OOR,type="l",col="blue",lwd=2)
lines(x=c(1:4),SOC,type="l",col="red",lwd=2)
lines(x=c(6:9),SOR,type="l",col="red",lwd=2)
abline(h=0, col="black", lty=2, lwd=1)

x=c(1:4)
arrows(x, OOC-sdev1, x, OOC+sdev1, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, SOC-sdev3, x, SOC+sdev3, length=0.05, angle=90, code=3,col="red",lwd=2)
x=c(6:9)
arrows(x, OOR-sdev2, x, OOR+sdev2, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, SOR-sdev4, x, SOR+sdev4, length=0.05, angle=90, code=3,col="red",lwd=2)


title(main="Cluster 1 Genes")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")

axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="black")
legend("topleft",legend = c("OreR;OreR","sm21:OreR"),fill = c("blue","red"),cex = .5)
abline(h=0,lty=2)

#### Figure 6C



#### Figure 6E
Cluster5.Genes=row.names(convert[convert$cluster==5&convert$mt==0,])
Cluster5.Genes=intersect(Cluster5.Genes,c(row.names(ORSR),row.names(OCSC)))
temp.genes=cpm
temp.genes[,1]=0
temp.genes[Cluster5.Genes,1]=1
genes=temp.genes[,1]
names(genes)=convert[row.names(temp.genes),"Symbol"]
genes=(genes[-grep("jus",names(genes))[1]])
genes=genes[na.omit(names(genes))]

table(genes)
pwf=nullp(genes,"dm3","geneSymbol")
head(pwf)

##KEGG term enrichment is KEGG
KEGG=goseq(pwf,gene2cat=kegg.sym)
KEGG$adjp=p.adjust(KEGG$over_represented_pvalue,method="BH")
row.names(KEGG)=substr(as.character(KEGG[,1]),6,13)
KEGG$Name=kegg.names[row.names(KEGG)]
##Remove Global and Overview Maps
KEGG=KEGG[substr(row.names(KEGG),1,6)!="dme011",]
Cluster5.KEGG.Table=KEGG[substr(row.names(KEGG),1,6)!="dme012",]
sigCluster5.KEGG.Table=Cluster5.KEGG.Table[Cluster5.KEGG.Table$adjp<.05,]

#### Figure 6D
inKegg=(kegg[,2])
names(inKegg)=kegg[,1]
kegg.genes=rep(0,length(Cluster5.Genes))
names(kegg.genes)=Cluster5.Genes
i=1
while(i<=nrow(sigCluster5.KEGG.Table)){
  temp=grep(sigCluster5.KEGG.Table[i,"category"],inKegg)
  kegg.genes[names(inKegg[temp])]=1
  i=i+1
}

kegg.genes=names(kegg.genes[kegg.genes==1])
kegg.genes=intersect(kegg.genes,Cluster5.Genes)

rawclusterdata=as.matrix(clusterdata[kegg.genes,])

OOC=c(mean(rawclusterdata[,1]),mean(rawclusterdata[,2]),mean(rawclusterdata[,3]),mean(rawclusterdata[,4]))
OOR=c(mean(rawclusterdata[,5]),mean(rawclusterdata[,6]),mean(rawclusterdata[,7]),mean(rawclusterdata[,8]))
SOC=c(mean(rawclusterdata[,9]),mean(rawclusterdata[,10]),mean(rawclusterdata[,11]),mean(rawclusterdata[,12]))
SOR=c(mean(rawclusterdata[,13]),mean(rawclusterdata[,14]),mean(rawclusterdata[,15]),mean(rawclusterdata[,16]))
timelabs=c("0H", "1H", "2H", "4H")

sdev1=c(sd(rawclusterdata[,1]),sd(rawclusterdata[,2]),sd(rawclusterdata[,3]),sd(rawclusterdata[,4]))
sdev2=c(sd(rawclusterdata[,5]),sd(rawclusterdata[,6]),sd(rawclusterdata[,7]),sd(rawclusterdata[,8]))
sdev3=c(sd(rawclusterdata[,9]),sd(rawclusterdata[,10]),sd(rawclusterdata[,11]),sd(rawclusterdata[,12]))
sdev4=c(sd(rawclusterdata[,13]),sd(rawclusterdata[,14]),sd(rawclusterdata[,15]),sd(rawclusterdata[,16]))

dev.new()
plot(OOC,col="blue",ylim=c((min(c(OOR,OOC,SOC,SOR))-max(c(sdev1,sdev2,sdev3,sdev4))),(max(c(OOR,OOC,SOC,SOR))+max(c(sdev1,sdev2,sdev3,sdev4)))),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,9),xlab="",type="l",lwd=2)
lines(x=c(6:9),OOR,type="l",col="blue",lwd=2)
lines(x=c(1:4),SOC,type="l",col="red",lwd=2)
lines(x=c(6:9),SOR,type="l",col="red",lwd=2)
abline(h=0, col="black", lty=2, lwd=1)

x=c(1:4)
arrows(x, OOC-sdev1, x, OOC+sdev1, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, SOC-sdev3, x, SOC+sdev3, length=0.05, angle=90, code=3,col="red",lwd=2)
x=c(6:9)
arrows(x, OOR-sdev2, x, OOR+sdev2, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, SOR-sdev4, x, SOR+sdev4, length=0.05, angle=90, code=3,col="red",lwd=2)


title(main="Cluster 5 Genes")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")

axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="black")
legend("topleft",legend = c("OreR;OreR","sm21:OreR"),fill = c("blue","red"),cex = .5)
abline(h=0,lty=2)
