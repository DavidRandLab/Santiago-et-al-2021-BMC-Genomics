library(parallel)
library(BiocParallel)
library(DESeq2)
library(compiler)
library(SummarizedExperiment)

##MDS distance matrix to exclude outlier replicate libraries
counts=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Full%20Count%20Table.csv",row.names=1)
groups=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Full%20Count%20Table%20Meta%20Data.csv",row.names=1)
counts=counts[,row.names(groups)]

##Generating filter vector called trim that removes the most outlying replicate from each treatment group
##Using distance between libraries in an MDS plot to determine which library 

##uses the edger function that is defined below (needs to be ran before using) to trim out the genes with low read counts
x<-DGEList(counts[,row.names(groups)])
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 

MDSdata=z$counts
MDSgroups=groups

lcpm=cpm(MDSdata,log=TRUE)
colors=factor(groups$Group)
label=groups$Group
levels(colors)=c(1:length(levels(colors)))
temp=plotMDS(lcpm, col=as.vector(colors),labels=label,gene.selection="pairwise",top=nrow(lcpm))

similarity=temp$distance.matrix
samples=unique(groups$Group)

similaritytable=groups[,c(6,5)]
similaritytable$rep1=c(rep(c(1,1,2),14))
similaritytable$rep2=c(rep(c(2,3,3),14))
similaritytable$Replicate=c(rep(c("1~2","1~3","2~3"),14))
similaritytable$distance=0
i=1
j=1
while(i<=length(samples)){
  coords=groups$Group==samples[i]
  temp=similarity[coords,coords]
  similaritytable$distance[j]=temp[2,1]
  j=j+1
  similaritytable$distance[j]=temp[3,1]
  j=j+1
  similaritytable$distance[j]=temp[3,2]
  j=j+1
  i=i+1
}

tempkeep=similaritytable[similaritytable$Group==samples[1],]
tempkeep=tempkeep[tempkeep$distance==min(tempkeep$distance),]
keep=tempkeep
i=2
while(i<=length(samples)){
  tempkeep=similaritytable[similaritytable$Group==samples[i],]
  tempkeep=tempkeep[tempkeep$distance==min(tempkeep$distance),]
  keep=rbind(keep,tempkeep)
  i=i+1
}

trim=c(rep(0,42))
i=1
while(i<15){
  trim=trim+(groups$Group==keep[i,1]&groups$Replicate==keep[i,3])
  trim=trim+(groups$Group==keep[i,1]&groups$Replicate==keep[i,4])
  i=i+1
}

trim=trim==1
names(trim)=groups$Group

################
################

##ImpulseDE2
i=1
while(i <= length(list.files("/Users/johncsantiago/Documents/ImpulseDE2/ImpulseDE2/R/"))){
  source(paste("/Users/johncsantiago/Documents/ImpulseDE2/ImpulseDE2/R/",list.files("/Users/johncsantiago/Documents/ImpulseDE2/ImpulseDE2/R/")[i],sep="")) 
  i=i+1
}


counts=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Full%20Count%20Table.csv",row.names=1)
groups=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Full%20Count%20Table%20Meta%20Data.csv",row.names=1)
counts=counts[,row.names(groups)]
counts=counts[,trim]
groups=groups[trim,]
order=c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)
counts=counts[,order]
groups=groups[order,]
genotype=substr(groups$Mito,1,1)=="O"
counts=counts[,genotype]
groups=groups[genotype,]

counts=as.matrix(counts[1:100,])

##Comparing Control to Rapamycin treated
counts=counts[,groups$Mito=="OreR"]
groups=groups[groups$Mito=="OreR",]
groups[,"Condition"]=c(rep("control",8),rep("case",8))
groups[,"Sample"]=colnames(counts)
groups[,"Time"]=c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4)
groups[,"Batch"]=c(rep(c(1,2),8))
row.names(groups)=groups$Sample
groups=groups[,c(1,7,4,9)]

##running impulse de2 without including transiently expressed genes
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = groups,
  boolCaseCtrl    = TRUE,
  vecConfounders  = NULL,
  scaNProc        = 1 )


##write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCountOreRCC_TrimmedLibraries_ImpulseDEResuults.csv")
oocc=read.csv("SingleCountOreRCC_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)







colnames(cpmdata)=groups[colnames(cpmdata), "Group"]

heatmaply(cpmdata[cluster5,c(1:4,7,8,11,12,15:18,21,22,25,26,1,2,5,6,9,10,13,14,15,16,19,20,23,24,27,28)], scale="row",Colv = FALSE)
heatmaply(cpmdata[cluster5,c(1,2,5,6,9,10,13,14,15,16,19,20,23,24,27,28)], scale="row",Colv = FALSE)
colnames(cpmdata)=groups[colnames(cpmdata), "Group"]
heatmaply(cpmdata[intersect(cluster5,row.names(cpmdata)),c(1:3,7:9,13:15,19:21,22:24,28:30,34:36,40:42)], scale="row",Colv = FALSE)


ywd="/Users/johncsantiago/Google Drive/My Drive/Santiago_RandLab_DigitalNotebook/FlyRefeedingPaper/"
setwd(ywd)

boxplot(as.numeric(cpmdata["FBgn0003651",])~groups$Group, main="svp")
cpmdata=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)

so1=read.csv("/Users/johncsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/GenomeBiologyPaper/Supplementary/Table S2J SO1HR.SO0HC.csv",row.names=1)

groups=groups
groups=as.matrix(groups)
groups=as.data.frame(groups)
x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 

so1.c5=so1[row.names(cluster5),]
so1.c5=so1.c5[order(so1.c5$FDR),]

i=1
barplot(height = as.numeric(cpmdata[row.names(so1.c5)[i],]),names.arg = groups$Group,las=2,col=c(rep("grey",14),"blue","blue","grey","grey","red","red",rep("grey",10)),main=row.names(so1.c5)[i])
row.names(so1.c5)[i]
i=i+1

##Cluster Line Graph
##gets the mean cluster data for line plots
clusterdata=read.csv("/Users/johncsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv",row.names=1)
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
##Cluster 1
par(mar=c(5,5,1,1),bg="white")
plot(allclusters[1,c(1:9)],col="red", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=c((1.5*min(na.omit(allclusters[1,]))),(1.75*max(na.omit(allclusters[1,])))),lwd=2)
title(main="Cluster 1")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(2.5),labels=c("Control"),padj=3,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=3,tick=F,col.axis="black")
lines(allclusters[1,c(11:19)],type="l",lty=2,col="red",lwd=2)
##legend("topleft", legend = c("OreR;OreR","sm21;OreR"), col = c("red"), lty= c(1,2), lwd = 1, cex=.7)
abline(h=0, col="black", lty=2, lwd=1)

##Cluster2
par(mar=c(5,5,1,1),bg="white")
plot(allclusters[2,c(1:9)],col="green3", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=c((1.5*min(na.omit(allclusters[2,]))),(1.5*max(na.omit(allclusters[2,])))),lwd=2)
title(main="Cluster 2")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(2.5),labels=c("Control"),padj=3,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=3,tick=F,col.axis="black")
lines(allclusters[2,c(11:19)],type="l",lty=2,col="green3",lwd=2)
##legend("topleft", legend = c("OreR;OreR","sm21;OreR"), col = c("green3"), lty= c(1,2), lwd = 1, cex=.7)
abline(h=0, col="black", lty=2, lwd=1)

##Cluster3
par(mar=c(5,5,1,1),bg="white")
plot(allclusters[3,c(1:9)],col="black", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=c((1.5*min(na.omit(allclusters[3,]))),(1.5*max(na.omit(allclusters[3,])))),lwd=2)
title(main="Cluster 3")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(2.5),labels=c("Control"),padj=3,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=3,tick=F,col.axis="black")
lines(allclusters[3,c(11:19)],type="l",lty=2,col="black",lwd=2)
##legend("topleft", legend = c("OreR;OreR","sm21;OreR"), col = c("black"), lty= c(1,2), lwd = 1, cex=.7)
abline(h=0, col="black", lty=2, lwd=1)

##Cluster4
par(mar=c(5,5,1,1),bg="white")
plot(allclusters[4,c(1:9)],col="gold3", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=c((1.5*min(na.omit(allclusters[4,]))),(1.5*max(na.omit(allclusters[4,])))),lwd=2)
title(main="Cluster 4")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(2.5),labels=c("Control"),padj=3,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=3,tick=F,col.axis="black")
lines(allclusters[4,c(11:19)],type="l",lty=2,col="gold3",lwd=2)
##legend("topleft", legend = c("OreR;OreR","sm21;OreR"), col = c("gold3","gold3"), lty= c(1,2), lwd = 1, cex=.7)
abline(h=0, col="black", lty=2, lwd=1)

##Cluster5
par(mar=c(5,5,1,1),bg="white")
plot(allclusters[5,c(1:9)],col="blue", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=c((1.5*min(na.omit(allclusters[5,]))),(1.5*max(na.omit(allclusters[5,])))),lwd=2)
title(main="Cluster 5")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(2.5),labels=c("Control"),padj=3,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=3,tick=F,col.axis="black")
lines(allclusters[5,c(11:19)],type="l",lty=2,col="blue",lwd=2)
##legend("topleft", legend = c("OreR;OreR","sm21;OreR"), col = c("blue","blue"), lty= c(1,2), lwd = 1, cex=.7)
abline(h=0, col="black", lty=2, lwd=1)





counts=read.csv("SingleCountReads.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)

genesinkegg=function(gois,keggs){
  
  convert=read.csv(("FBgnConversionTable.csv"),row.names=1)
  temp=convert
  temp$inKEGG=0
  kegg2=as.list(names(kegg))
  names(kegg2)=unlist(kegg)
  inkegg=as.vector(unlist(kegg2[names(kegg2)==paste("path:dme",keggs,sep="")]))
  info=(temp[inkegg,])
  i=1
  while(i<=length(inkegg)){
    temp$inKEGG[temp$Symbol==inkegg[i]]=1
    i=i+1
  }
  kgenes=temp[gois,]
  kgenes=kgenes[kgenes$inKEGG==1,]
  return(kgenes)
}


gs=intersect(row.names(ORSR), row.names(logFCdata[logFCdata[,17]==1,]))
gs=as.character(convert[gs,4])
##motif enrichment
setwd("/Users/johnsantiago/Downloads/")
data(motifAnnotations_dmel_v8)
summary(gs)
motifRankings <- importRankings("dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
setwd(ywd)
motifEnrichmentTable_wGenes <- cisTarget(gs, motifRankings, motifAnnot=motifAnnotations_dmel_v8)

motifs_AUC <- calcAUC(gs, motifRankings, nCores=1)
hist(motifs_AUC["geneSet",c(1:100)], main="hypoxia", xlab="AUC histogram",
     breaks=100, col="#ff000050", border="darkred")
nes3 <- (3*sd(motifs_AUC[1,])) + mean(motifs_AUC[1],])
abline(v=nes3, col="red")


enrtf=as.data.frame(motifEnrichmentTable_wGenes)
head(enrtf)
motifEnrichmentTable=motifEnrichmentTable_wGenes

signifMotifNames <- motifEnrichmentTable$motif[1:3]

incidenceMatrix <- getSignificantGenes(gs, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix



library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)


setwd("/Users/johnsantiago/Desktop/TF Analysis/")
write.csv(enrtf,"Cluster5_TFenrichment.csv")
tfbreakdown=function(tfvar){
  tfvar=tfvar[tfvar$TF_highConf!="",]
  tfvar$TF_highConf=as.character(tfvar$TF_highConf)
  i=1
  while(i<=nrow(tfvar)){
    tfvar$TF_highConf[i]=strsplit(as.character(tfvar$TF_highConf[i])," \\(")[[1]][1]
    i=i+1
  }
  
  temp=strsplit(as.character(tfvar$TF_highConf),"; ")
  i=1
  temp2=0
  while(i<=length(temp)){
    temp2=c(temp2,rep(i,length(temp[[i]])))
    i=i+1
  }
  tfvar=tfvar[temp2,]
  tfvar$TF_highConf=unlist(temp)
  
  temp=as.character(tfvar$enrichedGenes)
  temp=strsplit(temp,";")
  names(temp)=tfvar$TF_highConf
  namestf=1
  i=1
  while(i<=length(temp)){
    namestf=c(namestf,rep(names(temp)[i],length(temp[[i]])))
    i=i+1
  }
  namestf=namestf[-1]
  
  temp2=matrix(0,nrow=length(unlist(temp)),ncol=2)
  colnames(temp2)=c("Transcription Factor","Enriched Gene")
  temp2[,1]=namestf
  temp2[,2]=unlist(temp)
  return(temp2)
}


##Cluster 5 TF and KEGG Network Diagram
tfvar=read.csv("/Users/johnsantiago/Desktop/TF Analysis/Cluster5_TFenrichment.csv")
tftable5=tfbreakdown(tfvar)
tftable5=tftable5[,c(1,2,2)]
tftable5[,3]=0
i=1
while(i<=nrow(tftable5)){
  tftable5[i,3]=row.names(convert[convert$Symbol==tftable5[i,2],])
  i=i+1
}


gois=row.names(ORSR)
gois=setdiff(gois,row.names(convert[convert$mt==1,]))
gois=intersect(row.names(logFCdata)[logFCdata[,17]=="5"],gois)
keggs=c("03010","00190","03050","00071","00280","04146","00061","00640","04140")
i=1
setwd(ywd)
while(i<10){
  kgenes=genesinkegg(gois,keggs[i])
  kgenes$KEGG.Path=keggs[i]
  if(i==1){
    sigkgenes=kgenes
  }
  if(i>1){
    sigkgenes=rbind(sigkgenes,kgenes)
  }
  i=i+1
}
kgenes=sigkgenes[as.character(unique(sigkgenes[,1])),]

i=1
while(i<=nrow(kgenes)){
  tftable5[tftable5[,2]==kgenes[i,4],3]=1
  i=i+1
}

alltf=unique(tftable5[,1])
names(alltf)=alltf
i=1
while(i<=length(alltf)){
  alltf[i]=length(unique(tftable5[tftable5[,1]==alltf[i],3]))
  i=i+1
}

tftable=sigkgenes[,c(11,4)]
tftable$category="KEGG"
tftable5[,3]="TF"
tftable5=tftable5[tftable5[,1]=="Dref"|tftable5[,1]=="gt",]
colnames(tftable5)=colnames(tftable)
i=1
while(i<=length(unique(tftable[,2]))){
  tftable=rbind(tftable,tftable5[tftable5[,2]==(unique(tftable[,2])[i]),])
  i=i+1
}


im=matrix(0,nrow = length(unique(tftable[,1])),ncol= length(unique(tftable[,2])))
row.names(im)=unique(tftable[,1])
colnames(im)=unique(tftable[,2])
i=1
while(i<=nrow(im)){
  im[i,as.character(unique(tftable[tftable[,1]==row.names(im)[i],2]))]=1
  i=i+1
}

library(reshape2)
edges=melt(im)

edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
edgecolors=genes
edgecolors[c(1:length(genes))]="gray"
reds=unique(tftable5[tftable5[,1]=="Dref",2])
i=1
while(i<=length(reds)){
  edgecolors[genes==reds[i]]="red"
  i=i+1
}
blues=unique(tftable5[tftable5[,1]=="gt",2])
i=1
while(i<=length(blues)){
  edgecolors[genes==blues[i]]="blue"
  i=i+1
}
purples=intersect(reds,blues)
i=1
while(i<=length(purples)){
  edgecolors[genes==purples[i]]="purple"
  i=i+1
}
legendinfo=edgecolors
legendinfo[legendinfo=="red"]="Inferred Dref Target"
legendinfo[legendinfo=="blue"]="Inferred gt Target"
legendinfo[legendinfo=="purple"]="Inferred Dref and gt Target"
legendinfo[legendinfo=="gray"]="Not Dref or gt Target"
motifnames=c("Ribosome","Dref","giant","OXPHOS","Proteasome","FA Degradation","BCAA Degradation","Peroxisome","FA Biosynthesis","Propanoate Metabolism","Autophagy")
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifnames, rep("",length(genes))),    
                    title=c(rep("",length(motifs)), genes), 
                    shape=c(rep("star", length(motifs)), rep("dot", length(genes))),
                    color=c(rep("grey",11),edgecolors))
nodes$group=c("KEGG Pathway", "Dref","giant",rep("KEGG Pathway",8),legendinfo)

nodes <- data.frame(id=c(genes,motifs),   
                    label=c(rep("",length(genes)),motifnames),    
                    title=c(genes,motifnames), 
                    shape=c(rep("dot", length(genes)),rep("star", length(motifs))),
                    color=c(edgecolors,rep("grey",11)))
nodes$group=c(legendinfo,"KEGG Pathway", "Dref","giant",rep("KEGG Pathway",8))
##edges is 2 columns, c("from","to) which is c(hub, node connected to hub)
##nodes is 4 columns, c("id", "label", "titel", "groups"). 
  ##id corresponds to the connections made in edges. it is all unique hubs and associated nodes
  ##label is the always visible wording on the node
  ##title is the invisible wording that shows up for a node when its hovered over
  ##group is for groups assignment and legend


##takes a "tftable" and turns it into an edge table
im=matrix(0,nrow = length(unique(tftable[,1])),ncol= length(unique(tftable[,2])))
row.names(im)=unique(tftable[,1])
colnames(im)=unique(tftable[,2])
i=1
while(i<=nrow(im)){
  im[i,as.character(unique(tftable[tftable[,1]==row.names(im)[i],2]))]=1
  i=i+1
}
edges=melt(im)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")


##edge table is just the hub in 1 column and the associated nodes in the other
##example "OXPHOS" fill column 1 and all the OXPHOS genes in column 2
##can just paste these together
##edge data columns c("from","to")
##node data columns c("id","label","title","group")
##plotdata needs a group for each element in the legend
##plotdata columns c("group.label","fill.color","border.color","physics")

##makes edges
##Cluster 5 TF data
tfvar=read.csv("/Users/johncsantiago/Documents/Old Mac Files/Old Mac/TF Analysis/Cluster5_TFenrichment.csv")
tftable5=tfbreakdown(tfvar)



##gather fbgn for gene symbols used in TF analysis
symconvert=row.names(convert)
names(symconvert)=convert[,4]
orsrtf=unique(symconvert[tftable5[,2]])

temp=tftable5[,1]
names(temp)=tftable5[,2]
temp2=temp[unique(Dref[,2])]

gs=intersect(row.names(ORSR), row.names(logFCdata[logFCdata[,17]==5,]))
##gs=as.character(convert[gs,4])
orsrtf=gs

##Cluster 5 KEGG data that was used in TF analysis
OXPHOS=as.matrix(genesinkegg(orsrtf,"00190")[,c(1,4)])
OXPHOS[,1]="OXPHOS"
Proteasome=as.matrix(genesinkegg(orsrtf,"03050")[,c(1,4)])
Proteasome[,1]="Proteasome"
Ribosome=as.matrix(genesinkegg(orsrtf,"03010")[,c(1,4)])
Ribosome[,1]="Ribosome"
Propanoate.Metabolism=as.matrix(genesinkegg(orsrtf,"00640")[,c(1,4)])
Propanoate.Metabolism[,1]="Propanoate Metabolism"
BCAA.Degradation=as.matrix(genesinkegg(orsrtf,"00280")[,c(1,4)])
BCAA.Degradation[,1]="BCAA Degradation"
FA.Biosynthesis=as.matrix(genesinkegg(orsrtf,"00061")[,c(1,4)])
FA.Biosynthesis[,1]="FA Biosynthesis"
FA.Degradation=as.matrix(genesinkegg(orsrtf,"00071")[,c(1,4)])
FA.Degradation[,1]="FA Degradation"
Autophagy=as.matrix(genesinkegg(orsrtf,"04140")[,c(1,4)])
Autophagy[,1]="Autophagy"
Peroxisome=as.matrix(genesinkegg(orsrtf,"04146")[,c(1,4)])
Peroxisome[,1]="Peroxisome"  

##genes regulated by tf that will be nodes
Dref=as.matrix(unique(tftable5[tftable5[,1]=="Dref",c(1,2)]))
giant=as.matrix(unique(tftable5[tftable5[,1]=="gt",c(1,2)]))
giant[,1]="giant"
row.names(Dref)=Dref[,2]
Dref=Dref[intersect(Dref[,2],keggedges[,2]),]
row.names(giant)=giant[,2]
giant=giant[intersect(giant[,2],keggedges[,2]),]


hubs=list(Dref,giant,OXPHOS,Proteasome,Ribosome,Propanoate.Metabolism,BCAA.Degradation,FA.Biosynthesis,FA.Degradation,Autophagy,Peroxisome)

keggedges=edges=as.data.frame(do.call(rbind,hubs[c(3:11)]))
tfedges=edges=as.data.frame(do.call(rbind,hubs[c(1,2)]))

edges=as.data.frame(do.call(rbind,hubs))
colnames(edges)=c("from","to")

##specific for thesis presentation
edges=read.csv("/Users/johnsantiago/Desktop/Cluster5_KEGG_Edges.csv")
colnames(edges)=c("from","to")
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(genes,motifs),   
                    label=c(genes,motifs),    
                    title=c(genes,motifs),
                    group=c(genes,motifs))

visnet.cluster5=visNetwork(nodes, edges) %>%
  visNodes(font=list(color="black",size=80),mass=3.5) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE)


##hub node and group data
motifs <- unique(as.character(edges[,1]))
  motifs.label=motifs
  motifs.title=rep("",length(motifs))
  motifs.group=c("Dref","giant",rep("KEGG Pathway",9))
 head motifs.group.order=motifs.group[c(1,2,3)]
  motifs.fill.color=c("#03BDF6","#FB699A","#FEF861")
  motifs.border.color=c("#0350F6","#DB0003","#E56E00")
  motifs.physics=c("FALSE","FALSE","TRUE")

##gene node and group data    
genes <- unique(as.character(edges[,2]))
  genes.label=rep("",length(genes))
  genes.title=genes
  genes.group.names=c("Dref Target Gene","giant Target Gene","Dref and giant Target Gene","Not Dref or giant Target")
  genes.fill.color=c("#0350F6","#F8253D","#A73DFF","#FFC942")
  genes.border.color=c("black","black","black","black")
  genes.group=c(rep(genes.group.names[4],length(genes)))
  names(genes.group)=genes
  genes.group[Dref[,2]]=genes.group.names[1]
  genes.group[giant[,2]]=genes.group.names[2]
  genes.group[intersect(Dref[,2],giant[,2])]=genes.group.names[3]


##Assign info to node table  
nodes <- data.frame(id=c(genes,motifs),   
                    label=c(genes.label,motifs.label),    
                    title=c(genes.title,motifs.title),
                    group=c(genes.group,motifs.group))

##Plot the network
##Start the plot
visnet.cluster5=visNetwork(nodes, edges) %>%
  visNodes(font=list(color="black",size=80),mass=1.5) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
  visPhysics(solver="forceAtlas2Based",enabled=TRUE,timestep=.2)



##Hub Groups
i=1
while(i<=length(motifs.group.order)){
  ##visnet.cluster5=visGroups(visnet.cluster5, groupname = motifs.group.order[i], shape = "box", color = list(background = motifs.fill.color[i], border=motifs.border.color[i]),physics=motifs.physics[i])
  visnet.cluster5=visGroups(visnet.cluster5, groupname = motifs.group.order[i], shape = "box", color = list(background = motifs.fill.color[i], border=motifs.border.color[i]),physics=TRUE)
  i=i+1
}

##Target Node Groups
i=1
while(i<=length(genes.group.names)){
  visnet.cluster5=visGroups(visnet.cluster5, groupname = genes.group.names[i], shape = "dot", color = list(background = genes.fill.color[i], border=genes.border.color[i]))
  i=i+1
}

##Load network
visLegend(visnet.cluster5, main="Legend", position="right", ncol=1,width=.1) 

##Cluster 1 TF and KEGG Network Diagram
tfvar=read.csv("/Users/johnsantiago/Desktop/TF Analysis/Cluster1_TFenrichment.csv")
tftable1=tfbreakdown(tfvar)
tftable1=tftable1[,c(1,2,2)]
tftable1[,3]=0
gois=row.names(ORSR)
gois=setdiff(gois,row.names(convert[convert$mt==1,]))
gois=intersect(row.names(logFCdata)[logFCdata[,17]=="1"],gois)
keggs=c("00010","00190","00620","00230","04341","03050","00020","00052","04145")
i=1
setwd(ywd)
while(i<10){
  kgenes=genesinkegg(gois,keggs[i])
  kgenes$KEGG.Path=keggs[i]
  if(i==1){
    sigkgenes=kgenes
  }
  if(i>1){
    sigkgenes=rbind(sigkgenes,kgenes)
  }
  i=i+1
}
tftable=sigkgenes[,c(11,4)]
tftable$category="KEGG"
tftable1[,3]="TF"
tftable1=tftable1[tftable1[,1]=="Abd-B",]
colnames(tftable1)=colnames(tftable)
i=1
while(i<=length(unique(tftable[,2]))){
  tftable=rbind(tftable,tftable1[tftable1[,2]==(unique(tftable[,2])[i]),])
  i=i+1
}




im=matrix(0,nrow = length(unique(tftable[,1])),ncol= length(unique(tftable[,2])))
row.names(im)=unique(tftable[,1])
colnames(im)=unique(tftable[,2])
i=1
while(i<=nrow(im)){
  im[i,as.character(unique(tftable[tftable[,1]==row.names(im)[i],2]))]=1
  i=i+1
}
edges=melt(im)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
edgecolors=genes
edgecolors[c(1:length(genes))]="gray"
reds=unique(tftable1[tftable1[,1]=="Abd-B",2])
i=1
while(i<=length(reds)){
  edgecolors[genes==reds[i]]="red"
  i=i+1
}
abdorno=edgecolors
abdorno[abdorno=="red"]="Inferred Abd-B Target"
abdorno[abdorno=="gray"]="Not Abd-B Target"
motifnames=c("Glycolysis","Pyruvate Metabolism","TCA","Galactose Metabolism","Purine Metabolism","Abd-B","OXPHOS","Phagosome","Hedgehog Signaling","Proteasome")
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifnames, rep("",length(genes))),    
                    title=c(rep("",length(motifs)), genes), 
                    shape=c(rep("star", length(motifs)), rep("dot", length(genes))),
                    color=c(rep("grey",10),edgecolors))
nodes$group=c(rep("KEGG Pathway",5),"Abd-B",rep("KEGG Pathway",4),abdorno)

visnet.cluster1=visNetwork(nodes[,c(1:3,6)], edges[,c(1:2)]) %>%
  visNodes(font=list(color="black",size=60),mass=3.5) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE)
visnet.cluster1=visGroups(visnet.cluster1, groupname = "KEGG Pathway", shape = "box", color = list(background = "#FEF861", border="#E56E00"),size=70,physics=TRUE)

visnet.cluster1=visGroups(visnet.cluster1, groupname = "Abd-B", shape = "box", color = list(background = '#FB699A', border='#DB0003'),size=80,physics=FALSE,mass=10)
visnet.cluster1=visGroups(visnet.cluster1, groupname = "Inferred Abd-B Target", shape = "dot", color = list(background = '#F8253D', border="black"))
visnet.cluster1=visGroups(visnet.cluster1, groupname = "Not Abd-B Target", shape = "dot", color = list(background = '#FFC942', border="black"))
visLegend(visnet.cluster1, main="Legend", position="right", ncol=1,width=.1) 





##Cluster 1 and Cluster 5 OXPHOS Complex Distribution
oxphossigs=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/OxphosUnitLabeledORSR.csv",row.names=1)
oxphostable=(oxphossigs[,c(4,6,13)])
oxphostable[,2]=paste("C",oxphostable[,2],sep="")
i=1
while(i<=nrow(oxphostable)){
  oxphostable[i,3]=convert[row.names(oxphostable)[i],6]
  i=i+1
}
oxphostable=oxphostable[oxphostable[,3]!="2",]
oxphostable[oxphostable[,2]=="Cv",2]="C5"


im=matrix(0,nrow = length(c(unique(oxphostable[,2]),unique(oxphostable[,3]))),ncol= length(unique(oxphostable[,1])))
row.names(im)=c(unique(oxphostable[,2]),unique(oxphostable[,3]))
colnames(im)=unique(oxphostable[,1])
i=1
while(i<=length(unique(oxphostable[,2]))){
  im[i,as.character(unique(oxphostable[oxphostable[,2]==row.names(im)[i],1]))]=1
  i=i+1
}
while(i<=nrow(im)){
  im[i,as.character(unique(oxphostable[oxphostable[,3]==row.names(im)[i],1]))]=1
  i=i+1
}

edges=melt(im)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
edgecolors=genes
edgecolors[c(1:length(genes))]="gray"
reds=unique(oxphostable[oxphostable[,3]=="1",1])
i=1
while(i<=length(reds)){
  edgecolors[genes==reds[i]]="red"
  i=i+1
}
blues=unique(oxphostable[oxphostable[,3]=="5",1])
i=1
while(i<=length(blues)){
  edgecolors[genes==blues[i]]="blue"
  i=i+1
}

genecat=edgecolors
genecat[genecat=="red"]="Cluster 1"
genecat[genecat=="blue"]="Cluster 5"
genecat[genecat=="gray"]="Not Cluster 1 or Cluster 5"
motifnames=c("Complex II","Cluster 1","Complex III","Complex IV","Complex V","Complex I","Cluster 5")
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifnames, rep("",length(genes))),    
                    title=c(rep("",length(motifs)), genes), 
                    shape=c(rep("star", length(motifs)), rep("dot", length(genes))),
                    color=c(rep("grey",7),edgecolors))
nodes$group=c("OXPHOS Complex","Co-Expression Cluster",rep("OXPHOS Complex",4),"Co-Expression Cluster",genecat)
visnet.cluster1=visNetwork(nodes[,c(1:3,6)], edges[,c(1:2)]) %>%
  visNodes(font=list(color="black",size=40),size=10) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE)
visnet.cluster1=visGroups(visnet.cluster1, groupname = "OXPHOS Complex", shape = "star", color = list(background = "yellow", border="black"),size=30,physics=FALSE)
visnet.cluster1=visGroups(visnet.cluster1, groupname = "Co-Expression Cluster", shape = "star", color = list(background = "red", border="black"),size=30,physics=FALSE)
visnet.cluster1=visGroups(visnet.cluster1, groupname = "Cluster 1", shape = "dot", color = list(background = "blue", border="black"))
visnet.cluster1=visGroups(visnet.cluster1, groupname = "Cluster 5", shape = "square", color = list(background = "green", border="black"))
visLegend(visnet.cluster1, main="Legend", position="right", ncol=1,width=.1) 


##motif table
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
resultsSubset <- motifEnrichmentTable_wGenes_wLogo

##network diagram
library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

signifMotifNames <- motifEnrichmentTable_wGenes$motif

incidenceMatrix <- getSignificantGenes(gs, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

im=matrix(0,nrow = length(unique(tftable[,1])),ncol= length(unique(tftable[,2])))
row.names(im)=unique(tftable[,1])
colnames(im)=unique(tftable[,2])
dupedges=intersect(row.names(im),colnames(im))
i=1
while(i<=length(dupedges)){
  row.names(im)[row.names(im)==dupedges[i]]=paste(dupedges[i],"(TF)",sep=" ")
i=i+1
}



i=1
while(i<=nrow(im)){
  im[i,unique(tftable[tftable[,1]==row.names(im)[i],2])]=1
  i=i+1
}

imdref=im[c(3,6,6),]
imdref[3,]=0
imnotdref=im[-c(3,6),]
i=1
while(i<=nrow(imnotdref)){
  imdref[3,]=imdref[3,]+imnotdref[i,]
  i=i+1
}
row.names(imdref)[3]="other"
imdref[3,imdref[3,]>0]=1
edges=melt(imdref)
edges=melt(im[c(1:10),])
##nodes <- read.csv("/Users/johnsantiago/Downloads/sunbelt2019/Data files/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
##links <- read.csv("/Users/johnsantiago/Downloads/sunbelt2019/Data files/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
##nodes$shape="diamond"
##nodes$color="red"
##nodes$background="red"
##nodes$border="black"

edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c("gray","red","green", rep("black", length(genes))))
visNetwork(nodes[,c(1,5)], edges[,c(1:2)]) %>%
  visPhysics(stabilization = FALSE)

    visNodes(color= list(background="red",border="black")) %>%
  visEdges(color="black")

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
colors=as.matrix(edges)
colors[,2]="black"
colors[colors[,1]=="Dref",]=("blue")
colors[colors[,1]=="gt",]=("green")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
colors=rep("black",length(motifs))
colors[motifs=="Dref"]=("blue")
colors[motifs=="gt"]=("green")
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(colors, rep("red", length(genes))))
##nodes <- data.frame(id=c(motifs, genes),   
##                    label=c(motifs, genes),    
##                    title=c(motifs, genes), # tooltip 
##                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
##                    color=c(rep("black", length(motifs)), rep("red", length(genes))))

visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)


tffbgn=as.character(tfsum$TF_highConf)
i=1
while(i<=length(tffbgn)){
  if(substring(tffbgn[i],1,1)==" "){
    tffbgn[i]=substring(tffbgn[i],2,nchar(tffbgn[i]))
  }
  i=i+1
}
xx <- as.list(org.Dm.egSYMBOL2EG)

i=1
while(i<=length(tffbgn)){
  ##looks for the gene symbol for each entrez id then replaces the entrez id
  if(length(xx[[tffbgn[i]]])>0){
    tffbgn[i]=(xx[[tffbgn[i]]])
  }
  i=i+1
}

tfbreakdown=function(tfvar){
  tfvar=tfvar[tfvar$TF_highConf!="",]
  tfvar$TF_highConf=as.character(tfvar$TF_highConf)
  i=1
  while(i<=nrow(tfvar)){
    tfvar$TF_highConf[i]=strsplit(as.character(tfvar$TF_highConf[i])," \\(")[[1]][1]
    i=i+1
  }
  
  temp=strsplit(as.character(tfvar$TF_highConf),"; ")
  i=1
  temp2=0
  while(i<=length(temp)){
    temp2=c(temp2,rep(i,length(temp[[i]])))
    i=i+1
  }
  tfvar=tfvar[temp2,]
  tfvar$TF_highConf=unlist(temp)
  
  temp=as.character(tfvar$enrichedGenes)
  temp=strsplit(temp,";")
  names(temp)=tfvar$TF_highConf
  namestf=1
  i=1
  while(i<=length(temp)){
    namestf=c(namestf,rep(names(temp)[i],length(temp[[i]])))
    i=i+1
  }
  namestf=namestf[-1]
  
  temp2=matrix(0,nrow=length(unlist(temp)),ncol=2)
  colnames(temp2)=c("Transcription Factor","Enriched Gene")
  temp2[,1]=namestf
  temp2[,2]=unlist(temp)
  return(temp2)
}
tfvar=read.csv("/Users/johnsantiago/Desktop/TF Analysis/Cluster5_TFenrichment.csv")
tftable=tfbreakdown(tfvar)




tfsum$tffentrez=tffbgn
i="04350"
gik=as.list(org.Dm.egPATH2EG)
tfegs=matrix(unique(tfsum$tffentrez),ncol=1)
row.names(tfegs)=tfegs[,1]
blank=matrix("",nrow=nrow(tfegs),ncol=1)
tfegs=cbind(tfegs,blank)
i=1
j=2
while(i<=nrow(tfegs)){
  if(length(intersect(row.names(tfegs),gik[[i]]))>0){
    temp=intersect(row.names(tfegs),gik[[i]])
    tfegs[temp,j]=names(gik[i])
    tfegs=cbind(tfegs,blank)
    j=j+1
  }
  i=i+1
}

intersect(row.names(tfegs),gik[[i]])
i=i-1

# Save results as a data.frame
enrtf=as.data.frame(motifEnrichmentTable_wGenes)
head(enrtf)
setwd("/Users/johnsantiago/Desktop/TF Analysis/")
write.csv(enrtf,"OXPHOS_Cluster2_TranscriptionFactorEnrichement.csv")

##TF analysis
tfsum=read.csv("OXPHOS_Cluster2_TranscriptionFactorEnrichement.csv",row.names=1)
tfdesc=tfsum[tfsum$Description!="",]
##tf1=tfsum[tfsum[,1]==1,]


temp=unique(tftable[,1])
tftable=tftable[tftable[,1]=="ERR",]

temp2=convert[convert$Symbol==tftable[1,2],]
i=2
while(i<=length(temp)){
  temp3=convert[convert$Symbol==tftable[i,2],]
  temp2=rbind(temp2,temp3)
  i=i+1
}
write.csv(temp2,"temp.csv")





##FigureS3
cats=c("00190","03050","00010","00020","00280","00071","00030","00052","00640","00500","04145","00061","00620","00630","04146","00051","04141","04137","00750","00520","04213")
LETTERS[4]
catcount=1
gois=row.names(ORSR)
gois=setdiff(gois,row.names(convert[convert$mt==1,]))
keggs=cats[catcount]
hmtitle=paste("KEGGid",keggs,sep="")
pngname=paste("FigureS3",LETTERS[catcount],"_KEGGid",keggs,".png",sep="")
temp=convert[convert$mt==1,]
gois=setdiff(gois,row.names(temp))
kgenes=genesinkegg(gois,keggs)
heatmapdata=keggheatmapinfo(kgenes)
col.pan <- colorpanel(100, "black","white","red")

setwd("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/")
png(filename=pngname,pointsize=12,width=6000,height=6000,res=400)
##png(filename="OXPHOS.png",pointsize=12,width=5000,height=2000,res=400)
heatmapped=heatmap.2(as.matrix(t(heatmapdata[[2]])),Rowv=F,Colv=F,trace="none",col=col.pan,scale="col",key=TRUE,cexRow = .6,cexCol = .6,main=hmtitle)
dev.off()
setwd(ywd)

,RowSideColors = heatmapdata[[3]]

catcount=catcount+1
catcount


##TableS2 code
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")
singlecounts=singlecounts[,trim]
groups=groups[trim,]
groups=as.matrix(groups)
groups=as.data.frame(groups[c(1:28),])
x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 
design<-model.matrix(~0+groups$Group)
colnames(design) <- levels(groups$Group)
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)
compare = makeContrasts((SO4HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
G_X_E<-topTags(lrt, n=14003,adjust.method="BH", sort.by="PValue")
results=G_X_E$table
write.csv(results,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS2L_SO4HR-SO0HC.csv")

##TableS3
write.csv(oocontrol,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3A_OOC.csv")
write.csv(oorapa,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3B_OOR.csv")
write.csv(socontrol,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3C_SOC.csv")
write.csv(sorapa,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3D_SOR.csv")

write.csv(oocc,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3E_OOC-OOR.csv")
write.csv(socc,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3F_SOC-SOR.csv")

write.csv(oocsoc,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3G_OOC-SOC.csv")
write.csv(oorsor,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS3H_OOR-SOR.csv")


##TableS4-11
###################
#### GOseq analysis
###################
##gets the KEGG enrichment values from a geneset and writes to a pdf (also does GO enrichment but not written to a pdf)
##run 1 time before using goseq to generate the kegg correlation list
kegg=keggLink("pathway","dme")
temp=names(kegg)
temp2=convert[convert[,5]!="",]
i=1
while(i<=nrow(temp2)){
  if(length(temp[temp==temp2[i,5]])){
    temp[temp==temp2[i,5]]=as.character(temp2[i,4])
  }
  i=i+1
}
names(kegg)=temp
kegg=as.list(kegg[substr(names(kegg),1,3)!="dme"])

##sig genes
setwd(ywd)
##significant genes. FBgn IDs
genelist=row.names(ORSR)
genelist=setdiff(genelist,row.names(convert[convert$mt==1,]))
genelist=intersect(genelist,row.names(logFCdata)[logFCdata[,17]==5])
###the output file = csvname_KEGGEnrichmentAnalysis.csv
gocsvname="TableS10"
keggcsvname="TableS11"
##backgroud genes to be used. an impulseDE output file
data=oocontrol
data=(data[setdiff(row.names(data),row.names(convert[convert$mt==1,])),])
data=data[,c(1,2)]
data[,1]=as.character(convert[row.names(data),4])
data[data[,1]=="jus",][1,1]="jus2"
data[,2]=0
data[genelist,2]=1
genes=as.integer(data[,2])
names(genes)=data[,1]
table(genes)
pwf=nullp(genes,"dm3","geneSymbol")
head(pwf)

##GO term enrichment is GO.wall
GO.wall=goseq(pwf,"dm3","geneSymbol")

##KEGG term enrichment is KEGG
KEGG=goseq(pwf,gene2cat=kegg)
KEGG$adjp=p.adjust(KEGG$over_represented_pvalue,method="BH")
row.names(KEGG)=substr(as.character(KEGG[,1]),9,13)
write.csv(GO.wall,paste("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/",gocsvname,".csv",sep=""))
write.csv(KEGG,paste("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/",keggcsvname,".csv",sep=""))

##TableS12
write.csv(logFCdata,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS12.csv")


##heatmapdata from cluster 1 kegg genes
##c1k=heatmapdata[[1]]

##heatmapdata from cluster 5 kegg genes
##c5k=heatmapdata[[1]]
c1kt=(log(c1k,base=10))
c5kt=(log(c5k,base=10))
##c1kt=c1k
##c5kt=c5k

c1kt=as.numeric(as.matrix(c1kt))
c5kt=as.numeric(as.matrix(c5kt))
c1kt=matrix(c1kt,ncol=1)
c5kt=matrix(c5kt,ncol=1)
Cluster1KEGG=data.frame(c1kt)
Cluster5KEGG=data.frame(c5kt)
names(Cluster1KEGG)="Cluster 1 KEGG Genes"
names(Cluster5KEGG)="Cluster 5 KEGG Genes"

boxplot(c(Cluster1KEGG,Cluster5KEGG),ylab="Expression Across All Conditions (log10)",col=c("firebrick","dodgerblue"))

t.test((1.6*c1kt),c5kt,paired=F)
Cluster1KEGG=Cluster1KEGG[,c(1,1)]
Cluster5KEGG=Cluster5KEGG[,c(1,1)]
colnames(Cluster1KEGG)=c("Cluster","Expression")
colnames(Cluster5KEGG)=c("Cluster","Expression")
canovadata=rbind(Cluster1KEGG,Cluster5KEGG)
canova=aov(Expression~Cluster,data=canovadata)

##uses heatmapdata for pca
kgenes=convert[intersect(row.names(ORSR),row.names(logFCdata)[logFCdata[,17]=="1"]),]
kgenes[,"inKEGG"]=1
heatmapdata=keggheatmapinfo(kgenes)

pca <- prcomp(t(heatmapdata[[1]]), scale.=TRUE) 
gr <- factor(colnames(heatmapdata[[1]]))
color=factor(colnames(heatmapdata[[1]]))
shape=factor(colnames(heatmapdata[[1]]))
pcacols=c(1,2,3,2,3,2,3,4,5,6,5,6,5,6)
par(xpd=T,mfrow=c(1,1))
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA,components = c(1,2),show.ellipses = F,title = metatype[i],radius = 1.5,legend=NULL)

eigs <- pca$sdev^2
ve=(eigs / sum(eigs))*100
names(ve)=c("PC1","PC2","PC3")
barplot(ve,ylim=c(0,100))

e1=(pca[[5]])[,1]

e1plotnames=names(e1)[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]
e1plotdata=c(1:16)
e1plot[1]=sum(e1[1],e1[2])/2
e1plot[2]=sum(e1[3],e1[4])/2
e1plot[3]=sum(e1[5],e1[6])/2
e1plot[4]=sum(e1[7],e1[8])/2
e1plot[5]=sum(e1[9],e1[10])/2
e1plot[6]=sum(e1[11],e1[12])/2
e1plot[7]=sum(e1[13],e1[14])/2

e1plot[8]=sum(e1[15],e1[16])/2
e1plot[9]=sum(e1[17],e1[18])/2
e1plot[10]=sum(e1[19],e1[20])/2
e1plot[11]=sum(e1[21],e1[22])/2
e1plot[12]=sum(e1[23],e1[24])/2
e1plot[13]=sum(e1[25],e1[26])/2
e1plot[14]=sum(e1[27],e1[28])/2

e1plot[15]=sum(e1[29],e1[30])/2
e1plot[16]=sum(e1[31],e1[32])/2

names(e1plot)=e1plotnames
dev.off()
plot(as.numeric(e1plot[1:4]),type="l",col="blue",ylim=c(-50,50))
lines(as.numeric(e1plot[5:8]),col="red")
lines(as.numeric(e1plot[9:12]),col="blue",lty=2)
lines(as.numeric(e1plot[13:16]),col="red",lty=2)
abline(h=0,lty=2)

##cluster 5 plot using line plot data from kgenes
sdev1=c(sqrt(var(plotdata[,1])),sqrt(var(plotdata[,2])),sqrt(var(plotdata[,3])),sqrt(var(plotdata[,4])))
sdev2=c(sqrt(var(plotdata[,6])),sqrt(var(plotdata[,7])),sqrt(var(plotdata[,8])),sqrt(var(plotdata[,9])))
sdev3=c(sqrt(var(plotdata[,11])),sqrt(var(plotdata[,12])),sqrt(var(plotdata[,13])),sqrt(var(plotdata[,14])))
sdev4=c(sqrt(var(plotdata[,16])),sqrt(var(plotdata[,17])),sqrt(var(plotdata[,18])),sqrt(var(plotdata[,19])))

plot(line1,col="blue", ylim=c(-2,2.5),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,9),xlab="",type="line",lwd=2)
lines(x=c(6:9),line2[c(1:4)],type="l",col="blue",lwd=2)
lines(x=c(1:4),line3,type="l",col="red",lwd=2)
lines(x=c(6:9),line4,type="l",col="red",lwd=2)
abline(h=0, col="black", lty=2, lwd=1)

x=c(1:4)
arrows(x, line1-sdev1, x, line1+sdev1, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, line3-sdev3, x, line3+sdev3, length=0.05, angle=90, code=3,col="red",lwd=2)
x=c(6:9)
arrows(x, line2-sdev2, x, line2+sdev2, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, line4-sdev4, x, line4+sdev4, length=0.05, angle=90, code=3,col="red",lwd=2)


title(main="Cluster 5 Genes")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")

axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="black")
legend("topleft",legend = c("OreR;OreR","sm21:OreR"),fill = c("blue","red"),cex = .5)

##cluster 1 plot using line plot data from kgenes
sdev1=c(sqrt(var(plotdata[,1])),sqrt(var(plotdata[,2])),sqrt(var(plotdata[,3])),sqrt(var(plotdata[,4])))
sdev2=c(sqrt(var(plotdata[,6])),sqrt(var(plotdata[,7])),sqrt(var(plotdata[,8])),sqrt(var(plotdata[,9])))
sdev3=c(sqrt(var(plotdata[,11])),sqrt(var(plotdata[,12])),sqrt(var(plotdata[,13])),sqrt(var(plotdata[,14])))
sdev4=c(sqrt(var(plotdata[,16])),sqrt(var(plotdata[,17])),sqrt(var(plotdata[,18])),sqrt(var(plotdata[,19])))

plot(line1,col="blue", ylim=c(-2,2.5),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,9),xlab="",type="line",lwd=2)
lines(x=c(6:9),line2[c(1:4)],type="l",col="blue",lwd=2)
lines(x=c(1:4),line3,type="l",col="red",lwd=2)
lines(x=c(6:9),line4,type="l",col="red",lwd=2)
abline(h=0, col="black", lty=2, lwd=1)

x=c(1:4)
arrows(x, line1-sdev1, x, line1+sdev1, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, line3-sdev3, x, line3+sdev3, length=0.05, angle=90, code=3,col="red",lwd=2)
x=c(6:9)
arrows(x, line2-sdev2, x, line2+sdev2, length=0.05, angle=90, code=3,col="blue",lwd=2)
arrows(x, line4-sdev4, x, line4+sdev4, length=0.05, angle=90, code=3,col="red",lwd=2)


title(main="Cluster 1 KEGG Genes")
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")

axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="black")
legend("topleft",legend = c("OreR;OreR","sm21:OreR"),fill = c("blue","red"),cex = .5)



##taking another look at the noise between libraries
##Generating filter vector called trim that removes the most outlying replicate from each treatment group
##Using distance between libraries in an MDS plot to determine which library 
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")

##uses the edger function that is defined below (needs to be ran before using) to trim out the genes with low read counts
x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 

##generates a table showing distance between replicates in MDS space
##MDSdata=z$counts
MDSgroups=groups
lcpm=cpm(z,log=TRUE)
##lcpm=testesremoved
colors=as.factor(groups$Group)
label=groups$Group
levels(colors)=c(1:length(levels(colors)))
temp=plotMDS(lcpm, col=as.vector(colors),labels=label,gene.selection="pairwise",top=nrow(lcpm))
similarity=temp$distance.matrix
samples=unique(groups$Group)
similaritytable=groups[,c(5,4)]
similaritytable$rep1=c(rep(c(1,1,2),14))
similaritytable$rep2=c(rep(c(2,3,3),14))
similaritytable$Replicate=c(rep(c("1~2","1~3","2~3"),14))
similaritytable$distance=0
i=1
j=1
while(i<=length(samples)){
  coords=groups$Group==samples[i]
  temp=similarity[coords,coords]
  similaritytable$distance[j]=temp[2,1]
  j=j+1
  similaritytable$distance[j]=temp[3,1]
  j=j+1
  similaritytable$distance[j]=temp[3,2]
  j=j+1
  i=i+1
}

similaritytable
plot(similaritytable[,5])

library(scatterplot3d)
library(plot3D)
library("plot3Drgl")
x=1
y=2
z=3

scatter3D(cpm2[,x],cpm2[,y],cpm2[,z],cex=.3,xlim=c(0,4000),ylim=c(0,4000),zlim=c(0,4000),xlab=groups[x,6],ylab=groups[y,6],zlab=groups[z,6])
segments3D(0,0,0,4000,4000,4000,add=TRUE)
plotrgl()
rglwidget(width=500,height=500)




lcpm=cpm(z,log=T)
pca <- prcomp(t(lcpm[,c(1:3)]), scale.=TRUE) 
gr <- factor(colnames(lcpm)[c(1:3)])
pcacols=c("black","red","blue")
##pca3d(pca, group=gr,show.group.labels=T)
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA)

eigs <- pca$sdev^2
ve=(eigs / sum(eigs))*100
names(ve)=c("PC1","PC2","PC3")
barplot(ve,ylim=c(0,100))

tissuedata=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/FlyAtlas2Data.csv",row.names=1)
tissueid=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/FlyAtlasTissueIDKey.csv",row.names=2)
##tissueid=tissueid[,c(2:4)]
tissue1=row.names(lcpm)
##generates a table (exampleid) with tissue specific FPKM for each gene in the tissue1 list
##Also generates a table (enriched) which is the FPKM normalized to the whole fly FPKM. Method is used by FlyAtlas2
exampleid=as.matrix(tissueid)
enriched=exampleid
i=1
while(i<=length(tissue1)){
  gene=tissue1[i]
  example=tissuedata[tissuedata[,1]==gene,]
  
  exampleid=cbind(exampleid,0)
  colnames(exampleid)[length(colnames(exampleid))]=gene
  exampleid[as.character(example$TissueID),gene]=as.vector(example$FPKM)
  exampleid[,gene]=as.numeric(exampleid[,gene])
  
  enriched=cbind(enriched,0)
  colnames(enriched)[length(colnames(enriched))]=gene
  enriched[as.character(example$TissueID),gene]=as.vector(example$FPKM)
  
  
  enriched[enriched[,2]=="'Female'",gene]=as.numeric(enriched[enriched[,2]=="'Female'",gene])/as.numeric(enriched["200",gene])
  enriched[enriched[,2]=="'Male'",gene]=as.numeric(enriched[enriched[,2]=="'Male'",gene])/as.numeric(enriched["100",gene])
  enriched[enriched[,2]=="'Larval'",gene]=as.numeric(enriched[enriched[,2]=="'Larval'",gene])/as.numeric(enriched["300",gene])
  enriched[enriched[,gene]=="Nan"|enriched[,gene]=="Inf",gene]=0
  enriched[,gene]=as.numeric(enriched[,gene])
  i=i+1
}

maleaovdata=enriched[enriched[,3]=="'Male'",]
testesspec=(t(maleaovdata[c(1:3),-c(1:9)]))
colnames(testesspec)=c("non-testes","testes","non-testes/testes")
i=1
while(i<=nrow(testesspec)){
  testesspec[i,1]=sum(as.numeric(maleaovdata[c(1:13,15),row.names(testesspec)[i]]))
  testesspec[i,2]=as.numeric(maleaovdata[14,row.names(testesspec)[i]])
  testesspec[i,3]=as.numeric(testesspec[i,1])/as.numeric(testesspec[i,2])
  i=i+1
}

par(mfrow=c(2,3))
i=0
plot(x=lcpm[,((i*3)+1)],y=lcpm[,((i*3)+2)],cex=.1,xlab=colnames(lcpm)[((i*3)+1)],ylab=colnames(lcpm)[((i*3)+2)])
plot(x=lcpm[,((i*3)+1)],y=lcpm[,((i*3)+3)],cex=.1,xlab=colnames(lcpm)[((i*3)+1)],ylab=colnames(lcpm)[((i*3)+3)])
plot(x=lcpm[,((i*3)+2)],y=lcpm[,((i*3)+3)],cex=.1,xlab=colnames(lcpm)[((i*3)+2)],ylab=colnames(lcpm)[((i*3)+3)])

plot(x=lcpm[row.names(testesspec)[testesspec[,3]>1],((i*3)+1)],y=lcpm[row.names(testesspec)[testesspec[,3]>1],((i*3)+2)],cex=.1,xlab=colnames(lcpm)[((i*3)+1)],ylab=colnames(lcpm)[((i*3)+2)])
plot(x=lcpm[row.names(testesspec)[testesspec[,3]>1],((i*3)+1)],y=lcpm[row.names(testesspec)[testesspec[,3]>1],((i*3)+3)],cex=.1,xlab=colnames(lcpm)[((i*3)+1)],ylab=colnames(lcpm)[((i*3)+2)])
plot(x=lcpm[row.names(testesspec)[testesspec[,3]>1],((i*3)+2)],y=lcpm[row.names(testesspec)[testesspec[,3]>1],((i*3)+3)],cex=.1,xlab=colnames(lcpm)[((i*3)+2)],ylab=colnames(lcpm)[((i*3)+3)])
i=i+1


testesremoved=lcpm[row.names(testesspec)[testesspec[,3]>1],]


##re-run edger with testes genes removed
x<-DGEList(singlecounts[row.names(testesspec)[testesspec[,3]>1],])
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 


##Showing mean values of OXPHOS cluster 1 and 5 genes
##Gene symbol from cluster 5
gene5="ATPsynG"

##Gene symbol from cluster 1
gene1="ATPsynGL"

barplot(as.matrix(cluster5means[gene5,]),main="Cluster 5")
barplot(as.matrix(cluster1means[gene1,]),main="Cluster1")

##KEGG enrichment analysis for cluster 1 and cluster 5
###################
#### GOseq analysis
###################
##gets the KEGG enrichment values from a geneset and writes to a pdf (also does GO enrichment but not written to a pdf)
##run 1 time before using goseq to generate the kegg correlation list
kegg=keggLink("pathway","dme")
temp=names(kegg)
temp2=convert[convert[,5]!="",]
i=1
while(i<=nrow(temp2)){
  if(length(temp[temp==temp2[i,5]])){
    temp[temp==temp2[i,5]]=as.character(temp2[i,4])
  }
  i=i+1
}
names(kegg)=temp
kegg=as.list(kegg[substr(names(kegg),1,3)!="dme"])

##sig genes
setwd(ywd)
##significant genes. FBgn IDs
genelist=row.names(ORSR)
genelist=setdiff(genelist,row.names(convert[convert$mt==1,]))
genelist=intersect(row.names(logFCdata)[logFCdata[,17]=="1"],genelist)

##the output file = csvname_KEGGEnrichmentAnalysis.csv
csvname="/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/testesGenes"
##backgroud genes to be used. an impulseDE output file
data=oocontrol
data=(data[setdiff(row.names(data),row.names(convert[convert$mt==1,])),])
data=data[,c(1,2)]
data[,1]=as.character(convert[row.names(data),4])
data[data[,1]=="jus",][1,1]="jus2"
data[,2]=0
data[genelist,2]=1
genes=as.integer(data[,2])
names(genes)=data[,1]
table(genes)
pwf=nullp(genes,"dm3","geneSymbol")
head(pwf)

##GO term enrichment is GO.wall
GO.wall=goseq(pwf,"dm3","geneSymbol")
##write.csv(GO.wall,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/Cluster5_GOenrichment.csv")


##KEGG term enrichment is KEGG
KEGG=goseq(pwf,gene2cat=kegg)
KEGG$adjp=p.adjust(KEGG$over_represented_pvalue,method="BH")
row.names(KEGG)=substr(as.character(KEGG[,1]),9,13)
##write.csv(KEGG,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/testesGenes.csv")

##To rerun heatmaps
##for oxphos figures only
oxphossigs=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/OxphosUnitLabeledORSR.csv",row.names=1)

##list of FBgn
gois=row.names(ORSR)
gois=setdiff(gois,row.names(convert[convert$mt==1,]))
##gois=intersect(row.names(oxphossigs)[oxphossigs[,6]=="5"],gois)

gois=intersect(row.names(logFCdata)[logFCdata[,17]=="5"],gois)
##list of 1 or more keggs
##cluster 5 keggs
##oatkegg=c("00310","00190","03050","00071","00280","04146","00061","00640","04136")
##cluster 1 keggs
##oatkegg=c("00010","00190","00620","00230","04341","03050","00020","00052","04145")
##ok=1
keggs="00190"
hmtitle=""

##Generates a kgenes file for the given kegg category with excessive gene info
temp=convert[convert$mt==1,]
gois=setdiff(gois,row.names(temp))
if(length(keggs)==1){
kgenes=genesinkegg(gois,keggs)
}
if(length(keggs)>1){
  
##For more than 1 kegg category
i=1
while(i<11){
kgenes=genesinkegg(gois,keggs[i])
if(i==1){
  sigkgenes=kgenes
}
if(i>1){
sigkgenes=rbind(sigkgenes,kgenes)
}
i=i+1
}

kgenes=sigkgenes[as.character(unique(sigkgenes[,1])),]
}

i=1
while(i<=nrow(kgenes)){
  kgenes$cluster[i]=logFCdata[row.names(kgenes)[i],17]
  i=i+1
}
dim(kgenes)


##plots the kgenes into a heatmap
heatmapdata=keggheatmapinfo(kgenes)
##temp=plothmdata(heatmapdata)
col.pan <- colorpanel(100, "black","white","red")
heatmapped=heatmap.2(as.matrix(heatmapdata[[2]]),Rowv=F,Colv=F,RowSideColors = heatmapdata[[3]],trace="none",col=col.pan,scale="row",key=TRUE,cexRow = .8,cexCol = 1,main=hmtitle)
##heatmapped=heatmap.2(as.matrix(heatmapdata[[2]]),Rowv=F,Colv=F,RowSideColors = heatmapdata[[3]],trace="none",col=col.pan,scale="row",key=TRUE,cexRow = .4,cexCol = 1,main=hmtitle)
##,colsep=c(4,8,12),rowsep=heatmapdata[[4]])
##,sepwidth=c(.1,.1),sepcolor="blue",main=hmtitle)
##legend("left",      
##legend = unique(heatmapdata[[5]]),
##col = unique(heatmapdata[[3]]), 
##lty= 1,             
##lwd = 5,           
##cex=.7)

##ok=ok+1

##line graphs for kegg cats
##Uses the heatmap data of the kegg category organized by MBClusterSeq Clusters and makes line plots for each gene cluster
clusterdata=as.matrix(t(heatmapped$carpet))
allclusters=unique(heatmapdata[[5]])
i=1
while(i<=length(allclusters)){
if(length((row.names(clusterdata)[heatmapdata[[5]]==allclusters[i]]))>1){  
  usethese=(heatmapdata[[5]]==allclusters[i])
  usethese=usethese[c(length(usethese):1)]
plotdata=clusterdata[usethese,]
min=min(plotdata)
max=max(plotdata)
plotdata=plotdata[,c(1:4,NA,5:8,NA,9:12,NA,13:16)]

plot(x=c(1:19),plotdata[1,],type="l",col="grey", ylim=c(min,max),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,19),xlab="")
title(main=paste(hmtitle,"\n",allclusters[i]))
title(xlab="Samples Over Time",line=4) 
axis(1,at=c(1:4),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(6:9),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(11:14),labels=c(0,1,2,4),col.axis="Black")
axis(1,at=c(16:19),labels=c(0,1,2,4),col.axis="Black")

axis(1,at=c(2.5),labels=c("Control"),padj=2,tick=F,col.axis="blue")
axis(1,at=c(7.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="red")
axis(1,at=c(12.5),labels=c("Control"),padj=2,tick=F,col.axis="blue")
axis(1,at=c(17.5),labels=c("Rapamycin"),padj=2,tick=F,col.axis="red")

axis(1,at=c(5),labels=c("OreR;OreR"),padj=3.5,tick=F,col.axis="black")
axis(1,at=c(15),labels=c("sm21;OreR"),padj=3.5,tick=F,col.axis="black")


j=2
while(j<=nrow(plotdata)){
lines(plotdata[j,],type="l",col="grey",lwd=1)
j=j+1
}
line1=c(mean(plotdata[,1]),mean(plotdata[,2]),mean(plotdata[,3]),mean(plotdata[,4]))
line2=c(mean(plotdata[,6]),mean(plotdata[,7]),mean(plotdata[,8]),mean(plotdata[,9]))
line3=c(mean(plotdata[,11]),mean(plotdata[,12]),mean(plotdata[,13]),mean(plotdata[,14]))
line4=c(mean(plotdata[,16]),mean(plotdata[,17]),mean(plotdata[,18]),mean(plotdata[,19]))


  
lines(line1,type="l",col="blue",lwd=2)
lines(x=c(6:9),line2,type="l",col="red",lwd=2)
lines(x=c(11:14),line3,type="l",col="blue",lwd=2)
lines(x=c(16:19),line4,type="l",col="red",lwd=2)


abline(h=0, col="black", lty=2, lwd=1)
}
i=i+1
}




###################################################
##rerun mbclusterseq using CPM data not raw data...
###################################################
##get the EdgeR output for normalization and trimming
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")

singlecounts=singlecounts[,trim]
groups=groups[trim,]
edger=function(counts, groups){
  x<-DGEList(counts)
  class(x)
  geneid <- rownames(x)
  x$samples=cbind(x$samples,groups)
  atleast=(10000000/min(x$samples$lib.size))
  keep=rowSums(cpm(x)>=atleast)>=3
  y=x[keep, ,keep.lib.sizes=FALSE]
  z <- calcNormFactors(y, method = "TMM") 
  return(z)
}

##uses the edger function that is defined below (needs to be ran before using) to trim out the genes with low read counts
x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 

##create heatmaps and files
geneset=allsigs

counts=cpm(z$counts)[geneset,]

groups=z$samples
order=c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)

counts=counts[,order]
groups=groups[order,]
groups$Time=c(rep(c(0,0,1,1,2,2,4,4),4))
groups$Treat[c(9,10,25,26)]="Rapa"
groups$Group=paste(substr(groups$Sample,1,2),substr(groups$Time,1,1),"H",substr(groups$Treat,1,1),sep="")
GeneID=row.names(counts)  
Normalizer=groups$norm.factors
Treatment=groups$Group
Treatment=paste(substr(Treatment,1,2),substr(Treatment,5,5),substr(Treatment,3,4),sep="")

mydata=RNASeq.Data(counts,Normalizer,Treatment,GeneID) 

c0=KmeansPlus.RNASeq(mydata,nK=5)$centers

cls=Cluster.RNASeq(data=mydata,model="nbinom",centers=c0,method="DA")$cluster

tr=Hybrid.Tree(data=mydata,cluster0 =cls,model="nbinom") 

image=plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)

logFCdata=cbind(mydata[[6]],cls) 
colnames(logFCdata)=c(levels(as.factor(Treatment)),"Cluster")
##pdf("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_Heatmap.pdf")
plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)
##dev.off()

##dir.create("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_data")
setwd("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_data/")
##save(mydata,file="AllSignificantGenes_forcolors_MBClusterSeqOutput_mydata.RData")
##pdf("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_Heatmap.pdf")
plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)
##dev.off()
##write.csv(c0,"AllSignificantGenes_forcolors_MBClusterSeqOutput_c0.csv")
##write.csv(cls,"AllSignificantGenes_forcolors_MBClusterSeqOutput_cls.csv")
##write.csv(tr,"AllSignificantGenes_forcolors_MBClusterSeqOutput_tr.csv")
##write.csv(logFCdata,"AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv")
setwd(ywd)


##Cluster Line Graph
##gets the mean cluster data for line plots
clusterdata=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv",row.names=1)
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
par(mar=c(5,5,1,1),bg="transparent")
plot(allclusters[1,],col="red", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=limits,lwd=2)
##plot(x=c(1:5),clusterdata[i,c(1:5)],type="l",col="grey", ylim=c(min,max),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,10),xlab="")
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

##plots mean cluster values as boxplots connected by lines
cluster1$blank=NA
cluster1=cluster1[,c(1:4,18,5:8,18,9:12,18,13:16)]
timelabs=c("0H","1H","2H","4H")


boxplot(cluster1,outcex=0,ylim=c(-1.5,1),xaxt="n",xlab="")
lines(allclusters[1,],type="l",col="red",lwd=2)
title(xlab="Samples Over Time",line=4,cex=2) 
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(6:9),labels=timelabs,col.axis="Black")
axis(1,at=c(11:14),labels=timelabs,col.axis="Black")
axis(1,at=c(16:19),labels=timelabs,col.axis="Black")

axis(1,at=c(2.55),labels=c("OreR;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("OreR;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(17.5),labels=c("sm21;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
abline(h=0, col="black", lty=2, lwd=1)
##axis(1,at=c(5),labels=c("OreR;OreR"),padj=3,tick=F,col.axis="black")
##axis(1,at=c(15),labels=c("sm21;OreR"),padj=3,tick=F,col.axis="black")

cluster2$blank=NA
cluster2=cluster2[,c(1:4,18,5:8,18,9:12,18,13:16)]
timelabs=c("0H","1H","2H","4H")


boxplot(cluster2,outcex=0,ylim=c(-.75,.75),xaxt="n",xlab="")
lines(allclusters[2,],type="l",col="blue",lwd=2)
title(xlab="Samples Over Time",line=4,cex=2) 
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(6:9),labels=timelabs,col.axis="Black")
axis(1,at=c(11:14),labels=timelabs,col.axis="Black")
axis(1,at=c(16:19),labels=timelabs,col.axis="Black")

axis(1,at=c(2.55),labels=c("OreR;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("OreR;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(17.5),labels=c("sm21;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
abline(h=0, col="black", lty=2, lwd=1)


cluster3$blank=NA
cluster3=cluster3[,c(1:4,18,5:8,18,9:12,18,13:16)]
timelabs=c("0H","1H","2H","4H")


boxplot(cluster3,outcex=0,ylim=c(-.75,.75),xaxt="n",xlab="")
lines(allclusters[3,],type="l",col="black",lwd=2)
title(xlab="Samples Over Time",line=4,cex=2) 
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(6:9),labels=timelabs,col.axis="Black")
axis(1,at=c(11:14),labels=timelabs,col.axis="Black")
axis(1,at=c(16:19),labels=timelabs,col.axis="Black")

axis(1,at=c(2.55),labels=c("OreR;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("OreR;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(17.5),labels=c("sm21;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
abline(h=0, col="black", lty=2, lwd=1)


##########################
##Tissue specific heatmaps
##########################

##give a list of genes (FBgn ID) and a graph title
##tissue1=intersect(row.names(ORSR),row.names(logFCdata)[logFCdata[,17]==1|logFCdata[,17]=="1"|logFCdata[,17]=="1"])
tissue1=intersect(row.names(ORSR),row.names(logFCdata)[logFCdata[,17]==5])
graphtitle="All Sigs Cluster 1"
hmtitle="All Sigs Cluster 1"
##plots the tissue specific expression for the list of genes in a heatmap
##loads the formatted tissues specific data sets
tissuedata=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/FlyAtlas2Data.csv",row.names=1)
tissueid=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/FlyAtlasTissueIDKey.csv",row.names=2)
tissueid=tissueid[,c(2:4)]

##generates a table (exampleid) with tissue specific FPKM for each gene in the tissue1 list
##Also generates a table (enriched) which is the FPKM normalized to the whole fly FPKM. Method is used by FlyAtlas2
exampleid=as.matrix(tissueid)
enriched=exampleid
i=1
while(i<=length(tissue1)){
gene=tissue1[i]
example=tissuedata[tissuedata[,1]==gene,]

exampleid=cbind(exampleid,0)
colnames(exampleid)[length(colnames(exampleid))]=gene
exampleid[as.character(example$TissueID),gene]=as.vector(example$FPKM)
exampleid[,gene]=as.numeric(exampleid[,gene])

enriched=cbind(enriched,0)
colnames(enriched)[length(colnames(enriched))]=gene
enriched[as.character(example$TissueID),gene]=as.vector(example$FPKM)


enriched[enriched[,2]=="'Female'",gene]=as.numeric(enriched[enriched[,2]=="'Female'",gene])/as.numeric(enriched["200",gene])
enriched[enriched[,2]=="'Male'",gene]=as.numeric(enriched[enriched[,2]=="'Male'",gene])/as.numeric(enriched["100",gene])
enriched[enriched[,2]=="'Larval'",gene]=as.numeric(enriched[enriched[,2]=="'Larval'",gene])/as.numeric(enriched["300",gene])
enriched[enriched[,gene]=="Nan"|enriched[,gene]=="Inf",gene]=0
enriched[,gene]=as.numeric(enriched[,gene])
i=i+1
}

##Organizes exampleid data for heatmap (hmdata)
hmdata=matrix(as.numeric(exampleid[,-c(1:3)]),ncol=(ncol(exampleid)-3))
row.names(hmdata)=row.names(exampleid)
colnames(hmdata)=colnames(exampleid)[-c(1:3)]
hmdata=t(hmdata)
hmnames=(paste(exampleid[,1],exampleid[,2],exampleid[,3],sep=" "))
colnames(hmdata)=hmnames

displayorder=c(c(1:42)[exampleid[,2]=="'Both'"],c(1:42)[exampleid[,2]=="'Male'"],c(1:42)[exampleid[,2]=="'Female'"])


##plots hmdata
col.pan <- colorpanel(100, "black","white","red")
heatmap.2(hmdata[,displayorder],Rowv=T,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,main=graphtitle,cexCol = .5,cexRow=.7,colCol=c(rep("dark green",9),rep("blue",16),rep("red",17)))

##Organizes enriched data for heatmap (hmdata)
hmdata=matrix(as.numeric(enriched[,-c(1:3)]),ncol=(ncol(enriched)-3))
row.names(hmdata)=row.names(enriched)
colnames(hmdata)=colnames(enriched)[-c(1:3)]
hmdata=t(hmdata)
hmnames=(paste(substr(enriched[,1],2,2),substr(enriched[,2],2,2)," ",enriched[,3],sep=""))
colnames(hmdata)=hmnames

##plots the larval, male and female tissue expression data
displayorder=c(c(1:42)[exampleid[,2]=="'Both'"],c(1:42)[exampleid[,2]=="'Male'"],c(1:42)[exampleid[,2]=="'Female'"])
col.pan <- colorpanel(100, "black","white","red")
heatmap.2(hmdata[,displayorder],Rowv=T,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,main=graphtitle,cexCol = .5,cexRow=.5,colCol=c(rep("dark green",9),rep("blue",16),rep("red",17)))

##Plots only male and female tissues, not larval
displayorder=c(c(1:42)[exampleid[,2]=="'Male'"],c(1:42)[exampleid[,2]=="'Female'"])
heatmap.2(na.omit(hmdata[,displayorder]),Rowv=T,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,main=graphtitle,cexCol = .4,cexRow=.5,colCol=c(rep("blue",16),rep("red",17)))

##Plot with rows using symbol instead of FBgn
symboldevstages=hmdata
rowi=1
while(rowi<=nrow(symboldevstages)){
  row.names(symboldevstages)[rowi]=as.character(convert[row.names(symboldevstages)[rowi],4])
  rowi=rowi+1
}
heatmap.2(na.omit(symboldevstages[,displayorder]),Rowv=T,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,main=graphtitle,cexCol = .4,cexRow=.35,colCol=c(rep("blue",16),rep("red",17)))

##Plot MALES ONLY with rows using symbol instead of FBgn
symboldevstages=hmdata
rowi=1
while(rowi<=nrow(symboldevstages)){
  row.names(symboldevstages)[rowi]=as.character(convert[row.names(symboldevstages)[rowi],4])
  rowi=rowi+1
}
displayorder=c(c(1:42)[exampleid[,2]=="'Male'"])
heatmap.2(na.omit(symboldevstages[,displayorder]),Rowv=T,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,main=graphtitle,cexCol = .4,cexRow=.35)


##Preparing the raw FlyAtlas2 data to a usable format
##original data downloaded from: http://motif.gla.ac.uk/downloads/FlyAtlas2_19.10.15.sql
##fpkm=read.csv("/Users/johnsantiago/Desktop/FlyAtlas2Data.csv",header=F)

##organized=fpkm[,c(1:8)]
##i=1
##while((i*8)<=ncol(fpkm)){
  ##temporganized=organized
  ##tempfpkm=fpkm[,c(((8*i)+1):((8*i)+8))]
  ##colnames(tempfpkm)=colnames(temporganized)
  ##organized=rbind(temporganized,tempfpkm)
  ##print(i)
  ##print(dim(organized))
  ##i=i+1
##}

##colnames(organized)=c("FBgn","TissueID","FPKM","Rep1","Rep2","Rep3","SD","Status")
##write.csv(organized,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/FlyAtlas2Data.csv")

##ids=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/FlyAtlas2tissueid.csv",header=F)

##idtable=matrix(c(rep(0,378)),ncol=9)

##j=1
  ##while(j<=42){
    ##idtable[j,1]=as.character(ids[((j-1)*9)+1,1])
    ##idtable[j,2]=as.character(ids[((j-1)*9)+2,1])
    ##idtable[j,3]=as.character(ids[((j-1)*9)+3,1])
    ##idtable[j,4]=as.character(ids[((j-1)*9)+4,1])
    ##idtable[j,5]=as.character(ids[((j-1)*9)+5,1])
    ##idtable[j,6]=as.character(ids[((j-1)*9)+6,1])
    ##idtable[j,7]=as.character(ids[((j-1)*9)+7,1])
    ##idtable[j,8]=as.character(ids[((j-1)*9)+8,1])
    ##idtable[j,9]=as.character(ids[((j-1)*9)+9,1])
    ##j=j+1
  ##}

##write.csv(idtable,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/FlyAtlasTissueIDKey.csv")




bpds=na.omit(symboldevstages[,displayorder])
boxplot(bpds,main=hmtitle,cex.axis=.5,col="dodgerblue")

i=1
while(i<16){
tempbpds=bpds[,c(i,i)]
if(i==14){
  tempbpds[,1]="testes" 
}

if(i!=14){
  tempbpds[,1]="other" 
}

tempbpds[,1]=colnames(bpds)[i]
colnames(tempbpds)=c("sample","Expression")
if(i==1){
  bpdsaov=tempbpds
}
if(i>1){
  bpdsaov=rbind(bpdsaov,tempbpds)
}
i=i+1
}

##bpdsaov=as.table(bpdsaov)
##bpdsaov[,2]=as.numeric(bpdsaov[,2])
tissueaov=aov(bpdsaov[,2]~bpdsaov[,1])
tukey=TukeyHSD(tissueaov)

i=1
while(i<16){
  print(mean(bpds[,i]))
##  print(t.test(bpds[,14],bpds[,i]))
  i=i+1
}

######################################
##making developmental stages heatmaps
######################################

##create a list of genes (FBgn ID) for the heatmap
devstage1=row.names(kgenes[kgenes$cluster=="1",])
hmtitle="All Sigs Cluster 1"
devstage1=intersect(row.names(ORSR),row.names(logFCdata)[logFCdata[,17]==1])
##plots the development stage specific expression for the list of genes in a heatmap
##load the modencode data set downloadd from https://github.com/modENCODE-DCC/www/tree/master/html/docs/flyscores
devstages=read.csv("/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/devstages.csv",header=T,row.names=1)

##remove cell line data
devstagesonly=devstages[,c(1:4,30:53,55,57,59,54,56,58)]

##generates a heatmap of the development stage expression levels for genes in devstage1
hmdevstages=na.omit(devstagesonly[as.character(devstage1),])
col.pan <- colorpanel(100, "black","white","red")
heatmap.2(as.matrix(hmdevstages[,c(5:34)]),Rowv=T,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,main=hmtitle,cexRow = .7)

##Use symbols for row names instead of FBgn
symboldevstages=hmdevstages
rowi=1
while(rowi<=nrow(symboldevstages)){
  row.names(symboldevstages)[rowi]=as.character(convert[row.names(symboldevstages)[rowi],4])
  rowi=rowi+1
}
heatmap.2(as.matrix(symboldevstages[,c(29:34)]),Rowv=T,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,main=hmtitle,cexRow = .35)

##cluster5 dev stage data
##c5ds=symboldevstages[,c(29:34)]
##c1ds=symboldevstages[,c(29:34)]
c5males=c(c5ds[,1],c5ds[,2],c5ds[,3])
c5females=c(c5ds[,4],c5ds[,5],c5ds[,6])

c1males=c(c1ds[,1],c1ds[,2],c1ds[,3])
c1females=c(c1ds[,4],c1ds[,5],c1ds[,6])


boxplot(bpds,main=hmtitle,cex.axis=.5,col=c(rep("dodgerblue",3),rep("firebrick",3)))




##run 1 time before trying to plot kegg data
##generates kegg, a list with names as [symbol] and kegg category as the [[data]]
##kegg correlation list
kegg=keggLink("pathway","dme")
temp=names(kegg)
temp2=convert[convert[,5]!="",]
i=1
while(i<=nrow(temp2)){
  if(length(temp[temp==temp2[i,5]])){
    temp[temp==temp2[i,5]]=as.character(temp2[i,4])
  }
  i=i+1
}
names(kegg)=temp
kegg=as.list(kegg[substr(names(kegg),1,3)!="dme"])

genesinkegg=function(gois,keggs){
  
  convert=read.csv(("FBgnConversionTable.csv"),row.names=1)
  temp=convert
  temp$inKEGG=0
  kegg2=as.list(names(kegg))
  names(kegg2)=unlist(kegg)
  inkegg=as.vector(unlist(kegg2[names(kegg2)==paste("path:dme",keggs,sep="")]))
  info=(temp[inkegg,])
  i=1
  while(i<=length(inkegg)){
    temp$inKEGG[temp$Symbol==inkegg[i]]=1
    i=i+1
  }
  kgenes=temp[gois,]
  kgenes=kgenes[kgenes$inKEGG==1,]
  return(kgenes)
}

##takes the genesinkegg output and generates a heatmap of the individual librarys and a heatmap of the library means
keggheatmapinfo=function(kgenes){
  cpmtable=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
  groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)
  hmcolors=kgenes$cluster
  hmdata=cpmtable[row.names(kgenes),]
  row.names(hmdata)=kgenes$Symbol
  hmdata=hmdata[,c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)]
  groups=groups[c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28),]
  meanhmdata=hmdata[,c(1:14)]
  meanhmdata[,1]=(hmdata[,1]+hmdata[,2])/2
  meanhmdata[,2]=(hmdata[,3]+hmdata[,4])/2
  meanhmdata[,3]=(hmdata[,5]+hmdata[,6])/2
  meanhmdata[,4]=(hmdata[,7]+hmdata[,8])/2
  meanhmdata[,5]=(hmdata[,9]+hmdata[,10])/2
  meanhmdata[,6]=(hmdata[,11]+hmdata[,12])/2
  meanhmdata[,7]=(hmdata[,13]+hmdata[,14])/2
  meanhmdata[,8]=(hmdata[,15]+hmdata[,16])/2
  meanhmdata[,9]=(hmdata[,17]+hmdata[,18])/2
  meanhmdata[,10]=(hmdata[,19]+hmdata[,20])/2
  meanhmdata[,11]=(hmdata[,21]+hmdata[,22])/2
  meanhmdata[,12]=(hmdata[,23]+hmdata[,24])/2
  meanhmdata[,13]=(hmdata[,25]+hmdata[,26])/2
  meanhmdata[,14]=(hmdata[,27]+hmdata[,28])/2
  meanhmdata[,15]=(hmdata[,29]+hmdata[,30])/2
  meanhmdata[,16]=(hmdata[,31]+hmdata[,32])/2
  meanhmdata=meanhmdata[order(hmcolors),]
  hmdata=hmdata[order(hmcolors),]
  hmcolors=hmcolors[order(hmcolors)]
  cluster=paste("cluster",hmcolors)
  colnames(hmdata)=groups$Group
  colnames(meanhmdata)=groups$Group[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]
  basecolors=c("red","green","black","yellow","blue")
  hmcolors=basecolors[(hmcolors)]
  row.sep=c(1:(length(unique(cluster))-1))
  i=1
  while(i<length(unique(cluster))){
    row.sep[i]=length(cluster[cluster==unique(cluster)[i]])
    i=i+1
  }
  i=2
  while(i<=length(row.sep)){
    row.sep[i]=row.sep[i]+row.sep[i-1]
    i=i+1
  }
  
  heatmapdata=list(hmdata,meanhmdata,hmcolors,row.sep,cluster)
  names(heatmapdata)=c("hmdata","meanhmdata","hmcolors","row.sep","cluster")
  return(heatmapdata)
}


##turns the data into a heatmap
plothmdata=function(heatmapdata){
  col.pan <- colorpanel(100, "black","white","red")
  heatmap.2(as.matrix(heatmapdata[[2]]),Rowv=F,Colv=F,RowSideColors = heatmapdata[[3]],trace="none",col=col.pan,scale="row",key=TRUE,colsep=c(4,8,12),rowsep=heatmapdata[[4]],sepwidth=c(.1,.1),sepcolor="blue",main=hmtitle)
  ##legend("left",      
         ##legend = unique(heatmapdata[[5]]),
         ##col = unique(heatmapdata[[3]]), 
         ##lty= 1,             
         ##lwd = 5,           
         ##cex=.7)
}


##Prior to 01132020
##your working directory (ywd) is the path to access the R outputs folder
ywd="/Users/johncsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/FlyRefeedingPaper/"
setwd(ywd)

#############
####Libraries
#############
library(ImpulseDE2)
library(edgeR)
library(goseq)
library(gplots)
library(pca3d)
library(MBCluster.Seq)
library(org.Dm.eg.db)
library(RcisTarget)
library(stats)
library(KEGGREST)

###############
####Quick Files
###############
oocc=read.csv("SingleCountOreRCC_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
oorapa=read.csv("SingleCountOreRRapa_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
oocontrol=read.csv("SingleCountOreRControl_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
socc=read.csv("SingleCountsm21OreRCC_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
sorapa=read.csv("SingleCountsm21OreRRapa_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
socontrol=read.csv("SingleCountsm21OreRControl_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
oocsoc=read.csv("/Users/johncsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/FlyRefeedingPaper/SingleCount_OreROreRvssm21OreRControl_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
oorsor=read.csv("SingleCount_OreROreRvssm21OreRRapa_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)
cpmsc=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
OCSC=na.omit(oocsoc[oocsoc$padj<.05,])
ORSR=na.omit(oorsor[oorsor$padj<.05,])
OO=na.omit(oocc[oocc$padj<.05,])
SO=na.omit(socc[socc$padj<.05,])
SOR=na.omit(sorapa[sorapa$padj<.05,])
SOC=na.omit(socontrol[socontrol$padj<.05,])
OOR=na.omit(oorapa[oorapa$padj<.05,])
OOC=na.omit(oocontrol[oocontrol$padj<.05,])
allsigs=unique(c(row.names(OCSC),row.names(ORSR),row.names(OO),row.names(SO),row.names(OOC),row.names(OOR),row.names(SOC),row.names(SOR)))
trim=c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE,  TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)
convert=read.csv("FBgnConversionTable.csv",row.names=1)
cpmtable=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)
logFCdata=read.csv("/Users/johncsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/Chapter1Data/AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv",row.names=1)


countdata=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
colnames(countdata)[1:9]=paste(substring(colnames(countdata)[1:9],1,2),substring(colnames(countdata)[1:9],4,4),sep="")
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv",row.names=1)
colnames(countdata)=row.names(groups[colnames(countdata),])
x <- countdata
group <- factor(groups$Group)
y <- DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE] 
z <- calcNormFactors(y, method = "TMM") 
cpmdata=cpm(z)
design<-model.matrix(~0+group)
colnames(design) <- levels(group)
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)
compare = makeContrasts((SO1HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
G_X_E<-topTags(lrt,adjust.method="BH",n = nrow(z$counts), sort.by="PValue")
fullSO1R=G_X_E$table

trim=c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE,  TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)
countdata=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
colnames(countdata)[1:9]=paste(substring(colnames(countdata)[1:9],1,2),substring(colnames(countdata)[1:9],4,4),sep="")
countdata=countdata[,trim]
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv",row.names=1)
groups=groups[trim,]
colnames(countdata)=row.names(groups[colnames(countdata),])
x <- countdata
group <- factor(groups$Group)
y <- DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE] 
z <- calcNormFactors(y, method = "TMM") 
cpmdata=cpm(z)
design<-model.matrix(~0+group)
colnames(design) <- levels(group)
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)
compare = makeContrasts((SO1HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
G_X_E<-topTags(lrt,adjust.method="BH",n = nrow(z$counts), sort.by="PValue")
trimSO1R=G_X_E$table

trim=c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE,  TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)
countdata=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
colnames(countdata)[1:9]=paste(substring(colnames(countdata)[1:9],1,2),substring(colnames(countdata)[1:9],4,4),sep="")
countdata=countdata[,trim]
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv",row.names=1)
groups=groups[trim,]
colnames(countdata)=row.names(groups[colnames(countdata),])
countdata=countdata[,15:28]
groups=groups[15:28,]
x <- countdata
group <- factor(groups$Group)
y <- DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE] 
z <- calcNormFactors(y, method = "TMM") 
cpmdata=cpm(z)
design<-model.matrix(~0+group)
colnames(design) <- levels(group)
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)
compare = makeContrasts((SO1HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
G_X_E<-topTags(lrt,adjust.method="BH",n = nrow(z$counts), sort.by="PValue")
onlySO1R=G_X_E$table

sfull=fullSO1R[fullSO1R$FDR<.05,]
strim=trimSO1R[trimSO1R$FDR<.05,]
sonly=onlySO1R[onlySO1R$FDR<.05,]

dim(sfull)
dim(strim)
dim(sonly)


##Chapter 1 Figures

##Figure 2 A-D
##Runs EdgeR on the libraries. Then do individual comparisons for heatmap
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")

singlecounts=singlecounts[,trim]
groups=groups[trim,]
groups=as.matrix(groups)
groups=as.data.frame(groups[c(1:28),])


x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 
sc=z
cpmsc=cpm(z)
##write.csv(cpmsc,"SingleCountReads_NormalizedCPM.csv")
##write.csv(sc$counts,"SingleCountReads.csv")
##write.csv(groups,"SingleCountReads_Metadata.csv")


##To run an edger comparison
design<-model.matrix(~0+groups$Group)
colnames(design) <- levels(factor(groups$Group))
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)

##pdf("SO1HR_Volcano.pdf",width = 6,height = 6)
##title="SO1HR"
compare = makeContrasts((SO4HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
G_X_E<-topTags(lrt, n=14003,adjust.method="BH", sort.by="PValue")
##hist(G_X_E$table$PValue, breaks=100,main=(G_X_E$comparison))
results=G_X_E$table
write.csv(results,"/Users/johnsantiago/Google Drive File Stream/My Drive/Santiago_RandLab_DigitalNotebook/Thesis Files/Chapter 1/BMCGenomics/Supplementary/TableS2L_SO4HR-SO0HC.csv")

##VOLCANO!!
results=G_X_E$table
notsig=results
sig=results[results$FDR<=.05,]
sig$FDR=-log10(sig$FDR)
sig$dist=0
sig$dist=abs(sig$logFC)+(sig$FDR)
sig$dist=(sig$dist/10)*100
sig$dist[sig$dist>100]=100
col.pan <- colorpanel(100, "red","yellow","blue")

plot(x=notsig$logFC,y=-log10(notsig$FDR),cex=.3,main=title,xlab="log Fold Change",ylab="-log10 FDR",col="black",pch=19,xlim=c(-5,5),ylim=c(0,10))
##points(x=sig$logFC,y=(sig$FDR),pch=20,cex=.4,col=col.pan[sig$dist], bg=col.pan[sig$dist])
points(x=sig$logFC,y=(sig$FDR),pch=20,cex=.4,col="red")
abline(v=0,lty=2,col="black")
abline(h=(-log10(.05)),lty=2,col="red")
dev.off()


##Figure 2 E-G
##Figures for whole transcriptome trends section in draft
##Total DE genes within a condition bargraph 
setwd("/Users/johncsantiago/Google Drive File Stream/My Drive/Sanders Lab Files/Santiago Manuscripts/GenomeBiologyPaper/Revisions/Code for Revisions/")
pdf("Figure1C.pdf",height=6,width=6)
numgene=c(nrow(OOC),nrow(OOR),nrow(SOC),nrow(SOR))
mp=barplot(numgene,col =c(rep("#004C9966",2),rep("#FF000080",2)),border =c(rep("#000066",3),rep("#990000",3),NA,rep("#000066",3),rep("#990000",3)),ylab="Total Genes",ylim=c(0,5000),
           main="Total Genes Differentially Expressed\nIn Each Condition")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c(rep(c("Control","Rapamycin"),2)),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3,3.7),labels=c("OreR;OreR","sm21;OreR"),padj=.6,tick=F,cex.axis=1)
dev.off()

pdf("Figure1C.pdf",height=6,width=6)
numgene=c(nrow(OOC),nrow(OOR),nrow(SOC),nrow(SOR))
mp=barplot(numgene,col =c(rep("blue",2),rep("red",2)),ylab="Total Genes",ylim=c(0,5000),
           main="Total Genes Differentially Expressed\nIn Each Condition")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c(rep(c("Control","Rapamycin"),2)),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3,3.7),labels=c("OreR;OreR","sm21;OreR"),padj=.6,tick=F,cex.axis=1)
dev.off()

##Total DE genes between control and rapa for each genotype
pdf("Figure1D.pdf",height=6,width=4)
numgene=c(nrow(OO),nrow(SO))
mp=barplot(numgene,col =c(("#004C9966"),("#FF000080")),border =c(rep("#000066",3),rep("#990000",3),NA,rep("#000066",3),rep("#990000",3)),ylab="Total Genes",ylim=c(0,5000),
           main="Treatment Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c((c("OreR;OreR","sm21;OreR"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.35),labels=c("Control ~ Rapamycin"),padj=.6,tick=F,cex.axis=1)
dev.off()

pdf("Figure1D.pdf",height=6,width=4)
numgene=c(nrow(OO),nrow(SO))
mp=barplot(numgene,col =c(("blue"),("red")),ylab="Total Genes",ylim=c(0,5000),
           main="Treatment Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c((c("OreR;OreR","sm21;OreR"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.35),labels=c("Control ~ Rapamycin"),padj=.6,tick=F,cex.axis=1)
dev.off()

##Total DE genes between genotype for each treatment
pdf("Figure1E.pdf",height=6,width=4)
numgene=c(nrow(OCSC),nrow(ORSR))
mp=barplot(numgene,col =c(("#004C9966"),("#FF000080")),border =c(rep("#000066",3),rep("#990000",3),NA,rep("#000066",3),rep("#990000",3)),ylab="Total Genes",ylim=c(0,5000),
           main="Genotype Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c((c("Control","Rapamycin"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3),labels=c("OreR;OreR ~ sm21;OreR"),padj=.6,tick=F,cex.axis=1)
dev.off()

pdf("Figure1E.pdf",height=6,width=4)
numgene=c(nrow(OCSC),nrow(ORSR))
mp=barplot(numgene,col =c(("blue"),("red")),ylab="Total Genes",ylim=c(0,5000),
           main="Genotype Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c((c("Control","Rapamycin"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3),labels=c("OreR;OreR ~ sm21;OreR"),padj=.6,tick=F,cex.axis=1)
dev.off()

##Figure 2 H-J
##PCA data
cpmtable=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)


temp=c(row.names(OO),row.names(SO),row.names(OOR),row.names(OOC),row.names(SOR),row.names(SOC),row.names(OCSC),row.names(ORSR))
temp=unique(temp)
##temp=unique(row.names(cpmtable))
temp2=cpmtable[temp,]
temp3=temp2[,c(1,2,1,2,15,16,15,16)]
temp4=temp2[,c(3,4,5,6,17,18,19,20)]
row.names(temp4)=paste(row.names(temp4),"1hour",sep="")
colnames(temp4)=colnames(temp3)
temp3=rbind(temp3,temp4)
temp4=temp2[,c(7,8,9,10,21,22,23,24)]
row.names(temp4)=paste(row.names(temp4),"2hour",sep="")
colnames(temp4)=colnames(temp3)
temp3=rbind(temp3,temp4)
temp4=temp2[,c(11,12,13,14,25,26,27,28)]
row.names(temp4)=paste(row.names(temp4),"4hour",sep="")
colnames(temp4)=colnames(temp3)
temp3=rbind(temp3,temp4)
tempgroups=groups[c(1,2,1,2,15,16,15,16),]
tempgroups$Group=c("OOC","OOC","OOR","OOR","SOC","SOC","SOR","SOR")
tempgroups$Treat[c(3,4,7,8)]="Rapa"
tempgroups$Replicate=c(1,2,1,2,1,2,1,2)
row.names(tempgroups)=paste(tempgroups$Group,tempgroups$Replicate,sep="")
tempgroups=tempgroups[,-3]
colnames(temp3)=row.names(tempgroups)

combodata=t(temp3)
pca <- prcomp(combodata, scale.=TRUE) 
gr <- factor(tempgroups$Group)

pca3d(pca, group=gr,show.group.labels=T,palette=c("blue","red","green","black"),show.plane=F,show.centroids = T) 
pca2d(pca, group=gr, show.group.labels = T, palette=c("blue","red","green","black"))
pca2d(pca, components=c(2,3),group=gr, show.group.labels = T, palette=c("blue","red","green","black"),show.centroids = T)
pca2d(pca, components=c(1,3),group=gr, show.group.labels = T, palette=c("blue","red","green","black"))
pca2d(pca, components=c(3,4),group=gr, show.group.labels = T, palette=c("blue","red","green","black"))

##proportion explained by each pca
summary(pca)
##scree plot
plot(pca)

##Figure 3A

tr=read.csv("AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_tr.csv", row.names=1)
cls=read.csv("AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_cls.csv", row.names=1)
load(file="AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_mydata.RData")

##pdf("AllSignificantGenes_forcolors_MBClusterSeq_Heatmap.pdf", width=6, height=12)
plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)
##dev.off()

##Figure 4





##############################
##Fly Count Table Construction
##############################

##Generating filter vector called trim that removes the most outlying replicate from each treatment group
##Using distance between libraries in an MDS plot to determine which library 
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")
edger=function(counts, groups){
  x<-DGEList(counts)
  class(x)
  geneid <- rownames(x)
  x$samples=cbind(x$samples,groups)
  atleast=(10000000/min(x$samples$lib.size))
  keep=rowSums(cpm(x)>=atleast)>=3
  y=x[keep, ,keep.lib.sizes=FALSE]
  z <- calcNormFactors(y, method = "TMM") 
  return(z)
}

##uses the edger function that is defined below (needs to be ran before using) to trim out the genes with low read counts
x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 

MDSdata=z$counts
MDSgroups=groups

lcpm=cpm(MDSdata,log=TRUE)
colors=(groups$Group)
label=groups$Group
levels(colors)=c(1:length(levels(colors)))
temp=plotMDS(lcpm, col=as.vector(colors),labels=label,gene.selection="pairwise",top=nrow(lcpm))

similarity=temp$distance.matrix
samples=unique(groups$Group)

similaritytable=groups[,c(6,5)]
similaritytable$rep1=c(rep(c(1,1,2),14))
similaritytable$rep2=c(rep(c(2,3,3),14))
similaritytable$Replicate=c(rep(c("1~2","1~3","2~3"),14))
similaritytable$distance=0
i=1
j=1
while(i<=length(samples)){
  coords=groups$Group==samples[i]
  temp=similarity[coords,coords]
  similaritytable$distance[j]=temp[2,1]
  j=j+1
  similaritytable$distance[j]=temp[3,1]
  j=j+1
  similaritytable$distance[j]=temp[3,2]
  j=j+1
  i=i+1
}

tempkeep=similaritytable[similaritytable$Group==samples[1],]
tempkeep=tempkeep[tempkeep$distance==min(tempkeep$distance),]
keep=tempkeep
i=2
while(i<=length(samples)){
  tempkeep=similaritytable[similaritytable$Group==samples[i],]
  tempkeep=tempkeep[tempkeep$distance==min(tempkeep$distance),]
  keep=rbind(keep,tempkeep)
  i=i+1
}

trim=c(rep(0,42))
i=1
while(i<15){
  trim=trim+(groups$Group==keep[i,1]&groups$Replicate==keep[i,3])
  trim=trim+(groups$Group==keep[i,1]&groups$Replicate==keep[i,4])
  i=i+1
}

trim=trim==1

##Runs EdgeR on the libraries. Saves normalized CPM file and a count table both with low count reads removed
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")

singlecounts=singlecounts[,trim]
groups=groups[trim,]

x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 
sc=z
cpmsc=cpm(z)
##write.csv(cpmsc,"SingleCountReads_NormalizedCPM.csv")
##write.csv(sc$counts,"SingleCountReads.csv")
##write.csv(groups,"SingleCountReads_Metadata.csv")


##To run an edger comparison
design<-model.matrix(~0+groups$Group)
colnames(design) <- levels(groups$Group)
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)
compare = makeContrasts((SO4HR-OO4HR), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
G_X_E<-topTags(lrt, n=14003,adjust.method="BH", sort.by="PValue")
hist(G_X_E$table$PValue, breaks=100,main=(G_X_E$comparison))


##VOLCANO!!
results=G_X_E$table
notsig=results
sig=results[results$FDR<=.05,]
##pdf("SO1HC_Volcano.pdf",width = 6,height = 12)
plot(x=notsig$logFC,y=-log10(notsig$FDR),cex=.3,main="OreR;OreR\n1 Hour Control",xlab="log Fold Change",ylab="-log10 FDR",col="black",pch=19)
points(x=sig$logFC,y=-log10(sig$FDR),cex=.4,col="red")
abline(v=0,lty=2,col="black")
abline(v=-2,lty=2,col=)
cpmdata["FBgn0039475",]

##dev.off()
nrow(sig)


#######
##EdgeR overview of kegg metabolism over time
#######
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")

singlecounts=singlecounts[,trim]
groups=groups[trim,]

x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 
sc=z
cpmsc=cpm(z)



design<-model.matrix(~0+groups$Group)
colnames(design) <- levels(groups$Group)
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)

##OO1HC
compare = makeContrasts((OO1HC-OO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

edgerdata=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(edgerdata)=row.names(temp)
colnames(edgerdata)=c("rownames.OO1HC","OO1HC.FC","OO1HC.FDR")

##OO2HC
compare = makeContrasts((OO2HC-OO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by ="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.OO2HC","OO2HC.FC","OO2HC.FDR")
edgerdata=cbind(edgerdata,ted)

##OO4HC
compare = makeContrasts((OO4HC-OO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.OO4HC","OO4HC.FC","OO4HC.FDR")
edgerdata=cbind(edgerdata,ted)


##OO1HR
compare = makeContrasts((OO1HR-OO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.OO1HR","OO1HR.FC","OO1HR.FDR")
edgerdata=cbind(edgerdata,ted)

##OO2HR
compare = makeContrasts((OO2HR-OO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by ="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.OO2HR","OO2HR.FC","OO2HR.FDR")
edgerdata=cbind(edgerdata,ted)

##OO4HR
compare = makeContrasts((OO4HR-OO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.OO4HR","OO4HR.FC","OO4HR.FDR")
edgerdata=cbind(edgerdata,ted)


##SO1HC
compare = makeContrasts((SO1HC-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.SO1HC","SO1HC.FC","SO1HC.FDR")
edgerdata=cbind(edgerdata,ted)

##SO2HC
compare = makeContrasts((SO2HC-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by ="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.SO2HC","SO2HC.FC","SO2HC.FDR")
edgerdata=cbind(edgerdata,ted)

##SO4HC
compare = makeContrasts((SO4HC-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.SO4HC","SO4HC.FC","SO4HC.FDR")
edgerdata=cbind(edgerdata,ted)

##SO1HR
compare = makeContrasts((SO1HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.SO1HR","SO1HR.FC","SO1HR.FDR")
edgerdata=cbind(edgerdata,ted)

##SO2HR
compare = makeContrasts((SO2HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by ="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.SO2HR","SO2HR.FC","SO2HR.FDR")
edgerdata=cbind(edgerdata,ted)

##SO4HR
compare = makeContrasts((SO4HR-SO0HC), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
temp<-topTags(lrt, n=14003,adjust.method="BH", sort.by="none")
temp=temp$table[order(row.names(temp$table)),]

ted=cbind(row.names(temp),temp$logFC,temp$FDR)
row.names(ted)=row.names(temp)
colnames(ted)=c("rownames.SO4HR","SO4HR.FC","SO4HR.FDR")
edgerdata=cbind(edgerdata,ted)

##write.csv(edgerdata,"EdgeRTimepointFCandFDR.csv")

###############
#### ImpulseDE 
###############
##Fly Ore;OreR case-control
##Setup input libraries and metadata
counts=read.csv("SingleCountReads.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)
counts=counts[,c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)]
groups=groups[c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28),]

counts=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Full%20Count%20Table.csv",row.names=1)
groups=read.csv("https://raw.githubusercontent.com/DavidRandLab/Santiago-et-al-2021-BMC-Genomics/main/Data%20Files/Full%20Count%20Table%20Meta%20Data.csv",row.names=1)
counts=counts[,row.names(groups)]
counts=counts[,trim]
groups=groups[trim,]
order=c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)
counts=counts[,order]
groups=groups[order,]
genotype=substr(groups$Mito,1,1)=="O"
counts=counts[,genotype]
groups=groups[genotype,]

counts=as.matrix(counts)

counts=counts[,groups$Mito=="OreR"]
groups=groups[groups$Mito=="OreR",]
groups[,"Condition"]=c(rep("control",8),rep("case",8))
groups[,"Sample"]=colnames(counts)
groups[,"Time"]=c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4)
groups[,"Batch"]=c(rep(c(1,2),8))
row.names(groups)=groups$Sample
groups=groups[,c(1,7,4,9)]

##running impulse de2 without including transiently expressed genes
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = groups,
  boolCaseCtrl    = TRUE,
  vecConfounders  = NULL,
  scaNProc        = 1 )
write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCountOreRCC_TrimmedLibraries_ImpulseDEResuults.csv")
oocc=read.csv("SingleCountOreRCC_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)

##Fly Ore;OreR Rapa Case only
case=(groups$Condition)=="case"

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts[,case], 
  dfAnnotation    = groups[case,],
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCountOreRRapa_TrimmedLibraries_ImpulseDEResuults.csv")
oorapa=read.csv("SingleCountOreRRapa_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)

##Fly Ore;OreR Control Case only
control=(groups$Condition)=="control"
tempmat=counts[,control]
tempdf=groups[control,]

tempdf$Condition="case"

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = tempmat, 
  dfAnnotation    = tempdf,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCountOreRControl_TrimmedLibraries_ImpulseDEResuults.csv")
oocontrol=read.csv("SingleCountOreRControl_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)


##Fly sm21;OreR case-control
##Setup input libraries and metadata
counts=read.csv("SingleCountReads.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)
counts=counts[,c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)]
groups=groups[c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28),]
genotype=substr(groups$Mito,1,1)=="s"
counts=counts[,genotype]
groups=groups[genotype,]

counts=as.matrix(counts)

groups[,"Condition"]=c(rep("control",8),rep("case",8))
groups[,"Sample"]=colnames(counts)
groups[,"Time"]=c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4)
groups[,"Batch"]=c(rep(c(1,2),8))
row.names(groups)=groups$Sample
groups=groups[,c(1,7,4,8)]

##running impulse de2 without including transiently expressed genes
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = groups,
  boolCaseCtrl    = TRUE,
  vecConfounders  = NULL,
  scaNProc        = 1 )
write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCountsm21OreRCC_TrimmedLibraries_ImpulseDEResuults.csv")
socc=read.csv("SingleCountsm21OreRCC_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)

##Fly sm21;OreR Rapa Case only
case=(groups$Condition)=="case"

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts[,case], 
  dfAnnotation    = groups[case,],
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCountsm21OreRRapa_TrimmedLibraries_ImpulseDEResuults.csv")
sorapa=read.csv("SingleCountsm21OreRRapa_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)

##Fly sm21;OreR Control Case only
control=(groups$Condition)=="control"
tempmat=counts[,control]
tempdf=groups[control,]

tempdf$Condition="case"

objectImpulseDE2 <- runImpulseDE2(
  matCountData    = tempmat, 
  dfAnnotation    = tempdf,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCountsm21OreRControl_TrimmedLibraries_ImpulseDEResuults.csv")
socontrol=read.csv("SingleCountsm21OreRControl_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)

##Fly Ore;OreR vs sm21;OreR Control
##Setup input libraries and metadata
counts=read.csv("SingleCountReads.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)
counts=counts[,c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)]
groups=groups[c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28),]
groups$Treat[c(9,10,25,26)]="Rapa"

counts=counts[,groups$Treat=="Control"]
groups=groups[groups$Treat=="Control",]

counts=as.matrix(counts)

groups[,"Condition"]=c(rep("control",8),rep("case",8))
groups[,"Sample"]=colnames(counts)
groups[,"Time"]=c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4)
groups[,"Batch"]=c(rep(c(1,2),8))
row.names(groups)=groups$Sample
groups=groups[,c(1,7,4,8)]

##running impulse de2 without including transiently expressed genes
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = groups,
  boolCaseCtrl    = TRUE,
  vecConfounders  = NULL,
  scaNProc        = 1 )
write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCount_OreROreRvssm21OreRControl_TrimmedLibraries_ImpulseDEResuults.csv")
oocsoc=read.csv("SingleCount_OreROreRvssm21OreRControl_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)

##Fly Ore;OreR vs sm21;OreR Rapa
##Setup input libraries and metadata
counts=read.csv("SingleCountReads.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)

counts=counts[,c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)]
groups=groups[c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28),]
groups$Treat[c(9,10,25,26)]="Rapa"

counts=counts[,groups$Treat=="Rapa"]
groups=groups[groups$Treat=="Rapa",]

counts=as.matrix(counts)

groups[,"Condition"]=c(rep("control",8),rep("case",8))
groups[,"Sample"]=colnames(counts)
groups[,"Time"]=c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4)
groups[,"Batch"]=c(rep(c(1,2),8))
row.names(groups)=groups$Sample
groups=groups[,c(1,7,4,8)]

##running impulse de2 without including transiently expressed genes
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts, 
  dfAnnotation    = groups,
  boolCaseCtrl    = TRUE,
  vecConfounders  = NULL,
  scaNProc        = 1 )
write.csv(objectImpulseDE2$dfImpulseDE2Results,"SingleCount_OreROreRvssm21OreRRapa_TrimmedLibraries_ImpulseDEResuults.csv")
oorsor=read.csv("SingleCount_OreROreRvssm21OreRRapa_TrimmedLibraries_ImpulseDEResuults.csv",row.names=1)

##Figures for whole transcriptome trends section in draft
##Total DE genes within a condition bargraph 
numgene=c(nrow(OOC),nrow(OOR),nrow(SOC),nrow(SOR))
mp=barplot(numgene,col =c(rep("#004C9966",2),rep("#FF000080",2)),border =c(rep("#000066",3),rep("#990000",3),NA,rep("#000066",3),rep("#990000",3)),ylab="Total Genes",ylim=c(0,max(numgene,na.rm = T)),
main="Total Genes Differentially Expressed\nIn Each Condition")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c(rep(c("Control","Rapamycin"),2)),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3,3.7),labels=c("OreR;OreR","sm21;OreR"),padj=.6,tick=F,cex.axis=1)


##Total DE genes between control and rapa for each genotype
numgene=c(nrow(OO),nrow(SO))
mp=barplot(numgene,col =c(("#004C9966"),("#FF000080")),border =c(rep("#000066",3),rep("#990000",3),NA,rep("#000066",3),rep("#990000",3)),ylab="Total Genes",ylim=c(0,4000),
           main="Treatment Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c((c("OreR;OreR","sm21;OreR"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.35),labels=c("Control ~ Rapamycin"),padj=.6,tick=F,cex.axis=1)

##Total DE genes between genotype for each treatment
numgene=c(nrow(OCSC),nrow(ORSR))
mp=barplot(numgene,col =c(("#004C9966"),("#FF000080")),border =c(rep("#000066",3),rep("#990000",3),NA,rep("#000066",3),rep("#990000",3)),ylab="Total Genes",ylim=c(0,max(numgene,na.rm = T)),
           main="Genotype Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(mp),labels=c((c("Control","Rapamycin"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3),labels=c("OreR;OreR ~ sm21;OreR"),padj=.6,tick=F,cex.axis=1)

##PCA data
cpmtable=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)


temp=c(row.names(OO),row.names(SO),row.names(OOR),row.names(OOC),row.names(SOR),row.names(SOC),row.names(OCSC),row.names(ORSR))
temp=unique(temp)
##temp=unique(row.names(cpmtable))
temp2=cpmtable[temp,]
temp3=temp2[,c(1,2,1,2,15,16,15,16)]
temp4=temp2[,c(3,4,5,6,17,18,19,20)]
row.names(temp4)=paste(row.names(temp4),"1hour",sep="")
colnames(temp4)=colnames(temp3)
temp3=rbind(temp3,temp4)
temp4=temp2[,c(7,8,9,10,21,22,23,24)]
row.names(temp4)=paste(row.names(temp4),"2hour",sep="")
colnames(temp4)=colnames(temp3)
temp3=rbind(temp3,temp4)
temp4=temp2[,c(11,12,13,14,25,26,27,28)]
row.names(temp4)=paste(row.names(temp4),"4hour",sep="")
colnames(temp4)=colnames(temp3)
temp3=rbind(temp3,temp4)
tempgroups=groups[c(1,2,1,2,15,16,15,16),]
tempgroups$Group=c("OOC","OOC","OOR","OOR","SOC","SOC","SOR","SOR")
tempgroups$Treat[c(3,4,7,8)]="Rapa"
tempgroups$Replicate=c(1,2,1,2,1,2,1,2)
row.names(tempgroups)=paste(tempgroups$Group,tempgroups$Replicate,sep="")
tempgroups=tempgroups[,-3]
colnames(temp3)=row.names(tempgroups)

combodata=t(temp3)
pca <- prcomp(combodata, scale.=TRUE) 
gr <- factor(tempgroups$Group)

pca3d(pca, group=gr,show.group.labels=T,palette=c("blue","red","green","black"),show.plane=F,show.centroids = T) 
pca2d(pca, group=gr, show.group.labels = T, palette=c("blue","red","green","black"))
pca2d(pca, components=c(2,3),group=gr, show.group.labels = T, palette=c("blue","red","green","black"),show.centroids = T)
pca2d(pca, components=c(1,3),group=gr, show.group.labels = T, palette=c("blue","red","green","black"))
pca2d(pca, components=c(3,4),group=gr, show.group.labels = T, palette=c("blue","red","green","black"))

##proportion explained by each pca
summary(pca)
##scree plot
plot(pca)

#########################
##Master Conversion Table
#########################
##making a file with all conversions. fbgn, eg, symbol
convert=matrix(row.names(oocontrol),ncol=1)
convert=convert[,c(1,1,1)]
x=as.list(org.Dm.egFLYBASE2EG)

i=1
while(i<=nrow(convert)){
  if(length(x[[convert[i,2]]])>0){
    convert[i,2]=(x[[convert[i,2]]])
  }
  i=i+1
}
convert=convert[,c(1,3,2,3)]
colnames(convert)=c("original FBgn","updated FBgn","Entrez/NCBI gene ID","Symbol")
##write.csv(convert,"FBgnConversionTable.csv")


xx=as.list(org.Dm.egSYMBOL)
i=1
while(i<=nrow(convert)){
  if(length(xx[[convert[i,3]]])>0){
    convert[i,4]=(xx[[convert[i,3]]])
  }
  i=i+1
}
##write.csv(convert,"FBgnConversionTable.csv")

xxx=keggConv("dme","ncbi-geneid")
convert$'KEGG ID'=paste("ncbi-geneid:",convert$`Entrez/NCBI gene ID`,sep="")
i=1

convert[i,]
i
i=i+1
while(i<=nrow(convert)){
  if(length(xxx[[convert[i,5]]])>0){
    convert[i,5]=(xxx[[convert[i,5]]])
  }
  i=i+1
}
convert[substr(convert[,5],1,3)=="ncb",5]=""
##write.csv(convert,"FBgnConversionTable.csv")
convert=read.csv("FBgnConversionTable.csv",row.names=1)
convert$cluster=0
clusterdata=read.csv("AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv",row.names=1)
sigs=intersect(row.names(clusterdata),row.names(convert))
i=1
while(i<=length(sigs)){
  convert$cluster[row.names(convert)==sigs[i]]=clusterdata$cls[row.names(clusterdata)==sigs[i]]
  i=i+1
}
##write.csv(convert,"FBgnConversionTable.csv")

###################
#### GOseq analysis
###################
##gets the KEGG enrichment values from a geneset and writes to a pdf (also does GO enrichment but not written to a pdf)
##run 1 time before using goseq to generate the kegg correlation list
kegg=keggLink("pathway","dme")
temp=names(kegg)
temp2=convert[convert[,5]!="",]
i=1
while(i<=nrow(temp2)){
  if(length(temp[temp==temp2[i,5]])){
    temp[temp==temp2[i,5]]=as.character(temp2[i,4])
  }
  i=i+1
}
names(kegg)=temp
kegg=as.list(kegg[substr(names(kegg),1,3)!="dme"])

##sig genes
setwd(ywd)
##significant genes. FBgn IDs
genelist=row.names(tempsigs)
genelist=setdiff(genelist,row.names(convert[convert$mt==1,]))
##the output file = csvname_KEGGEnrichmentAnalysis.csv
csvname="tempsigs"
##backgroud genes to be used. an impulseDE output file
data=oocontrol
data=(data[setdiff(row.names(data),row.names(convert[convert$mt==1,])),])
data=data[,c(1,2)]
data[,1]=as.character(convert[row.names(data),4])
data[data[,1]=="jus",][1,1]="jus2"
data[,2]=0
data[genelist,2]=1
genes=as.integer(data[,2])
names(genes)=data[,1]
table(genes)
pwf=nullp(genes,"dm3","geneSymbol")
head(pwf)

##GO term enrichment is GO.wall
GO.wall=goseq(pwf,"dm3","geneSymbol")

##KEGG term enrichment is KEGG
KEGG=goseq(pwf,gene2cat=kegg)
KEGG$adjp=p.adjust(KEGG$over_represented_pvalue,method="BH")
row.names(KEGG)=substr(as.character(KEGG[,1]),9,13)
write.csv(KEGG,paste(csvname,"_KEGGEnrichmentAnalysis.csv"))



########################################
##KEGG Supervised heatmaps and pca plots
########################################
##run 1 time before trying to plot kegg data
##generates kegg, a list with names as [symbol] and kegg category as the [[data]]
##kegg correlation list
kegg=keggLink("pathway","dme")
temp=names(kegg)
temp2=convert[convert[,5]!="",]
i=1
while(i<=nrow(temp2)){
  if(length(temp[temp==temp2[i,5]])){
    temp[temp==temp2[i,5]]=as.character(temp2[i,4])
  }
  i=i+1
}
names(kegg)=temp
kegg=as.list(kegg[substr(names(kegg),1,3)!="dme"])

genesinkegg=function(gois,keggs){
  
convert=read.csv(("FBgnConversionTable.csv"),row.names=1)
temp=convert
temp$inKEGG=0
kegg2=as.list(names(kegg))
names(kegg2)=unlist(kegg)
inkegg=as.vector(unlist(kegg2[names(kegg2)==paste("path:dme",keggs,sep="")]))
info=(temp[inkegg,])
i=1
while(i<=length(inkegg)){
  temp$inKEGG[temp$Symbol==inkegg[i]]=1
  i=i+1
}
kgenes=temp[gois,]
kgenes=kgenes[kgenes$inKEGG==1,]
return(kgenes)
}

##takes the genesinkegg output and generates a heatmap of the individual librarys and a heatmap of the library means
keggheatmapinfo=function(kgenes){
cpmtable=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)
hmcolors=kgenes$cluster
hmdata=cpmtable[row.names(kgenes),]
row.names(hmdata)=kgenes$Symbol
hmdata=hmdata[,c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)]
groups=groups[c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28),]
meanhmdata=hmdata[,c(1:14)]
meanhmdata[,1]=(hmdata[,1]+hmdata[,2])/2
meanhmdata[,2]=(hmdata[,3]+hmdata[,4])/2
meanhmdata[,3]=(hmdata[,5]+hmdata[,6])/2
meanhmdata[,4]=(hmdata[,7]+hmdata[,8])/2
meanhmdata[,5]=(hmdata[,9]+hmdata[,10])/2
meanhmdata[,6]=(hmdata[,11]+hmdata[,12])/2
meanhmdata[,7]=(hmdata[,13]+hmdata[,14])/2
meanhmdata[,8]=(hmdata[,15]+hmdata[,16])/2
meanhmdata[,9]=(hmdata[,17]+hmdata[,18])/2
meanhmdata[,10]=(hmdata[,19]+hmdata[,20])/2
meanhmdata[,11]=(hmdata[,21]+hmdata[,22])/2
meanhmdata[,12]=(hmdata[,23]+hmdata[,24])/2
meanhmdata[,13]=(hmdata[,25]+hmdata[,26])/2
meanhmdata[,14]=(hmdata[,27]+hmdata[,28])/2
meanhmdata[,15]=(hmdata[,29]+hmdata[,30])/2
meanhmdata[,16]=(hmdata[,31]+hmdata[,32])/2
meanhmdata=meanhmdata[order(hmcolors),]
hmdata=hmdata[order(hmcolors),]
hmcolors=hmcolors[order(hmcolors)]
cluster=paste("cluster",hmcolors)
colnames(hmdata)=groups$Group
colnames(meanhmdata)=groups$Group[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]
basecolors=c("red","blue","black","yellow","green")
hmcolors=basecolors[(hmcolors)]
row.sep=c(1:(length(unique(cluster))-1))
i=1
while(i<length(unique(cluster))){
  row.sep[i]=length(cluster[cluster==unique(cluster)[i]])
  i=i+1
}
i=2
while(i<=length(row.sep)){
  row.sep[i]=row.sep[i]+row.sep[i-1]
  i=i+1
}

heatmapdata=list(hmdata,meanhmdata,hmcolors,row.sep,cluster)
names(heatmapdata)=c("hmdata","meanhmdata","hmcolors","row.sep","cluster")
return(heatmapdata)
}


##turns the data into a heatmap
plothmdata=function(heatmapdata){
  col.pan <- colorpanel(100, "black","white","red")
heatmap.2(as.matrix(heatmapdata[[2]]),Rowv=F,Colv=F,RowSideColors = heatmapdata[[3]],trace="none",col=col.pan,scale="row",key=TRUE,colsep=c(4,8,12),rowsep=heatmapdata[[4]],sepwidth=c(.1,.1),sepcolor="blue",main=hmtitle)
legend("left",      
       legend = unique(heatmapdata[[5]]),
       col = unique(heatmapdata[[3]]), 
       lty= 1,             
       lwd = 5,           
       cex=.7)
}

##list of FBgn
gois=row.names(ORSR)
temp=convert[convert$mt==1,]
gois=setdiff(gois,row.names(temp))
##list of keggs
keggs=c("00190")
##plot title
hmtitle=""

kgenes=genesinkegg(gois,keggs)
i=1
while(i<=nrow(kgenes)){
  kgenes$cluster[i]=logFCdata[row.names(kgenes)[i],17]
  i=i+1
}

##pie for oxphos cluster distribution
piedata=summary(as.factor(kgenes$cluster))
names(piedata)=c("Cluster 1","Cluster 2","Cluster 5")
piecols=c("red","green","blue")
pie(piedata,col=piecols)


heatmapdata=keggheatmapinfo(kgenes)
temp=plothmdata(heatmapdata)

##output for kegg website mapping
kgenes=genesinkegg(gois,keggs)
colors=c("red","blue","purple","yellow","green")
kgenes$colors=colors[kgenes$cluster]
temp=kgenes[,c(5,11)]
write.csv(temp,"temp.csv")

pdf("Figure4D.pdf")
set=c(9:16,25:32,7,8,23,24)
colset=c("black","blue","red","blue","red","blue","red","black","dark green","purple","dark green","purple","dark green","purple")
temp=t(heatmapdata[[1]])
##temp=temp[set,]
pca <- prcomp(temp, scale.=TRUE) 
gr <- factor(row.names(temp))
##pca3d(pca, group=gr,show.group.labels=T)
pca2d(pca, group=gr, show.group.labels = T, palette=colset)
dev.off()
plot(pca)

##heatmap for line graphs
cluster="cluster 2"
col.pan <- colorpanel(100, "black","white","red")
oglinedata=heatmap.2(as.matrix(heatmapdata[[2]]),Rowv=F,Colv=F,RowSideColors = heatmapdata[[3]],trace="none",col=col.pan,scale="row",key=TRUE,colsep=c(4,8,12),rowsep=heatmapdata[[4]],sepwidth=c(.1,.1),sepcolor="blue",main=hmtitle)
legend("left",      
       legend = unique(heatmapdata[[5]]),
       col = unique(heatmapdata[[3]]), 
       lty= 1,             
       lwd = 5,           
       cex=.7)

linedata=oglinedata$carpet
linedata=linedata[c(1:4,NA,5:8,NA,9:12,NA,13:16),]
linedata=linedata[,c(ncol(linedata):1)]
use=heatmapdata[[5]]==cluster
linedata=linedata[,use]
plot(linedata[,1],type="l",ylim=c(min(na.omit(linedata)),max(na.omit(linedata))))
i=1
while(i<=ncol(linedata)){
  lines(x=c(1:19),y=linedata[,i])
  i=i+1
}
abline(h=0,lty=2,col="red")

##output for kegg website mapping
kgenes=genesinkegg(gois,keggs)
colors=c("red","blue","purple","yellow","green")
kgenes$colors=colors[kgenes$cluster]
temp=kgenes[,c(5,11)]
write.csv(temp,"temp.csv")

##pca plots
##individual time points
cpmtable=read.csv("SingleCountReads_NormalizedCPM.csv",row.names=1)
groups=read.csv("SingleCountReads_Metadata.csv",row.names=1)
colnames(cpmtable)=groups$Group
cpmtable2=cpmtable[setdiff(row.names(cpmtable),allsigs),]
cpmtable2=cpmtable[row.names(ORSR),]
cpmtable2=cpmtable[allsigs,]
temp=t(cpmtable2)
pca <- prcomp(temp, scale.=TRUE) 
gr <- factor(row.names(temp))
pcacols=c("black",rep(c("blue","red"),3),"black",rep(c("dark green","purple"),3))
##pca3d(pca, group=gr,show.group.labels=T)
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA)
legend("left", legend = c("OOC","OOR","SOC","SOR"), col = c("red","blue"), lty= c(1,1,2,2), lwd = 1, cex=.7)
proportion=(summary(pca)$importance)[2,]*100
barplot(proportion,ylim=c(0,100),axes = F)
axis(2,at=c((0:20)*5),cex.axis=.5)


set=c(9:16,25:32,7,8,23,24)
par(bg=NA)
colset=c("black","red","red","blue","red","black","purple","purple","dark green","purple")
temp=t(heatmapdata[[1]])
temp=temp[set,]
pca <- prcomp(temp, scale.=TRUE) 
gr <- factor(row.names(temp))

##pca3d(pca, group=gr,show.group.labels=T)
pca2d(pca, group=gr, show.group.labels = T, palette=colset,bg=NA)
plot(pca)

actualgenelevels=heatmapdata$meanhmdata
actualgenelevels$Mean=0
i=1
while(i<=nrow(actualgenelevels)){
  actualgenelevels$Mean[i]=mean(as.numeric(actualgenelevels[i,c(1:16)]))
  i=i+1
}









pdf("ORSR_FattyAcidDegradationPathway_Heatmap.pdf",height = 14)
plothmdata(heatmapdata)
dev.off()


actualgenelevels[,16][actualgenelevels$Mean>600]=600
actualgenelevels$Mean[actualgenelevels$Mean>600]=600
col.pan <- colorpanel(1000, "blue","yellow","red")
heatmap.2((as.matrix(actualgenelevels[,c(16,17)])),Rowv="NULL",Colv="NULL",trace="none",col=col.pan,scale="none",key=TRUE,cexRow = .3)




##Plot function for personalizing
col.pan <- colorpanel(100, "black","white","red")
heatmap.2(as.matrix(heatmapdata[[2]]),Rowv=F,Colv=F,trace="none",col=col.pan,scale="row",key=TRUE,colsep=c(4,8,12),rowsep=heatmapdata[[4]],sepwidth=c(.1,.1),sepcolor="blue",main=hmtitle)
legend("left", legend = unique(cluster), col = unique(hmcolors), lty= 1, lwd = 5, cex=.7)





##pca plots
##individual time points
temp=t(heatmapdata[[1]])
pca <- prcomp(temp, scale.=TRUE) 
gr <- factor(colnames(heatmapdata[[1]]))
pcacols=c("black",rep(c("blue","red"),3),"black",rep(c("dark green","purple"),3))
##pca3d(pca, group=gr,show.group.labels=T)
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols)
plot(pca)
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,components = c(2,3))

##pca plots
##individual time points
temp=t(heatmapdata[[1]])
pca <- prcomp(temp, scale.=TRUE) 
gr <- factor(colnames(heatmapdata[[1]]))
pcacols=c("black",rep(c("blue","red"),3),"black",rep(c("dark green","purple"),3))
##pca3d(pca, group=gr,show.group.labels=T)
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols)
plot(pca)
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,components = c(2,3))

##condensed into genotype*treatment effect
temp=heatmapdata[[1]]
temp2=temp[c(1,2,9,10,17,18,25,26)]
temp3=temp[,c(3,4,11,12,19,20,27,28)]
colnames(temp3)=colnames(temp2)
row.names(temp3)=paste(row.names(temp3),".1hour",sep="")
temp2=rbind(temp2,temp3)
temp3=temp[,c(5,6,13,14,21,22,29,30)]
colnames(temp3)=colnames(temp2)
row.names(temp3)=paste(row.names(temp3),".2hour",sep="")
temp2=rbind(temp2,temp3)
temp3=temp[,c(7,8,15,16,23,24,31,32)]
colnames(temp3)=colnames(temp2)
row.names(temp3)=paste(row.names(temp3),".4hour",sep="")
temp2=rbind(temp2,temp3)
colnames(temp2)=c("OOC","OOC","OOR","OOR","SOC","SOC","SOR","SOR")


temp=t(temp2)
pca <- prcomp(temp, scale.=TRUE) 
gr <- factor(colnames(temp2))
pcacols=c("blue","red","green","purple")
pca3d(pca, group=gr,show.group.labels=T,palette=c("blue","red","green","black"),show.plane=F,show.centroids = T)
pca2d(pca, group=gr, show.group.labels = T, palette=pcacols)

##Metabolism (KEGG ID: 01100) time course image data
tckegg="01212"
##OOC
alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
oc1h=alltimes[alltimes$OO1HC.FDR<=.05,]
tempoc1h=alltimes[alltimes$OO1HC.FDR<=.05,]
oc1h=genesinkegg(row.names(oc1h),tckegg)
oc1h=oc1h[,c(1,5)]
oc1h$color="white"
tempoc1h=tempoc1h[row.names(oc1h),]
oc1h$color[(tempoc1h$OO1HC.FC<0)]="blue"
oc1h$color[(tempoc1h$OO1HC.FC>0)]="red"
write.csv(oc1h,paste(tckegg, "edgertimecoursekegg_oc1h.csv",sep=""))


alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
oc2h=alltimes[alltimes$OO2HC.FDR<=.05,]
tempoc2h=alltimes[alltimes$OO2HC.FDR<=.05,]
oc2h=genesinkegg(row.names(oc2h),tckegg)
oc2h=oc2h[,c(1,5)]
oc2h$color="white"
tempoc2h=tempoc2h[row.names(oc2h),]
oc2h$color[(tempoc2h$OO2HC.FC<0)]="blue"
oc2h$color[(tempoc2h$OO2HC.FC>0)]="red"
write.csv(oc2h,paste(tckegg, "edgertimecoursekegg_oc2h.csv",sep=""))

alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
oc4h=alltimes[alltimes$OO4HC.FDR<=.05,]
tempoc4h=alltimes[alltimes$OO4HC.FDR<=.05,]
oc4h=genesinkegg(row.names(oc4h),tckegg)
oc4h=oc4h[,c(1,5)]
oc4h$color="white"
tempoc4h=tempoc4h[row.names(oc4h),]
oc4h$color[(tempoc4h$OO4HC.FC<0)]="blue"
oc4h$color[(tempoc4h$OO4HC.FC>0)]="red"
write.csv(oc4h,paste(tckegg, "edgertimecoursekegg_oc4h.csv",sep=""))

##OOR
alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
or1h=alltimes[alltimes$OO1HR.FDR<=.05,]
tempor1h=alltimes[alltimes$OO1HR.FDR<=.05,]
or1h=genesinkegg(row.names(or1h),tckegg)
or1h=or1h[,c(1,5)]
or1h$color="white"
tempor1h=tempor1h[row.names(or1h),]
or1h$color[(tempor1h$OO1HR.FC<0)]="blue"
or1h$color[(tempor1h$OO1HR.FC>0)]="red"
write.csv(or1h,paste(tckegg, "edgertimecoursekegg_or1h.csv",sep=""))


alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
or2h=alltimes[alltimes$OO2HR.FDR<=.05,]
tempor2h=alltimes[alltimes$OO2HR.FDR<=.05,]
or2h=genesinkegg(row.names(or2h),tckegg)
or2h=or2h[,c(1,5)]
or2h$color="white"
tempor2h=tempor2h[row.names(or2h),]
or2h$color[(tempor2h$OO2HR.FC<0)]="blue"
or2h$color[(tempor2h$OO2HR.FC>0)]="red"
write.csv(or2h,paste(tckegg, "edgertimecoursekegg_or2h.csv",sep=""))

alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
or4h=alltimes[alltimes$OO4HR.FDR<=.05,]
tempor4h=alltimes[alltimes$OO4HR.FDR<=.05,]
or4h=genesinkegg(row.names(or4h),tckegg)
or4h=or4h[,c(1,5)]
or4h$color="white"
tempor4h=tempor4h[row.names(or4h),]
or4h$color[(tempor4h$OO4HR.FC<0)]="blue"
or4h$color[(tempor4h$OO4HR.FC>0)]="red"
write.csv(or4h,paste(tckegg, "edgertimecoursekegg_or4h.csv",sep=""))

##SOC
alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
sc1h=alltimes[alltimes$SO1HC.FDR<=.05,]
tempsc1h=alltimes[alltimes$SO1HC.FDR<=.05,]
sc1h=genesinkegg(row.names(sc1h),tckegg)
sc1h=sc1h[,c(1,5)]
sc1h$color="white"
tempsc1h=tempsc1h[row.names(sc1h),]
sc1h$color[(tempsc1h$SO1HC.FC<0)]="blue"
sc1h$color[(tempsc1h$SO1HC.FC>0)]="red"
write.csv(sc1h,paste(tckegg, "edgertimecoursekegg_sc1h.csv",sep=""))


alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
sc2h=alltimes[alltimes$SO2HC.FDR<=.05,]
tempsc2h=alltimes[alltimes$SO2HC.FDR<=.05,]
sc2h=genesinkegg(row.names(sc2h),tckegg)
sc2h=sc2h[,c(1,5)]
sc2h$color="white"
tempsc2h=tempsc2h[row.names(sc2h),]
sc2h$color[(tempsc2h$SO2HC.FC<0)]="blue"
sc2h$color[(tempsc2h$SO2HC.FC>0)]="red"
write.csv(sc2h,paste(tckegg, "edgertimecoursekegg_sc2h.csv",sep=""))

alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
sc4h=alltimes[alltimes$SO4HC.FDR<=.05,]
tempsc4h=alltimes[alltimes$SO4HC.FDR<=.05,]
sc4h=genesinkegg(row.names(sc4h),tckegg)
sc4h=sc4h[,c(1,5)]
sc4h$color="white"
tempsc4h=tempsc4h[row.names(sc4h),]
sc4h$color[(tempsc4h$SO4HC.FC<0)]="blue"
sc4h$color[(tempsc4h$SO4HC.FC>0)]="red"
write.csv(sc4h,paste(tckegg, "edgertimecoursekegg_sc4h.csv",sep=""))

##SOR
alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
sr1h=alltimes[alltimes$SO1HR.FDR<=.05,]
tempsr1h=alltimes[alltimes$SO1HR.FDR<=.05,]
sr1h=genesinkegg(row.names(sr1h),tckegg)
sr1h=sr1h[,c(1,5)]
sr1h$color="white"
tempsr1h=tempsr1h[row.names(sr1h),]
sr1h$color[(tempsr1h$SO1HR.FC<0)]="blue"
sr1h$color[(tempsr1h$SO1HR.FC>0)]="red"
write.csv(sr1h,paste(tckegg, "edgertimecoursekegg_sr1h.csv",sep=""))


alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
sr2h=alltimes[alltimes$SO2HR.FDR<=.05,]
tempsr2h=alltimes[alltimes$SO2HR.FDR<=.05,]
sr2h=genesinkegg(row.names(sr2h),tckegg)
sr2h=sr2h[,c(1,5)]
sr2h$color="white"
tempsr2h=tempsr2h[row.names(sr2h),]
sr2h$color[(tempsr2h$SO2HR.FC<0)]="blue"
sr2h$color[(tempsr2h$SO2HR.FC>0)]="red"
write.csv(sr2h,paste(tckegg, "edgertimecoursekegg_sr2h.csv",sep=""))

alltimes=read.csv("EdgeRTimepointFCandFDR.csv",row.names=1)
sr4h=alltimes[alltimes$SO4HR.FDR<=.05,]
tempsr4h=alltimes[alltimes$SO4HR.FDR<=.05,]
sr4h=genesinkegg(row.names(sr4h),tckegg)
sr4h=sr4h[,c(1,5)]
sr4h$color="white"
tempsr4h=tempsr4h[row.names(sr4h),]
sr4h$color[(tempsr4h$SO4HR.FC<0)]="blue"
sr4h$color[(tempsr4h$SO4HR.FC>0)]="red"
write.csv(sr4h,paste(tckegg, "edgertimecoursekegg_sr4h.csv",sep=""))


###################
#######MBClusterseq
###################
##get the EdgeR output for normalization and trimming
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")

singlecounts=singlecounts[,trim]
groups=groups[trim,]
edger=function(counts, groups){
  x<-DGEList(counts)
  class(x)
  geneid <- rownames(x)
  x$samples=cbind(x$samples,groups)
  atleast=(10000000/min(x$samples$lib.size))
  keep=rowSums(cpm(x)>=atleast)>=3
  y=x[keep, ,keep.lib.sizes=FALSE]
  z <- calcNormFactors(y, method = "TMM") 
  return(z)
}

##uses the edger function that is defined below (needs to be ran before using) to trim out the genes with low read counts
x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 

##create heatmaps and files
geneset=allsigs

counts=z$counts[geneset,]

groups=z$samples
order=c(1:4,7:8,11:12,1:2,5:6,9:10,13:14,15:18,21:22,25:26,15:16,19:20,23:24,27:28)

counts=counts[,order]
groups=groups[order,]
groups$Time=c(rep(c(0,0,1,1,2,2,4,4),4))
groups$Treat[c(9,10,25,26)]="Rapa"
groups$Group=paste(substr(groups$Sample,1,2),substr(groups$Time,1,1),"H",substr(groups$Treat,1,1),sep="")
GeneID=row.names(counts)  
Normalizer=groups$norm.factors
Treatment=groups$Group
Treatment=paste(substr(Treatment,1,2),substr(Treatment,5,5),substr(Treatment,3,4),sep="")

mydata=RNASeq.Data(counts,Normalizer,Treatment,GeneID) 

c0=KmeansPlus.RNASeq(mydata,nK=5)$centers

cls=Cluster.RNASeq(data=mydata,model="nbinom",centers=c0,method="DA")$cluster

tr=Hybrid.Tree(data=mydata,cluster0 =cls,model="nbinom") 

image=plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)

logFCdata=cbind(mydata[[6]],cls) 
colnames(logFCdata)=c(levels(as.factor(Treatment)),"Cluster")
##pdf("AllSignificantGenes_forcolors_MBClusterSeq_Heatmap.pdf")
plotHybrid.Tree(merge=tr,cluster=cls,logFC=mydata$logFC,tree.title=NULL,colorful=F)
##dev.off()

##dir.create("AllSignificantGenes_forcolors_MBClusterSeq_data")
setwd("AllSignificantGenes_forcolors_MBClusterSeq_data/")
##save(mydata,file="AllSignificantGenes_forcolors_MBClusterSeqOutput_mydata.RData")
##write.csv(c0,"AllSignificantGenes_forcolors_MBClusterSeqOutput_c0.csv")
##write.csv(cls,"AllSignificantGenes_forcolors_MBClusterSeqOutput_cls.csv")
##write.csv(tr,"AllSignificantGenes_forcolors_MBClusterSeqOutput_tr.csv")
##write.csv(logFCdata,"AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv")
setwd(ywd)

##Cluster Line Graph
##gets the mean cluster data for line plots
clusterdata=read.csv("AllSignificantGenes_forcolors_MBClusterSeq_data/AllSignificantGenes_forcolors_MBClusterSeqOutput_logFCdata.csv",row.names=1)
cluster1=clusterdata[clusterdata$cls==1,]
line1=c(rep(0,16))
i=1
while(i<17){
  line1[i]=mean(cluster1[,i])
  i=i+1
}

cluster2=clusterdata[clusterdata$cls==2,]
line2=c(rep(0,16))
i=1
while(i<17){
  line2[i]=mean(cluster2[,i])
  i=i+1
}
allclusters=rbind(line1,line2)

cluster3=clusterdata[clusterdata$cls==3,]
line3=c(rep(0,16))
i=1
while(i<17){
  line3[i]=mean(cluster3[,i])
  i=i+1
}
allclusters=rbind(allclusters,line3)
cluster4=clusterdata[clusterdata$cls==4,]
line4=c(rep(0,16))
i=1
while(i<17){
  line4[i]=mean(cluster4[,i])
  i=i+1
}
allclusters=rbind(allclusters,line4)
cluster5=clusterdata[clusterdata$cls==5,]
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

##Actually plots it
par(mar=c(5,5,1,1),bg="transparent")
plot(allclusters[1,],col="red", type="l",xaxt="n",xlab="",ylab="Mean Cluster Value",ylim=limits,lwd=2)
##plot(x=c(1:5),clusterdata[i,c(1:5)],type="l",col="grey", ylim=c(min,max),ylab="log Fold Change from the Mean",xaxt="n",xlim=c(1,10),xlab="")
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

lines(allclusters[2,],type="l",col="blue",lwd=2)
lines(allclusters[3,],type="l",col="black",lwd=2)
lines(allclusters[4,],type="l",col="yellow",lwd=2)
lines(allclusters[5,],type="l",col="green",lwd=2)


abline(h=0, col="black", lty=2, lwd=1)


cluster1$blank=NA
cluster1=cluster1[,c(1:4,18,5:8,18,9:12,18,13:16)]
timelabs=c("0H","1H","2H","4H")


boxplot(cluster1,outcex=0,ylim=c(-1.5,1),xaxt="n",xlab="")
lines(allclusters[1,],type="l",col="red",lwd=2)
title(xlab="Samples Over Time",line=4,cex=2) 
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(6:9),labels=timelabs,col.axis="Black")
axis(1,at=c(11:14),labels=timelabs,col.axis="Black")
axis(1,at=c(16:19),labels=timelabs,col.axis="Black")

axis(1,at=c(2.55),labels=c("OreR;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("OreR;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(17.5),labels=c("sm21;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
abline(h=0, col="black", lty=2, lwd=1)
##axis(1,at=c(5),labels=c("OreR;OreR"),padj=3,tick=F,col.axis="black")
##axis(1,at=c(15),labels=c("sm21;OreR"),padj=3,tick=F,col.axis="black")

cluster2$blank=NA
cluster2=cluster2[,c(1:4,18,5:8,18,9:12,18,13:16)]
timelabs=c("0H","1H","2H","4H")


boxplot(cluster2,outcex=0,ylim=c(-.75,.75),xaxt="n",xlab="")
lines(allclusters[2,],type="l",col="blue",lwd=2)
title(xlab="Samples Over Time",line=4,cex=2) 
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(6:9),labels=timelabs,col.axis="Black")
axis(1,at=c(11:14),labels=timelabs,col.axis="Black")
axis(1,at=c(16:19),labels=timelabs,col.axis="Black")

axis(1,at=c(2.55),labels=c("OreR;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("OreR;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(17.5),labels=c("sm21;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
abline(h=0, col="black", lty=2, lwd=1)


cluster3$blank=NA
cluster3=cluster3[,c(1:4,18,5:8,18,9:12,18,13:16)]
timelabs=c("0H","1H","2H","4H")


boxplot(cluster3,outcex=0,ylim=c(-.75,.75),xaxt="n",xlab="")
lines(allclusters[3,],type="l",col="black",lwd=2)
title(xlab="Samples Over Time",line=4,cex=2) 
axis(1,at=c(1:4),labels=timelabs,col.axis="Black")
axis(1,at=c(6:9),labels=timelabs,col.axis="Black")
axis(1,at=c(11:14),labels=timelabs,col.axis="Black")
axis(1,at=c(16:19),labels=timelabs,col.axis="Black")

axis(1,at=c(2.55),labels=c("OreR;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(7.5),labels=c("OreR;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(12.5),labels=c("sm21;OreR\nControl"),padj=1.5,tick=F,col.axis="black")
axis(1,at=c(17.5),labels=c("sm21;OreR\nRapa"),padj=1.5,tick=F,col.axis="black")
abline(h=0, col="black", lty=2, lwd=1)

#############################
####Transcription factor data
#############################

gois=row.names(ORSR)
keggs="00190"
kgenes=genesinkegg(gois,keggs)
gs=kgenes[kgenes$cluster==2,]
gs=as.character(gs[gs$mt==0,4])

gs=row.names(ORSR)
gs=intersect(row.names(ORSR), row.names(logFCdata[logFCdata[,17]==1,]))
gs=convert[gs,4]

setwd("/Users/johnsantiago/Downloads/")
data(motifAnnotations_dmel_v8)
summary(gs)

motifRankings <- importRankings("dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
  
setwd(ywd)
motifEnrichmentTable_wGenes <- cisTarget(gs, motifRankings, motifAnnot=motifAnnotations_dmel_v8)
  
# Save results as a data.frame
enrtf=as.data.frame(motifEnrichmentTable_wGenes)
head(enrtf)
setwd(ywd)
##write.csv(enrtf,"OXPHOS_Cluster2_TranscriptionFactorEnrichement.csv")

##TF analysis
tfsum=read.csv("OXPHOS_Cluster2_TranscriptionFactorEnrichement.csv")
tfdesc=tfsum[tfsum$Description!="",]
##tf1=tfsum[tfsum[,1]==1,]

tffbgn=as.character(tfsum$TF_highConf)
i=1
while(i<=length(tffbgn)){
  if(substring(tffbgn[i],1,1)==" "){
tffbgn[i]=substring(tffbgn[i],2,nchar(tffbgn[i]))
  }
  i=i+1
}
xx <- as.list(org.Dm.egSYMBOL2EG)

i=1
while(i<=length(tffbgn)){
  ##looks for the gene symbol for each entrez id then replaces the entrez id
  if(length(xx[[tffbgn[i]]])>0){
    tffbgn[i]=(xx[[tffbgn[i]]])
  }
  i=i+1
}

tfsum$tffentrez=tffbgn
i="04350"
gik=as.list(org.Dm.egPATH2EG)
tfegs=matrix(unique(tfsum$tffentrez),ncol=1)
row.names(tfegs)=tfegs[,1]
blank=matrix("",nrow=nrow(tfegs),ncol=1)
tfegs=cbind(tfegs,blank)
i=1
j=2
while(i<=nrow(tfegs)){
  if(length(intersect(row.names(tfegs),gik[[i]]))>0){
    temp=intersect(row.names(tfegs),gik[[i]])
    tfegs[temp,j]=names(gik[i])
    tfegs=cbind(tfegs,blank)
    j=j+1
  }
  i=i+1
}

intersect(row.names(tfegs),gik[[i]])
i=i-1

tfbreakdown=function(tfvar){
tfvar=tfvar[tfvar$TF_highConf!="",]
tfvar$TF_highConf=as.character(tfvar$TF_highConf)
i=1
while(i<=nrow(tfvar)){
tfvar$TF_highConf[i]=strsplit(as.character(tfvar$TF_highConf[i])," \\(")[[1]][1]
i=i+1
}

temp=strsplit(as.character(tfvar$TF_highConf),"; ")
i=1
temp2=0
while(i<=length(temp)){
  temp2=c(temp2,rep(i,length(temp[[i]])))
  i=i+1
}
tfvar=tfvar[temp2,]
tfvar$TF_highConf=unlist(temp)
  
temp=as.character(tfvar$enrichedGenes)
temp=strsplit(temp,";")
names(temp)=tfvar$TF_highConf
namestf=1
i=1
while(i<=length(temp)){
  namestf=c(namestf,rep(names(temp)[i],length(temp[[i]])))
  i=i+1
}
namestf=namestf[-1]

temp2=matrix(0,nrow=length(unlist(temp)),ncol=2)
colnames(temp2)=c("Transcription Factor","Enriched Gene")
temp2[,1]=namestf
temp2[,2]=unlist(temp)
return(temp2)
}

tfvar=read.csv("Cluster2_TranscriptionFactorEnrichement.csv")
tftable=tfbreakdown(tfvar)
temp=unique(tftable[,1])
tftable=tftable[tftable[,1]=="ERR",]

temp2=convert[convert$Symbol==tftable[1,2],]
i=2
while(i<=length(temp)){
  temp3=convert[convert$Symbol==tftable[i,2],]
  temp2=rbind(temp2,temp3)
  i=i+1
}
write.csv(temp2,"temp.csv")

#####################################
####  Manahttan Plot of Gene Location
#####################################

xx2=as.list(org.Dm.egCHRLOC)
convert$Chromosome=0
convert$Location=0

convert$Chromosome[1]=names(xx2[[as.character(convert[1,3])]])[1]
convert$Location[1]=xx2[[as.character(convert[1,3])]][1]
i=2
while(i<=nrow(convert)){
  if(length(names(xx2[[as.character(convert[i,3])]])[1])>0){
  convert$Chromosome[i]=names(xx2[[as.character(convert[i,3])]])[1]
  }
  if(length(xx2[[as.character(convert[i,3])]][1])>0){
  convert$Location[i]=xx2[[as.character(convert[i,3])]][1]
  }
  i=i+1
}

##write.csv(convert,"FBgnConversionTable.csv")

temp=convert[row.names(ORSR),]
temp=temp[temp$Chromosome!=0,]


##Runs EdgeR on the libraries. Saves normalized CPM file and a count table both with low count reads removed
singlecounts=read.csv("CombinedCountTable_SingleReadOnly.csv",row.names=1)
groups=read.csv("Fly_Combined_Count_Table_Metadata.csv")

singlecounts=singlecounts[,trim]
groups=groups[trim,]

x<-DGEList(singlecounts)
class(x)
geneid <- rownames(x)
x$samples=cbind(x$samples,groups)

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast=(10000000/min(x$samples$lib.size))
keep=rowSums(cpm(x)>=atleast)>=3
y=x[keep, ,keep.lib.sizes=FALSE]
z <- calcNormFactors(y, method = "TMM") 
sc=z
cpmsc=cpm(z)
##write.csv(cpmsc,"SingleCountReads_NormalizedCPM.csv")
##write.csv(sc$counts,"SingleCountReads.csv")
##write.csv(groups,"SingleCountReads_Metadata.csv")


##To run an edger comparison
design<-model.matrix(~0+groups$Group)
colnames(design) <- levels(groups$Group)
z = estimateGLMCommonDisp(z,design, verbose=T)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)

compare = makeContrasts((SO1HR-OO1HR), levels=design)
lrt <- glmLRT(fit,contrast=as.vector(compare))		
G_X_E<-topTags(lrt, n=14003,adjust.method="BH", sort.by="PValue")
hist(G_X_E$table$PValue, breaks=100,main=(G_X_E$comparison))


##VOLCANO!!
results=G_X_E$table
notsig=results
sig=results[results$FDR<=.05,]
nrow(sig)
##pdf("SO1HR_Manhattans.pdf",width = 6,height = 12)
locate=convert[row.names(sig),]
locate=cbind(locate,sig)
xgenes=locate[locate$Chromosome=="X",]
alllocate=convert[row.names(notsig),]
alllocate=cbind(alllocate,notsig)
alllocate=alllocate[alllocate$Chromosome=="X",]

plot(x=abs(alllocate$Location),y=alllocate$logFC,cex=.4,col="black",main="X")
points(x=abs(xgenes$Location),y=xgenes$logFC,cex=.4,col="red")

hist(x=abs(xgenes$Location),breaks = 1000)
locate=convert[row.names(sig),]
locate=cbind(locate,sig)
xgenes=locate[locate$Chromosome=="Y",]
alllocate=convert[row.names(notsig),]
alllocate=cbind(alllocate,notsig)
alllocate=alllocate[alllocate$Chromosome=="Y",]

plot(x=abs(alllocate$Location),y=alllocate$logFC,cex=.4,col="black",main="Y")
points(x=abs(xgenes$Location),y=xgenes$logFC,cex=.4,col="red")
abline(h = 0,lty=2,col="black")
hist(x=abs(xgenes$Location),breaks = 1000)
locate=convert[row.names(sig),]
locate=cbind(locate,sig)
xgenes=locate[locate$Chromosome=="4",]
alllocate=convert[row.names(notsig),]
alllocate=cbind(alllocate,notsig)
alllocate=alllocate[alllocate$Chromosome=="4",]

plot(x=abs(alllocate$Location),y=alllocate$logFC,cex=.4,col="black",main="4")
points(x=abs(xgenes$Location),y=xgenes$logFC,cex=.4,col="red")
abline(h = 0,lty=2,col="black")
hist(x=abs(xgenes$Location),breaks = 1000)
locate=convert[row.names(sig),]
locate=cbind(locate,sig)
xgenes=locate[locate$Chromosome=="3R",]
alllocate=convert[row.names(notsig),]
alllocate=cbind(alllocate,notsig)
alllocate=alllocate[alllocate$Chromosome=="3R",]

plot(x=abs(alllocate$Location),y=alllocate$logFC,cex=.4,col="black",main="3R")
points(x=abs(xgenes$Location),y=xgenes$logFC,cex=.4,col="red")
abline(h = 0,lty=2,col="black")
hist(x=abs(xgenes$Location),breaks = 1000)
locate=convert[row.names(sig),]
locate=cbind(locate,sig)
xgenes=locate[locate$Chromosome=="3L",]
alllocate=convert[row.names(notsig),]
alllocate=cbind(alllocate,notsig)
alllocate=alllocate[alllocate$Chromosome=="3L",]

plot(x=abs(alllocate$Location),y=alllocate$logFC,cex=.4,col="black",main="3L")
points(x=abs(xgenes$Location),y=xgenes$logFC,cex=.4,col="red")
abline(h = 0,lty=2,col="black")
hist(x=abs(xgenes$Location),breaks = 1000)
locate=convert[row.names(sig),]
locate=cbind(locate,sig)
xgenes=locate[locate$Chromosome=="2R",]
alllocate=convert[row.names(notsig),]
alllocate=cbind(alllocate,notsig)
alllocate=alllocate[alllocate$Chromosome=="2R",]

plot(x=abs(alllocate$Location),y=alllocate$logFC,cex=.4,col="black",main="2R")
points(x=abs(xgenes$Location),y=xgenes$logFC,cex=.4,col="red")
abline(h = 0,lty=2,col="black")
hist(x=abs(xgenes$Location),breaks = 1000)
locate=convert[row.names(sig),]
locate=cbind(locate,sig)
xgenes=locate[locate$Chromosome=="2L",]
alllocate=convert[row.names(notsig),]
alllocate=cbind(alllocate,notsig)
alllocate=alllocate[alllocate$Chromosome=="2L",]

plot(x=abs(alllocate$Location),y=alllocate$logFC,cex=.4,col="black",main="2L")
points(x=abs(xgenes$Location),y=xgenes$logFC,cex=.4,col="red")
abline(h = 0,lty=2,col="black")
hist(x=abs(xgenes$Location),breaks = 1000)
dev.off()
