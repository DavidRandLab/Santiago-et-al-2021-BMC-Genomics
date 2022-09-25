library(visNetwork)
library(goseq)
library(RcisTarget)
library(dplyr)
library(KEGGREST)

##Only download these files 1 time since they are very large
##download.file("https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc8nr/gene_based/dm6-5kb-upstream-full-tx-11species.mc8nr.genes_vs_motifs.rankings.feather","/dm6-5kb-upstream-full-tx-11species.mc8nr.genes_vs_motifs.rankings.feather")
##download.file("https://resources.aertslab.org/cistarget/motif2tf/motifs-v8-nr.flybase-m0.001-o0.0.tbl","/motifs-v8-nr.flybase-m0.001-o0.0.tbl")

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
kegg.names=unique(substr((kegg[,2]),6,13))
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
##RcisTarget to get enriched TFs
data(motifAnnotations_dmel_v8)
motifRankings = importRankings("dm6-5kb-upstream-full-tx-11species.mc8nr.genes_vs_motifs.rankings.feather")

gs=convert[kegg.genes,"Symbol"]

motifEnrichmentTable_wGenes <- cisTarget(gs, motifRankings, motifAnnot=motifAnnotations_dmel_v8)

motifs_AUC <- calcAUC(gs, motifRankings, nCores=1)

##RcisTarget output to setup edge table for visnetwork
tfvar = read.csv("/Users/johncsantiago/Documents/GitHub/Santiago-et-al-2021-BMC-Genomics/Data Files/Cluster 1 TF Analysis.csv",row.names=1)


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
  
tftable=tfbreakdown(tfvar)
Cluster1.TF.edges=tftable[tftable[,1]=="Abd-B",2]
Cluster1.TF.edges=intersect(Cluster1.TF.edges,convert[kegg.genes,"Symbol"])

from=rep("Abd-B",length(Cluster1.TF.edges))
to=as.matrix(Cluster1.TF.edges)
  
i=1
while(i<=length(kegg.genes)){
  temp=intersect(kegg[grep(kegg.genes[i],(kegg[,1])),2],sigCluster1.KEGG.Table$category)
  to=c(to,convert[rep(kegg.genes[i],length(temp)),"Symbol"])
  from=c(from,as.vector(kegg.names[substr(temp,6,13)]))
  i=i+1
}

edges=data.frame("from"=from,"to"=to)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
edgecolors=genes
edgecolors=rep("gray",length(genes))
edgecolors[grep("Abd-B",edges[,1])]="red"
abdorno=edgecolors
abdorno[abdorno=="red"]="Inferred Abd-B Target"
abdorno[abdorno=="gray"]="Not Abd-B Target"

nodes <- data.frame(id=c(genes, motifs),   
                    label=c(genes,motifs),   
                    title=c(genes,motifs),
                    shape=c(rep("dot", length(genes)), rep("star", length(motifs))),
                    color=c(rep("grey",9),edgecolors))
nodes$group=c(abdorno,"Abd-B",rep("KEGG Pathway",8))

visnet.cluster1=visNetwork(nodes[,c(1:3,6)], edges[,c(1:2)]) %>%
  visNodes(font=list(color="black",size=60),mass=3.5) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE)
visnet.cluster1=visGroups(visnet.cluster1, groupname = "KEGG Pathway", shape = "box", color = list(background = "#FEF861", border="#E56E00"),size=70,physics=TRUE)

visnet.cluster1=visGroups(visnet.cluster1, groupname = "Abd-B", shape = "box", color = list(background = '#FB699A', border='#DB0003'),size=80,physics=FALSE,mass=10)
visnet.cluster1=visGroups(visnet.cluster1, groupname = "Inferred Abd-B Target", shape = "dot", color = list(background = '#F8253D', border="black"))
visnet.cluster1=visGroups(visnet.cluster1, groupname = "Not Abd-B Target", shape = "dot", color = list(background = '#FFC942', border="black"))
visLegend(visnet.cluster1, main="Legend", position="right", ncol=1,width=.1) 


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


#### Figure 6F
##RcisTarget to get enriched TFs
data(motifAnnotations_dmel_v8)
motifRankings = importRankings("dm6-5kb-upstream-full-tx-11species.mc8nr.genes_vs_motifs.rankings.feather")

gs=convert[kegg.genes,"Symbol"]

motifEnrichmentTable_wGenes <- cisTarget(gs, motifRankings, motifAnnot=motifAnnotations_dmel_v8)

motifs_AUC <- calcAUC(gs, motifRankings, nCores=1)

##Network analysis using RcisTarget output
tfvar = read.csv("/Users/johncsantiago/Documents/GitHub/Santiago-et-al-2021-BMC-Genomics/Data Files/Cluster 5 TF Analysis.csv",row.names=1)

tftable=tfbreakdown(tfvar)

Cluster5.Dref.edges=tftable[tftable[,1]=="Dref",2]
Cluster5.Dref.edges=intersect(Cluster5.Dref.edges,convert[kegg.genes,"Symbol"])

Cluster5.giant.edges=tftable[tftable[,1]=="gt",2]
Cluster5.giant.edges=intersect(Cluster5.giant.edges,convert[kegg.genes,"Symbol"])

from=c(rep("Dref",length(Cluster5.Dref.edges)),rep("gt",length(Cluster5.giant.edges)))
to=as.matrix(c(Cluster5.Dref.edges,Cluster5.giant.edges))

i=1
while(i<=length(kegg.genes)){
  temp=intersect(kegg[grep(kegg.genes[i],(kegg)[,1]),2],sigCluster5.KEGG.Table$category)
  to=c(to,convert[rep(kegg.genes[i],length(temp)),"Symbol"])
  from=c(from,as.vector(kegg.names[substr(temp,6,13)]))
  i=i+1
}

edges=data.frame("from"=from,"to"=to)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
edgecolors=c(rep("gray",length(genes)))
names(edgecolors)=genes
edgecolors[edges[grep("Dref",edges[,1]),2]]="blue"
edgecolors[edges[grep("gt",edges[,1]),2]]="red"
edgecolors[intersect(edges[grep("Dref",edges[,1]),2],edges[grep("gt",edges[,1]),2])]="purple"


TForno=edgecolors
TForno[edgecolors=="gray"]="Not a Dref or gt Target Gene"
TForno[edgecolors=="blue"]="Dref Target Gene"
TForno[edgecolors=="red"]="gt Target Gene"
TForno[edgecolors=="purple"]="Dref and gt Target Gene"

nodes <- data.frame(id=c(genes, motifs),   
                    label=c(genes,motifs),   
                    title=c(genes,motifs),
                    shape=c(rep("dot", length(genes)), rep("star", length(motifs))),
                    color=c(edgecolors,rep("grey",length(motifs))))
nodes$group=c(TForno,"Dref","gt",rep("KEGG Pathway",9))

visnet.cluster5=visNetwork(nodes[,c(1:3,6)], edges[,c(1:2)]) %>%
  visNodes(font=list(color="black",size=60),mass=3.5) %>%
  visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE)

visnet.cluster5=visGroups(visnet.cluster5, groupname = "KEGG Pathway", shape = "box", color = list(background = "#FEF861", border="#E56E00"),size=70,physics=TRUE)
visnet.cluster5=visGroups(visnet.cluster5, groupname = "gt", shape = "box", color = list(background = '#FB699A', border='#DB0003'),size=80,physics=FALSE,mass=10)
visnet.cluster5=visGroups(visnet.cluster5, groupname = "Dref", shape = "box", color = list(background = '#03BDF6', border='#0350F6'),size=80,physics=FALSE,mass=10)

visnet.cluster5=visGroups(visnet.cluster5, groupname = "gt Target Gene", shape = "dot", color = list(background = '#F8253D', border="black"))
visnet.cluster5=visGroups(visnet.cluster5, groupname = "Dref Target Gene", shape = "dot", color = list(background = 'dodgerblue', border="black"))
visnet.cluster5=visGroups(visnet.cluster5, groupname = "Dref and gt Target Gene", shape = "dot", color = list(background = 'mediumorchid', border="black"))
visnet.cluster5=visGroups(visnet.cluster5, groupname = "Not a Dref or gt Target Gene", shape = "dot", color = list(background = '#FFC942', border="black"))

visLegend(visnet.cluster5, main="Legend", position="right", ncol=1,width=.1) 