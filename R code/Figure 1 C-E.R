
##ImpulseDE2 output tables with significant genes only
OCSC=na.omit(oocsoc[oocsoc$padj<.05,])
ORSR=na.omit(oorsor[oorsor$padj<.05,])
OO=na.omit(oocc[oocc$padj<.05,])
SO=na.omit(socc[socc$padj<.05,])
SOR=na.omit(sorapa[sorapa$padj<.05,])
SOC=na.omit(socontrol[socontrol$padj<.05,])
OOR=na.omit(oorapa[oorapa$padj<.05,])
OOC=na.omit(oocontrol[oocontrol$padj<.05,])

##Figure 1C
##Total DE genes differentially expressed in each condition
dev.new(height=7.7, width=6.6, noRStudioGD = TRUE, units = "inch")
numgene=c(nrow(OOC),nrow(OOR),nrow(SOC),nrow(SOR))
barplot(numgene,col =c(rep("blue",2),rep("red",2)),ylab="Total Genes",ylim=c(0,5000),
        main="Total Genes Differentially Expressed\nIn Each Condition")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(.65,1.9,3.1,4.3),labels=c(rep(c("Control","Rapamycin"),2)),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3,3.7),labels=c("OreR;OreR","sm21;OreR"),padj=.6,tick=F,cex.axis=1)

##Figure 1D
##Total DE genes between control and rapa for each genotype
dev.new(height=7.7, width=6.6, noRStudioGD = TRUE, units = "inch")
numgene=c(nrow(OO),nrow(SO))
barplot(numgene,col =c(("blue"),("red")),ylab="Total Genes",ylim=c(0,5000),
        main="Treatment Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(.65,1.9),labels=c((c("OreR;OreR","sm21;OreR"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.35),labels=c("Control ~ Rapamycin"),padj=.6,tick=F,cex.axis=1)

##Figure 1E
##Total DE genes between OreR;OreR and sm21;OreR for each treatment
dev.new(height=7.7, width=6.6, noRStudioGD = TRUE, units = "inch")
numgene=c(nrow(OCSC),nrow(ORSR))
barplot(numgene,col =c(("blue"),("red")),ylab="Total Genes",ylim=c(0,5000),
        main="Genotype Effect")
axis(1,seq(-1,17),lwd.ticks=0,labels=rep("",19))
axis(1,at=c(.65,1.9),labels=c((c("Control","Rapamycin"))),col.axis="Black",tick=F,cex=.5,padj=-1.5)
axis(1,at=c(1.3),labels=c("OreR;OreR ~ sm21;OreR"),padj=.6,tick=F,cex.axis=1)
