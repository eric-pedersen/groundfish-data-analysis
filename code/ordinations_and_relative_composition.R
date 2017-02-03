#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
library(plyr)
library(dplyr)
library(vegan)
library(bootstrap)
library(RColorBrewer)


#Loading data and functions ####
source("code/functions.R")
DFO_Com = read.csv("data/DFO_Com.csv",stringsAsFactors = F)
DFO_Dataset = read.csv("data/DFO_Dataset.csv",stringsAsFactors = F)
load("data/year_geom_means.Rdata")

Div<-DFO_Dataset$DIV
Year<-factor(DFO_Dataset$Year,ordered=T)

top4_sp = c("GADUS_MORHUA","SEBASTES_MENTELLA",
            "REINHARDTIUS_HIPPOGLOSSOIDES", "HIPPOGLOSSOIDES_PLATESSOIDES")





#NMDS Calculations####
#taxonomy
# Calculate Bray Curtis dissimilarity, perform NMDS, record site scores
# All spp
diss.tax<-vegdist(decostand(Year_Geom_Means,"total"),"bray")
mds.tax<-metaMDS(diss.tax) #3D : k=3
clust.tax<-scores(mds.tax,display="sites",scaling=1)

# Minus 4 spp
diss.taxmin4<-vegdist(decostand(Year_Geom_Means_rare,"total"),"bray")
mds.taxmin4<-metaMDS(diss.taxmin4) #3D : k=3
clust.taxmin4<-scores(mds.taxmin4,display="sites",scaling=1)

#Calculating relative percentage of the community for each species
Percent_Com<-Year_Geom_Means/rowSums(Year_Geom_Means)


#Figure 2####
ColV2<-c("navy","gold", "mediumaquamarine")


colV<-c(rep(ColV2[1],9),rep(ColV2[2],5),rep(ColV2[3],23))
colVseg<-c(rep(ColV2[1],9),rep(ColV2[2],4),NA,rep(ColV2[3],23))



axis.V<-1.1
label.V<-1.2

#Colour Vector
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)
top4_palette = c("#13ABDA", "#FFD578", "#D17634","#34A576")


pdf("Figures/Fig. 2.pdf", height=6,width=6)
par(mar=c(4,4,1,1),oma=c(2,2,2,2),las=1)
layout(matrix(c(1,2,3,3),2,2, byrow=T))
par(pty='s')
# Taxonomy
ordiplot(mds.tax,type="n",xlab="",ylab="NMDS 2",main="All species", cex.main=1,axes=F,frame=F,cex.lab=label.V, cex.axis=axis.V)
points(clust.tax[,1],clust.tax[,2],pch=19,col= colV,cex=1)
s<-seq(33);s<-s[-length(s)]
segments(clust.tax[s,1],clust.tax[s,2],clust.tax[s+1,1],clust.tax[s+1,2],col=colVseg)
segments(clust.tax[14,1],clust.tax[14,2],clust.tax[15,1],clust.tax[15,2], col=colV[14],lty=2)
text(clust.tax[c(1,10,15,33),1]+c(0,0,0,0),clust.tax[c(1,10,15,33),2]+c(-0.04,0.05,0.04,0.04),rownames(Year_Geom_Means)[c(1,10,15,33)],cex=1,col=colV[c(1,10,15,33)])
axis(1,at=c(-0.3,0,0.3))
axis(2,at=c(-0.3,0,0.3))
legend("bottomright",paste("S = ",round(mds.tax$stress,digits=2)), bty='n', cex=1.2)
legend("topleft", "a", bty='n', cex=1.8,inset=c(-0.15,0.04))


ordiplot(mds.taxmin4,type="n",xlab="",ylab=NA,xlim=c(-0.4,0.4),ylim=c(-0.3,0.45),
         main="Non-commercial species",cex.main=1,axes=F,frame=F,cex.lab=label.V,
         cex.axis=axis.V)
points(clust.taxmin4[,1],clust.taxmin4[,2],pch=19,col=colV,cex=1)
segments(clust.taxmin4[s,1],clust.taxmin4[s,2],clust.taxmin4[s+1,1],
         clust.taxmin4[s+1,2],col=colVseg)
segments(clust.taxmin4[14,1],clust.taxmin4[14,2],
         clust.taxmin4[15,1],clust.taxmin4[15,2], 
         col=colV[14],lty=2)
axis(1,at=c(-0.3,0,0.3))
axis(2,at=c(-0.3,0,0.3))
text(clust.taxmin4[c(1,10,15,33),1]+c(0,0.07,0,-0.03), clust.taxmin4[c(1,10,15,33),2]+c(-0.08,-0.05,0.05,-0.05),rownames(Year_Geom_Means)[c(1,10,15,33)],cex=1,col=colV[c(1,10,15,33)])
legend("bottomright",paste("S = ",round(mds.taxmin4 $stress,digits=2)), bty='n', cex=1.2)
legend("topleft", "b", bty='n', cex=1.8,inset=c(-0.15,0.04))

par(pty='m')
plot(Percent_Com$GADUS_MORHUA*100~c(1981:2013),type='n',lwd=2, xaxs="i",yaxs="i",
     ylab="% of total biomass", xlab="Year", 
     ylim=c(0,100),col=top4_palette[1],cex.lab=label.V,cex.axis=axis.V)
polygon(c(1981,1981:1994,1994),
        c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$SEBASTES_MENTELLA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[1:14]*100,0), col=top4_palette[4])
polygon(c(1981,1981:1994,1994),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[1:14]*100,0), col=top4_palette[3])
polygon(c(1981,1981:1994,1994),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA)[1:14]*100,0), col=top4_palette[2])
polygon(c(1981,1981:1994,1994),c(0,Percent_Com$GADUS_MORHUA[1:14]*100,0), col=top4_palette[1])
polygon(c(1995,1995:2013,2013),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$SEBASTES_MENTELLA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[15:33]*100,0), col=top4_palette[4])
polygon(c(1995,1995:2013,2013),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[15:33]*100,0), col=top4_palette[3])
polygon(c(1995,1995:2013,2013),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA)[15:33]*100,0), col=top4_palette[2])
polygon(c(1995,1995:2013,2013),c(0,Percent_Com$GADUS_MORHUA[15:33]*100,0), col=top4_palette[1])
abline(v=Eras, lwd=1,lty=2, col=8)
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
legend("topleft", "c", bty='n', cex=1.8,inset=c(-0.035,0))
legend("topright",
       legend = c("cod","halibut","plaice","redfish"),
       lwd=2,col=top4_palette,
       box.col="#FF003300",inset=0.005,cex=0.85,ncol=2)
dev.off()


