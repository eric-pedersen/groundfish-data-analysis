#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
library(plyr)
library(dplyr)
library(synchrony)
library(bootstrap)
library(plotrix)
library(RColorBrewer)

#Loading data and functions ####
source("code/functions.R")
DFO_Com = read.csv("data/DFO_Com.csv",stringsAsFactors = F)
DFO_Dataset = read.csv("data/DFO_Dataset.csv",stringsAsFactors = F)

Div<-DFO_Dataset$DIV
Year<-factor(DFO_Dataset$Year,ordered=T)

top4_sp = c("GADUS_MORHUA","SEBASTES_MENTELLA",
            "REINHARDTIUS_HIPPOGLOSSOIDES", "HIPPOGLOSSOIDES_PLATESSOIDES")


#Total Biomass####

Total_Biomass<- ddply(DFO_Com,.variables=.(Year),
                      .fun=function(x){
                        time_series = rowSums(x)
                        jack<-jackknife(time_series,CalcZeroInfGeomDens)
                        return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
                      })

Total_Biomass_rare<- ddply(DFO_Com[,!names(DFO_Com)%in%top4_sp],
                           .variables=.(Year),
                           .fun=function(x){
                             time_series = rowSums(x)
                             jack<-jackknife(time_series,CalcZeroInfGeomDens)
                             return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
                           })


#Calculating yearly geometric means for synchrony calculations
Year_Geom_Means<-data.frame(matrix(NA,length(unique(Year)),
                                   ncol(DFO_Com),
                                   dimnames=list(levels(Year),
                                                 colnames(DFO_Com))))
Year_Geom_Means_SE<-Year_Geom_Means
for(i in 1:ncol(DFO_Com)){
  hold<- ddply(DFO_Com,.variables=.(Year),
               .fun=function(x){
                 time_series = x[,i]
                 jack<-jackknife(time_series,CalcZeroInfGeomDens)
                 return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
               })
  Year_Geom_Means[,i]<-hold$Bmass
  Year_Geom_Means_SE[,i]<-hold$jack.se
}


#mean density and subsets of mean densities across years
Year_Geom_Means_all<-Year_Geom_Means
Year_Geom_Means_rare<-Year_Geom_Means[,!names(Year_Geom_Means)%in% top4_sp]


#Synchrony calculations  ####
twindow = 5
sync_allsp = community.sync.window(Year_Geom_Means[rownames(Year_Geom_Means)<1995,], twin=twindow)
sync_allsp2 = community.sync.window(Year_Geom_Means[rownames(Year_Geom_Means)>1994,], twin=twindow)
sync_rare=community.sync.window(Year_Geom_Means_rare[rownames(Year_Geom_Means)<1995,], twin=twindow)
sync_rare2=community.sync.window(Year_Geom_Means_rare[rownames(Year_Geom_Means)>1994,], twin=twindow)

alpha=0.05 





#Figure 1#####

axis.V<-1.1
label.V<-1.2

#Colour Vector
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)

pdf("figures/Fig. 1 - 5yr.pdf", height=7,width=7)
par(mfrow=c(3,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)
plot(Total_Biomass$Bmass[1:14]~c(1981:1994),type='l',lty=1, lwd=2, ylab="Biomass (kg/tow)", xlab=NA,ylim=c(0,210),cex.lab=label.V,cex.axis=axis.V,xaxt='n',xlim=c(1981,2013))
lines(Total_Biomass$Bmass[15:33]~c(1995:2013),lty=1, lwd=2)
plotCI(c(1981:2013),Total_Biomass$Bmass,uiw=Total_Biomass$jack.se,add=T,pch=NA,sfrac=0)
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
plotCI(c(1981:2013),Total_Biomass_rare$Bmass,uiw=Total_Biomass_rare$jack.se,add=T,pch=NA,sfrac=0,lty=3)
lines(c(1981:1994),Total_Biomass_rare$Bmass[1:14],type='l',lty=3, lwd=2, ylim=c(0,1))
lines(c(1995:2013),Total_Biomass_rare$Bmass[15:33],type='l',lty=3, lwd=2, ylim=c(0,1))
abline(v=Eras, lwd=1,lty=2, col=8)
legend("topright",c("total community","remaining species"),lty=c(1,3),lwd=2,bty='o',bg="white",box.col="#FF003300",inset=0.005,cex=1.1)
legend("topleft", "a",bty='n', cex=1.8, inset=c(-0.04,-0.025))



plot(Year_Geom_Means$GADUS_MORHUA[1:14]~c(1981:1994),type='l',lwd=2, ylab="Biomass (kg/tow)", xlab=NA,col=ColV[1],cex.lab=label.V,cex.axis=axis.V,xaxt='n',ylim=c(0,33),xlim=c(1981,2013))
lines(Year_Geom_Means$GADUS_MORHUA[15:33]~c(1995:2013),type='l',lwd=2,col=ColV[1])
plotCI(c(1981:2013),Year_Geom_Means$GADUS_MORHUA,uiw=Year_Geom_Means_SE$GADUS_MORHUA,add=T, pch=NA,sfrac=0, col=ColV[1])
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
lines(Year_Geom_Means$REINHARDTIUS_HIPPOGLOSSOIDES[1:14]~c(1981:1994),col=ColV[2],lwd=2)
lines(Year_Geom_Means$REINHARDTIUS_HIPPOGLOSSOIDES[15:33]~c(1995:2013),col=ColV[2],lwd=2)
plotCI(c(1981:2013),Year_Geom_Means$REINHARDTIUS_HIPPOGLOSSOIDES,uiw=Year_Geom_Means_SE$REINHARDTIUS_HIPPOGLOSSOIDES,add=T, pch=NA,sfrac=0, col=ColV[2])
lines(Year_Geom_Means$HIPPOGLOSSOIDES_PLATESSOIDES[1:14]~c(1981:1994),col=ColV[3],lwd=2)
lines(Year_Geom_Means$HIPPOGLOSSOIDES_PLATESSOIDES[15:33]~c(1995:2013),col=ColV[3],lwd=2)
plotCI(c(1981:2013),Year_Geom_Means$HIPPOGLOSSOIDES_PLATESSOIDES,uiw=Year_Geom_Means_SE$HIPPOGLOSSOIDES_PLATESSOIDES,add=T, pch=NA,sfrac=0, col=ColV[3])
lines(Year_Geom_Means$SEBASTES_MENTELLA[1:14]~c(1981:1994),col=ColV[4],lwd=2)
lines(Year_Geom_Means$SEBASTES_MENTELLA[15:33]~c(1995:2013),col=ColV[4],lwd=2)
plotCI(c(1981:2013),Year_Geom_Means$SEBASTES_MENTELLA,uiw=Year_Geom_Means_SE$SEBASTES_MENTELLA,add=T, pch=NA,sfrac=0, col=ColV[4])
abline(v=Eras, lwd=1,lty=2, col=8)
legend("topright",legend = c("cod","halibut","plaice","redfish"),lwd=2,col=c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),box.col="#FF003300",inset=0.005,cex=1.1)
legend("topleft", "b",bty='n', cex=1.8, inset=c(-0.04,-0.025))

phase<-cbind(Total_Biomass$Bmass[5:33]-Total_Biomass$Bmass[1:29],1985:2013)
grad.V<-rev(brewer.pal(9,"RdBu"))
grad.V[5]<-"lightgrey"

plot(sync_allsp[,"end"], sync_allsp[,"sync"], t="n",
     xlab="", ylab="Community synchrony", ylim=c(0, 1), lwd = 4, cex.lab= label.V, cex.axis=axis.V, xlim=c(1981,2013))
#polygon(x = c(1975,1975,2015,2015),y=c(-1,2,2,-1),col="grey")
segments(x0 =sync_allsp[-nrow(sync_allsp),"end"],y0 =sync_allsp[-nrow(sync_allsp),"sync"],x1 =sync_allsp[-1,"end"],y1 =sync_allsp[-1,"sync"],col= grad.V[(((phase[,1]/max(abs(phase[,1])))+1)*5)+1],lwd=3)
segments(x0 =sync_allsp2[-nrow(sync_allsp2),"end"],y0 =sync_allsp2[-nrow(sync_allsp2),"sync"],x1 =sync_allsp2[-1,"end"],y1 =sync_allsp2[-1,"sync"],col= grad.V[(((phase[11:29,1]/max(abs(phase)))+1)*5)+1],lwd=3)
points(sync_allsp[,"end"], sync_allsp[,"sync"], col="black", pch=21, 
       bg=ifelse(sync_allsp[,"pval"] < alpha, 'black', 'white'))
points(sync_allsp2[,"end"], sync_allsp2[,"sync"], col="black", pch=21, 
       bg=ifelse(sync_allsp2[,"pval"] < alpha, 'black', 'white'))
abline(v=Eras, lwd=1,lty=2, col=8)


phase<-cbind(Total_Biomass_rare$Bmass[5:33]-Total_Biomass_rare $Bmass[1:29],1985:2013)
segments(x0 = sync_rare[-nrow(sync_rare),"end"],y0 = sync_rare[-nrow(sync_rare),"sync"],x1 = sync_rare[-1,"end"],y1 = sync_rare[-1,"sync"],col= grad.V[(((phase[,1]/max(abs(phase[,1])))+1)*5)+1],lwd=3, lty=3)
segments(x0 = sync_rare2[-nrow(sync_rare2),"end"],y0 = sync_rare2[-nrow(sync_rare2),"sync"],x1 = sync_rare2[-1,"end"],y1 = sync_rare2[-1,"sync"],col= grad.V[(((phase[11:29,1]/max(abs(phase)))+1)*5)+1],lwd=3,lty=3)
points(sync_rare[,"end"], sync_rare[,"sync"], col="black", pch=21, 
       bg=ifelse(sync_rare[,"pval"] < alpha, 'black', 'white'))
points(sync_rare2[,"end"], sync_rare2[,"sync"], col="black", pch=21, 
       bg=ifelse(sync_rare2[,"pval"] < alpha, 'black', 'white'))
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
color.legend(1980.5,0.1,1982,0.8,legend=c("-","0","+"),rect.col= grad.V,gradient = "y",align="rb")
mtext("Year",side=1,outer=T,line=1,  adj=0.525,cex=1)
legend("topleft", "c",bty='n', cex=1.8, inset=c(-0.04,-0.025))
dev.off()