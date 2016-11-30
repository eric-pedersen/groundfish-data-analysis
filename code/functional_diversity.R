
#Clear workspace
rm(list=ls(all=TRUE))

#Loading data and packages ####
library(plyr)
library(dplyr)
library(bootstrap)
library(RColorBrewer)
library(plotrix)
require(FD)

source("code/functions.R")

DFO_Com = read.csv("data/DFO_Com.csv",stringsAsFactors = F)
DFO_Dataset = read.csv("data/DFO_Dataset.csv",stringsAsFactors = F)
load("data/year_geom_means.Rdata")
top4_sp = c("GADUS_MORHUA","SEBASTES_MENTELLA",
            "REINHARDTIUS_HIPPOGLOSSOIDES", "HIPPOGLOSSOIDES_PLATESSOIDES")

#Calculate Functional Groups####
Trait_Match<-read.csv("data/Heike_Traits.csv",row.names=1)
Traits<-Trait_Match[,c("vertical.position","Food.Items","double.time","length",
                       "trophic.level","Aggregation")]

Traits<-Trait_Match[match(names(DFO_Com),Trait_Match$DFO_clean_2014_name),
                    c("vertical.position","Food.Items","double.time","length",
                      "trophic.level","Aggregation")]

top4_index = match(top4_sp, names(DFO_Com))
row.names(Traits)<-names(DFO_Com)

#total community
Func_Div<-dbFD_batch(Traits,Year_Geom_Means,w.abun=T,stand.x=T,
               calc.FGR=T,calc.FRic=F,
               stand.FRic=F,scale.RaoQ=F,calc.CWM=F,cut_type = "G",cut_val = 11)



#remaining species
Func_Div_rare<-dbFD_batch(Traits[-top4_index,],
                          Year_Geom_Means[,-top4_index],
                          w.abun=T,stand.x=T,calc.FGR=T,
                          calc.FRic=F,stand.FRic=F,scale.RaoQ=F,calc.CWM=F,
                          cut_type = "G",cut_val = 11)

Year = DFO_Dataset$Year
# #jackknife FDis
FDis_jack = data_frame(Year = 1981:2013, FDis = 0,jack.se=0 )
for(y in 1981:2013){
  hold1<-DFO_Com[Year==y,]
  hold3<-matrix(1,nrow(hold1)+1,ncol(DFO_Com))
  colnames(hold3)<-colnames(DFO_Com)
  for(j in 1:nrow(hold1)){
    hold2<-hold1[-j,]
    for(i in 1:ncol(DFO_Com)){
      hold3[j,i]<-CalcZeroInfGeomDens(hold2[,i])
    }}
  A<-dbFD(Traits,hold3,w.abun=T,stand.x=T,calc.FGR=F,calc.FRic=F,
          stand.FRic=F,scale.RaoQ=F,calc.CWM=F)$FDis[-(j+1)]
  FDis_jack[y-1980,2]<-mean(A)
  FDis_jack[y-1977,3]<-sd(A)/sqrt(length(A))}


FDis_jack_rare<-FDis_jack

for(y in 1981:2013){
  hold1<-DFO_Com[,-top4_index][Year==y,]
  hold3<-matrix(1,nrow(hold1)+1,ncol(hold1))
  colnames(hold3)<-colnames(hold1)
  for(j in 1:nrow(hold1)){
    hold2<-hold1[-j,]
    for(i in 1:ncol(hold1)){
      hold3[j,i]<-CalcZeroInfGeomDens(hold2[,i])
    }}
  A<-dbFD(Traits[-top4_index,],hold3,w.abun=T,stand.x=T,calc.FGR=F,
          calc.FRic=F,stand.FRic=F,scale.RaoQ=F,calc.CWM=F)$FDis[-(j+1)]
  FDis_jack_rare[y-1980,2]<-mean(A)
  FDis_jack_rare[y-1980,3]<-sd(A)/sqrt(length(A))}

#Colour Vector
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)

pdf("figures/Fig. 3.pdf",width=6, height=5)
plot(Func_Div$FDis[1:14]~c(1981:1994), type='l', lwd=2, lty=1, 
     ylab="Functional dispersion (FDis)",
     xlab= NA,
     cex.lab=1,cex.axis=1,ylim=c(0.25,0.4),xlim=c(1981,2013))

lines(Func_Div$FDis[15:33]~c(1995:2013), type='l', lwd=2, lty=1)
lines(Func_Div_rare$FDis[1:14]~c(1981:1994), type='l', lwd=2, lty=3)
lines(Func_Div_rare$FDis[15:33]~c(1995:2013), type='l', lwd=2, lty=3)
#lines(decostand(diversity(Year_Geom_Means),method="range")~c(1978:2012),lwd=2,lty=1, col="grey90")
abline(v=Eras, lwd=1,lty=2, col=8)
axis(1,seq(1980,2010,by=5),cex=1,labels=F,tick=T)
legend("bottomright",c("total community","non-commercial species"),
       lty=c(1,3),lwd=2,bty='o',bg=NA,box.col="#FF003300",
       inset=0.005,cex=1.1)
dev.off()
