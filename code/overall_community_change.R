#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
library(plyr)
library(dplyr)
library(tidyr)
library(vegan)
library(RColorBrewer)



#Loading data and functions ####
source("code/functions.R")
load("data/year_geom_means.Rdata")
voronoi_data= read.csv("data/voronoi_data.csv", stringsAsFactors = F)
DFO_Dataset = read.csv("data/DFO_Dataset.csv", stringsAsFactors = F)
Func_Div = read.csv("data/func_div_data.csv",stringsAsFactors = F)

biomass = DFO_Dataset %>%
  group_by(Year) %>%
  dplyr::select(ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)%>%
  mutate(obs = 1:n())%>%
  gather(species, abundance, ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)%>%
  group_by(Year, obs) %>%
  summarize(total = sum(abundance))%>%
  group_by(Year) %>%
  summarize(Bmass = CalcZeroInfGeomDens(total))




polygons = unique(voronoi_data$polygon)
diss.tax<-as.matrix(vegdist(decostand(Year_Geom_Means,"total"),"bray"))
spatial_BC <- matrix(NA, nrow=length(1981:2013), 
                     ncol=length(polygons))
for(i in 1:length(polygons)){
  current_voronoi <- subset(voronoi_data,polygon==polygons[i]&Year>1980)
  current_dist <- vegdist(current_voronoi[,8:37])
  current_dist <- as.matrix(current_dist)
  spatial_BC[current_voronoi$Year-1980,i] <- current_dist[,1]
  rm(current_dist, current_voronoi)
}

Ref_com_1981<-(1-decostand(diss.tax[,1],"range"))*100

Ref_bmass_1981<-biomass$Bmass/biomass$Bmass[1]
#Ref_bmass_1981<-Ref_bmass_1981[-1]
Ref_bmass_1981<-(Ref_bmass_1981-min(Ref_bmass_1981))/(Ref_bmass_1981-min(Ref_bmass_1981))[1]*100

Ref_cod_1981<-Year_Geom_Means$GADUS_MORHUA/Year_Geom_Means$GADUS_MORHUA[1]
#Ref_cod_1981 <-Ref_cod_1981[-1]
Ref_cod_1981 <-(Ref_cod_1981-min(Ref_cod_1981))/(Ref_cod_1981-min(Ref_cod_1981))[1]*100

Ref_fdis_1981 = Func_Div$FDis/Func_Div$FDis[1]
Ref_fdis_1981 = (Ref_fdis_1981 - min(Ref_fdis_1981))/(Ref_fdis_1981-min(Ref_fdis_1981))[1]*100


#Figure 4####

axis.V<-1.1
label.V<-1.2

#Colour Vector
cod_col ="#FF4040"
Eras<-c(1990,1995)


pdf("Figures/Fig. 5.pdf",height=6,width=8)
par(las=1, mfrow=c(1,1),pty='m')
plot(Ref_bmass_1981[1:14]~c(1981:1994), type='l', lwd=2, ylab="Scaled similarity to 1981", 
     xlab="Year",pch=19, xlim=c(1980,2014), ylim=c(0,120))
lines(Ref_bmass_1981[15:33]~c(1995:2013), type='l',lwd=2, col=1, pch=19)
lines(Ref_com_1981[1:14]~c(1981:1994), type='l',lwd=2, col="olivedrab3", pch=19)
lines(Ref_com_1981[15:33]~c(1995:2013), type='l',lwd=2, col="olivedrab3", pch=19)
lines(Ref_cod_1981[1:14]~c(1981:1994), type='l', lwd=2,col=cod_col, pch=19)
lines(Ref_cod_1981[15:33]~c(1995:2013), type='l', lwd=2,col=cod_col, pch=19)
lines(Ref_fdis_1981[1:14]~c(1981:1994), type='l', lwd=2,col="deepskyblue", pch=19)
lines(Ref_fdis_1981[15:33]~c(1995:2013), type='l', lwd=2,col="deepskyblue", pch=19)
text(x=2014.2,y=Ref_bmass_1981[33],labels=round(Ref_bmass_1981,digits=0)[33])
text(x=2014.2,y= Ref_cod_1981[33],labels=round(Ref_cod_1981,digits=0)[33], col=cod_col)
text(x=2014.2,y=Ref_com_1981[33],labels=round(Ref_com_1981,digits=0)[33], col= "olivedrab3")
text(x=2014.2,y=Ref_fdis_1981[33],labels=round(Ref_fdis_1981,digits=0)[33], col= "deepskyblue")
abline(v=c(1990,1995), lty=2)
legend("bottomleft",legend=c("community biomass","cod biomass","community composition", "functional diversity"), 
       lwd=2,col=c(1,cod_col,"olivedrab3", "deepskyblue"), bty='n')
dev.off()

