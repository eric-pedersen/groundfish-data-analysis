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

Ref_bmass_1981<-biomass$Bmass/biomass$Bmass[1]
#Ref_bmass_1981<-Ref_bmass_1981[-1]
Ref_bmass_1981<-(Ref_bmass_1981-min(Ref_bmass_1981))/(Ref_bmass_1981-min(Ref_bmass_1981))[1]*100

Ref_cod_1981<-Year_Geom_Means$GADUS_MORHUA/Year_Geom_Means$GADUS_MORHUA[1]
#Ref_cod_1981 <-Ref_cod_1981[-1]
Ref_cod_1981 <-(Ref_cod_1981-min(Ref_cod_1981))/(Ref_cod_1981-min(Ref_cod_1981))[1]*100


#Figure 4####

axis.V<-1.1
label.V<-1.2

#Colour Vector
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)


pdf("Figures/Fig. 4.pdf",height=5,width=7)
par(las=1, mfrow=c(1,1),pty='m')
plot(Ref_bmass_1981[1:14]~c(1981:1994), type='l', lwd=2, ylab="Scaled similarity to 1981", xlab="Year",pch=19, xlim=c(1980,2014), ylim=c(0,120))
lines(Ref_bmass_1981[15:33]~c(1995:2013), type='l',lwd=2, col=1, pch=19)
lines(Ref_com_1981[1:14]~c(1981:1994), type='l',lwd=2, col="dodgerblue3", pch=19)
lines(Ref_com_1981[15:33]~c(1995:2013), type='l',lwd=2, col="dodgerblue3", pch=19)
lines(Ref_cod_1981[1:14]~c(1981:1994), type='l', lwd=2,col=ColV[1], pch=19)
lines(Ref_cod_1981[15:33]~c(1995:2013), type='l', lwd=2,col=ColV[1], pch=19)
text(x=2014.2,y=Ref_bmass_1981[33],labels=round(Ref_bmass_1981,digits=0)[33])
text(x=2014.2,y= Ref_cod_1981[33],labels=round(Ref_cod_1981,digits=0)[33], col=ColV[1])
text(x=2014.2,y=Ref_com_1981[33],labels=round(Ref_com_1981,digits=0)[33], col= "dodgerblue3")
abline(v=c(1990,1995), lty=2)
legend("topright",legend=c("community biomass","cod biomass","community composition"), lwd=2,col=c(1,ColV[1],"dodgerblue3"
), bty='n')
dev.off()

