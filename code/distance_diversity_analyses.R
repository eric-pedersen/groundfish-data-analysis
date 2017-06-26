#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
library(plyr)
library(reshape2)
library(dplyr)
library(tidyr)
library(vegan)
library(bootstrap)
library(RColorBrewer)
library(SoDA)
library(stringr)
library(plotrix)

#Loading data and functions ####
source("code/functions.R")
voronoi_data =read.csv("data/voronoi_data.csv",stringsAsFactors = F)

#Options for diversity-distance plots ####
com_dist = "bray" #uses Bray-Curtis dissimilarity
geodesic_dist = T #Use geodesic distances for calculating inter-polygon distance
boot_errors = T #Should jackknife se values be calculated?
include_variance= T
include_mean = T
min_years_in_group = 30 
#specifies how many years a polygon has to have observations in it for it to be used in the analysis

top4_sp_names = c("GADUS_MORHUA","SEBASTES_MENTELLA",
            "REINHARDTIUS_HIPPOGLOSSOIDES", "HIPPOGLOSSOIDES_PLATESSOIDES")
top4_sp = match(top4_sp_names, names(voronoi_data))
nontop_sp = (8:37)[!(8:37)%in%top4_sp]


group_counts = plyr::count(subset(voronoi_data,Year>1980), .(polygon))
group_counts = group_counts$polygon[group_counts$freq>=min_years_in_group]
voronoi_data = voronoi_data[voronoi_data$polygon%in%group_counts,]


#Calculating spatial, depth, and community distances
mean_dist_data = ddply(voronoi_data, .(Year), 
                       function(x){
                         dist_data = CalcSpaceDepthDists(x[,c("LAT_DEC","LONG_DEC","Depth")],
                                                         use_geodesic=geodesic_dist)
                         dist_data$sp_dist = as.vector(vegdist(decostand(x[,8:37],method="total"),
                                                               method=com_dist))
                         dist_data$sp_nontop_dist = as.vector(vegdist(decostand(x[,nontop_sp],method="total"),
                                                                      method=com_dist))
                         mean_dist_labels = outer(x$polygon,x$polygon,FUN=paste)
                         dist_data$pairs = mean_dist_labels[lower.tri(mean_dist_labels,diag=F)]
                         return(dist_data)
                       })


mean_dist_data = transform(mean_dist_data, space_dist_scl = space_trans(space_dist),
                           depth_dist_scl = depth_trans(depth_dist),
                           sp_dist_scl = comm_trans(sp_dist,com_dist=com_dist),
                           sp_nontop_dist_scl = comm_trans(sp_nontop_dist,com_dist=com_dist)
)

pairs = str_split_fixed(mean_dist_data$pairs,pattern=" ",n=2)
mean_dist_data$site_1  = as.numeric(pairs[,1])
mean_dist_data$site_2  = as.numeric(pairs[,2])
mean_dist_data = melt(mean_dist_data,measure.vars=c("sp_dist_scl","sp_nontop_dist_scl"),
                      variable.name="community_subset",value.name = "distance_scl")

#Calculating regression models 
div_dist_models = ddply(mean_dist_data,.(Year,community_subset), 
                        .fun=CalculateRSquareValues)

mean_com_var = mean_dist_data %>%
  group_by(Year, community_subset) %>% 
  summarize(com_var = var(distance_scl))

if(!include_mean){
  div_dist_models = subset(div_dist_models, coef!="Mean community distance")
}
if(!include_variance){
  div_dist_models = subset(div_dist_models, coef!="variance")
}


# Calculating average diversity ####
mean_diversity_data = voronoi_data %>%
  select(Year,polygon, ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)%>%
  gather(species,abundance,-Year,-polygon)%>%
  group_by(Year,polygon) %>%
  mutate(prop = abundance/sum(abundance),
         abundance_rare = ifelse(species%in%top4_sp_names, 
                                 NA, abundance),
         prop_rare = abundance_rare/sum(abundance_rare,na.rm = T))%>%
  summarise(diversity= sum(-prop*ifelse(prop==0, 0, log(prop))),
            diversity= exp(diversity),
            diversity_rare = sum(-prop_rare*ifelse(prop_rare==0, 0, log(prop_rare)),na.rm = T),
            diversity_rare = exp(diversity_rare))%>%
  group_by(Year)%>%
  summarise(div_sd = sd(diversity,na.rm = T)/sqrt(n()), 
            div_mean = mean(diversity,na.rm = T),
            div_rare_sd = sd(diversity_rare, na.rm=T)/sqrt(n()-4),
            div_rare_mean = mean(diversity_rare, na.rm=T))

#Figure 5####
ComV<-c("sp_dist_scl","sp_nontop_dist_scl")

axis.V<-1.1
label.V<-1.2

#Colour Vector
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)

pdf("Figures/Fig. 6.pdf", height=6,width=6)
par(mfrow=c(2,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)
plot(c(1981:2013),div_dist_models[div_dist_models$community_subset=="sp_dist_scl" & div_dist_models$coef=="Spatial Distance","R2_value"],type='n', ylab=expression(paste("Distance R"^"2")), ylim=c(0,0.3), xlab=NA,xaxt='n', cex.lab=label.V,cex.axis=axis.V, xlim=c(1981,2013))
for(i in 1:3){
  plotCI(c(1981:2013),
         div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","R2_value"],
         uiw=div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","se"],
         pch=NA,slty=i,sfrac=0,add=T)
  lines(c(1981:1994),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","R2_value"][1:14], lty=i, lwd=2)
  lines(c(1995:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","R2_value"][15:33], lty=i, lwd=2)
}
abline(v=Eras, lwd=1,lty=2, col=8)
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
legend("topleft", "A",bty='n', cex=1.8, inset=c(-0.04,-0.025))
legend(x="topright", legend=c("Total community","Noncommercial species"),
       col="black", lty=1:3,lwd=3,bty='n',box.col="#FF003300", cex=0.9)

plot(c(1981:2013),div_dist_models[div_dist_models$community_subset=="sp_dist_scl" & div_dist_models$coef=="Depth","R2_value"],type='n', ylab=expression(paste("Depth R"^"2")), ylim=c(0,0.5), cex.lab=label.V,cex.axis=axis.V, xlim=c(1980,2013),xlab="")
for(i in 1:3){
  plotCI(c(1981:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Depth","R2_value"],uiw=div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Depth","se"],pch=NA,slty=i,sfrac=0,add=T)
  lines(c(1981:1994),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Depth","R2_value"][1:14], lty=i, lwd=2)
  lines(c(1995:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Depth","R2_value"][15:33], lty=i, lwd=2)
}
abline(v=Eras, lwd=1,lty=2, col=8)
mtext("Year",side=1,outer=T,line=1,  adj=0.525,cex=1,padj=0)
legend("topleft", "B",bty='n', cex=1.8, inset=c(-0.03,-0.025))
dev.off()

png("Figures/Fig. S4.png", height=9,width=6,units="in",res = 800 )
par(mfrow=c(3,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)
plot(c(1981:2013),mean_diversity_data$div_mean,type='n',xaxt='n',ylim = c(2.5,4.75), 
     ylab="Average diversity\n(effective number of species)", 
     cex.lab=label.V,cex.axis=axis.V,xlab="", xlim=c(1980,2013))
plotCI(c(1981:2013),mean_diversity_data$div_mean,
       uiw=mean_diversity_data$div_sd,pch=NA,slty=1,sfrac=0,add=T)
lines(c(1981:1994),mean_diversity_data$div_mean[1:14], lty=1, lwd=2)
lines(c(1995:2013),mean_diversity_data$div_mean[15:33], lty=1, lwd=2)

plotCI(c(1981:2013),mean_diversity_data$div_rare_mean,
       uiw=mean_diversity_data$div_rare_sd,pch=NA,slty=2,sfrac=0,add=T)
lines(c(1981:1994),mean_diversity_data$div_rare_mean[1:14], lty=2, lwd=2)
lines(c(1995:2013),mean_diversity_data$div_rare_mean[15:33], lty=2, lwd=2)
legend("topleft", "(a)",bty='n', cex=1.8, inset=c(-0.03,-0.025))
abline(v=Eras, lwd=1,lty=2, col=8)


plot(c(1981:2013),div_dist_models[div_dist_models$community_subset=="sp_dist_scl" & div_dist_models$coef=="Mean community distance","R2_value"],type='n',ylim = c(0.3,2.4), ylab="mean dissimilarity", cex.lab=label.V,cex.axis=axis.V, xlim=c(1980,2013),xaxt='n',xlab="")
for(i in 1:3){
  plotCI(c(1981:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Mean community distance","R2_value"],uiw=div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Mean community distance","se"],pch=NA,slty=i,sfrac=0,add=T)
  lines(c(1981:1994),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Mean community distance","R2_value"][1:14], lty=i, lwd=2)
  lines(c(1995:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Mean community distance","R2_value"][15:33], lty=i, lwd=2)
}
abline(v=Eras, lwd=1,lty=2, col=8)
mtext("Year",side=1,outer=T,line=1,  adj=0.525,cex=1,padj=0)
legend("topleft", "(b)",bty='n', cex=1.8, inset=c(-0.03,-0.025))


plot(c(1981:2013),div_dist_models[div_dist_models$community_subset=="sp_dist_scl" & div_dist_models$coef=="variance","R2_value"],type='n',ylim = c(0.9,4.55), ylab="total variance", cex.lab=label.V,cex.axis=axis.V,xlab="", xlim=c(1980,2013))
for(i in 1:3){
  plotCI(c(1981:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="variance","R2_value"],uiw=div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="variance","se"],pch=NA,slty=i,sfrac=0,add=T)
  lines(c(1981:1994),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="variance","R2_value"][1:14], lty=i, lwd=2)
  lines(c(1995:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="variance","R2_value"][15:33], lty=i, lwd=2)
}
abline(v=Eras, lwd=1,lty=2, col=8)
mtext("Year",side=1,outer=T,line=1,  adj=0.525,cex=1,padj=0)
legend("topleft", "(c)",bty='n', cex=1.8, inset=c(-0.03,-0.025))
legend(x="topright", legend=c("Total community","Noncommercial species"),
       col="black", lty=1:3,lwd=3,bty='n',box.col="#FF003300", cex=0.9)

dev.off()

