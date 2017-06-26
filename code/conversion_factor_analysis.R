#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
library(plyr)
library(FD)
library(dplyr)
library(tidyr)
library(class) #required for k-nearest neighbours matching
library(stringr)
library(mgcv)
library(ggplot2)
library(bootstrap)
library(RColorBrewer)
library(plotrix)
library(vegan)

source("code/functions.R")
#Loading data and functions ####
DFO_Dataset_Init<-read.csv("data/DFO_SURVEYS_text_cleaner.csv")#DFO Dataset

#Add 2013 data to the older dataset
DFO_Dataset_2013<-read.csv("data/DFO_RV_Surveys_kg-per-tow_fish+shrimp+crab_Spring_Fall_2013.csv")
DFO_Dataset<-rbind.fill(DFO_Dataset_Init,DFO_Dataset_2013)[,1:ncol(DFO_Dataset_Init)]

#remove inverts
DFO_Dataset <- DFO_Dataset[,colnames(DFO_Dataset) !="CHIONOECETES_OPILIO"]
DFO_Dataset <- DFO_Dataset[,colnames(DFO_Dataset) !="PANDALUS_BOREALIS"]
DFO_Dataset <- DFO_Dataset[,colnames(DFO_Dataset) !="PANDALUS_MONTAGUI"]
DFO_Dataset <- DFO_Dataset[,colnames(DFO_Dataset) !="PANDALUS_PROPINQUUS"]


#remove inshore samples
DFO_Dataset<-DFO_Dataset[DFO_Dataset$Strata_Type!="InshoreNew",]
DFO_Dataset<-DFO_Dataset[DFO_Dataset$Strata_Type!="InshoreOld",]
DFO_Dataset$Strata_Type<-droplevels(DFO_Dataset$Strata_Type)

#remove spring samples
DFO_Dataset<-DFO_Dataset[DFO_Dataset$Season!="Spring",]
DFO_Dataset$Season<-droplevels(DFO_Dataset$Season)

#remove 3N and 3O
DFO_Dataset<-DFO_Dataset[DFO_Dataset$DIV!="3N",]
DFO_Dataset<-DFO_Dataset[DFO_Dataset$DIV!="3O",]
DFO_Dataset<-DFO_Dataset[DFO_Dataset$DIV!="2H",]
DFO_Dataset$DIV<-droplevels(DFO_Dataset$DIV)

#Set all NA values for species abundances equal to zero
DFO_Dataset[,31:ncol(DFO_Dataset)][is.na(DFO_Dataset[,31:ncol(DFO_Dataset)])] = 0

#remove pre 1981 data####
DFO_Dataset<-DFO_Dataset[DFO_Dataset$Year>=1981,]


#List of column numbers of the comunity matrix, and names of species identified
#as gear-sensitive.
gear_sensitive_species<-c(1,2,3,7,8,9,10,11,12,14,15,16,20,22,25,28,30,31,32,33,
                          34,39,40,42,47,54,55,56,57)
gear_sensitive_species_names = names(DFO_Dataset[,31:89])[gear_sensitive_species]

# Remove gear-sensitive species from the data
DFO_Dataset = DFO_Dataset[,!names(DFO_Dataset)%in% gear_sensitive_species_names]


# this function takes a specific data frame and an id vector specificing which
# points should be in the pre- or post- data set. It then matches data points
# in the post data set with those in the pre data set, finding the nearest neighbour
# based on scaled latitude, longitude, and log-transformed depth. 
# It then transforms both data sets into long form (with one column for species 
# and one for density) then joins the two data sets based on the match.
find_matched_species = function(full_data, test_data_id){
  stopifnot(is.logical(test_data_id)) 
  stopifnot(length(test_data_id)==nrow(full_data))
  test_data = filter(full_data, test_data_id>0)
  train_data = filter(full_data, test_data_id==0)
  train_data$id = 1:nrow(train_data)
  
  test_data$id = knn(select(train_data,lat_scale,long_scale,depth_scale),
                     select(test_data,lat_scale, long_scale,depth_scale),
                     train_data$id)
  test_data$id = as.numeric(test_data$id)
  
  
  train_data_long = gather(train_data, species,train_density, 
                           ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)
  
  test_data_long = gather(test_data,species,test_density, 
                          ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)
  test_data_long = select(test_data_long, -lat_scale, -long_scale, -depth_scale,-Year)
  test_data_long = left_join(test_data_long, train_data_long)
  
  # creates a new conversion factor based on the ratio of the before and after
  # densities in matched observations
  test_data_long = mutate(test_data_long, 
                         conv_factor = log((test_density+0.01)/(train_density+0.01)),
                         species =factor(species)
  )
  return(test_data_long)
}

# scaling latitude, longitude, and depth in the DFO biomass dataset
# and filtering only for those years immediately before and after the gear change
gearchange_data = DFO_Dataset %>%
  mutate(lat_scale = as.vector(scale(LAT_DEC,center = T,scale=T)),
         long_scale = as.vector(scale(LONG_DEC, center=T,scale = T)),
         depth_scale = as.vector(scale(log10(Depth), center=T,scale=T)))%>%
  filter(Year %in%(1993:1996)) %>% 
  select(Year,lat_scale,long_scale,depth_scale, ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)

prechange_data = filter(gearchange_data,Year<1995)
postchange_data = filter(gearchange_data, Year>=1995)

# Creating three data sets: one for points matched before and after
# the gear change, one for points immediately before the gear change and one for points
# immediately after the gear change
gearchange_conv = find_matched_species(gearchange_data, gearchange_data$Year>1994)
prechange_conv = find_matched_species(prechange_data, prechange_data$Year==1994)
postchange_conv = find_matched_species(postchange_data, postchange_data$Year==1996)                                      

gearchange_conv = gearchange_conv %>%
  filter(test_density>0|train_density>0)

# Mixed effect models of the gear change and the before and after changes
conv_gearchange_model = gam(conv_factor~s(species,bs="re"),data= gearchange_conv)
conv_prechange_model = gam(conv_factor~s(species,bs="re"),data= prechange_conv)
conv_postchange_model = gam(conv_factor~s(species,bs="re"),data= postchange_conv)

conv_fit_data = data.frame(species= unique(gearchange_conv$species))
conv_fit_gearchange = predict(conv_gearchange_model, conv_fit_data,se.fit=T)
conv_fit_prechange  = predict(conv_prechange_model, conv_fit_data, se.fit=T)
conv_fit_postchange  = predict(conv_postchange_model, conv_fit_data, se.fit=T)

conv_fit_data$gearchange_fit = as.vector(conv_fit_gearchange$fit)
conv_fit_data$gearchange_se = as.vector(conv_fit_gearchange$se.fit)

conv_fit_data$prechange_fit = as.vector(conv_fit_prechange$fit)
conv_fit_data$prechange_se = as.vector(conv_fit_prechange$se.fit)

conv_fit_data$postchange_fit = as.vector(conv_fit_postchange$fit)
conv_fit_data$postchange_se = as.vector(conv_fit_postchange$se.fit)

# This adjusts the gear change by the average of the rates of range before and
# after the gear change
conv_fit_data = conv_fit_data %>%
  mutate(total_fit = gearchange_fit - (prechange_fit+postchange_fit)/2,
         total_se = sqrt(gearchange_se^2 + (1/2^2)*(prechange_se^2 +postchange_se^2)))

# This arranges the data by corrected conversion factor
conv_fit_data = conv_fit_data %>%
  mutate(species_label = str_to_title(species),
         species_label = str_replace(species_label,"_", " "))%>%
  arrange(gearchange_fit)%>%
  mutate(species_label = factor(species_label, levels = species_label),
         top4 = ifelse(species_label%in% c("Sebastes mentella", "Gadus morhua",
                                     "Reinhardtius hippoglossoides",
                                     "Hippoglossoides platessoides"), 1, 0))

conv_gearchange_plot = ggplot(aes(x=species_label, y=exp(gearchange_fit),
                                  color= factor(top4)), data=conv_fit_data)+
  geom_point()+
  scale_color_manual(values = c("black","red"))+
  geom_hline(yintercept = 1)+
  geom_linerange(aes(ymin = exp(gearchange_fit-2*gearchange_se), 
                     ymax= exp(gearchange_fit+2*gearchange_se)))+
  coord_flip()+
  labs(y= "Estimated conversion factor (q)")+
  theme_bw()+
  theme(axis.text.y = element_text(face="italic"), legend.position = "none",
        axis.title.y = element_blank())+
  scale_y_log10()


conv_prechange_plot = ggplot(aes(x=species, y=exp(prechange_fit)), data=conv_fit_data)+
  geom_point()+
  geom_hline(yintercept = 1)+
  geom_linerange(aes(ymin = exp(prechange_fit-2*prechange_se), 
                     ymax= exp(prechange_fit+2*prechange_se)))+
  coord_flip()+
  labs(y= "conversion factor")+
  theme_bw() 



conv_postchange_plot = ggplot(aes(x=species, y=exp(postchange_fit)), data=conv_fit_data)+
  geom_point()+
  geom_hline(yintercept = 1)+
  geom_linerange(aes(ymin = exp(postchange_fit-2*postchange_se), 
                     ymax= exp(postchange_fit+2*postchange_se)))+
  coord_flip()+
  labs(y= "conversion factor")+
  theme_bw() 

conv_agg_plot = ggplot(aes(x=species, y=exp(total_fit)),
                       data=conv_fit_data)+
  geom_point()+
  geom_hline(yintercept = 1)+
  geom_linerange(aes(ymin = exp(total_fit-2*total_se), 
                     ymax= exp(total_fit+2*total_se)))+
  coord_flip()+
  labs(y= "conversion factor")+
  theme_bw() 

conversion_factors = conv_fit_data %>%
  select(species, gearchange_fit)%>%
  mutate(conversion = 1/exp(gearchange_fit))%>%
  select(-gearchange_fit)



write.csv(conversion_factors ,"data/conversion_factors.csv",row.names = F)
ggsave(filename = "figures/fig. S5.png",
       plot = conv_gearchange_plot, width = 6, height = 8, units="in", dpi=400)

load("data/year_geom_means.Rdata")

Div<-DFO_Dataset$DIV
Year<-factor(DFO_Dataset$Year,ordered=T)
DFO_Dataset_conv = DFO_Dataset
Year_Geom_Means_Conv = Year_Geom_Means
Year_Geom_Means_rare_Conv = Year_Geom_Means_rare
top4_sp = c("GADUS_MORHUA","SEBASTES_MENTELLA",
            "REINHARDTIUS_HIPPOGLOSSOIDES", "HIPPOGLOSSOIDES_PLATESSOIDES")


# re-creating figures with conversion factors
for(i in conversion_factors$species){
  if(i %in% names(DFO_Dataset_conv)){
    convert =  conversion_factors$conversion[conversion_factors$species==i]
    DFO_Dataset_conv[DFO_Dataset_conv$Year>1994, names(DFO_Dataset_conv)==i] = convert*DFO_Dataset_conv[DFO_Dataset_conv$Year>1994, names(DFO_Dataset_conv)==i]
    Year_Geom_Means_Conv[(1981:2013)>1994, names(Year_Geom_Means_Conv)==i] = convert*Year_Geom_Means_Conv[(1981:2013)>1994, names(Year_Geom_Means_Conv)==i] 
    if(!i%in% top4_sp){
      Year_Geom_Means_rare_Conv[(1981:2013)>1994, names(Year_Geom_Means_rare_Conv)==i] = convert*Year_Geom_Means_rare_Conv[(1981:2013)>1994, names(Year_Geom_Means_rare_Conv)==i] 
    }
  }
}

DFO_Com_conv = DFO_Dataset_conv %>%
  select(ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)



#Calculating overall biomass taking conversion factors into account ####
Total_Biomass<- ddply(DFO_Com_conv,.variables=.(Year),
                      .fun=function(x){
                        time_series = rowSums(x)
                        jack<-jackknife(time_series,CalcZeroInfGeomDens)
                        return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
                      })


Total_Biomass_rare<- ddply(DFO_Com_conv[,!names(DFO_Com_conv)%in%top4_sp],
                           .variables=.(Year),
                           .fun=function(x){
                             time_series = rowSums(x)
                             jack<-jackknife(time_series,CalcZeroInfGeomDens)
                             return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
                           })



# Calculating new NMDS scores for all and rare taxa with conversion factors ####
#NMDS Calculations####
#taxonomy
# Calculate Bray Curtis dissimilarity, perform NMDS, record site scores
# All spp
diss.tax<-vegdist(decostand(Year_Geom_Means_Conv,"total"),"bray")
mds.tax<-metaMDS(diss.tax) #3D : k=3
clust.tax<-scores(mds.tax,display="sites",scaling=1)

# Minus 4 spp
diss.taxmin4<-vegdist(decostand(Year_Geom_Means_rare_Conv,"total"),"bray")
mds.taxmin4<-metaMDS(diss.taxmin4) #3D : k=3
clust.taxmin4<-scores(mds.taxmin4,display="sites",scaling=1)


#Functional diversity analysis
Trait_Match<-read.csv("data/Heike_Traits.csv",row.names=1)
Traits<-Trait_Match[,c("vertical.position","Food.Items","double.time","length",
                       "trophic.level","Aggregation")]

Traits<-Trait_Match[match(names(DFO_Com_conv),Trait_Match$DFO_clean_2014_name),
                    c("vertical.position","Food.Items","double.time","length",
                      "trophic.level","Aggregation")]

top4_index = match(top4_sp, names(DFO_Com_conv))
row.names(Traits)<-names(DFO_Com_conv)


# Calculate functional diversities ####
#total community
Func_Div<-dbFD_batch(Traits,Year_Geom_Means_Conv,w.abun=T,stand.x=T,
                     calc.FGR=T,calc.FRic=F,
                     stand.FRic=F,scale.RaoQ=F,calc.CWM=F,cut_type = "G",cut_val = 11)



#remaining species
Func_Div_rare<-dbFD_batch(Traits[-top4_index,],
                          Year_Geom_Means_Conv[,-top4_index],
                          w.abun=T,stand.x=T,calc.FGR=T,
                          calc.FRic=F,stand.FRic=F,scale.RaoQ=F,calc.CWM=F,
                          cut_type = "G",cut_val = 11)


#Figure 1#####

axis.V<-1.1
label.V<-1.2

#Colour Vector
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)
top4_palette = c("#13ABDA", "#FFD578", "#D17634","#34A576")

#Figure 2####
ColV2<-c("navy","gold", "mediumaquamarine")
colV<-c(rep(ColV2[1],9),rep(ColV2[2],5),rep(ColV2[3],23))
colVseg<-c(rep(ColV2[1],9),rep(ColV2[2],4),NA,rep(ColV2[3],23))


png("figures/fig. S6.png", height=8,width=12,units="in", res=400)
par(mar=c(4,4,1,1),oma=c(2,2,2,2),las=1)
layout(matrix(c(1,1,4,4,1,1,4,4,2,2,4,4,2,2,5,5,3,3,5,5,3,3,5,5),6,4, byrow=T))
plot(Total_Biomass$Bmass[1:14]~c(1981:1994),type='l',lty=1, lwd=2, ylab="Biomass (kg/tow)", xlab=NA,ylim=c(0,210),cex.lab=label.V,cex.axis=axis.V,xaxt='n',xlim=c(1981,2013))
lines(Total_Biomass$Bmass[15:33]~c(1995:2013),lty=1, lwd=2)
plotCI(c(1981:2013),Total_Biomass$Bmass,uiw=Total_Biomass$jack.se,add=T,pch=NA,sfrac=0)
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
plotCI(c(1981:2013),Total_Biomass_rare$Bmass,uiw=Total_Biomass_rare$jack.se,add=T,pch=NA,sfrac=0,lty=3)
lines(c(1981:1994),Total_Biomass_rare$Bmass[1:14],type='l',lty=3, lwd=2, ylim=c(0,1))
lines(c(1995:2013),Total_Biomass_rare$Bmass[15:33],type='l',lty=3, lwd=2, ylim=c(0,1))
abline(v=Eras, lwd=1,lty=2, col=8)
legend("top",c("total community","non-commercial species"),
       lty=c(1,3),lwd=2,bty='o',bg="white",box.col="#FF003300",inset=0.005,cex=1.1)
legend("topleft", "(a)",bty='n', cex=1.8, inset=c(-0.04,-0.025))



plot(Year_Geom_Means_Conv$GADUS_MORHUA[1:14]~c(1981:1994),type='l',lwd=2, ylab="Biomass (kg/tow)", xlab=NA,col=top4_palette[1],cex.lab=label.V,cex.axis=axis.V,xaxt='n',ylim=c(0,33),xlim=c(1981,2013))
lines(Year_Geom_Means_Conv$GADUS_MORHUA[15:33]~c(1995:2013),type='l',lwd=2,col=top4_palette[1])
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
lines(Year_Geom_Means_Conv$REINHARDTIUS_HIPPOGLOSSOIDES[1:14]~c(1981:1994),col=top4_palette[2],lwd=2)
lines(Year_Geom_Means_Conv$REINHARDTIUS_HIPPOGLOSSOIDES[15:33]~c(1995:2013),col=top4_palette[2],lwd=2)
lines(Year_Geom_Means_Conv$HIPPOGLOSSOIDES_PLATESSOIDES[1:14]~c(1981:1994),col=top4_palette[3],lwd=2)
lines(Year_Geom_Means_Conv$HIPPOGLOSSOIDES_PLATESSOIDES[15:33]~c(1995:2013),col=top4_palette[3],lwd=2)
lines(Year_Geom_Means_Conv$SEBASTES_MENTELLA[1:14]~c(1981:1994),col=top4_palette[4],lwd=2)
lines(Year_Geom_Means_Conv$SEBASTES_MENTELLA[15:33]~c(1995:2013),col=top4_palette[4],lwd=2)
abline(v=Eras, lwd=1,lty=2, col=8)
legend("topright",legend = c("cod","halibut","plaice","redfish"),lwd=2,col=top4_palette ,box.col="#FF003300",inset=0.005,cex=1.1)
legend("topleft", "(b)",bty='n', cex=1.8, inset=c(-0.04,-0.025))


plot(Func_Div$FDis[1:14]~c(1981:1994),type='l',lty=1, lwd=2, ylab="Functional diversity", 
     xlab="Year",ylim=c(0.25,0.4),cex.lab=label.V,cex.axis=axis.V,xlim=c(1981,2013))
lines(Func_Div$FDis[15:33]~c(1995:2013),lty=1, lwd=2)
lines(c(1981:1994),Func_Div_rare$FDis[1:14],type='l',lty=3, lwd=2, ylim=c(0,1))
lines(c(1995:2013),Func_Div_rare$FDis[15:33],type='l',lty=3, lwd=2, ylim=c(0,1))
abline(v=Eras, lwd=1,lty=2, col=8)
legend("topleft", "(c)",bty='n', cex=1.8, inset=c(-0.04,-0.025))



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
legend("topleft", "(d)", bty='n', cex=1.8,inset=c(0,0.04))


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
legend("topleft", "(e)", bty='n', cex=1.8,inset=c(0,0.04))


dev.off()