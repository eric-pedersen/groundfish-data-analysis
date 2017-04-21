#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
library(plyr)
library(dplyr)
library(tidyr)
library(class) #required for k-nearest neighbours matching
library(stringr)
library(mgcv)
library(ggplot2)

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
  labs(y= "conversion factor")+
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
ggsave(filename = "figures/S0 Gearchange conversion estimates.png",
       plot = conv_gearchange_plot, width = 6, height = 8, units="in", dpi=400)
