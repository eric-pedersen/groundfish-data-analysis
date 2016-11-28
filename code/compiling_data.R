#This file will load each individual data set and clean and compile them
#into intermediate data files for analyses


# This code librarys a clean groundfish community data set added to the data folder. 
# Data are available by request from Fisheries and Oceans Canada (pierre.pepin@dfo-mpo.gc.ca) 

#Calculate Year Means####
#Clear workspace
rm(list=ls(all=TRUE))

source("code/functions.R")

#Packages####
library(plyr)
library(dplyr)
library(sp)




#Load and Preprocess Dataset####
#Main Dataset
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

#Subset community Matrix, selecting only species data
DFO_Com<-DFO_Dataset[,31:ncol(DFO_Dataset)]



#Loads previously calculated Voronoi polygon centroids ####
load("data/voronoi_shapes.Rdat")
voronoi_centroids = as.data.frame(getSpPPolygonsLabptSlots(voronoi_shapes))
names(voronoi_centroids) = c("LONG_DEC","LAT_DEC")
voronoi_centroids$polygon = voronoi_shapes$voroLatt

#removes species that are sensitive to the gear change
voronoi_data = DFO_Dataset[,!names(DFO_Dataset)%in%gear_sensitive_species_names]
voronoi_data$polygon= over(SpatialPoints(cbind(voronoi_data$LONG_DEC, voronoi_data$LAT_DEC)), voronoi_shapes)$voroLatt
voronoi_data = voronoi_data[!is.na(voronoi_data$polygon),]
voronoi_data[,c("LONG_DEC","LAT_DEC")] =voronoi_centroids[match(voronoi_data$polygon, 
                                                                voronoi_centroids$polygon),
                                                          c("LONG_DEC","LAT_DEC")]

# Calculates mean depth of all the trawls in each voronoi polygon                                                         
voronoi_data = voronoi_data %>%
  group_by(polygon) %>%
  mutate(Depth = mean(Depth))%>%
  group_by(polygon, Year) %>%
  mutate(Temp_at_fishing = mean(Temp_at_fishing))%>%
  mutate_each(funs = funs(CalcZeroInfGeomDens),
              ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)%>%
  summarize_each(funs = funs(first), 
                 LAT_DEC:Temp_at_fishing, 
                 DIV, ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)

write.csv(DFO_Com, "data/DFO_Com.csv",row.names = F)
write.csv(DFO_Dataset, "data/DFO_Dataset.csv",row.names = F)
write.csv(voronoi_data, "data/voronoi_data.csv",row.names = F)


