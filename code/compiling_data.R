#This file will load each individual data set and clean and compile them
#into intermediate data files for analyses


# This code librarys a clean groundfish community data set added to the data folder. 
# Data are available by request from Fisheries and Oceans Canada (pierre.pepin@dfo-mpo.gc.ca) 

#Clear workspace
rm(list=ls(all=TRUE))

source("code/functions.R")

#Packages####
library(plyr)
library(dplyr)
library(sp)
library(bootstrap)




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


Year<-factor(DFO_Dataset$Year,ordered=T)

#Calculating yearly geometric means for synchrony and ordination calculations####
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


#Transforming data from  clustering maps
clust_map_data = DFO_Dataset
clust_map_data[,c("polygon","LONG_DEC","LAT_DEC")]= over(SpatialPoints(cbind(clust_map_data$LONG_DEC, clust_map_data$LAT_DEC)), voronoi_shapes)
clust_map_data = clust_map_data %>% 
  group_by(polygon) %>% 
  mutate(Depth=mean(Depth,na.rm = T))%>% 
  ungroup()

year_groups = c("1981-1984","1985-1989","1990-1994","1995-2001","2002-2006", "2007-2013")

#classifies years into 6 sets of years.
clust_map_data$Year_group = with(clust_map_data,
                                 ifelse(Year<1985, year_groups[1],
                                        ifelse(Year<1990, year_groups[2],
                                               ifelse(Year<1995, year_groups[3],
                                                      ifelse(Year<2002, year_groups[4],
                                                             ifelse(Year<2007, year_groups[5],
                                                                    year_groups[6]
                                                             ))))))

clust_map_data = clust_map_data[!rowSums(clust_map_data[,31:60])==0,]


# Calculates geom-mean density for each species in each year group in each polygon
clust_map_data = ddply(clust_map_data, .(polygon, Year_group),function(x){
  out_data = x[1,]
  out_data[1,31:60] = apply(x[,31:60],MARGIN=2, 
                            FUN=CalcZeroInfGeomDens)
  #out_data[1,31:60] = out_data[1,31:60]/sum(out_data[1,31:60])
  return(out_data)
})


write.csv(DFO_Com, "data/DFO_Com.csv",row.names = F)
write.csv(DFO_Dataset, "data/DFO_Dataset.csv",row.names = F)
write.csv(voronoi_data, "data/voronoi_data.csv",row.names = F)
write.csv(clust_map_data, "data/clust_map_data.csv",row.names = F)
save(Year_Geom_Means,Year_Geom_Means_all,Year_Geom_Means_rare,Year_Geom_Means_SE,
     file = "data/year_geom_means.Rdata")


