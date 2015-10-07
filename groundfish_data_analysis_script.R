#Figs for Groundfish Working Group 
#Paper 1
#March 2014

#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
require(plyr)
require(dplyr)
require(FD)
require(picante)
require(plotrix)
require(RColorBrewer)
require(scatterplot3d)
require(synchrony)
require(foreign)
require(vegan)
require(SoDA)
require(stringr)
require(reshape2)
require(bootstrap)
require(ggplot2)
require(mgcv)

require(sp)
require(grid)
require(rgdal)
require(maptools)
require(rgeos)
require(ggdendro)
require(mapproj)

#Load and Preprocess Dataset####
#Main Dataset
DFO_Dataset_Init<-read.csv("~/Dropbox/Groundfish working group shared files/data/DFO_SURVEYS_text_cleaner.csv")#DFO Dataset
DFO_Dataset_2013<-read.csv("~/Dropbox/Groundfish working group shared files/data/DFO_RV_Surveys_kg-per-tow_fish+shrimp+crab_Spring_Fall_2013.csv")

DFO_Dataset<-rbind.fill(DFO_Dataset_Init,DFO_Dataset_2013)[,1:ncol(DFO_Dataset_Init)]



#Colour Vector
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)

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

#Set all NA values equal to zero
DFO_Dataset[,31:ncol(DFO_Dataset)][is.na(DFO_Dataset[,31:ncol(DFO_Dataset)])] = 0

#remove pre 1981 data####
DFO_Dataset<-DFO_Dataset[DFO_Dataset$Year>=1981,]

#Subset community Matrix
DFO_Com<-DFO_Dataset[,31:ncol(DFO_Dataset)]
DFO_Com[is.na(DFO_Com)]<-0

Div<-DFO_Dataset$DIV
Year<-factor(DFO_Dataset$Year,ordered=T)

#Calculate Year Means####
#function to calculate geometric mean abundance
CalcZeroInfGeomDens = function(x){
  prob_obs= sum(x>0)/length(x)
  geom_mean_dens = ifelse(prob_obs>0,exp(mean(log(x[x>0]))),0)
  return(geom_mean_dens*prob_obs)
}

#geometric mean yearly biomass for each species
Year_Geom_Means<-data.frame(matrix(NA,length(unique(Year)),ncol(DFO_Com),dimnames=list(levels(Year),colnames(DFO_Com))))
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


# Plots species by the ratio of densities before and after the gear change
# to assess sensitivity to the gear change.
ratio_col<-rep(1,length(names(Year_Geom_Means)))
ratio_col[gear_sensitive_species]<-2
ratio<-colMeans(Year_Geom_Means[15:20,])/colMeans(Year_Geom_Means[1:14,])                                                    
plot(ratio, pch=20,col=ratio_col,type='p',ylim=c(0,40))
abline(h=7, lty=2)
abline(h=1, lty=2)
text(ratio,col=ratio_col)


#List of column numbers of the comunity matrix, and names of species identified as gear-sensitive.
gear_sensitive_species<-c(1,2,3,7,8,9,10,11,12,14,15,16,20,22,25,28,30,31,32,33,34,39,40,42,47,54,55,56,57)
gear_sensitive_species_names = names(DFO_Dataset[,31:89])[gear_sensitive_species]

#mean density and subsets of mean densities across years
Year_Geom_Means_all<-Year_Geom_Means
Year_Geom_Means<-Year_Geom_Means[,-gear_sensitive_species]
Year_Geom_Means_4<-Year_Geom_Means[,c("GADUS_MORHUA","REINHARDTIUS_HIPPOGLOSSOIDES","HIPPOGLOSSOIDES_PLATESSOIDES","SEBASTES_MENTELLA")]
Year_Geom_Means_rare<-Year_Geom_Means[,-which(names(Year_Geom_Means)=="GADUS_MORHUA" |names(Year_Geom_Means)== "REINHARDTIUS_HIPPOGLOSSOIDES" | names(Year_Geom_Means)=="HIPPOGLOSSOIDES_PLATESSOIDES" | names(Year_Geom_Means)=="SEBASTES_MENTELLA")]


#Total Biomass####
DFO_Com_All<-DFO_Com
DFO_Com<-DFO_Com[,-gear_sensitive_species]
Total_Biomass<- ddply(DFO_Com,.variables=.(Year),
                      .fun=function(x){
                        time_series = rowSums(x)
                        jack<-jackknife(time_series,CalcZeroInfGeomDens)
                        return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
                      })
Total_Biomass_4<- ddply(DFO_Com[,c("GADUS_MORHUA","REINHARDTIUS_HIPPOGLOSSOIDES","HIPPOGLOSSOIDES_PLATESSOIDES","SEBASTES_MENTELLA")],.variables=.(Year),
                        .fun=function(x){
                          time_series = rowSums(x)
                          jack<-jackknife(time_series,CalcZeroInfGeomDens)
                          return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
                        })
Total_Biomass_rare<- ddply(DFO_Com[,-which(names(DFO_Com)=="GADUS_MORHUA" |names(DFO_Com)== "REINHARDTIUS_HIPPOGLOSSOIDES" | names(DFO_Com)=="HIPPOGLOSSOIDES_PLATESSOIDES" | names(DFO_Com)=="SEBASTES_MENTELLA")],.variables=.(Year),
                           .fun=function(x){
                             time_series = rowSums(x)
                             jack<-jackknife(time_series,CalcZeroInfGeomDens)
                             return(data.frame(Bmass=mean(jack$jack.values),jack.se=jack$jack.se))
                           })

						   
#NMDS Calculations####
#taxonomy
# Calculate Bray Curtis dissimilarity, perform NMDS, record site scores
# All spp
diss.tax<-vegdist(decostand(Year_Geom_Means,"total"),"bray")
mds.tax<-metaMDS(diss.tax) #3D : k=3
clust.tax<-scores(mds.tax,display="sites",scaling=1)

# 4 spp
diss.tax4<-vegdist(decostand(Year_Geom_Means_4,"total"),"bray")
mds.tax4<-metaMDS(diss.tax4) #3D : k=3
clust.tax4<-scores(mds.tax4,display="sites",scaling=1)

# Minus 4 spp
diss.taxmin4<-vegdist(decostand(Year_Geom_Means_rare,"total"),"bray")
mds.taxmin4<-metaMDS(diss.taxmin4) #3D : k=3
clust.taxmin4<-scores(mds.taxmin4,display="sites",scaling=1)



#Community Synchrony Calculations####
# code by Tarik
#modified by Patrick - April 15
community.sync.window <- function (d, twin=5, nrands=999) {
  years = as.numeric(rownames(d))
  starts=seq(from=1, to=length(years)-twin+1, by=1)
  ends=starts+twin-1
  mid=starts+(twin-1)/2
  
  results=matrix(nrow=length(starts), ncol=6, NA)
  for (i in 1:length(starts)) {
    locs=which(rownames(d) %in% years[starts[i]:ends[i]])
    c = community.sync(d[locs,], nrands=nrands)
    results[i,] = c(years[starts[i]], years[ends[i]],years[mid[i]], c$meancorr, c$obs, c$pval)
  }
  colnames(results)=c("start", "end","mid", "meancorr", "sync", "pval")
  results=as.data.frame(results)
  return(results)
}

#Synchrony calculations. 
twindow = 5
sync_allsp = community.sync.window(Year_Geom_Means[rownames(Year_Geom_Means)<1995,], twin=twindow)
sync_allsp2 = community.sync.window(Year_Geom_Means[rownames(Year_Geom_Means)>1994,], twin=twindow)
sync_common =community.sync.window(Year_Geom_Means_4[rownames(Year_Geom_Means)<1995,], twin=twindow)
sync_common2 =community.sync.window(Year_Geom_Means_4[rownames(Year_Geom_Means)>1994,], twin=twindow)
sync_rare=community.sync.window(Year_Geom_Means_rare[rownames(Year_Geom_Means)<1995,], twin=twindow)
sync_rare2=community.sync.window(Year_Geom_Means_rare[rownames(Year_Geom_Means)>1994,], twin=twindow)
alpha=0.05 


#Loading Voronoi polygon data ####
# This section is based on loading map layers from the data/map_files folder. 
# adjust paths to fit. 

#Loads previously calculated Voronoi polygon centroids
load("~/Dropbox/Groundfish working group shared files/data/map_files/meeting 2/voronoi_shapes.Rdat")
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
voronoi_data = ddply(voronoi_data,.(polygon), function(x) {
  x$Depth  = mean(x$Depth)
  return(x)
})

voronoi_data = ddply(voronoi_data, .(polygon, Year), 
                     function(x){
                       out_data = x[1,]
                       out_data$Temp_at_fishing = mean(x$Temp_at_fishing)
                       out_data[1,31:60] = apply(x[,31:60],MARGIN=2, 
                                                 FUN=CalcZeroInfGeomDens)
                       return(out_data)
                     })

                     
top4_sp = c("GADUS_MORHUA","SEBASTES_MENTELLA",
            "REINHARDTIUS_HIPPOGLOSSOIDES", "HIPPOGLOSSOIDES_PLATESSOIDES")
top4_sp = match(top4_sp, names(voronoi_data))
nontop_sp = (31:60)[!(31:60)%in%top4_sp]

group_counts = plyr::count(subset(voronoi_data,Year>1980), .(polygon))
group_counts = group_counts$polygon[group_counts$freq>=min_years_in_group]
voronoi_data = voronoi_data[voronoi_data$polygon%in%group_counts,]


#Diversity by Distance Calculations####
#code by Eric
#Functions
CalcGeoDistMat = function(lat, long){
  # This function calculates the geodesic distance matrix between all points in a set of 
  # lat/long coordinates.
  stopifnot(is.numeric(lat)&is.numeric(long))
  n_data= length(lat)
  lat_data = expand.grid(lat,lat)
  long_data = expand.grid(long,long)
  dist_calc = geoDist(lat1=lat_data[,1],lat2=lat_data[,2],
                      lon1=long_data[,1],lon2=long_data[,2],NAOK=T)
  dist_mat = matrix(dist_calc,nrow=n_data)
  return(dist_mat)
}

CalcSpaceDepthDists = function(x,use_geodesic){
  #Calculates in a single function both depth and distance dissimilarity between polygons
  if(use_geodesic){
    spatial_dist = CalcGeoDistMat(x$LAT_DEC,x$LONG_DEC)
    spatial_dist = spatial_dist[lower.tri(spatial_dist,diag=F)]
  }else{
    spatial_dist = as.vector(dist(cbind(x$LAT_DEC,x$LONG_DEC)))
  }
  depth_dist = as.vector(dist(x$Depth))
  return(data.frame(space_dist=spatial_dist, depth_dist=depth_dist))
}

EstimateDistR2 = function(x){
  # Estimates the partial R^2 explained of species dissimilarity
  # by differences in distance and depth (after transformations)
  distance_scl = x$distance_scl
  boot_data = x
  boot_data = subset(boot_data, !is.na(distance_scl))
  distance_scl = distance_scl[!is.na(distance_scl)]
  var_decomp = varpart(Y=distance_scl, ~space_dist_scl,~depth_dist_scl,
                       data=boot_data)
  return(var_decomp$part$indfract$Adj.R.square[c(1,3)])
}

CalculateRSquareValues = function(x){
  #Wrapper function for EstimateDistR2 that also calculates
  #jackknife se values for the estimates
  var_frac = EstimateDistR2(x)
  if(boot_errors){
    unique_sites = unique(c(x$site_1,x$site_2))
    n_sites=length(unique_sites)
    jack_R_sqr = matrix(rep(0,times=n_sites*2),ncol=2)
    jack_mean = rep(0,times=n_sites)
    jack_var =rep(0,times=n_sites)
    for(i in 1:n_sites){
      indexes = x$site_1!=i&x$site_2!=i
      current_data = x[indexes,]
      jack_mean[i] = mean(current_data$distance_scl,na.rm=T)
      jack_var[i] =var(current_data$distance_scl,na.rm=T)
      jack_R_sqr[i,]  =  EstimateDistR2(current_data)
    }
    jack_R_sqr_est = colMeans(jack_R_sqr)
    jack_mean_est = mean(jack_mean)
    jack_var_est = mean(jack_var)
    jack_R_sqr_se = sqrt((n_sites-1)*(colMeans(jack_R_sqr^2)-jack_R_sqr_est^2))
    jack_mean_se = sqrt((n_sites-1)*var(jack_mean))
    jack_var_se = sqrt((n_sites-1)*var(jack_var))
  }
  model_slopes = coef(lm(distance_scl~space_dist_scl+depth_dist_scl,
                         data=x))
  output_data = data.frame(coef = c("Mean community distance","variance",
                                    "Spatial Distance","Depth"),
                           R2_value = c(mean(x$distance_scl,na.rm=T),
                                        var(x$distance_scl,na.rm=T), 
                                        var_frac),
                           slope_value = c(mean(x$distance_scl,na.rm=T),
                                           var(x$distance_scl,na.rm=T),
                                           model_slopes[2:3]))
  output_data$coef= factor(output_data$coef, 
                           levels= c("Mean community distance",
                                     "variance","Spatial Distance",
                                     "Depth"))
  if(boot_errors){
    
    R2_value = c(jack_mean_est,jack_var_est, jack_R_sqr_est)
    output_data$se = c(jack_mean_se,jack_mean_se,jack_R_sqr_se)
  }
  return(output_data)
}

R2NameLabeller =function(variable, value){
  #Function to aid in plotting. Parses various expressions into caption labels.
  dist_labels= list("Mean community distance"=expression("Mean scaled community distance"), 
                    "variance" = expression("Variance of community distance"),
                    "Spatial Distance"= expression(paste('Partial ',R^2, ' explained by distance')),
                    "Depth"= expression(paste('Partial ',R^2, ' explained by depth differential'))
  )
  R2_labels = list(
    "R2_value" = "Total community",
    "R2_top4"="four most abundant",
    "R2_nontop"="Remaining species"
  )
  if(variable=="coef"){
    return(dist_labels[value])
  }else{
    return(R2_labels[value])
  }
}


#Options for diversity-distance plots
com_dist = "bray" #uses Bray-Curtis dissimilarity
geodesic_dist = T #Use geodesic distances for calculating inter-polygon distance
use_voronoi_data =T
boot_errors = T #Should jackknife se values be calculated?
include_variance= F
include_mean = F
min_years_in_group = 30 #specifies how many years a polygon has to have observations in it for it to be used in the analysis

space_trans = function(x) {
  #Returns centered, scaled, log-10 transformation of spatial distances
  return(scale(log10(x)))
}
depth_trans = function(x) {
  #Returns centered, scaled depth dissimilarities
  return(scale(x,center=T,scale=T))
}

comm_trans = function(x, logt=T,scale_to_max=T,com_dist="euclidean", pad_val =0.01,
                      center_val=F,scale_val = F) {
  # a modified logistic transformation for species dissimilarity. Transforms a value
  # ranging from 0 to some maximum value into a range from -inf to +inf. 
  if(com_dist=="euclidean"){
    max_val = sqrt(2)
  }else if(com_dist=="bray"){
    max_val = 1 
  }else stop("community distance ",com_dist, " is not a recognized distance currently") 
  new_x = x
  if(scale_to_max) new_x = (new_x+pad_val)/(max_val+pad_val-new_x)
  if(logt) new_x =log(new_x)
  new_x= scale(new_x,center=center_val,scale=scale_val)
  return(new_x)
}



#Calculating spatial, depth, and community distances
mean_dist_data = ddply(voronoi_data, .(Year), 
                       function(x){
                         dist_data = CalcSpaceDepthDists(x[,c("LAT_DEC","LONG_DEC","Depth")],
                                                         use_geodesic=geodesic_dist)
                         dist_data$sp_dist = as.vector(vegdist(decostand(x[,31:60],method="total"),
                                                               method=com_dist))
                         dist_data$sp_top4_dist= as.vector(vegdist(decostand(x[,top4_sp],method="total"),
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
                           sp_top4_dist_scl = comm_trans(sp_top4_dist,com_dist=com_dist),
                           sp_nontop_dist_scl = comm_trans(sp_nontop_dist,com_dist=com_dist)
)

pairs = str_split_fixed(mean_dist_data$pairs,pattern=" ",n=2)
mean_dist_data$site_1  = as.numeric(pairs[,1])
mean_dist_data$site_2  = as.numeric(pairs[,2])
mean_dist_data = melt(mean_dist_data,measure.vars=c("sp_dist_scl","sp_top4_dist_scl","sp_nontop_dist_scl"),
                      variable.name="community_subset",value.name = "distance_scl")

#Calculating regression models 
div_dist_models = ddply(mean_dist_data,.(Year,community_subset), 
                        .fun=CalculateRSquareValues)

if(!include_mean){
  div_dist_models = subset(div_dist_models, coef!="Mean community distance")
}
if(!include_variance){
  div_dist_models = subset(div_dist_models, coef!="variance")
}

ComV<-c("sp_dist_scl","sp_top4_dist_scl","sp_nontop_dist_scl")



#Spatial Cluster Analysis####
PlotMultipleGgplotObjs <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

clust_map_data = subset(DFO_Dataset[,!names(DFO_Dataset)%in%gear_sensitive_species_names], Year>1980)
clust_map_data[,c("polygon","LONG_DEC","LAT_DEC")]= over(SpatialPoints(cbind(clust_map_data$LONG_DEC, clust_map_data$LAT_DEC)), voronoi_shapes)
clust_map_data = clust_map_data %>% group_by(polygon) %>% mutate(Depth=mean(Depth,na.rm = T))%>% ungroup()

clust_map_data_all_years<-clust_map_data

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

# Calculates geom-mean density for each species in each year group in each polygon
clust_map_data = ddply(clust_map_data, .(polygon, Year_group),function(x){
  out_data = x[1,]
  out_data[1,31:60] = apply(x[,31:60],MARGIN=2, 
                            FUN=CalcZeroInfGeomDens)
  #out_data[1,31:60] = out_data[1,31:60]/sum(out_data[1,31:60])
  return(out_data)
})
#removing one row with all zero data
clust_map_data = clust_map_data[!rowSums(clust_map_data[,31:60])==0,]


# Clusters the data into 7 groupings, as well as clustering the top 4 and remaining species
top4_sp = c("GADUS_MORHUA","SEBASTES_MENTELLA","REINHARDTIUS_HIPPOGLOSSOIDES", "HIPPOGLOSSOIDES_PLATESSOIDES")
top4_sp = match(top4_sp, names(clust_map_data))
nontop_sp = (31:60)[!(31:60)%in%top4_sp]
n_cluster= 7
clust_tree = hclust(vegdist(decostand(clust_map_data[,31:60],method = "total")))
top4_voronoi_dissim = vegdist(decostand(clust_map_data[,top4_sp],method = "total"))
top4_voronoi_dissim[which(is.nan(top4_voronoi_dissim))]=0
clust_tree_top4 = hclust(top4_voronoi_dissim)
remain_dissim=vegdist(decostand(clust_map_data[,nontop_sp],method = "total"))
remain_dissim[which(is.nan(remain_dissim))]=0
clust_tree_remain = hclust(remain_dissim)


clust_map_data$cluster = cutree(clust_tree,k = n_cluster)
clust_map_data$cluster = factor(ifelse(clust_map_data$cluster==8,7, clust_map_data$cluster ) )
clust_map_data$cluster_top4 = (cutree(clust_tree_top4,k = 4))


clust_map_data$cluster_remain = (cutree(clust_tree_remain,k = 7))
#clust_map_data$cluster_remain = mapvalues(clust_map_data$cluster_remain,
#                                          1:4, 
#                                          c(1,2,3,4))
clust_map_data$cluster_remain =factor(clust_map_data$cluster_remain)
clust_map_data$cluster_top4 =factor(clust_map_data$cluster_top4)

clust_map_data$cod_percent<-clust_map_data$GADUS_MORHUA/rowSums(clust_map_data[,31:60])


dendro_plot_data =dendro_data(clust_tree) 
dendro_top4_plot_data = dendro_data(clust_tree_top4)
dendro_remain_plot_data = dendro_data(clust_tree_remain)


#### Turning shape files into plotable polygon data ####
polygon_base_data = fortify(voronoi_shapes)
names(polygon_base_data)[1:2] = c("LONG_DEC","LAT_DEC")
polygon_base_data$polygon = as.numeric(polygon_base_data$id)+1

polygon_base_data = polygon_base_data[!polygon_base_data$polygon==171,] #polygon 171 and 170 are duplicates here
n_poly_pts = nrow(polygon_base_data)
polygon_map_data = polygon_base_data[rep(1:n_poly_pts, 
                                         times=length(year_groups)),]
polygon_map_data$Year_group = rep(year_groups,each=n_poly_pts)  

polygon_map_data = left_join(polygon_map_data, 
                             clust_map_data[,c("polygon","Year_group","cluster",
                                               "cluster_top4","cluster_remain", "Depth","cod_percent")])



#### Building a depth-voronoi polygon dataset ####
depth_voronoi_map_data = left_join(polygon_base_data,clust_map_data[,c("polygon","cluster", "Depth")])
depth_voronoi_map_data = depth_voronoi_map_data%>%
  group_by(polygon,order) %>% 
  summarise(Depth=mean(Depth,na.rm = T), LAT_DEC=mean(LAT_DEC),LONG_DEC=mean(LONG_DEC))
depth_voronoi_map_data$Depth = with(depth_voronoi_map_data,
                                    trunc(Depth/150)+1)
depth_voronoi_map_data$Depth = factor(with(depth_voronoi_map_data,
                                           ifelse(Depth==10,9,Depth)))

#### Building aggregate composition data for each cluster ####
## Averaging over each cluster, and setting each row to sum to one,then taking  column means
cluster_composition = ddply(clust_map_data, .(cluster),
                            function(x){
                              x[,31:60] = decostand(x[,31:60],method ="total") 
                              out_data = data.frame(species = c(names(x)[top4_sp]),
                                                    proportion = rep(0,times=4))
                              for(i in 1:4){
                                out_data[i,2] = mean(x[,top4_sp[i]])
                              }
                              return(out_data)
                            })
cluster_composition$species = factor(cluster_composition$species, 
                                     levels=c("GADUS_MORHUA","REINHARDTIUS_HIPPOGLOSSOIDES",
                                              "HIPPOGLOSSOIDES_PLATESSOIDES",
                                              "SEBASTES_MENTELLA"),
                                     labels=c("G. Morhua", "R. Hippoglossoides","H. Platessoides",
                                              "S. Mentella"))
cluster_order = ddply(cluster_composition,.(cluster), function(x){
  proportion = sum(x$proportion)
  return(data.frame(proportion=proportion))
})
cluster_order = cluster_order$cluster[order(cluster_order$proportion,
                                            decreasing = T)]
cluster_composition$cluster =factor(cluster_composition$cluster,
                                    levels = cluster_order)

polygon_outline = unionSpatialPolygons(voronoi_shapes,
                                       IDs = rep(1, nrow(coordinates(voronoi_shapes))))
polygon_outline = fortify(polygon_outline)  
names(polygon_outline)[1:2] = c("LONG_DEC","LAT_DEC")


cluster_palette=c("#1F78B4","#E31A1C","#A6CEE3",'#984ea3',"#33A02C","#b2df8a",'#fdbf6f')


cluster_palette_top4 = c("#1F78B4", "#E31A1C", "#33A02C", "#555555" )
lat_long_labels = list(scale_x_continuous(breaks = c(-55,-50), 
                                          labels = c(expression(55*degree*0*minute*0*second*'W'),
                                                     expression(50*degree*0*minute*0*second*'W'))),
                       scale_y_continuous(breaks = c(50,55), 
                                          labels = c(expression(50*degree*0*minute*0*second*'N'),
                                                     expression(55*degree*0*minute*0*second*'N'))))

polygon_plot = ggplot(aes(x=LONG_DEC, y=LAT_DEC),
                      data= polygon_map_data)+
  geom_polygon(data=polygon_outline, fill=NA, colour="black",
               aes(group=piece,order=order))+
  coord_map(projection="albers",lat0=min(polygon_map_data$LAT_DEC),
            lat1= max(polygon_map_data$LAT_DEC))+
  facet_grid(.~Year_group)+
  lat_long_labels+
  scale_size_continuous("Depth (m)",range = c(0.25,0.01))+
  scale_colour_gradient("Depth (m)",high= "#093f9c", low = "#d6d6ff",guide="none")+
  theme_bw()+
  theme(text=element_text(size=15),panel.grid.minor=element_blank(),legend.position="bottom",
        axis.text = element_text(size=10),
        axis.title=element_blank())

polygon_total_plot =polygon_plot + 
  scale_fill_manual(values= cluster_palette,guide="none",na.value = "#aaaaaa")+
  geom_polygon(aes(group= polygon,fill= cluster,order=order))

polygon_top4_plot = polygon_plot + 
  scale_fill_manual(values= cluster_palette_top4,guide="none",na.value = "#aaaaaa")+
  geom_polygon(aes(group= polygon,fill= cluster_top4,order=order))

polygon_remain_plot =polygon_plot + 
  scale_fill_manual(values= cluster_palette,guide="none")+
  geom_polygon(aes(group= polygon,fill= cluster_remain,order=order))


polygon_cod_plot=polygon_plot+
  scale_fill_continuous(low="#ffffff",high="#e41a1c",guide="none",na.value="#aaaaaa")+
  geom_polygon(aes(group= polygon,fill= cod_percent,order=order))



#Depth map with Voronoi polygons
voronoi_depth_classes = paste(seq(0,1200, by=150),seq(149,1349, by=150),sep="-")
voronoi_depth_classes[9] = "1200+"
depth_voronoi_plot = ggplot(aes(x=LONG_DEC, y=LAT_DEC),data=depth_voronoi_map_data)+
  geom_polygon(aes(group= polygon,fill= Depth,order=order))+
  geom_polygon(data=polygon_outline, fill=NA, colour="black",
               aes(group=piece,order=order))+
  
  coord_map(projection="albers",lat0=min(polygon_map_data$LAT_DEC),
            lat1= max(polygon_map_data$LAT_DEC))+
  theme_bw()+
  lat_long_labels+
  scale_fill_brewer("Depth (m)",palette="Blues",breaks=factor(1:9), 
                    labels=voronoi_depth_classes)+
  theme(text=element_text(size=15),panel.grid.minor=element_blank(),
        legend.position="right", 
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.text = element_text(size=10),
        axis.title=element_blank())

comp_plot = ggplot(aes(x=cluster, y=proportion),
                   data=cluster_composition)+
  geom_bar(stat="identity", aes(fill=species,order= species))+
  scale_fill_brewer("Species",palette="Set1")+
  annotate(x=factor(1:7),y = rep(-0.1,times=7),
           colour=cluster_palette,
           geom="point",size=10,shape=18)+
  coord_cartesian(ylim=c(-0.2,1))+
  scale_y_continuous("Average % of community",breaks = c(0,0.25,0.5,0.75,1))+
  scale_x_discrete("Community cluster")+
  theme_bw()+
  theme(text= element_text(size=15),panel.grid.major.x = element_blank(),axis.text.x=element_blank(),
        axis.ticks.x  = element_blank(),legend.direction="horizontal",
        legend.position = c(0.5,0.9))

dendro_plot_base = list(
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)),
  theme_bw(),
  theme(legend.position="none",panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text=element_blank(),axis.ticks=element_blank(),
        axis.title= element_blank()))


dendro_total_plot = ggplot(segment(dendro_plot_data)) +
  geom_tile(data=label(dendro_plot_data),size=0.5,
            aes(label=label, x=x, y=-0.025, 
                fill= clust_map_data$cluster[as.numeric(as.character(dendro_plot_data$labels$label))]),
            height=0.05,width=1)+
  scale_fill_manual(values=cluster_palette)+
  dendro_plot_base

dendro_top4_plot = ggplot(segment(dendro_top4_plot_data)) +
  geom_tile(data=label(dendro_top4_plot_data),size=0.5,
            aes(label=label, x=x, y=-0.025, 
                fill= clust_map_data$cluster_top4[as.numeric(as.character(dendro_top4_plot_data$labels$label))]),
            height=0.05,width=1)+
  scale_fill_manual(values=cluster_palette_top4)+
  dendro_plot_base+
  annotate("text",x=-100,y=0.5,label="Top 4 species",angle=90,size=7)

dendro_remain_plot = ggplot(segment(dendro_remain_plot_data)) +
  geom_tile(data=label(dendro_remain_plot_data),size=0.5,
            aes(label=label, x=x, y=-0.025, 
                fill= clust_map_data$cluster_remain[as.numeric(as.character(dendro_remain_plot_data$labels$label))]),
            height=0.05,width=1)+
  scale_fill_manual(values=cluster_palette)+
  dendro_plot_base+
  annotate("text",x=-100,y=0.5,label="Remaining species",angle=90,size=7)


#Figures####
axis.V<-1.1
label.V<-1.2

#Figure 1#####
pdf("Fig. 1 - 5yr.pdf", height=7,width=7)
par(mfrow=c(3,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)
plot(Total_Biomass$Bmass[1:14]~c(1981:1994),type='l',lty=1, lwd=2, ylab="Biomass (kg/tow)", xlab=NA,ylim=c(0,210),cex.lab=label.V,cex.axis=axis.V,xaxt='n',xlim=c(1981,2013))
lines(Total_Biomass$Bmass[15:33]~c(1995:2013),lty=1, lwd=2)
plotCI(c(1981:2013),Total_Biomass$Bmass,uiw=Total_Biomass$jack.se,add=T,pch=NA,sfrac=0)
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
plotCI(c(1981:2013),Total_Biomass_4$Bmass,uiw=Total_Biomass_4$jack.se,add=T,pch=NA,sfrac=0,lty=2)
lines(c(1981:1994),Total_Biomass_4$Bmass[1:14],type='l',lty=2, lwd=2, ylim=c(0,1))
lines(c(1995:2013),Total_Biomass_4$Bmass[15:33],type='l',lty=2, lwd=2, ylim=c(0,1))
plotCI(c(1981:2013),Total_Biomass_rare$Bmass,uiw=Total_Biomass_rare$jack.se,add=T,pch=NA,sfrac=0,lty=3)
lines(c(1981:1994),Total_Biomass_rare$Bmass[1:14],type='l',lty=3, lwd=2, ylim=c(0,1))
lines(c(1995:2013),Total_Biomass_rare$Bmass[15:33],type='l',lty=3, lwd=2, ylim=c(0,1))
abline(v=Eras, lwd=1,lty=2, col=8)
legend("topright",c("total community","4 commercial species","remaining species"),lty=c(1,2,3),lwd=2,bty='o',bg="white",box.col="#FF003300",inset=0.005,cex=1.1)
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

phase<-cbind(Total_Biomass_4$Bmass[5:33]-Total_Biomass_4$Bmass[1:29],1985:2013)
segments(x0 = sync_common[-nrow(sync_common),"end"],y0 = sync_common[-nrow(sync_common),"sync"],x1 = sync_common[-1,"end"],y1 = sync_common[-1,"sync"],col= grad.V[(((phase[,1]/max(abs(phase[,1])))+1)*5)+1],lwd=3, lty=2)
segments(x0 = sync_common2[-nrow(sync_common2),"end"],y0 = sync_common2[-nrow(sync_common2),"sync"],x1 = sync_common2[-1,"end"],y1 = sync_common2[-1,"sync"],col= grad.V[(((phase[11:29]/max(abs(phase)))+1)*5)+1],lwd=3,lty=2)
points(sync_common[,"end"], sync_common[,"sync"], col="black", pch=21, 
       bg=ifelse(sync_common[,"pval"] < alpha, 'black', 'white'))
points(sync_common2[,"end"], sync_common2[,"sync"], col="black", pch=21, 
       bg=ifelse(sync_allsp2[,"pval"] < alpha, 'black', 'white'))
#lines(sync_rare[,"end"], sync_rare[,"sync"], col="black", lty = 1, lwd = 3)

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

#Figure 3####
ColV2<-c("#253494","#fc4e2a","#4eb3d3")
ColV2<-c(brewer.pal(3,"Dark2")[1:2],"#e425a7")

colV<-c(rep(ColV2[1],9),rep(ColV2[2],5),rep(ColV2[3],23))
colVseg<-c(rep(ColV2[1],9),rep(ColV2[2],4),NA,rep(ColV2[3],23))

Percent_Com<-Year_Geom_Means/rowSums(Year_Geom_Means)

pdf("Fig. 2.pdf", height=5,width=7)
par(mar=c(4,4,1,1),oma=c(2,2,2,2),las=1)
layout(matrix(c(1,2,3,4,4,4),2,3, byrow=T))
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

ordiplot(mds.tax4,type="n",xlab="NMDS1",ylab=NA,xlim=c(-0.4,0.3),main="4 commercial species",cex.main=1,axes=F,frame=F,cex.lab=label.V, cex.axis=axis.V)
points(clust.tax4[,1],clust.tax4[,2],pch=19,col=colV,cex=1)
segments(clust.tax4[s,1],clust.tax4[s,2],clust.tax4[s+1,1],clust.tax4[s+1,2],col=colVseg)
segments(clust.tax4[14,1],clust.tax4[14,2],clust.tax4[15,1],clust.tax4[15,2], col=colV[14],lty=2)
text(clust.tax4[c(1,10,15,33),1]+c(0.03,0,0,0.02), clust.tax4[c(1,10,15,33),2]+c(0.04,-0.04,0.04,-0.04),rownames(Year_Geom_Means)[c(1,10,15,33)],cex=1,col=colV[c(1,10,15,33)])
axis(1,at=c(-0.3,0,0.3))
axis(2,at=c(-0.3,0,0.3))
legend("bottomright",paste("S = ",round(mds.tax4 $stress,digits=2)), bty='n', cex=1.2)
legend("topleft", "b", bty='n', cex=1.8,inset=c(-0.15,0.04))

ordiplot(mds.taxmin4,type="n",xlab="",ylab=NA,xlim=c(-0.4,0.4),ylim=c(-0.3,0.45),main="Remaining species",cex.main=1,axes=F,frame=F,cex.lab=label.V, cex.axis=axis.V)
points(clust.taxmin4[,1],clust.taxmin4[,2],pch=19,col=colV,cex=1)
segments(clust.taxmin4[s,1],clust.taxmin4[s,2],clust.taxmin4[s+1,1],clust.taxmin4[s+1,2],col=colVseg)
segments(clust.taxmin4[14,1],clust.taxmin4[14,2],clust.taxmin4[15,1],clust.taxmin4[15,2], col=colV[14],lty=2)
axis(1,at=c(-0.3,0,0.3))
axis(2,at=c(-0.3,0,0.3))
text(clust.taxmin4[c(1,10,15,33),1]+c(0,0.07,0,-0.03), clust.taxmin4[c(1,10,15,33),2]+c(-0.08,-0.05,0.05,-0.05),rownames(Year_Geom_Means)[c(1,10,15,33)],cex=1,col=colV[c(1,10,15,33)])
legend("bottomright",paste("S = ",round(mds.taxmin4 $stress,digits=2)), bty='n', cex=1.2)
legend("topleft", "c", bty='n', cex=1.8,inset=c(-0.15,0.04))

par(pty='m')
plot(Percent_Com$GADUS_MORHUA*100~c(1981:2013),type='n',lwd=2, ylab="% of total biomass", xlab="Year", ylim=c(0,105),col=ColV[1],cex.lab=label.V,cex.axis=axis.V)
polygon(c(1981,1981:1994,1994),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$SEBASTES_MENTELLA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[1:14]*100,0), col=ColV[4])
polygon(c(1981,1981:1994,1994),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[1:14]*100,0), col=ColV[3])
polygon(c(1981,1981:1994,1994),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA)[1:14]*100,0), col=ColV[2])
polygon(c(1981,1981:1994,1994),c(0,Percent_Com$GADUS_MORHUA[1:14]*100,0), col=ColV[1])
polygon(c(1995,1995:2013,2013),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$SEBASTES_MENTELLA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[15:33]*100,0), col=ColV[4])
polygon(c(1995,1995:2013,2013),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA+Percent_Com$HIPPOGLOSSOIDES_PLATESSOIDES)[15:33]*100,0), col=ColV[3])
polygon(c(1995,1995:2013,2013),c(0,(Percent_Com$REINHARDTIUS_HIPPOGLOSSOIDES+Percent_Com$GADUS_MORHUA)[15:33]*100,0), col=ColV[2])
polygon(c(1995,1995:2013,2013),c(0,Percent_Com$GADUS_MORHUA[15:33]*100,0), col=ColV[1])
abline(v=Eras, lwd=1,lty=2, col=8)
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
legend("topleft", "d", bty='n', cex=1.8,inset=c(-0.035,0))
legend("topright",legend = c("cod","halibut","plaice","redfish"),lwd=2,col=c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),box.col="#FF003300",inset=0.005,cex=0.85)
dev.off()


#Figure 3####
pdf("Fig. 3.pdf",width = 12,height = 10)
{
PlotMultipleGgplotObjs(polygon_cod_plot,polygon_total_plot,depth_voronoi_plot, comp_plot, 
                       layout=matrix(c(1,1,1,2,2,2,3,4,4),nrow = 3,byrow = T))
}
dev.off()

#Figure 4####
polygons = unique(voronoi_data$polygon)
diss.tax<-as.matrix(vegdist(decostand(Year_Geom_Means,"total"),"bray"))
spatial_BC <- matrix(NA, nrow=length(1981:2013), 
                     ncol=length(polygons))
for(i in 1:length(polygons)){
  current_voronoi <- subset(voronoi_data,polygon==polygons[i]&Year>1980)
  current_dist <- vegdist(current_voronoi[,31:60])
  current_dist <- as.matrix(current_dist)
  spatial_BC[current_voronoi$Year-1980,i] <- current_dist[,1]
  rm(current_dist, current_voronoi)
}

Ref_bmass_1981<-Total_Biomass$Bmass/Total_Biomass$Bmass[1]
#Ref_bmass_1981<-Ref_bmass_1981[-1]
Ref_bmass_1981<-(Ref_bmass_1981-min(Ref_bmass_1981))/(Ref_bmass_1981-min(Ref_bmass_1981))[1]*100

Ref_cod_1981<-Year_Geom_Means$GADUS_MORHUA/Year_Geom_Means$GADUS_MORHUA[1]
#Ref_cod_1981 <-Ref_cod_1981[-1]
Ref_cod_1981 <-(Ref_cod_1981-min(Ref_cod_1981))/(Ref_cod_1981-min(Ref_cod_1981))[1]*100


Ref_com_1981<-(1-decostand(diss.tax[,1],"range"))*100

Ref_spatial_1981<-(1-spatial_BC[-1,])*100
#Ref_spatial_1981<-apply(Ref_spatial_1981,2,function(x){x/x[1]})
matplot(Ref_spatial_1981,type='l', lty=1,col=1)

pdf("Fig. 4.pdf",height=5,width=7)
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



#Fig S1####
pdf("~/Dropbox/Groundfish working group shared files/Manuscript 1/Figures/Fig. S1.pdf",width = 15,height = 10)
{
PlotMultipleGgplotObjs(polygon_total_plot ,
                       polygon_top4_plot, 
                       polygon_remain_plot,
                       layout=matrix(c(1,2,3),nrow = 3,byrow = T))
}
dev.off()


#Figure S2####
pdf("Fig. S2.pdf", height=7,width=7)
par(mfrow=c(3,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)
plot(c(1981:2013),div_dist_models[div_dist_models$community_subset=="sp_dist_scl" & div_dist_models$coef=="Spatial Distance","R2_value"],type='n', ylab=expression(paste("Distance R"^"2")), ylim=c(0,0.3), xlab=NA,xaxt='n', cex.lab=label.V,cex.axis=axis.V, xlim=c(1980,2013))
for(i in 1:3){
  plotCI(c(1981:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","R2_value"],uiw=div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","se"],pch=NA,slty=i,sfrac=0,add=T)
  lines(c(1981:1994),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","R2_value"][1:14], lty=i, lwd=2)
  lines(c(1995:2013),div_dist_models[div_dist_models$community_subset==ComV[i] & div_dist_models$coef=="Spatial Distance","R2_value"][15:33], lty=i, lwd=2)
}
abline(v=Eras, lwd=1,lty=2, col=8)
axis(1,seq(1980,2010,by=5),cex=1.4,labels=F,tick=T)
legend("topleft", "A",bty='n', cex=1.8, inset=c(-0.04,-0.025))
legend(x="topright", legend=c("Total community", "4 large commercial species", "Remaining species"),
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



