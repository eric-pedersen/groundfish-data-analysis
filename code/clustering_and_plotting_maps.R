#Clear workspace
rm(list=ls(all=TRUE))

#Packages####
library(plyr)
library(reshape2)
library(dplyr)
library(vegan)
library(bootstrap)
library(RColorBrewer)
library(ggplot2)
library(ggdendro)
library(sp)
library(maptools)
library(rgeos)
library(ggpolypath)
library(stringr)



#Loading data and functions ####
source("code/functions.R")
clust_map_data=read.csv("data/clust_map_data.csv",stringsAsFactors = F)
load("data/voronoi_shapes.Rdat")
year_groups = c("1981-1984","1985-1989","1990-1994","1995-2001","2002-2006", "2007-2013")



# Clusters the data into 7 groupings, as well as clustering the top 4 and remaining species
top4_sp = c("GADUS_MORHUA","SEBASTES_MENTELLA","REINHARDTIUS_HIPPOGLOSSOIDES", 
            "HIPPOGLOSSOIDES_PLATESSOIDES")
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

#Merge polygons into contiguous regions with the same clusters ####
#This is based off code from http://gis.stackexchange.com/questions/63577/joining-polygons-in-r
year_group_polygons = list()
n_vals = length(unique(paste(polygon_map_data$Year_group, polygon_map_data$polygon)))
year_group_data = data.frame(cluster=rep(0, times=n_vals),
                             Year_group = rep(0, times= n_vals),
                             polygon = rep(0, times=n_vals))
counter= 1
for(i in unique(polygon_map_data$polygon)){
  for(j in unique(polygon_map_data$Year_group)){
    current_data =polygon_map_data%>%
      filter(Year_group==j, polygon==i)%>%
      select(LONG_DEC,LAT_DEC,cluster,Year_group,polygon)
    current_id = paste(i,j, sep = "_")
    year_group_polygons[[counter]] = Polygon(coords=current_data[,1:2],hole = F)
    year_group_polygons[[counter]] = Polygons(list(year_group_polygons[[counter]]),current_id)
    year_group_data[counter,] = current_data[1,3:5]
    counter = counter+1
  }
}
rownames(year_group_data) = paste(year_group_data$polygon, 
                                  year_group_data$Year_group,sep="_")


polygon_map_total_data = SpatialPolygonsDataFrame(SpatialPolygons(year_group_polygons), 
                             year_group_data)
#clears up some issues with the shapes of the polygons having intersections
#based on code from http://gis.stackexchange.com/questions/163445/r-solution-for-topologyexception-input-geom-1-is-invalid-self-intersection-er
polygon_map_total_data = gBuffer(polygon_map_total_data,byid = T,width = 0) 
polygon_map_total_data = unionSpatialPolygons(polygon_map_total_data, 
                                        paste(year_group_data$cluster, 
                                              year_group_data$Year_group, sep= " "))

polygon_map_total_data= fortify(polygon_map_total_data) %>%
  mutate(cluster = str_split_fixed(id," ",n = 2)[,1],
         Year_group =str_split_fixed(id," ",n = 2)[,2],
         LONG_DEC= long, LAT_DEC= lat,
         cluster= ifelse(cluster=="NA", NA, cluster))



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

cluster_palette = c("#FFCF4C", "#13ABDA", "#384A81", "#34A576", "#D17634","#6E2152","#A1D664") 

species_palette = c("#13ABDA", "#FFD578", "#D17634","#34A576")
cluster_palette_top4 = c("#1F78B4", "#FF4040", "#33A02C", "#555555" )
lat_long_labels = list(scale_x_continuous(breaks = c(-55,-50), 
                                          labels = c(expression(55*degree*0*minute*0*second*'W'),
                                                     expression(50*degree*0*minute*0*second*'W'))),
                       scale_y_continuous(breaks = c(50,55), 
                                          labels = c(expression(50*degree*0*minute*0*second*'N'),
                                                     expression(55*degree*0*minute*0*second*'N'))))


polygon_plot = ggplot(aes(x=LONG_DEC, y=LAT_DEC),
                      data= polygon_map_data)+
  geom_polygon(fill=NA, colour="black",aes(group=group))+
  coord_map(projection="albers",lat0=min(polygon_map_data$LAT_DEC),
            lat1= max(polygon_map_data$LAT_DEC))+
  facet_grid(.~Year_group)+
  lat_long_labels+
  scale_size_continuous("Depth (m)",range = c(0.25,0.01))+
  scale_colour_gradient("Depth (m)",high= "#093f9c", low = "#d6d6ff",
                        guide="none")+
  theme_bw()+
  theme(text=element_text(size=15),panel.grid.minor=element_blank(),
        legend.position="bottom",
        axis.text = element_text(size=10),
        axis.title=element_blank())


polygon_total_plot =polygon_plot + 
  scale_fill_manual(values= cluster_palette,guide="none",na.value = "#aaaaaa")+
  geom_polypath(data=polygon_map_total_data,
               aes(group= group,fill= cluster),
               col="black",size=0.5)

polygon_top4_plot = polygon_plot + 
  scale_fill_manual(values= cluster_palette_top4,guide="none",
                    na.value = "#aaaaaa")+
  geom_polygon(aes(group= polygon,fill= cluster_top4,order=order))

polygon_remain_plot =polygon_plot + 
  scale_fill_manual(values= cluster_palette,guide="none")+
  geom_polygon(aes(group= polygon,fill= cluster_remain,order=order))


polygon_cod_plot=polygon_plot+
  scale_fill_continuous(low="#ffffff",high="#13ABDA",guide="none",
                        na.value="#aaaaaa")+
  geom_polygon(aes(group= polygon,fill= cod_percent,order=order))



#Depth map with Voronoi polygons
voronoi_depth_classes = paste(seq(0,1200, by=150),seq(149,1349, by=150),sep="-")
voronoi_depth_classes[9] = "1200+"
depth_voronoi_plot = ggplot(aes(x=LONG_DEC, y=LAT_DEC),
                            data=depth_voronoi_map_data)+
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
  scale_fill_manual(values= species_palette)+
  annotate(x=factor(1:7),y = rep(-0.1,times=7),
           colour=cluster_palette,
           geom="point",size=10,shape=18)+
  coord_cartesian(ylim=c(-0.2,1))+
  scale_y_continuous("Average % of community",breaks = c(0,0.25,0.5,0.75,1))+
  scale_x_discrete("Community cluster")+
  theme_bw()+
  theme(text= element_text(size=15),panel.grid.major.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x  = element_blank(),legend.direction="horizontal",
        legend.position = c(0.5,0.9))



#Figure 4####
pdf("figures/Fig. 4.pdf",width = 12,height = 7)
{
  PlotMultipleGgplotObjs(polygon_total_plot,depth_voronoi_plot, comp_plot, 
                         layout=matrix(c(1,1,1,2,3,3),nrow = 2,byrow = T))
}
dev.off()

#Figure S3 ####
pdf("figures/Fig. S3.pdf", width= 12, height =3.33)
print(polygon_cod_plot)
dev.off()
