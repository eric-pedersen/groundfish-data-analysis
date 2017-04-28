
#Clear workspace
rm(list=ls(all=TRUE))

#Loading data and packages ####
library(plyr)
library(dplyr)
library(bootstrap)
library(RColorBrewer)
library(plotrix)
require(FD)
library(ggplot2)
library(viridis)
library(tidyr)


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


# Calculate functional diversities ####
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

# Calculating weighted mean trait data for total community ####
vertical_total_data = sum_traits_by_biomass("vertical.position",Traits,
                                      Year_Geom_Means,discrete = T )
vertical_total_data$trait_value = factor(vertical_total_data$trait_value, 
                                   levels =c("bathypelagic","bathydemersal",
                                             "demersal","benthopelagic",
                                             "pelagic-oceanic"))
food_total_data  = sum_traits_by_biomass("Food.Items",Traits,
                                   Year_Geom_Means,discrete = T )
food_total_data$trait_value = factor(food_total_data$trait_value,
                               levels = c("Small Benthivore","Medium Benthivore",
                                          "Large Benthivore", "PlankPiscivore",
                                          "Piscivore"))
agg_total_data  = sum_traits_by_biomass("Aggregation",Traits,
                                   Year_Geom_Means,discrete = T )
agg_total_data$trait_value = factor(agg_total_data$trait_value,
                               levels = c("irregular","none",
                                          "rare", "schools",
                                          "shoal"))

double_total_data = sum_traits_by_biomass("double.time",Traits,
                                    Year_Geom_Means,discrete = F )
length_total_data = sum_traits_by_biomass("length",Traits,
                                    Year_Geom_Means,discrete = F )
trophic_total_data = sum_traits_by_biomass("trophic.level",Traits,
                                    Year_Geom_Means,discrete = F )

# Calculating weighted mean trait data for rare species community ####
vertical_rare_data = sum_traits_by_biomass("vertical.position",Traits,
                                            Year_Geom_Means_rare,discrete = T )
vertical_rare_data$trait_value = factor(vertical_rare_data$trait_value, 
                                         levels =c("bathypelagic","bathydemersal",
                                                   "demersal","benthopelagic",
                                                   "pelagic-oceanic"))
food_rare_data  = sum_traits_by_biomass("Food.Items",Traits,
                                         Year_Geom_Means_rare,discrete = T )
food_rare_data$trait_value = factor(food_rare_data$trait_value,
                                     levels = c("Small Benthivore","Medium Benthivore",
                                                "Large Benthivore", "PlankPiscivore",
                                                "Piscivore"))
agg_rare_data  = sum_traits_by_biomass("Aggregation",Traits,
                                        Year_Geom_Means_rare,discrete = T )
agg_rare_data$trait_value = factor(agg_rare_data$trait_value,
                                    levels = c("irregular","none",
                                               "rare", "schools",
                                               "shoal"))

double_rare_data = sum_traits_by_biomass("double.time",Traits,
                                          Year_Geom_Means_rare,discrete = F )
length_rare_data = sum_traits_by_biomass("length",Traits,
                                          Year_Geom_Means_rare,discrete = F )
trophic_rare_data = sum_traits_by_biomass("trophic.level",Traits,
                                           Year_Geom_Means_rare,discrete = F )


# Creating plots ####
#Colours and eras used for plotting 
ColV<-brewer.pal(9,"Set1")
Eras<-c(1990,1995)

pdf("figures/Fig. 4.pdf",width=6, height=5)
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

func_div_data = data.frame(Year = 1981:2013, FDis = Func_Div$FDis, 
                           FDis_rare= Func_Div_rare$FDis)
#save the functional dispersion code so it doesn't have to be rerun all the time
write.csv(func_div_data,"data/func_div_data.csv",row.names = F)


#Plotting mean trait values for the total community ####
discrete_plot = list(aes(Year, proportion, color=trait_value,
                         group = paste(trait_value,Year<1995)),
                     geom_line(size=2), theme_bw(12),
                     geom_vline(xintercept = Eras,linetype=2),
                     scale_y_continuous(limits = c(0,1),expand=c(0,0)),
                     theme(legend.position = "bottom"),
                     guides(col = guide_legend(nrow = 3,byrow = T)))

continuous_plot = list(aes(Year, value,
                           group = Year<1995),
                       geom_line(size=2), theme_bw(12),
                       geom_vline(xintercept = Eras,linetype=2))

agg_total_plot = ggplot(data= agg_total_data) + 
  discrete_plot+
  scale_color_viridis("Aggregation",discrete = T)+
  annotate(geom="text",label = "A", x = 1982, y=0.9,size=5)

food_total_plot = ggplot(data= food_total_data) + 
  discrete_plot+
  scale_color_viridis("Food\nniche",discrete = T,option = "inferno")+
  annotate(geom="text",label = "B", x = 1982, y=0.9,size=5)

vertical_total_plot = ggplot(data= vertical_total_data) + 
  discrete_plot+
  scale_color_viridis("Vertical\nposition",discrete = T,option = "plasma")+
  annotate(geom="text",label = "C", x = 1982, y=0.9,size=5)


double_total_plot =  ggplot(data= double_total_data) + 
  continuous_plot+
  scale_y_continuous("Mean doubling time (years)")+
  annotate(geom="text",label = "D", x = 1982, y=9.5,size=5)

trophic_total_plot =  ggplot(data= trophic_total_data) + 
  continuous_plot+
  scale_y_continuous("Mean trophic level")+
  annotate(geom="text",label = "E", x = 1982, y=4.15,size=5)

length_total_plot =  ggplot(data= length_total_data) + 
  continuous_plot+
  scale_y_continuous("Mean maximum body length (cm)")+
  annotate(geom="text",label = "F", x = 1982, y=140,size=5)

pdf("figures/Fig. S2.pdf",width=12, height=6)

PlotMultipleGgplotObjs(agg_total_plot, food_total_plot, vertical_total_plot, double_total_plot,trophic_total_plot,length_total_plot,
                       layout = matrix(c(1:3,1:3, 1:3, 4:6,4:6),ncol = 3,byrow = T))
dev.off()


#Plotting mean trait values for the rare species community ####

agg_rare_plot = ggplot(data= agg_rare_data) + 
  discrete_plot+
  scale_color_viridis("Aggregation",discrete = T)+
  annotate(geom="text",label = "A", x = 1982, y=0.9,size=5)

food_rare_plot = ggplot(data= food_rare_data) + 
  discrete_plot+
  scale_color_viridis("Food\nniche",discrete = T,option = "inferno")+
  annotate(geom="text",label = "B", x = 1982, y=0.9,size=5)

vertical_rare_plot = ggplot(data= vertical_rare_data) + 
  discrete_plot+
  scale_color_viridis("Vertical\nposition",discrete = T,option = "plasma")+
  annotate(geom="text",label = "C", x = 1982, y=0.9,size=5)


double_rare_plot =  ggplot(data= double_rare_data) + 
  continuous_plot+
  scale_y_continuous("Mean doubling time (years)", 
                     limits = range(double_total_data$value))+
  annotate(geom="text",label = "D", x = 1982, y=9.5,size=5)

trophic_rare_plot =  ggplot(data= trophic_rare_data) + 
  continuous_plot+
  scale_y_continuous("Mean trophic level")+
  annotate(geom="text",label = "E", x = 1982, y=3.95,size=5)

length_rare_plot =  ggplot(data= length_rare_data) + 
  continuous_plot+
  scale_y_continuous("Mean maximum body length (cm)", 
                     limits = c(77.66, 145))+
  annotate(geom="text",label = "F", x = 1982, y=140,size=5)

pdf("figures/Fig. S3.pdf",width=12, height=6)

PlotMultipleGgplotObjs(agg_rare_plot, food_rare_plot, vertical_rare_plot, double_rare_plot,trophic_rare_plot,length_rare_plot,
                       layout = matrix(c(1:3,1:3, 1:3, 4:6,4:6),ncol = 3,byrow = T))
dev.off()