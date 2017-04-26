library(assertthat)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(vegan)
library(lme4)
library(ggplot2)
library(broom)
library(RcppRoll)

rm(list = ls())
#fishing data was obtained from "https://www.nafo.int/Data/STATLANT" (the NAFO
# STATLANT 21A Data Extraction Tool). Search parameters were for all species and 
# all countries' catches for all years in divisions 2J, 3K, and 3L.
catch_data = read_csv("data/STATLANT21A_Extraction.csv",skip=10)
catch_data = catch_data %>%
  transmute(year =Year, country=Country, species =Species, catch= `Catch ('000 Kg)`)%>%
  filter(year>1980)%>%
  group_by(species,year)%>%
  summarise(catch = sum(catch))%>%
  group_by(species)%>%
  filter(!all(catch<100))
  
benthic_fish = c("Greenland Halibut", "Atlantic Halibut" ,
                 "Winter Flounder", "Skates (ns)", 
                 "American Plaice", "Atlantic Redfishes (ns)", "Atlantic Cod",
                 "Witch Flounder", "Yellowtail Flounder", "American Angler",
                 "Beaked Redfish(deep-water)", 
                 "Dogfishes (ns)","Wolffishes (catfish) (ns)","Atlantic Wolffish", 
                 "Finfishes (ns)", "Flatfishes (ns)","Greenland Cod","Groundfishes (ns)",
                 "Haddock","Lampfishes (ns)","Pollock (saithe)","Yellowtail Flounder")

                 
pelagic_fish = c("Roughhead Grenadier", "Roundnose Grenadier","Atlantic Herring",
                 "Atlantic Mackerel","Capelin","Red Hake",
                 "Atlantic Salmon","Baird's Slickhead","Chars (ns)", "Large Sharks (ns)",
                 "Northern Bluefin Tuna", "Porbeagle","Swordfish","White Hake")
                 
invertebrates = c("(northern) Shortfin Squid", "American Lobster", "Whelks (ns)",
                  "Atlantic Rock Crab", "Marine Crabs (ns)","Bay Scallop",
                  "Northern Prawn","Aesop Shrimp","Sea Urchin" ,"Queen Crab",
                  "Icelandic Scallop","Marine Molluscs (ns)", "Pink (=pandalid) Shrimps",
                  "Sea Scallop","Surf Clam" )

#test if all species are in one of these lists
assert_that(all(c(benthic_fish,pelagic_fish, invertebrates)%in%unique(catch_data$species)))
assert_that(all(unique(catch_data$species)%in%c(benthic_fish,pelagic_fish, invertebrates)))

catch_summary = catch_data %>%
  mutate(type = ifelse(species%in%benthic_fish, "benthic",
                       ifelse(species%in%pelagic_fish, "pelagic", "invertebrate")))%>%
  group_by(year,type)%>%
  summarise(catch = sum(catch))


# Catch effort data ####
# data is from NAFO catch effort time series 21B, downloaded from: https://www.nafo.int/Data/Catch-Statistics
effort_1980 = read_csv("data/NAFO effort data/NAFO21B-80-89.txt")
effort_1990 = read_csv("data/NAFO effort data/NAFO21B-90-99.txt")
effort_2000 = read_csv("data/NAFO effort data/NAFO21B-2000-09.txt")%>%
  rename(Catches = Month_NK) #fixing an inconsistency between time series column names
effort_2010 = read_csv("data/NAFO effort data/NAFO21B-2010-14.csv")%>%
  rename(Catches = Month_NK,GearCode = Gear,Divcode = AreaCode,Code = SpeciesEffort)
assert_that(all(names(effort_1980)==names(effort_1990)))
assert_that(all(names(effort_1980)==names(effort_2000)))
assert_that(all(names(effort_1980)==names(effort_2010)))

effort_all = rbind(effort_1980,effort_1990, effort_2000,effort_2010)%>%
  rename(year = Year)%>%
  filter(Divcode %in%c(23,31,32),#filter only for 2J,3K,3L
         between(year, 1981,2013), #filter for only those years in the data set
         Code == 2 #filter only effort measured in days
         )%>%
  rowwise()%>%
  mutate(tonne_mean = recode(Tonnage,`1` = 25/2, `2` = 25, `3`= (50+150)/2,
                             `4`= (150+500)/2, `5`= (500+1000)/2, `6`= (1000+2000)/2,
                             `7` = 3000, `0`=0),#recodes as the midpoint of each tonnage category
         tonne_mean = ifelse(tonne_mean==0, NA, tonne_mean), #recodes zeros as missing tonnage data
         tonne_mean = ifelse(Country%in%c(2,3)&Tonnage==2,(25 +50)/2, tonne_mean),# Canada has a special code for tonnage 2; it divides the smallest tonnages into 2 sizes 
         total_days = sum(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep,Oct, Nov,Dec),
         effort = tonne_mean*total_days,
         type = ifelse(MainSpecies%in% c(1,2,3,6,10, 11,12,13,14,16,19,22, 29,54,59),
                        "benthic",
                        ifelse(MainSpecies%in% c(20,21,30,31,41, 42, 49, 52,53,99),
                               "pelagic", "invertebrate"))
         )

effort_summary = effort_all %>%
  group_by(year, type)%>%
  summarize(effort= sum(effort,na.rm = T))


fishing_summary = effort_summary %>%
  left_join(catch_summary)
fishing_summary_long = fishing_summary%>%
  gather(driver, value,-type,-year)

fishing_totals = fishing_summary %>%
  group_by(year) %>%
  summarize(effort= sum(effort), catch = sum(catch))



# Climate data ####
dfo_data = read_csv("data/DFO_Dataset.csv")

temperature_data = dfo_data %>%
  filter(Season!="Spring", Month %in%c(10,11,12),
         !is.na(Depth), !is.na(Temp_at_fishing))%>%
  select(Year, Month,Depth, Temp_at_fishing)%>%
  rename(year=Year, month=Month)%>%
  group_by(year)%>%
  summarise(temperature = mean(Temp_at_fishing,na.rm = T))%>%
  ungroup()%>%
  mutate(temp_mean = ifelse(year<1986,cummean(temperature),
                            roll_meanr(temperature,n = 5,align = "right")))

# AMO time series from: https://www.esrl.noaa.gov/psd/data/timeseries/AMO/
# I downloaded the "AMO smoothed, short (1948 to present)" series
amo_data = read_fwf("data/amon.us.data",
                    col_positions = fwf_positions(start = c(2,seq(9,108,by=9)),
                                                  end   = c(5,seq(14,113, by=9)),
                    col_names =c("year", 1:12)),skip = 1,n_max = 66)%>%
  filter(year<=2013)%>%
  gather(month, amo,-year)%>%
  mutate(month=as.numeric(month))%>%
  arrange(year,month)%>%
  group_by(year)%>%
  summarise(amo =mean(amo))%>%
  mutate(amo_mean = roll_mean(amo,n=5, align = "right",fill = NA))

ao_data = read_csv("data/ao_data.csv",skip = 1)%>%
  mutate(year = str_sub(Date, 1,4),
         year = as.numeric(year),
         month = as.numeric(str_sub(Date,5,-1)))%>%
  filter(year<=2013)%>%
  arrange(year,month)%>%
  group_by(year)%>%
  summarize(ao = mean(Value))%>%
  mutate(ao_mean = roll_mean(ao,n=5, align = "right",fill = NA))


# NAO time series downloaded from ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/nao_index.tim
# found via this website: http://www.cpc.ncep.noaa.gov/data/teledoc/nao.shtml
nao_data =read_fwf("data/nao_index.tim",skip=9, 
                   col_positions = fwf_positions(start=c(1,8,12),
                                                 end = c(4,9,16),
                                                 col_names =c("year","month", "nao")))%>%
  filter(year<=2013)%>%
  arrange(year, month)%>%
  group_by(year)%>%
  summarise(nao = mean(nao))%>%
  mutate(nao_mean = roll_mean(nao,n=5, align = "right",fill = NA))


climate_data = nao_data %>%
  left_join(ao_data)%>%
  left_join(amo_data)%>%
  left_join(temperature_data)%>%
  filter(between(year, 1981,2013))%>%
  select(-nao, -ao,-amo,-temperature)

climate_princomp = princomp(select(climate_data,-year),cor = T)


climate_data$pc1 = climate_princomp$scores[,1]
climate_data$pc2 = climate_princomp$scores[,2]

climate_data_long =climate_data%>%
  select(-pc1,-pc2)%>%
  gather(index, value,-year)%>%
  mutate(index = recode(index,nao_mean = "NAO (5-year avg)",
                        ao_mean = "AO (5-year avg)",
                        amo_mean = "AMO (5-year avg)",
                        temp_mean = "mean bottom\ntemperature (C, 5-year avg)"),
         index = factor(index, levels = c("NAO (5-year avg)","AO (5-year avg)",
                                          "AMO (5-year avg)","mean bottom\ntemperature (C, 5-year avg)")))


climate_plot = ggplot(data=climate_data_long, aes(x= year, value))+
  geom_point() +
  geom_line()+ 
  facet_grid(index~., scales="free_y")+
  theme_bw()

climate_ord_plot = ggplot(data= climate_data, aes(x=pc1,y=pc2,label=year))+
  geom_path()+
  geom_text()+
  labs(x = "Climate principal component 1",
       y = "Climate principal component 2")+
  theme_bw()

fishing_plot = ggplot(data= fishing_summary_long, aes(x = year, y= value,color=type))+
  geom_point() +
  geom_line()+
  scale_color_brewer(palette="Set1")+
  facet_grid(driver~.,scales = "free_y")+
  theme_bw()


ggsave("figures/climate_ts.pdf",climate_plot, width= 6, height=8) 
ggsave("figures/climate_ord.pdf",climate_ord_plot, width= 6, height=6) 
ggsave("figures/fishing_ts.pdf", fishing_plot, width= 6, height=6)

write.csv(fishing_summary, "data/fishing_effort_data.csv",row.names = F)
write.csv(climate_data, "data/climate_data.csv",row.names = F)

driver_data = data_frame(year=climate_data$year,
                         `Benthic fishing effort\n(megatonne-days)`=fishing_summary$effort[fishing_summary$type=="benthic"]/1e6,
                         `Bottom temperature\n(5-year running mean)`=climate_data$temp_mean,
                         `Climate\n(first principle component)` = climate_data$pc1
                    )

driver_data = driver_data%>%
  gather(index, value, -year)

driver_labels = driver_data %>%
  group_by(index)%>%
  summarise(n=n())%>%
  mutate(label =LETTERS[1:3],
         year = 1981,
         value = c(15, 2.5,3))

driver_plot = ggplot(driver_data,aes(year, value)) + 
  geom_line()+
  geom_text(data =driver_labels,aes(label = label))+
  facet_grid(index~.,scales = "free_y",switch = "y")+
  theme_bw()+
  labs(y=NULL)+
  theme(strip.placement = "outside",strip.background = element_blank(),panel.grid = element_blank()
          )
  


ggsave(filename = "figures/Fig. XX.pdf", driver_plot,height=8,width=5)



