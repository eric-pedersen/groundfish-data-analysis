library(assertthat)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(vegan)
library(lme4)
library(ggplot2)
library(broom)

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
  summarise(temperature = mean(Temp_at_fishing,na.rm = T))

# AMO time series from: https://www.esrl.noaa.gov/psd/data/timeseries/AMO/
# I downloaded the "AMO smoothed, short (1948 to present)" series
amo_data = read_fwf("data/amon.us.data",
                    col_positions = fwf_positions(start = c(2,seq(9,108,by=9)),
                                                  end   = c(5,seq(14,113, by=9)),
                    col_names =c("year", 1:12)),skip = 1,n_max = 66)%>%
  filter(between(year, 1971, 2013))%>%
  gather(month, amo,-year)%>%
  mutate(month=as.numeric(month))%>%
  arrange(year,month)%>%
  mutate(amo_sum = cumsum(amo))%>%
  group_by(year)%>%
  summarise(amo =mean(amo),amo_sum = last(amo_sum))

ao_data = read_csv("data/ao_data.csv",skip = 1)%>%
  mutate(year = str_sub(Date, 1,4),
         year = as.numeric(year),
         month = as.numeric(str_sub(Date,5,-1)))%>%
  filter(between(year,1971,2013))%>%
  arrange(year,month)%>%
  mutate(ao_sum = cumsum(Value))%>%
  group_by(year)%>%
  summarize(ao = mean(Value),
            ao_sum = last(ao_sum))

# NAO time series downloaded from ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/nao_index.tim
# found via this website: http://www.cpc.ncep.noaa.gov/data/teledoc/nao.shtml
nao_data =read_fwf("data/nao_index.tim",skip=9, 
                   col_positions = fwf_positions(start=c(1,8,12),
                                                 end = c(4,9,16),
                                                 col_names =c("year","month", "nao")))%>%
  filter(between(year, 1971,2013))%>%
  arrange(year, month)%>%
  mutate(nao_sum = cumsum(nao))%>%
  group_by(year)%>%
  summarise(nao = mean(nao),nao_sum = last(nao_sum))

climate_data = nao_data %>%
  left_join(ao_data)%>%
  left_join(amo_data)%>%
  left_join(temperature_data)%>%
  filter(between(year, 1981,2013))

climate_princomp = princomp(select(climate_data,-year),cor = T)

climate_data$pc1 = climate_princomp$scores[,1]
climate_data$pc2 = climate_princomp$scores[,2]

climate_data_long =climate_data%>%
  gather(index, value,-year)


# Loading community composition data averaged over years
load("data/year_geom_means.Rdata")

community_mean_data = Year_Geom_Means
community_mean_data$Total = rowSums(community_mean_data)
community_mean_data$year = 1981:2013

community_driver_data  = community_mean_data %>%
  left_join(climate_data)%>%
  ungroup()%>%
  gather(species, biomass, ANARHICHAS_DENTICULATUS:UROPHYCIS_TENUIS)%>%
  group_by(species)%>%
  mutate(biomass_l = log10(biomass+0.01))%>%
  mutate_all(.funs = funs(lag =lag))%>%
  mutate(growth =biomass_l - biomass_l_lag)



model_autocor = lmer(growth~biomass_l_lag+
                 (1+biomass_l_lag|species),
               data=community_driver_data,na.action = na.exclude)

model_full = lmer(growth~biomass_l_lag+nao +ao+amo+temperature+
                 (1+biomass_l_lag+nao+ao+amo+temperature|species),
               data=community_driver_data,na.action = na.exclude)

driver_model_coefs = tidy(model_full, effect= "ran_modes")%>%
  rename(species= level)%>%
  filter(!term=="(Intercept)")%>%
  mutate(species =str_to_title(species),
         species = str_replace(species, "_", " "))
         
species_order = driver_model_coefs%>%
  filter(term =="biomass_l_lag")%>%
  arrange(estimate)

driver_model_coefs = driver_model_coefs%>%
  mutate(species = factor(species, levels = species_order$species))

ggplot(aes(species, estimate),data=driver_model_coefs)+
  geom_point()+
  facet_grid(.~term,scales = "free_x") +
  coord_flip()+
  geom_linerange(aes(ymin = estimate-2*std.error, ymax =estimate+2*std.error))+
  geom_hline(yintercept = 0, lty=2)

model_predict = community_driver_data%>%
  ungroup()%>%
  mutate(fit = 10^as.vector(predict(model_autocor)+biomass_l_lag))%>%
  select(year, species, fit)%>%
  arrange(year)%>%
  spread(species, fit)%>%
  select(-year)
