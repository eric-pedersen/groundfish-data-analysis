library(assertthat)
library(dplyr)
library(readr)
#fishing data was obtained from "https://www.nafo.int/Data/STATLANT" (the NAFO
# STATLANT 21A Data Extraction Tool). Search parameters were for all species and 
# all countries' catches for all years in divisions 2J, 3K, and 3L.
fishing_data = read_csv("data/STATLANT21A_Extraction.csv",skip=10)
fishing_data = fishing_data %>%
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
assert_that(all(c(benthic_fish,pelagic_fish, invertebrates)%in%unique(fishing_data$species)))
assert_that(all(unique(fishing_data$species)%in%c(benthic_fish,pelagic_fish, invertebrates)))

fishing_summary = fishing_data %>%
  mutate(type = ifelse(species%in%benthic_fish, "benthic",
                       ifelse(species%in%pelagic_fish, "pelagic", "invertebrate")))%>%
  group_by(year,type)%>%
  summarise(catch = sum(catch))


# Catch effort data ####
# data is from NAFO catch effort time series 21B, downloaded from: https://www.nafo.int/Data/Catch-Statistics
nafo_effort_1980 = read_csv("data/NAFO effort data/NAFO21B-80-89.txt")
nafo_effort_1990 = read_csv("data/NAFO effort data/NAFO21B-90-99.txt")
nafo_effort_2000 = read_csv("data/NAFO effort data/NAFO21B-2000-09.txt")%>%
  rename(Catches = Month_NK) #fixing an inconsistency between time series column names
nafo_effort_2010 = read_csv("data/NAFO effort data/NAFO21B-2010-14.csv")%>%
  rename(Catches = Month_NK,GearCode = Gear,Divcode = AreaCode,Code = SpeciesEffort)
assert_that(all(names(nafo_effort_1980)==names(nafo_effort_1990)))
assert_that(all(names(nafo_effort_1980)==names(nafo_effort_2000)))
assert_that(all(names(nafo_effort_1980)==names(nafo_effort_2010)))

nafo_effort_all = rbind(nafo_effort_1980,nafo_effort_1990, nafo_effort_2000,nafo_effort_2010)%>%
  filter(Divcode %in%c(23,31,32),#filter only for 2J,3K,3L
         between(Year, 1981,2013), #filter for only those years in the data set
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
                        "benthic_fish",
                        ifelse(MainSpecies%in% c(20,21,30,31,41, 42, 49, 52,53,99),
                               "pelagic_fish", "invertebrates"))
         )

nafo_effort_summary = nafo_effort_all %>%
  group_by(Year, type)%>%
  summarize(effort= sum(effort,na.rm = T))

