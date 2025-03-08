
library(surveyjoin)
library(dplyr)
library(ggplot2)

cache_data()
load_sql_data()

data_version()

get_metadata()
get_shapefiles()
get_rawdata()

get_species()

#example single species data for region and coast
d <- get_data(common = "sablefish", years = 2013:2023, regions = "afsc")

d <- get_data(common = "arrowtooth flounder", years = 2003:2018)

#
data("spp_dictionary")

#
d3 <- get_data(common = c("Arrowtooth Flounder", "Pacific Ocean Perch", "Walleye Pollock", "Pacific Cod", "Pacific Halibut", "Rex Sole", "Flathead Sole", "Northern Rock Sole", "English Sole", "Sablefish", "Dover Sole", "Yelloweye Rockfish", "Shortraker Rockfish", "Lingcod"), years = 2003:2024, surveys = c("Aleutian Islands", "Gulf of Alaska", "NWFSC.Combo", "NWFSC.Shelf", "NWFSC.Hypoxia", "Triennial", "SYN QCS", "SYN HS", "SYN WCVI", "SYN WCHG"))

saveRDS(d3, "dfa_data_origplus_fullcoast_Feb2025.RDS")


