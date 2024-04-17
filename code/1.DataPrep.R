#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0. SET UP: pacakges adn create metadata ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 0.1. METADATE: Create meta data file to Summarise our data----
# 0.1. CREATE TEMPLATES
create_spice() # creates template from directory

# 0.2. EDIT TEMPLATES
#attributes.csv: This is where most of the user data entry will take place. 
# For each variable, its name, units, and a written description are filled in.
edit_attributes()

# access.csv: Includes a row for each file that was read in, 
# and documents the name of each file and its format.
edit_access()

# creators.csv: One row for each creator, and gives their affiliation,
# contact email, ORCID, etc
edit_creators() #

#biblio.csv: Citation information about the project, 
# as much or as little data as possible can be included, 
# but if things like bounding box coordinates are not included, 
# then when the website is generated there will not be a bounding box map generated
edit_biblio() 

#---- 0.2. Import packages----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages("HDInterval")
library("HDInterval")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("dplyr")
library(dplyr)
library(ggpubr)
library(ggplot2)
#rstan_options(auto_write = TRUE)
library(tidyr) #fill is part of tidyr
library(lme4)
library(car)
library(loo)
library(wesanderson) # for color palette
library(ggthemes) 
library(grid)
setwd("~/Documents/Projects/Facilitation_gradient")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Data from Caracoles, Donana National Park, Spain----
# Main collector Oscar Godoy and Nacho Bartomeus
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#---- 1.1. Import the competitive data ----
competition_spain <- read.csv("data/competition_spain.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))

species.neigh <- names(competition_spain)
species.neigh  <- species.neigh[!species.neigh %in% c("day","month","year","plot","subplot","focal","fruit","comments",
                                                      "observers", "site", "treatment")]

competition.to.keep <- competition_spain %>%  
  dplyr::select(all_of(c(species.neigh)))%>%
  mutate_at( species.neigh, as.numeric) %>%
  colSums(na.rm=T)
length(competition.to.keep[competition.to.keep == 0]) # check if any is 0

# change the format of the rows to numeric 
competition_spain[species.neigh] <- sapply(competition_spain[species.neigh],as.numeric)

# change na values to 0
competition_spain[is.na(competition_spain)] <- 0

focal.levels <- levels(as.factor(competition_spain$focal))

#---- 1.2. Read Teasorus ----
plant_code_spain <- read.csv( "data/plant_code_spain.csv",
                        header = T, stringsAsFactors = F, sep=",",
                        na.strings = c("","NA"))

competition.long.neigh_spain <- competition_spain %>%
  gather(any_of(plant_code_spain$code.plant), key="code.plant",value="abundance") %>%
  left_join(plant_code_spain) 

competition.long.focal_spain <- competition.long.neigh_spain %>%
  filter(focal == code.plant_spain) %>%
  rename("conspecific"="abundance") %>%
  dplyr::select(focal,conspecific,fruit,seed,plot,subplot,year) 

competition.long_spain <- competition.long.neigh_spain %>%
  aggregate(abundance ~  focal + fruit + seed + family + plot + subplot + year, sum )%>%
  spread(family, abundance) %>% # change this ti have grass /forb or native/exotic
  right_join(competition.long.focal_spain) 

head(competition.long_spain )  


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Data from Perenjori, WA, Australia----
# Main collector Trace Martin, Courtney Taylor, Margie Mayfield
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2.1. Import the competitive data ----

competition_aus <- read.csv("data/competition_aus_2016.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))

head(competition_aus)

species.neigh <- names(competition_aus)
species.neigh  <- species.neigh[!species.neigh %in% c("day","month","year","plot","subplot","seed","focal.ID",
                                                      "focal","fruit","comments",
                                                      "observers", "site", "treatment")]

competition.to.keep <- competition_aus %>%  
  dplyr::select(all_of(c(species.neigh)))%>%
  mutate_at( species.neigh, as.numeric) %>%
  colSums(na.rm=T)
length(competition.to.keep[competition.to.keep == 0]) # check if any is 0

# change the format of the rows to numeric 
competition_aus[species.neigh] <- sapply(competition_aus[species.neigh],as.numeric)

# change na values to 0
competition_aus[is.na(competition_aus)] <- 0

focal.levels <- levels(as.factor(competition_aus$focal))

#---- 2.2. Read Teasorus ----
plant_code_aus <- read.csv("data/plant_code_aus.csv",
                            header = T, stringsAsFactors = F, sep=",",
                            na.strings = c("","NA"))

to.remove.sp <- species.neigh[!species.neigh %in% plant_code_aus$code]

competition.long.neigh_aus <- competition_aus%>%
  select(!all_of(to.remove.sp)) %>%
  gather(any_of(plant_code_aus$code), key="code",value="abundance") %>%
  right_join(plant_code_aus) 
str(competition.long.neigh_aus)

competition.long.focal_aus <- competition.long.neigh_aus%>%
  filter(focal == final.code) %>%
  rename("conspecific"="abundance") %>%
  dplyr::select(focal,conspecific,seed,plot) 

competition.long_aus <- competition.long.neigh_aus %>%
  aggregate(abundance ~  focal +  seed + family + plot, sum )%>%
  spread(family, abundance) %>% 
  right_join(competition.long.focal_aus) 

head(competition.long_aus)  

