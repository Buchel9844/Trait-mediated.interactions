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
library(plotly)
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
#---- Wainwright 2014 preparation
plant_code_aus <- read.csv("data/plant_code_aus.csv",
                           header = T, stringsAsFactors = F, sep=",",
                           na.strings = c("","NA"))

Wainwright2014.seed <-  read.csv("data/aus_rawdata/Wainwright2014_seed_aus.csv" )
nrow(Wainwright2014.seed)
Wainwright2014.seed <- Wainwright2014.seed %>%
  filter(!is.na(seed)) %>%
  filter( trt.water=="D" & reserve == "Perenjori") %>%
  mutate(year = 2014,
         plot = str_remove_all(plot, '\\.')) %>%
  unite("plot",c(plot,quadrant),sep="_",remove=T) %>%
  mutate(focal=case_when(focal == "A" ~ "ARCA",
                         focal == "T" ~ "TRCY",
                         focal == "W" ~ "WAAC",
                         focal == "H" ~ "HYGA")) %>%
  aggregate(seed  ~ plot + focal + year ,mean) 
str(Wainwright2014.seed)  

Wainwright2014 <-  read.csv("data/aus_rawdata/Wainwright2014_com_aus.csv" )
names(Wainwright2014) <- tolower(names(Wainwright2014))
Wainwright2014<- Wainwright2014 %>%
  filter( trt.water=="D" & reserve == "perenjori") %>%
  mutate(year = 2014) %>%
  unite("plot",c(plot,quadrant),sep="_",remove=T) %>%
  separate(col=species, into=c("genus","species")) %>%
  filter(!is.na(genus)) %>%
  right_join(plant_code_aus[,c("genus","species","final.code")], 
             by=c("genus","species"),
             relationship = "many-to-many") %>%
  filter(!is.na(count)) %>%
  mutate(focal=case_when(focal == "a" ~ "ARCA",
                         focal == "t" ~ "TRCY",
                         focal == "w" ~ "WAAC",
                         focal == "h" ~ "HYGA")) %>%
  aggregate(count ~ plot + focal + final.code + year,mean) %>%
  spread(final.code,count) %>%
  right_join(Wainwright2014.seed, relationship = "many-to-many") %>%
  filter(!is.na(seed))


Wainwright2014[is.na(Wainwright2014)] <- 0
view(Wainwright2014) 
table(Wainwright2014$focal)

#---- Martyn 2016 preparation
Martyn2016 <-  read.csv("data/aus_rawdata/Martyn2016_aus.csv" )
Martyn2016<- Martyn2016 %>%
  select(-"focal.ID")
str(Martyn2016)
table(Martyn2016$focal)

#---- Pastore 2017 preparation
Pastore2017_com <-  read.csv("data/aus_rawdata/Pastore2017_com_aus.csv" )
Pastore2017_seed <-  read.csv("data/aus_rawdata/Pastore2017_seed_aus.csv" )

Pastore2017_com <- Pastore2017_com %>% 
  filter(!is.na(.[,"ARCA"])) %>%
  mutate(plot.id.focal = paste0(.[,"plot"],.[,"id.focal"],.[,"focal"])) 
Pastore2017_com <- Pastore2017_com[!duplicated(Pastore2017_com$plot.id.focal),]

head(Pastore2017_seed)

#---- Merge
competition_aus <- bind_rows(Martyn2016,Wainwright2014)


write.csv(competition_aus,
          "data/competition_aus.csv")


#---- 2.1. Import the competitive data ----

competition_aus <- read.csv("data/competition_aus.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))
competition_aus$year <- 2016
head(competition_aus)

view(Wainwright2014)



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

focal.levels_aus <- levels(as.factor(competition_aus$focal))

ggplotly(ggplot(competition_aus, aes(x=seed,color=focal)) +geom_density() + xlim(0,250))

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

