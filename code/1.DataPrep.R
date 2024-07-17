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
library(pals) # for lots of colors
setwd("~/Documents/Projects/Facilitation_gradient")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Data from Caracoles, Donana National Park, Spain----
# Main collector Oscar Godoy and Nacho Bartomeus
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1.1. Read Teasorus ----
plant_code_spain <- read.csv( "data/plant_code_spain.csv",
                              header = T, stringsAsFactors = F, sep=",",
                              na.strings = c("","NA"))



#---- 2.2. Identify most abundance species ----
 
abundance_spain <- read.csv("data/abundance_spain.csv",
                            header = T,stringsAsFactors = F, sep=",",
                            na.strings=c("","NA"))

abundance_spain.summary <- abundance_spain %>%
  aggregate(individuals ~ year + species, sum)

species.list.to.keep <- c("BEMA","CETE","CHFU",
                          "CHMI","HOMA","LEMA","MESU","PAIN","PLCO",
                          "POMA","POMO","SASO","SCLA","SPRU")
# Amaranthaceae, Gentianaceae, Asteraceae, 
# Poaceae, Fabaceae,Plantaginaceae, Amaranthaceae, Caryophyllaceae
final.species.list <- c("BEMA","CETE","CH.sp",
                        "HOMA","LEMA","ME.sp","PAIN","PLCO",
                        "PO.sp","SASO","SCLA","SPRU","rare")

# regroup POMA and POMO under Polypogon
# regroup MEEL and MESU under Melilotus ? 
# regroup CHFU and CHMI under Chamaemelum? Or put CHMI under rare 
species.list.to.rare <- c("COSQ","ACHI","ANAR","FRPU","LYTR","MEEL",
                          "MEPO","PUPA","RAPE","SOAS","SUSP")

competition.long.neigh_spain.sum <- competition.long.neigh_spain %>%
  aggregate(fruit ~ year + focal, length)
# PLCO and SCLA missing 2018
# PUPA has good data

#---- 2.3. Import competition data -----
competition.spain <- read.csv("data/spain_rawdata/competition.csv")
# maybe can add 2022 and 2023 later
str(competition.spain)
species.neigh <- names(competition.spain)
species.neigh  <- species.neigh[!species.neigh %in% c("day","month","year","plot","subplot","focal","fruit","comments",
                                                      "observers", "site", "treatment")]
competition.to.keep <- competition.spain %>%  
  dplyr::select(all_of(c(species.neigh)))%>%
  mutate_at( species.neigh, as.numeric) %>%
  colSums(na.rm=T)
length(competition.to.keep[competition.to.keep == 0]) # check if any is 0

# change the format of the rows to numeric 
competition.spain[species.neigh] <- sapply(competition.spain[species.neigh],as.numeric)

# change na values to 0
competition.spain[is.na(competition.spain)] <- 0


competition.spain_long <- competition.spain %>%
  gather(any_of(plant_code_spain$code.plant), key="code.plant",
         value="abundance") %>%
  left_join(plant_code_spain, by="code.plant") %>%
  left_join(plant_code_spain %>%
              select(code.plant,code.analysis) %>%
              rename(focal="code.plant",focal.analysis="code.analysis"),
            by="focal")  %>%
  select(day,month,year,plot,subplot,focal.analysis,fruit,seed,code.analysis,
         abundance) %>%
  aggregate(abundance ~ day + month + year + plot + subplot + focal.analysis +
              fruit + seed + code.analysis, sum) %>%
  spread(code.analysis, abundance)

comp.plot <- competition.spain_long %>%
  aggregate(fruit ~ focal.analysis + year, length) %>%
  ggplot(aes(y=fruit,x=as.factor(year),fill=as.factor(focal.analysis))) +
  geom_bar(stat="identity", position="dodge",
           color="black") +
  theme_bw() +
  scale_fill_manual(values=unname(kelly())) +
  labs(y="number of observations", 
       x="year",fill="focal",
       title="Number of observation of each focal for the competitive experirement in Spain") 
comp.plot

ggsave(paste0("figures/supp/observations.spain.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       comp.plot)

#---- 2.2. Seed production -----

numb.seed.spain <- read.csv("data/spain_rawdata/number_seed.csv")
str(numb.seed.spain)
numb.seed.spain <- numb.seed.spain %>%
  left_join(plant_code_spain %>%
              select(code.plant,code.analysis) %>%
              rename(focal="code.plant",focal.analysis="code.analysis"),
            by="focal") %>%
  mutate(fruits.panicle = ifelse(is.na(fruits.panicle),
                                 1,
                                 fruits.panicle ),
         fruits.panicle = ifelse(fruits.panicle==seeds.fruit,
                                 1,
                                 fruits.panicle )) %>%
  mutate(total.seeds = fruits.panicle*seeds.fruit) %>%
  filter(total.seeds < 3000)

numb.seed.spain.plot <- numb.seed.spain %>%
  ggplot(aes(y=total.seeds,x=as.factor(focal.analysis),
             fill=as.factor(focal.analysis))) +
  geom_boxplot(color="black") +
  labs(x="focal",y="seeds per individual",color="focal",
       title="Number of seeds per fruit of individuals of each focal for the competitive experiment in Spain") +
  scale_fill_manual(values=unname(kelly())) +
  theme_bw() 
numb.seed.spain.plot

ggsave(paste0("figures/supp/Seed.per.flower.spain.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       numb.seed.spain.plot)

#---- 2.3. Join seed and interactions -----

set.seed(1616)

for(i in c(1:nrow(competition.spain_long))){
  year.i <- competition.spain_long[i,"year"] 
  focal.i <- competition.spain_long[i,"focal.analysis"]
  plot.i <- competition.spain_long[i,"plot"] 
  subplot.i <- competition.spain_long[i,"subplot"] 
  seeds.fruit.i <- NA

  
  seeds.fruit.i <- abs(rnorm(1,mean = mean(numb.seed.spain[which(numb.seed.spain$year %in% year.i &
                                                numb.seed.spain$focal.analysis == focal.i &
                                                numb.seed.spain$plot == plot.i &
                                                numb.seed.spain$subplot == subplot.i),
                                  "total.seeds"], na.rm=T),
                         sd = sd(numb.seed.spain[which(numb.seed.spain$year %in% year.i &
                                                           numb.seed.spain$focal.analysis == focal.i &
                                                           numb.seed.spain$plot == plot.i &
                                                         numb.seed.spain$subplot == subplot.i),
                                                   "total.seeds"], na.rm=T)))
  
  if(is.na(seeds.fruit.i)){
    seeds.fruit.i <- abs(rnorm(1,mean(numb.seed.spain[which(numb.seed.spain$year  %in% year.i &
                                                        numb.seed.spain$focal.analysis == focal.i &
                                                        numb.seed.spain$plot == plot.i),
                                                "total.seeds"], na.rm=T),
                               sd = sd(numb.seed.spain[which(numb.seed.spain$year  %in% year.i &
                                                         numb.seed.spain$focal.analysis == focal.i &
                                                           numb.seed.spain$plot == plot.i),
                                                 "total.seeds"], na.rm=T)))
    
    
  }
  if(is.na(seeds.fruit.i)){
    seeds.fruit.i <- abs(rnorm(1,mean(numb.seed.spain[which(numb.seed.spain$year %in% year.i &
                                                        numb.seed.spain$focal.analysis == focal.i),
                                                "total.seeds"], na.rm=T),
                               sd = sd(numb.seed.spain[which(numb.seed.spain$year %in% year.i &
                                                               numb.seed.spain$focal.analysis == focal.i ),
                                                 "total.seeds"], na.rm=T)))
    
  }
  if(is.na(seeds.fruit.i)){
    seeds.fruit.i <- abs(rnorm(1,mean(numb.seed.spain[which(numb.seed.spain$focal.analysis == focal.i),
                                                      "total.seeds"], na.rm=T),
                               sd = sd(numb.seed.spain[which(numb.seed.spain$focal.analysis == focal.i ),
                                                       "total.seeds"], na.rm=T)))
    
  }
  
  competition.spain_long$seed[i] <-  seeds.fruit.i*competition.spain_long[i,"fruit"]
}

str(competition.spain_long)
focal.comp.seed.plot <- competition.spain_long %>%
  filter(seed < 5000) %>%
  ggplot(aes(y=seed,x=as.factor(focal.analysis),
             fill=as.factor(focal.analysis))) +
  geom_boxplot(color="black") +
  labs(color="focal",
       fill="",
       x="focal",
       y="total seeds peer individual") + 
  scale_fill_manual(values=unname(kelly())) +
  theme_bw()  
focal.comp.seed.plot 

ggsave(paste0("figures/supp/Seed.per.ind.spain.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       focal.comp.seed.plot)

write.csv(competition.spain_long,
          file=paste0(home.dic,"results/competition.spain_long.csv"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Data from Perenjori, WA, Australia----
# Main collector Trace Martin, Courtney Taylor, Margie Mayfield
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Wainwright 2014 preparation ----
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
  left_join(plant_code_aus[,c("genus","species","final.code")], 
             by=c("genus","species"),
             relationship = "many-to-many") %>%
  filter(!is.na(count)) %>%
  mutate(focal=case_when(focal == "a" ~ "ARCA",
                         focal == "t" ~ "TRCY",
                         focal == "w" ~ "WAAC",
                         focal == "h" ~ "HYGA")) %>%
  aggregate(count ~ plot + focal + final.code + year,mean) %>%
  spread(final.code,count) %>%
  left_join(Wainwright2014.seed, relationship = "many-to-many") %>%
  filter(!is.na(seed))


Wainwright2014[is.na(Wainwright2014)] <- 0
view(Wainwright2014) 
head(Wainwright2014) 
table(Wainwright2014$focal)

#---- Martyn 2016 preparation----
Martyn2016 <-  read.csv("data/aus_rawdata/Martyn2016_aus.csv" )
Martyn2016<- Martyn2016 %>%
  select(-"focal.ID") %>%
  gather(any_of(plant_code_aus$code), key="code", value="count") %>%
  select(plot,seed,focal,code,count) %>%
  left_join(plant_code_aus[,c("code","genus",
                               "species","final.code")], 
             by=c("code"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + focal + final.code + seed,sum) %>%
  spread(final.code,count) 

#str(Martyn2016)
#table(Martyn2016$focal)
#head(Martyn2016)


#---- Pastore 2017 preparation----
Pastore2017_com <-  read.csv("data/aus_rawdata/Pastore2017_com_aus.csv" )
Pastore2017_seed <-  read.csv("data/aus_rawdata/Pastore2017_seed_aus.csv" )

Pastore2017 <- Pastore2017_com %>% 
  left_join(Pastore2017_seed, by=c("plot","focal","id.focal"),
             relationship = "many-to-many") %>%
  filter(!is.na(seed)) %>%
  gather(any_of(plant_code_aus$code), key="code", value="count") %>%
  select(plot,seed,focal,code,count,id.focal) %>%
  left_join(plant_code_aus[,c("code","genus",
                               "species","final.code")], 
             by=c("code"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + id.focal + focal + final.code + seed,sum) %>%
  spread(final.code,count) 

#view(Pastore2017)
#head(Pastore2017)
#table(Pastore2017$focal)

#---- Sevenello 2022 preparation----

Sevenello2022_com <- read.csv("data/aus_rawdata/Sevenello2022_com_aus.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))
Sevenello2022_seed <- read.csv("data/aus_rawdata/Sevenello2022_seed_aus.csv",
                         header = T,stringsAsFactors = F, sep=",",
                         na.strings=c("","NA"))
Sevenello2022 <- Sevenello2022_com %>% 
  left_join(Sevenello2022_seed,relationship = "many-to-many") %>%
  dplyr::filter(treatment %in% c("OP")) %>%
  filter(!is.na(seed)) %>%
  gather(any_of(plant_code_aus$code), key="code", value="count") %>%
  select(plot,seed,focal,code,count) %>%
  left_join(plant_code_aus[,c("code","genus",
                               "species","final.code")], 
             by=c("code"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + focal + final.code + seed,sum) %>%
  spread(final.code,count) 
  
#str(Sevenello2022)
#head(Sevenello2022)
#table(Sevenello2022$focal)


#---- Taylor 2023 preparation----
load("data/aus_rawdata/Taylor2023_com_aus.RData")
Taylor2023_com <- bind_rows(neighbourhoods[["BOCL"]] %>% as.data.frame(),
                            neighbourhoods[["BOOP"]] %>% as.data.frame(),
                            neighbourhoods[["CACL"]] %>% as.data.frame(),
                            neighbourhoods[["CAOP"]] %>% as.data.frame(),
                            neighbourhoods[["PECL"]] %>% as.data.frame(),
                            neighbourhoods[["PEOP"]] %>% as.data.frame())
str(Taylor2023_com)
Taylor2023_seed <- read.csv("data/aus_rawdata/Taylor2023_seed_aus.csv",
                               header = T,stringsAsFactors = F, sep=",",
                               na.strings=c("","NA"))
str(Taylor2023_seed ) 
Taylor2023 <- Taylor2023_com %>% 
  rename("focal" ="SPECIES" ) %>%
  left_join(Taylor2023_seed,relationship = "many-to-many") %>%
  rename("plot" = "unique_id") %>%
  filter(!is.na(seed)) %>%
  gather(any_of(c(plant_code_aus$name)[!is.na(c(plant_code_aus$name))]), 
         key="species", value="count") %>%
  select(plot,seed,focal,code,count) %>%
  left_join(plant_code_aus[,c("code","genus",
                               "species","final.code")], 
             by=c("code"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + focal + final.code + seed,sum) %>%
  spread(final.code,count) %>%
  left_join(plant_code_aus[,c("name","final.code")] %>%
               rename("focal" ="name" ) %>%
               rename("final.focal" ="final.code" ), 
             by=c("focal"),
             relationship = "many-to-many") %>%
  select(-"focal") %>%
  rename("focal" ="final.focal" )

#str(Taylor2023)
#head(Taylor2023)
#table(Taylor2023$focal)

# summary table

summary_table_aus <- bind_rows(as.data.frame(table(Martyn2016$focal)) %>%
  mutate(year=2016),
  as.data.frame(table(Wainwright2014$focal)) %>%
    mutate(year=2014),
  as.data.frame(table(Pastore2017$focal)) %>%
    mutate(year=2017),
  as.data.frame(table(Taylor2023$focal)) %>%
    mutate(year=2023),
  as.data.frame(table(Sevenello2022$focal)) %>%
    mutate(year=2022)
    ) %>%
  rename("focal"="Var1") %>%
  mutate(focal = case_when(focal="VERO"~"GORO",
                           T~focal))
ggplot(summary_table_aus[which(summary_table_aus$Freq>50),],
       aes(y=Freq,x=year,color=focal,fill=focal)) +
  geom_bar(stat="identity",position="dodge") +
  geom_text(aes(label=focal),
            position=position_dodge(width=1),
           #vjust=0.3,
            angle=90)

#---- Merge
competition_aus <- bind_rows(Martyn2016,Wainwright2014)


write.csv(competition_aus,
          "data/competition_aus.csv")

