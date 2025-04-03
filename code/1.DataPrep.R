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
library(tidyverse)
#rstan_options(auto_write = TRUE)
library(tidyr) #fill is part of tidyr
library(lme4)
library(car)
library(loo)
library(wesanderson) # for color palette
library(ggthemes) 
library(grid)
library(pals) # for lots of colors
library(FactoMineR)
#install.packages("factoextra")
library(factoextra)
library(vegan)
setwd("~/Documents/Projects/Facilitation_gradient")


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
# 1. Data from Caracoles, Donana National Park, Spain----
# Main collector Oscar Godoy and Nacho Bartomeus
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.1. Read Teasorus ----
plant_code_spain <- read.csv( "data/spain_rawdata/plant_code_spain.csv",
                              header = T, stringsAsFactors = F, sep=",",
                              na.strings = c("","NA"))


#---- 2.2. Identify most abundance species ----
 
abundance_spain <- read.csv("data/spain_rawdata/abundance.csv",
                            header = T,stringsAsFactors = F, sep=",",
                            na.strings=c("","NA"))

abundance_spain.summary <- abundance_spain %>%
  aggregate(individuals ~ year + species, sum) %>%
  spread(species,individuals)

species.list.to.keep <- c("BEMA","CETE","CHFU",
                          "CHMI","HOMA","LEMA","MESU","PAIN","PLCO",
                          "POMA","POMO","SASO","SCLA","SPRU")
# Amaranthaceae, Gentianaceae, Asteraceae, 
# Poaceae, Fabaceae,Plantaginaceae, Amaranthaceae, Caryophyllaceae
final.species.list.spain <- c("BEMA","CETE","CHFU",
                        "HOMA","LEMA","ME.sp","PAIN","PLCO",
                        "PO.sp","SASO","SCLA","SPRU")
head(abundance_spain)
zscore.fct <- function(x){
  (x-mean(x))/var(x)
}
abundance_spain.summary <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  filter( code.analysis %in% final.species.list.spain ) %>%
  aggregate(individuals ~ code.analysis + year + plot + subplot , sum) %>%
  mutate(individuals =individuals/10000, # number of individuals per cm squarre
         com_id = paste(plot,subplot,sep="_")) %>%
  rename("species"="code.analysis") %>%
  group_by(year) %>%
  mutate(count.zscore = zscore.fct(individuals)) %>%
  ungroup()

# regroup POMA and POMO under Polypogon
# regroup MEEL and MESU under Melilotus ? 
# regroup CHFU and CHMI under Chamaemelum? Or put CHMI under rare 
species.list.to.rare <- c("COSQ","ACHI","ANAR","FRPU","LYTR","MEEL",
                          "MEPO","PUPA","RAPE","SOAS","SUSP")

# PLCO and SCLA missing 2018
# PUPA has good data
abundance.plot <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  filter( code.analysis %in% final.species.list.spain ) %>%
  aggregate(individuals ~ code.analysis + year + plot + subplot , sum) %>%
  mutate(individuals =individuals/100,
         com_id = paste(plot,subplot,sep="_")) %>%
  ggplot(aes(y=individuals,x=as.character(year),
             fill=as.factor(code.analysis),
             group=as.factor(code.analysis),
             color=as.factor(code.analysis))) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",size=2) +
  stat_summary(fun.y = mean,
               geom = "line",size=1) +
  theme_bw() +
  scale_y_log10() +
  scale_color_manual(values=unname(kelly())) +
  scale_fill_manual(values=unname(kelly())) +
  labs(y="number of individual per centimeter", 
       x="year",fill="species",color="species",
       title="Number of individuals of each focal for abundance observations in spain")  +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=28),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=20, angle=66, hjust=1),
         axis.text.y= element_text(size=20),
         axis.title.x= element_text(size=24),
         axis.title.y= element_text(size=24),
         title=element_text(size=16))
abundance.plot

ggsave(paste0("figures/abundance.spain.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       abundance.plot)

str(abundance_spain.summary)
abundance_spain.summary %>%
  #filter()%>%
  ggplot(aes(y=CHFU, x=SCLA)) +
  geom_point() +
  geom_smooth() +
  coord_cartesian(xlim=c(0,100),ylim=c(0,100))


abundance_spain.summary.df <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  filter( code.analysis %in% final.species.list.spain ) %>%
  aggregate(individuals ~ code.analysis + year + plot + subplot +species +family+Native.sp, sum) %>%
  group_by(code.analysis,species,family,Native.sp) %>%
  summarise(count = n(),
            mean = mean(individuals, na.rm = TRUE)/4, 
            sd = sd(individuals, na.rm = TRUE)/4)
write.csv(abundance_spain.summary.df,
          "results/abundance_spain.summary.df.csv")
head(abundance_spain.summary.df)
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
              dplyr::select(code.plant,code.analysis) %>%
              rename(focal="code.plant",focal.analysis="code.analysis"),
            by="focal")  %>%
  dplyr::select(day,month,year,plot,subplot,focal.analysis,fruit,seed,code.analysis,
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

#---- 2.4. Seed production -----

numb.seed.spain <- read.csv("data/spain_rawdata/number_seed.csv")
str(numb.seed.spain)
numb.seed.spain <- numb.seed.spain %>%
  left_join(plant_code_spain %>%
              dplyr::select(code.plant,code.analysis) %>%
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

#---- 2.5. Join seed and interactions -----

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
#---- 2.6. Seed germination and survival -----

seed_germination_spain <- read.csv(paste0("data/spain_rawdata/germination_2015.csv"),
                                   header = T,stringsAsFactors = F, sep=",",
                                   na.strings=c("","NA")) %>%
  bind_rows(read.csv(paste0("data/spain_rawdata/germination_2021.csv"),
                     header = T,stringsAsFactors = F, sep=",",
                     na.strings=c("","NA"))  %>%
              tidyr::gather(any_of(plant_code_spain$code.plant),
                            key="focal.code",value="germination"))%>%
  bind_rows(read.csv(paste0("data/spain_rawdata/germination_2020.csv"),
                     header = T,stringsAsFactors = F, sep=",",
                     na.strings=c("","NA")) %>%
                          tidyr::gather(any_of(plant_code_spain$code.plant),
                                        key="focal.code",value="germination")
  ) %>%
  dplyr::select(year,focal.code,focal,germination,survival) %>%
  rename("code.plant" ="focal.code") %>%
  left_join(plant_code_spain %>%
              dplyr::select(code.plant,code.analysis)) %>%
  mutate(germination = case_when(year ==2015 ~ germination/100,
                                 T ~ germination/50)) %>%
  group_by(code.analysis) %>% 
  mutate(g.mean = mean(germination, na.rm = T),
         g.sd = sd(germination, na.rm = T)) %>%
  dplyr::select(code.analysis,g.mean,g.sd) %>%
  unique() %>%
  ungroup() %>%
  as.data.frame()

#view(seed_germination_spain)

seed_survival_spain <- read.csv(paste0("data/spain_rawdata/germination_2015.csv"),
                                header = T,stringsAsFactors = F, sep=",",
                                na.strings=c("","NA")) %>%
  bind_rows(read.csv(paste0("data/spain_rawdata/seed_survival_2020.csv"),
                     header = T,stringsAsFactors = F, sep=",",
                     na.strings=c("","NA"))  %>%
              tidyr::gather(any_of(plant_code_spain$code.plant),
                            key="focal.code",value="survival")) %>%
  dplyr::select(year,focal.code,focal,germination,survival) %>%
  rename("code.plant" ="focal.code") %>%
  left_join(plant_code_spain, by="code.plant") %>%
  mutate(survival = case_when(year ==2015 ~ survival/100,
                                 T ~ survival/10)) %>%
  group_by(code.analysis) %>% 
  mutate(s.mean = mean(survival, na.rm = T),
         s.sd = sd(survival, na.rm = T)) %>%
  dplyr::select(code.analysis,s.mean,s.sd) %>%
  unique() %>%
  ungroup()%>%
  as.data.frame()

#---- 3.0 Trait----
#---- 3.1 Export and exploration Trait----
Flower_traits_spain <- read.csv("data/spain_rawdata/spain_flower_trait.csv",
                               header = T, stringsAsFactors = F, sep=",",
                               na.strings = c("","NA")) %>%
  dplyr::select(species,code.plant,code.analysis,flowerwidth) %>% 
  left_join(plant_code_spain %>% dplyr::select(code.plant,code.analysis)) 

Seed_traits_spain <- read.csv( "data/spain_rawdata/spain_seed_mass.txt",
                                header = T, stringsAsFactors = F, sep="\t",
                                na.strings = c("","NA")) %>%
  rename("seed.mass"="mean.seed.mass..100.seeds.aprox.",
         "code.plant"="code")%>%
  dplyr::select(species,code.plant,seed.mass) %>% 
  left_join(plant_code_spain %>% dplyr::select(code.plant,code.analysis)) 

plant_traits_spain <- read.csv("data/spain_rawdata/spain_trait_df.csv",
                               header = T, stringsAsFactors = F, sep=",",
                               na.strings = c("","NA")) %>%
  rename("code.plant"="code") %>%
  left_join(plant_code_spain %>% dplyr::select(code.plant,code.analysis),
            by="code.plant") %>%
  full_join(Flower_traits_spain) %>%
  full_join(Seed_traits_spain ) %>%
  filter( code.analysis %in% final.species.list.spain ) %>%
  dplyr::select(-c(code.plant,species)) %>%
  gather(c("C13","seed.mass", "flowerwidth","C.N", "CS","DR","heigh",
    "TDMr","LAI","SLA","SRL","SRA"),key="trait",value="value") %>%
  aggregate(value ~ trait + code.analysis, mean) %>%
  spread(trait,value) %>%
  column_to_rownames("code.analysis") %>%
  rename("Leaf C to N ratio" ="C.N",
         "C13 water use efficiency"="C13",
         "Canopy shape"="CS",
         "Root diameter"="DR",
         "Stem height"="heigh",
         "Leaf area index"="LAI",
         "Flower width"="flowerwidth",
         "Seed mass"= "seed.mass",
         #"Leaf area"="LeafArea",
         #"Leaf nitrogen cc"="N15",
         "Root mass density"="TDMr") 
  view(plant_traits_spain )
  #---- 3.1.1 Data specific to PAIN----
  # Data from https://uol.de/en/landeco/research/leda/data-files
  #The LEDA Traitbase: A database of life-history traits of Northwest European flora. Journal of Ecology 96: 1266-1274.
  
  LS <- read_delim("data/spain_rawdata/PAIN_traits/LA.csv", delim = ";", escape_double = FALSE, 
                   trim_ws = TRUE, skip = 3) %>%
    rename("species"="SBS name") %>%
    dplyr::filter(species %in% plant_code_spain$species[which(plant_code_spain$code.analysis %in% final.species.list.spain)] |
                    str_detect(species, "^Parapholis")) %>%
    dplyr::select(c("species","single value [mm^2]","mean LS [mm^2]")) %>% as.data.frame() %>%
    rename("value"="single value [mm^2]",
           "mean.value"="mean LS [mm^2]") %>% 
    mutate(mean.value = coalesce(mean.value, value)) %>%
    mutate(trait = "Leaf area")%>%
    group_by( species,trait) %>%
    summarise(value.trait = mean(mean.value)/100)
  
  
  SLA <- read_delim("data/spain_rawdata/PAIN_traits/SLA.csv", delim = ";", escape_double = FALSE, 
                    trim_ws = TRUE, skip = 4) %>%
    rename("species"="SBS name") %>%
    dplyr::filter(species %in% plant_code_spain$species[which(plant_code_spain$code.analysis %in% final.species.list.spain)] |
                    str_detect(species, "^Parapholis")) %>%
    dplyr::select(c("species","single value [mm^2/mg]","mean SLA [mm^2/mg]")) %>% as.data.frame() %>%
    rename("value"="single value [mm^2/mg]",
           "mean.value"="mean SLA [mm^2/mg]") %>% 
    mutate(mean.value = coalesce(mean.value, value)) %>%
    mutate(trait = "SLA") %>%
    group_by( species,trait) %>%
    summarise(value.trait = mean(mean.value)*100)
  
  height <- read_delim("data/spain_rawdata/PAIN_traits/height.csv", delim = ";", escape_double = FALSE, 
                       trim_ws = TRUE, skip = 4) %>%
    rename("species"="SBS name") %>%
    dplyr::filter(species %in% plant_code_spain$species[which(plant_code_spain$code.analysis %in% final.species.list.spain)] |
                    str_detect(species, "^Parapholis")) %>%
    dplyr::select(c("species","single value [m]","mean CH [m]")) %>% as.data.frame() %>%
    rename("value"="single value [m]",
           "mean.value"="mean CH [m]") %>% 
    mutate(mean.value = coalesce(mean.value, value)) %>%
    mutate(trait = "Stem height") %>%
    group_by( species,trait) %>%
    summarise(value.trait = mean(mean.value)*100)
  
  
  LMDC <- read_delim("data/spain_rawdata/PAIN_traits/LMDC.csv", delim = ";", escape_double = FALSE, 
                     trim_ws = TRUE, skip = 4) %>%
    rename("species"="SBS name") %>%
    dplyr::filter(species %in% plant_code_spain$species[which(plant_code_spain$code.analysis %in% final.species.list.spain)] |
                    str_detect(species, "^Parapholis")) %>%
    dplyr::select(c("species","single value [mg/g]","mean LMDC [mg/g]")) %>% as.data.frame() %>%
    rename("value"="single value [mg/g]",
           "mean.value"="mean LMDC [mg/g]") %>% 
    mutate(mean.value = coalesce(mean.value, value)) %>%
    mutate(trait = "LeafDryMatter") %>%
    group_by(species,trait) %>%
    summarise(value.trait = mean(mean.value)/100) 
  
  PAIN_traits <- bind_rows(SLA,height) %>% # Seed mas and leaf dry matter,LS
    dplyr::filter(str_detect(species, "^Parapholis")) %>%
    aggregate(value.trait ~ trait, function(x) mean(x,na.rm=T)) %>%
    spread(trait,value.trait)  %>% 
    mutate(code="PAIN") %>%
    column_to_rownames("code")
  
  plant_traits_spain["PAIN","SLA"] <-PAIN_traits$SLA
  plant_traits_spain["PAIN","Stem height"] <-PAIN_traits$'Stem height'
  view(plant_traits_spain)
#---- 4. Save data SPAIN ----
competition.spain_long <- competition.spain_long[,-c(1:2)] %>%
  rename("focal"="focal.analysis") 

clean.data.spain = list(species_spain = final.species.list.spain,
                      competition_spain =competition.spain_long,
                      abundance_spain.summary=abundance_spain.summary,
                      seed_germination_spain =seed_germination_spain,
                      seed_survival_spain = seed_survival_spain,
                      plant_traits =plant_traits_spain,
                      trait.dist_spain.df = trait.dist_spain.df)
#load("data/clean.data.spain.RData")
clean.data.spain$plant_traits <- plant_traits_spain
## clean.data.spain$trait.dist_spain.df  <- trait.dist_spain.df 
#clean.data.spain$abundance_spain.summary <- abundance_spain.summary
save(clean.data.spain,
     file="data/clean.data.spain.RData")
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
# 3. Data from Perenjori, WA, Australia----
# Main collector Trace Martin, Courtney Taylor, Margie Mayfield
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plant_code_aus <- read.csv("data/aus_rawdata/plant_code_aus.csv",
                           header = T, stringsAsFactors = F, sep=",",
                           na.strings = c("","NA"))

species.list.to.keep.aus <- levels(as.factor(plant_code_aus$final.code))

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.0. Identify most abundance species ----
load("/Users/lisabuche/Documents/Projects/Perenjori/results/community_id_df.csv.gz")
abundance_aus <- community_id_df

abundance_aus.summary.year <- abundance_aus %>%
  mutate(count=as.numeric((count/(scale.area))*625)) %>%
  filter(!stringr::str_detect(id.plot, 'BO_|CA_')) %>%
  dplyr::select(count,year,family) %>%
  aggregate(count~ family, mean)

abundance_aus.summary.year <- abundance_aus %>%
  mutate(count=as.numeric((count/(scale.area))*625)) %>%
  filter(!stringr::str_detect(id.plot, 'BO_|CA_')) %>%
  dplyr::select(count,year,final.code) %>%
  filter( final.code %in% species.list.to.keep.aus) %>%
  aggregate(count~ year + final.code, mean) %>%
  spread(final.code,count) %>%
  mutate(year=as.character(year)) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), function(x) sum(x,na.rm = T)),
                      across(where(is.character), ~"Total")))

names(abundance_aus.summary.year)[abundance_aus.summary.year[10,] > 10]
view(abundance_aus.summary.year)
view(summary_table_aus)
names(summary_table_aus)
#Zac, do you think it is possible that in perenjori, people have been mistaking Goodenia pusilliflora for GOCY? cause GOPU is present before 2020 and not after - but trace collected more than 100 data point and fecundity on it - which is weird bc
# only keep the one that have more than 100 data obs as focal in summary_table_aus
# only keep the species that have maximum of two years without abundance data
final.species.list.aus <- c("ARCA","GOBE","GORO","HYGL",
                            "LARO","MIMY","PEAI","PLDE",
                            "POAR", "POLE","TRCY","WAAC")

# Abundance clean data
abundance_aus.summary <- abundance_aus %>%
  filter(!stringr::str_detect(id.plot, 'BO_|CA_')) %>% # remove reserve outside Perenjory in 2023
  dplyr::select(count,year,final.code,id.plot,collector,scale.width,scale.area ) %>%
  filter( final.code %in% final.species.list.aus) %>%
  rename("species"="final.code") %>%
  rename("individuals" ="count") %>%
  rename("com_id" ="id.plot")%>%
  aggregate(individuals ~ year + species + com_id +collector +scale.width+scale.area , sum) %>%
  mutate(individuals =individuals/scale.area ) %>% # number of individuals per cm^2
  group_by(year) %>%
  mutate(count.zscore = zscore.fct(individuals)) %>%
  ungroup()

#abundance_aus_med <-  abundance_aus.preclean  %>%
# aggregate(individuals ~ species + year, function(x) median(x[!x==0 & !is.na(x)]))

# VISUALISATION
library(pals)
color.palette <- unname(kelly())[1:length(final.species.list.aus)]

abundance_aus_plot <- ggplot() +
  stat_summary(data=abundance_aus.clean ,
               aes(x=as.character(year), y = individuals,
                   group=as.factor(species),
                   color=as.factor(species)),
               fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",size=2) +
  stat_summary(data=abundance_aus.clean ,
               aes(x=as.character(year), y = individuals,
                   group=as.factor(species),
                   color=as.factor(species)),
               fun.y = mean,
               geom = "line",size=1) +
  scale_x_discrete("year",limits=c("2010","2011","2012-2013","2014","2015","2016",
                                   "2017","2018","2019","2020","2021","2022","2023")) +
  scale_y_log10()+
  labs(color="species",y="Mean number of \nindividuals per centimeter",
       title="Density over time of annual plants in Perenjory region") +
  scale_color_manual(values=color.palette) +
  theme_bw() +
  
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=28),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=20, angle=66, hjust=1),
         axis.text.y= element_text(size=20),
         axis.title.x= element_text(size=24),
         axis.title.y= element_text(size=24),
         title=element_text(size=16))

abundance_aus_plot

library(plotly)
plotly::ggplotly(abundance_aus_plot)
ggsave(abundance_aus_plot,
       file="figures/abundance_aus_plot.pdf")

abundance_aus.summary.df <- abundance_aus %>%
  mutate(individuals=as.numeric(count/(scale.width))) %>%
  filter(!stringr::str_detect(id.plot, 'BO_|CA_')) %>% # remove reserve outside Perenjory in 2023
  filter( final.code %in% final.species.list.aus) %>%
  group_by(final.code,species_id,family) %>%
  summarise(count = n(),
            mean = mean(individuals, na.rm = TRUE)*25, 
            sd = sd(individuals, na.rm = TRUE)*25)
head(abundance_aus.summary.df)
write.csv(abundance_aus.summary.df,
          "results/abundance_aus.summary.df.csv")

#---- 3.1. Wainwright 2014 preparation ----

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
  filter(!is.na(seed)) %>%
  mutate(plot=as.character(plot)) %>%
  dplyr::select(any_of(c("plot","focal","year","seed",final.species.list.aus ))) %>%
  mutate(scale=25)


Wainwright2014[is.na(Wainwright2014)] <- 0
str(Wainwright2014) 
head(Wainwright2014) 
table(Wainwright2014$focal)

#---- 3.2. Wainwright 2015 - used in Stouffer2017 Cyclic paper ----
#load("~/Documents/Projects/Original/Stouffer2017_Cyclic/data.Rdata")
# to reattribute other to actual name - intern fils of Claire Wainwright:
Wainwright2015 <-  read.csv("data/aus_rawdata/Wainwright2015_comp_aus.csv",sep=",") %>%
  left_join(read.csv("data/aus_rawdata/Wainwright2015_seed_aus.csv")%>%
              filter(competition =="Highcomp"),
            by=c("treatment","block","plot.ID","plot.letter","focal"),
            relationship = "many-to-many") %>%
  bind_rows(read.csv("data/aus_rawdata/Wainwright2015_seed_aus.csv",sep=",") %>%
            filter(competition =="Solo") %>%
            mutate(density=0,neigh=focal)) %>%
  filter(treatment =="Open") %>%
  mutate(scale="30",
         year=2015) %>%
  rename("plot"="plot.ID")%>%
  left_join(plant_code_aus[,c("genus","species","final.code")] %>%
              mutate(focal =paste0(substr(genus,1L,1L),".",species)) ,
            by=c("focal"),
            relationship = "many-to-many")%>%
  rename("final.code.focal"="final.code") %>%
  left_join(plant_code_aus[,c("genus","species","final.code")] %>%
              mutate(neigh =paste0(substr(genus,1L,1L),".",species)) ,
            by=c("neigh"),
            relationship = "many-to-many") %>%
  rename("final.code.neigh"="final.code") %>%
  mutate(density= ifelse(is.na(density), 0, density)) %>%
  stats::aggregate(density ~ year + plot + competition + final.code.focal + final.code.neigh +scale+ seed,sum) %>%
  rename("focal"="final.code.focal",
         "neigh"="final.code.neigh") %>%
  spread(neigh,density) %>%
  dplyr::select(any_of(c("plot","focal","year",
                  "seed",final.species.list.aus ))) %>%
  mutate(scale=30)

Wainwright2015[is.na(Wainwright2015)] <- 0
str(Wainwright2015) 
head(Wainwright2015) 
table(Wainwright2015$focal)


#---- 3.3. Martyn 2016 preparation----
Martyn2016 <-  read.csv("data/aus_rawdata/Martyn2016_aus.csv" )
Martyn2016 <- Martyn2016 %>%
  dplyr::select(-"focal.ID") %>%
  gather(any_of(plant_code_aus$code), key="code", value="count") %>%
  dplyr::select(plot,seed,focal,code,count) %>%
  left_join(plant_code_aus[,c("code","genus",
                               "species","final.code")], 
             by=c("code"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + focal + final.code + seed,sum) %>%
  spread(final.code,count) %>%
  mutate( year=2016)%>%
  mutate(plot=as.character(plot)) %>%
  dplyr::select(any_of(c("plot","focal","year","seed",final.species.list.aus ))) %>%
  mutate(scale=15)

levels(as.factor(Martyn2016$focal))
#str(Martyn2016)
#table(Martyn2016$focal)
#head(Martyn2016)


#---- 3.4. Pastore 2017 preparation----
Pastore2017_com <-  read.csv("data/aus_rawdata/Pastore2017_com_aus.csv" )
Pastore2017_seed <-  read.csv("data/aus_rawdata/Pastore2017_seed_aus.csv" )

Pastore2017 <- Pastore2017_com %>% 
  left_join(Pastore2017_seed, by=c("plot","focal","id.focal"),
             relationship = "many-to-many") %>%
  filter(!is.na(seed)) %>%
  gather(any_of(plant_code_aus$code), key="code", value="count") %>%
  dplyr::select(plot,seed,focal,code,count,id.focal) %>%
  left_join(plant_code_aus[,c("code","genus",
                               "species","final.code")], 
             by=c("code"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + id.focal + focal + final.code + seed,sum) %>%
  spread(final.code,count) %>%
  mutate( year=2017)%>%
  mutate(plot=as.character(plot)) %>%
  dplyr::select(any_of(c("plot","focal","year","seed",final.species.list.aus )))%>%
  mutate(scale=15)


#view(Pastore2017)
#str(Pastore2017)
#table(Pastore2017$focal)

#---- 3.5. Sevenello 2022 preparation----

Sevenello2022_com <- read.csv("data/aus_rawdata/Sevenello2022_com_aus.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))
Sevenello2022_seed <- read.csv("data/aus_rawdata/Sevenello2022_seed_aus.csv",
                         header = T,stringsAsFactors = F, sep=",",
                         na.strings=c("","NA"))
view(Sevenello2022_seed )
Sevenello2022 <- Sevenello2022_com %>% 
  left_join(Sevenello2022_seed,
              relationship = "many-to-many") %>%
  dplyr::filter(treatment %in% c("OP")) %>%
  filter(!is.na(seed)) %>%
  gather(any_of(plant_code_aus$code), key="code", value="count") %>%
  dplyr::select(plot,seed,focal,code,count) %>%
  left_join(plant_code_aus[,c("code","genus",
                               "species","final.code")], 
             by=c("code"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + focal + final.code + seed,sum) %>%
  spread(final.code,count) %>%
  mutate( year=2022)%>%
  mutate(plot=as.character(plot)) %>%
  dplyr::select(any_of(c("plot","focal","year","seed",final.species.list.aus )))%>%
  mutate(scale=15)

  
#str(Sevenello2022)
#head(Sevenello2022)
#table(Sevenello2022$focal)


#---- 3.6. Taylor 2023 preparation----p
load("data/aus_rawdata/Taylor2023_com_aus.RData")
Taylor2023_com <- bind_rows(neighbourhoods[["BOCL"]] %>% as.data.frame(),
                            neighbourhoods[["BOOP"]] %>% as.data.frame(),
                            neighbourhoods[["CACL"]] %>% as.data.frame(),
                            neighbourhoods[["CAOP"]] %>% as.data.frame(),
                            neighbourhoods[["PECL"]] %>% as.data.frame(),
                            neighbourhoods[["PEOP"]] %>% as.data.frame())
view(Taylor2023_com)
Taylor2023_seed <- read.csv("data/aus_rawdata/Taylor2023_seed_aus.csv",
                               header = T,stringsAsFactors = F, sep=",",
                               na.strings=c("","NA"))
view(Taylor2023_seed ) 

Taylor2023 <- Taylor2023_com %>% 
  rename("focal" ="SPECIES" ) %>%
  left_join(Taylor2023_seed,relationship = "many-to-many") %>%
  rename("plot" = "unique_id") %>%
  filter(!is.na(seed)) %>%
  gather(any_of(c(plant_code_aus$name)[!is.na(c(plant_code_aus$name))]), 
         key="species", value="count") %>%
  dplyr::select(plot,seed,code,count,species) %>%
  rename("focal" ="code",
         "name"="species") %>%
  left_join(plant_code_aus[,c("code","genus","name",
                               "species","final.code")], 
             by=c("name"),
             relationship = "many-to-many") %>%
  aggregate(count ~ plot + focal + final.code + seed,sum) %>%
  spread(final.code,count) %>%
  left_join(plant_code_aus[,c("code","final.code")] %>%
               rename("focal" ="code" ) %>%
               rename("final.focal" ="final.code" ), 
             by=c("focal"),
             relationship = "many-to-many") %>%
  dplyr::select(-"focal") %>%
  rename("focal" ="final.focal" ) %>%
  mutate( year=2023)%>%
  mutate(plot=as.character(plot)) %>%
  dplyr::select(any_of(c("plot","focal","year","seed",final.species.list.aus )))%>%
  mutate(scale=15)

view((Taylor2023))
#str(Taylor2023)
#head(Taylor2023)
#table(Taylor2023$focal)

#---- 3.7. summary table ----

summary_table_aus <- bind_rows(as.data.frame(table(Martyn2016$focal)) %>%
  mutate(year="2016"),
  as.data.frame(table(Wainwright2014$focal)) %>%
    mutate(year="2014"),
  as.data.frame(table(Wainwright2015$focal)) %>%
    mutate(year="2015"),
  as.data.frame(table(Pastore2017$focal)) %>%
    mutate(year="2017"),
  as.data.frame(table(Taylor2023$focal)) %>%
    mutate(year="2023"),
  as.data.frame(table(Sevenello2022$focal)) %>%
    mutate(year="2022")
    ) %>%
  rename("focal"="Var1") %>%
  mutate(focal = case_when(focal=="VERO"~"GORO",
                           focal=="POGN"~"POAR",
                           focal=="VECY"~"GOCY",
                           T~focal)) %>%
  filter(focal %in% final.species.list.aus)  %>%
  spread(focal, Freq) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), function(x) sum(x,na.rm = T)),
                      across(where(is.character), ~"Total")))
view(summary_table_aus)


ggplot(summary_table_aus[which(summary_table_aus$year=="Total"),]%>%
         gather(final.species.list.aus, key="focal",value="Freq"),
       aes(y=Freq,x=focal)) +
  geom_bar(stat="identity",position="dodge") +
  theme_bw()


#---- 3.8. Merge ----
competition_aus <- bind_rows(Martyn2016,Wainwright2014,Wainwright2015,
                             Pastore2017,Taylor2023,
                             Sevenello2022) 


head(competition_aus)
str(competition_aus)
names(competition_aus)


#----3.9. Seed survival and germination----
View(plant_code_aus)
seed_germination_aus <- read.csv(paste0("data/aus_rawdata/Wainwright2015_seedgermss_aus.csv"),
                                   header = T,stringsAsFactors = F, sep=",",
                                   na.strings=c("","NA")) 
  #dplyr::select(code.plant,year,species,germination,survival) %>%
  #left_join(plant_code_aus %>%
    #          dplyr::select(code.plant))

#----4.0. Traits list----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----4.1 Regroup trait dataset----
#Web traits
library(austraits) 
austraits <- load_austraits("10.5281/zenodo.11188867")

aus.names <- plant_code_aus %>%
  filter(code %in% final.species.list.aus) %>%
  select(name) %>%
  unlist() %>%
  as.vector()
trait.of.interest.aus <- c("plant_height","seed_dry_mass","seed_germination",
                           "root_diameter","root_distribution_coefficient")
web.traits <- austraits$traits %>%
  filter(taxon_name  %in% c(aus.names,"Velleia rosea","Podolepis canescence")) %>%
  filter(trait_name %in% trait.of.interest.aus) %>%
  mutate(value = as.numeric(value)) %>% 
  aggregate(value ~ taxon_name + trait_name + unit,mean) %>%
  dplyr::select(taxon_name,trait_name,value) %>%
  spread(trait_name,value) %>%
  mutate(taxon_name=case_when(taxon_name == "Velleia rosea"~"Goodenia rosea",
                              T~taxon_name)) %>%
  rename("species"="taxon_name")

view(web.traits)
# internal dataset
OG.plant_traits_aus <- read.csv("data/aus_rawdata/Traits_database.csv",
                             header = T, stringsAsFactors = F, sep=",",
                             na.strings = c("","NA"))
plant_germ_aus <- read.csv("~/Documents/Projects/Perenjori/data/seed_rates.csv",
                           header = T, stringsAsFactors = F, sep=",",
                           na.strings = c("","NA"))
ManuelSevenello.traits <- read.csv("data/aus_rawdata/ManuelSevenello_traits.csv")
WinnieSiu.traits <- read.csv("data/aus_rawdata/WinnieSiu_traits.csv") %>%
  gather(root.biomass,srl,key="traits",value="value") %>%
  aggregate(value ~ traits + species + final.code, median) %>%
  spread(traits,value)


aus_traits_df <- OG.plant_traits_aus %>%
  left_join(plant_code_aus[,c("name","code","final.code")] %>%
              rename("species" ="name" ), 
            by=c("species"),
            relationship = "many-to-many") %>%
  filter(final.code %in% final.species.list.aus ) %>%
  dplyr::select(final.code,species,aboverground.biomass.mg, 
                height.mm,mean.seed.mass.mg,sla.mm2.mg,
                width.longest.mm, width.90.from.longest.mm,
                C13.plant.permill,
                total.root.length.cm,number.of.root.tips,
                root.volume.less.than.0.5mm.diameter.mm3,
                total.root.volume.cm3) %>%
  group_by(final.code,species) %>%
  summarise_all(list(~ median(.x, na.rm = TRUE))) %>%
  left_join(plant_germ_aus[,c("name","code","germ","surv")] %>%
              rename("species" ="name",
                     "final.code"="code")) %>%
  left_join(web.traits) %>%
  left_join(ManuelSevenello.traits) %>%
  left_join(WinnieSiu.traits) %>%
  left_join(clean.data.aus$abundance_aus.summary %>%
              dplyr::select(species,individuals) %>%
              group_by(species) %>%
              summarise(Abundance = mean(individuals, na.rm=T)) %>% 
              rename("final.code"="species")) %>%
  left_join(clean.data.aus$competition_aus %>%
              dplyr::select(seed,focal) %>%
              group_by(focal) %>%
              summarise(Fecundity = mean(seed, na.rm=T)) %>% 
              rename("final.code"="focal")) 

view(aus_traits_df)

write.csv(aus_traits_df,
          "data/aus_traits_df.csv")
aus_traits_df <- read.csv("data/aus_traits_df.csv")
#rename and clean data set for analysis
plant_traits_aus <- aus_traits_df %>%
  mutate(mean.seed.mass.mg = case_when(is.na(mean.seed.mass.mg)~seed.dry.mass,
                                       T~mean.seed.mass.mg)) %>%
  mutate(CanopyArea = pi*width.longest.mm*width.90.from.longest.mm,
         sla.mm2.mg = 10*sla.mm2.mg, # to have cm2/g
         srl= 1000*srl, # to have cm2/g
         height.mm = height.mm/10, # to have cm
         Root.mass.density =root.biomass/(total.root.volume.cm3*1000)) %>% # convert to mg/cm3
  dplyr::select(final.code,
                flower.size.numb,
                Fecundity,
                #surv,
                mean.seed.mass.mg,
                height.mm,
                sla.mm2.mg,
                #width.longest.mm, 
                #width.90.from.longest.mm,
                CanopyArea,
                total.root.length.cm,
                number.of.root.tips,
                C13.plant.permill,
                srl,
                Root.mass.density 
                #root.volume.less.than.0.5mm.diameter.mm3
                )  %>%
  rename(#"Root volume"="root.volume.less.than.0.5mm.diameter.mm3",
         "SLA"="sla.mm2.mg",
         "SRL"="srl",
         "Root length"="total.root.length.cm",
         "Canopy shape" = "CanopyArea",
         "C13 water use efficiency"="C13.plant.permill",
         #"Canopy width"="width.longest.mm",
         #"Canopy width 90deg"="width.90.from.longest.mm",
         "Stem height"="height.mm",
         "Seed mass"= "mean.seed.mass.mg",
         "Root tips"="number.of.root.tips",
         "Root mass density"="Root.mass.density",
         "Flower width"="flower.size.numb" ,
         "Mean fecundity" ="Fecundity") %>%
  #dplyr::filter(!final.code %in% c("ARCA","PEAI")) %>%
  column_to_rownames("final.code") 
view(plant_traits_aus)

# ---- 5. Save data AUS ----
clean.data.aus = list(seed_germination_aus=seed_germination_aus,
                      species_aus = final.species.list.aus,
                      competition_aus =competition_aus,
                      abundance_aus.summary=abundance_aus.clean,
                      #abundance_aus.preclean = abundance_aus.preclean,
                      aus_above_grouping = aus_above_grouping,
                      aus_pol_grouping =aus_pol_grouping,
                      aus_below_grouping = aus_below_grouping,
                      plant_pol.traits_aus.dist = plant_pol.traits_aus.dist,
                      trait.dist_aus.df  =trait.dist_aus.df,
                      plant_traits = plant_traits_aus)
load("data/clean.data.aus.RData")
#clean.data.aus$abundance_aus.summary <- abundance_aus.summary
#clean.data.aus$aus_above_grouping <-aus_above_grouping
#clean.data.aus$aus_pol_grouping <-aus_pol_grouping
#clean.data.aus$aus_below_grouping <-aus_below_grouping
#clean.data.aus$trait.dist_aus.df <-trait.dist_aus.df
clean.data.aus$plant_traits <-plant_traits_aus
save(clean.data.aus,
     file="data/clean.data.aus.RData")
