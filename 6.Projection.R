
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 1.1. Import packages ----
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
library(ggridges)

home.dic <- "/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
year.int = "All"    
Code.focal = "LEMA"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Projection species abundance ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spain_env_pdsi<- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv"))

spain_env_pdsi_med <- spain_env_pdsi %>%
  dplyr::filter(month >3 & month < 8) %>%
  aggregate(spain_pdsi ~ year, median) 

year.levels <- names(parameter$df_lambda_sd)

#---- 2.1 Lambda ----
list.lambda <- list()
for(Code.focal in species.spain){
  year.levels <- RawData[[Code.focal]]$year.levels
  list.lambda[[Code.focal]] <- Parameters[[Code.focal]]$df_lambda_sd %>%
    mutate_at(year.levels, ~ rowSums(cbind(., Parameters[[Code.focal]]$df_lambda_mean))) %>%
    gather(key="year", value="lambda_sd") %>%
    mutate(year=as.numeric(year)) %>%
    left_join(spain_env_pdsi_med, by="year") %>%
    mutate(env.cat = case_when(spain_pdsi > 2 ~"wet",
                               spain_pdsi < - 2 ~"dry",
                               T~ "variable"))
  
  
}
# visualisation of the dry/wet/variable intrinsici growth rate
Code.focal <- "LEMA"
year.levels <- RawData[[Code.focal]]$year.levels
Parameters[[Code.focal]]$df_lambda_sd %>%
  mutate_at(year.levels, ~ rowSums(cbind(., Parameters[[Code.focal]]$df_lambda_mean))) %>%
  gather(key="year", value="lambda_sd") %>%
  mutate(year=as.numeric(year)) %>%
  left_join(spain_env_pdsi_med, by="year") %>%
  mutate(env.cat = case_when(spain_pdsi > 2 ~"wet",
                             spain_pdsi < - 2 ~"dry",
                             T~ "variable")) %>%
  #filter(env.cat =="wet") %>%
  ggplot(aes(x=lambda_sd,y=env.cat,fill=spain_pdsi)) +
  geom_density_ridges() +
  geom_vline(xintercept=median(Parameters[[Code.focal]]$df_lambda_mean[,1]))  
  
head( list.lambda[[Code.focal]] )

#---- 2.2 Alpha----

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
list.alpha <- list()

for(Code.focal in species.spain){
  df_alpha_generic_param = Parameters[[Code.focal]]$df_alpha_generic_param
  
  Sp.names = colnames(RawData[[Code.focal]]$SpMatrix)
  year.levels <- RawData[[Code.focal]]$year.levels
  print(Code.focal)
  param.neigh <- NULL
  for( neigh in Sp.names){
    for(year.int in year.levels){
      param.neigh.n <- NULL
      print(year.int)
      alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                             neigh]
      if(RawData[[Code.focal]]$Inclusion_alpha_initial[year.int,neigh]>0){
        alpha_initial = alpha_initial + Parameters[[Code.focal]]$df_alpha_initial[,paste(year.int,neigh,sep="_")]
      }
      
      alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                           neigh]
      if(RawData[[Code.focal]]$Inclusion_alpha_slope[year.int,neigh]>0){
        alpha_slope = alpha_slope + Parameters[[Code.focal]]$df_alpha_slope[,paste(year.int,neigh,sep="_")]
      }
      
      alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                       neigh]
      if(RawData[[Code.focal]]$Inclusion_c[year.int,neigh]>0){
        alpha_c = alpha_c + Parameters[[Code.focal]]$df_alpha_c[,paste(year.int,neigh,sep="_")]
      }
      
      
      param.neigh.n <- data.frame(neigh = neigh, 
                                year = year.int,
                                alpha_initial = alpha_initial,
                                alpha_slope = alpha_slope,
                                alpha_c=  alpha_c,
                                N_opt_mean = Parameters[[Code.focal]]$df_N_opt_mean[,neigh],
                                focal=Code.focal)
      param.neigh <- bind_rows(param.neigh,param.neigh.n )
    }
  }
  list.alpha[[Code.focal]] <-  param.neigh
}
#---- 2.3 Germination ----
seed_germination_spain <- read.csv(paste0("data/spain_rawdata/seed_germination.csv"),
                            header = T,stringsAsFactors = F, sep=",",
                            na.strings=c("","NA"))

seed_germination_spain <- seed_germination_spain %>%
  mutate(code.plant=toupper(focal.code)) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  mutate(germination = case_when(year ==2015 ~ germination/100,
                                 T ~ germination/50)) %>%
  #aggregate(germination ~ focal.code + year,sd)
  group_by(code.analysis) %>% 
  mutate(g.mean = mean(germination, na.rm = T),
         g.sd = sd(germination, na.rm = T)) %>%
  select(code.analysis,g.mean,g.sd) %>%
  unique() %>%
  ungroup() %>%
  as.data.frame()

#view(seed_germination_spain)

seed_survival_spain <- read.csv(paste0("data/spain_rawdata/seed_survival_2020.csv"),
                                   header = T,stringsAsFactors = F, sep=",",
                                   na.strings=c("","NA"))

seed_survival_spain <- seed_survival_spain%>%
  gather(any_of(plant_code_spain$code.plant),
         key="code.plant",value="survival") %>%
  left_join(plant_code_spain, by="code.plant") %>%
  mutate(survival = survival/10) %>%
  group_by(code.analysis) %>% 
  mutate(s.mean = mean(survival, na.rm = T),
         s.sd = sd(survival, na.rm = T)) %>%
  select(code.analysis,s.mean,s.sd) %>%
  unique() %>%
  ungroup()%>%
  as.data.frame()

#---- 2.4 Projection ----
abundance_spain <- read.csv(paste0("data/abundance_spain.csv"),
                            header = T,stringsAsFactors = F, sep=",",
                            na.strings=c("","NA"))
plant_code_spain <- read.csv(paste0( "data/plant_code_spain.csv"),
                             header = T, stringsAsFactors = F, sep=",",
                             na.strings = c("","NA"))

abundance_spain_short <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  filter(year==2021) %>%
  aggregate(individuals ~code.analysis  , mean) %>% 
  spread(code.analysis,individuals) 

abundance_spain_rare <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  filter(code.analysis=="rare") %>% 
  aggregate(individuals ~ subplot + year + code.plant, sum) %>% 
  aggregate(individuals ~ subplot + year , mean) %>% 
  ggplot(aes(y=individuals ,x=as.factor(year))) +
  geom_boxplot()
  
abundance_spain_rare <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  filter(code.analysis=="ME.sp") %>% 
  aggregate(individuals ~ subplot + year + code.plant, sum) %>% 
  ggplot(aes(y=individuals ,x=as.factor(year))) +
  geom_boxplot()

#---- 3. Projection ----
time.step <- 10 
dry.year <- c("2019","2020","2021")
var.year <- c("2015","2016")
wet.year <- c("2018","2017")


g <-  data.frame(matrix(0.7,ncol=length(species.spain),
                         nrow=1)) 
names(g) <- species.spain
s <-  data.frame(matrix(0.7,ncol=length(species.spain),
                        nrow=1))
names(s) <- species.spain

list.projection <- list()
for(env.cat.int in c("dry","wet","variable")){
  print(env.cat.int)
df.projection <- NULL
 for(iteration in 1:100){
   print(iteration)
   g.df <-seed_germination_spain %>% 
     rowwise() %>%
     mutate(g.est = abs(rnorm(1, mean = g.mean,sd=g.sd))) %>% 
     select(code.analysis, g.est) %>%
     spread(code.analysis, g.est)
   
   s.df <- seed_survival_spain %>% 
     rowwise() %>%
     mutate(s.est = abs(rnorm(1, mean = s.mean,sd=s.sd))) %>%
     select(code.analysis, s.est) %>%
     spread(code.analysis, s.est)
   
  df.projection.n <- ceiling(abundance_spain_short/g.df) # extracting seed bank from abundance and germination 
  df.projection.n[2:(time.step+1),] <- NA
  df.projection.n$time <- 1:(time.step+1) # first abundance is the abundance in 2021
  


  for(time in 1:time.step){
      density.vec <- df.projection.n[time,]
      # determine the density of each focal in the next time step
      
    for(Code.focal in species.spain){ 
      # sample the germination
      g.focal <-  g.df[,Code.focal]
      s.focal <-  s.df[,Code.focal]
       
      # sample the intrinsic growth rate
    lambda.vec <- list.lambda[[Code.focal]] %>%
      filter(.[,"env.cat"] == env.cat.int) 
    number.vec <- sample(1:nrow(lambda.vec),1)
    lambda.val = lambda.vec[number.vec,"lambda_sd"]
    year.int = lambda.vec[number.vec,"year"]
   # determine the interaction effect
    alpha.vec <- list.alpha[[Code.focal]] %>%
    filter(.[,"year"] == year.int)
    
    Sp.names = levels(as.factor(alpha.vec$neigh))
    alpha.impact <- NULL
    # effect of neighbours on focal
      for(neigh.int in species.spain){
      density.neigh <- df.projection.n[time,neigh.int]
    
      if(neigh.int %in% Sp.names){
      alpha.vec.n <- alpha.vec %>%
        filter(neigh == neigh.int)
      }else{
      alpha.vec.n <- alpha.vec %>%
        filter(neigh == "others")
      }

      alpha.impact.n <- data.frame(Code.focal=Code.focal,
                              neigh =neigh.int,
                              density.neigh = density.neigh,
                              g.neigh = abs(rnorm(1,mean=seed_germination_spain[which(seed_germination_spain$code.analysis==neigh.int),"g.mean"],
                                              sd=seed_germination_spain[which(seed_germination_spain$code.analysis==neigh.int),"g.sd"])),
                              effect = mean(alpha_function4(alpha.vec.n$alpha_initial,
                                          alpha.vec.n$alpha_slope,
                                          alpha.vec.n$alpha_c,
                                          density.neigh,
                                          alpha.vec.n$N_opt_mean))) %>%
      mutate(total.impact = density.neigh*g.neigh*effect)
    
      alpha.impact <- bind_rows( alpha.impact, alpha.impact.n)
      
        }
  if(df.projection.n[time,Code.focal] < 1|
     is.infinite(df.projection.n[time,Code.focal])|
     is.na(df.projection.n[time,Code.focal])){
    df.projection.n[time,Code.focal] = 0
  }
    fecundity.n <- exp(lambda.val + sum(alpha.impact$total.impact))
    abundance.n <- ((1-g.focal)*s.focal + fecundity.n*g.focal)*df.projection.n[time,Code.focal]
    df.projection.n[time+1,Code.focal] <- abundance.n 
      }
    }
  df.projection.n$iteration <- iteration
  df.projection <- bind_rows(df.projection,df.projection.n)
  }
  list.projection[[env.cat.int]] <-  df.projection
}
 
nsafe_colorblind_palette <- c("black","#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)

plot.projection <- list()
for(env.cat.int in c("dry","wet","variable")){
plot.projection[[env.cat.int]] <- list.projection[[env.cat.int]] %>%
  gather(all_of(species.spain), key="species",value="density") %>%
  filter(density < 10000) %>%
  ggplot(aes(x=time,y=density,
             group=as.factor(species),
             color=as.factor(species))) + 
  stat_summary(fun = mean,
               fun.min = function(x) quantile(x,0.025), 
               fun.max = function(x) quantile(x,0.975), 
               geom = "pointrange",size=2) +
  stat_summary(fun = mean,
               geom = "line",size=1) +
  labs(color="Groups",
       y="Averaged number of individuals \n in 1meter squarred plot",
       x="year",
       title="Density over time of annual plants in Caracoles") +
  #coord_cartesian( xlim = NULL, ylim = c(0,500),expand = TRUE, default = FALSE, clip = "on") +
  scale_color_manual(values=safe_colorblind_palette) +
  #scale_y_log10() +
  #scale_y_continuous(limits=c(0,10000)) +
  coord_cartesian(ylim=c(0,1000)) +
  theme_bw() +
  guides(fill=guide_legend(nrow = 1,
                           direction="horizontal",
                           byrow = TRUE,
                           title.hjust = 0.1)) +
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
}
ggarrange(plotlist = plot.projection,
          nrow=3,
          common.legend = T,
          legend="bottom")
plot.projection[["wet"]]
plot.projection[["dry"]]
plot.projection[["variable"]]

