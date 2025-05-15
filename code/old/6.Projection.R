
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
#---- 2.0 Load data ----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))

spain_env_pdsi <- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv")) %>%
  group_by(year) %>% 
  mutate(pdsi.mean = mean(spain_pdsi, na.rm = T),
         pdsi.sd = sd(spain_pdsi, na.rm = T)) %>%
  select(year,pdsi.mean,pdsi.sd) %>%
  unique() %>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(year = as.integer(year))

head(spain_env_pdsi)
#---- 2.1 Lambda ----
source(paste0(home.dic,"code/PopProjection_toolbox.R"))

list_lambda <- list()

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  for(Code.focal in Code.focal.list){ #focal.levels
    
    load(file= paste0(home.dic,
                      "results/Parameters_",Code.focal,"_",country,".Rdata"))
    
    df_alpha_generic_param = parameter$df_alpha_generic_param
    
    list_alpha_mean[[paste0(country,"_",Code.focal)]] <- df_alpha_generic_param %>%
      group_by(parameter) %>% 
      summarise(across(any_of( Code.focal.list),mean)) %>%
      bind_rows(as.data.frame(colMeans(parameter$df_N_opt))%>%
                  rownames_to_column("neigh") %>%
                  spread(neigh,
                         "colMeans(parameter$df_N_opt)") %>%
                  mutate(parameter="N_opt")) %>%
      column_to_rownames("parameter")
    
    list_alpha_sd[[paste0(country,"_",Code.focal)]] <- df_alpha_generic_param %>%
      group_by(parameter) %>% 
      summarise(across(any_of( Code.focal.list),sd)) %>%
      bind_rows(as.data.frame(sapply(parameter$df_N_opt,sd))%>%
                  rownames_to_column("neigh") %>%
                  spread(neigh,
                         "sapply(parameter$df_N_opt, sd)") %>%
                  mutate(parameter="N_opt"))%>%
      column_to_rownames("parameter")
    
  }
}
#---- 2.2 Alpha----


source(paste0(home.dic,"code/PopProjection_toolbox.R"))

list_alpha_mean <- list()
list_alpha_sd <- list()

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  for(Code.focal in Code.focal.list){ #focal.levels
    
    load(file= paste0(home.dic,
                      "results/Parameters_",Code.focal,"_",country,".Rdata"))
    
    df_alpha_generic_param = parameter$df_alpha_generic_param
    
    list_alpha_mean[[paste0(country,"_",Code.focal)]] <- df_alpha_generic_param %>%
      group_by(parameter) %>% 
      summarise(across(any_of( Code.focal.list),mean)) %>%
      bind_rows(as.data.frame(colMeans(parameter$df_N_opt))%>%
                  rownames_to_column("neigh") %>%
                  spread(neigh,
                         "colMeans(parameter$df_N_opt)") %>%
                  mutate(parameter="N_opt")) %>%
      column_to_rownames("parameter")
    
    list_alpha_sd[[paste0(country,"_",Code.focal)]] <- df_alpha_generic_param %>%
      group_by(parameter) %>% 
      summarise(across(any_of( Code.focal.list),sd)) %>%
      bind_rows(as.data.frame(sapply(parameter$df_N_opt,sd))%>%
                  rownames_to_column("neigh") %>%
                  spread(neigh,
                         "sapply(parameter$df_N_opt, sd)") %>%
                  mutate(parameter="N_opt"))%>%
      column_to_rownames("parameter")
    
  }
}



#---- 3. Projection ----
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_short_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  sgerm.df <- get(paste0("clean.data.",country))[[paste0("seed_germination_",country)]]
  ssurv.df <- get(paste0("clean.data.",country))[[paste0("seed_survival_",country)]]
  env_pdsi <- read.csv(paste0(home.dic,"results/",country,"_env_pdsi.csv")) %>%
    group_by(year) %>% 
    mutate(PDSI.mean = mean(get(paste0(country,"_pdsi")), na.rm = T),
           PDSI.sd = sd(get(paste0(country,"_pdsi")),na.rm = T)) %>%
    select(year,PDSI.mean,PDSI.sd) %>%
    unique() %>%
    ungroup() %>%
    as.data.frame()
  for(env.cat.int in c("dry","wet","variable")){
  for(Code.focal in Code.focal.list){ #focal.levels
    load(file =paste0(home.dic,"results/Parameters_",Code.focal,"_",country,".Rdata"))
    #assign(paste0("parameter_",Code.focal),
    #       parameter)
    Parameters[[paste(country,"_",Code.focal)]] <- parameter
  }
}


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
 
safe_colorblind_palette <- c("black","#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
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
       #y="Averaged number of individuals \n in 1meter squarred plot",
       x="year",
       title=env.cat.int) +
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

