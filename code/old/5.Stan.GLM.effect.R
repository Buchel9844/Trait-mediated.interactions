
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
library(brms)
library(loo)
library(rstanarm)
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
library(wesanderson) # for color palette
library(ggthemes) 
library(grid)

home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
year.int = "All"    
Code.focal = "LEMA"
source(paste0(home.dic,"code/PopProjection_toolbox.R"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Visualisation Species interactions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load(file= paste0(home.dic,"results/Parameters_",Code.focal,"_All.Rdata"))
load(file= paste0(home.dic,"results/Inclusion",Code.focal,"_",year.int,".Rdata"))
load(paste0(home.dic,"results/test.sigmoid.rData"))

year.veclevels <- factor(c(Inclusion$DataVec$year))
levels(year.veclevels) <- list("2015"="1","2016"="2","2017"="3",
                               "2018"="4","2019"="5","2020"="6","2021"="7")
year.veclevels <- as.vector(year.veclevels)

SpMatrix <- Inclusion$DataVec$SpMatrix

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 4. Model of interaction according to PDSI and neighbourhood ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SpMatrix.year <- as.data.frame(SpMatrix) %>%
  mutate(year = as.numeric(year.veclevels)) 

interaction.effect <- NULL
for( n in levels(as.factor(test.sigmoid$neigh))){
  for( y in levels(as.factor(SpMatrix.year$year))){
    print(paste(n," ",y))
    test.sigmoid.n <- test.sigmoid %>% 
      dplyr::filter( neigh == n  & year == y & density==0) %>%
      select(!"density") %>%
      unique()
    
    SpMatrix.year.n <- SpMatrix.year %>% 
      dplyr::select(all_of(n),year) %>%
      dplyr::filter( year == as.numeric(y)) %>%
      unique()
    
    interaction.effect.out <- alpha_function4(rep(test.sigmoid.n$alpha_initial,each=nrow(SpMatrix.year.n)),
                                              rep(test.sigmoid.n$alpha_slope, each=nrow(SpMatrix.year.n)),
                                              rep(test.sigmoid.n$alpha_c, each=nrow(SpMatrix.year.n)),
                                              rep(as.numeric(SpMatrix.year.n[,n]),
                                                  times=nrow(test.sigmoid.n)),
                                              rep(test.sigmoid.n$N_opt_mean, each=nrow(SpMatrix.year.n)))
    
    interaction.effect.n.y <-  data.frame(effect.raw =  interaction.effect.out,
                                          density = rep(as.numeric(SpMatrix.year.n[,n]),
                                                        times=nrow(test.sigmoid.n)))  %>%
      mutate(effect.on.lambda = effect.raw * density,
             neigh=n,
             year=y)
    
    interaction.effect <-  bind_rows(interaction.effect , interaction.effect.n.y )
  }
}

spain_env_pdsi  <- read.csv(paste0(home.dic,
                                   "results/spain_env_pdsi.csv"))

spain_env_pdsi_med <- spain_env_pdsi %>%
  filter(month >3 & month < 8) %>%
  aggregate(spain_pdsi ~ year, median) 

interaction.effect.env <- interaction.effect %>%
  mutate(year=as.numeric(year)) %>%
  left_join(spain_env_pdsi_med, by="year") %>%
  mutate(density=scale(density)) %>%
  unique()

interaction.effect.env %>%
  ggplot(aes(x=effect.on.lambda)) +
  geom_density()+
  facet_wrap(.~as.factor(neigh),scale="free")

interaction_stan.df <- NULL
for( n in levels(as.factor(interaction.effect.env$neigh))){
  print(n)
  interaction_stan_0 <- stan_glm(effect.on.lambda ~  1 + density,
                                 data = interaction.effect.env%>%dplyr::filter(neigh==n),
                                 family = gaussian(),
                                 prior_intercept = normal(0, 1), chains = 3, 
                                 iter = 1000, warmup = 500)
  interaction_stan_psi <- stan_glm(effect.on.lambda ~  1 + density + spain_pdsi,
                                   data = interaction.effect.env%>%dplyr::filter(neigh==n),
                                   family = gaussian(),
                                   prior_intercept = normal(0, 1), chains = 3, 
                                   iter = 1000, warmup = 500)
  interaction_stan_int <- stan_glm(effect.on.lambda ~  1 + density*spain_pdsi,
                                   data = interaction.effect.env%>%dplyr::filter(neigh==n),
                                   family = gaussian(),
                                   prior_intercept = normal(0, 1), chains = 3, 
                                   iter = 1000, warmup = 500)
  
  loo.0 <- loo(interaction_stan_0, cores = 2)
  loo.psi <- loo(interaction_stan_psi, cores = 2)
  loo.int <- loo(interaction_stan_int, cores = 2)
  loo.comp <- loo_compare(loo.0 , loo.psi, loo.int) %>%
    as.data.frame() %>%
    rownames_to_column(var = "model")
  
  
  interaction_stan.df.n <- bind_rows(as.data.frame(interaction_stan_0$stan_summary) %>%
                                       mutate(model="interaction_stan_0"),
                                     as.data.frame(interaction_stan_psi$stan_summary)%>%
                                       mutate(model="interaction_stan_psi")) %>%
    bind_rows(as.data.frame(interaction_stan_int$stan_summary)%>%
                mutate(model="interaction_stan_int")) %>%
    rownames_to_column(var = "parameter") %>% 
    left_join(loo.comp, by="model") %>%
    mutate(neigh=n)
  
  interaction_stan.df <- bind_rows(interaction_stan.df,
                                   interaction_stan.df.n)
  
}

write.csv(interaction_stan.df,
          paste0(home.dic,"results/interaction_stan.df.csv"))

interaction_stan.df<- read.csv(paste0(home.dic,"results/interaction_stan.df.csv"))
interaction_stan.df <- interaction_stan.df[,-1]
interaction_stan.df <- interaction_stan.df[,-1]

names(interaction_stan.df)
view(interaction_stan.df)
interaction_stan.df <- interaction_stan.df %>%
  as.data.frame()%>%
  setNames(c("parameter","mean","se_mean","sd","Q.2.5",
             "Q.10","Q.25","Q.50","Q.75","Q.90","Q.97.5",
             "n_eff","Rhat","model","elpd_diff",
             "se_diff","elpd_loo","se_elpd_loo","p_loo","se_p_loo",
             "looic","se_looic","neigh")) %>%
  mutate(parameter = rep(c("intercept","density","sigma","mean_ppd","log_post",
                           "intercept","density","spain_pdsi","sigma","mean_ppd","log_post",
                           "intercept","density","spain_pdsi","density:spain_pdsi","sigma","mean_ppd","log_post"),
                         times=nlevels(as.factor(interaction_stan.df$neigh))))


plot.glm.spdi <- interaction_stan.df %>%
  filter(parameter %in% c("intercept","density:spain_pdsi","spain_pdsi","density")) %>%
  mutate(parameter = case_when(parameter=="density:spain_pdsi" ~ "density:pdsi",
                               parameter=="spain_pdsi" ~ "pdsi",
                               T~parameter)) %>%
  mutate(parameter = factor(parameter,
                            levels=c("intercept","density","pdsi","density:pdsi"))) %>%
  mutate(model = case_when(model =="interaction_stan_0" ~ "effect ~ density",
                           model =="interaction_stan_psi" ~ "effect ~ density + pdsi",
                           model =="interaction_stan_int" ~ "effect ~ density + pdsi + density:pdsi")) %>%
  ggplot(aes(y=mean,x=as.factor(parameter),
             color=model)) +
  geom_errorbar(aes(y=mean,
                    ymin=Q.2.5,
                    ymax=Q.97.5)) +
  geom_hline(yintercept=0) +
  labs(x="impact of parameters",
       y="Effect on the intrinsic growth rate of LEMA",
       color="Stan glm fitting,\n from simple to complexe:",
       subtitle = "Mean effect, with 95% interval",
       title = " Distingeling the effect of density and environement (pdsi) on a focal species: LEMA")+
  facet_wrap(.~as.factor(neigh),
             scale="free",nrow=2) +
  theme_bw() +
  theme(legend.key.size = unit(1.2, 'cm'),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=10),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12),
        #axis.ticks.x=element_blank(),
        axis.text.x= element_text(size=10, 
                                  angle=66,
                                  vjust=0.95,
                                  hjust=1),
        axis.text.y= element_text(size=10),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        title=element_text(size=12))


plot.glm.spdi


ggsave(plot.glm.spdi,
       file=paste0(home.dic,"figures/plot.glm.spdi.png"))

library(tidybayes)


