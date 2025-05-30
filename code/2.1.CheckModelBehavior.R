#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(cli)
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
#setwd("/home/lbuche/Eco_Bayesian/chapt3")
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2. Check Model----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
country.list <- c("aus","spain") #aus
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))

bayes_RMDS <- function(post.draws,
                       var_name,
                       obs.var){
  y <- post.draws[[var_name]]
  x <-  obs.var
  n <- length(x)
  rmds <- c()
  for(i in 1:nrow(y)){
    rmds[i] <- sqrt(sum((x - y[i,])^2)/n)
  }
  data.frame(mean.fec=mean(x),
             sd.fec= sd(x),
             mean.rmds = mean(rmds),
             min.rmds = min(rmds),
             max.rmds =max(rmds),
             Sim.Bias = (mean(x)-mean(y))^2)
}
ModelCheck.df <- NULL
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  for(Code.focal in Code.focal.list){ #focal.levels
    print(paste(country,Code.focal))
    if(Code.focal =="BEMA") next
    load(paste0(project.dic,"results/Modelfit_",Code.focal,"_",country,".rds"))
    
    
    ModelfitPosteriors <- rstan::extract(Modelfit)
    
    
    load(file= paste0(project.dic,"results/ModelfitPosteriors",Code.focal,"_",country,".Rdata"))
    
    load(paste0(project.dic,"results/Parameters_",Code.focal,"_",country,".Rdata"))
    options(mc.cores =  parallel::detectCores())
    loo.fit <- rstan::loo(Modelfit,pars ="F_sim")
    loo.fit.df <- as.data.frame(loo.fit$estimates)
    DataVec <- parameter[["DataVec"]]
    
    mc <- data.frame(focal =Code.focal , 
                     country=country,
                     n.obs =  DataVec$N,
                     n.parameter  = (1+1 + dim(ModelfitPosteriors$c)[2]*4), # lambda_mean, lambda_sd, number of neighbours X 4 ( number of paramters to estimate for species interaactions)
                     Rhat = max(summary(Modelfit)$summary[,"Rhat"],na.rm =T),
                     Neff = min(summary(Modelfit)$summary[,"n_eff"],na.rm = T)) %>%
      bind_cols( bayes_RMDS(ModelfitPosteriors,
                            var_name = 'F_hat',
                            DataVec$Fecundity)) %>%
      mutate(perc.K = sum(loo.fit$diagnostics$pareto_k > min(1-(1/log10(DataVec$N)),0.7))/
               length(loo.fit$diagnostics$pareto_k),
             p_loo = paste0(round(loo.fit.df["p_loo","Estimate"],digits = 1),"+-",
                            round(loo.fit.df["p_loo","SE"],digits = 1))) %>%
      mutate(p_loo.check = case_when(loo.fit.df["p_loo","Estimate"]-loo.fit.df["p_loo","SE"] < n.parameter ~ "Good specification",
                                     T~"misspecification")) %>%
      mutate(elpd_loo = loo.fit.df["elpd_loo","Estimate"],
             elpd_loo_se = loo.fit.df["elpd_loo","SE"])
    
    ModelCheck.df <- bind_rows(ModelCheck.df ,  mc)
    write.csv(ModelCheck.df,
              file=paste0(home.dic,"results/ModelCheck.df.csv"))
  }
}
