#Inversed modelling to predict lambda 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages("HDInterval")
library("HDInterval")
#install.packages("tidyverse")
library(cli)
library(tidyverse)
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
setwd("/home/lbuche/Eco_Bayesian/chapt3")
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"



#home.dic <- "/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
#year.int = "All"    
#Code.focal = "LEMA"
country.list <- c("aus","spain")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Projection species abundance ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))


#---- 2.1 Alpha----


source(paste0(home.dic,"code/PopProjection_toolbox.R"))

list_alpha_mean <- list()
list_alpha_sd <- list()
country.list="spain"
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


#---- 3. Extraction of realised intrinsic growth rate ----
country.list="spain"

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_short_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  sgerm.df <- get(paste0("clean.data.",country))[[paste0("seed_germination_",country)]]
  ssurv.df <- get(paste0("clean.data.",country))[[paste0("seed_survival_",country)]]
  for(Code.focal in Code.focal.list){ #focal.levels
    year.levels =  levels(as.factor(abundance_short_df$year))
    year.int = nlevels(as.factor(abundance_short_df$year))
    com_id =  levels(as.factor(abundance_short_df$com_id))
    
    growth.ratio.df <- NULL
    for(y in 1: year.int){
      for(i in  com_id){
        if(y ==7) next
        growth.ratio.df <- bind_rows(growth.ratio.df,
                                     data.frame(country=country,
                                                focal = Code.focal,
                                                year = as.integer(year.levels[y]),
                                                com_id=i,
                                                growth.ratio=  abundance_short_df[which( abundance_short_df$year== as.integer(year.levels[y+1]) &
                                                                                           abundance_short_df$com_id ==i),"LEMA"] / 
                                                  abundance_short_df[which( abundance_short_df$year== as.integer(year.levels[y]) &
                                                                              abundance_short_df$com_id ==i),"LEMA"]))
      }
    }
    
    abundance.GR.df.sp <-  abundance_short_df %>%
      select(-c("plot","subplot")) %>% 
      left_join(growth.ratio.df,by=c("year", 'com_id')) %>%
      filter(!is.na(growth.ratio) & !is.nan(growth.ratio) & !is.infinite(growth.ratio)) %>%
      filter(!growth.ratio==0)
    
    
    view(abundance.GR.df.sp)
    DataVec <- list(N= nrow(abundance.GR.df.sp),
                    S= length(Code.focal.list),
                    SpAbundance = abundance.GR.df.sp %>%
                      select(all_of(Code.focal.list)),
                    growth_ratio = abundance.GR.df.sp$growth.ratio ,
                    lambda_mean_obs = mean(abundance.GR.df.sp$growth.ratio),
                    year =  as.integer(factor(abundance.GR.df.sp$year,unique(abundance.GR.df.sp$year))),
                    Y= nlevels(as.factor(abundance.GR.df.sp$year)),
                    g_mean = sgerm.df$g.mean[which(sgerm.df$code.analysis == Code.focal)],
                    g_sd = sgerm.df$g.sd[which(sgerm.df$code.analysis == Code.focal)],
                    seed_s_mean = ssurv.df$s.mean[which(ssurv.df$code.analysis == Code.focal)],
                    seed_s_sd = ssurv.df$s.sd[which(ssurv.df$code.analysis == Code.focal)],
                    alpha_initial_mean = list_alpha_mean[[paste0(country,"_",Code.focal)]]["alpha_initial",],
                    alpha_slope_mean = list_alpha_mean[[paste0(country,"_",Code.focal)]]["alpha_slope",],
                    c_mean = list_alpha_mean[[paste0(country,"_",Code.focal)]]["c",],
                    N_opt_mean = list_alpha_mean[[paste0(country,"_",Code.focal)]]["N_opt",],
                    alpha_initial_sd = list_alpha_sd[[paste0(country,"_",Code.focal)]]["alpha_initial",],
                    alpha_slope_sd = list_alpha_sd[[paste0(country,"_",Code.focal)]]["alpha_slope",],
                    c_sd = list_alpha_sd[[paste0(country,"_",Code.focal)]]["c",],
                    N_opt_sd = list_alpha_sd[[paste0(country,"_",Code.focal)]]["N_opt",]
    )
    list.init <- function(...)list(N_opt= array(as.numeric(sapply(abundance.GR.df.sp%>%
                                                                    select(all_of(Code.focal.list)),median),
                                                           dim = DataVec$S)),
                                   lambda_mean= array(as.numeric( mean(abundance.GR.df.sp$growth.ratio),
                                                                  dim = 1))
    ) 
    
    
    Inv.Modelfit <- stan(file = paste0(home.dic,"code/Abundance_fit.stan") ,
                         #fit= PrelimFit, 
                         data = DataVec,
                         init ="random", # all initial values are 0 
                         control=list(max_treedepth=15),
                         warmup = 500,
                         iter = 1000, 
                         init_r = 2,
                         chains = 4,
                         seed= 1616) 
    
    save(file= paste0(project.dic,"results/Inv.Modelfit_",country,"_",Code.focal,".rds"),
         Inv.Modelfit)
    
    #load(paste0(home.dic,"results/Inv.Modelfit_",Code.focal,"_",year.int,".rds"))
    
    Inv.ModelfitPosteriors <- rstan::extract(Inv.Modelfit)
    
    save(file= paste0(project.dic,"results/Inv.ModelfitPosteriors",country,"_",Code.focal,".Rdata"),
         Inv.ModelfitPosteriors)
    #load(file= paste0(project.dic,"results/Inv.ModelfitPosteriors",country,"_",Code.focal,".Rdata"))
    
    
    Inv.Modelfit_loo <- rstan::loo(Inv.Modelfit,pars ="GR")
    
    save(file= paste0(project.dic,"results/Inv.ModelfitLOO",country,"_",Code.focal,".Rdata"),
         Inv.Modelfit_loo)
    
    print("Final Fit done")
    
    #---- 3.3. Final fit posterior check and behavior checks---- 
    
    ##### Diagnostic plots and post prediction 
    pdf(paste0(home.dic,"figures/Inv.Modelfit_",country,"_",Code.focal,".pdf"))
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    #source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    stan_post_pred_check(Inv.ModelfitPosteriors,"GR",
                         abundance.GR.df.sp$growth.ratio)
    hist(log(abundance.GR.df.sp$growth.ratio), breaks=150)
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(Inv.Modelfit)$summary[,"Rhat"],
         main = paste("Inversed Fit: Histogram of Rhat for",
                      Code.focal," and ",country))
    hist(summary(Inv.Modelfit)$summary[,"n_eff"],
         main = paste("Inversed Fit: Histogram of Neff for",
                      Code.focal," and ",country))
    
    # plot the corresponding graphs
    param <- c("lambda_mean","lambda_sd",
               "N_opt","disp_dev")
    
    trace <- stan_trace(Inv.Modelfit, pars=param,
                        inc_warmup = TRUE)
    print(trace)
    dens <- stan_dens(Inv.Modelfit, 
                      pars=param)
    print(dens)
    splot <- stan_plot(Inv.Modelfit, 
                       pars=param)
    print(splot)
    sampler_params <- get_sampler_params(Inv.Modelfit, inc_warmup = TRUE)
    summary(do.call(rbind, sampler_params), digits = 2)
    #pairs(Prelimfit, pars = param)
    
    dev.off()
    
    #---- 3.3. Extract Estimate----
    
    df_lambda_mean <- Inv.ModelfitPosteriors$lambda_mean %>% 
      as.data.frame() %>%
      setNames(Code.focal)
    
    df_lambda_sd <- Inv.ModelfitPosteriors$lambda_sd %>% 
      as.data.frame() %>%
      setNames(levels(as.factor(abundance.GR.df.sp$year)))
    
    write.csv(bind_cols(df_lambda_mean,df_lambda_sd),
              file=paste0(home.dic,"results/IntGR_",country,"_",Code.focal,".csv"))
  }
}
