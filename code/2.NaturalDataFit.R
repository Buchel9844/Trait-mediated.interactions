
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

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
home.dic <- ""
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"

#---- 1.2. Import the competitive data ----
#setwd("~/Documents/Projects/Facilitation_gradient")
#source(paste0(home.dic,"code/1.DataPrep.R"))
# see competition.spain_long
#write.csv(competition.spain_long,
#          file=paste0(home.dic,"results/competition.spain_long.csv"))

competition.spain_long <- read.csv(paste0(home.dic,"results/competition.spain_long.csv"))      
species.spain <- levels(as.factor(competition.spain_long$focal.analysis))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 2. Run the model for each focal----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run.diagnostic =1
run.fit = 1
#grouping ="family"
#year.int = "All"    
#Code.focal = "LEMA"

load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))

for(country in c("aus","spain")){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species.",country)]]
  for(Code.focal in Code.focal.list){ #focal.levels
    print(paste(Country,Code.focal))
    
    competition_focal_neigh <- get(paste0("clean.data.",country))[[paste0("competition_",country)]] %>%
      filter(focal == Code.focal) %>%
      gather(any_of(Code.focal.list),
             key="code.analysis",value="abundance") %>%
      aggregate(abundance ~ code.analysis + year, sum) %>%
      spread(code.analysis, abundance)
    
    competition.spain_focal[competition.spain_focal==0] <- NA
    names.to.keep <- names(competition_focal_neigh)[
      colSums(competition_focal_neigh)>10 & # more than 10 individuals found across all sample in the neighbourhood of the focal
        !is.na(colSums(competition_focal_neigh)) 
      & names(competition_focal_neigh) %in% Code.focal.list] 
    
    # filter not working so doing it that hard way 
    
    competition.to.keep <- get(paste0("clean.data.",country))[[paste0("competition_",country)]] %>%
      filter(focal == Code.focal) %>%
      gather(any_of(names.to.keep),
             key="code.analysis",value="abundance")  %>%
      aggregate(abundance ~ year + plot +  focal + seed + code.analysis, sum) %>%
      spread(code.analysis, abundance) %>%
      dplyr::filter(seed < 5000)
    
    competition.to.keep <-  competition.to.keep[!is.na(competition.to.keep$seed),]
    competition.to.keep[is.na(competition.to.keep)] <- 0
    # Next continue to extract the data needed to run the model. 
    N <- as.integer(nrow(competition.to.keep))
    Fecundity <- round(competition.to.keep$seed) 
    year.vec <- as.integer(factor(competition.to.keep$year,unique(competition.to.keep$year))) # vector of levels
    
    #---- 2. ABUDANCE MATRIX----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #---- 2.1. Interaction (direct) matrix of plant with COMP ----
    
    # Now calculate the total number of plant species to use for the model, discounting
    #       any species columns with 0 abundance. Save a vector of the species names
    #       corresponding to each column for easy matching later.
    AllSpNames <- names(competition.to.keep)[!names(competition.to.keep) %in% c("focal","year","day","month",                
                                                                                "seed","fruit","subplot","plot")]
    AllSpAbunds <- competition.to.keep %>% 
      dplyr::select(all_of(c(AllSpNames)))%>%
      mutate_at( AllSpNames, as.numeric)
    
    SpTotals <- colSums(AllSpAbunds)
    
    SpToKeep <- SpTotals > 0
    sum(SpTotals)
    NamesSpToKeep <- names(SpTotals[SpToKeep])
    Stotal <- sum(SpToKeep)
    S <- Stotal
    SpMatrix <- NULL
    SpMatrix <- matrix(NA, nrow = N, ncol = S)
    i <- 1
    for(s in 1:length(AllSpNames)){
      if(SpToKeep[s] == 1){
        SpMatrix[,i] <- AllSpAbunds[,s]
        i <- i + 1
      }else{next}
    }
    
    SpNames <- c(AllSpNames[SpToKeep])
    colnames(SpMatrix) <-  SpNames
    Intra <- ifelse(SpNames == Code.focal, 1, 0)
    #SpMatrix[,which(SpNames ==specific.focal)] <- SpMatrix[,which(SpNames == specific.focal)] - SpMatrix[,"conspecific"]  # remove conspecific obs from focal family
    
    # max fecundity
    Nmax <- c(SpMatrix[which.max(Fecundity),])
    
    # Upper bound intrinsic fecundity
    FMax <- ceiling(log(median(Fecundity)))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #---- 3. BAYES FIT----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##---- 3.1. Set up summary interactions df and parameters ---- 
    
     DataVec <- list(N=N, 
                    S=S,
                    FMax = FMax,
                    Nmax=Nmax,
                    year= year.vec,
                    Y = nlevels(as.factor(year.vec)),
                    Fecundity = Fecundity,
                    SpMatrix =SpMatrix,
                    Intra=Intra,
                    run_estimation=1,
                    tau0 = 1,
                    slab_scale = sqrt(2),
                    slab_df = 4
    )
    
    
    ##---- 3.2. Run Final fit ----
    # Now run a fianl fit of the model to assess parameter 
    print("Final fit begins")
    
    set.seed(1616)
    std.error <- function(x) sd(x)/sqrt(length(x))
    
    list.init <- function(...)list(#lambda_mean= array(as.numeric(mean(log(Fecundity)), 
      #                               dim = 1)),
      N_opt_mean= array(as.numeric(DataVec$Nmax, 
                                   dim = DataVec$S)))
    
    if(run.fit == 1){
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores()) # to use the core at disposition 
      Modelfit <- stan(file = paste0(home.dic,"code/Ricker_fit.stan") ,
                       #fit= PrelimFit, 
                       data = DataVec,
                       init =list.init, # all initial values are 0 
                       control=list(max_treedepth=15),
                       warmup = 500,
                       iter = 1000, 
                       init_r = 2,
                       chains = 4,
                       seed= 1616) 
      
      save(file= paste0(home.dic,"results/Modelfit_",Code.focal,"_",year.int,".rds"),
           Modelfit)
    }
    
    load(paste0(home.dic,"results/Modelfit_",Code.focal,"_",year.int,".rds"))
    
    ModelfitPosteriors <- rstan::extract(Modelfit)
    
    save(file= paste0(home.dic,"results/ModelfitPosteriors",Code.focal,"_",year.int,".Rdata"),
         ModelfitPosteriors)
    load(file= paste0(home.dic,"results/ModelfitPosteriors",Code.focal,"_",year.int,".Rdata"))
    
    
    Modelfit_loo <- rstan::loo(Modelfit,pars ="F_sim")
    
    save(file= paste0(home.dic,"results/ModelfitLOO",Code.focal,"_",year.int,".Rdata"),
         Modelfit_loo)
    
    print("Final Fit done")
    
    #---- 3.3. Final fit posterior check and behavior checks---- 
    
    ##### Diagnostic plots and post prediction 
    pdf(paste0(home.dic,"figures/Modelfit_",Code.focal,"_",year.int,".pdf"))
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    #source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    stan_post_pred_check(ModelfitPosteriors,"F_hat",Fecundity,
                         paste0("results/PostFec_Modelfit_",Code.focal,"_",year.int,".csv.gz"),
                         limx=max(Fecundity)+100) 
    
    
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(Modelfit)$summary[,"Rhat"],
         main = paste("Final Fit: Histogram of Rhat for",
                      Code.focal," and ",year.int))
    hist(summary(Modelfit)$summary[,"n_eff"],
         main = paste("Final Fit: Histogram of Neff for",
                      Code.focal," and ",year.int))
    
    # plot the corresponding graphs
    param <- c("alpha_initial","alpha_slope","c",
               "lambda_mean","lambda_sd",
               "N_opt_i")
    
    trace <- stan_trace(Modelfit, pars=param,
                        inc_warmup = TRUE)
    print(trace)
    dens <- stan_dens(Modelfit, 
                      pars=param)
    print(dens)
    splot <- stan_plot(Modelfit, 
                       pars=param)
    print(splot)
    sampler_params <- get_sampler_params(Modelfit, inc_warmup = TRUE)
    summary(do.call(rbind, sampler_params), digits = 2)
    #pairs(Prelimfit, pars = param)
    
    dev.off()
    #---- 3.3. Extract Estimate----
    #---- Generic parameters---
    
    df_parameter_n <-NULL
    
    df_lambda_mean <- ModelfitPosteriors$lambda_mean %>% 
      as.data.frame() %>%
      setNames(Code.focal)
    
    df_lambda_sd <- ModelfitPosteriors$lambda_sd %>% 
      as.data.frame() %>%
      setNames(levels(as.factor(competition.to.keep$year)))
    
    
    df_N_opt <- ModelfitPosteriors$N_opt %>% 
      as.data.frame() %>%
      setNames(SpNames)
    
    
    generic_name <- c("alpha_initial","alpha_slope","c")
    generic_name <- generic_name[generic_name %in%  names(ModelfitPosteriors)]
    df_alpha_generic_param <- NULL
    
    for ( n in generic_name ){
      df_parameter_hat_n <- ModelfitPosteriors[[n]] %>% 
        as.data.frame() %>%
        setNames(SpNames) %>%
        mutate(parameter = n)
      
      df_alpha_generic_param  <- bind_rows(df_alpha_generic_param , df_parameter_hat_n)
    }
   
    parameter = list(df_N_opt = df_N_opt, 
                     df_lambda_mean = df_lambda_mean, 
                     df_lambda_sd = df_lambda_sd,
                     df_alpha_generic_param = df_alpha_generic_param
                     )
    
    save(file= paste0(home.dic,"results/Parameters_",Code.focal,"_",year.int,".Rdata"),
         parameter)
    
  }
}



