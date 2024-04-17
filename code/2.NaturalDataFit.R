
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

source("code/1.DataPrep.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 2. Run the model for each focal----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run.diagnostic =1
run.stan = 1
function.grouping = 1
grouping ="family"
year.levels = levels(as.factor(competition.long$year))
year.int = "All"    

for(Code.focal in "LEMA"){ #focal.levels
  for(year.int in c("All")){ #year.levels

  print(paste(Code.focal,year.int))
  
 # data for the focal
  if(year.int == "All"){
    SpDataFocal <- competition.long[which(competition.long$focal == Code.focal),] 
    
  if(length(levels(as.factor(  SpDataFocal$year)))>1){
    "All good for year"
  }
    }else{
  SpDataFocal <- competition.long[which(competition.long$focal == Code.focal & 
                                           competition.long$year == year.int),] 
  
  if(lenght(levels(as.factor(  SpDataFocal$year))) == 1 &
     length(levels(as.factor(  SpDataFocal$focal))) == 1){
    "All good for year and focal"
  }
  }
  
  names(SpDataFocal) <- tolower(names(SpDataFocal))
  SpDataFocal <-  SpDataFocal[!is.na(SpDataFocal$seed),]
  SpDataFocal[is.na(SpDataFocal)] <- 0
  # Next continue to extract the data needed to run the model. 
  N <- as.integer(nrow(SpDataFocal))
  Fecundity <- as.integer(SpDataFocal$seed) 
  year.vec <- as.integer(factor(SpDataFocal$year,unique(SpDataFocal$year))) # vector of levels
  Y <- length(levels(as.factor(  SpDataFocal$year)))
  if(grouping =="family"){
  specific.focal <- tolower(plant_code$family[which(plant_code$code.plant==Code.focal)])
  }else{
  specific.focal <- tolower(SpDataFocal$functional.group[1])
  }
  #---- 2. ABUDANCE MATRIX----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 2.1. Interaction (direct) matrix of plant with COMP ----
  
  # Now calculate the total number of plant species to use for the model, discounting
  #       any species columns with 0 abundance. Save a vector of the species names
  #       corresponding to each column for easy matching later.
  AllSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("focal","year","day","month",                
                                                              "seed","fruit","subplot","plot")]
  AllSpAbunds <- SpDataFocal %>% 
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
      #if(SpTotals[s] < 5){
      #  SpMatrix[,S+1] <- SpMatrix[,S+1] + AllSpAbunds[,s]
     # }else{
     SpMatrix[,i] <- AllSpAbunds[,s]
     i <- i + 1
    }else{next}
  }
   
   SpNames <- c(AllSpNames[SpToKeep])
   colnames(SpMatrix) <-  SpNames
   Intra <- ifelse(SpNames == "conspecific", 1, 0)
   SpMatrix[,which(SpNames ==specific.focal)] <- SpMatrix[,which(SpNames == specific.focal)] - SpMatrix[,"conspecific"]  # remove conspecific obs from focal family
  
  # max fecundity
  Nmax <- c(SpMatrix[which.max(Fecundity),])

  # Upper bound intrinsic fecundity
  U <- ceiling(log(max(Fecundity)))
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. BAYES FIT----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##---- 3.1. Set up summary interactions df and parameters ---- 
  
  run_estimation <- 1

  DataVec <- list(N=N, 
                  S=S,
                  year= year.vec,
                  Y = Y,
                  U= U,
                  Nmax=Nmax,
                  Fecundity=Fecundity,
                  SpMatrix =SpMatrix,
                  Intra=Intra,
                  run_estimation=1,
                  tau0 = 1,
                  slab_scale = sqrt(2),
                  slab_df = 4
                  )
  
  
  ##---- 3.2. Run  final fit ----
  # Now run a fianl fit of the model to assess parameter 
  print("Preliminary Sparse fit begins")
  
  #install.packages("codetools")
  library("codetools")
  options(mc.cores = parallel::detectCores())
 
  if( run.stan == 1){                               
  Prelimfit <- stan(file = "code/Preliminary_fit.stan", 
                   data = DataVec,
                   init =  "random",
                   warmup= 500,
                   iter = 1000, 
                   init_r = 1,
                   chains = 3,
                   control=list(max_treedepth=15),
                   seed= 1644)
  
  save(file= paste0("results/stan/Prelimfit_",Code.focal,"_",year.int,".rds"),
       Prelimfit)
  }

    
  #load(paste0("results/stan/FinalFit_",Code.focal,"_",alpha.function,".rds"))
  
 load(paste0("results/stan/Prelimfit_",Code.focal,"_",year.int,".rds"))

 PrelimfitPosteriors <- rstan::extract(Prelimfit)
  
  print("Final Fit done")
  
  #---- 3.3. Final fit posterior check and behavior checks---- 

  ##### Diagnostic plots and post prediction 
  pdf(paste0("figures/stan/Prelimfit_",Code.focal,"_",year.int,".pdf"))
  # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
  source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
  # check the distribution of Rhats and effective sample sizes 
  ##### Posterior check
  stan_post_pred_check(PrelimfitPosteriors,"F_hat",Fecundity,
                       paste0("results/stan/PostFec_",Code.focal,"_",year.int,".csv.gz")) 

  
  # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
  hist(summary(Prelimfit)$summary[,"Rhat"],
       main = paste("Finat Fit: Histogram of Rhat for",
                    Code.focal," and ",year.int))
  hist(summary(Prelimfit)$summary[,"n_eff"],
       main = paste("Finat Fit: Histogram of Neff for",
                    Code.focal," and ",year.int))
  
  # plot the corresponding graphs
  param <- c("Rho","alpha_initial","alpha_slope","c",
             "alpha_initial_intra","alpha_slope_intra",
             "c_intra","lambdas",
             "alpha_initial_hat_tilde[1,1]",
             "initial_hat_shrinkage[1,1]",
             "alpha_slope_hat_tilde[1,1]",
             "slope_hat_shrinkage[1,1]",
             "c_hat_tilde[1,1]",
             "c_hat_shrinkage[1,1]")
  
  trace <- stan_trace(Prelimfit, pars=param,
                      inc_warmup = TRUE)
  print(trace)
  dens <- stan_dens(FinalFit, 
                    pars=param)
  print(dens)
  splot <- stan_plot(FinalFit, 
                     pars=param)
  print(splot)
  sampler_params <- get_sampler_params(FinalFit, inc_warmup = TRUE)
  summary(do.call(rbind, sampler_params), digits = 2)
  pairs(FinalFit, pars = param)
  
  dev.off()
  #---- 3.3. Extract Inclusion ----
  

  Inclusion_alpha_initial <- matrix(data = 0, nrow = Y, 
                                    ncol = S,
                                    dimnames=list(c(""),
                                                  SpNames))
  
  Inclusion_alpha_slope <- matrix(data = 0, nrow = Y, 
                                    ncol = S,
                                    dimnames=list(c(""),
                                                  SpNames))
  
  Inclusion_c <- matrix(data = 0, nrow = Y, 
                                  ncol = S,
                                  dimnames=list(c(""),
                                                SpNames))
  alpha_initial_ys <- c()
  alpha_slope_ys <- c()
  alpha_c_ys <- c()
  
  IntLevel <- 0.2 #0.5 usually, 0.75 for Waitzia, shade
  for(y in 1:Y){
    for(s in 1:length(SpNames)){
   # ALPHA INITIAL
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
      alpha_initial_ys <- HDInterval::hdi(PrelimPosteriors$alpha_initial_hat[y,s],
                               credMass = IntLevel)
    if(alpha_initial_ys[1] > 0 | alpha_initial_ys[2] < 0){
      Inclusion_alpha_initial[y,s] <- 1
    }
  # ALPHA SLOPE
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    alpha_slope_ys <- HDInterval::hdi(PrelimPosteriors$alpha_slope_hat[y,s],
                                        credMass = IntLevel)
    if(alpha_slope_ys[1] > 0 | alpha_slope_ys[2] < 0){
      Inclusion_alpha_slope[y,s] <- 1
    }
  # ALPHA C
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    alpha_c_ys <- HDInterval::hdi(PrelimPosteriors$c_hat[y,s],
                                      credMass = IntLevel)
    if(alpha_c_ys[1] > 0 | alpha_c_ys[2] < 0){
      Inclusion_c[y,s] <- 1
      }
    }
  }
  
  Inclusion <- as.data.frame(bind.rows(bind.rows(Inclusion_alpha_initial,Inclusion_alpha_slope),
            Inclusion_c)) %>%
      colnames(SpNames) %>%
      mutate( parameter = rep(c("alpha_initial","alpha_slope","c"),times=3),
              year = rep(year.vec,times=3)) 
  write.csv(Inclusion,
            paste0("results/Inclusion",Code.focal,"_",year.int))
  
  }
}

#---- 3.4. Extract Parameters ----

assign(paste0("Parameters_",Code.focal,"_",alpha.function),
       list(DataVec = DataVec,
            alpha_value= PrelimPosteriors$alpha_value,      
            lambda =  FinalPosteriors$lambda_ei,
            alpha_slope= FinalPosteriors$alpha_slope_ei,
            alpha_init =FinalPosteriors$alpha_init,
            alpha_c = FinalPosteriors$c_ei
       ))

save(list =paste0("Parameters_",Code.focal,"_",alpha.function),
     file = paste0("results/stan/Parameters_",Code.focal,"_",alpha.function,".RData"))

#---- 3.5. Extraction fecundity---

assign(paste0("Fsim_",Code.focal,"_",alpha.function),
       list(DataVec = DataVec,
            F_sim=FinalPosteriors$F_sim
       ))

save(list =paste0("Fsim_",Code.focal,"_",alpha.function),
     file = paste0("results/stan/Fsim_",Code.focal,"_",alpha.function,".RData"))


