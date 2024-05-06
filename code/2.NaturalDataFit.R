
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
competition <- read.csv(paste0(home.dic,"data/competition_spain.csv"),
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))

species.neigh <- names(competition)
species.neigh  <- species.neigh[!species.neigh %in% c("day","month","year","plot","subplot","focal","fruit","comments",
                                                      "observers", "site", "treatment")]

competition.to.keep <- competition %>%  
  dplyr::select(all_of(c(species.neigh)))%>%
  mutate_at( species.neigh, as.numeric) %>%
  colSums(na.rm=T)
length(competition.to.keep[competition.to.keep == 0]) # check if any is 0

# change the format of the rows to numeric 
competition[species.neigh] <- sapply(competition[species.neigh],as.numeric)

# change na values to 0
competition[is.na(competition)] <- 0

focal.levels <- levels(as.factor(competition$focal))

#---- 1.3. Make Teasorus ----
plant_code <- read.csv(paste0(home.dic, "data/plant_code_spain.csv"),
                       header = T, stringsAsFactors = F, sep=",",
                       na.strings = c("","NA"))

competition.long.neigh <- competition %>%
  gather(any_of(plant_code$code.plant), key="code.plant",value="abundance") %>%
  left_join(plant_code) 

competition.long.focal <- competition.long.neigh %>%
  filter(focal == code.plant) %>%
  rename("conspecific"="abundance") %>%
  dplyr::select(focal,conspecific,fruit,seed,plot,subplot,year) 

competition.long <- competition.long.neigh %>%
  aggregate(abundance ~  focal + fruit + seed + family + plot + subplot + year, sum )%>%
  spread(family, abundance) %>% # change this ti have grass /forb or native/exotic
  right_join(competition.long.focal,
             multiple = "all") 

head(competition.long ) 

ggplot(competition.long[which(competition.long$focal =="LEMA"),],aes(y=seed, x=year,group=year) ) + 
  geom_boxplot()

ggplot(competition.long[which(competition.long$focal =="LEMA"),],aes(x=seed) ) + 
  geom_density()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 2. Run the model for each focal----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run.diagnostic =1
run.prelimfit = 0
run.finalfit = 0
grouping ="family"
year.int = "All"    
Code.focal = "LEMA"

for(Code.focal in c("LEMA")){ #focal.levels
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
    Y <- length(levels(as.factor(SpDataFocal$year)))
    year.levels = levels(as.factor(SpDataFocal$year))
    year.veclevels <- SpDataFocal$year # vector of levels
    
    
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
    
    
    DataVec <- list(N=N, 
                    S=S,
                    year= year.vec,
                    year.veclevels = year.veclevels,
                    year.levels =year.levels ,
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
    
    
    ##---- 3.2. Run Preliminary fit ----
    # Now run a fianl fit of the model to assess parameter 
    print("Preliminary Sparse fit begins")
    
    #install.packages("codetools")
    library("codetools")
    options(mc.cores = parallel::detectCores())
    
    list.init <- function(...)list(lambda_mean= array(as.numeric(max(log(Fecundity)), 
                                                                 dim = 1)),
                                   N_opt_mean= array(as.numeric(DataVec$Nmax, 
                                                                dim = DataVec$S)))
    
    if( run.prelimfit == 1){                               
      Prelimfit <- stan(file = paste0(home.dic,"code/chapt3/Preliminary_fit.stan"), 
                        data = DataVec,
                        init =  list.init,
                        warmup= 500,
                        iter = 1000, 
                        init_r = 1,
                        chains = 3,
                        control=list(max_treedepth=15),
                        seed= 1644)
      
      save(file= paste0(home.dic,"results/chapt3/Prelimfit_",Code.focal,"_",year.int,".rds"),
           Prelimfit)
    }
    
    load(paste0(home.dic,"results/chapt3/Prelimfit_",Code.focal,"_",year.int,".rds"))
    
    PrelimfitPosteriors <- rstan::extract(Prelimfit)
    
    save(file= paste0(home.dic,"results/chapt3/PrelimfitPosteriors",Code.focal,"_",year.int,".Rdata"),
         PrelimfitPosteriors)
    load(file= paste0(home.dic,"results/chapt3/PrelimfitPosteriors",Code.focal,"_",year.int,".Rdata"))
    
    print("Preliminary Fit done")
    
    #---- 3.3. Preliminary fit posterior check and behavior checks---- 
    
    ##### Diagnostic plots and post prediction 
    pdf(paste0(home.dic,"figure/chapt3/Prelimfit_",Code.focal,"_",year.int,".pdf"))
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    #source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    stan_post_pred_check(PrelimfitPosteriors,"F_hat",Fecundity,
                         paste0("results/chapt3/PostFec_",Code.focal,"_",year.int,".csv.gz"),
                         limx=max(Fecundity)+100) 
    
    
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(Prelimfit)$summary[,"Rhat"],
         main = paste("Prelim Fit: Histogram of Rhat for",
                      Code.focal," and ",year.int))
    hist(summary(Prelimfit)$summary[,"n_eff"],
         main = paste("Prelim Fit: Histogram of Neff for",
                      Code.focal," and ",year.int))
    
    # plot the corresponding graphs
    param <- c("alpha_initial","alpha_slope","c",
               "lambda_mean","N_opt_mean[1]",
               "alpha_initial_hat_tilde[1,1]",
               "initial_hat_shrinkage[1,1]",
               "alpha_slope_hat_tilde[1,1]",
               "slope_hat_shrinkage[1,1]",
               "c_hat_tilde[1,1]",
               "c_hat_shrinkage[1,1]")
    
    trace <- stan_trace(Prelimfit, pars=param,
                        inc_warmup = TRUE)
    print(trace)
    dens <- stan_dens(Prelimfit, 
                      pars=param)
    print(dens)
    splot <- stan_plot(Prelimfit, 
                       pars=param)
    print(splot)
    sampler_params <- get_sampler_params(Prelimfit, inc_warmup = TRUE)
    summary(do.call(rbind, sampler_params), digits = 2)
    #pairs(Prelimfit, pars = param)
    
    dev.off()
    #---- 3.3. Extract Inclusion ----
    
    
    Inclusion_alpha_initial <- matrix(data = rep(0,times=Y*S), nrow = Y, 
                                      ncol = S,
                                      dimnames = list(year.levels,
                                                      SpNames))
    
    Inclusion_alpha_slope <- matrix(data = 0, nrow = Y, 
                                    ncol = S,
                                    dimnames = list(year.levels,
                                                    SpNames))
    
    Inclusion_c <- matrix(data = 0, nrow = Y, 
                          ncol = S,
                          dimnames = list(year.levels,
                                          SpNames))
    alpha_initial_ys <- c()
    alpha_slope_ys <- c()
    alpha_c_ys <- c()
    
    IntLevel <- 0.1 #0.5 usually, 0.75 for Waitzia, shade
    for(y in 1:Y){
      for(s in 1:length(SpNames)){
        # ALPHA INITIAL
        # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
        alpha_initial_ys <- HDInterval::hdi(PrelimfitPosteriors$alpha_initial_hat[,y,s],
                                            credMass = IntLevel)
        if(alpha_initial_ys[1] > 0 | alpha_initial_ys[2] < 0){
          Inclusion_alpha_initial[y,s] <- 1
        }
        # ALPHA SLOPE
        # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
        alpha_slope_ys <- HDInterval::hdi(PrelimfitPosteriors$alpha_slope_hat[,y,s],
                                          credMass = IntLevel)
        if(alpha_slope_ys[1] > 0 | alpha_slope_ys[2] < 0){
          Inclusion_alpha_slope[y,s] <- 1
        }
        # ALPHA C
        # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
        alpha_c_ys <- HDInterval::hdi(PrelimfitPosteriors$c_hat[,y,s],
                                      credMass = IntLevel)
        if(alpha_c_ys[1] > 0 | alpha_c_ys[2] < 0){
          Inclusion_c[y,s] <- 1
        }
      }
    }
    
    Inclusion <- list(DataVec=DataVec, 
                      Inclusion_alpha_initial = Inclusion_alpha_initial,
                      Inclusion_alpha_slope = Inclusion_alpha_slope,
                      Inclusion_c = Inclusion_c
    )
    
    
    
    save(file= paste0(home.dic,
                      "results/chapt3/Inclusion",Code.focal,"_",year.int,".Rdata"),
         Inclusion)
    
    DataVec.final <- append(DataVec, 
                            list(Inclusion_alpha_initial = Inclusion_alpha_initial,
                                 Inclusion_alpha_slope = Inclusion_alpha_slope,
                                 Inclusion_c = Inclusion_c))
    
    save(file= paste0(home.dic,
                      "results/chapt3/Inclusion",
                      Code.focal,"_",year.int,".Rdata"),
         DataVec.final)
    
    ##---- 3.4. Run Final fit ----
    # Now run a fianl fit of the model to assess parameter 
    print("Final fit begins")
    
    
    set.seed(1616)
    std.error <- function(x) sd(x)/sqrt(length(x))
    
    list.init <- function(...)list(lambda_mean= array(as.numeric( abs(mean(PrelimfitPosteriors$lambda_mean))), dim = 1),
                                   N_opt_mean= array(as.numeric(abs(mean(PrelimfitPosteriors$N_opt_mean))), dim = DataVec$S),
                                   alpha_initial = array(as.numeric( colMeans(PrelimfitPosteriors$alpha_initial)), 
                                                         dim = DataVec$S),
                                   alpha_slope = array(as.numeric(colMeans(PrelimfitPosteriors$alpha_slope)), 
                                                       dim = DataVec$S),
                                   c = array(as.numeric(colMeans(PrelimfitPosteriors$c)),
                                             dim = DataVec$S),
                                   disp_dev = array(as.numeric(rnorm(1,mean = mean(PrelimfitPosteriors$disp_dev),
                                                                     sd= std.error(PrelimfitPosteriors$disp_dev))), 
                                                    dim = 1))
    
    if(run.finalfit == 1){
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores()) # to use the core at disposition 
      Finalfit <- stan(file = paste0(home.dic,"code/chapt3/Final_fit.stan") ,
                       #fit= PrelimFit, 
                       data = DataVec.final,
                       init =list.init, # all initial values are 0 
                       control=list(max_treedepth=15),
                       warmup = 500,
                       iter = 1000, 
                       init_r = 2,
                       chains = 4,
                       seed= 1616) 
      
      save(file= paste0(home.dic,"results/chapt3/Finalfit_",Code.focal,"_",year.int,".rds"),
           Finalfit)
    }
    
    load(paste0(home.dic,"results/chapt3/Finalfit_",Code.focal,"_",year.int,".rds"))
    
    FinalfitPosteriors <- rstan::extract(Finalfit)
    
    save(file= paste0(home.dic,"results/chapt3/FinalFitPosteriors",Code.focal,"_",year.int,".Rdata"),
         FinalfitPosteriors)
    load(file= paste0(home.dic,"results/chapt3/FinalFitPosteriors",Code.focal,"_",year.int,".Rdata"))
    
    
    print("Final Fit done")
    
    #---- 3.3. Final fit posterior check and behavior checks---- 
    
    ##### Diagnostic plots and post prediction 
    pdf(paste0(home.dic,"figure/chapt3/FinalFit_",Code.focal,"_",year.int,".pdf"))
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    #source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    stan_post_pred_check(FinalfitPosteriors,"F_hat",Fecundity,
                         paste0("results/chapt3/PostFec_Finalfit_",Code.focal,"_",year.int,".csv.gz"),
                         limx=max(Fecundity)+100) 
    
    
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(Finalfit)$summary[,"Rhat"],
         main = paste("Final Fit: Histogram of Rhat for",
                      Code.focal," and ",year.int))
    hist(summary(Finalfit)$summary[,"n_eff"],
         main = paste("Final Fit: Histogram of Neff for",
                      Code.focal," and ",year.int))
    
    # plot the corresponding graphs
    param <- c("alpha_initial[1]","alpha_slope[1]","c[1]",
               "lambda_mean[1]","lambda_sd[1]",
               "N_opt_mean[1]",
               "alpha_initial_hat[1,1]",
               "c_hat[1,1]",
               "alpha_slope_hat[1,1]")
    
    trace <- stan_trace(Finalfit, pars=param,
                        inc_warmup = TRUE)
    print(trace)
    dens <- stan_dens(Finalfit, 
                      pars=param)
    print(dens)
    splot <- stan_plot(Finalfit, 
                       pars=param)
    print(splot)
    sampler_params <- get_sampler_params(Finalfit, inc_warmup = TRUE)
    summary(do.call(rbind, sampler_params), digits = 2)
    #pairs(Prelimfit, pars = param)
    
    dev.off()
    #---- 3.3. Extract Estimate----
    #---- Generic parameters---
    
    df_parameter_n <-NULL
    
    df_lambda_mean <- FinalfitPosteriors$lambda_mean %>% 
      as.data.frame() %>%
      setNames(Code.focal)
    
    df_lambda_sd <- FinalfitPosteriors$lambda_sd %>% 
      as.data.frame() %>%
      setNames(year.levels)
    
    
    df_N_opt_mean <- FinalfitPosteriors$N_opt_mean %>% 
      as.data.frame() %>%
      setNames(Sp.names)
    
    
    generic_name <- c("alpha_initial","alpha_slope","c")
    generic_name <- generic_name[generic_name %in%  names(FinalfitPosteriors)]
    df_alpha_generic_param <- NULL
    
    for ( n in generic_name ){
      df_parameter_hat_n <- FinalfitPosteriors[[n]] %>% 
        as.data.frame() %>%
        setNames(SpNames) %>%
        mutate(parameter = n)
      
      df_alpha_generic_param  <- bind_rows(df_alpha_generic_param , df_parameter_hat_n)
    }
    
    # Species specific parameters
    df_alpha_initial <- matrix(nrow=2000, ncol=nrow(which(Inclusion_alpha_initial==1,arr.ind=T)))
    colnames_alpha_initial <- c()
    for( n in 1:nrow(which(Inclusion_alpha_initial==1,arr.ind=T))){
      coord.n <- which(Inclusion_alpha_initial==1,arr.ind=T)[n,]
      
      df_alpha_initial[,n] <- c(FinalfitPosteriors$alpha_initial_hat[,coord.n[1],coord.n[2]])
      colnames_alpha_initial[n] <- paste(year.levels[coord.n[1]],SpNames[coord.n[2]],sep="_")
      
    }
    colnames(df_alpha_initial) <- colnames_alpha_initial 
    
    
    df_alpha_slope <- matrix(nrow=2000, ncol=nrow(which(Inclusion_alpha_slope==1,arr.ind=T)))
    colnames_alpha_slope <- c()
    for( n in 1:nrow(which(Inclusion_alpha_slope ==1,arr.ind=T))){
      coord.n <- which(Inclusion_alpha_slope==1,arr.ind=T)[n,]
      
      df_alpha_slope[,n] <- FinalfitPosteriors$alpha_slope_hat[,coord.n[1],coord.n[2]]
      colnames_alpha_slope[n] <- paste(year.levels[coord.n[1]],SpNames[coord.n[2]],sep="_")
      
    }
    colnames(df_alpha_slope) <- colnames_alpha_slope
    
    
    df_alpha_c <- matrix(nrow=2000, ncol=nrow(which(Inclusion_c==1,arr.ind=T)))
    colnames_alpha_c <- c()
    for( n in 1:nrow(which(Inclusion_c==1,arr.ind=T))){
      coord.n <- which(Inclusion_c==1,arr.ind=T)[n,]
      
      df_alpha_c[,n] <- FinalfitPosteriors$c_hat[,coord.n[1],coord.n[2]]
      colnames_alpha_c [n] <- paste(year.levels[coord.n[1]],SpNames[coord.n[2]],sep="_")
      
    }
    colnames(df_alpha_c) <- colnames_alpha_c
    parameter = list(df_N_opt_mean = df_N_opt_mean, 
                     df_lambda_mean = df_lambda_mean, 
                     df_lambda_sd = df_lambda_sd,
                     df_alpha_generic_param = df_alpha_generic_param, 
                     df_alpha_c=df_alpha_c,
                     df_alpha_slope=df_alpha_slope,
                     df_alpha_initial =  df_alpha_initial)
    
    save(file= paste0(home.dic,"results/chapt3/Parameters_",Code.focal,"_",year.int,".Rdata"),
         parameter)
    
  }
}



