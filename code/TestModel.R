#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import packages----
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Test model with empty data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2.1 Create empty data ----

DataVec <- list(N=200, 
                S=4,
                U= 3,
                Y=2, # number of years
                year = rep(c(1,2),each=100), # column for the year
                Nmax=matrix(0,ncol=4,nrow=2),
                Fecundity=rep(0,200),
                SpMatrix =matrix(0,ncol=4,nrow=200),
                Intra=c(0,0,1,0),
                run_estimation=1,
                tau0 = 1,
                slab_scale = sqrt(2),
                slab_df = 4
)
#---- 2.2. Run preliminary test ----

test.prelim <- stan(file = "code/Preliminary_fit.stan", 
     data = DataVec,
     init =  "random",
     warmup= 750,
     iter = 1500, 
     init_r = 1,
     chains = 3,
     control=list(max_treedepth=15),
     seed= 1644)

save(file= paste0("results/stan/Prelim_nodata.rds"),
     test.prelim)

test.prelim.posteriors <- rstan::extract(test.prelim)

#---- 3.3. Preliminary posterior check and behavior checks---- 

##### Diagnostic plots and post prediction 
pdf(paste0("figures/stan/Prelim_nodata.pdf"))
# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
# check the distribution of Rhats and effective sample sizes 
##### Posterior check
stan_post_pred_check(test.prelim.posteriors,"F_hat",Fecundity,
                     paste0("results/stan/PostFec_Prelim_Nodata.csv.gz")) 

#log_lik_2 <- loo::extract_log_lik(FinalFit, 
#                                  parameter_name = "F_sim", 
#                                  merge_chains = F)

#r_eff <- loo::relative_eff(exp(log_lik_2), cores = 2) 

# loo_1 <- loo::loo(log_lik_2 , r_eff = r_eff, cores = 2)

# print(loo_1)

# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(test.prelim)$summary[,"Rhat"],
     main = paste("Finat Fit: Histogram of Rhat for",
                  Code.focal," and ",year.int))
hist(summary(test.prelim)$summary[,"n_eff"],
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
trace <- stan_trace(test.prelim, 
                    pars=param ,
                    inc_warmup = TRUE)
print(trace)
dens <- stan_dens(test.prelim, 
                  pars=param )
print(dens)
splot <- stan_plot(test.prelim, 
                   pars=param)
print(splot)
sampler_params <- get_sampler_params(test.prelim, 
                                     inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
pairs(test.prelim, pars = c("Rho","alpha_initial","alpha_slope","c",
                         "alpha_initial_intra","alpha_slope_intra",
                         "c_intra","lambdas"))

dev.off()

