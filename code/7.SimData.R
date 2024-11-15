
#Inversed modelling to predict lambda 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(cli)
library(rstan)
#install.packages("HDInterval")
library("HDInterval")
#install.packages("tidyverse")
library(tidyverse)
library(tibble)
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
library(truncnorm)
#setwd("/home/lbuche/Eco_Bayesian")
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"

home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2.Test model with simulated data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alpha_function4  <- function(Amin, Aslopes,c,N,N0){
  e = exp(Aslopes*(N-N0)) # c is stretching the graph horizontally 
  #if((N-N0) >10 ){
  #  e = exp(-Aslopes*(10))
  #}
  a = c*(1-e) #stretching the graph vertically
  d = Amin
  alpha = (a/(1 + e)) + d
  
  return(alpha)
}
Ricker_solution_NatData <- function(state,
                                    pars,
                                    neigh.vec,
                                    N,
                                    numb_plot,
                                    numb_year,
                                    widthplot =25) {
  #g_low <- pars$g_low[1] # germination rate 
  #g_high<- pars$g_high[1] # germination rate 
  g_init<- pars$g_init[1] # germination rate 
  s <- pars$s[1] #seed survival
  beta <- pars$beta
  beta_g <- pars$beta_g
  lambda.mean <- pars$lambda[1] # intrinsic growth rate
  df <- state
  PDSI <- rep(truncnorm::rtruncnorm(numb_year,0,1,a=-1,b=1),each= numb_plot) #scaled values
  g <- c()
  for(y in 1:(N)){
    g[y] <- g_init + beta_g*exp(PDSI[y])
  }
  for(y in 1:(N-numb_plot)){
    lambda <- lambda.mean + beta*exp(PDSI[y])
    Nt1 <- c()
    Fec <- c()
    Nt <-  df[y,]  # species i densities
    Interaction_effect <- c()
    for(n in 1:neigh.vec){
      Nmax <- pars$Nmax[n] # density at which fecundity is max - effect of neighbors is 0
      a_initial <- pars$a_init[n] # which int.function
      a_slope <- pars$a_slope[n]  # which int.function
      c <- pars$c[n]  # which int.function
      
      a <- alpha_function4(a_initial, a_slope,c,Nt[n]*(15/widthplot), Nmax)
      
      Interaction_effect[n] <- as.numeric(a*Nt[n]*(15/widthplot))
    }
    Fec = exp(sum(Interaction_effect))
    
    Nt1 <- (1-g[y])*s*Nt$seed.total.conspecific+ Fec*lambda*Nt$conspecific
    if(Nt1 <0 ){
      Nt1 = 1
    }
    df$conspecific[y+numb_plot] <-Nt1*g[y+numb_plot] #rpois(1,lambda=Nt1*g)
    
    df$seed.total.conspecific[y+numb_plot] <- Nt1
    
    df$PDSI[y] <- PDSI[y]
    df$lambda[y] <- lambda
    df$g[y] <- g[y]
    df$Fec[y] <- log(Fec)
  }
  
  return(df)
}
parameter.df <- NULL
#i = 1
#numb_plot =10
#numb_year=7
for(i in 1:100){
  # ### Select set values ###
  S <- 5
  alpha_initial_mean = truncnorm::rtruncnorm(S, mean=0,sd=0.2,a=-1,b=1)
  alpha_slope_mean = -abs(truncnorm::rtruncnorm(S, mean=0.2,sd=0.2,a=-1,b=0))
  c_mean = -abs(truncnorm::rtruncnorm(S, mean=0,sd=0.1,a=-1,b=0))
  N_opt_mean = floor(abs(rnorm(S, mean=0,sd=2)))
  N_opt_mean[S]<- 0
  alpha_initial_mean[S] <- -0.2
  sign.beta <- sample(c(-1,1),1)
  pars <- list(a_initial=alpha_initial_mean,
               a_slope=alpha_slope_mean,
               c=c_mean,
               Nmax =N_opt_mean,
               g_init =truncnorm::rtruncnorm(1, mean=0.5,sd=0.2,a=0.3,b=0.9),
               beta_g = truncnorm::rtruncnorm(1, mean=0.1,sd=0.1,a=0.001,b=0.5)*sign.beta,
               s=truncnorm::rtruncnorm(1, mean=0.5,sd=0.2,a=0.3,b=0.7),
               lambda = truncnorm::rtruncnorm(1, mean=20,sd=10,a=15,b=30),
               beta=truncnorm::rtruncnorm(1, mean=2,sd=2,a=0.5,b=5)*sign.beta)
  for(numb_plot in c(10,30,100)){
    for(numb_year in c(7,10,15)){
      print(paste(i,numb_plot,numb_year))
      # ### CREATE SIMULATED DATA ###
      # set up an empty matrix with columns for each neighbours
      # create a matrix of observations for each focal species
      S_obs <- numb_plot* (numb_year+100)
      K_Nmat <- matrix(ncol =  S_obs)
      SpAbundance <- matrix(data = 0, nrow = S_obs, ncol = S-1)
      
      for(s in 1:S-1){ 
        SpAbundance[,s] <- rpois(S_obs, lambda = abs(rnorm(1,mean=0,sd=5)))
      }
      
      state = SpAbundance %>%
        as.data.frame() %>%
        mutate(conspecific = c(rep(10,times=numb_plot),
                               rep(0,S_obs-numb_plot)),
               seed.total.conspecific = 0,
               time = rep(1: (numb_year+100),each=numb_plot),
               plot = rep(1: numb_plot,times=(numb_year+100)))
      
      
      SpAbundance_focal <- Ricker_solution_NatData(state,
                                                   pars,
                                                   neigh.vec = S,
                                                   N = S_obs,
                                                   numb_plot=numb_plot,
                                                   numb_year=numb_year+100,
                                                   widthplot =25)
      SpAbundance_focal <- SpAbundance_focal %>%
        filter(time > 100) %>%
        mutate(time = time - 100)
      #View(SpAbundance_focal)
      #head(SpAbundance_focal)
      #ggplot(SpAbundance_focal, aes(x=as.factor(time),y=conspecific)) + geom_boxplot()
      #plot(density(SpAbundance_focal$conspecific[which(SpAbundance_focal$time <10)],na.rm=T))
      #plot(density(SpAbundance_focal$seed.total.conspecific,na.rm=T))
      
      #SpAbundance_focal %>%
      # ggplot(aes(y=conspecific, x=exp(PDSI))) + geom_point() +
      # geom_smooth()
      SpAbundance_focal.glm <- glm(lambda~PDSI,data=SpAbundance_focal)
      
      sp.name<- c(paste0("V",1:(S-1)),"conspecific")
      DataVec <- list(N=  numb_plot* numb_year,
                      S= S,
                      plot_width =25,
                      ratio_scale = 15/25,
                      SpAbundance =   SpAbundance_focal %>%
                        dplyr::select(all_of(sp.name)),
                      SpAbundance_med = SpAbundance_focal %>%
                        group_by(time) %>%
                        summarise_at(sp.name, function(x) median(x[!x==0 & !is.na(x)]))%>%
                        fill(sp.name, .direction=c("downup"))%>%
                        dplyr::select(conspecific),
                      # obs_spI = round(as.vector(unlist(abundance_df_analysis[,Code.focal]))),
                      obs_spI = as.vector(unlist(round(SpAbundance_focal[,"conspecific"]))),
                      lambda_mean_prior =mean(SpAbundance_focal$lambda) + sd(SpAbundance_focal$lambda),
                      lambda_sd_prior =10,
                      focal_pos = which(sp.name == "conspecific"), 
                      numb_year= numb_year,
                      numb_plot= numb_plot,
                      year =  SpAbundance_focal$time,
                      plot =  SpAbundance_focal$plot,
                      PDSI_mean = c(SpAbundance_focal %>%
                                      dplyr::select(PDSI,time) %>% unique() %>% dplyr::select(PDSI) %>%unlist()),
                      g_vec = c(SpAbundance_focal %>%
                                  dplyr::select(g,time) %>% unique() %>% dplyr::select(g) %>%unlist()),
                      seed_s =0.5,
                      alpha_initial_mean = alpha_initial_mean,
                      alpha_slope_mean = alpha_slope_mean,
                      c_mean = c_mean ,
                      N_opt_mean = N_opt_mean,
                      beta_prior =SpAbundance_focal.glm$coefficients[2],
                      beta_up = 30,#-median(SpAbundance_focal$lambda)/min(SpAbundance_focal$PDSI),
                      beta_low = -30)#-median(SpAbundance_focal$lambda)/max(SpAbundance_focal$PDSI))
      
      list.init <- function(...)list(beta= array(as.numeric(SpAbundance_focal.glm$coefficients[2]),
                                                 dim = 1),
                                     g_beta= array(as.numeric(0.1*sign.beta),
                                                   dim = 1),
                                     lambda_mean =array(as.numeric(DataVec$lambda_mean_prior),
                                                        dim = 1),
                                     g_init =array(as.numeric(0.5),
                                                   dim = 1),
                                     seed_s =array(as.numeric(0.5),
                                                   dim = 1))
      # g_init = array(as.numeric(sgerm.df$g.mean[which(sgerm.df$code.analysis == Code.focal)]),
      #             dim = 1))
      
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores()) # to use the core at disposition 
      
      Inv.Modelfit <- stan(file = paste0(home.dic,"code/Abundance_fit.stan") ,
                           #fit= PrelimFit, 
                           data = DataVec,
                           #init ="random", # all initial values are 0 
                           init=list.init,
                           control=list(max_treedepth=15),
                           warmup = 500,
                           iter = 1000, 
                           init_r = 2,
                           chains = 3)  
      Inv.ModelfitPosteriors <- rstan::extract(Inv.Modelfit)
      if(i==1){
        save(file= paste0(project.dic,"results/Inv.Modelfit_test_dataVec_",i,"_",numb_year,"_",numb_plot,".RData"),
             DataVec )
        save(file= paste0(project.dic,"results/Inv.Modelfit_test_",i,"_",numb_year,"_",numb_plot,".rds"),
             Inv.Modelfit)
        #load(paste0(project.dic,"results/Inv.Modelfit_test.rds"))
        #load(paste0(project.dic,"results/Inv.Modelfit_",country,"_",Code.focal,".rds"))
        source(paste0(home.dic,"code/stan_modelcheck_rem.R")) # call the functions to check diagnistic plots
        pdf(paste0(project.dic,"figures/Inv.Modelfit_test",i,"_",numb_year,"_",numb_plot,".pdf"))
        
        hist(Inv.ModelfitPosteriors$seedtotal,
             breaks = 150)
        
        hist( SpAbundance_focal$seed.total.conspecific,breaks = 150)
        # check the distribution of Rhats and effective sample sizes 
        ##### Posterior check
        stan_post_pred_check_pois_glow_ghigh(Inv.ModelfitPosteriors,"seedtotal",
                                             DataVec$obs_spI,
                                             gvec = DataVec$g_vec,
                                             numb_plot,
                                             numb_year) 
        
        param <- c("lambda_mean","beta","g_init","g_beta","seed_s")
        
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
        #pairs(Inv.Modelfit, pars = param)
        
        dev.off()
      }
      
      parameter.df.i <- left_join(as.data.frame(pars) %>%
                                    dplyr::select(s,g_init,beta_g,lambda,beta) %>%
                                    unique() %>%
                                    gather(key="parameter",value="set.value"),
                                  
                                  as.data.frame(matrix(c(quantile(Inv.ModelfitPosteriors$seed_s,c(0.025,0.5,0.975)),
                                                         quantile(Inv.ModelfitPosteriors$g_init,c(0.025,0.5,0.975)),
                                                         quantile(Inv.ModelfitPosteriors$g_beta,c(0.025,0.5,0.975)),
                                                         quantile(Inv.ModelfitPosteriors$lambda_mean,c(0.025,0.5,0.975)),
                                                         quantile(Inv.ModelfitPosteriors$beta,c(0.025,0.5,0.975))),
                                                       ncol=3, byrow=T,dimnames =list(c(),
                                                                                      c("Q2.5","Q50","Q97.5")))) %>%
                                    mutate(parameter = c("s","g_init","beta_g","lambda","beta"))) %>%
        mutate(sim = i,
               numb_year=numb_year,
               numb_plot=numb_plot)
      
      parameter.df <- bind_rows(parameter.df,parameter.df.i)
      
      write.csv(parameter.df,
                paste0(project.dic,"results/parameter.df.csv"))  
    }
  }
}
parameter.df <- read.csv(paste0(project.dic,"results/parameter.df.csv"))  %>%
  mutate(parameter.2 = case_when((parameter=="s")~ "seed survival",
                                 (parameter=="lambda")~ "intrinsic performance",
                                 (parameter=="beta")~ "Effect of exp(PDSI)\n on the intrinsic performance",
                                 (parameter=="g_init")~ "seed germination",
                                 (parameter=="beta_g")~ "Effect of exp(PDSI)\n on germination",
                                 T~ parameter))

parameter.plot <- parameter.df %>%
  mutate(diff = Q50  -set.value) %>%
  mutate(parameter.2 = factor(parameter.2,
                              levels=c("seed survival","intrinsic performance","Effect of exp(PDSI)\n on the intrinsic performance",
                                       "seed germination","Effect of exp(PDSI)\n on germination"))) %>%
  ggplot(aes(x=diff,y=1)) +
  geom_density_ridges(quantile_lines=T,
                      quantiles = c(0.025,0.5, 0.975)) +
  geom_vline(xintercept=0, color="red") +
  facet_wrap(.~parameter.2, scale="free") +
  labs(y= "",
       x="Difference between the median of the estimated parameter and its corresponding set value") +
  theme_bw() +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=16),
         legend.text=element_text(size=12),
         legend.title=element_text(size=12),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=12, angle=0, hjust=1),
         axis.text.y= element_text(size=12, angle=20, hjust=1),
         axis.title.x= element_text(size=16),
         axis.title.y= element_text(size=16),
         title=element_text(size=16))
parameter.plot

ggsave(parameter.plot,
       file=paste0(home.dic,"figures/parameter.plot.pdf"))


parameter.plot.n.obs <- parameter.df %>%
  mutate(diff = Q50  -set.value) %>%
  mutate(parameter.2 = factor(parameter.2,
                              levels=c("seed survival","intrinsic performance","Effect of exp(PDSI)\n on the intrinsic performance",
                                       "seed germination","Effect of exp(PDSI)\n on germination"))) %>%
  ggplot(aes(x=diff, # already removed the parameter to recover
             y=as.factor(numb_plot))) +
  facet_grid(.~parameter.2, scale="free") +
  geom_density_ridges(quantile_lines = TRUE,
                      scale=1,
                      quantiles = c(0.025,0.5, 0.975)) +
  geom_vline(xintercept=0, color="red") +
  labs(y= "Number of plots/observations within each year",
       x="Difference between the median of the estimated parameter and its corresponding set value") +
  theme_bw() +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=14),
         legend.text=element_text(size=12),
         legend.title=element_text(size=12),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=12, angle=0, hjust=1),
         axis.text.y= element_text(size=12, angle=20, hjust=1),
         axis.title.x= element_text(size=16),
         axis.title.y= element_text(size=16),
         title=element_text(size=16),
         plot.margin = unit(c(1,1,1,1), "cm"))
parameter.plot.n.obs
ggsave(parameter.plot.n.obs,
       file=paste0(home.dic,"figures/parameter.plot.n.obs.pdf"))

parameter.plot.year <- parameter.df %>%
  mutate(diff = Q50  -set.value) %>%
  mutate(parameter.2 = factor(parameter.2,
                              levels=c("seed survival","intrinsic performance","Effect of exp(PDSI)\n on the intrinsic performance",
                                       "seed germination","Effect of exp(PDSI)\n on germination"))) %>%
  ggplot(aes(x=diff, # already removed the parameter to recover
             y=as.factor(numb_year))) +
  facet_grid(.~parameter.2, scale="free") +
  geom_density_ridges(quantile_lines = TRUE,
                      scale=1,
                      quantiles = c(0.025,0.5, 0.975)) +
  geom_vline(xintercept=0, color="red") +
  labs(y= "Number of year",
       x="Difference between the median of the estimated parameter and its corresponding set value") +
  theme_bw() +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=16),
         legend.text=element_text(size=12),
         legend.title=element_text(size=12),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=12, angle=0, hjust=1),
         axis.text.y= element_text(size=12, angle=20, hjust=1),
         axis.title.x= element_text(size=16),
         axis.title.y= element_text(size=16),
         title=element_text(size=16),
         plot.margin = unit(c(1,1,1,1), "cm"))
parameter.plot.year

ggsave(parameter.plot.year,
       file=paste0(home.dic,"figures/parameter.plot.year.pdf"))
# probability

library(dplyr)
head(parameter.df)
parameter.prob.df <- parameter.df %>%
  mutate(prob.95 = case_when((set.value<Q97.5 & 
                                set.value>Q2.5)~ 1,
                             T ~0)) 
parameter.prob.df %>%
  aggregate(prob.95 ~ parameter , function(x) round((sum(x)/900)*100,digits=2)) %>%
  mutate(numb_year="all levels",
         numb_plot ="all levels") %>%
  bind_rows(parameter.prob.df %>%
              mutate(numb_plot = as.character(numb_plot)) %>%
              aggregate(prob.95 ~ parameter + numb_plot, function(x) round((sum(x)/300)*100,digits=2)) %>%
              mutate(numb_year="all levels")) %>%
  bind_rows(parameter.prob.df %>%
              mutate(numb_year= as.character(numb_year)) %>%
              aggregate(prob.95 ~ parameter + numb_year, function(x) round((sum(x)/300)*100,digits=2)) %>%
              mutate(numb_plot="all levels")) %>%
  spread(parameter,prob.95) %>%
  dplyr::select(numb_year,numb_plot,s,lambda,beta,g_init,beta_g)

write.csv(parameter.prob.df,
          file=paste0(home.dic,"results/parameter.prob.df.csv"))
head(parameter.prob.df)
view(parameter.prob.df)
