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

plant_code_spain <- read.csv( "data/plant_code_spain.csv",
                              header = T, stringsAsFactors = F, sep=",",
                              na.strings = c("","NA"))


#---- 2.1 Alpha----
year.int = "All"    

source(paste0(home.dic,"code/PopProjection_toolbox.R"))

list_alpha_mean <- list()
list_alpha_sd <- list()

#for(Code.focal in species.spain){
  load(file= paste0(home.dic,
                    "results/Parameters_",Code.focal,"_",year.int,".Rdata"))
  
  df_alpha_generic_param = parameter$df_alpha_generic_param

  list_alpha_mean[[Code.focal]] <- df_alpha_generic_param %>%
      group_by(parameter) %>% 
      summarise(across(all_of(species.spain),mean)) %>%
      bind_rows(as.data.frame(colMeans(parameter$df_N_opt))%>%
                  rownames_to_column("neigh") %>%
                  spread(neigh,
                         "colMeans(parameter$df_N_opt)") %>%
                  mutate(parameter="N_opt")) %>%
    column_to_rownames("parameter")
  
  list_alpha_sd[[Code.focal]] <- df_alpha_generic_param %>%
    group_by(parameter) %>% 
    summarise(across(all_of(species.spain),sd)) %>%
    bind_rows(as.data.frame(sapply(parameter$df_N_opt,sd))%>%
                rownames_to_column("neigh") %>%
                spread(neigh,
                       "sapply(parameter$df_N_opt, sd)") %>%
                mutate(parameter="N_opt"))%>%
    column_to_rownames("parameter")

 #   }

#---- 2.2 Germination ----
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

#---- 2.3 Abundance over time ----
abundance_spain <- read.csv(paste0("data/abundance_spain.csv"),
                            header = T,stringsAsFactors = F, sep=",",
                            na.strings=c("","NA"))
plant_code_spain <- read.csv(paste0( "data/plant_code_spain.csv"),
                             header = T, stringsAsFactors = F, sep=",",
                             na.strings = c("","NA"))

abundance_spain_short <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  aggregate(individuals ~code.analysis + year + plot + subplot , mean) %>%
  mutate(com_id = paste(plot,subplot,sep="_")) %>%
  spread(code.analysis,individuals)
str(abundance_spain_short)

abundance_spain_short  %>% 
  ggplot(aes(y=individuals ,x=as.factor(year),
             color=as.factor(code.analysis))) +
  geom_boxplot()


abundance_spain_rare <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  filter(code.analysis=="rare") %>% 
  aggregate(individuals ~ subplot + year + code.plant, sum) %>% 
  aggregate(individuals ~ subplot + year , mean) %>% 
  ggplot(aes(y=individuals ,x=as.factor(year))) +
  geom_boxplot()

#---- 3. Extraction of realised intrinsic growth rate ----

focal = "LEMA"
year.levels = levels(as.factor(abundance_spain_short$year))
year.int = nlevels(as.factor(abundance_spain_short$year))
growth.ratio.df <- NULL
for(y in 1:nlevels(as.factor(abundance_spain_short$year))){
  for(i in levels(as.factor(abundance_spain_short$com_id))){
 if(y ==7) next
    growth.ratio.df <- bind_rows(growth.ratio.df,
    data.frame(focal = focal,
               year = as.integer(year.levels[y]),
               com_id=i,
               growth.ratio= abundance_spain_short[which(abundance_spain_short$year== as.integer(year.levels[y+1]) &
                                                   abundance_spain_short$com_id ==i),"LEMA"] / 
                 abundance_spain_short[which(abundance_spain_short$year== as.integer(year.levels[y]) &
                                    abundance_spain_short$com_id ==i),"LEMA"]))
  }
}

abundance_spain_LEMA <- abundance_spain_short %>%
  select(-c("plot","subplot")) %>% 
  left_join(growth.ratio.df,by=c("year", 'com_id')) %>%
  filter(!is.na(growth.ratio) & !is.nan(growth.ratio) & !is.infinite(growth.ratio))


view(abundance_spain_LEMA)
DataVec <- list(N= nrow(abundance_spain_LEMA),
                S= length(species.spain),
                SpAbundance = abundance_spain_LEMA %>%
                select(all_of(species.spain)),
                growth_ratio = abundance_spain_LEMA$growth.ratio ,
                lambda_mean_obs = mean(abundance_spain_LEMA$growth.ratio),
                year =  as.integer(factor(abundance_spain_LEMA$year,unique(abundance_spain_LEMA$year))),
                Y= nlevels(as.factor(abundance_spain_LEMA$year)),
                g_mean = seed_germination_spain$g.mean[which(seed_germination_spain$code.analysis == focal)],
                g_sd = seed_germination_spain$g.sd[which(seed_germination_spain$code.analysis == focal)],
                seed_s_mean = seed_survival_spain$s.mean[which(seed_survival_spain$code.analysis == focal)],
                seed_s_sd = seed_survival_spain$s.sd[which(seed_survival_spain$code.analysis == focal)],
                alpha_initial_mean = list_alpha_mean[[focal]]["alpha_initial",],
                alpha_slope_mean = list_alpha_mean[[focal]]["alpha_slope",],
                c_mean = list_alpha_mean[[focal]]["c",],
                N_opt_mean = list_alpha_mean[[focal]]["N_opt",],
                alpha_initial_sd = list_alpha_sd[[focal]]["alpha_initial",],
                alpha_slope_sd = list_alpha_sd[[focal]]["alpha_slope",],
                c_sd = list_alpha_sd[[focal]]["c",],
                N_opt_sd = list_alpha_sd[[focal]]["N_opt",]
                )
list.init <- function(...)list(N_opt= array(as.numeric(sapply(abundance_spain_LEMA %>%
                                 select(all_of(species.spain)),median),
                                 dim = DataVec$S)),
                               lambda_mean= array(as.numeric( mean(abundance_spain_LEMA$growth.ratio),
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
                 chains = 1,
                 seed= 1616) 

save(file= paste0(home.dic,"results/Inv.Modelfit_",Code.focal,".rds"),
     Inv.Modelfit)

#load(paste0(home.dic,"results/Inv.Modelfit_",Code.focal,"_",year.int,".rds"))

Inv.ModelfitPosteriors <- rstan::extract(Inv.Modelfit)

save(file= paste0(home.dic,"results/Inv.ModelfitPosteriors",Code.focal,".Rdata"),
     Inv.ModelfitPosteriors)
load(file= paste0(home.dic,"results/Inv.ModelfitPosteriors",Code.focal,".Rdata"))


Inv.Modelfit_loo <- rstan::loo(Inv.Modelfit,pars ="F_sim")

save(file= paste0(home.dic,"results/Inv.ModelfitLOO",Code.focal,".Rdata"),
     Inv.Modelfit_loo)

print("Final Fit done")

#---- 3.3. Final fit posterior check and behavior checks---- 

##### Diagnostic plots and post prediction 
pdf(paste0(home.dic,"figures/Inv.Modelfit_",Code.focal,".pdf"))
# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
#source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots

# check the distribution of Rhats and effective sample sizes 
##### Posterior check
stan_post_pred_check(Inv.ModelfitPosteriors,"GR",
                     abundance_spain_LEMA$growth.ratio) 
hist(abundance_spain_LEMA$growth.ratio,breaks = 150)
# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(Inv.Modelfit)$summary[,"Rhat"],
     main = paste("Inversed Fit: Histogram of Rhat for",
                  Code.focal," and ",year.int))
hist(summary(Inv.Modelfit)$summary[,"n_eff"],
     main = paste("Inversed Fit: Histogram of Neff for",
                  Code.focal," and ",year.int))

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
  setNames(levels(as.factor(abundance_spain_LEMA$year)))

write.csv(bind_cols(df_lambda_mean,df_lambda_sd),
          file=paste0(home.dic,"results/IntGR_",Code.focal,".csv"))
