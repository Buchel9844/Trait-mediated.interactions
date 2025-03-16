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



#---- 2.Illustration -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ModelCheck.df <- read.csv(paste0(home.dic,"results/ModelCheck.df.csv"))


ggsave(
  plot=ggarrange(ggplot(Model.check, aes(sd.fec,mean.rmds,color=focal)) + 
                   geom_jitter(size=3,width = 3, height = 3,
                               aes(shape=as.factor(year)))+
                   geom_abline()+
                   labs(y="Mean RMDS across chains",
                        x="Standard deviation of observed fecundity",
                        shape="year") +
                   theme_bw(), 
                 ggplot(Model.check, aes(elpd_loo,mean.rmds,color=focal)) +
                   scale_x_reverse() +
                   geom_jitter(size=3,width = 3, height = 3,
                               aes(shape=as.factor(year)))+
                   labs(y="Mean RMDS across chains",
                        x="ELPD_loo",
                        shape="year") +
                   theme_bw(),
                 nrow=1,
                 labels = c("a.","b."),
                 common.legend = T,
                 legend="bottom"
  ),
  filename=paste0(home.dic,"figure/ModelBehaviorRMDS.linear.pdf"))


df_param_all <- read.csv(paste0(home.dic,"results/Chapt1_Parameters_values.csv"))

Trophic.mean.df <- read.csv(paste0(home.dic,"results/Trophic.mean.df.csv"))

abund.mean.df <- read.csv(paste0(home.dic,"results/abund.mean.df.csv"))



ModelCheck <- read.csv(paste0(home.dic,"results/ModelCheck.csv")) %>%
  rename("complexity.plant"=complexity)
y.lim <- data.frame(focal=rep(c("CETE","CHFU","HOMA","LEMA"),each=3),
                    year=rep(c(2019,2020,2021),times=4),
                    upp = c(1500,1500,1500,
                            1000,300,750,
                            100,700,700,
                            1000,300,500)
)
library(ggridges)
range.plot.rmds <- param.fec.pred %>%
  inner_join(y.lim) %>%
  filter(exp(effect.total.max) < upp &
           exp(effect.total.min) < upp &
           exp(effect.total.mean) < upp) %>%
  ggplot() +
  geom_density_ridges(aes(y=complexity.plant,x=exp(effect.total.min)),
                      alpha=0.2,fill="blue",scale=1) +
  geom_density_ridges(aes(y=complexity.plant,x=exp(effect.total.max)),
                      alpha=0.2,fill="red",scale=1) +
  geom_density_ridges(aes(y=complexity.plant,x=exp(effect.total.mean)),
                      alpha=0.2,fill="orange",scale=1) +
  #stat_summary(aes(y=complexity.plant,x=exp(effect.total.mean)),
  #            color="orange",
  #           fun = "median", geom = "point", size = 2) +
  #stat_summary(aes(y=complexity.plant,x=exp(effect.total.max)),
  #            fun = "median", geom = "point", 
  #           size = 2,color="red",alpha=0.5) +
  #stat_summary(aes(y=complexity.plant,x=exp(effect.total.min)),
  #         fun = "median", geom = "point",
  #         size = 2,color="blue",alpha=0.5) +
  
  geom_errorbarh(data=ModelCheck,height = .2,
                 aes(y=complexity.plant,
                     xmin=min.rmds,
                     xmax=max.rmds),
                 color="black",size=1,alpha=0.6,
                 position=position_dodge(0.3)) +
  geom_point(data=ModelCheck,size=2,shape=19,color="green",
             aes(y=complexity.plant,
                 x=mean.rmds),
             color="black",
             position=position_dodge(0.3)) +
  theme_bw() +
  labs(y="grouping factor",
       x="Predicted fecundity")+
  #scale_x_continuous(aes(limits=c(0,upp))) +
  facet_wrap(focal ~ year,ncol=3,
             scale="free") +
  guides(shape=guide_legend("RMDS",
                            override.aes = list(size = 10,color="black"),
                            nrow = 1,
                            direction="horizontal",
                            byrow = T,
                            title.hjust = 0.1))+
  theme(strip.placement = "outside",
        legend.key.size = unit(1, 'cm'),
        title =element_text(size=12),
        axis.text.x= element_text(size=16),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        strip.text = element_text(size=12),
        panel.border = element_rect(color = "black", 
                                    fill = NA))
range.plot.rmds

ggsave(paste0(home.dic,"figure/range.plot.rmds.pdf"),
       #dpi="retina",
       width = 30,
       height = 20,
       units = c("cm"),
       range.plot.rmds
)

