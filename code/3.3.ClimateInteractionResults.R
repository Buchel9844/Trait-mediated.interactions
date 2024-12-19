#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 0. SET UP: Import packages----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#packageurl <- "http://cran.r-project.org/src/contrib/Archive/cli/cli_3.6.3.tar.gz"
#install.packages(packageurl, repos=NULL, type="source",depedency=T)

library(cli)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages("HDInterval")
library("HDInterval")
#install.packages("tidyverse")
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
library(cowplot)
library(pals)
library(igraph)
library(statnet)
library(intergraph)
library(ggraph)
#setwd("/home/lbuche/Eco_Bayesian/chapt3")
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.Climate related interactions---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1 Along PDSI----
Realised.Int.PDSI.plotlist <- list()
country  ="aus"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  
  year.levels <-levels(as.factor(Realised.Int.Year.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  env_pdsi <- read.csv(paste0(home.dic,"results/",country,"_env_pdsi.csv")) 
  year.levels.all <-levels(as.factor(env_pdsi$year))
  
  env_pdsi <- read.csv(paste0(home.dic,"results/",country,"_env_pdsi.csv")) %>%
    mutate(year.thermique= as.numeric(factor(year))) %>% 
    mutate(year.thermique.2 = case_when(month > 9 ~ year.thermique + 1,
                                        T ~ year.thermique)) %>%
    group_by(year.thermique.2) %>% 
    mutate(PDSI.mean = mean(get(paste0(country,"_pdsi")), na.rm = T),
           PDSI.sd = sd(get(paste0(country,"_pdsi")),na.rm = T),
           Precip = sum(prec)) %>%
    ungroup() %>%
    mutate(year.thermique =year.levels.all[year.thermique.2]) %>%
    dplyr::select(year.thermique,PDSI.mean,PDSI.sd,Precip) %>%
    unique() %>%
    filter(year.thermique %in%  year.levels) %>%
    mutate(Precip.extrem = Precip - median(Precip)) %>%
    rename("year"="year.thermique") %>%
    mutate(PDSI.max = max(PDSI.mean)) %>%
    mutate_at("PDSI.mean",function(x) (x/ max(abs(x)))) %>% # scaling 
    as.data.frame()
  
  Realised.Int.sum <- Realised.Int.Year.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal,year) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.std= length(realised.effect.std[realised.effect.std<0]),
              RE.total.std = length(realised.effect.std),# high mean high competition
              RE.mean.std = mean(realised.effect.std, na.rm = T),
              RE.sd.std = sd(realised.effect.std,na.rm = T),
              RE.Q1.std = quantile(realised.effect.std, c(0.1), na.rm = T),
              RE.Q5.std = quantile(realised.effect.std, c(0.5), na.rm = T),
              RE.Q9.std = quantile(realised.effect.std,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T),
              RE.neg.sig.std= length(sigmoid[sigmoid.std<0]),
              RE.total.sig.std = length(sigmoid.std),# high mean high competition
              RE.mean.sig.std = mean(sigmoid.std, na.rm = T),
              RE.sd.sig.std = sd(sigmoid.std,na.rm = T),
              RE.Q1.sig.std = quantile(sigmoid.std, c(0.1), na.rm = T),
              RE.Q5.sig.std = quantile(sigmoid.std, c(0.5), na.rm = T),
              RE.Q9.sig.std = quantile(sigmoid.std,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.std = RE.neg.std/RE.total.std,
           RE.ratio.sig = RE.neg.sig/RE.total.sig,
           RE.ratio.sig.std = RE.neg.sig.std/RE.total.sig.std) %>%
    as.data.frame() %>%
    left_join(env_pdsi%>%mutate(year= as.numeric(year)))
  head(Realised.Int.sum)
  
  Realised.Int.PDSI.plotlist[[country]] <- Realised.Int.sum %>%
    #filter(Precip.extrem > -200) %>%
    filter(!focal == neigh) %>%
    #mutate(Precip.extrem.scaled = scale(Precip.extrem)) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.50 ~ "Comp",
                              RE.ratio.sig<0.50 ~ "Fac",
                              T ~ "Neutral")) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(y =RE.Q5.sig,x =Precip.extrem)) +
    geom_point(aes(fill=RE.ratio.sig),position="jitter",
               shape=21,
               alpha=0.2,
               size=2,
               stroke = 2) +
    geom_smooth(aes( group=as.factor(PosNeg),color=as.factor(PosNeg)),
                method="loess",se = T,size=2) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous")) +
    scale_color_manual(values=rev(c(wes_palette("Zissou1", 
                                                2, 
                                                type = "continuous")))) +
    #coord_cartesian(ylim=c(-0.5,0.5),clip="on") +
    #coord_cartesian(ylim=c(0.5,1.5),clip="on") +
    theme_bw()
  
  Realised.Int.PDSI.plotlist[[country]]
  
}
Realised.Int.PDSI.plotlist[["aus"]]
Realised.Int.PDSI.plotlist[["spain"]]