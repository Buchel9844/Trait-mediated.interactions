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
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- ""
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.Climate related interactions---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1 Along PDSI----
Realised.Int.PDSI.plotlist <- list()
country  ="aus"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  
  year.levels <-levels(as.factor(Realised.Int.Obs.list[[country]]$year))
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
  ggplot(env_pdsi) + geom_point(aes(y=PDSI.mean, x=year)) + theme_bw()
  
  Realised.Int.sum <- Realised.Int.Obs.list[[country]] %>%
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
    
  #view(  Median_weighted_df)
  Realised.Int.PDSI.plotlist[[paste0(country,"_sigmoid")]] <- Realised.Int.Obs.list[[country]] %>%
    left_join(env_pdsi%>%mutate(year= as.numeric(year))) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(y =sigmoid,x =Precip.extrem)) +
    geom_point(aes(fill=sigmoid),position="jitter",
               shape=21,
               alpha=0.2,
               size=2,
               stroke = 0) +
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous")),
                         limits=c(-max(Realised.Int.Obs.list[[country]]$sigmoid),max(Realised.Int.Obs.list[[country]]$sigmoid))) +
   geom_pointrange(data= Realised.Int.Obs.list[[country]] %>%
                    left_join(env_pdsi%>%mutate(year= as.numeric(year))) %>%
                     dplyr::filter(!sigmoid == 0) %>%
                    mutate(PosNeg = case_when(sigmoid<0 ~ "Comp",
                                              sigmoid>0 ~ "Fac",
                                              T ~ "Neutral"))  %>%
                    group_by(Precip.extrem,year,PosNeg ) %>%
                    summarise(sig.median = median(sigmoid,na.rm=T),
                              sig.q1= median(sigmoid,na.rm=T) - mad(sigmoid,na.rm=T),
                              sig.q9= median(sigmoid,na.rm=T) + mad(sigmoid,na.rm=T)),
                   aes(y=sig.median,
                       ymin=sig.q1,
                       ymax=sig.q9,
                       color=PosNeg )) +
    geom_line(data= Realised.Int.Obs.list[[country]] %>%
                      left_join(env_pdsi%>%mutate(year= as.numeric(year))) %>%
                      dplyr::filter(!sigmoid == 0) %>%
                      mutate(PosNeg = case_when(sigmoid<0 ~ "Comp",
                                                sigmoid>0 ~ "Fac",
                                                T ~ "Neutral"))  %>%
                      group_by(Precip.extrem,year,PosNeg ) %>%
                      summarise(sig.median = median(sigmoid,na.rm=T)),
                    aes(y=sig.median,
                        color=PosNeg )) +
   scale_color_manual(values=rev(c(wes_palette("Zissou1", 
                                                2, 
                                                type = "continuous")))) +
    theme_bw() +
    labs(x="Precipitation extrem (positive= wet)",
         y="Interaction strength",
         color="Median \ninteraction",
         fill="Strength of \ninteractions") + 
    theme(plot.margin = unit(c(1,0,0,0),"cm"),
          axis.title =element_text( size = 16), 
          legend.position="bottom",
          legend.key.size = unit(1, 'cm'),
          legend.title =element_text(size=16),
          legend.text =element_text(size=16),
          axis.text =element_text(size = 14))
  
  
  Median_weighted_df <- Realised.Int.sum %>%
    #filter(!focal == neigh) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.50 ~ "Comp",
                              RE.ratio.sig<0.50 ~ "Fac",
                              T ~ "Neutral")) %>%
    dplyr::filter(!RE.Q5.sig == 0) %>%
    group_by(Precip.extrem,year,PosNeg ) %>%
    mutate(RE.ratio.sig = abs(RE.ratio.sig -0.5)) %>%
    summarise(weighted.median = matrixStats::weightedMedian(RE.Q5.sig,RE.ratio.sig,na.rm=T),
              weighted.mad = matrixStats::weightedMad(RE.Q5.sig,RE.ratio.sig,na.rm=T)) %>%
    mutate(weightedQ1= weighted.median-weighted.mad,
           weightedQ9= weighted.median+weighted.mad) 
  
  Realised.Int.PDSI.plotlist[[country]] <- Realised.Int.sum %>%
    #filter(Precip.extrem > -200) %>%
    #filter(!focal == neigh) %>%
    #mutate(Precip.extrem.scaled = scale(Precip.extrem)) %>%
    dplyr::filter(!RE.Q5.sig == 0) %>%
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
    geom_pointrange(data=  Median_weighted_df,
                    aes(y=weighted.median,
                        x =Precip.extrem,
                        ymin=weightedQ1,
                        ymax=weightedQ9,
                        color=PosNeg),size=1) +
    geom_line(data=  Median_weighted_df,
                    aes(y=weighted.median,x =Precip.extrem,
                      color=PosNeg)) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous")) +
    scale_color_manual(values=rev(c(wes_palette("Zissou1", 
                                                2, 
                                                type = "continuous")))) +
    theme_bw() +
    labs(x="Precipitation extrem (positive= wet)",
         y="Median interaction strength",
         color="Weighted median \ninteraction",
         fill="Ratio of positive \nto competitive \ninteractions") + 
    theme(plot.margin = unit(c(1,0,0,0),"cm"),
          axis.title =element_text( size = 16), 
          legend.position="bottom",
          legend.key.size = unit(1, 'cm'),
          legend.title =element_text(size=16),
          legend.text =element_text(size=16),
          axis.text =element_text(size = 14))
  
  Realised.Int.PDSI.plotlist[[country]]
  
}
ggarrange(Realised.Int.PDSI.plotlist[["spain"]],
          Realised.Int.PDSI.plotlist[["aus"]],common.legend = T,
          nrow=2, legend = "bottom",labels=c("a.Spain","b.Australia"),
          font.label = list(size = 20, color = "black", face = "bold", family = NULL),
          label.x = 0,
          label.y = 1.0)
# Interactions.precipitation.pdf

ggarrange(Realised.Int.PDSI.plotlist[["spain_sigmoid"]],
          Realised.Int.PDSI.plotlist[["aus_sigmoid"]],common.legend = T,
          nrow=2, legend = "bottom",labels=c("a.Spain","b.Australia"),
          font.label = list(size = 20, color = "black", face = "bold", family = NULL),
          label.x = 0,
          label.y = 1.0)
#Interactions.precipitation.pdf

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.Looking at Interactions across time for answers----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1. Make data df ----
country = "spain"
Cool.theory.trait.year.df <- list()
density.quantile.name <- c("intercept","low","medium","high")
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
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
    mutate(Precip.extrem.scaled = scale(Precip.extrem)) %>%
    rename("year"="year.thermique") %>%
    mutate(PDSI.max = max(PDSI.mean)) %>%
    mutate_at("PDSI.mean",function(x) (x/ max(abs(x)))) %>% # scaling 
    mutate(year=as.numeric(year)) %>%
    as.data.frame()
  
  # make trait data frame
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] %>%
    dplyr::select(-c("Mean fecundity"))
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=trait.dist,
             scaled.trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list))) %>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>% 
                  rownames_to_column("focal") %>% rename("receiver.trait"=i) %>%
                  mutate(receiver.trait = scale(receiver.trait)))%>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>%
                  rownames_to_column("neigh") %>% rename("emitter.trait"=i)%>%
                  mutate(emitter.trait = scale(emitter.trait)))
    
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
    
    
  }
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  
  # Make data frame with trait and INTRA specific interactions
  trait.value.intra.year.df <- Realised.Obs.Year.list[[country]] %>%
    dplyr::filter(neigh ==focal) %>% # only INTRA
    left_join(specific.trait.dist,
              relationship ="many-to-many")
  Intra.trait.year.df <- NULL
  #trait.i = "SLA"
  for( trait.i in names(trait.df)){
    trait.intra.year.df.i <-  trait.value.intra.year.df  %>%
      dplyr::filter(trait==trait.i) %>%
      dplyr::select(year,neigh,trait,sigmoid,
                    emitter.trait,receiver.trait) %>%
      mutate(year=as.numeric(year)) %>%
      left_join(env_pdsi %>% dplyr::select(year,Precip.extrem.scaled))
    
    glm.intra.trait.year.i <- glmmTMB(sigmoid ~ 1 + receiver.trait + Precip.extrem.scaled +
                                        receiver.trait:Precip.extrem.scaled,
                                      trait.intra.year.df.i ,
                                      family="gaussian")
    #summary(glm.intra.trait.year.i)
    
    Intra.trait.year.df.i <- as.data.frame(confint( glm.intra.trait.year.i )) %>%
      rownames_to_column("parameters") %>%
      rename("Q2.5"="2.5 %",
             "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    Intra.trait.year.df <- bind_rows(Intra.trait.year.df,Intra.trait.year.df.i)
    
  }
  # Make data frame with trait and inter specific interactions
  trait.dist.year.df <- Realised.Obs.Year.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% # only INTRA
    left_join(specific.trait.dist,
              relationship ="many-to-many")
  Inter.trait.year.df <- NULL
  for( trait.i in names(trait.df)){
    trait.dist.year.df.i <-  trait.dist.year.df %>%
      dplyr::filter(trait==trait.i) %>%
      dplyr::select(neigh,focal,year,trait,sigmoid,emitter.trait,
                    receiver.trait,scaled.trait.dist)%>%
      mutate(year=as.numeric(year)) %>%
      left_join(env_pdsi %>% dplyr::select(year,Precip.extrem.scaled))
    
    glm.inter.trait.year.i <- glmmTMB(sigmoid ~ Precip.extrem.scaled +scaled.trait.dist+
                                        scaled.trait.dist*Precip.extrem.scaled +
                                        (1|focal) + (1|neigh),
                                      trait.dist.year.df.i,
                                      family="gaussian")
    
    Inter.trait.precip.df.i <- as.data.frame(confint( glm.inter.trait.year.i )) %>%
      rownames_to_column("parameters") %>%
      rename("Q2.5"="2.5 %",
             "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    Inter.trait.year.df <- bind_rows( Inter.trait.year.df, Inter.trait.precip.df.i)
  }
  
  Cool.theory.trait.year.df[[country]] <- list(
    trait.dist.year.df=trait.dist.year.df,
    Inter.trait.year.df=Inter.trait.year.df,
    Intra.trait.year.df =Intra.trait.year.df,
    env_pdsi=env_pdsi)
}
Cool.theory.trait.year.df[[country]]$Inter.trait.year.df
#---- 2.2. More plrect plot ----
Rect.glm.theory.trait.year.plotlist <- list()
country="aus" 
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
  if(country=="aus"){
    Inter.trait.year.df <- Cool.theory.trait.year.df[[country]]$Inter.trait.year.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA")))
    trait.levels <- c("SRL","Root tips","Root mass density","Root length",
                      "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                      "Canopy shape","Stem height","SLA")
  }
  if(country=="spain"){
    Inter.trait.year.df<- Cool.theory.trait.year.df[[country]]$Inter.trait.year.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA")))
    trait.levels <- c("SRL","Root diameter","Root mass density","SRA",
                      "Mean fecundity","C13 water use efficiency",
                      "Leaf C to N ratio","Leaf area index",
                      "Canopy shape","Stem height","SLA")
  }
  dummy.col <- c("SRL"="#4E79A7FF","SRA"="#76B7B2FF" ,"Root length"="#A4BED5FF","Root tips"="#512DA8FF" ,
                 "Root diameter"="#F28E2BFF" , "Root mass density"="#ED645AFF",
                 "Flower width"= "#FF9DA7FF" ,"Seed mass"="#B276B2FF",
                 "C13 water use efficiency"="#9C755FFF",
                 "Leaf C to N ratio"= "#BCBD22FF" ,
                 "Leaf area index"="#D4E157FF" ,"Canopy shape"="#72874EFF",
                 "SLA"="#59A14FFF","Stem height"="#FED789FF",
                 "intercept"="black")
  
  Inter.trait.year.df.i  <- Inter.trait.year.df %>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
    dplyr::filter(!signif  ==""|parameters=="(Intercept)") %>%
    mutate(parameters = case_when(parameters=="Precip.extrem.scaled" ~ "Extreme Precipitation regime\non INTERspecific interactions",
                                  parameters=="scaled.trait.dist" ~ "Distance trait",
                                  parameters=="Precip.extrem.scaled:scaled.trait.dist"~"Distancetrait.Precip",
                                  parameters=="(Intercept)"~"intercept",
                                  T~parameters) ,#"Dissimilarity between interactors on interspecific interaction according to precipitation regime"),
           parameter.impacted ="INTER interactions")
  head(  Inter.trait.year.df.i)
  Intra.trait.year.df.i  <- Cool.theory.trait.year.df[[country]]$Intra.trait.year.df%>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
    #mutate(trait=addline_format(trait)) %>%
    dplyr::filter(!signif  =="") %>%
    mutate(parameters = case_when(parameters=="Precip.extrem.scaled" ~ "Extreme Precipitation regime\n on INTRAspecific interactions",
                                  parameters=="receiver.trait" ~ "Receiver trait ",
                                  parameters=="(Intercept)"~"intercept",
                                  parameters=="receiver.trait:Precip.extrem.scaled" ~ "Receiver trait x Precip on INTRA"),
           parameter.impacted ="INTRA interactions")
  head(  Intra.trait.year.df.i)
  
  
  Cool.theory.trait.year.df[[country]]$trait.dist.year.df %>%
    left_join(Cool.theory.trait.year.df[[country]]$env_pdsi %>% dplyr::select(year,Precip.extrem.scaled)) %>%
    filter(trait=="SRA") %>%
    as.data.frame() %>%
    mutate(scaled.trait.dist=as.vector(scaled.trait.dist),
           Precip.extrem.scaled=as.vector(Precip.extrem.scaled)) %>%
    ggplot(aes(y=sigmoid,
               x=scaled.trait.dist,
               color=Precip.extrem.scaled,
               group=Precip.extrem.scaled)) +
    geom_jitter() +
    geom_smooth(method="lm")  
  
  Cool.theory.trait.year.df[[country]]$trait.dist.year.df %>%
    left_join(Cool.theory.trait.year.df[[country]]$env_pdsi %>% dplyr::select(year,Precip.extrem.scaled)) %>%
    filter(trait=="SLA") %>%
    as.data.frame() %>%
    mutate(scaled.trait.dist=as.vector(scaled.trait.dist),
           Precip.extrem.scaled=as.vector(Precip.extrem.scaled)) %>%
    ggplot(aes(y=as.factor(Precip.extrem.scaled),
               x=as.factor(scaled.trait.dist),
               color=sigmoid,
               fill=sigmoid)) + 
    #geom_point(shape=15, size=5) +
    geom_tile()+
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                   3, 
                                                   type = "continuous")),
                         limits=c(-0.3,0.3))+
    scale_color_gradientn(colours = rev(wes_palette("Zissou1", 
                                                    3, 
                                                    type = "continuous")),
                          limits=c(-0.3,0.3))+
    theme_clean()
  
  Inter.trait.year.plot.i <- Inter.trait.year.df.i %>%
    ggplot(aes(y=parameters,
               x=Estimate,
               color=trait)) +
    geom_pointrange(aes(xmin=  Q2.5,
                        xmax=Q97.5 ),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.3)) +
    scale_color_manual(values=dummy.col) +
    labs(x="coefficient",y="",
         color="") +
    guides(color = guide_legend(title.position = "top",
                                nrow=2))+ 
    geom_vline(xintercept=0,color="black") +
    theme_bw() +
    theme(legend.position="bottom",
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=16),
          axis.text = element_text(size=14),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
  Inter.trait.year.plot.i
  Intra.trait.year.plot.i<- Intra.trait.year.df.i%>%
    ggplot(aes(y=parameters,
               x=Estimate,
               color=trait)) +
    geom_pointrange(aes(xmin=  Q2.5,
                        xmax=Q97.5 ),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.3)) +
    scale_color_manual(values=dummy.col) +
    labs(x="coefficient",y="",
         color="") +
    guides(color = guide_legend(title.position = "top",
                                nrow=2))+ 
    geom_vline(xintercept=0,color="black") +
    theme_bw() +
    theme(legend.position="bottom",
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=16),
          axis.text = element_text(size=14),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
  Intra.trait.year.plot.i
  Rect.glm.theory.trait.year.plotlist[[country]] <- ggarrange(Inter.trait.year.plot.i,
                                                              Intra.trait.year.plot.i,
                                                              nrow=2,legend="bottom",label.x = 0.2,
                                                              align=c("v"),heights=c(1,1),
                                                              common.legend = T, 
                                                              labels=c("Inter","Intra"),
                                                              font.label = list(size = 20, color = "black", 
                                                                                face = "bold", family = NULL))
  Rect.glm.theory.trait.year.plotlist[[country]]
}
Rect.glm.theory.trait.year.plotlist[[paste0("spain")]]
Rect.glm.theory.trait.year.plotlist[[paste0("aus")]]
plot_grid(ggarrange(Rect.glm.theory.trait.year.plotlist[[paste0("spain")]],
                    Rect.glm.theory.trait.year.plotlist[[paste0("aus")]],
                    labels=c("a.Spain","b. Australia"),
                    font.label = list(size = 20, color = "black", 
                                      face = "bold", family = NULL),
                    legend="bottom",label.x = 0.2,
                    align=c("v"),heights=c(1,1),
                    common.legend = T,
                    ncol=2),
          ggpubr::get_legend(legend.plot),
          ncol = 1,
          rel_heights =c(1,0.2),
          nrow=2, labels=c("",""))
#figures/GLM.precip.trait.pdf
# 3D graphs
trait.dist.year.df.i <- Cool.theory.trait.year.df[[country]]$trait.dist.year.df %>%
  left_join(Cool.theory.trait.year.df[[country]]$env_pdsi %>% dplyr::select(year,Precip.extrem.scaled)) %>%
  filter(trait=="SLA") %>%
  as.data.frame() %>%
  mutate(z=sigmoid,
         x=as.vector(scaled.trait.dist),
         y=as.vector(Precip.extrem.scaled))
head(trait.dist.year.df.i)
str(trait.dist.year.df.i)
library(plotly)
summary(wet.df$y)
wet.df <- trait.dist.year.df.i[which(trait.dist.year.df.i$Precip.extrem.scaled >1),]
str(wet.df)
plot_ly(z = ~xtabs(z ~ x + y, 
                   data = wet.df )) %>% add_surface()
plot_ly() %>% 
  add_trace(data = wet.df,  
            x=wet.df$x, y=wet.df$y, z=wet.df$z, type="mesh3d" )

plot_ly(x=  as.vector(trait.dist.year.df.i$scaled.trait.dist), 
        y=  as.vector(trait.dist.year.df.i$Precip.extrem.scaled), 
        z=  trait.dist.year.df.i$sigmoid)%>% add_surface()