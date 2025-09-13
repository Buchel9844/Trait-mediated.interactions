
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
#install.packages("ggpubr")
library(ggpubr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("lme4")
library(lme4)
#install.packages("vegan")
library(vegan)
#install.packages("wesanderson")
library(wesanderson) # for color palette
#install.packages("ggthemes")
library(ggthemes) 
#install.packages("grid")
library(grid)
#install.packages("lmtest")
library(lmtest)
library(cowplot)
library(PerformanceAnalytics)
library(vegan)
library(brms)
library(ggridges)
library(tidybayes)
library(ggdist)
library(bayestestR) # for bayesian stats like p_direction
#setwd("/home/lbuche/Eco_Bayesian/chapt3")

project.dic <- ""
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"

#remove.packages("TMB")
#install.packages("/Users/lisabuche/Downloads/TMB_1.9.15.tar.gz",
#                 type = 'source')
#install.packages("TMB", version='1.9.15')

#remotes::install_github("glmmTMB/glmmTMB/glmmTMB")
library(glmmTMB)

#---- 0.1. Import results----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))
country.list <- c("aus","spain")
# parameter from models
load(paste0(home.dic,"results/parameters_alpha.RData")) 
# Raw sigmoid
load(paste0(project.dic,"results/Param.sigm.df.RData"))
# realised interactions for each year
load(paste0(project.dic,"results/Theoretical.Int.list.RData")) 

# PDSI
env_pdsi_aus <- read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) 
env_pdsi_spain <- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv")) 
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.Looking at trait-mediated interactions----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#----1.0. Get phylo distance----
library(taxize)
library(ape)
plant_code_aus <- read.csv("data/aus_rawdata/plant_code_aus.csv",
                           header = T, stringsAsFactors = F, sep=",",
                           na.strings = c("","NA"))%>%
  mutate(code.analysis=final.code,
         species=name)%>%
  mutate(species = case_when(species=="Podolepis aristata"~ "Podolepis canescens",
                             T~species))
plant_code_spain <- read.csv( "data/spain_rawdata/plant_code_spain.csv",
                              header = T, stringsAsFactors = F, sep=",",
                              na.strings = c("","NA")) %>%
  filter(!species %in% c("Melilotus sulcatus","Polypogon monspeliensis"))
taxize_dist.list <- list()
country="spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  # make trait data frame
  species_name <-  get(paste0("plant_code_",country)) %>%
    dplyr::filter(code.analysis %in% Code.focal.list) %>%
    select(code.analysis,species) %>%
    unique()

  taxize_class <-classification( species_name$species,
                                 db="ncbi")
  taxize_tree <- class2tree( taxize_class, check = TRUE)
  #plot(taxize_tree)
  taxize_dist.df <-   taxize_tree$distmat %>%
    as.matrix() %>%
    as.data.frame()%>%
    gather(.,key="neigh.name",value="phylo.dist") %>%
    mutate(phylo.dist=phylo.dist/max(phylo.dist),
           focal.name=rep(species_name$species,
                          times=length(species_name$species))) %>%
    left_join(species_name %>%
                dplyr::select(species,code.analysis) %>%
                rename("focal.name"="species",
                       "focal"="code.analysis"))%>%
    left_join(species_name %>%
                dplyr::select(species,code.analysis) %>%
                rename("neigh.name"="species",
                       "neigh"="code.analysis"))
  view(taxize_dist.df)
  
  taxize_dist.list[[country]]<- taxize_dist.df
}

save(taxize_dist.list,
     file="results/taxize_dist.list.RData")
load(file="results/taxize_dist.list.RData")

#---- 1.1. Make data df ----
Cool.theory.trait.df <- list()
density.quantile.name <- c("intercept","low","medium","high")
rsq <- function (x, y) cor(x, y,use="complete.obs") ^ 2
zscore <- function (x) (x - mean(x,na.rm=T))/sd(x,na.rm=T)

Cool.theory.trait.df[["spain"]] <- list(
  trait.dist.df=NULL,
  trait.value.lambda.df = NULL, 
  trait.value.intra.df=NULL,
  Inter.trait.df=NULL,
  Lambda.trait.df=NULL,
  Intra.trait.df =NULL)

Cool.theory.trait.df[["aus"]] <- list(
  trait.dist.df=NULL,
  trait.value.lambda.df = NULL, 
  trait.value.intra.df=NULL,
  Inter.trait.df=NULL,
  Lambda.trait.df=NULL,
  Intra.trait.df =NULL)

#---- 1.1.1 Make trait distant data.frame ----

for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      bind_cols(outer(trait.df[,i], trait.df[,i], '/') %>%
                  as.data.frame() %>%
                  gather(.,key="focal",value="trait.ratio")) %>%
      mutate(trait.dist=trait.dist,
             scaled.trait.dist = as.vector(scale(trait.dist)),
             scaled.trait.ratio = as.vector(scale(trait.ratio)),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list))) %>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>% 
                  rownames_to_column("focal") %>% rename("receiver.trait"=i) %>%
                  mutate(receiver.trait.scaled = as.vector(zscore(receiver.trait)),
                         receiver.trait.log.scaled = as.vector(zscore(log(receiver.trait)))))%>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>%
                  rownames_to_column("neigh") %>% rename("emitter.trait"=i)%>%
                  mutate(emitter.trait.scaled = as.vector(zscore(emitter.trait)),
                         emitter.trait.log.scaled = as.vector(zscore(log(emitter.trait)))))
    
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
    
    
  }
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      dplyr::filter(trait %in% c("SRL","Root tips","Root mass density","Root length",
                      "C13 water use efficiency","Flower width","Seed mass",
                      "Stem height","SLA"))%>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "C13 water use efficiency","Flower width","Seed mass",
                                            "Stem height","SLA"))) 
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      dplyr::filter(trait %in% c("SRL","Root diameter","Root mass density","SRA",
                      "C13 water use efficiency","Leaf C to N ratio",#"Leaf area index",
                      "Flower width","Seed mass","Stem height","SLA"))%>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                        "C13 water use efficiency","Leaf C to N ratio",#"Leaf area index",
                                            "Flower width","Seed mass","Stem height","SLA"))) 
  }
  Cool.theory.trait.df[[country]]$trait.dist.df <- specific.trait.dist
}

#---- 1.1.2 Run lambda regression----

for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] %>%
    dplyr::select(-any_of(c("Canopy shape","Mean fecundity","Leaf area index")))
  specific.trait.dist <- Cool.theory.trait.df[[country]]$trait.dist.df
  # Make data frame with trait and Lambda
  trait.value.lambda.df <- Theoretical.Int.list[[country]] %>%
    dplyr::select(focal,lambda,density.quantile) %>% # only lambda
    unique() %>%
    left_join(specific.trait.dist %>%
                dplyr::select(focal,receiver.trait.scaled,
                              receiver.trait,trait) %>% unique(),
              relationship ="many-to-many")
  Lambda.trait.df <- NULL
  glm.lambda.trait.summary <- NULL
  #trait.i = "Flower width"
  #n ="low"
  for( trait.i in names(trait.df)){
    for(n in c("low","high")){
      trait.lambda.df.i <-   trait.value.lambda.df  %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(density.quantile==n) %>%
        mutate(focal=as.factor(focal),
               lambda=zscore(log(lambda)))
      if(i %in% c("C13 water use efficiency")){
        trait.lambda.df.i <- trait.lambda.df.i %>%
          dplyr::select(-receiver.trait.scaled) %>%
          dplyr::rename("receiver.trait.scaled" ="receiver.trait.log.scaled") %>%
          mutate(receiver.trait= log(receiver.trait))
      }
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores())
      glm.lambda.trait.outlier.i <- brm(lambda ~ 1 + receiver.trait.scaled + (1|focal) ,
                                    trait.lambda.df.i %>%
                                      dplyr::filter(receiver.trait.scaled>-3 &
                                                      receiver.trait.scaled<3)%>%
                                      mutate(receiver.trait.scaled = zscore(receiver.trait)),
                                    family="gaussian",
                                    chains = 4,
                                    iter = 2000,
                                    warmup = floor(2000/2))
      #summary( glm.lambda.trait.outlier.i)
     # plot( glm.lambda.trait.outlier.i)
      
      Lambda.trait.df.i <- as.data.frame(glm.lambda.trait.outlier.i,
                  variable =c("b_Intercept","b_receiver.trait.scaled")) %>%
        rename("Intercept"="b_Intercept",
               "Receiver.trait.scaled"="b_receiver.trait.scaled" ) %>%
        mutate(trait=trait.i,
               density.quantile=n,
               rhat= max(rhat(glm.lambda.trait.outlier.i)),
               num.div = rstan::get_num_divergent(glm.lambda.trait.outlier.i$fit))
      
      Lambda.trait.df <- bind_rows( Lambda.trait.df, Lambda.trait.df.i)
      
    }
  }
  Cool.theory.trait.df[[country]]$trait.value.lambda.df <- trait.value.lambda.df
  Cool.theory.trait.df[[country]]$Lambda.trait.df <- Lambda.trait.df

}
Lambda.trait.df <- Cool.theory.trait.df[[country]]$Lambda.trait.df
save(Cool.theory.trait.df,
     file="results/Cool.theory.trait.df.Rdata")
#---- 1.1.3 Run INTRA regression----

for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] %>%
    dplyr::select(-any_of(c("Canopy shape","Mean fecundity","Leaf area index")))
  specific.trait.dist <- Cool.theory.trait.df[[country]]$trait.dist.df
  # Make data frame with trait and INTRA specific interactions
  trait.value.intra.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(neigh==focal) %>%
    left_join(specific.trait.dist,
              relationship ="many-to-many",
              )
  Intra.trait.df <- NULL
  glm.intra.trait.summary <- NULL
  #trait.i = "SLA"
  #trait.i = "Flower width"
  for( trait.i in names(trait.df)){
    print(trait.i)
    for(n in c("low","high")){#density.quantile.name){
      trait.intra.df.i <-  trait.value.intra.df  %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(density.quantile==n) %>%
        dplyr::select(focal,density.quantile,trait,theoretical.effect,
                      emitter.trait,receiver.trait,
                      emitter.trait.scaled,receiver.trait.scaled,
                      emitter.trait.log.scaled,receiver.trait.log.scaled) %>%
        mutate(focal=as.factor(focal))
      if(i %in% c("C13 water use efficiency")){
        trait.intra.df.i <-  trait.intra.df.i %>%
          dplyr::select(-receiver.trait.scaled) %>%
          dplyr::rename("receiver.trait.scaled" ="receiver.trait.log.scaled") %>%
          mutate(receiver.trait= log(receiver.trait))
      }
      
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores())
      glm.intra.trait.outlier.i <- brm(theoretical.effect ~ 1 + receiver.trait.scaled  + (1|focal),
                                   trait.intra.df.i %>%
                                     dplyr::filter(receiver.trait.scaled>-3 &
                                                      receiver.trait.scaled<3)%>%
                                     mutate(receiver.trait.scaled = zscore(receiver.trait)),
                                   family="gaussian",
                                   chains = 4,
                                   iter = 2000,
                                   warmup = floor(2000/2))
     
      Intra.trait.df.i <- as.data.frame(glm.intra.trait.outlier.i,
                                        variable =c("b_Intercept","b_receiver.trait.scaled")) %>%
        rename("Intercept"="b_Intercept",
               "Receiver.trait.scaled"="b_receiver.trait.scaled" ) %>%
        mutate(trait=trait.i,
               density.quantile=n,
               rhat= max(rhat(glm.intra.trait.outlier.i)),
               num.div = rstan::get_num_divergent(glm.intra.trait.outlier.i$fit))
      
      
      Intra.trait.df <- bind_rows(Intra.trait.df,Intra.trait.df.i)
      
    }
  }
  Cool.theory.trait.df[[country]]$trait.value.intra.df <-  trait.value.intra.df
  Cool.theory.trait.df[[country]]$Intra.trait.df  <- Intra.trait.df 
}
save(Cool.theory.trait.df,
     file="results/Cool.theory.trait.df.Rdata")

#---- 1.1.3 Run INTER regression----
country="spain"
for( country in country.list){
    Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
    trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] %>%
      dplyr::select(-any_of(c("Canopy shape","Mean fecundity","Leaf area index")))
    specific.trait.dist <- Cool.theory.trait.df[[country]]$trait.dist.df
  # Make data frame with trait and inter specific interactions
    taxize_dist.df <- taxize_dist.list[[country]]
  
  trait.dist.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    left_join(as.data.frame(specific.trait.dist),
              relationship ="many-to-many") %>%
    left_join(  taxize_dist.df %>%
                  dplyr::select(phylo.dist,focal,neigh))
  Inter.trait.df <- NULL

 # Inter.trait.df <- Cool.theory.trait.df[[country]]$Inter.trait.df  
  
    # trait.i ="SRL"
  # n ="low"
    for( trait.i in names(trait.df)){ #names(trait.df)){
      for(n in c("low","high")){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(density.quantile==n) %>%
        dplyr::select(neigh,focal,density.quantile,trait,
                      theoretical.effect,
                      emitter.trait,emitter.trait.scaled,emitter.trait.log.scaled,
                      receiver.trait,receiver.trait.scaled,receiver.trait.log.scaled,
                      trait.dist,scaled.trait.dist,trait.ratio,scaled.trait.ratio,
                      phylo.dist) %>%
        mutate(focal=as.factor(focal),
               neigh=as.factor(neigh),
               phylo.dist=abs(phylo.dist))
      
      if(i %in% c("C13 water use efficiency")){
        trait.dist.df.i <-    trait.dist.df.i %>%
          dplyr::select(-c('receiver.trait.scaled','emitter.trait.scaled','scaled.trait.dist','trait.dist')) %>%
          dplyr::rename("receiver.trait.scaled" ="receiver.trait.log.scaled",
                        "emitter.trait.scaled" ="emitter.trait.log.scaled",
                        "scaled.trait.dist" ="scaled.trait.ratio",
                        "trait.dist"="trait.ratio") %>%
          mutate(receiver.trait= log(receiver.trait),
                 emitter.trait= log(emitter.trait))
      }
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores())
  
      glm.inter.trait.outlier.i <- brm(theoretical.effect ~  receiver.trait.scaled  + emitter.trait.scaled  +   (1|focal) + (1|neigh) + (1|phylo.dist), 
                                                trait.dist.df.i %>%
                                             dplyr::filter(receiver.trait.scaled > -3 &
                                                             receiver.trait.scaled < 3 &
                                                             emitter.trait.scaled >-3 &
                                                             emitter.trait.scaled < 3)%>%
                                             mutate(receiver.trait.scaled = zscore(receiver.trait),
                                                    emitter.trait.scaled = zscore(emitter.trait)),
                                          family="gaussian",
                                           chains = 4,
                                           iter = 2000,
                                           prior =  set_prior("normal(0,0.1)",
                                                                  class = "b", 
                                                                  lb = -1, ub=1),
                                           warmup = floor(2000/2))
    

      #summary(glm.inter.trait.outlier.i)
      #plot(glm.inter.trait.outlier.i)
      Inter.trait.df.i <- as.data.frame(glm.inter.trait.outlier.i,
                                        variable =c("b_Intercept",
                                                    "b_receiver.trait.scaled",
                                                    "b_emitter.trait.scaled")) %>%
        rename("Intercept"="b_Intercept",
               "Receiver.trait.scaled"="b_receiver.trait.scaled",
               "Emitter.trait.scaled" = "b_emitter.trait.scaled") %>%
        mutate(trait=trait.i,
               density.quantile=n,
               rhat= max(rhat(glm.inter.trait.outlier.i)),
               num.div = rstan::get_num_divergent(glm.inter.trait.outlier.i$fit))
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores())
      glm.inter.trait.dist.i <- brm(theoretical.effect ~  scaled.trait.dist + (1|focal) + (1|neigh) + (1|phylo.dist), 
                                       trait.dist.df.i,
                                       family="gaussian",
                                       chains = 4,
                                       iter = 2000,
                                       prior =  set_prior("normal(0,0.1)",
                                                          class = "b", 
                                                          lb = -1, ub=1),
                                       warmup = floor(2000/2))
      
      Inter.trait.df.dist.i <- as.data.frame(glm.inter.trait.dist.i,
                                        variable =c("b_Intercept",
                                                    "b_scaled.trait.dist")) %>%
        rename("Intercept.dist"="b_Intercept",
               "Scaled.trait.dist"="b_scaled.trait.dist") %>%
        mutate(rhat.dist= max(rhat(glm.inter.trait.dist.i)),
               num.div.dist = rstan::get_num_divergent(glm.inter.trait.dist.i$fit))
      
      Inter.trait.df <- bind_rows(Inter.trait.df, bind_cols(Inter.trait.df.i,
                                                            Inter.trait.df.dist.i))
      Cool.theory.trait.df[[country]]$Inter.trait.df  <- Inter.trait.df
    #  save(Cool.theory.trait.df,
       #    file="results/Cool.theory.trait.df.Rdata")
      
    }
  }
  
  Cool.theory.trait.df[[country]]$trait.dist.df  <-  trait.dist.df 
  Cool.theory.trait.df[[country]]$Inter.trait.df  <- Inter.trait.df
  save(Cool.theory.trait.df,
       file="results/Cool.theory.trait.df.Rdata")
}
save(Cool.theory.trait.df,
     file="results/Cool.theory.trait.df.Rdata")

load("results/Cool.theory.trait.df.Rdata")

view(Cool.theory.trait.df[[country]]$Inter.trait.df)
str(Cool.theory.trait.df[[country]])
ggplot()+
  geom_point(aes(x=trait.dist.df$receiver.trait[which(trait.dist.df$trait=="SLA")],
                 y=trait.dist.df$receiver.trait[which(trait.dist.df$trait=="C13 water use efficiency")]))

#---- 1.2. Make detailed graphs ----
Cool.detailed.theory.trait.plotlist <- list()
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
density.quantile.name <- c("low","high")
country="aus"
dummy.col <- c("SRL"="#4E79A7FF","SRA"="#76B7B2FF" ,"Root length"="#A4BED5FF","Root tips"="#512DA8FF" ,
               "Root diameter"="#F28E2BFF" , "Root mass density"="#ED645AFF",
               "Flower width"= "#FF9DA7FF" ,"Seed mass"="#B276B2FF",
               "C13 water use efficiency"="#9C755FFF",
               "Leaf C to N ratio"= "#BCBD22FF" ,
               "Leaf area index"="#D4E157FF" ,"Canopy shape"="#72874EFF",
               "SLA"="#59A14FFF","Stem height"="#FED789FF",
               "intercept"="grey80")

dummy.names <- c("SRL"="Specific root length",
                 "SRA"="Specific root area" ,
                 "Root length"= "Root length","Root tips"="Root tips",
                 "Root diameter"="Root diameter" , "Root mass density"="Root tissue density",
                 "Flower width"= "Flower width" ,"Seed mass"="Seed mass",
                 "C13 water use efficiency"="Water use efficiency",
                 "Leaf C to N ratio"= "Nitrogen use efficiency" ,
                 "Leaf area index"="Leaf area index" ,"Canopy shape"="Canopy shape",
                 "SLA"="Specific leaf area","Stem height"="Stem height")

country="spain"

n="low"
for( country in country.list){
  for(n in c("low","high")){
    
    Inter.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      gather(any_of(c("Scaled.trait.dist",
               "Emitter.trait.scaled","Receiver.trait.scaled")),
             key="trait.param",value="trait.coeff")  %>%
      dplyr::filter(!is.na(trait.coeff)) %>%
      mutate(trait= factor(trait ,levels=names(dummy.names) )) %>%
      mutate(trait.param=tolower(trait.param))
          

    Intra.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      rename("trait.coeff"="Receiver.trait.scaled") %>%
      mutate(trait= factor(trait ,levels=names(dummy.names) ))
    
    
    Lambda.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      rename("trait.coeff"="Receiver.trait.scaled") %>%
      mutate(trait= factor(trait ,levels=names(dummy.names)))
    
    Inter.trait.df.i <- Cool.theory.trait.df[[country]]$trait.dist.df %>%
      dplyr::select(theoretical.effect,receiver.trait.scaled,trait,emitter.trait.scaled,scaled.trait.dist) %>%
      gather(receiver.trait.scaled,emitter.trait.scaled,scaled.trait.dist,
             key="trait.param", value="trait.value") %>%
      rename("raw.value"="theoretical.effect") %>%
      mutate(parameter="INTER") #%>% 
      #left_join(Inter.trait.df.long.i %>%
      #            mutate(trait.param=tolower(trait.param))) %>%
      #dplyr::filter(!is.na(trait.coeff))
    
    Intra.trait.df.i <- Cool.theory.trait.df[[country]]$trait.value.intra.df %>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(theoretical.effect,receiver.trait.scaled,trait) %>%
      rename("raw.value"="theoretical.effect",
             "trait.value" ="receiver.trait.scaled")%>%
      left_join(Intra.trait.df.long.i) %>%
      mutate(parameter="INTRA",
             trait.param="receiver.trait.scaled") %>% 
      dplyr::filter(!is.na(trait.coeff)) 
    
    Lambda.trait.df.i <- Cool.theory.trait.df[[country]]$trait.value.lambda.df %>%
      dplyr::filter( trait %in% levels(as.factor(Lambda.trait.df.long.i$trait))) %>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(lambda,receiver.trait.scaled,trait) %>%
      rename("raw.value"="lambda",
             "trait.value" ="receiver.trait.scaled")%>%
      left_join(Lambda.trait.df.long.i) %>%
      mutate(parameter="lambda",
             trait.param="receiver.trait.scaled",
             raw.value=zscore(log(raw.value)))%>% 
      dplyr::filter(!is.na(trait.coeff))
   
    
    pdf(file = paste0("figures/GLM/GLM_details_",country,"_",n,"_density.pdf"),
        onefile = TRUE)
    
    Cool.detailed.theory.trait.plotlist[[paste0(country,"_Inter")]] <- ggplot() +
      geom_point(data=Inter.trait.df.i %>% dplyr::select(raw.value,trait.value,
                                                         trait,trait.param) %>%unique(), 
                 aes(y=raw.value,
                     x=trait.value,
                     color=raw.value),
                 shape=16,
                 size=1,alpha=1) + 
      geom_abline(data=Inter.trait.df.long.i %>%
                    group_by(trait,trait.param) %>%
                    summarise(trait.coeff=median(trait.coeff),
                              Intercept=median(Intercept)),
                  aes(slope=trait.coeff, intercept=Intercept),
                  size=1.5) +
      geom_ribbon(data=Inter.trait.df.long.i %>%
                    group_by(trait,trait.param) %>%
                    summarise(trait.coeff.median= median(trait.coeff),
                              trait.coeff.min = quantile(trait.coeff,0.05),
                              trait.coeff.max = quantile(trait.coeff,0.95),
                              Intercept=median(Intercept))%>%
                    ungroup() %>%
                    cross_join(data.frame(trait.value = c(-3:3))) %>%
                    mutate(y= Intercept + trait.coeff.median*trait.value,
                           ymin= Intercept + trait.coeff.min*trait.value,
                           ymax= Intercept + trait.coeff.max*trait.value), 
                  aes(y=y, x=trait.value,ymin=ymin,ymax=ymax),
                  linetype=2, alpha=0.1) +
      facet_grid(addline_format(trait) ~ trait.param,scale="free") +
      geom_hline(yintercept=0,color="grey12",linetype="dashed") +
      labs(x="Trait value or absolute value of the difference",
           y=paste0("Interspecific interactions with ", n ,"density of emitter individuals"))+
      scale_color_gradientn(colours = rev(wes_palette("Zissou1", 
                                                     101, 
                                                     type = "continuous")))+
      theme_few() +
      theme(strip.text.y = element_text(angle=0),
            plot.background = element_rect(color="white",fill="white"),
            panel.background = element_rect(color="white",fill="white"),
            panel.grid.major.x = element_line(color="gray90"))
    print(Cool.detailed.theory.trait.plotlist[[paste0(country,"_Inter")]])
    
    Cool.detailed.theory.trait.plotlist[[paste0(country,"_Intra")]] <- ggplot() +
      geom_point(data=Intra.trait.df.i %>% dplyr::select(raw.value,trait.value,trait) %>%unique(), 
                 aes(y=raw.value,
                     x=trait.value,
                     color=raw.value),
                 shape=16,
                 size=3,alpha=1) + 
      geom_abline(data=Intra.trait.df.long.i%>%
                    group_by(trait) %>%
                    summarise(trait.coeff=median(trait.coeff),
                              Intercept=median(Intercept)),
                  aes(slope=trait.coeff, intercept=Intercept),
                  size=1.5) +
      geom_ribbon(data=Intra.trait.df.long.i%>%
                    group_by(trait) %>%
                    summarise(trait.coeff.median= median(trait.coeff,na.rm=T),
                              trait.coeff.min = quantile(trait.coeff,0.05),
                              trait.coeff.max= quantile(trait.coeff,0.95),
                              Intercept=median(Intercept,na.rm=T)) %>%
                    ungroup() %>%
                   cross_join(data.frame(trait.value = c(-3:3))) %>%
                  mutate(y= Intercept + trait.coeff.median*trait.value,
                         ymin= Intercept + trait.coeff.min*trait.value,
                         ymax= Intercept + trait.coeff.max*trait.value), 
                  aes(y=y, x=trait.value,ymin=ymin,ymax=ymax),
                  size=1.5,alpha=0.2) +
      facet_grid(addline_format(trait) ~ .,scale="free") +
      geom_hline(yintercept=0,color="grey12",linetype="dashed") +
      labs(x="Trait value of focal species",
           color="Response coefficient",
           y=paste0("Intraspecific interaction with ", n ," density of individuals"))+
      scale_color_gradientn(colours = rev(wes_palette("Zissou1", 
                                                      101, 
                                                      type = "continuous")))+
      theme_few() +
      theme(strip.text.y = element_text(angle=0),
            plot.background = element_rect(color="white",fill="white"),
            panel.background = element_rect(color="white",fill="white"),
            panel.grid.major.x = element_line(color="gray90"))
    print(Cool.detailed.theory.trait.plotlist[[paste0(country,"_Intra")]])
    
    Cool.detailed.theory.trait.plotlist[[paste0(country,"_Lambda")]] <- ggplot() +
      facet_grid(addline_format(trait) ~ .,scale="free") +
      geom_point(data=Lambda.trait.df.i %>% dplyr::select(raw.value,trait.value,trait) %>%unique(), 
                 aes(y=raw.value, x=trait.value,color=raw.value),
                 shape=16,size=3,alpha=1) + 
      geom_abline(data=Lambda.trait.df.long.i %>%
                    group_by(trait) %>%
                    summarise(trait.coeff=median(trait.coeff,na.rm=T), 
                              Intercept=median(Intercept,na.rm=T)),
                  aes(slope=trait.coeff, intercept=Intercept),
                  size=1.5) +
      geom_ribbon(data=Lambda.trait.df.long.i%>%
                    group_by(trait) %>%
                    summarise(trait.coeff.median= median(trait.coeff,na.rm=T),
                              trait.coeff.min = quantile(trait.coeff,0.05),
                              trait.coeff.max= quantile(trait.coeff,0.95),
                              Intercept=median(Intercept,na.rm=T)) %>%
                    ungroup() %>%
                    cross_join(data.frame(trait.value = c(-3:3))) %>%
                    mutate(y= Intercept + trait.coeff.median*trait.value,
                           ymin= Intercept + trait.coeff.min*trait.value,
                           ymax= Intercept + trait.coeff.max*trait.value), 
                  aes(y=y, x=trait.value,ymin=ymin,ymax=ymax),
                  size=1.5,alpha=0.2) +
      labs(x="Trait value of focal species",
           y=paste0("Intrinsic fecundity")) +
      scale_color_gradientn(colours = rev(wes_palette("Moonrise3", 
                                                      101, 
                                                      type = "continuous"))) +
      theme_few() +
      theme(strip.text.y = element_text(angle=0),
            plot.background = element_rect(color="white",fill="white"),
            panel.background = element_rect(color="white",fill="white"),
            panel.grid.major.x = element_line(color="gray90"))
    print(Cool.detailed.theory.trait.plotlist[[paste0(country,"_Lambda")]])
    dev.off()
    
  
   
   
    Cool.detailed.theory.trait.plotlist[[paste(country,"_",n)]] <- plot_grid(ggplot() +
                geom_point(data=Inter.trait.df.i %>%
                             mutate(trait.param.label = case_when(trait.param == "emitter.trait.scaled" ~ "Trait value of emitter",
                                                                  trait.param == "receiver.trait.scaled" ~ "Trait value of focal/receiver",
                                                                  trait.param == "scaled.trait.dist" ~ "Trait value of focal - Trait value of emitter")),
                           aes(y=raw.value,
                               x=trait.value,
                               color=trait),
                           shape=16,
                           size=3,alpha=0.2) +
      geom_abline(data=Inter.trait.df.long.i %>%
                    group_by(trait,trait.param) %>%
                    summarise(trait.coeff=median(trait.coeff),
                              Intercept=median(Intercept)) %>%
                    mutate(trait.param.label = case_when(trait.param == "emitter.trait.scaled" ~ "Trait value of emitter",
                                                         trait.param == "receiver.trait.scaled" ~ "Trait value of focal/receiver",
                                                         trait.param == "scaled.trait.dist" ~ "Trait value of focal - Trait value of emitter")),
                  aes(slope=trait.coeff, intercept=Intercept,
                      color=trait),
                  size=1.5) +
      facet_wrap(. ~ trait.param.label,scale="free",strip.position = "bottom") +
      geom_hline(yintercept=0,color="grey12",linetype="dashed") +
      labs(title= paste0("INTERspecific interactions with ", n ," density of emitter individuals"),
           x="",
           y=paste0("Interaction effect on the focal/receiver's fecundity"))+
      scale_color_manual(values= dummy.col ) +                       
      theme_few() +
      theme(strip.placement = "outside",
            legend.position="none",
            plot.title = element_text(size=20),
            axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            strip.text = element_text(size=16),
            panel.grid.major.x = element_line(color="gray90")),
     ggarrange( ggplot() +
                   geom_point(data=Intra.trait.df.i%>%
                                mutate(trait=as.factor(trait)),
                              aes(y=raw.value,
                                  x=trait.value,
                                  color=trait),
                              shape=16,
                              size=3,alpha=0.2) +
                   geom_abline(data=Intra.trait.df.long.i%>%
                                 group_by(trait) %>%
                                 summarise(trait.coeff=median(trait.coeff),
                                           Intercept=median(Intercept))%>%
                                 mutate(trait=as.factor(trait)),
                    aes(slope=trait.coeff, intercept=Intercept,
                        color=trait),
                    size=1.5) +
        geom_hline(yintercept=0,color="grey12",linetype="dashed") +
        labs(title=paste0( "INTRAspecific interactions with \n", n ," density of focal individuals"),
             x="Trait value of focal species",
             y=paste0("Interaction effect on the focal's fecundity"))+
          scale_color_manual(values= dummy.col ) +                       
          theme_few() +theme(legend.position="none",
                             plot.title = element_text(size=20),
                             axis.text=element_text(size=16),
                axis.title=element_text(size=16),
                strip.text = element_text(size=16),
                panel.grid.major.x = element_line(color="gray90")),
      ggplot() +
        geom_point(data=Lambda.trait.df.i %>%
                     mutate(trait=as.factor(trait)),
                   aes(y=raw.value,
                       x=trait.value,color=trait),
                   shape=16, 
                   size=3,alpha=0.5) +
        geom_abline(data=Lambda.trait.df.long.i%>%
                      group_by(trait) %>%
                      summarise(trait.coeff=median(trait.coeff),
                                Intercept=median(Intercept)) %>%
                      mutate(trait=as.factor(trait)),
                    aes(slope=trait.coeff, intercept=Intercept,
                        color=trait),
                    size=1.5) +
        labs(title= "Intrinsic fecundity\n",
             x="Trait value of focal species",
             y=paste0("Intrinsic fecundity of the focal")) +
        scale_color_manual(values= dummy.col ) +                       
        theme_few() +theme(legend.position="none",
                           axis.text=element_text(size=16),
                           plot.title = element_text(size=20),
              axis.title=element_text(size=16),
              strip.text = element_text(size=16),
              panel.grid.major.x = element_line(color="gray90")),
      ncol=2,
      common.legend = T,legend="bottom"),
      #ggpubr::get_legend(legend.plot),
     ncol = 1,
     rel_heights =c(1,0.8),
     labels = '')
    Cool.detailed.theory.trait.plotlist[[paste(country,"_",n)]]
    ggsave(
      Cool.detailed.theory.trait.plotlist[[paste(country,"_",n)]],
      file=paste("figures/GLM/CombinedGLM_details_",country,"_",n,".pdf"),
      width=13,
      height=12,
      unit="in")
  }
}

#---- 1.3. Make summary of glm plots ----
#---- 1.3.1 Make summary tables to do plots----
summary.table.for.plot.glm <- list()
for( country in country.list){

  Inter.trait.df <- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
    mutate(trait = factor(trait, levels=names(dummy.names)))

    df.i <- Inter.trait.df %>%
      gather(c("Receiver.trait.scaled","Emitter.trait.scaled","Scaled.trait.dist","Intercept","Intercept.dist"),
             key="parameters",
             value="estimate") %>%
    mutate(parameters  = case_when(parameters =="Receiver.trait.scaled" ~ "Focal trait",
                                   parameters =="Emitter.trait.scaled" ~ "Emitter trait",
                                   parameters =="Scaled.trait.dist" ~ "Focal trait -\nEmitter trait",
                                   parameters =="Intercept" ~ "intercept",
                                   parameters =="Intercept.dist" ~ "intercept.dist",
                                   T ~ parameters)) %>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))
  
  Intra.trait.df.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
    gather(c("Receiver.trait.scaled","Intercept"),
           key="parameters",
           value="estimate") %>%
    mutate(parameters  = case_when( parameters =="Intercept" ~ "intercept",
                                    parameters =="Receiver.trait.scaled" ~ "Focal trait",
                                    T ~ parameters))%>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))
  
  Lambda.trait.df.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
    gather(c("Receiver.trait.scaled","Intercept"),
           key="parameters",
           value="estimate") %>%
    mutate(parameters  = case_when(parameters =="Receiver.trait.scaled" ~ "Focal trait",
                                   parameters =="Intercept" ~ "intercept",
                                   T ~ parameters))%>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))

  summary.table.for.plot.glm[[country]] <- list(Lambda.trait.df.i=Lambda.trait.df.i,
                                               Intra.trait.df.i=Intra.trait.df.i,
                                               df.i=df.i)
  
}

dummy.col <- c("SRL"="#D55E00","SRA"=  "#F28E2BFF",
               "Root length"= "#F0B878FF","Root tips"= "#C07838FF",
               "Root diameter"= "#DAA51BFF", "Root mass density"="#FFD94AFF",
               "Floret width"= "#B276B2FF","Seed mass"="#FF9DA7FF",
               "C13 water use efficiency"="#24796CFF",#"#9C755FFF",
               "Leaf C to N ratio"= "#72874EFF",#"#BCBD22FF" ,
               "Leaf area index"="#D4E157FF" ,"Canopy shape"="#009E73",
               "SLA"="#99C945FF","Stem height"="#385028FF",
               "intercept"="grey80") ##32A251FF"
dummy.names <- c("SRL"="Specific root length",
                 "SRA"="Specific root area" ,
                 "Root length"= "Root length","Root tips"="Root tips",
                 "Root diameter"="Root diameter" , "Root mass density"="Root tissue density",
                 "Floret width"= "Flower width" ,"Seed mass"="Seed mass",
                 "C13 water use efficiency"="Water use efficiency",
                 "Leaf C to N ratio"= "Nitrogen use efficiency" ,
                 "Leaf area index"="Leaf area index" ,"Canopy shape"="Canopy shape",
                 "SLA"="Specific leaf area","Stem height"="Stem height")

#---- 1.4. FIG 1 - Make graph for main text -INTER - INTRA - horyzontal----
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

trait.to.keep <- c("C13 water use efficiency",
                     "SLA",
                     "Stem height",
                     "Root mass density", 
                     "Floret width",
                     "SRL")
library(purrr)
library(broom)
library(overlapping)
parameter.n="Focal trait"

probability.of.superiority.y.overlapping <- function(x,y){
  out.p.o.s <- NULL
  for(i in 1000){
  vector.of.difference <- data.frame(x.random=sample(x,size=4000, 
                                                     replace =F),
                                     y.random=sample(y,size=4000, 
                                                     replace =F)) %>% 
    mutate(diff =x.random-y.random,
      diff.prob = case_when(diff<0 ~ 1,
                                 diff > 0 ~ 0,
                                 diff==0 ~ 0.5))
  out.p.o.s.n <- data.frame(boot=i,
             probability.of.superiority = mean(vector.of.difference$diff.prob))
  out.p.o.s <- bind_rows(out.p.o.s,out.p.o.s.n)
  }
  out <- boot.overlap( list(x=x,y=y),
                B = 1000, pairsOverlap = FALSE)
  
  return(as.data.frame(out$OVboot_stats) %>% 
         mutate(probability.of.superiority = mean(out.p.o.s.n$probability.of.superiority),
                SE.probability.of.superiority = sd(out.p.o.s.n$probability.of.superiority)/sqrt(1000)))
}
result.ks.test <- NULL
set.seed(16)
for(parameter.n in  c("Focal trait","Emitter trait","Focal trait -\nEmitter trait")) {
#perform Kolmogorov-Smirnov test
D1 <- summary.table.for.plot.glm[["aus"]]$df.i %>%
  dplyr::filter(parameters %in% parameter.n &
                  density.quantile  %in% c("low") &
                  trait %in% trait.to.keep) %>%
  dplyr::select(trait,estimate) %>%
  group_by(trait) %>% mutate(id = row_number()) %>%
  ungroup() %>%
  tidyr::spread(key=trait,value=estimate) %>%
  as.data.frame()

D2 <- summary.table.for.plot.glm[["spain"]]$df.i %>%
  dplyr::filter(parameters %in% parameter.n &
                  density.quantile  %in% c("low") &
                  trait %in% trait.to.keep) %>%
  dplyr::select(trait,estimate) %>%
  group_by(trait) %>% mutate(id = row_number()) %>%
  ungroup() %>%
  tidyr::spread(key=trait,value=estimate)  %>%
  as.data.frame
names(D2)
# Loop through each column
result.ks.test.n <- colnames(D2)[colnames(D2) %in% names(dummy.names)] %>%
  set_names() %>% 
   map(~probability.of.superiority.y.overlapping(D1[, .x], D2[, .x])) |> 
  list_rbind() %>%
  mutate(trait=colnames(D2)[colnames(D2) %in% names(dummy.names)],
         test="Australia > Spain")
  # apply `ks.test` function for each column pair
  #map(~ ks.test(D1[, .x], D2[, .x])) %>%
  # extract test results using `tidy` then bind them together by rows
  #map_dfr(., broom::tidy, .id = "parameter")
result.ks.test <- bind_rows(result.ks.test,
                            result.ks.test.n %>% mutate(parameters =parameter.n))

}

names(result.ks.test)
result.ks.test <- result.ks.test %>%
  rename("Estimate % of overlap"="estOV",
         "Probabitlity of superiority"="probability.of.superiority",
         "direction of superiority"="test") %>%
  dplyr::select("parameters","trait",
                "Estimate % of overlap","bias","se",
                "Probabitlity of superiority","direction of superiority") %>%
  mutate(trait.names=dummy.names[trait]) 

view(result.ks.test)
write.csv(result.ks.test,
          "results/df.meta.analysis.distribution.csv")
result.ks.test <- read.csv("results/df.meta.analysis.distribution.csv")


# 

# test mass distribtuion of each distribution 
mass.distribution.df <- summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(country="Spain")) %>%
    mutate(model="inter",
           pointshape = "Heterospecific interactions") %>%
    dplyr::filter(parameters %in% c("Focal trait","Emitter trait","Focal trait -\nEmitter trait")) %>%
    dplyr::filter(density.quantile  %in% c("low")) %>%
    group_by(country, parameters,trait) %>%
    summarise(estimate.positive = length(estimate[estimate>0])/length(estimate),
              estimate.neg = length(estimate[estimate<0])/length(estimate),
              estimate.median = median(estimate),
              estimate.probabilitydirection = p_direction(estimate,method = "direct",null = 0,
                as_p = FALSE, remove_na = TRUE)$pd) %>%
    #dplyr::filter(trait %in% trait.to.keep ) %>%
    ungroup() %>% 
    dplyr::filter(trait %in% trait.to.keep ) %>%
    mutate(trait.names=dummy.names[trait]) %>%
    mutate(significance.pos = case_when((estimate.positive < 0.1 & estimate.positive > 0.05)~"˙", #paste0(round(estimate.positive,digits=3)," *"),
                                        (estimate.positive < 0.05 & estimate.positive> 0.001) ~"*", #paste0(round(estimate.positive,digits=3)," **"),
                                             estimate.positive < 0.001 ~"**", #paste0(round(estimate.positive,digits=3)," ***"),
                                             T ~ ""),
           significance.neg = case_when((estimate.neg < 0.1 & estimate.neg > 0.05) ~"˙", #paste0(round(estimate.neg,digits=3)," *"),
                                             (estimate.neg < 0.05 & estimate.neg > 0.001)~"*", #paste0(round(estimate.neg,digits=3)," **"),
                                             estimate.neg < 0.001 ~"**", #paste0(round(estimate.neg,digits=3)," ***"),
                                             T ~ ""),
           
           significance = case_when(((estimate.positive < 0.1 & estimate.positive > 0.05)|
                                      (estimate.neg < 0.1 & estimate.neg > 0.05))~"˙",
                                    ((estimate.positive < 0.05 & estimate.positive> 0.001)|
                                       (estimate.neg < 0.05 & estimate.neg > 0.001))~"*",
                                    ((estimate.positive < 0.001)|
                                       (estimate.neg < 0.001))~"**",
                                    T ~"")) 

mass.distribution.df <- mass.distribution.df%>%
  mutate(trait.names=dummy.names[trait])  %>%
  mutate(estimate.positive = round(estimate.positive*100,digits=3),
         estimate.neg= round(estimate.neg*100,digits=3),
         estimate.probabilitydirection = round(estimate.probabilitydirection*100,digits=3))
view(mass.distribution.df)
write.csv(mass.distribution.df,
          "results/df.mass.distribution.csv")

  #view(  mass.distribution.df)
data.subpanel.label <- data.frame(country=c(rep("Australia",18),rep("Spain",18)),
             text.label = c("a.","b.","c.",
                            "d.","e.","f.",
                            "g.","h.","i.",
                            "j.","k.","l.",
                            "m.","n.","o.",
                            "p.","q.","r.", rep("",18)),
             trait.names = rep(rep(c("Specific leaf area","Water use efficiency",
                                     "Flower width","Stem height",
                                     "Specific root length","Root tissue density"),each=3),times=2),
             parameters = rep(c("Focal trait","Emitter trait","Focal trait -\nEmitter trait"),times=12))
  
  
  inter.plot <- summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(country="Spain")) %>%
    mutate(model="inter",
           pointshape = "Heterospecific interactions") %>%
    dplyr::filter(trait %in% trait.to.keep ) %>%
    dplyr::filter(parameters %in% c("Focal trait","Emitter trait","Focal trait -\nEmitter trait")) %>%
    dplyr::filter(density.quantile  %in% c("low")) %>%
    mutate(trait.names=dummy.names[trait]) %>%
    right_join( mass.distribution.df, by=c("country", "parameters","trait","trait.names")) %>%
    mutate(y_numb =case_when(country=="Spain" ~0.4, T~0),
           y_trait=((as.numeric(trait))+y_numb)) %>%
    left_join(result.ks.test, by=c("parameters","trait.names")) %>%
    mutate(rect.color=case_when(`Estimate % of overlap`<0.5 ~ "Opposite",
                                 T~"Consistent")) %>%
    ggplot(aes(y=factor(country,levels=c("Spain","Australia")),
               x=estimate,
               fill=stat(x))) + 
    geom_density_ridges_gradient(scale=c(0.8),
                                 rel_min_height = 0.005) +
    scale_fill_gradientn("Effect direction",
                         colors = c(wes_palette("Zissou1",2,type = "continuous")[2],
     "white",wes_palette("Zissou1",2,type = "continuous")[1]),
      limits=c(-0.07,0.07),
     breaks=c(-0.07,0,0.07),
     labels=c("Negatively correlated with facilitation","Neutral",
              "Positively correlated with facilitation"))+
    facet_grid(factor(addline_format(trait.names),
                      addline_format(c("Specific leaf area","Water use efficiency",
                                       "Flower width","Stem height",
                                       "Specific root length","Root tissue density"))) ~ factor(parameters, 
                                                                             c("Focal trait","Emitter trait","Focal trait -\nEmitter trait"))) +
    geom_vline(xintercept=0) + 
    scale_x_continuous(breaks=c(-0.06,0,0.06)) +
    coord_cartesian(xlim=c(-0.07,0.07)) +
    geom_point(data=mass.distribution.df,
               aes(x=0.04,#estimate.median,
                   y=country,
                   shape=factor(significance,levels=c("˙","*",""))),
               position=position_nudge(y= .4,
                                       x= 0.015),
               colour="black", 
               size=10) +
    geom_text(data=data.subpanel.label,
              aes(x=-0.06,#estimate.median,
                  y=country,
                  label=text.label),
              fontface ="bold",
              position=position_nudge(y= .62),
              colour="black", 
              size=5) +
    scale_shape_manual("Direction significance",
                       values=c("˙","*",""),
                       labels=c("less than 5% on one side",
                                "less than 10% on one side","")) +
    geom_rect(mapping=aes(xmin=-0.07,xmax=0.07,#estimate.median,
                  ymin=0.75,ymax=2.85,
                  color=factor(rect.color,
                               levels=c("Opposite","Consistent",""))),
              fill=NA) +
    scale_color_manual("Country comparison",
                       values=c("#F0E442","white"),
                       labels=c("Less than 50% overlap",""),
                       breaks=c("Opposite","Consistent")) +
    labs(y="",
         x="Per capita effect size on heterospecific interactions") +
    guides(fill= guide_legend(title.position = "top",
                              nrow=3,order=1,
                              override.aes = list(shape="")),
           color = guide_legend(title.position = "top",
                                nrow=2,
                                order=2),
           shape= guide_legend(title.position = "top",
                               nrow=2,
                               order=3)) +
    theme_bw() +
    theme(legend.position="bottom",
          strip.text.x = element_text(size=16),
          strip.text.y = element_text(size=16,angle=0),
          strip.background = element_rect(fill="grey95",color="white",linewidth=4),
          panel.spacing = unit(3, "lines"),
          legend.title =element_text(size=20),
          legend.text =element_text(size=18),
          axis.title=element_text(size=20),
          axis.text.y = element_text(size=18),
          axis.text.x = element_text(size=18,angle=90),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = NA,
                                       fill=NA),
          panel.grid.major.x = element_blank(), #element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank(),
          panel.spacing.x=unit(1, "lines"),
          panel.spacing.y=unit(0, "lines"),
          plot.margin=unit(c(1.5,1,0.5,1),"cm"))
  
inter.plot

#figures/main/TraitEffectSize_distribution.pdf

#figures/main/GLM.traits.INTER.bis.pdf
#----Intra---
mass.distribution.intra.df <- summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
              mutate(country="Spain")) %>%
  mutate(model="intra",
         pointshape = "Conspecific interactions") %>%
  dplyr::filter(parameters %in% c("Focal trait")) %>%
  dplyr::filter(density.quantile  %in% c("low")) %>%
  group_by(country, parameters,trait) %>%
  summarise(estimate.positive = length(estimate[estimate>0])/length(estimate),
            estimate.neg = length(estimate[estimate<0])/length(estimate),
            estimate.median = median(estimate)) %>%
  #dplyr::filter(trait %in% trait.to.keep ) %>%
  ungroup() %>% 
  dplyr::filter(trait %in% trait.to.keep ) %>%
  mutate(trait.names=dummy.names[trait]) %>%
  mutate(significance.pos = case_when((estimate.positive < 0.1 & estimate.positive > 0.05)~"˙", #paste0(round(estimate.positive,digits=3)," *"),
                                      (estimate.positive < 0.05 & estimate.positive> 0.001) ~"*", #paste0(round(estimate.positive,digits=3)," **"),
                                      estimate.positive < 0.001 ~"**", #paste0(round(estimate.positive,digits=3)," ***"),
                                      T ~ ""),
         significance.neg = case_when((estimate.neg < 0.1 & estimate.neg > 0.05) ~"˙", #paste0(round(estimate.neg,digits=3)," *"),
                                      (estimate.neg < 0.05 & estimate.neg > 0.001)~"*", #paste0(round(estimate.neg,digits=3)," **"),
                                      estimate.neg < 0.001 ~"**", #paste0(round(estimate.neg,digits=3)," ***"),
                                      T ~ ""),
         
         significance = case_when(((estimate.positive < 0.1 & estimate.positive > 0.05)|
                                     (estimate.neg < 0.1 & estimate.neg > 0.05))~"˙",
                                  ((estimate.positive < 0.05 & estimate.positive> 0.001)|
                                     (estimate.neg < 0.05 & estimate.neg > 0.001))~"*",
                                  ((estimate.positive < 0.001)|
                                     (estimate.neg < 0.001))~"**",
                                  T ~"")) 

view(  mass.distribution.intra.df)
intra.plot <- summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
              mutate(country="Spain")) %>%
  mutate(model="intra",
         pointshape = "Conspecific interactions") %>%
  dplyr::filter(parameters %in% c("Focal trait")) %>%
  dplyr::filter(density.quantile  %in% c("low")) %>%
  left_join( mass.distribution.intra.df %>%
                select(country,parameters,trait,significance), 
              by=c("country", "parameters","trait")) %>%
  dplyr::filter(trait %in% trait.to.keep ) %>%
  mutate(trait.names=dummy.names[trait]) %>%
  mutate(rect.color=case_when((trait=="Root mass density" & parameters=="Focal trait") ~ "Consistent",
                              T~"")) %>%
  ggplot(aes(y=country,
             x=estimate,
             fill=stat(x))) + 
  geom_density_ridges_gradient(scale=c(0.8),
                               rel_min_height = 0.005) +
  scale_fill_gradientn("Direction of effect",
                       colors = c(wes_palette("Zissou1",2,type = "continuous")[2],
                                  "white",wes_palette("Zissou1",2,type = "continuous")[1]),
                       limits=c(-0.02,0.02),
                       breaks=c(-0.02,0,0.02),
                       labels=c("Negatively correlated with facilitation","Neutre",
                                "Positively correlated with facilitation"))+
  facet_grid(factor(addline_format(trait.names),
                    addline_format(c("Flower width","Stem height",
                                     "Specific leaf area","Water use efficiency",
                                     "Specific root length","Root tissue density"))) ~ factor(parameters, 
                                                                                              c("Focal trait"))) +
  geom_vline(xintercept=0) + 
  scale_x_continuous(breaks=c(-0.02,0,0.02)) +
  coord_cartesian(xlim=c(-0.02,0.02)) +
  geom_text(data=mass.distribution.intra.df,
            aes(x=0.01,#estimate.median,
                y=country,
                label=significance),
            fontface ="bold",
            position=position_nudge(y= .4,
                                    x= 0.005),
            colour="black", 
            size=10) +
  labs(y="",
       x="Per capita effect size on conspecific interactions",
       shape="Country",
       color="Functional trait") +
  guides(fill= guide_legend(title.position = "top",
                            nrow=3)) +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16,angle=0),
        strip.background = element_rect(fill="grey95",color="white",linewidth=4),
        panel.spacing = unit(3, "lines"),
        legend.title =element_text(size=20),
        legend.text =element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=90),
        panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
        panel.border = element_rect( color = NA,
                                     fill=NA),
        panel.grid.major.x = element_blank(), #element_line(colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0, "lines"),
        plot.margin=unit(c(1.5,1,0.5,1),"cm"))

intra.plot

#figures/main/TraitEffectSizeIntra_distribution.pdf
ggsave(GLM.traits.INTER.INTRA,
       width=13.48,
       height=10.48,
       unit="in",
       file="figures/main/GLM.traits.INTER.png")

library(grid)
library(pBrackets) 
#figures/GLM.traits.INTER.pdf
#---- 1.5. FIG R2 - Make graph for main text -INTER - DIAGONAL ----
#---- 1.5.1 stat ----

library(kableExtra)
library(knitr) 
library(scales)
library(HDInterval)
test.list <- list()
Hdi.obs <- function(x,Q){
  interval <- HDInterval::hdi(x, credMass = Q)
  qvec <- x[ which( x> interval[1] & x< interval[2])]

}
outer_product_by_group <- function(data, group_var, x_var, y_var, fun = "-") {
  # Split data by group
  split_data <- split(data, data[[group_var]])
  
  # Apply outer product to each group
  outer_products <- lapply(split_data, function(group_data) {
    v <- outer(group_data[[x_var]], group_data[[y_var]], fun)
    length(v[which(v>0)])/length(v)
  })
  
  return(outer_products)
}

inter_diag_df <- summary.table.for.plot.glm[["aus"]]$df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
              mutate(country="Spain")) %>%
  dplyr::filter(density.quantile  %in% c("low","high")) 


for( parameters.n in c("Focal trait","Emitter trait",
                       "Focal trait -\nEmitter trait")){
inter_diag_df.i <-  inter_diag_df %>%
  dplyr::filter(parameters ==parameters.n) %>%
  group_by(trait,country,density.quantile) %>%
  reframe(estimate=Hdi.obs(abs(estimate),0.8)) %>%
  group_by(trait,country,density.quantile) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = density.quantile, 
            values_from = c('estimate'), names_sep="")  %>%
  mutate(group=paste0(trait,"_",country)) 

#view(inter_diag_df.i )
#str(inter_diag_df.i )
 
result.df <- outer_product_by_group(inter_diag_df.i,"group",
                         "low", "high","-")
result.df <- as.data.frame(result.df) %>%
            gather(key="group",value='percentage.low.vs.high') %>%
            separate(group, into=c("trait","country"),sep="_")
result.df <- bind_rows(data.frame(trait="All",country="both",
                       'percentage.low.vs.high' = mean(result.df[,"percentage.low.vs.high"]),
                       sd = round(sd(result.df[,"percentage.low.vs.high"])*100,digits=2)),
                       result.df) %>%
  mutate(percentage.low.vs.high = round(percentage.low.vs.high*100,digits=2))


#test <- brms::hypothesis(inter_diag_df.i,"low > high")
#tab = test$hypothesis %>% select(-Star)
#a = map_chr(tab, ~ifelse(class(.x)=="numeric", "r","l"))

#tab = tab %>% 
#  mutate(across(where(is.numeric), ~comma(., accuracy=0.01))) %>% 
 # rename_all(~gsub("\\.", " ", .))

test.list[[parameters.n]] <- list(#test = test,
                                 #tab=tab,
                                 result.df=result.df)

}


result.df.all <- left_join(test.list[["Focal trait"]]$result.df %>% rename("Focal.trait.perc.low.vs.high" = 'percentage.low.vs.high',
                                                                           "Focal.trait.sd"="sd"),
                           left_join(test.list[["Emitter trait"]]$result.df%>% rename("Neigh.trait.perc.low.vs.high" = 'percentage.low.vs.high',
                                                                      "Neigh.trait.sd"="sd"),
                                      test.list[["Focal trait -\nEmitter trait"]]$result.df%>% rename("Dist.trait.perc.low.vs.high" = 'percentage.low.vs.high',
                                                                                     "Dist.trait.sd"="sd")))
result.df.all

write.csv(result.df.all,
          file="results/Percentage.diagonal.csv")

result.df.all <-read.csv(file="results/Percentage.diagonal.csv")

#---- 1.5.2 Fig----

Cool.glm.theory.trait.plotlist <- list()

for( parameters.n in c("Focal trait","Emitter trait","Focal trait -\nEmitter trait")){
Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag",parameters.n)]]  <- summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(country="Spain")) %>%
    dplyr::filter(density.quantile  %in% c("low","high")) %>%
    #dplyr::filter(trait %in% trait.to.keep ) %>%
    dplyr::filter(parameters==c(parameters.n)) %>%
    dplyr::select(trait,country,density.quantile,estimate) %>%
    group_by(trait,country,density.quantile) %>%
    #mutate(ID = row_number()) %>%
    summarise(median.est=median(abs(estimate)),
              low.est = HDInterval::hdi(abs(estimate),0.8)[1],
              up.est = HDInterval::hdi(abs(estimate),0.8)[2]) %>%
    ungroup() %>%
    mutate(trait=factor(trait,
                      levels=rev(names( dummy.col)))) %>%
    pivot_wider(names_from = density.quantile, 
                values_from = c('median.est',"low.est","up.est"), names_sep=".") %>%
    ggplot(aes(y=median.est.high,
               x=median.est.low,
               shape = country,
               color=as.factor(trait),
               group=as.factor(trait) )) + 
    #shape=density.quantile)) +
    geom_point(size=6) +
    geom_pointrange(aes(ymin=low.est.high,
                         ymax=up.est.high),alpha=0.6,size=1) +
   geom_errorbarh(aes(xmin=low.est.low,
                      xmax=up.est.low),alpha=0.6,size=1)+
    geom_abline(intercept=0,slope=1) +
    scale_shape_manual(values=c(16:17)) +
    scale_color_manual(values=dummy.col,
                       labels=dummy.names)+
    scale_y_continuous(breaks=c(0,0.01,0.03),labels=c(0,0.01,0.03)) +
    scale_x_continuous(breaks=c(0,0.01,0.03),labels=c(0,0.01,0.03)) +
    labs(x="Effect size at low neighbor density",
         y="Effect size at high neighbor density",
         shape="Community",
         color="Functional traits") +
    annotate(geom = "polygon", x = c(Inf, -Inf, -Inf), 
           y = c(Inf, -Inf, Inf), fill = "grey70", alpha = 0.1 )+
    annotate(geom = "text",
             label=paste0(test.list[[parameters.n]]$result.df$percentage.low.vs.high[1]," % \u00B1 ", test.list[[parameters.n]]$result.df$sd[1]),
           x=0.024,y=0.001,size=6,angle=0)+
    coord_cartesian( xlim = c(0,0.03), ylim=c(0,0.03),
                   expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           shape = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          plot.margin = unit(c(1,ifelse(country=="aus",
                                        0,3),0,
                               ifelse(country=="aus",
                                      3,0)),"cm"),
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=20),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major =  element_blank(), #element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())

}

Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag","Focal trait")]]
Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag","Emitter trait")]]
Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag","Focal trait -\nEmitter trait")]]
  

ggarrange(Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag","Focal trait")]],
          Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag","Emitter trait")]],
          Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag","Focal trait -\nEmitter trait")]],
          common.legend = T, legend = "bottom",
          label.x = c(-0.02,-0.04,-0.18),
          label.y = 1.01,align="h",
          font.label = list(size = 24, color = "black", 
                            face = "bold", family = NULL),
           labels=c("a. Focal trait", "b. Emitter trait",
                    "c. Focal trait - Emitter trait"),
           ncol=3)
   #figures/main/Oblique.INTER.pdf 
  
#---- 1.6. LAST figure manuscript -----
#---- 1.6.2. FIG R3 - Make graph for main text -INTER - INTRA ----
 
 plot_inter_intra  <-   summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(model="inter") %>%
    bind_rows( summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i%>%
                 mutate(model="intra"))   %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(model="inter") %>%
                bind_rows( summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
                             mutate(model="intra")) %>%
                mutate(country="Spain")) %>%
    dplyr::filter(density.quantile  %in% c("low") &
                    !trait %in%   trait.to.remove) %>%
    dplyr::filter(parameters==c("Focal trait"))  %>%
    dplyr::select(trait,country,estimate,model)  %>%
    group_by(trait,country,model) %>%
   group_by(trait,country,model) %>%
   summarise(median.est=median(estimate),
             low.est = HDInterval::hdi(estimate,0.8)[1],
             up.est = HDInterval::hdi(estimate,0.8)[2]) %>%
   ungroup() %>%
   mutate(trait=factor(trait,
                       levels=rev(names( dummy.col)))) %>%
   pivot_wider(names_from = model, 
               values_from = c('median.est',"low.est","up.est"), names_sep=".") %>%
   ggplot(aes(x=median.est.intra,
              y=median.est.inter)) + 
    geom_point( aes(shape = country,
                    color=trait),
                size=8,alpha=1) +
   geom_errorbar(aes(color=trait,
                     group=trait,
                     ymin=low.est.inter,
                        ymax=up.est.inter),
                    size=1,alpha=0.3) +
   geom_errorbarh(aes(color=trait,
                      group=trait,
                      xmin=low.est.intra,
                       xmax=up.est.intra),
                   size=1,alpha=0.3, height=0) +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    #geom_abline(intercept=0,slope=1) +
    scale_shape_manual(values=c(16,17)) +
    scale_color_manual(values=dummy.col,
                       labels=dummy.names)+
    scale_size(range=c(1,2))+
    scale_y_continuous(breaks=seq(-0.03,0.03,0.01)) +
    scale_x_continuous(breaks=seq(-0.01,0.01,0.005)) +
    labs(x="Effect size on conspecific interactions",
         y="Effect size on heterospecific interactions",
         shape="Community",
         color="Functional traits") +
    annotate(geom="segment",x=0, xend = +0.01 , 
                     y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom="segment",x=0, xend =0 , 
                     y=0, yend = 0.03, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom="segment",x=0, xend = -0.01 , 
                     y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom="segment",x=0, xend =0 , 
                     y=0, yend = -0.03, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom = "label",label="Increased   competition \nfrom   others",
             x=0.00,y=-0.025,size=5,angle=0,
             fill= scales::alpha("grey95", .8),
             color=scales::alpha("black"))+
    annotate(geom = "label",label="Increased   facilitation \nfrom   others",
             x=0.00,y=0.025,size=5,angle=0,
             fill= scales::alpha("grey95", .8),
             color=scales::alpha("black"))+
    annotate(geom = "label",label="Increased \nself-competition",
             x=-0.007,y=0.00,size=5,angle=0,
             fill= scales::alpha("grey95", .8),
             color=scales::alpha("black"))+
    annotate(geom = "label",label="Increased \nself-facilitation",
             x=0.007,y=0.00,size=5,angle=0,
             fill= scales::alpha("grey95", .8),
             color=scales::alpha("black"))+
    coord_cartesian( ylim = c(-0.03,0.03), 
                     xlim=c(-0.01,0.01),
                     expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           shape = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          plot.margin = unit(c(1,1,1,1),"cm"),
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=18),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major = element_blank(), #element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())
  
  
  plot_inter_intra
  #figures/main/Oblique.INTER.INTRA.pdf 
#---- 1.6.3. FIG R3 -Make graph for main text -INTER - LAMBDA ----
  
  trait.to.remove <- c(#"Root diameter",
                       #"Flower width","Seed mass",
                       "Leaf area index",
                       "Canopy shape")
  
  plot_inter_lambda  <- summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(model="inter") %>%
    bind_rows( summary.table.for.plot.glm[["aus"]]$Lambda.trait.df.i %>%
                 mutate(model="lambda")) %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(model="inter") %>%
                bind_rows( summary.table.for.plot.glm[["spain"]]$Lambda.trait.df.i %>%
                             mutate(model="lambda")) %>%
                mutate(country="Spain")) %>%
    dplyr::filter(density.quantile  %in% c("low") &
                    !trait %in%   trait.to.remove) %>%
    dplyr::filter(parameters==c("Focal trait"))  %>%
    dplyr::select(trait,country,estimate,model)  %>%
    group_by(trait,country,model) %>%
    #mutate(ID = row_number()) %>%
    summarise(median.est=median(estimate),
              low.est = HDInterval::hdi(estimate,0.8)[1],
              up.est = HDInterval::hdi(estimate,0.8)[2]) %>%
    ungroup() %>%
    mutate(trait=factor(trait,
                        levels=rev(names( dummy.col)))) %>%
    pivot_wider(names_from = model, 
                values_from = c('median.est',"low.est","up.est"), names_sep=".") %>%
    ggplot(aes(x=median.est.lambda,
               y=median.est.inter)) + 
    geom_errorbar(aes(ymin=low.est.inter,
                      ymax=up.est.inter,
                      group=as.factor(trait),
                      color=as.factor(trait)),
                  size=1,alpha=0.3) +
    geom_errorbarh(aes(xmin=low.est.lambda,
                       xmax=up.est.lambda,
                       group=as.factor(trait),
                       color=as.factor(trait)),
                   size=1,alpha=0.3, height=0) +
    geom_point(aes(group=as.factor(trait),shape = country, color=as.factor(trait)),
               size=8,alpha=1) +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    scale_shape_manual(values=c(16:17)) +
    scale_color_manual(values=dummy.col,
                       labels=dummy.names)+
    scale_y_continuous(breaks=seq(-0.03,0.03,0.01)) +
    #scale_x_continuous(breaks=seq(-0.1,0.1,0.05)) +
    labs(x="Effect size on intrinsic fitness",
         y="Effect size on heterospecific interactions",
         shape="Community",
         color="Functional traits") +
    annotate(geom="segment",x=0, xend = +1 , y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
  annotate(geom="segment",x=0, xend =0 ,  y=0, yend = 0.03, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
  annotate(geom="segment",x=0, xend = -1 , y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
  annotate(geom="segment",x=0, xend =0 ,  y=0, yend = -0.03, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom = "label",label="Increased   competition \nfrom   others",
             x=0.00,y=-0.025,size=5,angle=0,
             fill= scales::alpha("grey95", .8),
             color=scales::alpha("black"))+
    annotate(geom = "label",label="Increased   facilitation \nfrom   others",
             x=0.00,y=0.025,size=5,angle=0,
             fill= scales::alpha("grey95", .8),
             color=scales::alpha("black"))+
    annotate(geom = "label",label="Increased \nintrinsic fecundity",
             x=0.68,y=0,size=5,angle=0,
             fill= scales::alpha("grey95", .8))+
    annotate(geom = "label",label="Decreased \nintrinsic fecundity",
             x=-0.68,y=0,size=5,angle=0,
             fill= scales::alpha("grey95", .8))+
    coord_cartesian( ylim = c(-0.03,0.03),
                     xlim=c(-1,1),
                     expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           shape = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          plot.margin = unit(c(1,1,1,1),"cm"),
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=19),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major = element_blank(), #element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())
  
  plot_inter_lambda
  #figures/main/Oblique.INTER.INTRA.LAMBDA.pdf 
#---- 1.6.4. FIG R3 -Make graph for main text -INTRA - LAMBDA ----
  
  trait.to.remove <- c(#"Root diameter",
    #"Flower width","Seed mass",
    "Leaf area index",
    "Canopy shape")
  
  plot_intra_lambda  <- summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i%>%
    mutate(model="intra") %>%
    bind_rows( summary.table.for.plot.glm[["aus"]]$Lambda.trait.df.i %>%
                 mutate(model="lambda")) %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
                mutate(model="intra") %>%
                bind_rows( summary.table.for.plot.glm[["spain"]]$Lambda.trait.df.i %>%
                             mutate(model="lambda")) %>%
                mutate(country="Spain")) %>%
    dplyr::filter(density.quantile  %in% c("low") &
                    !trait %in%   trait.to.remove) %>%
    dplyr::filter(parameters==c("Focal trait"))  %>%
    dplyr::select(trait,country,estimate,model)  %>%
    group_by(trait,country,model) %>%
    #mutate(ID = row_number()) %>%
    summarise(median.est=median(estimate),
              low.est = HDInterval::hdi(estimate,0.8)[1],
              up.est = HDInterval::hdi(estimate,0.8)[2]) %>%
    ungroup() %>%
    mutate(trait=factor(trait,
                        levels=rev(names( dummy.col)))) %>%
    pivot_wider(names_from = model, 
                values_from = c('median.est',"low.est","up.est"), names_sep=".") %>%
    ggplot(aes(x=median.est.lambda,
               y=median.est.intra)) + 
    geom_errorbar(aes(ymin=low.est.intra,
                      ymax=up.est.intra,
                      group=as.factor(trait),
                      color=as.factor(trait)),
                  size=1,alpha=0.3) +
    geom_errorbarh(aes(xmin=low.est.lambda,
                       xmax=up.est.lambda,
                       group=as.factor(trait),
                       color=as.factor(trait)),
                   size=1,alpha=0.3, height=0) +
    geom_point(aes(group=as.factor(trait),shape = country, color=as.factor(trait)),
               size=8,alpha=1) +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    scale_shape_manual(values=c(16:17)) +
    scale_color_manual(values=dummy.col,
                       labels=dummy.names)+
    #scale_y_continuous(breaks=seq(-0.1,0.1,0.025)) +
    #scale_x_continuous(breaks=seq(-0.1,0.1,0.05)) +
    labs(x="Effect size on intrinsic fitness",
         y="Effect size on conspecific interactions",
         shape="Community",
         color="Functional traits") +
    annotate(geom="segment",x=0, xend = +1 , y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom="segment",x=0, xend =0 , y=0, yend = 0.01, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom="segment",x=0, xend = -1 , y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom="segment",x=0, xend =0 , 
                     y=0, yend = -0.010, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom="label",label="Increased \n self -competition",
              x=0.04,y=-0.0082,size=5,angle=0, 
              fill= scales::alpha("grey95", .8))+
    annotate(geom="label",label="Increased \n self-facilitation",
             x=0.025,y=0.0082,size=5,angle=0,
             fill= scales::alpha("grey95", .8))+
    annotate(geom="label",label="Increased \nintrinsic fecundity",
             x=0.68,y=0,size=5,angle=0, 
             fill= scales::alpha("grey95", .8))+
    annotate(geom="label",label="Decreased \nintrinsic fecundity",
             x=-0.68,y=0,size=5,angle=0, 
             fill= scales::alpha("grey95", .8))+
    coord_cartesian( ylim = c(-0.01,0.01),
                     xlim=c(-1,1),
                     expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           shape = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          plot.margin = unit(c(2,1,1,1),"cm"),
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=19),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major = element_blank(), #element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())
  plot_intra_lambda
#---- 1.6.4. FIG R3 -Make graph for main text -theoric ----
  
  plot_theory <- summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i%>%
    mutate(model="intra") %>%
    bind_rows( summary.table.for.plot.glm[["aus"]]$Lambda.trait.df.i %>%
                 mutate(model="lambda")) %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
                mutate(model="intra") %>%
                bind_rows( summary.table.for.plot.glm[["spain"]]$Lambda.trait.df.i %>%
                             mutate(model="lambda")) %>%
                mutate(country="Spain")) %>%
    dplyr::filter(density.quantile  %in% c("low") &
                    !trait %in%   trait.to.remove) %>%
    dplyr::filter(parameters==c("Focal trait"))  %>%
    dplyr::select(trait,country,estimate,model)  %>%
    group_by(trait,country,model) %>%
    #mutate(ID = row_number()) %>%
    summarise(median.est=median(estimate),
              low.est = HDInterval::hdi(estimate,0.8)[1],
              up.est = HDInterval::hdi(estimate,0.8)[2]) %>%
    ungroup() %>%
    mutate(trait=factor(trait,
                        levels=rev(names( dummy.col)))) %>%
    pivot_wider(names_from = model, 
                values_from = c('median.est',"low.est","up.est"), names_sep=".") %>%
    ggplot(aes(x=median.est.lambda,
               y=median.est.intra)) + 
    geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    scale_shape_manual(values=c(16:17)) +
    scale_color_manual(values=dummy.col,
                       labels=dummy.names)+
    scale_y_continuous(breaks=c(0)) +
    scale_x_continuous(breaks=c(0)) +
    labs(x="Effect size of traits on a parameter",
         y="Effect size of traits on a parameter",
         shape="Community",
         color="Functional traits") +
    geom_segment(aes(x=0, xend = +1 , 
                     y=0, yend = 0), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(aes(x=0, xend =0 , 
                     y=0, yend = 1), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(aes(x=0, xend = -1 , 
                     y=0, yend = 0), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(aes(x=0, xend =0 , 
                     y=0, yend = -1), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom = "text",label="Both parameters have \na positive impact on fitness",
             x=0.5,y=0.5,size=6,angle=0)+
    annotate(geom = "text",label="Both parameters have \na negative impact on fitness",
             x=-0.5,y=-0.5,size=6,angle=0)+
    annotate(geom = "text",label="Trade-off \nbetween both parameters",
             x=0.5,y=-0.5,size=6,angle=0)+
    annotate(geom = "text",label="Trade-off \nbetween both parameters",
             x=-0.5,y=0.5,size=6,angle=0)+
    coord_cartesian( ylim = c(-1,1), xlim=c(-1,1),
                     expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           shape = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          plot.margin = unit(c(2,1,1,1),"cm"),
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=19),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major = element_blank(), #element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())
  #plot_theory
  
  #figures/main/Oblique.INTER.INTRA.LAMBDA.pdf 
  ggarrange(plot_theory,  
            plot_intra_lambda,
            plot_inter_intra,
            plot_inter_lambda,
            ncol=2,
            nrow=2,
            common.legend = T, legend="bottom",
            label.x= c(-0.28,-0.28,
                       -0.28,-0.28),
            label.y= c(1.01,1.01,
                       1.02,1.02),
            align="hv",
            font.label = list(size = 22, color = "black", 
                              face = "bold", family = NULL),
            labels=c("a. Theoretical map of impact on fitness ",
                     "b. Conspecific vs intrinsic fecundity",
                     "c. Heterospecific vs conspecific effect",
                     "d. Heterospecific vs intrinsic fecundity" ))
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2. Supp figures----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#---- 2.1. Observed Trait Density figures----

dummy.names <- list("SRL"="SRL (cm/g)","SRA"="SRA (cm^2/g)" ,"Root length"="Root length \n(cm)",
                 "Root tips"="Root tips" ,
                 "Root diameter"="Root diameter \n(mm)" ,
                 "Root mass density"="Root mass \ndensity \n(mg/mm3)",
                 "Flower width"= "Flower width \n(mm)" ,"Seed mass"="Seed mass \n(mg)",
                 "C13 water use efficiency"="C13 water \nuse efficiency \n(per ml)",
                 "Leaf C to N ratio"= "Leaf C to N ratio" ,
                 "Leaf area index"="Leaf area index" ,
                 "Canopy shape"="Canopy shape",
                 "SLA"="SLA (cm^2/g)","Stem height"="Stem height \n(cm)")

trait_labeller <- function(variable,value){
  return(dummy.names[value])
}

plant_code_spain <- read.csv( "data/spain_rawdata/plant_code_spain.csv",
                              header = T, stringsAsFactors = F, sep=",",
                              na.strings = c("","NA"))
plant_code_aus <- read.csv("data/aus_rawdata/plant_code_aus.csv",
                           header = T, stringsAsFactors = F, sep=",",
                           na.strings = c("","NA"))

library(ggridges)

Cool.theory.trait.df[["spain"]]$trait.dist.df %>%
  dplyr::select(focal,trait,receiver.trait) %>%
    unique() %>%
    left_join(read.csv(paste0("data/spain_rawdata/plant_code_spain.csv"),
                       header = T, stringsAsFactors = F, sep=",",
                       na.strings = c("","NA")) %>% 
                dplyr::select(family,code.analysis) %>%
                unique() %>%
                rename("focal"="code.analysis")) %>%
  mutate(country="Spain") %>%
    bind_rows(Cool.theory.trait.df[["aus"]]$trait.dist.df %>%
                dplyr::select(focal,trait,receiver.trait) %>%
                unique() %>%
                left_join(read.csv(paste0("data/aus_rawdata/plant_code_aus.csv"),
                                   header = T, stringsAsFactors = F, sep=",",
                                   na.strings = c("","NA")) %>% 
                            dplyr::select(family,final.code) %>%
                            unique() %>%
                            filter(!is.na(final.code)) %>%
                            rename("focal"="final.code") %>%
                             mutate(country="Australia"))) %>%
  mutate(ynum= as.numeric(as.factor(country))) %>%
    ggplot() +
    ggridges::geom_density_ridges2(aes(y=country,
                                       x=receiver.trait),scale=1,
                                   fill="grey80",alpha=0.8) +
    geom_point(aes(y=ynum,
                   x=receiver.trait,
                   color=family),alpha=0.9,
               position=position_dodge2(width =0.16,
                                          preserve = "total",
                                          padding = 0.1,
                                          reverse = FALSE),
               shape = 15, size = 4) +

    facet_wrap(.~factor(trait,levels=names(dummy.names)),
               nrow=3,scale="free_x",
               labeller=trait_labeller) +
    scale_color_manual(values =c("#FBE183FF",   "#FE9B00FF","#D8443CFF",
                                 "#E6A2A6FF","#9F5691FF","#633372FF",
                                 "#1F6E9CFF", "#2B9B81FF","#92C051FF")) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(shape=15,
                                                     size=10,
                                                     alpha = 1))) +
   labs(y="", x="observed values") + 
    theme(legend.position="bottom",
          strip.text =element_text(size=16),
          legend.text =element_text(size=16),
          legend.title =element_text(size=20),
          axis.text=element_text(size=16),
          plot.title = element_text(size=20),
          axis.title=element_text(size=16))
#figures/supp/Traits.density.pdf


Cool.theory.trait.df[["spain"]]$trait.dist.df %>%
  dplyr::select(focal,trait,receiver.trait.scaled) %>%
  unique() %>%
  left_join(read.csv(paste0("data/spain_rawdata/plant_code_spain.csv"),
                     header = T, stringsAsFactors = F, sep=",",
                     na.strings = c("","NA")) %>% 
              dplyr::select(family,code.analysis) %>%
              unique() %>%
              rename("focal"="code.analysis")) %>%
  mutate(country="Spain") %>%
  bind_rows(Cool.theory.trait.df[["aus"]]$trait.dist.df %>%
              dplyr::select(focal,trait,receiver.trait.scaled) %>%
              unique() %>%
              left_join(read.csv(paste0("data/aus_rawdata/plant_code_aus.csv"),
                                 header = T, stringsAsFactors = F, sep=",",
                                 na.strings = c("","NA")) %>% 
                          dplyr::select(family,final.code) %>%
                          unique() %>%
                          filter(!is.na(final.code)) %>%
                          rename("focal"="final.code") %>%
                          mutate(country="Australia"))) %>%
  ggplot() +
  ggridges::geom_density_ridges2(aes(y=trait,
                                     x=receiver.trait.scaled),scale=1,
                                 fill="grey80",alpha=0.8) +
  geom_point(aes(y=trait,
                 x=receiver.trait.scaled,
                 color=family),alpha=0.9,
             position=position_dodge2(width =0.16,
                                      preserve = "total",
                                      padding = 0.1,
                                      reverse = FALSE),
             shape = 15, size = 4) +
  facet_wrap(.~country,scale="free_x",nrow=1) +
  scale_color_manual(values =c("#FBE183FF",   "#FE9B00FF","#D8443CFF",
                               "#E6A2A6FF","#9F5691FF","#633372FF",
                               "#1F6E9CFF", "#2B9B81FF","#92C051FF")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape=15,
                                                   size=10,
                                                   alpha = 1))) +
  labs(y="functional traits scaled", x="scaled values") + 
  theme(legend.position="bottom",
        strip.text =element_text(size=16),
        legend.text =element_text(size=16),
        legend.title =element_text(size=20),
        axis.text=element_text(size=16),
        plot.title = element_text(size=20),
        axis.title=element_text(size=16))

#figures/supp/Traits.density.scaled.pdf

#---- 2.2. Make fig 3 for all traits ----

intercept.df <- summary.table.for.plot.glm[[country]]$df.i %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    dplyr::mutate(trait = case_when(parameters =="intercept" ~ "intercept",
                                    T~trait)) %>%
    dplyr::mutate(parameters = case_when(parameters =="intercept" ~ NA,
                                         T~parameters)) %>%
    fill(parameters,.direction = c("up")) %>%
    filter(trait =="intercept") %>%
    mutate(pointshape = "Intercept of interaction strength")
  


mass.distribution.df <- summary.table.for.plot.glm[["aus"]]$df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
              mutate(country="Spain")) %>%
  mutate(model="inter",
         pointshape = "Heterospecific interactions") %>%
  dplyr::filter(parameters %in% c("Focal trait","Emitter trait","Focal trait -\nEmitter trait")) %>%
  dplyr::filter(density.quantile  %in% c("low")) %>%
  group_by(country, parameters,trait) %>%
  summarise(estimate.positive = length(estimate[estimate>0])/length(estimate),
            estimate.neg = length(estimate[estimate<0])/length(estimate),
            estimate.median = median(estimate)) %>%
  #dplyr::filter(trait %in% trait.to.keep ) %>%
  ungroup() %>% 
  mutate(trait.names=dummy.names[trait]) %>%
  mutate(significance.pos = case_when((estimate.positive < 0.1 & estimate.positive > 0.05)~"˙", #paste0(round(estimate.positive,digits=3)," *"),
                                      (estimate.positive < 0.05 & estimate.positive> 0.001) ~"*", #paste0(round(estimate.positive,digits=3)," **"),
                                      estimate.positive < 0.001 ~"**", #paste0(round(estimate.positive,digits=3)," ***"),
                                      T ~ ""),
         significance.neg = case_when((estimate.neg < 0.1 & estimate.neg > 0.05) ~"˙", #paste0(round(estimate.neg,digits=3)," *"),
                                      (estimate.neg < 0.05 & estimate.neg > 0.001)~"*", #paste0(round(estimate.neg,digits=3)," **"),
                                      estimate.neg < 0.001 ~"**", #paste0(round(estimate.neg,digits=3)," ***"),
                                      T ~ ""),
         
         significance = case_when(((estimate.positive < 0.1 & estimate.positive > 0.05)|
                                     (estimate.neg < 0.1 & estimate.neg > 0.05))~"˙",
                                  ((estimate.positive < 0.05 & estimate.positive> 0.001)|
                                     (estimate.neg < 0.05 & estimate.neg > 0.001))~"*",
                                  ((estimate.positive < 0.001)|
                                     (estimate.neg < 0.001))~"**",
                                  T ~"")) 
inter.plot <- summary.table.for.plot.glm[["aus"]]$df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
              mutate(country="Spain")) %>%
  mutate(model="inter",
         pointshape = "Heterospecific interactions") %>%
  dplyr::filter(parameters %in% c("Focal trait","Emitter trait","Focal trait -\nEmitter trait")) %>%
  dplyr::filter(density.quantile  %in% c("low")) %>%
  right_join( mass.distribution.df, by=c("country", "parameters","trait")) %>%
  mutate(trait.names=dummy.names[trait]) %>%
  mutate(y_numb =case_when(country=="Spain" ~0.4, T~0),
         y_trait=((as.numeric(trait))+y_numb)) %>%
  mutate(rect.color=case_when((trait=="Flower width" & parameters=="Focal trait")~ "Consistent",
                              (trait=="Flower width" & parameters=="Emitter trait")~ "Consistent",
                              (trait=="Stem height" & parameters=="Focal trait")~ "Consistent",
                              (trait=="SLA" & !parameters=="Focal trait") ~ "Consistent",
                              (trait=="C13 water use efficiency" & parameters=="Focal trait") ~ "Opposite",
                              (trait=="C13 water use efficiency" & parameters=="Emitter trait") ~ "Consistent",
                              (trait=="Root mass density" & parameters=="Focal trait") ~ "Opposite",
                              T~"")) %>%
  ggplot(aes(y=country,
             x=estimate,
             fill=stat(x))) + 
  geom_density_ridges_gradient(scale=c(0.8),
                               rel_min_height = 0.005) +
  scale_fill_gradientn("Effect direction",
                       colors = c(wes_palette("Zissou1",2,type = "continuous")[2],
                                  "white",wes_palette("Zissou1",2,type = "continuous")[1]),
                       limits=c(-0.07,0.07),
                       breaks=c(-0.07,0,0.07),
                       labels=c("Negatively correlated with facilitation","Neutral",
                                "Positively correlated with facilitation"))+
  facet_grid(factor(addline_format(trait.names)) ~ factor(parameters, 
                                                          c("Focal trait","Emitter trait","Focal trait -\nEmitter trait"))) +
  geom_vline(xintercept=0) + 
  scale_x_continuous(breaks=c(-0.06,0,0.06)) +
  coord_cartesian(xlim=c(-0.07,0.07)) +
  geom_text(data=mass.distribution.df,
            aes(x=0.04,#estimate.median,
                y=country,
                label=significance),
            fontface ="bold",
            position=position_nudge(y= .4,
                                    x= 0.015),
            colour="black", 
            size=10) +
  labs(y="",
       x="Per capita effect size on heterospecific interactions") +
  guides(fill= guide_legend(title.position = "top",
                            nrow=3)) +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16,angle=0),
        strip.background = element_rect(fill="grey95",color="white",linewidth=4),
        panel.spacing = unit(3, "lines"),
        legend.title =element_text(size=20),
        legend.text =element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=90),
        panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
        panel.border = element_rect( color = NA,
                                     fill=NA),
        panel.grid.major.x = element_blank(), #element_line(colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0, "lines"),
        plot.margin=unit(c(1.5,1,0.5,1),"cm"))

inter.plot
  #figures/supp/All_TraitEffectSize_distribution
mass.distribution.intra.df <- summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
              mutate(country="Spain")) %>%
  mutate(model="intra",
         pointshape = "Conspecific interactions") %>%
  dplyr::filter(parameters %in% c("Focal trait")) %>%
  dplyr::filter(density.quantile  %in% c("low")) %>%
  group_by(country, parameters,trait) %>%
  summarise(estimate.positive = length(estimate[estimate>0])/length(estimate),
            estimate.neg = length(estimate[estimate<0])/length(estimate),
            estimate.median = median(estimate)) %>%
  #dplyr::filter(trait %in% trait.to.keep ) %>%
  ungroup() %>% 
  mutate(trait.names=dummy.names[trait]) %>%
  mutate(significance.pos = case_when((estimate.positive < 0.1 & estimate.positive > 0.05)~"˙", #paste0(round(estimate.positive,digits=3)," *"),
                                      (estimate.positive < 0.05 & estimate.positive> 0.001) ~"*", #paste0(round(estimate.positive,digits=3)," **"),
                                      estimate.positive < 0.001 ~"**", #paste0(round(estimate.positive,digits=3)," ***"),
                                      T ~ ""),
         significance.neg = case_when((estimate.neg < 0.1 & estimate.neg > 0.05) ~"˙", #paste0(round(estimate.neg,digits=3)," *"),
                                      (estimate.neg < 0.05 & estimate.neg > 0.001)~"*", #paste0(round(estimate.neg,digits=3)," **"),
                                      estimate.neg < 0.001 ~"**", #paste0(round(estimate.neg,digits=3)," ***"),
                                      T ~ ""),
         
         significance = case_when(((estimate.positive < 0.1 & estimate.positive > 0.05)|
                                     (estimate.neg < 0.1 & estimate.neg > 0.05))~"˙",
                                  ((estimate.positive < 0.05 & estimate.positive> 0.001)|
                                     (estimate.neg < 0.05 & estimate.neg > 0.001))~"*",
                                  ((estimate.positive < 0.001)|
                                     (estimate.neg < 0.001))~"**",
                                  T ~"")) 

intra.plot <- summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
              mutate(country="Spain")) %>%
  mutate(model="intra",
         pointshape = "Conspecific interactions") %>%
  dplyr::filter(parameters %in% c("Focal trait")) %>%
  dplyr::filter(density.quantile  %in% c("low")) %>%
  left_join( mass.distribution.intra.df %>%
               select(country,parameters,trait,significance), 
             by=c("country", "parameters","trait")) %>%
  mutate(trait.names=dummy.names[trait]) %>%
  mutate(rect.color=case_when((trait=="Root mass density" & parameters=="Focal trait") ~ "Consistent",
                              T~"")) %>%
  ggplot(aes(y=country,
             x=estimate,
             fill=stat(x))) + 
  geom_density_ridges_gradient(scale=c(0.8),
                               rel_min_height = 0.005) +
  scale_fill_gradientn("Direction of effect",
                       colors = c(wes_palette("Zissou1",2,type = "continuous")[2],
                                  "white",wes_palette("Zissou1",2,type = "continuous")[1]),
                       limits=c(-0.02,0.02),
                       breaks=c(-0.02,0,0.02),
                       labels=c("Negatively correlated with facilitation","Neutre",
                                "Positively correlated with facilitation"))+
  facet_grid(factor(addline_format(trait.names)) ~ factor(parameters, 
                                                                                              c("Focal trait"))) +
  geom_vline(xintercept=0) + 
  scale_x_continuous(breaks=c(-0.02,0,0.02)) +
  coord_cartesian(xlim=c(-0.02,0.02)) +
  geom_text(data=mass.distribution.intra.df,
            aes(x=0.01,#estimate.median,
                y=country,
                label=significance),
            fontface ="bold",
            position=position_nudge(y= .4,
                                    x= 0.005),
            colour="black", 
            size=10) +
  labs(y="",
       x="Per capita effect size on conspecific interactions",
       shape="Country",
       color="Functional trait") +
  guides(fill= guide_legend(title.position = "top",
                            nrow=3)) +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16,angle=0),
        strip.background = element_rect(fill="grey95",color="white",linewidth=4),
        panel.spacing = unit(3, "lines"),
        legend.title =element_text(size=20),
        legend.text =element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=90),
        panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
        panel.border = element_rect( color = NA,
                                     fill=NA),
        panel.grid.major.x = element_blank(), #element_line(colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0, "lines"),
        plot.margin=unit(c(1.5,1,0.5,1),"cm"))

intra.plot
#figures/supp/All_TraitEffectSizeIntra_distribution

lambda.plot <- summary.table.for.plot.glm[["aus"]]$Lambda.trait.df.i %>%
  mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$Lambda.trait.df.i %>%
              mutate(country="Spain")) %>%
    dplyr::filter(!parameters == "intercept") %>%
    dplyr::filter(density.quantile  %in% c("low")) %>%
    mutate(trait.names=dummy.names[trait]) %>%
   ggplot(aes(y=country,
             x=estimate,
             fill=stat(x))) + 
  geom_density_ridges_gradient(scale=c(0.8),
                               rel_min_height = 0.005) +
  scale_fill_gradientn("Direction of effect",
                       colors = c(wes_palette("Zissou1",2,type = "continuous")[2],
                                  "white",wes_palette("Zissou1",2,type = "continuous")[1]),
                       limits=c(-3,3),
                       breaks=c(-3,0,3),
                       labels=c("Negatively correlated with facilitation","Neutre",
                                "Positively correlated with facilitation"))+
  facet_grid(factor(addline_format(trait.names)) ~ factor(parameters, 
                                                          c("Focal trait"))) +
  geom_vline(xintercept=0) + 
  #scale_x_continuous(breaks=c(-0.02,0,0.02)) +
  #coord_cartesian(xlim=c(-0.02,0.02)) +
  labs(y="",
       x="Per capita effect size on conspecific interactions",
       shape="Country",
       color="Functional trait") +
  guides(fill= guide_legend(title.position = "top",
                            nrow=3)) +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=16,angle=0),
        strip.background = element_rect(fill="grey95",color="white",linewidth=4),
        panel.spacing = unit(3, "lines"),
        legend.title =element_text(size=20),
        legend.text =element_text(size=18),
        axis.title=element_text(size=20),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=90),
        panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
        panel.border = element_rect( color = NA,
                                     fill=NA),
        panel.grid.major.x = element_blank(), #element_line(colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0, "lines"),
        plot.margin=unit(c(1.5,1,0.5,1),"cm"))
  
  lambda.plot
  
  #figures/supp/All_TraitEffectSizeLambda_distribution


#---- 2.3. Trait correlation----
country="aus"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] 

  pdf(paste0("figures/pca/",country,"cor.traits.pdf"))
chart.Correlation(trait.df, histogram=TRUE, pch=19) 
#figures/pca/AusCor.traits.pdf 
#figures/pca/SpainCor.traits.pdf 

  dev.off()

}
#----2.4. Trait coefficeint with intercept----
for( country in country.list){
  for(n in c("low")){
    
    Inter.trait.df.i  <- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      mutate(trait= factor(trait ,levels=names(dummy.names) )) %>%
      group_by(trait,density.quantile,rhat,num.div,rhat.dist) %>%
      summarise(median.intercept.model.1 = median(Intercept),
                CI.intercept.model.1 = paste0("[",paste(round(quantile(Intercept,c(0.1,0.9)),digits=3),
                                                        collapse = ";"),"]"),
                median.emitter.trait.model.1 = median(Emitter.trait.scaled),
                CI.emitter.trait.model.1 = paste0("[",paste(round(quantile(Emitter.trait.scaled,c(0.1,0.9)),digits=3),
                                                            collapse = ";"),"]"),
                median.receiver.trait.model.1 = median(Receiver.trait.scaled),
                CI.receiver.trait.model.1 = paste0("[", paste(round(quantile(Receiver.trait.scaled,c(0.1,0.9)),digits=3),
                                                              collapse = ";"),"]"),
                median.intercept.model.2= median(Intercept.dist),
                CI.intercept.model.2 = paste0("[",paste(round(quantile(Intercept.dist,c(0.1,0.9)),digits=3),
                                                        collapse = ";"),"]"),
                median.dist.trait.model.2= median(Scaled.trait.dist),
                CI.dist.trait.model.2 =paste0("[",paste(round(quantile(Scaled.trait.dist,c(0.1,0.9)),digits=3),
                                                        collapse = ";"),"]"))%>%
      mutate(interactions="Heterospecific")
    
    Intra.trait.df.i <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      mutate(trait= factor(trait ,levels=names(dummy.names) )) %>%
      group_by(trait,density.quantile,rhat,num.div) %>%
      summarise(median.intercept = median(Intercept),
                CI.intercept = paste0("[",paste(round(quantile(Intercept,c(0.1,0.9)),digits=3),
                                                collapse = ";"),"]"),
                median.receiver.trait = median(Receiver.trait.scaled),
                CI.receiver.trait = paste0("[",paste(round(quantile(Receiver.trait.scaled,c(0.1,0.9)),digits=3),
                                                     collapse = ";"),"]"),) %>%
      mutate(interactions="Conspecific")
    
    write.csv(Inter.trait.df.i,
              file=paste0("results/Inter.trait.coeff.",country,".csv"))
    write.csv(Intra.trait.df.i,
              file=paste0("results/Intra.trait.coeff.",country,".csv"))
    
  }
}

#---- 2.5. Species observations fecundity and in neighbours----
# for median of individuals 
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  competition_df <- get(paste0("clean.data.",country))[[paste0("competition_",country)]] 
  if(country=="aus"){
    competition_df <- competition_df %>%
      gather(any_of(Code.focal.list),
             key="species",value="individuals") %>%
      mutate(individuals = (individuals/scale)*15) %>%
      dplyr::rename("com_id"="plot") %>%
      aggregate(individuals ~ year + com_id + species + focal, sum)  %>%
      spread(species, individuals) 
    
    competition_df %>%
      group_by(focal) %>%
      summarise()
    
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3. Other trait graphs to explore----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#---- 3.1. Make figure for preso----

dummy.col <- c("C13 water use efficiency"="#539E59FF","SLA"="#F5EC54FF",
                 "Root mass density"="#B25D91FF",#"Root length"="#CB87B4FF" ,
                 "SRL"="#F28E2BFF")
dummy.names <- c("C13 water use efficiency"="Water use efficiency",
                 "SLA"= "Specific leaf area",
               "Root mass density"="Root tissue density",#"Root length"="#CB87B4FF" ,
               "SRL"="Specific root length")
dummy.col <- c("SRL"="#4E79A7FF","SRA"="#76B7B2FF" ,"Root length"="#A4BED5FF","Root tips"="#512DA8FF" ,
               "Root diameter"="#F28E2BFF" , "Root mass density"="#ED645AFF",
               "Flower width"= "#FF9DA7FF" ,"Seed mass"="#B276B2FF",
               "C13 water use efficiency"="#9C755FFF",
               "Leaf C to N ratio"= "#BCBD22FF" ,
               "Leaf area index"="#D4E157FF" ,"Canopy shape"="#72874EFF",
               "SLA"="#59A14FFF","Stem height"="#FED789FF",
               "intercept"="grey80")
dummy.names <- c("SRL"="Specific root length",
                 "SRA"="Specific root area" ,
                 "Root length"= "Root length","Root tips"="Root tips",
               "Root diameter"="Root diameter" , "Root mass density"="Root tissue density",
               "Flower width"= "Flower width" ,"Seed mass"="Seed mass",
               "C13 water use efficiency"="Water use efficiency",
               "Leaf C to N ratio"= "Nitrogen use efficiency" ,
               "Leaf area index"="Leaf area index" ,"Canopy shape"="Canopy shape",
               "SLA"="Specific leaf area","Stem height"="Stem height")

library(paletteer)
paletteer_d("LaCroixColoR::Berry")
paletteer_d("beyonce::X22")
trait.to.keep <- names(dummy.col)
legend.plot <- ggplot(data=Cool.theory.trait.df[["aus"]]$Intra.trait.df %>%
                        bind_rows(Cool.theory.trait.df[["spain"]]$Intra.trait.df) %>%
                        dplyr::filter(outlier=="included") %>%
                        #mutate(pointshape = rep(c("Trait effect on interaction",
                        #                         "Intercept of interaction strength"),
                        #                       each=120)) %>%
                        dplyr::filter(!parameters =="intercept") %>%
                        dplyr::filter(trait %in% trait.to.keep) %>%
                        #bind_rows(intercept.df) %>%
                        filter(density.quantile  %in% c("high","intercept")) %>%
                        mutate( density.quantile = factor(density.quantile ,
                                                          levels=c("intercept","high"))) %>%
                        mutate(trait=factor(trait,
                                            levels=rev(names( dummy.col)))),
                      aes(y=Estimate,
                          x=parameters,
                          color=trait#,
                          #shape = pointshape
                      ))+
  #shape=density.quantile )) + 
  geom_blank() + 
  geom_pointrange(aes(xmin=Q2.5,
                      xmax=Q97.5),
                  size=2,alpha=0.8) +
  #scale_shape_manual("",values= 17 )  +
  scale_color_manual(values= dummy.col,
                     labels=dummy.names)  +
  guides(shape= guide_legend(nrow=2,byrow=TRUE,
                            override.aes = list(alpha = 1,
                                             color=c("grey80"))),
   color= guide_legend(nrow=2,byrow=TRUE)) +
  labs(x="estimate",
       y=paste0(""),
       shape="Model estimates",
       color="Functional trait") +
  theme_few() +
  theme(legend.position="bottom",
        legend.key.size = unit(1, 'cm'),
        legend.title.position = "top",
        legend.title =element_text(size=20),
        legend.text =element_text(size=20),
        axis.text = element_text(size=18))
legend.plot


inter.plot <- summary.table.for.plot.glm[["spain"]]$df.i %>%
  mutate(country="Spain") %>%
  bind_rows(summary.table.for.plot.glm[["aus"]]$df.i %>%
              mutate(country="Australia")) %>%
  dplyr::filter(!parameters == "intercept") %>%
  dplyr::filter( outlier== "removed") %>%
  dplyr::filter(density.quantile  %in% c("low")) %>%
  dplyr::filter(trait %in% trait.to.keep )%>%
  mutate(Q2.5 = case_when(Q2.5 < -0.04 ~ -0.04,
                          T ~ Q2.5),
         Q97.5 = case_when(Q97.5 > 0.04 ~ 0.04,
                           T ~ Q97.5)) %>%
  mutate(y_numb= case_when(parameters =="Focal trait" ~ 3,
                           parameters =="Emitter trait"~2,
                           parameters =="Focal trait -\nEmitter trait" ~1)) %>%
  mutate(y_numb = case_when(model=="intra" ~ 1,
                            T~ y_numb)) %>%
  mutate(trait=factor(trait, 
                      rev(names(dummy.col)[names(dummy.col) %in% trait.to.keep]))) %>%
  mutate(y_trait=((as.numeric(trait)-2)*0.14 +y_numb)) %>%
  mutate(density.quantile.text = case_when(density.quantile =="low" ~"Low Neighbor's Density")) %>%
  ggplot(aes(y=y_trait,
             x=Estimate,
             color=as.factor(trait),
             group=as.factor(trait))) + 
  geom_linerange(aes(xmin=Q2.5,
                     xmax=Q97.5),
                 size=1,alpha=0.7) +
  geom_pointrange(aes(xmin=Q2.5,
                      xmax=Q97.5),
                  size=2,alpha=0.7) +
  geom_text(aes(y=y_trait, x=Estimate,
                group=as.factor(trait),label=signif),
            size=10,color="black") +
  scale_color_manual(values=dummy.col)+
  scale_y_continuous(breaks=c(1:3),
                     labels=rev(c("Focal trait",
                                  "Neighbor trait",
                                  "Focal trait -\nNeighbor trait"))) +
  annotate(geom = "text",label="Interaction \nmore positive",
           x=0.02,y=0.4,size=6,angle=0)+
  geom_segment(data=NULL,x=0, xend = +0.03 , 
               y=0.01, yend = 0.01, size=1,
               arrow = arrow(length = unit(0.6,"cm")),
               color="black")+
  annotate(geom = "text",label="Interaction \nmore negative",
           x=-0.02,y=0.4,size=6,angle=0)+
  geom_segment(data=NULL,x=0, xend = -0.03 , 
               y=0.01, yend = 0.01, size=1,
               arrow = arrow(length = unit(0.6,"cm")),
               color="black")+
  theme_bw() +
  coord_cartesian(xlim=c(-0.04,0.04),
                  ylim=c(0,3.6),
                  clip = "off",expand=F)+
  facet_wrap(.~ country,ncol=2,scale="free") +
  geom_vline(xintercept=0) + 
  labs(x="Effect size",
       y=paste0(""),
       shape="Model estimates",
       color="Functional trait") +
  guides(color = guide_legend(title.position = "top",
                              nrow=2),
         shape = guide_legend(title.position = "top",
                              nrow=2)) +
  theme(legend.position="none",
        strip.text = element_text(size=20),
        strip.background = element_rect(fill="grey98",color="black"),
        panel.spacing = unit(3, "lines"),
        legend.title =element_text(size=20),
        legend.text =element_text(size=18),
        axis.title=element_text(size=20),
        axis.text = element_text(size=18),
        panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
        panel.border = element_rect( color = "grey60"),
        panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,0,1),"cm"))

inter.plot

#-----3.2. Make plot trait and coefficient ----
Trait.coeff.plotlist <- list()
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
density.quantile.name <- c("low","high")
country="aus"
dummy.col <- c("SRL"="#4E79A7FF","SRA"="#76B7B2FF" ,"Root length"="#A4BED5FF","Root tips"="#512DA8FF" ,
               "Root diameter"="#F28E2BFF" , "Root mass density"="#ED645AFF",
               "Flower width"= "#FF9DA7FF" ,"Seed mass"="#B276B2FF",
               "C13 water use efficiency"="#9C755FFF",
               "Leaf C to N ratio"= "#BCBD22FF" ,
               "Leaf area index"="#D4E157FF" ,"Canopy shape"="#72874EFF",
               "SLA"="#59A14FFF","Stem height"="#FED789FF",
               "intercept"="grey80")

dummy.names <- c("SRL"="Specific root length",
                 "SRA"="Specific root area" ,
                 "Root length"= "Root length","Root tips"="Root tips",
                 "Root diameter"="Root diameter" , "Root mass density"="Root tissue density",
                 "Flower width"= "Flower width" ,"Seed mass"="Seed mass",
                 "C13 water use efficiency"="Water use efficiency",
                 "Leaf C to N ratio"= "Nitrogen use efficiency" ,
                 "Leaf area index"="Leaf area index" ,"Canopy shape"="Canopy shape",
                 "SLA"="Specific leaf area","Stem height"="Stem height")

country="spain"
n="low"
for( country in country.list){
  for(n in c("low")){
    
    Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
    trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] %>%
      dplyr::select(-any_of(c("Canopy shape","Mean fecundity","Leaf area index")))
    specific.trait.dist <- Cool.theory.trait.df[[country]]$trait.dist.df
    # Make data frame with trait and INTRA specific interactions
    taxize_dist.df <- taxize_dist.list[[country]]
    
    trait.coeff.df <- Theoretical.Int.list[[country]] %>%
      dplyr::filter(!neigh==focal & 
                      density.quantile  ==n) %>%
      dplyr::select(neigh,focal,country,lambda,density.quantile, theoretical.effect) %>%
      left_join(specific.trait.dist, 
                relationship ="many-to-many",
                by=c("neigh","focal")) %>%
      left_join(  taxize_dist.df %>%
                    dplyr::select(phylo.dist,focal,neigh),
                  by=c("neigh","focal"))
    Trait.coeff.plotlist[[paste0(country,".df")]]<-  trait.coeff.df
  }
}

for(traitagent in c("receiver.trait",
                     "emitter.trait",
                     "trait.dist")){
  print(traitagent)
  trait.coeff.plot <- bind_rows(Trait.coeff.plotlist[["aus.df"]],
              Trait.coeff.plotlist[["spain.df"]]) %>%
      mutate(sigmoid.exp = (exp(theoretical.effect)-1)) %>%
      ggplot(aes_string(x=traitagent,y="sigmoid.exp" ,
                 color="sigmoid.exp")) +
      geom_point() +
      facet_grid(country ~addline_format(trait),
                 scale="free_x") +
      scale_color_gradientn(colours = rev(wes_palette("Zissou1", 
                                                      101, 
                                                      type = "continuous")))+
    
       labs(x=traitagent,
            y="percentage of intrinsic lambda modifid",
            color="percentage of effect (y-axis)")+
      theme_few() +
      theme(legend.position ="bottom",
            strip.text.y = element_text(angle=0),
            axis.text.x = element_text(angle=90),
            plot.background = element_rect(color="white",fill="white"),
            panel.background = element_rect(color="white",fill="white"),
            panel.grid.major.x = element_line(color="gray90"))
  
  plot(trait.coeff.plot)
  Trait.coeff.plotlist[[ paste0(traitagent,".plot")]] <- trait.coeff.plot
  
}

ggarrange(Trait.coeff.plotlist[["receiver.trait.plot"]],
          Trait.coeff.plotlist[["emitter.trait.plot"]],
          Trait.coeff.plotlist[["trait.dist.plot"]],
          ncol=1,
          common.legend = T,
          legend="bottom")
#figures/supp/Traitvalues_sigmoid.pdf

for(traitagent in c("receiver.trait.scaled",
                    "emitter.trait.scaled",
                    "scaled.trait.dist")){
  print(traitagent)
  trait.coeff.plot <- bind_rows(Trait.coeff.plotlist[["aus.df"]],
                                Trait.coeff.plotlist[["spain.df"]]) %>%
    mutate(sigmoid.exp = (exp(theoretical.effect)-1)) %>%
    ggplot(aes_string(x=traitagent,y="sigmoid.exp" ,
                      color="sigmoid.exp")) +
    geom_point() +
    facet_grid(country ~addline_format(trait),
               scale="free_x") +
    scale_color_gradientn(colours = rev(wes_palette("Zissou1", 
                                                    101, 
                                                    type = "continuous")))+
    
    labs(x=traitagent,
         y="percentage of intrinsic lambda modifid",
         color="percentage of effect (y-axis)")+
    theme_few() +
    theme(legend.position ="bottom",
          strip.text.y = element_text(angle=0),
          axis.text.x = element_text(angle=90),
          plot.background = element_rect(color="white",fill="white"),
          panel.background = element_rect(color="white",fill="white"),
          panel.grid.major.x = element_line(color="gray90"))
  
  plot(trait.coeff.plot)
  Trait.coeff.plotlist[[ paste0(traitagent,".plot")]] <- trait.coeff.plot
  
}

ggarrange(Trait.coeff.plotlist[["receiver.trait.scaled.plot"]],
          Trait.coeff.plotlist[["emitter.trait.scaled.plot"]],
          Trait.coeff.plotlist[["scaled.trait.dist.plot"]],
          ncol=1,
          common.legend = T,
          legend="bottom")
#figures/supp/Trait.scaled.values_sigmoid.pdf
