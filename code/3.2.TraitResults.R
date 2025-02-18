
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
#setwd("/home/lbuche/Eco_Bayesian/chapt3")
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
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
Param.sigm.df <- list()
Param.sigm.df.aus <- read.csv(paste0(project.dic,"results/Param.sigmoid.aus.csv.gz"))
Param.sigm.df.spain <- read.csv(paste0(project.dic,"results/Param.sigmoid.spain.csv.gz"))
Param.sigm.df$spain <- Param.sigm.df.spain
Param.sigm.df$aus<- Param.sigm.df.aus
save(Param.sigm.df,
     file=paste0(project.dic,"results/Param.sigm.df.RData"))
load(paste0(project.dic,"results/Param.sigm.df.RData"))
# realised interactions for each year
load(paste0(project.dic,"results/Theoretical.Int.list.RData")) 


# realised interactions across time
load(paste0(project.dic,"results/Realised.Int.list.RData"))

# realised interactions for each year
load(paste0(project.dic,"results/Realised.Int.Year.list.RData")) 

load(paste0(project.dic,"results/Realised.Obs.Year.list.RData")) 
# PDSI
env_pdsi_aus <- read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) 
env_pdsi_spain <- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv")) 
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.Looking at theory interactions----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.1. Make data df ----
country = "aus" 
Cool.theory.trait.df <- list()
density.quantile.name <- c("intercept","low","medium","high")
rsq <- function (x, y) cor(x, y,use="complete.obs") ^ 2
zscore <- function (x) (x - mean(x,na.rm=T))/sd(x,na.rm=T)

for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  # make trait data frame
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] %>%
    dplyr::select(-c("Mean fecundity"))
  library(vegan)
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
  #trait.i = "SLA"
  for( trait.i in names(trait.df)){
    for(n in density.quantile.name){
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
      
      glm.lambda.trait.i <- glmmTMB(lambda ~ 1 + receiver.trait.scaled + (1|focal) ,
                                trait.lambda.df.i,
                                family="gaussian")
      glm.lambda.trait.outlier.i <- glmmTMB(lambda ~ 1 + receiver.trait.scaled + (1|focal) ,
                                    trait.lambda.df.i %>%
                                      dplyr::filter(receiver.trait.scaled>-3 &
                                                      receiver.trait.scaled<3)%>%
                                      mutate(receiver.trait.scaled = zscore(receiver.trait)),
                                    family="gaussian")
      
      
      Lambda.trait.df.i <- as.data.frame(confint(glm.lambda.trait.i)) %>%
        mutate(outlier= "included")%>%
        rownames_to_column("parameters") %>%
        bind_rows(as.data.frame(confint(glm.lambda.trait.outlier.i))%>%
                    rownames_to_column("parameters")%>%
                    mutate(outlier= "removed")) %>%
        rename("Q2.5"="2.5 %",
               "Q97.5"="97.5 %") %>%
        mutate(model="trait",
               trait=trait.i,
               density.quantile=n,
               signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                                (Q2.5>0 & Q97.5>0 )~"*",
                                T~""))
      
      Lambda.trait.df <- bind_rows( Lambda.trait.df, Lambda.trait.df.i)
      
    }
  }
  # Make data frame with trait and INTRA specific interactions
  trait.value.intra.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(neigh==focal)%>%
    left_join(specific.trait.dist,
              relationship ="many-to-many")
  Intra.trait.df <- NULL
  glm.intra.trait.summary <- NULL
  #trait.i = "SLA"
  #trait.i = "C13 water use efficiency"
  for( trait.i in names(trait.df)){
    for(n in density.quantile.name){
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
      
      glm.intra.trait.i <- glmmTMB(theoretical.effect ~ 1 + receiver.trait.scaled  + (1|focal),
                                    trait.intra.df.i,
                                    family="gaussian")
      #remove outlier
      glm.intra.trait.outlier.i <- glmmTMB(theoretical.effect ~ 1 + receiver.trait.scaled  + (1|focal),
                                   trait.intra.df.i %>%
                                     dplyr::filter(receiver.trait.scaled>-3 &
                                                      receiver.trait.scaled<3)%>%
                                     mutate(receiver.trait.scaled = zscore(receiver.trait)),
                                   family="gaussian")
      summary(glm.intra.trait.i)
      summary(glm.intra.trait.outlier.i)
      Intra.trait.df.i <- as.data.frame(confint(glm.intra.trait.i)) %>%
        mutate(outlier= "included")%>%
        rownames_to_column("parameters") %>%
        bind_rows(as.data.frame(confint(glm.intra.trait.outlier.i))%>%
                    rownames_to_column("parameters")%>%
                    mutate(outlier= "removed")) %>%
        rename("Q2.5"="2.5 %",
               "Q97.5"="97.5 %") %>%
        mutate(model="trait",
               trait=trait.i,
               density.quantile=n,
               signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                                (Q2.5>0 & Q97.5>0 )~"*",
                                T~""))
      
      
      Intra.trait.df <- bind_rows(Intra.trait.df,Intra.trait.df.i)
      
    }
  }
 
  # Make data frame with trait and inter specific interactions
  trait.dist.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    left_join(as.data.frame(specific.trait.dist),
              relationship ="many-to-many")
  Inter.trait.df <- NULL
  glm.inter.trait.summary<- NULL
  # trait.i ="SRL"
  # n ="intercept"
  for( trait.i in names(trait.df)){
    for(n in density.quantile.name){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(density.quantile==n) %>%
        dplyr::select(neigh,focal,density.quantile,trait,
                      theoretical.effect,
                      emitter.trait,emitter.trait.scaled,emitter.trait.log.scaled,
                      receiver.trait,receiver.trait.scaled,receiver.trait.log.scaled,
                      trait.dist,scaled.trait.dist,trait.ratio,scaled.trait.ratio) %>%
        mutate(focal=as.factor(focal),
               neigh=as.factor(neigh))
      
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
      
      glm.inter.trait.i.dist <- glmmTMB(theoretical.effect ~  scaled.trait.dist  + (1|focal) + (1|neigh), 
                               trait.dist.df.i,
                               family="gaussian")
      glm.inter.trait.i <- glmmTMB(theoretical.effect ~  receiver.trait.scaled  + emitter.trait.scaled  + (1|focal) + (1|neigh), 
                                   trait.dist.df.i,
                                   family="gaussian")
  
      glm.inter.trait.outlier.i<- glmmTMB(theoretical.effect ~  receiver.trait.scaled  + emitter.trait.scaled  + (1|focal) + (1|neigh), 
                                                trait.dist.df.i %>%
                                             dplyr::filter(receiver.trait.scaled > -3 &
                                                             receiver.trait.scaled < 3 &
                                                             emitter.trait.scaled >-3 &
                                                             emitter.trait.scaled < 3)%>%
                                             mutate(receiver.trait.scaled = zscore(receiver.trait),
                                                    emitter.trait.scaled = zscore(emitter.trait)),
                                           family="gaussian")
      #summary(glm.inter.trait.i)
      #summary(  glm.inter.trait.outlier.i )
      glm.inter.trait.outlier.i.dist  <-  glmmTMB(theoretical.effect ~  scaled.trait.dist  + (1|focal) + (1|neigh), 
                                            trait.dist.df.i %>%
                                             dplyr::filter(scaled.trait.dist>-3 &
                                                             scaled.trait.dist<3)%>%
                                             mutate(scaled.trait.dist= zscore(trait.dist)),
                                           family="gaussian")
      
      Inter.trait.df.i <- as.data.frame(confint(glm.inter.trait.i.dist)) %>%
        rownames_to_column("parameters") %>%
        mutate(model = "distance")%>%
      bind_rows(as.data.frame(confint(glm.inter.trait.i))%>%
                    mutate(model = "trait")%>%
                  rownames_to_column("parameters")) %>%
        rename("Q2.5"="2.5 %",
               "Q97.5"="97.5 %") %>%
        mutate(outlier= "included",
               trait=trait.i,
               density.quantile=n,
               signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                                (Q2.5>0 & Q97.5>0 )~"*",
                                T~""))
      
      Inter.trait.df.outlier.i <- as.data.frame(confint(glm.inter.trait.outlier.i.dist)) %>%
        rownames_to_column("parameters") %>%
        mutate(model = "distance")%>%
        bind_rows(as.data.frame(confint(glm.inter.trait.outlier.i))%>%
                    mutate(model = "trait")%>%
                    rownames_to_column("parameters")) %>%
        rename("Q2.5"="2.5 %",
               "Q97.5"="97.5 %") %>%
        mutate(outlier= "removed",
               trait=trait.i,
               density.quantile=n,
               signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                                (Q2.5>0 & Q97.5>0 )~"*",
                                T~""))
      
      Inter.trait.df <- bind_rows(Inter.trait.df, bind_rows(Inter.trait.df.i,
                                                            Inter.trait.df.outlier.i))
      
    }
  }
  
  Cool.theory.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    trait.value.lambda.df =  trait.value.lambda.df, 
    trait.value.intra.df=trait.value.intra.df,
    Inter.trait.df=Inter.trait.df,
    Lambda.trait.df=Lambda.trait.df,
    Intra.trait.df =Intra.trait.df)
}
view(Cool.theory.trait.df[[country]]$Inter.trait.df)
trait.dist.df <- Cool.theory.trait.df[["aus"]]$trait.dist.df
trait.dist.df %>%
  dplyr::filter(trait %in% c("SLA","C13 water use efficiency")) %>%
  dplyr::select(focal,trait,receiver.trait) %>%
  spread(trait,receiver.trait)
ggplot()+
  geom_point(aes(x=trait.dist.df$receiver.trait[which(trait.dist.df$trait=="SLA")],
                 y=trait.dist.df$receiver.trait[which(trait.dist.df$trait=="C13 water use efficiency")]))

#---- 1.2. Make detailed graphs ----
Cool.detailed.theory.trait.plotlist <- list()
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
density.quantile.name <- c("intercept","low","medium","high")
country="aus"
for( country in country.list){
  for(n in density.quantile.name){
    dummy.col <- c("SRL"="#4E79A7FF","SRA"="#76B7B2FF" ,"Root length"="#A4BED5FF","Root tips"="#512DA8FF" ,
                   "Root diameter"="#F28E2BFF" , "Root mass density"="#ED645AFF",
                   "Flower width"= "#FF9DA7FF" ,"Seed mass"="#B276B2FF",
                   "C13 water use efficiency"="#9C755FFF",
                   "Leaf C to N ratio"= "#BCBD22FF" ,
                   "Leaf area index"="#D4E157FF" ,"Canopy shape"="#72874EFF",
                   "SLA"="#59A14FFF","Stem height"="#FED789FF",
                   "intercept"="grey80")
    if(country=="aus"){
      trait.levels <- c("SRL","Root tips","Root mass density","Root length",
                        "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                        "Canopy shape","Stem height","SLA")
      
    }
    if(country=="spain"){
      trait.levels <- c("SRL","Root diameter","Root mass density","SRA",
                                              "Mean fecundity","C13 water use efficiency",
                        "Leaf C to N ratio","Leaf area index",
                                              "Canopy shape","Stem height","SLA")
    }

    
    Inter.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      dplyr::filter(outlier== "included") %>%
      dplyr::filter(!signif  =="" | parameters == "(Intercept)") %>%
      dplyr::select(parameters,Estimate, trait,model) %>%
      dplyr::filter(parameters  %in% c("emitter.trait.scaled",
                                       "receiver.trait.scaled",
                                       "scaled.trait.dist","(Intercept)")) %>%
      spread(parameters,Estimate) %>%
      rename("Intercept"="(Intercept)") %>%
      gather(any_of(c("scaled.trait.dist",
               "emitter.trait.scaled","receiver.trait.scaled")),
             key="trait.param",value="trait.coeff")  %>%
      dplyr::filter(!is.na(trait.coeff))
          

    Intra.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      dplyr::filter( outlier== "included") %>%
     dplyr::filter(!signif  =="" | parameters %in% c("(Intercept)")) %>%
      dplyr::select(parameters,Estimate, trait) %>%
      spread(parameters,Estimate) %>% 
      dplyr::filter(!is.na(receiver.trait.scaled)) %>%
      rename("trait.coeff"="receiver.trait.scaled",
             "Intercept"="(Intercept)") %>%
      mutate(trait= factor(trait ,levels=trait.levels )) 
    
    
    Lambda.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
      dplyr::filter( outlier== "included") %>%
       dplyr::filter(!signif  =="" | parameters %in% c("(Intercept)")) %>%
      dplyr::select(parameters,Estimate, trait) %>%
      spread(parameters,Estimate) %>% 
      dplyr::filter(!is.na(receiver.trait.scaled)) %>%
      rename("trait.coeff"="receiver.trait.scaled",
             "Intercept"="(Intercept)") %>%
      mutate(trait= factor(trait ,levels=trait.levels )) 
    
    Inter.trait.df.i <- Cool.theory.trait.df[[country]]$trait.dist.df %>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(theoretical.effect,receiver.trait.scaled,trait,emitter.trait.scaled,scaled.trait.dist) %>%
      gather(receiver.trait.scaled,emitter.trait.scaled,scaled.trait.dist,
             key="trait.param", value="trait.value") %>%
      rename("raw.value"="theoretical.effect") %>%
      mutate(parameter="INTER") %>% 
      left_join(Inter.trait.df.long.i) %>%
      dplyr::filter(!is.na(trait.coeff))
    
    Intra.trait.df.i <- Cool.theory.trait.df[[country]]$trait.value.intra.df%>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(theoretical.effect,receiver.trait.scaled,trait) %>%
      rename("raw.value"="theoretical.effect",
             "trait.value" ="receiver.trait.scaled")%>%
      mutate(parameter="INTRA",
             trait.param="receiver.trait.scaled")%>% 
      left_join(Intra.trait.df.long.i) %>%
      dplyr::filter(!is.na(trait.coeff)) 
    
    Lambda.trait.df.i <- Cool.theory.trait.df[[country]]$trait.value.lambda.df %>%
      dplyr::filter( trait %in% levels(as.factor(Lambda.trait.df.long.i$trait))) %>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(lambda,receiver.trait.scaled,trait) %>%
      rename("raw.value"="lambda",
             "trait.value" ="receiver.trait.scaled")%>%
      mutate(parameter="lambda",
             trait.param="receiver.trait.scaled")%>% 
      left_join(Lambda.trait.df.long.i) %>%
      dplyr::filter(!is.na(trait.coeff))
   
    
    pdf(file = paste0("figures/GLM/GLM_details_",country,"_",n,"_density.pdf"),
        onefile = TRUE)
    
    Cool.detailed.theory.trait.plotlist[[paste0(country,"_Inter")]] <- ggplot() +
      geom_point(data=Inter.trait.df.i,
                 aes(y=raw.value,
                     x=trait.value,
                     color=raw.value),
                 shape=16,
                 size=1,alpha=1) + 
      geom_abline(data=Inter.trait.df.long.i,
                  aes(slope=trait.coeff, intercept=Intercept),
                  size=1.5) +
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
      geom_point(data=Intra.trait.df.i,
                 aes(y=raw.value,
                     x=trait.value,
                     color=raw.value),
                 shape=16,
                 size=3,alpha=1) + 
      geom_abline(data=Intra.trait.df.long.i,
                  aes(slope=trait.coeff, intercept=Intercept),
                  size=1.5) +
      facet_grid(addline_format(trait) ~ .,scale="free") +
      geom_hline(yintercept=0,color="grey12",linetype="dashed") +
      labs(x="Trait value of focal species",
           y=paste0("Intraspecific interaction with ", n ,"density of individuals"))+
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
      geom_point(data=Lambda.trait.df.i,
                 aes(y=raw.value,
                     x=trait.value,
                     color=raw.value),
                 shape=16,
                 size=3,alpha=1) + 
      geom_abline(data=Lambda.trait.df.long.i,
                  aes(slope=trait.coeff, intercept=Intercept),
                  size=1.5) +
      facet_grid(addline_format(trait) ~ .,scale="free") +
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
    
  
   legend.plot <- ggplot(data=Cool.theory.trait.df[[country]]$Intra.trait.df %>%
                           dplyr::filter( outlier== "included") %>%
                           mutate(trait=factor(trait,
                                               levels=trait.levels)),
                         aes(y=Estimate,
                             x=parameters,
                             color=trait)) + geom_blank() + geom_line(size=3) +
     scale_color_manual(values= dummy.col )  +
     theme_few() +
     theme(legend.position="bottom",
           legend.key.size = unit(1, 'cm'),
           legend.title.position = "top",
           legend.title =element_text(size=20),
           legend.text =element_text(size=16))
   legend.plot
     
   
   
    Cool.detailed.theory.trait.plotlist[[paste(country,"_",n)]] <- plot_grid(ggplot() +
                geom_point(data=Inter.trait.df.i %>%
                             mutate(trait.param.label = case_when(trait.param == "emitter.trait.scaled" ~ "Trait value of emitter",
                                                                  trait.param == "receiver.trait.scaled" ~ "Trait value of focal/receiver",
                                                                  trait.param == "scaled.trait.dist" ~ "Trait value of focal - Trait value of emitter")),
                           aes(y=raw.value,
                               x=trait.value,
                               color=trait),
                           shape=16,
                           size=3,alpha=0.2)+
      geom_abline(data=Inter.trait.df.long.i%>%
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
      common.legend = T,legend="none"),
      ggpubr::get_legend(legend.plot),
     ncol = 1,
     rel_heights =c(1,0.6,0.2),
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
Cool.detailed.theory.trait.plotlist[["spain_intercept"]] # figures/Theory.int.intercept.spain.pdf
Cool.detailed.theory.trait.plotlist[[paste0("aus")]]

#---- 1.3. Make summary of glm plots ----
#---- 1.3.1 Make summary tables to do plots----
summary.table.for.plot.glm <- list()
for( country in country.list){
  if(country=="aus"){
    Inter.trait.df <- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA")))
    trait.levels <- c("SRL","Root tips","Root mass density","Root length",
                      "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                      "Canopy shape","Stem height","SLA")
  }
  if(country=="spain"){
    Inter.trait.df<- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
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
                 "intercept"="grey80")
  
  
    df.i <- Inter.trait.df %>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
    #dplyr::filter(!signif  ==""|parameters=="(Intercept)") %>%
    mutate(parameters  = case_when(parameters =="receiver.trait.scaled" ~ "Focal trait",
                                   parameters =="emitter.trait.scaled" ~ "Emitter trait",
                                   parameters =="scaled.trait.dist" ~ "Focal trait -\nEmitter trait",
                                   parameters =="(Intercept)" ~ "intercept",
                                   T ~ parameters)) %>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))
  
  Intra.trait.df.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
    #dplyr::filter(!signif  =="" ) %>%
    mutate(parameters  = case_when( parameters =="(Intercept)" ~ "intercept",
                                    parameters =="receiver.trait.scaled" ~ "Focal trait",
                                    T ~ parameters))%>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))
  
  Lambda.trait.df.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
    # dplyr::filter(!signif  ==""|parameters=="(Intercept)") %>%
    mutate(parameters  = case_when(parameters =="receiver.trait.scaled" ~ "Focal trait",
                                   parameters =="(Intercept)" ~ "intercept",
                                   T ~ parameters))%>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))
  
  intercept.df <- df.i %>%
    dplyr::mutate(trait = case_when(parameters =="intercept" ~ "intercept",
                                    T~trait)) %>%
    dplyr::mutate(parameters = case_when(parameters =="intercept" ~ NA,
                                         T~parameters)) %>%
    fill(parameters,.direction = c("up")) %>%
    filter(trait =="intercept") %>%
    group_by(parameters,trait,density.quantile) %>%
    summarise(Q2.5=min(Q2.5),
              Q97.5=max(Q97.5),
              Estimate =median(Estimate))%>%
    ungroup()
  

  summary.table.for.plot.glm[[country]] <- list(intercept.df =intercept.df,
                                               Lambda.trait.df.i=Lambda.trait.df.i,
                                               Intra.trait.df.i=Intra.trait.df.i,
                                               df.i=df.i)
  
}
#---- 1.3.2 Make plots for lambda, intra and inter----
Cool.glm.theory.trait.plotlist <- list()
density.quantile.name <- c("intercept","low","medium","high")
country="aus"

for( country in country.list){
  
  Inter.plot.sum <- summary.table.for.plot.glm[[country]]$df.i %>%
    dplyr::filter( outlier== "included") %>%
    dplyr::filter(!parameters =="intercept") %>%
    #bind_rows(intercept.df ) %>%
    mutate(trait=factor(trait, levels=names(dummy.col))) %>%
    ggplot(aes(y=as.factor(parameters),
               x=Estimate,
              color=as.factor(trait))) +
    geom_pointrange(aes(xmin=Q2.5,
                        xmax=Q97.5),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.6)) +
    scale_color_manual(values=dummy.col)+
    theme_bw() +
    geom_vline(xintercept=0) + 
    labs(x="estimate",
         y=paste0(""),
         color="Density of emitter") +
    facet_wrap(.~density.quantile,ncol=4,nrow=1) +
    guides(color = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="none",
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=16),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
  Inter.plot.sum
  
  intercept.df <- summary.table.for.plot.glm[[country]]$Intra.trait.df.i %>%
    dplyr::filter( outlier== "included") %>%
    dplyr::mutate(trait = case_when(parameters =="intercept" ~ "intercept",
                                    T~trait)) %>%
    dplyr::mutate(parameters = case_when(parameters =="intercept" ~ NA,
                                         T~parameters)) %>%
    fill(parameters,.direction = c("up")) %>%
    filter(trait =="intercept") %>%
    group_by(parameters,trait,density.quantile) %>%
    summarise(Q2.5=min(Q2.5),
              Q97.5=max(Q97.5),
              Estimate =median(Estimate))%>%
    ungroup()
  intercept.df
  
  Intra.plot.sum <- summary.table.for.plot.glm[[country]]$Intra.trait.df.i %>%
    dplyr::filter( outlier== "included") %>%
    dplyr::filter(!parameters =="intercept") %>%
    #bind_rows(intercept.df ) %>%
    mutate(trait=factor(trait, levels=names(dummy.col))) %>%
    ggplot(aes(y=as.factor(parameters),
               x=Estimate,
               color=as.factor(trait))) +
    geom_pointrange(aes(xmin=Q2.5,
                        xmax=Q97.5),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.9)) +
    scale_color_manual(values=dummy.col)+
    theme_bw() +
    geom_vline(xintercept=0) + 
    labs(x="estimate",
         y=paste0(""),
         color="Density of emitter") +
    facet_wrap(.~density.quantile,ncol=4,nrow=1) +
    guides(color = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="none",
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=16),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
 
  intercept.df <- summary.table.for.plot.glm[[country]]$Lambda.trait.df.i %>%
    dplyr::mutate(trait = case_when(parameters =="intercept" ~ "intercept",
                                    T~trait)) %>%
    dplyr::mutate(parameters = case_when(parameters =="intercept" ~ NA,
                                         T~parameters)) %>%
    fill(parameters,.direction = c("up")) %>%
    filter(trait =="intercept") %>%
    group_by(parameters,trait,density.quantile) %>%
    summarise(Q2.5=min(Q2.5),
              Q97.5=max(Q97.5),
              Estimate =median(Estimate))%>%
    ungroup()

  
   Lambda.plot.sum <- summary.table.for.plot.glm[[country]]$Lambda.trait.df.i %>%
     dplyr::filter( outlier== "included") %>%
     dplyr::filter(!parameters =="intercept") %>%
     bind_rows(intercept.df ) %>%
     mutate(trait=factor(trait, levels=names(dummy.col))) %>%
    ggplot(aes(y=as.factor(parameters),
               x=Estimate,
               color=as.factor(trait))) +
    geom_pointrange(aes(xmin=Q2.5,
                        xmax=Q97.5),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.9)) +
    scale_color_manual(values=dummy.col)+
    theme_bw() +
    labs(x="estimate",
         y=paste0(""),
         color="Density of emitter") +
   #facet_wrap(.~density.quantile,ncol=4,nrow=1) +
    guides(color = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="none",
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=16),
          axis.text = element_text(size=14),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
   
   legend.plot <- ggplot(data=Cool.theory.trait.df[["aus"]]$Intra.trait.df %>%
                           bind_rows(Cool.theory.trait.df[["spain"]]$Intra.trait.df) %>%
                           dplyr::filter(!parameters =="intercept") %>%
                           bind_rows(intercept.df )%>%
                           mutate(trait=factor(trait,
                                               levels=rev(names( dummy.col)))),
                         aes(y=Estimate,
                             x=parameters,
                             color=trait)) + geom_blank() + 
     geom_pointrange(aes(xmin=Q2.5,
                         xmax=Q97.5),
                     size=2,alpha=0.8) +
     scale_color_manual(values= dummy.col )  +
     theme_few() +
     theme(legend.position="bottom",
           legend.key.size = unit(1, 'cm'),
           legend.title.position = "top",
           legend.title =element_text(size=20),
           legend.text =element_text(size=16),
           axis.text = element_text(size=18))

  Cool.glm.theory.trait.plotlist[[country]] <- plot_grid( Inter.plot.sum,
                                                          Intra.plot.sum,
                                                          Lambda.plot.sum,
                                                          axis="tblr",
                                                          align="v",
                                                          ggpubr::get_legend(legend.plot),
                                                          nrow=4,rel_heights = c(1,0.6,0.6,0.3))
}
#figures/GLM.traits.AUS.pdf
Cool.glm.theory.trait.plotlist[[paste0("spain")]]#figures/GLM.traits.SPAIN.pdf
Cool.glm.theory.trait.plotlist[[paste0("aus")]]#figures/GLM.traits.AUS.pdf
#---- 1.3. Make graph for main text - INTRA ----
dummy <- data.frame(country = c("spain","aus"), 
                    Estimate = c(-0.025),
                    )

for( country in country.list){
  trait.to.remove <- c("Root diameter",
    #"Flower width","Seed mass",
    "Leaf area index",
    "Canopy shape")
  Cool.glm.theory.trait.plotlist[[paste0(country,"_intra")]]  <- summary.table.for.plot.glm[[country]]$Intra.trait.df.i %>%
  
    dplyr::filter(density.quantile  %in% c("low")) %>%
    mutate(facet.name=case_when(density.quantile == "low" ~ "Intraspecific interaction at\nLow conspecific's density",
                                density.quantile == "high" ~ "High neighbor's density")) %>%
    bind_rows( summary.table.for.plot.glm[[country]]$Lambda.trait.df.i %>%
                 mutate(facet.name="Intrinsic fecundity") ) %>%
    dplyr::filter( outlier== "removed") %>%
    dplyr::filter(!trait %in% trait.to.remove)%>%
    dplyr::filter(!parameters =="intercept") %>%
    mutate(parameters= factor(parameters,
                              levels=c("Focal trait")))%>%
    mutate(facet.name= factor(facet.name,
                              levels=c("Intraspecific interaction at\nLow conspecific's density",
                                       #"High neighbor's density",
                                       "Intrinsic fecundity"))) %>%
    mutate(y_numb= case_when(parameters =="Focal trait" ~ 1)) %>%
    mutate(trait=factor(trait, 
                        names(dummy.col)[!names(dummy.col) %in%trait.to.remove])) %>%
    mutate(y_trait=((as.numeric(trait)-6)*0.01 +y_numb)) %>%
    ggplot(aes(y=y_trait,
               x=Estimate,
               color=as.factor(trait),
               group=as.factor(trait)#,
               #shape = pointshape 
    )) + 
    #shape=density.quantile)) +
    geom_linerange(aes(xmin=Q2.5,
                       xmax=Q97.5),
                   size=1,alpha=0.7) +
    geom_pointrange(aes(xmin=Q2.5,
                        xmax=Q97.5),
                    size=2) +
    geom_text(aes(y=y_trait, x=Estimate,
                  group=as.factor(trait),label=signif),
              size=10,color="black") +
    scale_color_manual(values=dummy.col)+
    scale_y_continuous(labels=rev(c("Focal trait")),
                       #limits=c(0.92,1.07),
                       expand = c(0.01,0.01),
                       breaks=0) +
    facet_wrap(.~facet.name, ncol=3, scale="free_x") +
    theme_bw() +
    geom_vline(xintercept=0) + 
    labs(x="Effect size",y="Focal trait",
         color="Functional trait") +
    #coord_cartesian(expand=F) +
    guides(color = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          strip.text = element_text(size=20),
          strip.background = element_rect(fill="grey99"),
          legend.title =element_text(size=20),
          legend.text =element_text(size=18),
          axis.text = element_text(size=18),
          axis.title =element_text(size=20),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank(),
          plot.margin=unit(c(1,0.5,0.5,1),"cm"),
          panel.spacing = unit(3, "lines"))
  Cool.glm.theory.trait.plotlist[[paste0(country,"_intra")]]

}

legend.plot <- ggplot(data=Cool.theory.trait.df[["aus"]]$Intra.trait.df %>%
                        bind_rows(Cool.theory.trait.df[["spain"]]$Intra.trait.df) %>%
                        dplyr::filter(outlier=="included") %>%
                        #mutate(pointshape = rep(c("Trait effect on interaction",
                        #                         "Intercept of interaction strength"),
                        #                       each=120)) %>%
                        dplyr::filter(!parameters =="intercept") %>%
                        dplyr::filter(!trait %in% c("Root diameter",
                                                    #"Flower width","Seed mass",
                                                    "Leaf area index",
                                                    "Canopy shape")) %>%
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
  scale_color_manual(values= dummy.col )  +
  guides(shape= guide_legend(nrow=2,byrow=TRUE,
                             override.aes = list(alpha = 1,
                                                 color=c("grey80"))),
         color= guide_legend(nrow=4,byrow=TRUE)) +
  labs(x="Effect size",
       y=paste0(""),
       shape="Model estimates",
       color="Functional trait") +
  theme_few() +
  theme(legend.position="bottom",
        legend.key.size = unit(1, 'cm'),
        legend.title.position = "top",
        legend.title =element_text(size=20),
        legend.text =element_text(size=18),
        axis.text = element_text(size=18),
        axis.title =element_text(size=20),
        plot.margin=unit(c(1,0,0,1),"cm"))
legend.plot

GLM.traits.INTRA <- plot_grid( ggarrange(Cool.glm.theory.trait.plotlist[[paste0("aus_intra")]],
                                         Cool.glm.theory.trait.plotlist[[paste0("spain_intra")]],
                                         common.legend = T, legend = "none",
                                         label.x = c(-0.08,-0.06),
                                         label.y = 1.01,align="v",
                                         font.label = list(size = 26, color = "black", 
                                                           face = "bold", family = NULL),
                                         ncol=1,labels=c("a. Australia","b. Spain","")),
                               ggpubr::get_legend(legend.plot),
                               ncol=1,
                               rel_heights =c(1,0.2),
                               nrow=2,legend="none")
GLM.traits.INTRA

#figures/GLM.traits.INTRA.pdf

#---- 1.5. Make graph for main text -INTRA - LOW TO Lambda ----

trait.to.remove <- c("Root diameter",
                     #"Flower width","Seed mass",
                     "Leaf area index",
                     "Canopy shape")
"Intraspecific interaction at\nLow conspecific's density"
view(Cool.glm.theory.trait.plotlist[[paste0(country,"_intra_diag")]])
Cool.glm.theory.trait.plotlist[[paste0(country,"_intra_diag")]]  <- summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i %>%
  dplyr::filter(density.quantile  %in% c("low") &
                  parameters == "Focal trait") %>%
    mutate(parameters="FT") %>%
    bind_rows( summary.table.for.plot.glm[["aus"]]$Lambda.trait.df.i %>%
                 dplyr::filter(density.quantile  %in% c("low") &
                                 parameters == "Focal trait") %>%
               mutate(parameters="IF")) %>%
    mutate(country="Australia") %>%
  bind_rows(summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
            dplyr::filter(density.quantile  %in% c("low")&
                            parameters == "Focal trait") %>%
            mutate(parameters="FT") %>%
              bind_rows( summary.table.for.plot.glm[["spain"]]$Lambda.trait.df.i %>%
                           dplyr::filter(density.quantile  %in% c("low") &
                                           parameters == "Focal trait") %>%
                         mutate(parameters="IF")) %>%
              mutate(country="Spain")) %>%
  dplyr::filter( outlier== "removed") %>%
  dplyr::select(-signif) %>%
  dplyr::filter(!trait %in% trait.to.remove) %>%
  mutate(trait=factor(trait, 
                      levels=rev(names(dummy.col)[!names(dummy.col) %in% trait.to.remove]))) %>%
  pivot_wider(names_from = parameters, 
              values_from = c('Estimate','Q2.5','Q97.5'), names_sep="") %>%
  ggplot(aes(y=EstimateIF,
             x=EstimateFT,
             color=as.factor(trait),
             group=as.factor(trait),
             shape = country )) +
  #geom_rect(aes(xmin=-Inf, xmax=0, ymin=0, ymax=+Inf), color=NA,
  #          fill='grey98', alpha=0.98) +
  #geom_rect(aes(xmin=0, xmax=Inf, ymin=-Inf, ymax=0), color=NA,
 #           fill='grey98', alpha=0.98) +
  #shape=density.quantile)) +
  geom_point( size=8,alpha=1) +
  geom_pointrange(aes(ymin=Q2.5IF,
                      ymax=Q97.5IF),
                  size=2,alpha=0.3) +
  geom_errorbarh(aes(xmin=Q2.5FT,
                     xmax=Q97.5FT),
                 size=1,alpha=0.3,
                 height=0) +
  geom_segment(aes(x=0, xend = +0.15 , 
                   y=0, yend = 0), size=1,
               arrow = arrow(length = unit(0.6,"cm")),
               color="black")+
  geom_segment(aes(x=0, xend =0 , 
                   y=-1, yend = 1), size=1,
               arrow = arrow(length = unit(0.6,"cm")),
               color="black")+
  geom_segment(aes(x=0, xend = -0.15 , 
                   y=0, yend = 0), size=1,
               arrow = arrow(length = unit(0.6,"cm")),
               color="black")+
  geom_segment(aes(x=0, xend =0 , 
                   y=1, yend = -1), size=1,
               arrow = arrow(length = unit(0.6,"cm")),
               color="black")+
  annotate(geom = "text",label="Increases intrinsic fecundity",
           x=-0.008,y=0.55,size=6,angle=90)+
  annotate(geom = "text",label="Decreases intrinsic fecundity",
           x=0.008,y=-0.55,size=6,angle=-90)+
  annotate(geom = "text",label="Increases self-competition",
           x=-0.065,y=-0.05,size=6,angle=0)+
  annotate(geom = "text",label="Increases self-facilitation",
           x=0.065,y=0.05,size=6,angle=0)+
  scale_shape_manual(values=c(16:17)) +
  scale_color_manual(values=dummy.col)+
  scale_x_continuous(breaks=c(-0.03,0,0.03),
                     labels=c(-0.03,0,0.03)) +
  labs(x="Effect size on interactions at Low neighbor's density",
       y="Effect size on intrinsic fecundity",
       shape="Community",
       color="Functional trait") +
  coord_cartesian( xlim = c(-.15,0.15), ylim = c(-1,1),
                 expand = F, default = FALSE, clip = "on") +
  theme_bw() +
  guides(color = guide_legend(title.position = "top",
                              ncol=2),
         shape = guide_legend(title.position = "top",
                              nrow=2)) +
  theme(legend.position="bottom",
        strip.text = element_text(size=18),
        legend.title =element_text(size=16),
        legend.text =element_text(size=14),
        axis.title=element_text(size=20),
        axis.text = element_text(size=18),
        panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_line( color = "grey80",linetype="dashed"),
        panel.grid.minor = element_blank())
Cool.glm.theory.trait.plotlist[[paste0(country,"_intra_diag")]]
#figures/main/Oblique.INTRA.pdf
Lambda.trait <- NULL
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  if(country=="spain"){
    competition_df <- get(paste0("clean.data.",country))[[paste0("competition_",country)]]  %>%
      gather(any_of(Code.focal.list),key="species",value="individuals") %>%
      mutate(individuals = case_when(focal==species ~ (individuals + 1),
                                     T~individuals)) %>%
      mutate(com_id = paste0(plot,"_", subplot)) %>%
      group_by(species) %>%
      summarise(low.abundance = as.vector(quantile(individuals[individuals>0],0.1,na.rm=T)),
                median.abundance = as.vector(quantile(individuals[individuals>0],0.5,na.rm=T)),
                high.abundance = as.vector(quantile(individuals[individuals>0],0.9,na.rm=T)))
  }else{
    competition_df <- get(paste0("clean.data.",country))[[paste0("competition_",country)]]  %>%
      gather(any_of(Code.focal.list),key="species",value="individuals") %>%
      mutate(individuals = case_when(focal==species ~ (individuals + 1),
                                     T~individuals)) %>%
      mutate(individuals = (individuals/scale)*15) %>%
      dplyr::rename("com_id"="plot") %>%
      group_by(species) %>%
      summarise(low.abundance = as.vector(quantile(individuals[individuals>0],0.1,na.rm=T)),
                median.abundance = as.vector(quantile(individuals[individuals>0],0.5,na.rm=T)),
                high.abundance = as.vector(quantile(individuals[individuals>0],0.9,na.rm=T)))
  }
  
alpha_intra.n <- Cool.theory.trait.df[[country]]$trait.value.intra.df %>%
    mutate(country=country) %>%
    dplyr::filter(density.quantile=="low") %>%
    dplyr::select(neigh,country,theoretical.effect) %>%
    rename("intra.effect"="theoretical.effect") %>%
    unique()

alpha_inter.n <- Cool.theory.trait.df[[country]]$trait.dist.df %>%
  mutate(country=country) %>%
  dplyr::filter(density.quantile=="low") %>%
  group_by(focal,country) %>%
  summarise(inter.effect.mean = mean(theoretical.effect,na.rm=T),
            inter.effect.sd = sd(theoretical.effect,na.rm=T))%>%
  ungroup()
Lambda.trait.n <- Cool.theory.trait.df[[country]]$trait.value.lambda.df %>%
  mutate(country=country) %>%
  dplyr::filter(density.quantile=="low") %>%
  dplyr::select(-c('density.quantile', 'receiver.trait.scaled')) %>%
  spread(trait,receiver.trait) %>%
  left_join(competition_df  %>%
              dplyr::rename("focal"="species")) %>%
  left_join(alpha_intra.n%>%
  dplyr::rename("focal"="neigh")) %>%
  left_join(alpha_inter.n) %>%
  mutate(scaled.lambda= as.vector(scale(lambda)))


Lambda.trait <- bind_rows(Lambda.trait,Lambda.trait.n)

}
view(Lambda.trait)
str(Lambda.trait)
Lambda.trait %>%
  ggplot(aes(x=theoretical.effect, y=scaled.lambda,size=median.abundance)) +
  geom_point() +
  theme_bw()+
  scale_size_continuous(breaks=c(1,2,5,10,14))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  theme(legend.position="bottom",
        strip.text = element_text(size=18),
        legend.title =element_text(size=16),
        legend.text =element_text(size=14),
        axis.title=element_text(size=20),
        axis.text = element_text(size=18),
        panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_line( color = "grey80",linetype="dashed"),
        panel.grid.minor = element_blank())

write.csv(Lambda.trait,
          file="results/Lambda.trait.abundance.csv")
View(Lambda.trait)
 

#---- 1.4. Make graph for main text -INTER ----

for( country in country.list){
  trait.to.remove <- c("Root diameter",
                       #"Flower width","Seed mass",
                       "Leaf area index",
                       "Canopy shape")
  
  intercept.df <- bind_rows(summary.table.for.plot.glm[[country]]$df.i) %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    #dplyr::filter(!parameters =="intercept") %>%
    dplyr::mutate(trait = case_when(parameters =="intercept" ~ "intercept",
                                    T~trait)) %>%
    dplyr::mutate(parameters = case_when(parameters =="intercept" ~ NA,
                                         T~parameters)) %>%
    fill(parameters,.direction = c("up")) %>%
    filter(trait =="intercept") %>%
    group_by(parameters,trait,density.quantile) %>%
    summarise(Q2.5=min(Q2.5),
              Q97.5=max(Q97.5),
              Estimate =median(Estimate))%>%
    ungroup() %>%
    mutate(pointshape = "Intercept of interaction strength")
  
  
  Cool.glm.theory.trait.plotlist[[paste0(country,"_inter")]]  <- bind_rows(summary.table.for.plot.glm[[country]]$df.i) %>%
    dplyr::filter( outlier== "removed") %>%
    dplyr::filter(density.quantile  %in% c("low","high")) %>%
    dplyr::filter(!trait %in% trait.to.remove )%>%
    mutate(pointshape = "Trait effect on interaction") %>%
    #bind_rows(intercept.df) %>%
    mutate(parameters= factor(parameters,
                              levels=c("Focal trait","Emitter trait",
                                       "Focal trait -\nEmitter trait")))%>%
    mutate(y_numb= case_when(parameters =="Focal trait" ~ 3,
                             parameters =="Emitter trait"~2,
                             parameters =="Focal trait -\nEmitter trait" ~1)) %>%
    mutate(trait=factor(trait, 
                        names(dummy.col)[!names(dummy.col) %in%trait.to.remove])) %>%
    mutate(y_trait=((as.numeric(trait)-6)*0.06 +y_numb)) %>%
    mutate(density.quantile.text = case_when(density.quantile =="low" ~"Low Neighbor's Density",
                                             density.quantile =="high" ~"High Neighbor's Density")) %>%
    mutate(density.quantile.text=factor(density.quantile.text,
                                        levels=c("Low Neighbor's Density","High Neighbor's Density"))) %>%
    ggplot(aes(y=y_trait,
               x=Estimate,
               color=as.factor(trait),
               group=as.factor(trait)#,
               #shape = pointshape 
               )) + 
    #shape=density.quantile)) +
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
    scale_shape_manual(values=c(16:17))+
    scale_y_continuous(breaks=c(1:3),
                       limits=c(0.5,3.5),
                       labels=rev(c("Focal trait",
                                    "Neighbor trait",
                                    "Focal trait -\nNeighbor trait"))) +
    theme_bw() +
    coord_cartesian(xlim=c(-0.04,0.04),clip = "on",expand=F)+
    facet_wrap(.~density.quantile.text,ncol=2) +
    geom_vline(xintercept=0) + 
    labs(x="Effect size",
         y=paste0(""),
         shape="Model estimates",
         color="Functional trait") +
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           shape = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          strip.text = element_text(size=20),
          strip.background = element_rect(fill="grey98"),
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
  Cool.glm.theory.trait.plotlist[[paste0(country,"_inter")]]
  

  
}

legend.plot <- ggplot(data=Cool.theory.trait.df[["aus"]]$Intra.trait.df %>%
                        bind_rows(Cool.theory.trait.df[["spain"]]$Intra.trait.df) %>%
                        dplyr::filter(outlier=="included") %>%
                        #mutate(pointshape = rep(c("Trait effect on interaction",
                         #                         "Intercept of interaction strength"),
                         #                       each=120)) %>%
                        dplyr::filter(!parameters =="intercept") %>%
                        dplyr::filter(!trait %in% c("Root diameter",
                                                    #"Flower width","Seed mass",
                                                    "Leaf area index",
                                                    "Canopy shape")) %>%
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
  scale_color_manual(values= dummy.col )  +
  guides(shape= guide_legend(nrow=2,byrow=TRUE,
                             override.aes = list(alpha = 1,
                                                 color=c("grey80"))),
         color= guide_legend(nrow=4,byrow=TRUE)) +
  labs(x="estimate",
       y=paste0(""),
       shape="Model estimates",
       color="Functional trait") +
  theme_few() +
  theme(legend.position="bottom",
        legend.key.size = unit(1, 'cm'),
        legend.title.position = "top",
        legend.title =element_text(size=20),
        legend.text =element_text(size=18),
        axis.text = element_text(size=18))
legend.plot

GLM.traits.INTER <- plot_grid( ggarrange(Cool.glm.theory.trait.plotlist[[paste0("aus_inter")]],
           Cool.glm.theory.trait.plotlist[[paste0("spain_inter")]],
           common.legend = T, legend = "none",
           label.x = c(-0.05,-.03),
           label.y = 1.01,align="v",
           font.label = list(size = 26, color = "black", 
                             face = "bold", family = NULL),
           ncol=1,labels=c("a. Australia","b. Spain","")),
           ggpubr::get_legend(legend.plot),
           ncol=1,
          rel_heights =c(1,0.2),
          nrow=2,legend="none")
GLM.traits.INTER
#figures/GLM.traits.INTER.pdf
ggsave(last_plot(),
       width=15.41,
       height=10.30,
       unit="in",
       file="figures/GLM.traits.INTER.pdf")

library(grid)
library(pBrackets) 
#figures/GLM.traits.INTER.pdf
#---- 1.5. Make graph for main text -INTER - LOW TO HIGH ----

trait.to.remove <- c("Root diameter",
                       #"Flower width","Seed mass",
                       "Leaf area index",
                       "Canopy shape")

Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag_focal")]]  <- summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(country="Spain")) %>%
    dplyr::filter( outlier== "removed") %>%
    dplyr::select(-signif) %>%
    dplyr::filter(density.quantile  %in% c("low","high")) %>%
    dplyr::filter(!trait %in% trait.to.remove )%>%
    dplyr::filter(parameters==c("Focal trait")) %>%
    mutate(trait=factor(trait, 
                      levels=rev(names(dummy.col)[!names(dummy.col) %in% trait.to.remove]))) %>%
    mutate(Q2.5 = case_when(Estimate<0 ~ -Q2.5,
                            T ~ Q2.5),
           Q97.5 = case_when(Estimate<0 ~ -Q97.5,
                             T ~ Q97.5)) %>%
    pivot_wider(names_from = density.quantile, 
                values_from = c('Estimate','Q2.5','Q97.5'), names_sep="") %>%
    ggplot(aes(y=abs(Estimatehigh),
               x=abs(Estimatelow),
               color=as.factor(trait),
               group=as.factor(trait),
               shape = country )) + 
    #shape=density.quantile)) +
    geom_point( size=8,alpha=1) +
    geom_pointrange(aes(ymin=Q2.5high,
                        ymax=Q97.5high),
                    size=2,alpha=0.3) +
    geom_errorbarh(aes(xmin=Q2.5low,
                        xmax=Q97.5low),
                    size=1,alpha=0.3,
                   height=0) +
    geom_abline(intercept=0,slope=1) +
    scale_shape_manual(values=c(16:17)) +
    scale_color_manual(values=dummy.col)+
    scale_y_continuous(breaks=c(0,0.01,0.03),labels=c(0,0.01,0.03)) +
    scale_x_continuous(breaks=c(0,0.01,0.03),labels=c(0,0.01,0.03)) +
    labs(x="Effect size at Low neighbor's density",
         y="Effect size at High neighbor's density",
         shape="Community",
         color="Functional trait") +
    annotate(geom = "polygon", x = c(Inf, -Inf, -Inf), 
           y = c(Inf, -Inf, Inf), fill = "grey70", alpha = 0.1 )+
    annotate(geom = "text",label="Effect size at Low > High",
           x=0.033,y=0.032,size=6,angle=45)+
  coord_cartesian( xlim = c(0,0.04), ylim=c(0,0.04),
                   expand = F, default = FALSE, clip = "on") +
  theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=3),
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
          panel.grid.major = element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())

  
  Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag_neigh")]]  <- summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(country="Australia") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(country="Spain")) %>%
    dplyr::filter( outlier== "removed") %>%
    dplyr::select(-signif) %>%
    dplyr::filter(density.quantile  %in% c("low","high")) %>%
    dplyr::filter(!trait %in% trait.to.remove )%>%
    dplyr::filter(parameters==c("Emitter trait")) %>%
    mutate(trait=factor(trait, 
                        levels=rev(names(dummy.col)[!names(dummy.col) %in% trait.to.remove]))) %>%
    mutate(Q2.5 = case_when(Estimate<0 ~ -Q2.5,
                            T ~ Q2.5),
           Q97.5 = case_when(Estimate<0 ~ -Q97.5,
                             T ~ Q97.5)) %>%
    pivot_wider(names_from = density.quantile, 
                values_from = c('Estimate','Q2.5','Q97.5'), names_sep="") %>%
    ggplot(aes(y=abs(Estimatehigh),
               x=abs(Estimatelow),
               color=as.factor(trait),
               group=as.factor(trait),
               shape = country )) + 
    #shape=density.quantile)) +
    geom_point( size=8,alpha=1) +
    geom_pointrange(aes(ymin=Q2.5high,
                        ymax=Q97.5high),
                    size=2,alpha=0.3) +
    geom_errorbarh(aes(xmin=Q2.5low,
                       xmax=Q97.5low),
                   size=1,alpha=0.3,
                   height=0) +
    geom_abline(intercept=0,slope=1) +
    scale_shape_manual(values=c(16:17)) +
    scale_color_manual(values=dummy.col)+
   scale_y_continuous(breaks=c(0,0.01,0.03),labels=c(0,0.01,0.03)) +
    scale_x_continuous(breaks=c(0,0.01,0.03),labels=c(0,0.01,0.03)) +
    labs(x="Effect size at Low neighbor's density",
         y="Effect size at High neighbor's density",
         shape="Community",
         color="Functional trait") +
    annotate(geom = "polygon", x = c(Inf, -Inf, -Inf), 
             y = c(Inf, -Inf, Inf), fill = "grey70", alpha = 0.1 )+
    annotate(geom = "text",label="Effect size at Low > High",
             x=0.033,y=0.032,size=6,angle=45)+
    coord_cartesian( xlim = c(0,0.04), ylim=c(0,0.04),
                     expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=3),
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
          panel.grid.major = element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())
  Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag_neigh")]]

  ggarrange(Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag_focal")]],
             Cool.glm.theory.trait.plotlist[[paste0(country,"_inter_diag_neigh")]],
            common.legend = T, legend = "bottom",
            label.x = c(-0.08,-0.12),
            label.y = 1.01,align="h",
            font.label = list(size = 26, color = "black", 
                              face = "bold", family = NULL),
             labels=c("a. Focal's trait", "b. Neighbor's trait"),
             ncol=2)
  #figures/main/Oblique.INTER.pdf 
#---- 1.6. Make graph for main text -INTER - INTRA ----
  
  trait.to.remove <- c("Root diameter",
                       #"Flower width","Seed mass",
                       "Leaf area index",
                       "Canopy shape")
 
  plot_inter_intra  <- summary.table.for.plot.glm[["aus"]]$df.i %>%
    mutate(model="inter") %>%
    bind_rows( summary.table.for.plot.glm[["aus"]]$Intra.trait.df.i %>%
                 mutate(model="intra")) %>%
    mutate(country="aus") %>%
    bind_rows(summary.table.for.plot.glm[["spain"]]$df.i %>%
                mutate(model="inter") %>%
                bind_rows( summary.table.for.plot.glm[["spain"]]$Intra.trait.df.i %>%
                             mutate(model="intra")) %>%
                mutate(country="spain")) %>%
    dplyr::filter( outlier== "removed") %>%
    dplyr::select(-signif) %>%
    dplyr::filter(density.quantile  %in% c("low")) %>%
    dplyr::filter(!trait %in% trait.to.remove )%>%
    dplyr::filter(parameters==c("Focal trait")) %>%
    mutate(trait=factor(trait, 
                        levels=rev(names(dummy.col)[!names(dummy.col) %in% trait.to.remove]))) %>%
    pivot_wider(names_from = model, 
                values_from = c('Estimate','Q2.5','Q97.5'), names_sep="") %>%
    ggplot() + 
    #geom_pointrange(data= Lambda.trait,
     #               aes(y=inter.effect.mean,
       #                 x=intra.effect,
        #                ymin=inter.effect.mean-inter.effect.sd,
          #              ymax=inter.effect.mean+inter.effect.sd,
           #             shape = country,
            #            size=median.abundance),
             ##       color="black",
              #      alpha=0.7) +
    #shape=density.quantile)) +
    geom_point(aes(y=Estimateinter,
                   x=Estimateintra,
                   color=as.factor(trait),
                   group=as.factor(trait),
                   shape = country ),
               size=8,alpha=1) +
    geom_pointrange(aes(y=Estimateinter,
                        x=Estimateintra,
                        color=as.factor(trait),
                        group=as.factor(trait),
                        shape = country,
                        ymin=Q2.5inter,
                        ymax=Q97.5inter),
                    size=2,alpha=0.3) +
    geom_errorbarh(aes(y=Estimateinter,
                       xmin=Q2.5intra,
                       xmax=Q97.5intra,
                       color=as.factor(trait)),
                   size=1,alpha=0.3,
                   height=0) +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    #geom_abline(intercept=0,slope=1) +
    scale_shape_manual(values=c(16,17)) +
    scale_color_manual(values=dummy.col)+
    scale_size(range=c(1,2))+
    scale_y_continuous(breaks=seq(-0.1,0.1,0.025)) +
    scale_x_continuous(breaks=seq(-0.1,0.1,0.025)) +
    labs(x="Effect size on INTRA",
         y="Effect size on INTER",
         shape="Community",
         color="Functional trait") +
    geom_segment(data=NULL,x=0, xend = +0.12 , 
                     y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(data=NULL,x=0, xend =0 , 
                     y=0, yend = 0.05, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(data=NULL,x=0, xend = -0.12 , 
                     y=0, yend = 0, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(data=NULL,x=0, xend =0 , 
                     y=0, yend = -0.05, size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom = "text",label="Increases competition \nfrom others",
             x=0.002,y=-0.045,size=6,angle=0)+
    annotate(geom = "text",label="Increases facilitation \nfrom others",
             x=0.002,y=0.045,size=6,angle=0)+
    annotate(geom = "text",label="Increases \nself-competition",
             x=-0.095,y=0.00,size=6,angle=0)+
    annotate(geom = "text",label="Increases \nself-facilitation",
             x=0.095,y=0.00,size=6,angle=0)+
    coord_cartesian( ylim = c(-0.05,0.05), xlim=c(-0.12,0.12),
                     expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=3),
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
          panel.grid.major = element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())
  
  
  plot_inter_intra
  #figures/main/Oblique.INTER.INTRA.pdf 
  #---- 1.7. Make graph for main text -INTER - lambda ----
  
  trait.to.remove <- c("Root diameter",
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
    dplyr::filter( outlier== "removed") %>%
    dplyr::select(-signif) %>%
    dplyr::filter(density.quantile  %in% c("low")) %>%
    dplyr::filter(!trait %in% trait.to.remove )%>%
    dplyr::filter(parameters==c("Focal trait")) %>%
    mutate(trait=factor(trait, 
                        levels=rev(names(dummy.col)[!names(dummy.col) %in% trait.to.remove]))) %>%
    pivot_wider(names_from = model, 
                values_from = c('Estimate','Q2.5','Q97.5'), names_sep="") %>%
    ggplot(aes(y=Estimateinter,
               x=Estimatelambda,
               color=as.factor(trait),
               group=as.factor(trait),
               shape = country )) + 
    #shape=density.quantile)) +
    geom_point( size=8,alpha=1) +
    geom_pointrange(aes(ymin=Q2.5inter,
                        ymax=Q97.5inter),
                    size=2,alpha=0.3) +
    geom_errorbarh(aes(xmin=Q2.5lambda,
                       xmax=Q97.5lambda),
                   size=1,alpha=0.3,
                   height=0) +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    #geom_abline(intercept=0,slope=1) +
    scale_shape_manual(values=c(16:17)) +
    scale_color_manual(values=dummy.col)+
    #scale_y_continuous(breaks=seq(-0.1,0.1,0.025)) +
    #scale_x_continuous(breaks=seq(-0.1,0.1,0.05)) +
    labs(x="Effect size on Lambda",
         y="Effect size on INTER",
         shape="Community",
         color="Functional trait") +
    geom_segment(aes(x=0, xend = +1.2 , 
                     y=0, yend = 0), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(aes(x=0, xend =0 , 
                     y=0, yend = 0.05), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(aes(x=0, xend = -1.2 , 
                     y=0, yend = 0), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    geom_segment(aes(x=0, xend =0 , 
                     y=0, yend = -0.05), size=1,
                 arrow = arrow(length = unit(0.6,"cm")),
                 color="black")+
    annotate(geom = "text",label="Increases competition \nfrom others",
             x=0.002,y=-0.045,size=6,angle=0)+
    annotate(geom = "text",label="Increases facilitation \nfrom others",
             x=0.002,y=0.045,size=6,angle=0)+
    annotate(geom = "text",label="Increases \nintrinsic fecundity",
             x=0.9,y=0,size=6,angle=0)+
    annotate(geom = "text",label="Decreases \nintrinsic fecundity",
             x=-0.9,y=0,size=6,angle=0)+
    coord_cartesian( ylim = c(-0.05,0.05), xlim=c(-1.2,1.2),
                     expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(color = guide_legend(title.position = "top",
                                nrow=3),
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
          panel.grid.major = element_line( color = "grey70",linetype="dashed"),
          panel.grid.minor = element_blank())
  
  
  plot_inter_lambda
  #figures/main/Oblique.INTER.INTRA.LAMBDA.pdf 
  ggarrange(plot_inter_intra,
            plot_inter_lambda,
            ncol=2,
            common.legend = T, legend="bottom",
            label.x= c(-0.28,-0.35),
            align="h",
            font.label = list(size = 20, color = "black", 
                              face = "bold", family = NULL),
            labels=c("a. Mitigation of competition by self-facilitation",
                            "b. Mitigation of competitione by high intrinsic fecunidty"))
#---- 1.6.Summary across density----
for( country in country.list){
  
  Cool.glm.theory.trait.plotlist[[paste0(country,"_summary")]]  <- bind_rows(summary.table.for.plot.glm[[country]]$df.i) %>%
    dplyr::filter(!parameters =="intercept")%>%
    mutate(signif = case_when( (round(Q2.5,digits =3) <=0 &
                                  round(Q97.5,digits =3)-0.002 <=0) ~ "*",
                               (round(Q2.5,digits =3)+0.002 >=0 &
                                  round(Q97.5,digits =3) >=0) ~ "*",
                               T ~ "")) %>%
    dplyr::filter(!signif =="")%>%
    dplyr::filter( outlier== "removed") %>% #included
    dplyr::filter(!trait %in% trait.to.remove)%>%
    mutate(parameters= factor(parameters,
                              levels=c("Focal trait","Emitter trait",
                                       "Focal trait -\nEmitter trait")))%>%
    mutate(y_numb= case_when(parameters =="Focal trait" ~ 3,
                             parameters =="Emitter trait"~2,
                             parameters =="Focal trait -\nEmitter trait" ~1)) %>%
    mutate(trait=factor(trait, 
                        names(dummy.col)[!names(dummy.col) %in% trait.to.remove])) %>%
    mutate(y_trait=((as.numeric(trait)-6)*0.05 +y_numb)) %>%
    ggplot(aes(y=y_trait,
               x=Estimate,
               color=as.factor(trait),
               group=as.factor(trait) )) + 
    #shape=density.quantile)) +
    geom_linerange(aes(xmin=Q2.5,
                       xmax=Q97.5),
                   size=1,alpha=0.7,
                   position=position_dodge2(width=0.5)) +
    geom_pointrange(aes(xmin=Q2.5,
                        xmax=Q97.5),
                    size=2,alpha=0.7,
                    position=position_dodge2(width=0.5)) +
    scale_color_manual(values=dummy.col)+
    scale_shape_manual(values=c(16:19))+
    scale_y_continuous(breaks=c(1:4),
                       labels=rev(c("Focal trait",
                                    "Neighbor trait",
                                    "Focal trait -\nNeighbor trait",
                                    "Focal trait"))) +
    theme_bw() +
    facet_wrap(.~density.quantile, ncol=4) +
    geom_vline(xintercept=0) + 
    labs(x="estimate",
         y=paste0(""),
         color="Density of emitter") +
    guides(color = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          plot.margin = unit(c(1,ifelse(country=="aus",
                                        0,3),0,
                               ifelse(country=="aus",
                                      3,0)),"cm"),
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=18),
          axis.title=element_text(size=18),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
  Cool.glm.theory.trait.plotlist[[paste0(country,"_summary")]]
  
}
legend.plot <- ggplot(data=Cool.theory.trait.df[["aus"]]$Intra.trait.df %>%
                        bind_rows(Cool.theory.trait.df[["spain"]]$Intra.trait.df) %>%
                        dplyr::filter(!parameters =="intercept") %>%
                        dplyr::filter(!signif =="") %>%
                        dplyr::filter(!trait %in% c("Root diameter",
                                                    "Flower width","Seed mass",
                                                    "Leaf area index",
                                                    "Canopy shape","Stem height")) %>%
                        #bind_rows(intercept.df) %>%
                        mutate(trait=factor(trait,
                                            levels=rev(names( dummy.col)))),
                      aes(y=Estimate,
                          x=parameters,
                          color=trait))+
  #shape=density.quantile )) + 
  geom_blank() + 
  geom_pointrange(aes(xmin=Q2.5,
                      xmax=Q97.5),
                  size=2,alpha=0.8) +
  scale_color_manual("",values= dummy.col ) +
  theme_few() +
  theme(legend.position="bottom",
        legend.key.size = unit(1, 'cm'),
        legend.title.position = "top",
        legend.title =element_text(size=20),
        legend.text =element_text(size=16),
        axis.text = element_text(size=18))
legend.plot

GLM.traits.summary <- plot_grid( ggarrange(Cool.glm.theory.trait.plotlist[[paste0("aus_summary")]],
                                           Cool.glm.theory.trait.plotlist[[paste0("spain_summary")]],
                                           common.legend = T, legend = "none",
                                           label.x = c(0.18,0.2),
                                           label.y = 1.01,align="v",
                                           font.label = list(size = 26, color = "black", 
                                                             face = "bold", family = NULL),
                                           ncol=1,labels=c("a. Australia","b. Spain","")),
                                 ggpubr::get_legend(legend.plot),
                                 ncol=1,
                                 rel_heights =c(1,0.2),
                                 nrow=2,legend="none")
GLM.traits.summary

#figures/glm/GLM.traits.INTER.Intercept.to.High.pdf
#---- 1.7. Outliers ----
for( country in country.list){

  Cool.glm.theory.trait.plotlist[[paste0(country,"_removed.outlier")]]  <- bind_rows(summary.table.for.plot.glm[[country]]$df.i) %>%
    dplyr::filter(density.quantile  %in% c("low","high")) %>%
    dplyr::filter(!parameters =="intercept") %>%
    mutate(signif = case_when( (round(Q2.5,digits =3) <=0 &
                                  round(Q97.5,digits =3)-0.002 <=0) ~ "*",
                               (round(Q2.5,digits =3)+0.002 >=0 &
                                  round(Q97.5,digits =3) >=0) ~ "*",
                               T ~ "")) %>%
    #dplyr::filter(!signif =="")%>%
    dplyr::filter(!trait %in% c("Root diameter",
                                "Flower width","Seed mass",
                                "Leaf area index",
                                "Canopy shape")) %>%
    mutate(parameters= factor(parameters,
                              levels=c("Focal trait","Emitter trait",
                                       "Focal trait -\nEmitter trait")))%>%
    mutate(y_numb= case_when(parameters =="Focal trait" ~ 3,
                             parameters =="Emitter trait"~2,
                             parameters =="Focal trait -\nEmitter trait" ~1)) %>%
    mutate(trait=factor(trait, 
                        names(dummy.col))) %>%
    mutate(y_trait=((as.numeric(trait)-6)*0.035 +y_numb)) %>%
    ggplot(aes(y=y_trait,
               x=Estimate,
               color=as.factor(trait),
               shape=outlier))+
    #shape=density.quantile)) +
    geom_linerange(aes(xmin=Q2.5,
                       xmax=Q97.5),
                   size=1,alpha=0.7,
                   position=position_dodge2(width=0.5)) +
    geom_pointrange(aes(xmin=Q2.5,
                        xmax=Q97.5),
                    size=2,alpha=0.7,
                    position=position_dodge2(width=0.5)) +
    scale_color_manual(values=dummy.col)+
    scale_y_continuous(breaks=c(1:3),
                       labels=rev(c("Focal trait",
                                    "Neighbor trait",
                                    "Focal trait -\nNeighbor trait"))) +
    theme_bw() +
    facet_wrap(.~density.quantile, ncol=2) +
    geom_vline(xintercept=0) + 
    labs(x="estimate",
         y=paste0(""),
         color="Density of emitter") +
    guides(color = guide_legend(title.position = "top",
                                nrow=2)) +
    theme(legend.position="bottom",
          strip.text = element_text(size=18),
          legend.title =element_text(size=18),
          legend.text =element_text(size=18),
          axis.title=element_text(size=18),
          axis.text = element_text(size=18),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_rect( color = "grey60"),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
  Cool.glm.theory.trait.plotlist[[paste0(country,"_removed.outlier")]]
  
}

Cool.glm.theory.trait.plotlist[[paste0("aus_removed.outlier")]]#figures/glm/GLM.outliers.traits.AUS.pdf
Cool.glm.theory.trait.plotlist[[paste0("spain_removed.outlier")]]#figures/glm/GLM.outliers.traits.SPAIN.pdf



#---- 1.8. Supp figures----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 2.1. Observed Trait Density figures----

    trait.density.plot <- Cool.theory.trait.df[["spain"]]$trait.dist.df %>%
      bind_rows( Cool.theory.trait.df[["aus"]]$trait.dist.df) %>%
  ggplot(aes(group = country, color=country))

dummy.names <- list("SRL"="SRL (cm/g)","SRA"="SRA (cm^2/g)" ,"Root length"="Root length \n(cm)",
                 "Root tips"="Root tips" ,
                 "Root diameter"="Root diameter \n(mm)" ,
                 "Root mass density"="Root mass \ndensity \n(mg/mm3)",
                 "Flower width"= "Flower width \n(cm)" ,"Seed mass"="Seed mass \n(mg)",
                 "C13 water use efficiency"="C13 water \nuse efficiency \n(per ml)",
                 "Leaf C to N ratio"= "Leaf C to N ratio" ,
                 "Leaf area index"="Leaf area index" ,
                 "Canopy shape"="Canopy shape",
                 "SLA"="SLA (cm^2/g)","Stem height"="Stem height \n(cm)")

trait_labeller <- function(variable,value){
  return(dummy.names[value])
}

dummy.col <- c("SRL"="#4E79A7FF","SRA"="#76B7B2FF" ,"Root length"="#A4BED5FF","Root tips"="#512DA8FF" ,
               "Root diameter"="#F28E2BFF" , "Root mass density"="#ED645AFF",
               "Flower width"= "#FF9DA7FF" ,"Seed mass"="#B276B2FF",
               "C13 water use efficiency"="#9C755FFF",
               "Leaf C to N ratio"= "#BCBD22FF" ,
               "Leaf area index"="#D4E157FF" ,"Canopy shape"="#72874EFF",
               "SLA"="#59A14FFF","Stem height"="#FED789FF",
               "intercept"="black")
plant_code_spain <- read.csv( "data/spain_rawdata/plant_code_spain.csv",
                              header = T, stringsAsFactors = F, sep=",",
                              na.strings = c("","NA"))
plant_code_aus <- read.csv("data/aus_rawdata/plant_code_aus.csv",
                           header = T, stringsAsFactors = F, sep=",",
                           na.strings = c("","NA"))

library(ggridges)
get(paste0("clean.data.spain"))[["plant_traits"]] %>%
  dplyr::select(-c("Mean fecundity")) %>%
  rownames_to_column("species") %>%
  gather(any_of(names(dummy.col)),
         key="trait",value="trait.value") %>%
  mutate(country="spain") %>%
  left_join(plant_code_spain %>% dplyr::select(family,code.analysis) %>%
              unique() %>%
              rename("species"="code.analysis")) %>%
  bind_rows(get(paste0("clean.data.aus"))[["plant_traits"]] %>%
              dplyr::select(-c("Mean fecundity")) %>%
              rownames_to_column("species") %>%
              gather(any_of(names(dummy.col)),
                     key="trait",value="trait.value") %>%
              mutate(country="aus") %>%
              left_join(plant_code_aus %>% dplyr::select(family,final.code)%>%
                          unique() %>%
                          rename("species"="final.code"))) %>%
 
  mutate(ynum= case_when(trait %in% c("C13 water use efficiency",
                                      "Canopy shape","Root mass density",
                                      "SLA","SRL","Stem height") ~ as.numeric(as.factor(country)),
                                                     T~ 1)) %>%
  ggplot() +
  ggridges::geom_density_ridges2(aes(y=country,
                                     x=trait.value),
                                 fill="grey80",alpha=0.8) +
  geom_point(aes(y=ynum-0.01,
                 x=trait.value,
                 color=family),alpha=0.9,
             position=position_dodge2(  width =0.16,
                                        preserve = "total",
                                        padding = 0.1,
                                        reverse = FALSE),
             shape = 15, size = 4) +
  facet_wrap(trait~.,scale="free",nrow=2,
             labeller=trait_labeller) +
  scale_color_manual(values =c("#FBE183FF",   "#FE9B00FF","#D8443CFF",
                              "#E6A2A6FF","#9F5691FF","#633372FF",
                              "#1F6E9CFF", "#2B9B81FF","#92C051FF")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape=15,
                                                   size=10,
                                                   alpha = 1))) +
  theme(legend.position="bottom",
        strip.text =element_text(size=16),
        legend.text =element_text(size=16),
        axis.text=element_text(size=16),
        plot.title = element_text(size=20),
        axis.title=element_text(size=16))
library(paletteer)
#figures/Traits.density.pdf

#---- 2.2. Scaled Trait Density figures----

trait.density.df <- Cool.theory.trait.df[["spain"]]$trait.dist.df %>%
  bind_rows( Cool.theory.trait.df[["aus"]]$trait.dist.df) %>%
  ggplot(aes(group = country, color=country))
trait.density.df %>%
ggplot() +
  geom_density_ridges(aes(x=trait.emmiter))
+facet_grid(country~trait,nrow=2)
