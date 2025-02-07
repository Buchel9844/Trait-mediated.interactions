
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
      mutate(trait.dist=trait.dist,
             scaled.trait.dist=as.vector(scale(trait.dist)),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list))) %>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>% 
                  rownames_to_column("focal") %>% rename("receiver.trait"=i) %>%
                  mutate(receiver.trait = as.vector(scale(receiver.trait))))%>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>%
                  rownames_to_column("neigh") %>% rename("emitter.trait"=i)%>%
                  mutate(emitter.trait = as.vector(scale(emitter.trait))))
    
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
                dplyr::select(focal,receiver.trait,trait) %>% unique(),
              relationship ="many-to-many")
  Lambda.trait.df <- NULL
  glm.lambda.trait.summary <- NULL
  #trait.i = "SLA"
  for( trait.i in names(trait.df)){
    for(n in density.quantile.name){
      trait.lambda.df.i <-   trait.value.lambda.df  %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(density.quantile==n) %>%
        mutate(focal=as.factor(focal))
      
      glm.lambda.trait.i <- glmmTMB(lambda ~ 1 + receiver.trait + (1|focal) ,
                                trait.lambda.df.i,
                                family="gaussian")
      
      
      Lambda.trait.df.i <- as.data.frame(confint(glm.lambda.trait.i)) %>%
        mutate(model = "scaled.trait.dist") %>%
        rownames_to_column("parameters") %>%
        rename("Q2.5"="2.5 %",
               "Q97.5"="97.5 %") %>%
        mutate(trait=trait.i,
               density.quantile=n,
               signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                                (Q2.5>0 & Q97.5>0 )~"*",
                                T~""))
      
      Lambda.trait.df <- bind_rows( Lambda.trait.df, Lambda.trait.df.i)
      
    }
  }
  # Make data frae with trait and INTRA specific interactions
  trait.value.intra.df <- Theoretical.Int.list[[country]] %>%
    left_join(specific.trait.dist,
              relationship ="many-to-many")
  Intra.trait.df <- NULL
  glm.intra.trait.summary <- NULL
  #trait.i = "SLA"
  for( trait.i in names(trait.df)){
    for(n in density.quantile.name){
      trait.intra.df.i <-  trait.value.intra.df  %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(density.quantile==n) %>%
        dplyr::select(focal,density.quantile,trait,theoretical.effect,
                      emitter.trait,receiver.trait) %>%
        mutate(focal=as.factor(focal))
      
      glm.intra.trait.i <- glmmTMB(theoretical.effect ~ 1 + receiver.trait  + (1|focal),
                                    trait.intra.df.i,
                                    family="gaussian")
      
      
      Intra.trait.df.i <- as.data.frame(confint(glm.intra.trait.i)) %>%
        mutate(model = "scaled.trait.dist") %>%
        rownames_to_column("parameters") %>%
        rename("Q2.5"="2.5 %",
               "Q97.5"="97.5 %") %>%
        mutate(trait=trait.i,
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
                      theoretical.effect,emitter.trait,receiver.trait,scaled.trait.dist) %>%
        mutate(focal=as.factor(focal),
               neigh=as.factor(neigh))
      
      glm.inter.trait.i.dist <- glmmTMB(theoretical.effect ~  scaled.trait.dist  + (1|focal) + (1|neigh), 
                               trait.dist.df.i,
                               family="gaussian")
      glm.inter.trait.i <- glmmTMB(theoretical.effect ~  receiver.trait  + emitter.trait  + (1|focal) + (1|neigh), 
                                   trait.dist.df.i,
                                   family="gaussian")
  
    
      Inter.trait.df.i <- as.data.frame(confint(glm.inter.trait.i.dist)) %>%
        mutate(model = "scaled.trait.dist") %>%
        rownames_to_column("parameters") %>%
      bind_rows(as.data.frame(confint(glm.inter.trait.i))%>%
                    mutate(model = "trait")%>%
                  rownames_to_column("parameters")) %>%
        rename("Q2.5"="2.5 %",
               "Q97.5"="97.5 %") %>%
        mutate(trait=trait.i,
               density.quantile=n,
               signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                                (Q2.5>0 & Q97.5>0 )~"*",
                                T~""))
      
      Inter.trait.df <- bind_rows(Inter.trait.df,Inter.trait.df.i)
      
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
# for trade off - to add later
trait.dist.df.pivot <- trait.dist.df%>%
  mutate(trait = gsub(" ", ".", trait)) %>%
  dplyr::select(density.quantile,trait, scaled.trait.dist,theoretical.effect)%>%
  spread(trait, scaled.trait.dist)
trait.comb <- combn(gsub(" ", ".", names(trait.df)),2)
trait.comb <- apply(trait.comb,2,paste,collapse=":")
glm.int.all.trait <- glm(formula(paste("theoretical.effect ~ 1",
                                       paste0(gsub(" ", ".", names(trait.df)),
                                              collapse = "+"),
                                       paste0(trait.comb,collapse = "+"),
                                       sep="+")), 
                         data=trait.dist.df.pivot,
                         family='gaussian',
                         control=glm.control(maxit=1000))


Test.all.trait.df <- as.data.frame(summary(glm.int.all.trait)$coef) %>%
  rownames_to_column("parameters") %>%
  rename("estimate"="Estimate",
         "std.error"="Std. Error",
         "t.value"="t value",
         "pvalue"="Pr(>|t|)") %>%
  mutate(trait=trait.i,
         signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                          (pvalue<0.05 & pvalue >0.01 )~"**",
                          (pvalue<0.01)~"***",
                          T~""),
         RMSE = sqrt(mean(glm.int.all.trait$residuals^2)))

#---- 1.2. Make detailed graphs ----
Cool.detailed.theory.trait.plotlist <- list()
density.quantile.name <- c("intercept","low","medium","high")
# apply this later
lm.mod<-function(df){
  m1<-lm(n~year, data=df)
  m2<-lm(n~year+I(year^2), data=df)
  p <- ifelse(AIC(m1)<AIC(m2), "y~x", "y~poly(x, 2)")
  return(p) 
}
country="aus"
for( country in country.list){
  for(n in density.quantile.name){

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
      #dplyr::filter( best.model =="glm.inter.trait.i") %>%
      dplyr::filter(!signif  ==""   | parameters == "(Intercept)") %>%
      dplyr::select(parameters,Estimate, trait,model) %>%
      dplyr::filter(parameters  %in% c("emitter.trait","receiver.trait","scaled.trait.dist","(Intercept)")) %>%
      spread(parameters,Estimate) %>%
      rename("Intercept"="(Intercept)") %>%
      gather(any_of(c("scaled.trait.dist",
               "emitter.trait","receiver.trait")),
             key="trait.param",value="trait.coeff")  %>%
      dplyr::filter(!is.na(trait.coeff))
          

    Intra.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
     dplyr::filter(!signif  =="" | parameters %in% c("(Intercept)")) %>%
      dplyr::select(parameters,Estimate, trait) %>%
      spread(parameters,Estimate) %>% 
      dplyr::filter(!is.na(receiver.trait)) %>%
      rename("trait.coeff"="receiver.trait",
             "Intercept"="(Intercept)") %>%
      mutate(trait= factor(trait ,levels=trait.levels )) 
    
    
    Lambda.trait.df.long.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
      dplyr::filter(density.quantile %in% c(n)) %>%
       dplyr::filter(!signif  =="" | parameters %in% c("(Intercept)")) %>%
      dplyr::select(parameters,Estimate, trait) %>%
      spread(parameters,Estimate) %>% 
      dplyr::filter(!is.na(receiver.trait)) %>%
      rename("trait.coeff"="receiver.trait",
             "Intercept"="(Intercept)") %>%
      mutate(trait= factor(trait ,levels=trait.levels )) 
    
    Inter.trait.df.i <- Cool.theory.trait.df[[country]]$trait.dist.df %>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(theoretical.effect,receiver.trait,trait,emitter.trait,scaled.trait.dist) %>%
      gather(receiver.trait,emitter.trait,scaled.trait.dist,
             key="trait.param", value="trait.value") %>%
      rename("raw.value"="theoretical.effect") %>%
      mutate(parameter="INTER") %>% 
      left_join(Inter.trait.df.long.i) %>%
      dplyr::filter(!is.na(trait.coeff))
    
    Intra.trait.df.i <- Cool.theory.trait.df[[country]]$trait.value.intra.df%>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(theoretical.effect,receiver.trait,trait) %>%
      rename("raw.value"="theoretical.effect",
             "trait.value" ="receiver.trait")%>%
      mutate(parameter="INTRA",
             trait.param="receiver.trait")%>% 
      left_join(Intra.trait.df.long.i) %>%
      dplyr::filter(!is.na(trait.coeff)) 
    
    Lambda.trait.df.i <- Cool.theory.trait.df[[country]]$trait.value.lambda.df %>%
      dplyr::filter( trait %in% levels(as.factor(Lambda.trait.df.long.i$trait))) %>%
      dplyr::filter(density.quantile %in% c(n))  %>%
      dplyr::select(lambda,receiver.trait,trait) %>%
      rename("raw.value"="lambda",
             "trait.value" ="receiver.trait")%>%
      mutate(parameter="lambda",
             trait.param="receiver.trait")%>% 
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
                             mutate(trait.param.label = case_when(trait.param == "emitter.trait" ~ "Trait value of emitter",
                                                                  trait.param == "receiver.trait" ~ "Trait value of focal/receiver",
                                                                  trait.param == "scaled.trait.dist" ~ "Trait value of focal - Trait value of emitter")),
                           aes(y=raw.value,
                               x=trait.value,
                               color=trait),
                           shape=16,
                           size=3,alpha=0.2)+
      geom_abline(data=Inter.trait.df.long.i%>%
                    mutate(trait.param.label = case_when(trait.param == "emitter.trait" ~ "Trait value of emitter",
                                                         trait.param == "receiver.trait" ~ "Trait value of focal/receiver",
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
Cool.glm.theory.trait.plotlist <- list()
density.quantile.name <- c("intercept","low","medium","high")
country="aus"

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
    #dplyr::filter(density.quantile == n) %>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
    #dplyr::filter(!signif  ==""|parameters=="(Intercept)") %>%
    mutate(parameters  = case_when(parameters =="receiver.trait" ~ "Focal trait\non interspecific",
                                   parameters =="emitter.trait" ~ "Emitter trait\non interspecific",
                                   parameters =="scaled.trait.dist" ~ "Focal trait -\nEmitter trait\non interspecific",
                                   parameters =="(Intercept)" ~ "intercept",
                                   T ~ parameters)) %>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))
  
  Intra.trait.df.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
    #dplyr::filter(!signif  =="" ) %>%
    mutate(parameters  = case_when( parameters =="(Intercept)" ~ "intercept",
                                    parameters =="receiver.trait" ~ "Focal trait\non intraspecific",
                                   T ~ parameters))%>%
    mutate(density.quantile= factor(density.quantile,
                                    levels=c("intercept","low","medium","high")))
  
 Lambda.trait.df.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
    dplyr::filter(!parameters %in% c("Std.Dev.(Intercept)|focal","Std.Dev.(Intercept)|neigh")) %>%
   # dplyr::filter(!signif  ==""|parameters=="(Intercept)") %>%
   mutate(parameters  = case_when(parameters =="receiver.trait" ~ "Focal trait\non lambda",
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

 
  Inter.plot.sum <- df.i %>%
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
  
  intercept.df <- Intra.trait.df.i %>%
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
  
  Intra.plot.sum <- Intra.trait.df.i %>%
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
 
  intercept.df <- Lambda.trait.df.i %>%
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

  
   Lambda.plot.sum <- Lambda.trait.df.i %>%
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
                                                          #Lambda.plot.sum,
                                                          axis="tblr",
                                                          align="v",
                                                          ggpubr::get_legend(legend.plot),
                                                          nrow=4,rel_heights = c(1,0.6,0.6,0.3))
 # Cool.glm.theory.trait.plotlist[[country]]  
  
  # JUST INTERCEPT 
  intercept.df <- bind_rows(df.i,Intra.trait.df.i) %>%
    dplyr::filter(density.quantile %in% c("intercept")) %>%
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
    mutate(pointshape = "intercept")

 
  Cool.glm.theory.trait.plotlist[[paste0(country,"_intercept")]]  <- bind_rows(df.i,Intra.trait.df.i) %>%
    dplyr::filter(density.quantile  %in% c("intercept")) %>%
    dplyr::filter(!parameters =="intercept")%>%
    dplyr::filter(!trait %in% c("Root diameter",
                               "Flower width","Seed mass",
                               "Leaf area index",
                               "Canopy shape"))%>%
    mutate(pointshape = "estimate") %>%
    bind_rows(intercept.df) %>%
    mutate(parameters= factor(parameters,
                              levels=c("Focal trait\non interspecific","Emitter trait\non interspecific",
                                       "Focal trait -\nEmitter trait\non interspecific",
                                       "Focal trait\non intraspecific")))%>%
    mutate(y_numb= case_when(parameters =="Focal trait\non interspecific" ~ 4,
                             parameters =="Emitter trait\non interspecific"~3,
                             parameters =="Focal trait -\nEmitter trait\non interspecific" ~2,
                             parameters =="Focal trait\non intraspecific" ~ 1)) %>%
    mutate(trait=factor(trait, 
                        names(dummy.col))) %>%
    mutate(y_trait=((as.numeric(trait)-6)*0.04 +y_numb)) %>%
    ggplot(aes(y=y_trait,
               x=Estimate,
               color=as.factor(trait),
               group=as.factor(trait),
               shape = pointshape )) + 
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
    scale_y_continuous(breaks=c(1:4),
                     labels=rev(c("Focal trait",
                                  "Neighbor trait",
                                  "Focal trait -\nNeighbor trait",
                                  "Focal trait"))) +
    theme_bw() +
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
  Cool.glm.theory.trait.plotlist[[paste0(country,"_intercept")]]
}


Cool.glm.theory.trait.plotlist[[paste0("spain")]]#figures/GLM.traits.SPAIN.pdf
Cool.glm.theory.trait.plotlist[[paste0("aus")]]#figures/GLM.traits.AUS.pdf
legend.plot <- ggplot(data=Cool.theory.trait.df[["aus"]]$Intra.trait.df %>%
                        bind_rows(Cool.theory.trait.df[["spain"]]$Intra.trait.df) %>%
                        mutate(pointshape = rep(c("intercept"),each=240)) %>%
                        dplyr::filter(!parameters =="intercept") %>%
                        dplyr::filter(!trait %in% c("Root diameter",
                                                    "Flower width","Seed mass",
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
                          color=trait,
                          shape = pointshape))+
                          #shape=density.quantile )) + 
  geom_blank() + 
  geom_pointrange(aes(xmin=Q2.5,
                      xmax=Q97.5),
                  size=2,alpha=0.8) +
  scale_shape_manual("",values= 17 )  +
  scale_color_manual("",values= dummy.col )  +
  guides(shape= guide_legend(nrow=2,byrow=TRUE,
                             override.aes = list(alpha = 1,
                                                 color=c("grey80")))) +
  theme_few() +
  theme(legend.position="bottom",
        legend.key.size = unit(1, 'cm'),
        legend.title.position = "top",
        legend.title =element_text(size=20),
        legend.text =element_text(size=16),
        axis.text = element_text(size=18))
legend.plot

GLM.traits.Intercep <- plot_grid( ggarrange(Cool.glm.theory.trait.plotlist[[paste0("aus_intercept")]],
           Cool.glm.theory.trait.plotlist[[paste0("spain_intercept")]],
           common.legend = T, legend = "none",
           label.x = c(0.27,0.15),
           label.y = 1.01,
           font.label = list(size = 26, color = "black", 
                             face = "bold", family = NULL),
           ncol=2,labels=c("a. Australia","b. Spain","")),
           ggpubr::get_legend(legend.plot),
           ncol=1,
          rel_heights =c(1,0.2),
          nrow=2,legend="none")
GLM.traits.Intercep

GLM.traits.Intercep + annotation_custom(
  grob = textGrob(label = "Effect on \ninter- \ninteractions",
                  hjust = 0, gp = gpar(cex = 1.5),rot=0),
  ymin = 0.68,      # Vertical position of the textGrob
  ymax = 0.68,
  xmin = 0.002,         # Note: The grobs are positioned outside the plot area
  xmax = 0.002)+
  annotation_custom(grob = textGrob(label = "Effect on \nintra- \ninteractions",
                    hjust = 0, gp = gpar(cex = 1.5),rot=0),
    ymin = 0.3,      # Vertical position of the textGrob
    ymax = 0.3,
    xmin = 0.002,         # Note: The grobs are positioned outside the plot area
    xmax = 0.002)

bottom_x = 115
grid.brackets(bottom_x,420, bottom_x,40, lwd=2, col="black")
grid.brackets(bottom_x,570, bottom_x,425, lwd=2, col="black")

ggsave(last_plot(),
       width=15.41,
       height=10.30,
       unit="in",
       file="figures/GLM.traits.Intercept.pdf")

library(grid)
library(pBrackets) 
#figures/GLM.traits.Intercept.pdf

#---- 1.4.PCA axes----
PCA.trait.df <- list()
density.quantile.name <- c("intercept","low","medium","high")
data_normalized <- list()
res.traits <- list()
country="aus"
library(FactoMineR)
library(factoextra)
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
      mutate(trait.dist=trait.dist,
             scaled.trait.dist=as.vector(scale(trait.dist)),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list))) %>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>% 
                  rownames_to_column("focal") %>% dplyr::rename("receiver.trait"=i) %>%
                  mutate(receiver.trait = as.vector(scale(receiver.trait))))%>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>%
                  rownames_to_column("neigh") %>% dplyr::rename("emitter.trait"=i)%>%
                  mutate(emitter.trait = as.vector(scale(emitter.trait))))
    
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
    
    
  }
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","SLA","Stem height"))) 
    toremove <- c("Root diameter",
                  "Flower width","Seed mass",
                  "Leaf area index",
                  "Canopy shape","SRA","Root tips")
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","SLA","Stem height"))) 
    toremove <- c("Root diameter",
                  "Flower width","Seed mass",
                  "Leaf area index",
                  "Canopy shape","Root tips")
  }
  #country="aus"
  #country="spain"
 
  
  data_normalized[[country]] <- specific.trait.dist %>% 
    dplyr::select(trait,receiver.trait,focal) %>%
    dplyr::filter(!trait %in% toremove ) %>%
    unique()%>%
    spread(trait,receiver.trait) %>%
    column_to_rownames("focal")
  data_normalized[[country]] 

  #data_normalized[["spain"]] 
  #data_normalized[["aus"]] 
  #https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/
  ifelse(country=="aus",c(2,3,1),c(3,2,1))
  res.traits[[country]] <- MFA(data_normalized[[country]], 
                        group = if(country=="aus"){c(3,2,1)}else{c(2,3,1)},#(1,1,1,1,1,1,1,1,1,1), 
                        type = c("s","s","s"),#c("n","n","n","s","s","s","s","s","n","n"),
                        name.group = c("Root cooperation spectrum",
                                       "Drought sensitivity",
                                       "Plant size"),
                        graph = T)
  
  #FactoMineR::PCA(data_normalized[[country]])
  
  #country="aus"
  #country="spain"
  eig.val <- get_eigenvalue(res.traits[[country]])
  head(eig.val)
  pdf(paste0("figures/pca/",country,"PCAaxes.pdf"))
  chart.Correlation(  data_normalized[[country]], histogram=TRUE, pch=19)
  
  fviz_contrib(res.traits[[country]],
               choice = "quanti.var", axes = 1, top = 20,palette = "jco")
  fviz_contrib(res.traits[[country]], choice = "quanti.var",
               axes = 2, top = 20,palette = "jco")
  fviz_contrib(res.traits[[country]], choice = "quanti.var", axes = 3, top = 20,palette = "jco")
  fviz_mfa_var(res.traits[[country]], "quanti.var", palette = "jco", 
               col.var.sup = "violet", repel = TRUE)
  
  fviz_mfa_ind(res.traits[[country]], col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE)
  dev.off()
  # cluster analysis
  library(vegan)
  if(country=="spain"){
  PCA.trait.df[[country]]  <-  data.frame("BelowGroundStrategy"=res.traits[[country]]$ind$coord[,1],
                                          "PlantSize"=res.traits[[country]]$ind$coord[,2],
                                          "AboveGroundStrategy"=res.traits[[country]]$ind$coord[,3],
                                          "BelowGroundTrait"=data_normalized[[country]]$SRA,
                                          "PlantSizeTrait"=data_normalized[[country]]$`Stem height`,
                                          "AboveGroundTrait"=data_normalized[[country]]$SLA)
  }
  
  if(country=="aus"){
    PCA.trait.df[[country]]  <-  data.frame("PlantSize"=-1*res.traits[[country]]$ind$coord[,1],
                                            "AboveGroundStrategy"=-1*res.traits[[country]]$ind$coord[,2],
                                            "BelowGroundStrategy"=res.traits[[country]]$ind$coord[,3],
                                            "BelowGroundTrait"=data_normalized[[country]]$`Root length`,
                                            "PlantSizeTrait"=data_normalized[[country]]$`Stem height`,
                                            "AboveGroundTrait"=data_normalized[[country]]$`C13 water use efficiency`)
  }
}

#---- 1.5.PCA GLM- FOCAL----
glm.trait.pca.list <- list()
glm.trait.list <- list()
glm.new.data.list <- list()
glm.new.data.trait.list<- list()
for( country in country.list){
  
trait.pca.df <- Theoretical.Int.list[[country]] %>%
  dplyr::filter(!neigh ==focal) %>% 
  left_join(as.data.frame(PCA.trait.df[[country]] %>%
                            rownames_to_column("focal")),
            relationship ="many-to-many")

glm.trait.pca.df <- NULL
glm.trait.df <- NULL
glm.new.data.df <- NULL
glm.new.data.trait.df <- NULL
#glm.new.data.df.3D <- NULL
#glm.new.data.df.3D.i <- data.frame(AboveGroundStrategy = c(0,seq(min(trait.pca.df$AboveGroundStrategy),max(trait.pca.df$AboveGroundStrategy),by=(max(trait.pca.df$AboveGroundStrategy)-min(trait.pca.df$AboveGroundStrategy))/100)),
#                                   BelowGroundStrategy = c(0,seq(min(trait.pca.df$BelowGroundStrategy),max(trait.pca.df$BelowGroundStrategy),by=(max(trait.pca.df$BelowGroundStrategy)-min(trait.pca.df$BelowGroundStrategy))/100)))
#glm.new.data.df.3D.i <- expand.grid(glm.new.data.df.3D.i)
#glm.new.data.df.3D.i$PlantSize <- glm.new.data.df.3D.i$AboveGroundStrategy + glm.new.data.df.3D.i$BelowGroundStrategy
#glm.new.data.df.3D.i <- glm.new.data.df.3D.i %>%
#  dplyr::filter(PlantSize>min(trait.pca.df$PlantSize) &
 #          PlantSize < max(trait.pca.df$PlantSize))
#glm.new.data.df.i <- data.frame(PlantSize=c(0,seq(min(trait.pca.df$PlantSize),max(trait.pca.df$PlantSize),by=(max(trait.pca.df$PlantSize)-min(trait.pca.df$PlantSize))/100)),# interaction Size-Above
#                                AboveGroundStrategy = c(0,seq(min(trait.pca.df$AboveGroundStrategy),max(trait.pca.df$AboveGroundStrategy),by=(max(trait.pca.df$AboveGroundStrategy)-min(trait.pca.df$AboveGroundStrategy))/100)),
#                                BelowGroundStrategy = c(0,seq(min(trait.pca.df$BelowGroundStrategy),max(trait.pca.df$BelowGroundStrategy),by=(max(trait.pca.df$BelowGroundStrategy)-min(trait.pca.df$BelowGroundStrategy))/100)))

#glm.new.data.df.i <- expand.grid(glm.new.data.df.i)

glm.new.data.df.i <- data.frame(PlantSize=c(runif(50000,min = min(trait.pca.df$PlantSize), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                            #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                            runif(50000,min = 0, max = max(trait.pca.df$PlantSize))),# interaction Size-Above # quantile(trait.pca.df$PlantSize,0.66)
                                AboveGroundStrategy = runif(100000, min =min(trait.pca.df$AboveGroundStrategy),
                                                            max=max(trait.pca.df$AboveGroundStrategy)),
                                BelowGroundStrategy = runif(100000, min =min(trait.pca.df$BelowGroundStrategy),
                                                            max=max(trait.pca.df$BelowGroundStrategy)),
                                PlantSizeCat =rep(c("BelowAverage","AboveAverage"),each=50000),
                                BelowGroundTrait=runif(100000, min =min(trait.pca.df$BelowGroundTrait,na.rm=T),
                                                       max=max(trait.pca.df$BelowGroundTrait,na.rm=T)),
                                PlantSizeTrait=c(runif(50000,min = min(trait.pca.df$PlantSizeTrait,na.rm=T), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                                 #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                                 runif(50000,min = 0, max = max(trait.pca.df$PlantSizeTrait))),
                                AboveGroundTrait=runif(100000, min =min(trait.pca.df$AboveGroundTrait,na.rm=T),
                                                       max=max(trait.pca.df$AboveGroundTrait,na.rm=T))) #"AroundAverage",
glm.new.data.df.i <- glm.new.data.df.i %>%
  dplyr::filter(PlantSize>min(trait.pca.df$PlantSize) &
          PlantSize < max(trait.pca.df$PlantSize))


for(n in density.quantile.name){
    trait.pca.df.i <-  trait.pca.df %>%
      dplyr::filter(density.quantile==n) %>%
      dplyr::select(focal,theoretical.effect,density.quantile,
                    PlantSize,AboveGroundStrategy,BelowGroundStrategy,
                    BelowGroundTrait,PlantSizeTrait,AboveGroundTrait) %>%
      mutate(focal=as.factor(focal))
    
    glm.trait.pca.i <- glmmTMB(theoretical.effect ~  PlantSize*AboveGroundStrategy*BelowGroundStrategy, 
                                      trait.pca.df.i,
                                      family="gaussian")
    summary(glm.trait.pca.i)
    glm.trait.pca.df.i <- as.data.frame(confint(glm.trait.pca.i)) %>%
      mutate(density.quantile=n) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
             "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.pca.df <- bind_rows(glm.trait.pca.df,glm.trait.pca.df.i)

    
    glm.new.data.df.i$theoretical.effect.pca <- predict(glm.trait.pca.i,
                                                glm.new.data.df.i)
    
    glm.trait.i <- glmmTMB(theoretical.effect ~  PlantSizeTrait*AboveGroundTrait*BelowGroundTrait, 
                               trait.pca.df.i,
                               family="gaussian")
    summary(glm.trait.i)
    glm.trait.df.i <- as.data.frame(confint(glm.trait.i)) %>%
      mutate(density.quantile=n) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
                    "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.df <- bind_rows(glm.trait.df,glm.trait.df.i)
    
    
    glm.new.data.df.i$theoretical.effect <- predict(glm.trait.i,
                                                    glm.new.data.df.i)
    glm.new.data.df <- bind_rows( glm.new.data.df, glm.new.data.df.i%>%
                                    mutate(density.quantile=n))

    
}
glm.trait.pca.list[[country]] <-  glm.trait.pca.df
glm.trait.list[[country]] <-  glm.trait.df
glm.new.data.list[[country]] <-glm.new.data.df
glm.new.data.trait.list[[country]] <- glm.new.data.trait.df
}
view(glm.trait.pca.list[["aus"]])  
view(glm.trait.pca.list[["spain"]])

#---- 1.6.PCA plot glm-FOCAL----
# Above ground and below ground
glm.plot.list <- list()
# for focal species
for( country in country.list){
  if(country=="aus"){
    limits.vec = c(-0.4,0.4)
    ylabname = "13C, WUE"
    xlabname = "Root length"
  }else{limits.vec = c(-0.4,0.4)
  ylabname = "SLA - Inversed 13C, WUE"
  xlabname = "SRA"}
    focal.pca.df.i <- Theoretical.Int.list[[country]] %>%
      dplyr::filter(!neigh ==focal) %>% 
      aggregate(theoretical.effect~ focal, median) %>%
      left_join(as.data.frame(PCA.trait.df[[country]] %>%
                                rownames_to_column("focal")))%>%
      mutate(PlantSizeCat = case_when(PlantSizeTrait < 0 ~ "BelowAverage", #quantile(trait.pca.df$PlantSize,0.33)
                                      PlantSizeTrait > 0 ~ "AboveAverage", #quantile(trait.pca.df$PlantSize,0.66)
                                      T~"AroundAverage"))
    
    glm.new.data.df.i <- glm.new.data.list[[country]]%>%
      filter(theoretical.effect  > min(Theoretical.Int.list[[country]]$theoretical.effect) &
               theoretical.effect < max(Theoretical.Int.list[[country]]$theoretical.effect))
glm.plot.list[[country]] <-glm.new.data.df.i  %>%
  dplyr::filter(density.quantile %in%c("intercept"))%>%
  ggplot(aes(y=AboveGroundTrait, #Strategy,
             x=BelowGroundTrait, #Strategy,
             fill=theoretical.effect)) +
  geom_tile(width = 0.04,
              height = 0.04,
            alpha=0.8)+
  geom_point(data=focal.pca.df.i,
            aes(y=AboveGroundTrait, #Strategy,
                x=BelowGroundTrait, #Strategy,
                fill=theoretical.effect),
            shape=21,size=3) +
  geom_text(data=focal.pca.df.i,
             aes(y=AboveGroundTrait, #Strategy, 
                 x=BelowGroundTrait, #Strategy,
                 label=focal),
            size=4,hjust=0.5,vjust=-0.5) +
  scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                             101, 
                                             type = "continuous")),
                       limits=limits.vec)+
  facet_wrap(density.quantile~PlantSizeCat,nrow=1)+
  theme_clean()+
  labs(title=country,
       y =   ylabname, #"Above ground strategy, \nfrom less to more drought sensitive",
       x =   xlabname, #"Below ground strategy, \nfrom more to less cooperative",
       subtitle = "Plant height", #"Plant size, from high/fast to small/slow growing",
       fill="Predicted interaction effect on \nof focal species with that trait combinaison")+
  theme(legend.position="bottom",
        legend.key.size = unit(1, 'cm'),
        legend.background = element_blank(),
        plot.background = element_blank(),
        axis.text =element_text(size=20),
        axis.title =element_text(size=20),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title.position = "top",
        legend.title =element_text(size=20),
        legend.text =element_text(size=20))

glm.plot.list[[country]]

 
}

glm.plot.list[["spain"]] # figures/pca/Flat3D.focal.spain.pdf
glm.plot.list[["aus"]] # figures/pca/Flat3D.focal.aus.pdf

ggarrange(glm.plot.list[["spain"]],
          glm.plot.list[["aus"]],
          common.legend = T,legend="bottom",
          nrow=2, 
          labels=c("a. Spain","b. Australia"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
# figures/pca/Flat3D.focal.intercept.pdf
levels(as.factor(Chapt1_Parameters_values$parameter))
names(Chapt1_Parameters_values)
ParameterValuesForJurg <- Chapt1_Parameters_values %>%
  filter(parameter %in% c("Intraspecific","Plant - plant") &
           year== "2021" &
           complexity.plant =="family") %>%
  dplyr::rename("Iterations"="X" ) %>%
  dplyr::select(Iterations,focal,lambdas,parameter,estimate) %>%
  spread(parameter,estimate)
  view(ParameterValuesForJurg)
  write.csv(ParameterValuesForJurg,
            file="/Users/lisabuche/Downloads/ParameterValuesForJurg")
#---- 1.7.PCA  GLM-NEIGH----
glm.trait.pca.neigh.list <- list()
glm.trait.neigh.list <- list()
glm.new.data.neigh.list <- list()
glm.new.data.neigh.3D.list<- list()
for( country in country.list){
  
  trait.pca.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    left_join(as.data.frame(PCA.trait.df[[country]] %>%
                              rownames_to_column("neigh")),
              relationship ="many-to-many")
  
  
  glm.new.data.df.i <- data.frame(PlantSize=c(runif(50000,min = min(trait.pca.df$PlantSize), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                              #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                              runif(50000,min = 0, max = max(trait.pca.df$PlantSize))),# interaction Size-Above # quantile(trait.pca.df$PlantSize,0.66)
                                  AboveGroundStrategy = runif(100000, min =min(trait.pca.df$AboveGroundStrategy),
                                                              max=max(trait.pca.df$AboveGroundStrategy)),
                                  BelowGroundStrategy = runif(100000, min =min(trait.pca.df$BelowGroundStrategy),
                                                              max=max(trait.pca.df$BelowGroundStrategy)),
                                  PlantSizeCat =rep(c("BelowAverage","AboveAverage"),each=50000),
                                  BelowGroundTrait=runif(100000, min =min(trait.pca.df$BelowGroundTrait,na.rm=T),
                                                         max=max(trait.pca.df$BelowGroundTrait,na.rm=T)),
                                  PlantSizeTrait=c(runif(50000,min = min(trait.pca.df$PlantSizeTrait,na.rm=T), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                                   #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                                   runif(50000,min = 0, max = max(trait.pca.df$PlantSizeTrait))),
                                  AboveGroundTrait=runif(100000, min =min(trait.pca.df$AboveGroundTrait,na.rm=T),
                                                         max=max(trait.pca.df$AboveGroundTrait,na.rm=T))) #"AroundAverage",
  glm.new.data.df.i <- glm.new.data.df.i %>%
    dplyr::filter(PlantSize>min(trait.pca.df$PlantSize) &
                    PlantSize < max(trait.pca.df$PlantSize))
  
  glm.trait.df <- NULL
  glm.trait.pca.df<- NULL
  glm.new.data.df <- NULL
  for(n in density.quantile.name){
    trait.pca.df.i <-  trait.pca.df %>%
      dplyr::filter(density.quantile==n) %>%
      dplyr::select(neigh,theoretical.effect,density.quantile,
                    PlantSize,AboveGroundStrategy,BelowGroundStrategy,
                    PlantSizeTrait,AboveGroundTrait,BelowGroundTrait) %>%
      mutate(neigh=as.factor(neigh))
    
    glm.trait.pca.i <- glmmTMB(theoretical.effect ~  PlantSize*AboveGroundStrategy*BelowGroundStrategy, 
                               trait.pca.df.i,
                               family="gaussian")
    summary(glm.trait.pca.i)
    glm.trait.pca.df.i <- as.data.frame(confint(glm.trait.pca.i)) %>%
      mutate(density.quantile=n) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
                    "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.pca.df <- bind_rows(glm.trait.pca.df,glm.trait.pca.df.i)
    
    
    glm.new.data.df.i$theoretical.effect.pca <- predict(glm.trait.pca.i,
                                                    glm.new.data.df.i)
    glm.trait.i <- glmmTMB(theoretical.effect ~  PlantSizeTrait*AboveGroundTrait*BelowGroundTrait, 
                           trait.pca.df.i,
                           family="gaussian")
    summary(glm.trait.i)
    glm.trait.df.i <- as.data.frame(confint(glm.trait.i)) %>%
      mutate(density.quantile=n) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
                    "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.df <- bind_rows(glm.trait.df,glm.trait.df.i)
    
    
    glm.new.data.df.i$theoretical.effect <- predict(glm.trait.i,
                                                    glm.new.data.df.i)
    glm.new.data.df <- bind_rows( glm.new.data.df, glm.new.data.df.i%>%
                                    mutate(density.quantile=n))
  }
  glm.trait.pca.neigh.list[[country]] <-  glm.trait.pca.df
  glm.trait.neigh.list[[country]] <-  glm.trait.df
  glm.new.data.neigh.list[[country]] <-glm.new.data.df
}
#---- 1.8.PCA plot glm-NEIGH----
glm.plot.neigh.list <- list()
for( country in country.list){
  if(country=="aus"){
    limits.vec = c(-0.4,0.4)
    ylabname = "13C, WUE"
    xlabname = "Root length"
  }else{limits.vec = c(-0.4,0.4)
  ylabname = "SLA - Inversed 13C, WUE"
  xlabname = "SRA"}
  neigh.pca.df.i <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    aggregate(theoretical.effect~ neigh , median) %>%
    left_join(as.data.frame(PCA.trait.df[[country]] %>%
                              rownames_to_column("neigh")))%>%
    mutate(PlantSizeCat = case_when(PlantSizeTrait < 0 ~ "BelowAverage",
                                    PlantSizeTrait > 0 ~ "AboveAverage",
                                    T~"AroundAverage"))
  
  glm.new.data.df.i <- glm.new.data.neigh.list[[country]]%>%
    filter(theoretical.effect  > min(Theoretical.Int.list[[country]]$theoretical.effect) &
             theoretical.effect < max(Theoretical.Int.list[[country]]$theoretical.effect))
  
  glm.plot.neigh.list[[country]] <-glm.new.data.df.i  %>%
    dplyr::filter(density.quantile %in%c("intercept"))%>%
    ggplot(aes(y=AboveGroundTrait, #Strategy, 
               x=BelowGroundTrait, #Strategy, 
               fill=theoretical.effect)) +
    geom_tile(width = 0.04,
              height = 0.04,
              alpha=0.8)+
    geom_point(data= neigh.pca.df.i,
               aes(y=AboveGroundTrait, #Strategy, 
                   x=BelowGroundTrait, #Strategy, 
                   fill=theoretical.effect),
               shape=21,size=3) +
    geom_text(data= neigh.pca.df.i,
              aes(y=AboveGroundTrait, #Strategy,  
                  x=BelowGroundTrait, #Strategy, 
                  label=neigh),
              size=4,hjust=0.5,vjust=-0.5) +
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                   101, 
                                                   type = "continuous")),
                         limits=limits.vec)+
    facet_wrap(density.quantile~PlantSizeCat,nrow=1)+
    theme_clean()+
    labs(title=country,
         y = ylabname, #"Above ground strategy, \nfrom less to more drought sensitive",
         x =  xlabname, #"Below ground strategy, \nfrom more to less cooperative",
         subtitle = "Plant height",#"Plant size, from high/fast to small/slow growing",
         fill="Predicted interaction effect on \nof focal species with that trait combinaison")+
    theme(legend.position="bottom",
          legend.key.size = unit(2, 'cm'),
          legend.background = element_blank(),
          plot.background = element_blank(),
          axis.text =element_text(size=20),
          axis.title =element_text(size=20),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title.position = "top",
          legend.title =element_text(size=20),
          legend.text =element_text(size=20))
  
  glm.plot.neigh.list[[country]]
  
  
}
glm.plot.neigh.list[["aus"]]# figures/pca/Flat3D.neigh.aus.pdf
glm.plot.neigh.list[["spain"]]# figures/pca/Flat3D.neigh.spain.pdf

ggarrange(glm.plot.neigh.list[["spain"]],
          glm.plot.neigh.list[["aus"]],
          common.legend = T,legend="bottom",
          nrow=2, 
          #labels=c("a. Spain","b. Australia"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
# figures/pca/Flat3D.neigh.intercept.pdf

#----Not use----
data.3d <- glm.new.data.3D.list[[country]]%>%
  dplyr::filter(density.quantile=="intercept") %>%
  dplyr::select(-density.quantile) %>%
  as.matrix()%>%
  as.data.frame()%>%
  mutate_all(function(x) as.numeric(x))

str(data.3d)
ThreeDfig <- plot_ly(x=~PlantSize, y= ~AboveGroundStrategy, 
                     z=~BelowGroundStrategy, data=data.3d,
                     #type="scatter3d", mode="markers", 
                     #color=~theoretical.effect)
                     type='mesh3d',
                     flatshading=T,
                     intensity =~theoretical.effect,
                     colors=colorRamp(rev(wes_palette("Zissou1", 
                                                      round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*200)+1, 
                                                      type = "continuous")))) %>%
  add_trace(x=~PlantSize, y= ~AboveGroundStrategy, 
            z=~BelowGroundStrategy, data=focal.pca.df.i,
            color= ~theoretical.effect,mode = "markers",
            type = "scatter3d",
            marker = list(size = 5),
            colors=colorRamp(rev(wes_palette("Zissou1", 
                                             round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*200)+1, 
                                             type = "continuous"))))
ThreeDfig
library(scatterplot3d)
focal.pca.df.i <-  PCA.trait.df[[country]] %>%
  rownames_to_column("focal")

colors.gradient <- rev(wes_palette("Zissou1", 
                                   round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*200)+1, 
                                   type = "continuous"))
focal.colors.gradient <- colors.gradient[(round(focal.pca.df.i$theoretical.effect*100)+round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*100)+1)]
colors.gradient <- colors.gradient[(round(data.3d[,4]*100)+round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*100)+1)]


s3d <- scatterplot3d(data.3d[,1:3], angle = 60,
                     pch = 16, color=colors.gradient,
                     grid=TRUE,
                     main=country,
                     xlab = "Above ground strategy, from less to more drought sensitive",
                     ylab = "Below ground strategy, from more to less cooperative",
                     zlab = "Plant size, from small/slow to high/fast growing")

s3d$points3d(focal.pca.df.i$AboveGroundStrategy,
             focal.pca.df.i$BelowGroundStrategy,
             focal.pca.df.i$PlantSize, 
             col = focal.colors.gradient , type = "h", pch = 8)

text(s3d$xyz.convert(focal.pca.df.i$AboveGroundStrategy,
                     focal.pca.df.i$BelowGroundStrategy,
                     focal.pca.df.i$PlantSize)
     labels=focal.pca.df.i$focal,
     col='tomato', cex=2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Supp figures----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 2.1. Trait Density figures----

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

