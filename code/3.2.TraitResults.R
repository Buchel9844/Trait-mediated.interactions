
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
library(lavaan)
library(emmeans)
library(piecewiseSEM)
library(lmtest)
#setwd("/home/lbuche/Eco_Bayesian/chapt3")
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
project.dic <- ""
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"

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
        dplyr::filter(density.quantile==n) 
      
      rsq <- function (x, y) cor(x, y,use="complete.obs") ^ 2
      glm.lambda.trait.i <- glm(lambda ~ 1 + receiver.trait,
                                trait.lambda.df.i ,
                                family="gaussian")
      summary(glm.lambda.trait.i)
      #rsq(trait.lambda.df.i$lambda,trait.lambda.df.i$receiver.trait)
      
      trait.lambda.df.i<-  trait.lambda.df.i[complete.cases( trait.lambda.df.i),]
      glm.lambda.trait.zero.i <- glm(lambda~ 1,
                                     trait.lambda.df.i,
                                     na.action = na.omit,
                                     family="gaussian")
      
      
      glm.lambda.trait.summary.i <-  as.data.frame(lrtest( glm.lambda.trait.zero.i,
                                                           glm.lambda.trait.i)) %>%
        mutate(model.comp =c("null","null vs linear"),
               model=c("glm.lambda.trait.zero.i","glm.lambda.trait.line.i")) %>%
        dplyr::filter(!model.comp=="H0") %>%
        rename("pvalue"="Pr(>Chisq)") %>%
        mutate(signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~" "),
               trait=trait.i,
               density.quantile=n)  
      glm.intra.trait.summary <- bind_rows( glm.lambda.trait.summary, glm.lambda.trait.summary.i)
      
      if(glm.lambda.trait.summary.i$signif[which(glm.lambda.trait.summary.i$model.comp=="null vs linear")]==" "){
        model.to.keep <- "glm.lambda.trait.zero.i" 
      }else{
        model.to.keep <-  "glm.lambda.trait.i" 
      }
      
      Lambda.trait.df.i <- as.data.frame(summary(get(model.to.keep))$coef) %>%
        rownames_to_column("parameters") %>%
        rename("estimate"="Estimate",
               "std.error"="Std. Error",
               "t.value"="t value",
               "pvalue"="Pr(>|t|)") %>%
        mutate(trait=trait.i,
               best.model =model.to.keep,
               density.quantile=n,
               signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~""),
               RMSE = sqrt(mean(get(model.to.keep)$residuals^2)))
      
      Lambda.trait.df <- bind_rows( Lambda.trait.df, Lambda.trait.df.i)
      
    }
  }
  # Make data frae with trait and INTRA specific interactions
  trait.value.intra.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(neigh ==focal) %>% # only INTRA
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
        dplyr::select(neigh,density.quantile,trait,theoretical.effect,
                      emitter.trait,receiver.trait)
      
      rsq <- function (x, y) cor(x, y,use="complete.obs") ^ 2
      glm.intra.trait.i <- glm(theoretical.effect ~ 1 + receiver.trait,
                               trait.intra.df.i ,
                               family="gaussian")
      #summary(glm.intra.trait.i)
      
      trait.intra.df.i<- trait.intra.df.i[complete.cases(trait.intra.df.i),]
      glm.intra.trait.zero.i <- glm(theoretical.effect ~ 1,
                                    trait.intra.df.i,
                                    na.action = na.omit,
                                    family="gaussian")
      
      
      glm.intra.trait.summary.i <-  as.data.frame(lrtest( glm.intra.trait.zero.i,
                                                          glm.intra.trait.i)) %>%
        mutate(model.comp =c("null","null vs linear"),
               model=c("glm.inter.trait.zero.i","glm.inter.trait.line.i")) %>%
        dplyr::filter(!model.comp=="H0") %>%
        rename("pvalue"="Pr(>Chisq)") %>%
        mutate(signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~" "),
               trait=trait.i,
               density.quantile=n)  
      glm.intra.trait.summary <- bind_rows( glm.intra.trait.summary, glm.intra.trait.summary.i)
      
      if(glm.intra.trait.summary.i$signif[which(glm.intra.trait.summary.i$model.comp=="null vs linear")]==" "){
        model.to.keep <- "glm.intra.trait.zero.i" 
      }else{
        model.to.keep <-  "glm.intra.trait.i" 
      }
      
      Intra.trait.df.i <- as.data.frame(summary(get(model.to.keep))$coef) %>%
        rownames_to_column("parameters") %>%
        rename("estimate"="Estimate",
               "std.error"="Std. Error",
               "t.value"="t value",
               "pvalue"="Pr(>|t|)") %>%
        mutate(trait=trait.i,
               best.model =model.to.keep,
               density.quantile=n,
               signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~""),
               RMSE = sqrt(mean(get(model.to.keep)$residuals^2)))
      
      Intra.trait.df <- bind_rows(Intra.trait.df,Intra.trait.df.i)
      
    }
  }
  # Make data frame with trait and inter specific interactions
  trait.dist.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% # only inter
    left_join(specific.trait.dist,
              relationship ="many-to-many")
  Inter.trait.df <- NULL
  glm.inter.trait.summary<- NULL
  for( trait.i in names(trait.df)){
    for(n in density.quantile.name){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(density.quantile==n) %>%
        dplyr::select(neigh,density.quantile,trait,theoretical.effect,emitter.trait,receiver.trait,scaled.trait.dist)
      
      glm.inter.trait.i <- glm(theoretical.effect ~ emitter.trait + receiver.trait + abs(scaled.trait.dist), trait.dist.df.i,
                               family="gaussian")
      #summary(glm.inter.trait.i)
      Inter.trait.df.i <- as.data.frame(summary(glm.inter.trait.i)$coef) %>%
        rownames_to_column("parameters") %>%
        rename("estimate"="Estimate",
               "std.error"="Std. Error",
               "t.value"="t value",
               "pvalue"="Pr(>|t|)") %>%
        mutate(trait=trait.i,
               signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~""))
      
      variable.to.keep <- Inter.trait.df.i$parameters[which(!Inter.trait.df.i$signif=="")]
      variable.to.keep <- variable.to.keep[!variable.to.keep == "(Intercept)"]
      if(!length(variable.to.keep) == 0){
        trait.dist.df.i <- trait.dist.df.i[complete.cases(trait.dist.df.i),]
        glm.inter.trait.zero.i <- glm(theoretical.effect ~ 1, trait.dist.df.i,
                                      na.action = na.omit,
                                      family="gaussian")
        
        glm.inter.trait.line.i <- glm(formula(paste("theoretical.effect ~ ",paste0(variable.to.keep, collapse="+"))), 
                                      trait.dist.df.i,
                                      na.action = na.omit,
                                      family="gaussian")
        
        glm.inter.trait.quad.i <- glm(formula(paste0("theoretical.effect ~ ",paste0(variable.to.keep, collapse="+"),"+",
                                                     paste0(paste0("I(",variable.to.keep,"^2)"), collapse="+"))),
                                      trait.dist.df.i,
                                      family="gaussian")
        glm.inter.trait.summary.i <- bind_rows(
          as.data.frame(lrtest( glm.inter.trait.zero.i,glm.inter.trait.line.i)) %>%
            mutate(model.comp =c("null","null vs linear"),
                   model=c("glm.inter.trait.zero.i","glm.inter.trait.line.i")),
          as.data.frame(lrtest(  glm.inter.trait.zero.i,glm.inter.trait.quad.i))%>%
            mutate(model.comp =c("H0","null vs quadratic"),
                   model=c("glm.inter.trait.zero.i","glm.inter.trait.quad.i")),
          as.data.frame(lrtest(  glm.inter.trait.line.i,glm.inter.trait.quad.i))%>%
            mutate(model.comp =c("H0","linear vs quadratic"),
                   model=c("glm.inter.trait.line","glm.inter.trait.quad.i"))) %>%
          dplyr::filter(!model.comp=="H0") %>%
          rename("pvalue"="Pr(>Chisq)") %>%
          mutate(signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                  (pvalue<0.05 & pvalue >0.01 )~"**",
                                  (pvalue<0.01)~"***",
                                  T~" "),
                 trait=trait.i,
                 density.quantile=n)  
        glm.inter.trait.summary <- bind_rows( glm.inter.trait.summary,
                                              glm.inter.trait.summary.i)
        
        if(glm.inter.trait.summary.i$signif[which(glm.inter.trait.summary.i$model.comp=="null vs linear")]==" "){
          model.to.keep <- "glm.inter.trait.zero.i" 
        }else{model.to.keep <- "glm.inter.trait.i"
        }
        best.model = model.to.keep }else{
          model.to.keep <- "glm.inter.trait.i"
          best.model <- "no significant variables"
        }
      Inter.trait.df.i <- as.data.frame(summary(get(model.to.keep))$coef) %>%
        rownames_to_column("parameters") %>%
        rename("estimate"="Estimate",
               "std.error"="Std. Error",
               "t.value"="t value",
               "pvalue"="Pr(>|t|)") %>%
        mutate(trait=trait.i,
               best.model =best.model,
               density.quantile=n,
               signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~""),
               RMSE = sqrt(mean(get(model.to.keep)$residuals^2)))
      
      Inter.trait.df <- bind_rows(Inter.trait.df,Inter.trait.df.i)
      
    }
  }
  
  Cool.theory.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    Inter.trait.df=Inter.trait.df,
    Lambda.trait.df=Lambda.trait.df,
    Intra.trait.df =Intra.trait.df,
    glm.inter.trait.summary=glm.inter.trait.summary)
}
Cool.theory.trait.df[[country]]$Inter.trait.df
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

for( country in country.list){
  for(n in density.quantile.name){
    
    # Change to have quadratic too 
    
    Cool.detailed.theory.trait.plotlist[[paste0(country,"_",n)]] <- Cool.theory.trait.df[[country]]$trait.dist.df %>%
      dplyr::filter(density.quantile == n) %>%
      ggplot(aes(y=theoretical.effect,
                 x=scaled.trait.dist,
                 fill=theoretical.effect)) +
      geom_point(shape=21,color="white",
                 size=1,alpha=1) + 
      geom_hline(yintercept=0,color="black") +
      geom_smooth(method="lm",color="black",se=F) + 
      geom_smooth(method="glm",formula="y ~poly(x, 2)",color="grey",se=F) + 
      facet_wrap(trait~.,scale="free") + 
      labs(x="Focal trait value - Neigh trait value",
           y=paste0("Theoretical interaction at ",n," density"))+
      scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                     101, 
                                                     type = "continuous")))+
      theme_bw() 
    Cool.detailed.theory.trait.plotlist[[paste0(country,"_",n)]]
    
  }
}
Cool.detailed.theory.trait.plotlist[[paste0("spain_intercept")]] # figures/Theory.int.intercept.spain.pdf
Cool.detailed.theory.trait.plotlist[[paste0("aus_intercept")]]

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
  }
  if(country=="spain"){
    Inter.trait.df<- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  
  df.i <- Inter.trait.df %>%
    #dplyr::filter(density.quantile == n) %>%
    dplyr::filter(!parameters %in% c("(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    dplyr::filter(!signif  =="") %>%
    #mutate(parameters=factor(parameters,
    #                               levels=c("receiver.trait","emitter.trait","I(receiver.trait^2)","I(emitter.trait^2)"))) %>%
    mutate(n = as.numeric(as.factor(parameters))) %>%
    mutate(traitnum  = as.numeric(as.factor(trait))) %>%
    # Add scaling factor by n group 
    mutate(y_num = traitnum + case_when(n == 1 ~  0,
                                        n == 2 ~ -0.1,
                                        n == 3 ~   0.1,
                                        n == 4 ~   0.2,
                                        n == 5 ~   -0.2)) %>%
    mutate(density.quantile=factor(density.quantile,
                                   levels=density.quantile.name))
  
  Cool.glm.theory.trait.plotlist[[country]] <- df.i %>%
    ggplot(aes(y=y_num,
               x=estimate,
               shape=as.factor(parameters),
               color=density.quantile)) +
    geom_pointrange(aes(xmin=estimate - std.error,
                        xmax=estimate + std.error),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.2)) +
    scale_y_continuous(breaks = unique(df.i$traitnum), 
                       labels = unique(df.i$trait),
                       limits=c(0,11.4)) +
    scale_color_colorblind()+
    scale_shape_manual(values=c(15:18)) +
    theme_bw() +
    geom_vline(xintercept=0) + 
    labs(x="estimate",
         y=paste0(""),
         shape="parameter",
         color="Density of emitter") +
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           shape = guide_legend(title.position = "top",
                                nrow=2)) +
    annotate("text", x = -0.022, y=0.2, 
             label = "Increase in competition",
             size=4) + 
    geom_segment(aes(x = 0, y = 0, xend = 0.015, yend = 0),
                 arrow = arrow(length = unit(0.4, "cm")), color="black")+
    annotate("text", x = 0.022, y=0.2, 
             label = "Increase in facilitation",
             size=4) + 
    geom_segment(aes(x = 0, y = 0, xend = -0.015, yend = 0),
                 arrow = arrow(length = unit(0.4, "cm")),color="black") +
    
    theme(legend.position="bottom",
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=16),
          axis.text = element_text(size=14),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
  
}
Cool.glm.theory.trait.plotlist[[paste0("spain")]]
Cool.glm.theory.trait.plotlist[[paste0("aus")]]
ggarrange(Cool.glm.theory.trait.plotlist[[paste0("aus")]],
          Cool.glm.theory.trait.plotlist[[paste0("spain")]],
          nrow=1,legend="bottom",label.x = 0.2,
          align=c("h"),
          common.legend = T, 
          labels=c("a. Australia","b. Spain"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
#figures/GLM.traits.pdf
#---- 1.5. Make rec of intercept/lambda and intra ----
Rect.glm.theory.trait.plotlist <- list()
density.quantile.name <- c("intercept","low","medium","high")
country="aus"
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
  if(country=="aus"){
    Inter.trait.df <- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  if(country=="spain"){
    Inter.trait.df<- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  Inter.trait.df.i  <- Inter.trait.df %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    dplyr::filter(!parameters %in% c("(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    #dplyr::filter(!signif  =="") %>%
    dplyr::filter( best.model =="glm.inter.trait.i") %>%
    mutate(trait=addline_format(trait)) %>%
    mutate(parameters = case_when(parameters=="receiver.trait" ~ "Tolerance of receiver from interspecific interaction",
                                  parameters=="emitter.trait"~"Effect of emitter on interspecific interaction",
                                  parameters=="abs(scaled.trait.dist)"~"Dissimilarity between interactors on interspecific interaction"),
           parameter.impacted ="INTER interactions")
  
  Intra.trait.df.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    dplyr::filter(!parameters %in% c("(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    mutate(trait=addline_format(trait)) %>%
    #dplyr::filter(!signif  =="") %>%
    mutate(parameters = case_when(parameters=="receiver.trait" ~ "Effect on intraspecific interactions"),
           parameter.impacted ="INTRA interactions")
  
  Lambda.trait.df.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    dplyr::filter(!parameters %in% c("(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    mutate(trait=addline_format(trait)) %>%
    #dplyr::filter(!signif  =="") %>%
    mutate(parameters = case_when(parameters=="receiver.trait" ~ "Effect on intrinsic fecundity"),
           estimate=estimate/10000,
           std.error=std.error/10000,
           parameter.impacted ="Intrinsic fitness")
  
  
  Rect.glm.theory.trait.plotlist[[country]] <- bind_rows(Inter.trait.df.i,
                                                         Intra.trait.df.i,
                                                         Lambda.trait.df.i) %>%
    mutate(parameters =factor(parameters,
                              levels=c("Effect on intrinsic fecundity",
                                       "Effect on intraspecific interactions",
                                       "Dissimilarity between interactors on interspecific interaction",
                                       "Tolerance of receiver from interspecific interaction",
                                       "Effect of emitter on interspecific interaction"))) %>%
    ggplot(aes(y=trait,
               x=estimate,
               color=parameters)) +
    geom_pointrange(aes(xmin=estimate - std.error,
                        xmax=estimate + std.error),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.3)) +
    scale_color_colorblind()+
    labs(x="coefficient",y="",
         color="") +
    guides(color = guide_legend(title.position = "top",
                                nrow=5))+ 
    geom_vline(xintercept=0,color="black") +
    annotate("text", x = -0.05, y=0.4, 
             label = "Increase competition/\nDecrease paramater",
             size=6) + 
    geom_segment(data=NULL,aes(x = 0, y = 0, xend = 0.1, yend = 0),
                 arrow = arrow(length = unit(0.4, "cm")), color="black")+
    annotate("text", x = 0.05, y=0.4, 
             label = "Increase facilitation/\nIncrease parameter",
             size=6) + 
    geom_segment(data=NULL,aes(x = 0, y = 0, xend = -0.1, yend = 0),
                 arrow = arrow(length = unit(0.4, "cm")),color="black") +
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
  
  Rect.glm.theory.trait.plotlist[[country]]
}
Rect.glm.theory.trait.plotlist[[paste0("spain")]]
Rect.glm.theory.trait.plotlist[[paste0("aus")]]
ggarrange(Rect.glm.theory.trait.plotlist[[paste0("spain")]],
          Rect.glm.theory.trait.plotlist[[paste0("aus")]],
          nrow=1,legend="bottom",label.x = 0.2,
          align=c("h"),
          common.legend = T, 
          labels=c("a.Spain","b. Australia"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
#figures/GLM.traits.pdf
#---- 1.6. Make circle of intercept ----
Circ.glm.theory.trait.plotlist<- list()
density.quantile.name <- c("intercept","low","medium","high")
country="aus"
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
  if(country=="aus"){
    Inter.trait.df <- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  if(country=="spain"){
    Inter.trait.df<- Cool.theory.trait.df[[country]]$Inter.trait.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  Inter.trait.df.i  <- Inter.trait.df %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    dplyr::filter(!parameters %in% c("(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    #dplyr::filter(!signif  =="") %>%
    dplyr::filter( best.model =="glm.inter.trait.i") %>%
    mutate(trait=addline_format(trait)) %>%
    mutate(parameters = case_when(parameters=="receiver.trait" ~ "Tolerance of receiver from interspecific interaction",
                                  parameters=="emitter.trait"~"Effect of emitter on interspecific interaction",
                                  parameters=="abs(scaled.trait.dist)"~"Dissimilarity between interactors on interspecific interaction"),
           parameter.impacted ="INTER interactions")
  
  Intra.trait.df.i  <- Cool.theory.trait.df[[country]]$Intra.trait.df %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    dplyr::filter(!parameters %in% c("(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    mutate(trait=addline_format(trait)) %>%
    mutate(parameters = case_when(parameters=="receiver.trait" ~ "Effect on intraspecific interactions"),
           parameter.impacted ="INTRA interactions")
  
  Lambda.trait.df.i  <- Cool.theory.trait.df[[country]]$Lambda.trait.df %>%
    dplyr::filter(density.quantile %in% c("low")) %>%
    dplyr::filter(!parameters %in% c("(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    mutate(trait=addline_format(trait)) %>%
    mutate(parameters = case_when(parameters=="receiver.trait" ~ "Effect on intrinsic fecundity"),
           estimate=estimate/10000,
           std.error=std.error/10000,
           parameter.impacted ="Intrinsic fitness")
  
  
  Circ.glm.theory.trait.plotlist[[country]]<- bind_rows(Inter.trait.df.i,
                                                        Intra.trait.df.i,
                                                        Lambda.trait.df.i) %>%
    mutate(parameters =factor(parameters,
                              levels=c("Effect on intrinsic fecundity",
                                       "Effect on intraspecific interactions",
                                       "Dissimilarity between interactors on interspecific interaction",
                                       "Tolerance of receiver from interspecific interaction",
                                       "Effect of emitter on interspecific interaction"))) %>%
    ggplot(aes(x=trait,
               y=estimate,
               color=parameters)) +
    geom_pointrange(aes(ymin=estimate - std.error,
                        ymax=estimate + std.error),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.2)) +
    guides(color = guide_legend(title.position = "top",
                                nrow=3))+ 
    annotate("text",x =2.5, y = -.02, label = "increase in \ncompetition",
             angle =0,color = "gray12",size = 4)+ # spain: 14, aus: 6
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 0.05),
                 arrow = arrow(length = unit(0.5, "cm")),color = "gray12")+
    annotate("text",x =2.5, y = 0.02, label = "increase in \nfacilitation",#"trait value",
             angle = 0,color = "gray12",size = 4)+
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = -0.05),
                 arrow = arrow(length = unit(0.5, "cm")),color = "gray12")+
    scale_color_colorblind()+
    geom_hline(yintercept=0,color="black") +
    scale_y_continuous(breaks=c(0),expand = c(0, 0)) + 
    labs(color="")+
    coord_polar() +
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title =element_text(size=20),
          legend.text =element_text(size=20),
          axis.text.x=element_text(color = "gray12", size = 16,
                                   face=ifelse(Inter.trait.df$trait %in% Inter.trait.df$trait[Inter.trait.df$p.value<0.10],
                                               "bold","plain")),
          plot.margin = unit(c(0,0,0,0),"cm"))
  Circ.glm.theory.trait.plotlist[[country]]
}
Circ.glm.theory.trait.plotlist[[paste0("spain")]]
Circ.glm.theory.trait.plotlist[[paste0("aus")]]
ggarrange(Circ.glm.theory.trait.plotlist[[paste0("spain")]],
          Circ.glm.theory.trait.plotlist[[paste0("aus")]],
          nrow=1,legend="bottom",label.x = 0.2,
          align=c("h"),
          common.legend = T, 
          labels=c("a.Spain","b. Australia"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
#figures/GLM.circle.traits.pdf


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
  glm.intra.trait.year.summary <- NULL
  #trait.i = "SLA"
  for( trait.i in names(trait.df)){
      trait.intra.year.df.i <-  trait.value.intra.year.df  %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::select(year,neigh,trait,sigmoid,
                      emitter.trait,receiver.trait) %>%
        mutate(year=as.numeric(year)) %>%
        left_join(env_pdsi %>% dplyr::select(year,Precip.extrem.scaled))
      
      glm.intra.trait.year.i <- glm(sigmoid ~ 1 + abs(Precip.extrem.scaled)+ receiver.trait:abs(Precip.extrem.scaled),
                               trait.intra.year.df.i ,
                               family="gaussian")
      #summary(glm.intra.trait.year.i)
      
      trait.intra.year.df.i<- trait.intra.year.df.i[complete.cases(trait.intra.year.df.i),]
      glm.intra.trait.year.zero.i <- glm(sigmoid ~ 1,
                                    trait.intra.year.df.i,
                                    na.action = na.omit,
                                    family="gaussian")
      
      
      glm.intra.trait.year.summary.i <-  as.data.frame(lrtest( glm.intra.trait.year.zero.i,
                                                          glm.intra.trait.year.i)) %>%
        mutate(model.comp =c("null","null vs linear"),
               model=c("glm.inter.trait.year.zero.i","glm.inter.trait.year.line.i")) %>%
        dplyr::filter(!model.comp=="H0") %>%
        rename("pvalue"="Pr(>Chisq)") %>%
        mutate(signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~" "),
               trait=trait.i)  
      glm.intra.trait.year.summary <- bind_rows( glm.intra.trait.year.summary, glm.intra.trait.year.summary.i)
      
      if(glm.intra.trait.year.summary.i$signif[which(glm.intra.trait.year.summary.i$model.comp=="null vs linear")]==" "){
        model.to.keep <- "glm.intra.trait.year.zero.i" 
      }else{
        model.to.keep <-  "glm.intra.trait.year.i" 
      }
      
      Intra.trait.year.df.i <- as.data.frame(summary(get(model.to.keep))$coef) %>%
        rownames_to_column("parameters") %>%
        rename("estimate"="Estimate",
               "std.error"="Std. Error",
               "t.value"="t value",
               "pvalue"="Pr(>|t|)") %>%
        mutate(trait=trait.i,
               best.model =model.to.keep,
               signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~""),
               RMSE = sqrt(mean(get(model.to.keep)$residuals^2)))
      
      Intra.trait.year.df <- bind_rows(Intra.trait.year.df,Intra.trait.year.df.i)
      
  }
  # Make data frame with trait and inter specific interactions
  trait.dist.year.df <- Realised.Obs.Year.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% # only INTRA
    left_join(specific.trait.dist,
              relationship ="many-to-many")
  Inter.trait.year.df <- NULL
  glm.inter.trait.year.summary<- NULL
  for( trait.i in names(trait.df)){
      trait.dist.year.df.i <-  trait.dist.year.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::select(neigh,year,trait,sigmoid,emitter.trait,
                      receiver.trait,scaled.trait.dist)%>%
        mutate(year=as.numeric(year)) %>%
        left_join(env_pdsi %>% dplyr::select(year,Precip.extrem.scaled))
      
       #rsq(trait.dist.year.df.i$theoretical.effect,trait.dist.year.df.i$trait.dist)
      glm.inter.trait.year.i <- glm(sigmoid ~ abs(Precip.extrem.scaled) + emitter.trait:abs(Precip.extrem.scaled) + receiver.trait:abs(Precip.extrem.scaled) + abs(scaled.trait.dist):abs(Precip.extrem.scaled), 
                               trait.dist.year.df.i,
                               family="gaussian")
      summary(glm.inter.trait.year.i)
      Inter.trait.year.df.i <- as.data.frame(summary(glm.inter.trait.year.i)$coef) %>%
        rownames_to_column("parameters") %>%
        rename("estimate"="Estimate",
               "std.error"="Std. Error",
               "t.value"="t value",
               "pvalue"="Pr(>|t|)") %>%
        mutate(trait=trait.i,
               signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~""))
      
      variable.to.keep <- Inter.trait.year.df.i$parameters[which(!Inter.trait.year.df.i$signif=="")]
      variable.to.keep <- variable.to.keep[!variable.to.keep == "(Intercept)"]
      if(!length(variable.to.keep) == 0){
        trait.dist.year.df.i <- trait.dist.year.df.i[complete.cases(trait.dist.year.df.i),]
        glm.inter.trait.zero.year.i <- glm(sigmoid ~ 1, 
                                      trait.dist.year.df.i,
                                      na.action = na.omit,
                                      family="gaussian")
        
        glm.inter.trait.line.year.i <- glm(formula(paste("sigmoid~ ",paste0(variable.to.keep, collapse="+"))), 
                                      trait.dist.year.df.i,
                                      na.action = na.omit,
                                      family="gaussian")
        
        glm.inter.trait.year.summary.i <- as.data.frame(lrtest( glm.inter.trait.zero.year.i,glm.inter.trait.line.year.i)) %>%
            mutate(model.comp =c("null","null vs linear"),
                   model=c("glm.inter.trait.zero.year.i","glm.inter.trait.line.year.i")) %>%
          dplyr::filter(!model.comp=="H0") %>%
          rename("pvalue"="Pr(>Chisq)") %>%
          mutate(signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                  (pvalue<0.05 & pvalue >0.01 )~"**",
                                  (pvalue<0.01)~"***",
                                  T~" "),
                 trait=trait.i)  
        glm.inter.trait.year.summary <- bind_rows( glm.inter.trait.year.summary, glm.inter.trait.year.summary.i)
        
        if(glm.inter.trait.year.summary.i$signif[which(glm.inter.trait.year.summary.i$model.comp=="null vs linear")]==" "){
          model.to.keep <- "glm.inter.trait.zero.year.i"
        }else{model.to.keep <- "glm.inter.trait.year.i"
        }
        best.model = model.to.keep }else{
          model.to.keep <- "glm.inter.trait.year.i"
          best.model <- "no significant variables"
        }
      Inter.trait.year.df.i <- as.data.frame(summary(get(model.to.keep))$coef) %>%
        rownames_to_column("parameters") %>%
        rename("estimate"="Estimate",
               "std.error"="Std. Error",
               "t.value"="t value",
               "pvalue"="Pr(>|t|)") %>%
        mutate(trait=trait.i,
               best.model =best.model,
               signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                                (pvalue<0.05 & pvalue >0.01 )~"**",
                                (pvalue<0.01)~"***",
                                T~""),
               RMSE = sqrt(mean(get(model.to.keep)$residuals^2)))
      
      Inter.trait.year.df <- bind_rows(Inter.trait.year.df,Inter.trait.year.df.i)
  }
  
  Cool.theory.trait.year.df[[country]] <- list(
    trait.dist.year.df=trait.dist.year.df,
    Inter.trait.year.df=Inter.trait.year.df,
    Intra.trait.year.df =Intra.trait.year.df,
    glm.inter.trait.year.summary=glm.inter.trait.year.summary)
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
  }
  if(country=="spain"){
    Inter.trait.year.df<- Cool.theory.trait.year.df[[country]]$Inter.trait.year.df %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) 
  }
  Inter.trait.year.df.i  <- Inter.trait.year.df %>%
    dplyr::filter(!parameters %in% c("abs(Precip.extrem.scaled)","(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    dplyr::filter(!signif  =="") %>%
    dplyr::filter( best.model =="glm.inter.trait.year.i") %>%
    mutate(trait=addline_format(trait)) %>%
    mutate(parameters = case_when(parameters=="abs(Precip.extrem.scaled)" ~ "Extreme Precipitation regime (either very dry or very warm) on INTERspecific interactions",
                                  parameters=="abs(Precip.extrem.scaled):emitter.trait"~"Effect of emitter on interspecific interaction according to precipitation regime",
                                  parameters=="abs(Precip.extrem.scaled):receiver.trait" ~ "Tolerance of receiver from interspecific interaction according to precipitation regime",
                                  parameters=="abs(Precip.extrem.scaled):abs(scaled.trait.dist)"~"Dissimilarity between interactors on interspecific interaction according to precipitation regime"),
           parameter.impacted ="INTER interactions")
  
  Intra.trait.year.df.i  <- Cool.theory.trait.year.df[[country]]$Intra.trait.year.df %>%
    dplyr::filter(!parameters %in% c("abs(Precip.extrem.scaled)","(Intercept)","I(receiver.trait^2)","I(emitter.trait^2)")) %>%
    mutate(trait=addline_format(trait)) %>%
    dplyr::filter(!signif  =="") %>%
    mutate(parameters = case_when(parameters=="abs(Precip.extrem.scaled)" ~ "Extreme Precipitation regime (either very dry or very warm) on INTRAspecific interactions",
                                  parameters=="abs(Precip.extrem.scaled):receiver.trait" ~ "Effect on intraspecific interactions according to precipitation regime"),
           parameter.impacted ="INTRA interactions")
  
  Rect.glm.theory.trait.year.plotlist[[country]] <- bind_rows(Inter.trait.year.df.i,
                                                         Intra.trait.year.df.i) %>%
    ggplot(aes(y=trait,
               x=estimate,
               color=parameters)) +
    geom_pointrange(aes(xmin=estimate - std.error,
                        xmax=estimate + std.error),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.3)) +
    scale_color_colorblind()+
    labs(x="coefficient",y="",
         color="") +
    guides(color = guide_legend(title.position = "top",
                                nrow=5))+ 
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
  
  Rect.glm.theory.trait.year.plotlist[[country]]
}
Rect.glm.theory.trait.year.plotlist[[paste0("spain")]]
Rect.glm.theory.trait.year.plotlist[[paste0("aus")]]
ggarrange(Rect.glm.theory.trait.year.plotlist[[paste0("spain")]],
          Rect.glm.theory.trait.year.plotlist[[paste0("aus")]],
          nrow=1,legend="bottom",label.x = 0.2,
          align=c("h"),
          common.legend = T, 
          labels=c("a.Spain","b. Australia"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
#figures/GLM.precip.trait.pdf
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.Looking at INTER Interactions----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.1. Make data df ----
country = "aus" 
Cool.trait.df <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  Focal.group.df <- Realised.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.50 ~ "Receives Fac",
                                 sigmoid > 0.50 ~ "Receives Comp",
                                 T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid")
  
  Neigh.group.df <- Realised.Int.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh, function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.50 ~ "Gives Fac",
                                 sigmoid > 0.50 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")
  
  
  
  if(country=="aus"){
    Neigh.group.df <-Neigh.group.df %>%
      mutate(neigh.org = case_when( sigmoid.neigh < 0.50 ~ "Gives Fac",
                                    sigmoid.neigh > 0.50 ~ "Gives Comp",
                                    T ~ "Gives both"))
    
    Focal.group.df <-   Focal.group.df %>%
      mutate(focal.org = case_when( sigmoid.focal< 0.5 ~ "Receives Fac",
                                    sigmoid.focal > 0.50 ~ "Receives Comp",
                                    T ~ "Receives both"))
  }
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
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
             focal= rep(rownames(trait.df),times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root mass density","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Stem height","Canopy shape") ~ "2.Aboveground",
                                         trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Water use efficiency","Canopy shape","Stem height","Leaf C to N ratio","Leaf area",
                                                      "Leaf nitrogen cc") ~ "2.Aboveground",
                                         trait %in% c("Mean fecundity") ~ "3.Reproduction"))
  }
  trait.dist.df <- Realised.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    group_by(neigh,focal) %>%
    summarise(Sigm.Q5= median(sigmoid),
              Sigm.neg= length(sigmoid[sigmoid<0]),
              Sigm.total = length(sigmoid)#,# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(Sigm.ratio= Sigm.neg/Sigm.total) %>% 
    left_join(specific.trait.dist %>% 
                left_join(Focal.group.df%>% dplyr::select(focal, focal.org,sigmoid.focal), multiple="all") %>%
                left_join(Neigh.group.df %>% dplyr::select(neigh, neigh.org,sigmoid.neigh), multiple="all"),
              relationship ="many-to-many") %>%
    mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                 Sigm.ratio< 0.5 ~"Fac",
                                 Sigm.ratio== 0.5 ~"Neutre"))
  
  sum.trait.dist.median.df <-  trait.dist.df %>%
    group_by(compORfac,trait) %>%
    mutate(Sigm.ratio = abs(Sigm.ratio -0.5)) %>%
    summarise(weighted.median = matrixStats::weightedMedian(scaled.trait.dist,Sigm.ratio,na.rm=T),
              weighted.mad = matrixStats::weightedMad(scaled.trait.dist,Sigm.ratio,na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T)) %>%
    mutate(weightedQ1= weighted.median-weighted.mad,
           weightedQ9= weighted.median+weighted.mad)
  
  
  Cool.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    specific.trait.dist=specific.trait.dist,
    sum.trait.dist.median.df=sum.trait.dist.median.df,
    Focal.group.df=Focal.group.df,
    Neigh.group.df=Neigh.group.df,
    Sym.group.df=Sym.group.df)
}

#---- 3.2. Detailed plot of dist and median interaction ----
Cool.detailed.trait.plotlist <- list()
for( country in country.list){
  
  Cool.detailed.trait.plotlist[[paste0(country)]] <- ggplot() +
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
               aes(y=scaled.trait.dist,
                   x=Sigm.Q5,
                   group=as.factor(compORfac),
                   fill=Sigm.ratio,
                   color=as.factor(compORfac)),
               shape=21,
               size=1,
               alpha=0.2) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df,
                    aes(y=weighted.median,
                        x=median.sigmoid,
                        xmin=Q10.sigmoid,
                        xmax=Q90.sigmoid,
                        color=as.factor(compORfac)),
                    shape=15,
                    size=1) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df,
                    aes(y=weighted.median,
                        x=median.sigmoid,
                        ymin=weightedQ1,
                        ymax=weightedQ9,
                        color=as.factor(compORfac)),
                    shape=15,
                    size=1) +
    facet_wrap(.~trait,scale="free") + 
    labs(y="Focal trait value - Neigh trait value",
         x="Median sigmoid value given by Neigh")+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual(values = c("#F21A00","#3B9AB2","#EBCC2A"))+
    theme_bw() 
  
}
Cool.detailed.trait.plotlist[[paste0("spain")]]
Cool.detailed.trait.plotlist[[paste0("aus")]]

#---- 3.3. Stat Test ----
Test.trait.list <- list()
for( country in country.list){
  trait.dist.df <- Cool.trait.df[[country]]$trait.dist.df
  sum.trait.dist.median.df <- Cool.trait.df[[country]]$sum.trait.dist.median.df
  Test.neigh.trait.df <- NULL
  Test.focal.trait.df <- NULL
  Test.trait.df <- NULL
  for( trait.i in levels(as.factor(trait.dist.df$trait))){
    trait.dist.df.i <-  trait.dist.df %>%
      dplyr::filter(trait==trait.i) %>%
      mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                   Sigm.ratio< 0.5 ~"Fac",
                                   Sigm.ratio== 0.5 ~"Neutre")) %>%
      dplyr::filter(!compORfac =="Neutre") %>%
      dplyr::select(trait,compORfac,trait.dist,Sigm.ratio)
    
    median.comp.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Comp" &
                                                                sum.trait.dist.median.df$trait==trait.i]
    median.fac.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Fac" &
                                                               sum.trait.dist.median.df$trait==trait.i]
    
    if(median.comp.n > median.fac.n){test.n ="greater"}else{test.n="less"}
    #test.n ="two.sided"
    test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$compORfac,
                        paired = F,
                        alternative = test.n)
    Test.trait.df.n <- data.frame(trait=trait.i,
                                  p.value = test$p.value,
                                  stat= test$statistic,
                                  alternative=test.n
    )
    Test.trait.df <- bind_rows(Test.trait.df,Test.trait.df.n)
    for( i in levels(as.factor(trait.dist.df$focal.org))){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(!focal.org==i )  %>%
        dplyr::select(trait,focal.org,trait.dist)
      if(nlevels(as.factor(trait.dist.df.i$focal.org))<2)next
      #https://www.r-bloggers.com/2020/06/wilcoxon-test-in-r-how-to-compare-2-groups-under-the-non-normality-assumption-2/
      test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$focal.org,
                          paired = F) # I am not sure if they should be paired or not 
      Test.focal.trait.df.n <- data.frame(trait=trait.i,
                                          focal.org.1 = levels(as.factor(trait.dist.df.i$focal.org))[1],
                                          focal.org.2= levels(as.factor(trait.dist.df.i$focal.org))[2],
                                          p.value = test$p.value,
                                          stat= test$statistic
      )
      Test.focal.trait.df <- bind_rows(Test.focal.trait.df,Test.focal.trait.df.n)
    }
    
    for( i in levels(as.factor(trait.dist.df$neigh.org))){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(!neigh.org==i )  %>%
        dplyr::select(trait,neigh.org,trait.dist)
      if(nlevels(as.factor(trait.dist.df.i$neigh.org))<2)next
      test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$neigh.org,
                          paired = F) # I am not sure if they should be paired or not 
      Test.neigh.trait.df.n <- data.frame(trait=trait.i,
                                          neigh.org.1 = levels(as.factor(trait.dist.df.i$neigh.org))[1],
                                          neigh.org.2= levels(as.factor(trait.dist.df.i$neigh.org))[2],
                                          p.value = test$p.value,
                                          stat= test$statistic)
      Test.neigh.trait.df <- bind_rows(Test.neigh.trait.df,Test.neigh.trait.df.n)
    }
  }
  Test.trait.list[[country]] <- list( Test.focal.trait.df= Test.focal.trait.df,
                                      Test.neigh.trait.df=Test.neigh.trait.df,
                                      Test.trait.df=Test.trait.df)
}
Test.trait.list[["spain"]]$Test.trait.df
Test.trait.list[["aus"]]$Test.trait.df
#---- 3.4. Summary plot with dist - rect ----
Cool.rect.trait.plotlist <- list()
country="aus"
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
for( country in country.list){
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  
  Test.trait.df <-   Test.trait.list[[country]]$Test.trait.df
  Cool.rect.trait.plotlist[[paste0(country,"_Pairwise_interactions")]] <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                              Sigm.ratio< 0.5 ~"Fac",
                                              Sigm.ratio== 0.5 ~"Neutre")) %>%
                 dplyr::filter(!compORfac =="Neutre") %>%
                 mutate(trait = factor(trait, levels = levels(specific.trait.dist$trait))),
               aes(x=scaled.trait.dist,  
                   y=trait,
                   #alpha=Sigm.ratio,
                   #fill=compORfac, # focal.org,
                   color=compORfac), # focal.org,
               position = position_dodge(width = 0.7),
               shape=16,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df %>%
                      dplyr::filter(!compORfac =="Neutre") %>%
                      mutate(trait = factor(trait, 
                                            levels = levels(specific.trait.dist$trait))),
                    position = position_dodge(width = 0.7),
                    aes(x=weighted.median, # median.trait.dist,
                        y=trait,
                        xmin=weightedQ1,#Q10.trait.dist,
                        xmax=weightedQ9,#Q90.trait.dist,
                        color=as.factor(compORfac)),
                    shape=16,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1)) +
    annotate("text", x = 2, y=0.2, 
             label = "Focal has higher \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    annotate("text", x = -2, y=0.2, 
             label = "Focal has lower \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Pairwise interaction mainly")+#"Species grouping based on received interactions")+
    scale_color_manual(values = c("#F21A00","#3B9AB2"),
                       labels=c("Competitive","Facilitative"))+ # "#EBCC2A",
    scale_alpha_continuous(range=c(0.1,1)) +
    theme_bw()  +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          legend.title =element_text(size=16),
          legend.text =element_text(size=14),
          axis.title=element_text(size=14),
          axis.text.x=element_text(color = "gray12", size = 16,
                                   face=ifelse(Test.trait.df$trait %in% Test.trait.df$trait[Test.trait.df$p.value<0.10],
                                               "bold","plain")))
  
}
Cool.rect.trait.plotlist[["aus_Pairwise_interactions"]]
Cool.rect.trait.plotlist[["spain_Pairwise_interactions"]]
#---- 3.5. Summary plot with dist - circ ----
Cool.circ.trait.plotlist <- list()
country="aus"
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  
  Test.trait.df <-   Test.trait.list[[country]]$Test.trait.df %>%
    mutate(labels.pvalue = case_when((p.value <= 0.1 & p.value > 0.05 )~ "*",
                                     (p.value <= 0.05 & p.value > 0.01 )~ "**",
                                     p.value <= 0.01 ~ "***",
                                     T~ ""))
  
  Cool.circ.trait.plotlist[[paste0(country,"_Pairwise_interactions")]] <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                              Sigm.ratio< 0.5 ~"Fac",
                                              Sigm.ratio== 0.5 ~"Neutre")) %>%
                 dplyr::filter(!compORfac =="Neutre") %>%
                 mutate(trait = factor(trait, levels = levels(specific.trait.dist$trait))),
               #mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(y=scaled.trait.dist,  
                   x=trait,
                   #alpha=Sigm.ratio,
                   #fill=compORfac, # focal.org,
                   color=compORfac), # focal.org,
               position = position_dodge(width = 0.7),
               shape=16,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df %>%
                      dplyr::filter(!compORfac =="Neutre"), 
                    position = position_dodge(width = 0.7),
                    aes(y=weighted.median, # median.trait.dist,
                        x=trait,
                        ymin=weightedQ1,#Q10.trait.dist,
                        ymax=weightedQ9,#Q90.trait.dist,
                        color=as.factor(compORfac)),
                    shape=16,
                    size=1) +
    #geom_text(data=Test.trait.df,aes(x=trait,y=3,label=labels.pvalue),size=12)+
    geom_hline(yintercept=0,color="black") +
    scale_y_continuous(breaks=c(0),expand = c(0, 0)) + 
    scale_x_discrete(labels=addline_format(paste(Test.trait.df$trait,Test.trait.df$labels.pvalue))) +
    labs(color="Pairwise interaction mainly")+
    scale_color_manual(values = c("#F21A00","#3B9AB2"),
                       labels=c("Competitive","Facilitative"))+ # "#EBCC2A",
    coord_polar() +
    annotate("text",x =3.4, y = 1.3, label = "Focal>Neigh",
             angle =-8,color = "gray12",size = 4.5)+ # spain: 14, aus: 6
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 3),
                 arrow = arrow(length = unit(0.5, "cm")))+
    annotate("text",x =3.64, y = -1.3, label = "Focal<Neigh",#"trait value",
             angle = -8,color = "gray12",size = 4.5)+
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = -3),
                 arrow = arrow(length = unit(0.5, "cm")))+
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           fill = guide_legend(title.position = "top",
                               nrow=2)) + 
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title =element_text(size=20),
          legend.text =element_text(size=20),
          axis.text.x=element_text(color = "gray12", size = 16,
                                   face=ifelse(Test.trait.df$trait %in% Test.trait.df$trait[Test.trait.df$p.value<0.10],
                                               "bold","plain")),
          plot.margin = unit(c(0,0,0,0),"cm"))
  Cool.circ.trait.plotlist[[paste0(country,"_Pairwise_interactions")]]
  
}
Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]] #figures/Cool.circ.pairwise.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]] # figures/Cool.circ.pairwise.trait.aus.pdf

Cool.circ.trait.plotlist[["Pairwise_interactions"]] <-
  ggarrange(Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]],
            Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]],
            nrow=1,legend="bottom",
            common.legend = T, labels=c("a. Spain","b. Australia"),
            font.label = list(size = 20, color = "black", face = "bold", family = NULL))

Cool.circ.trait.plotlist["Pairwise_interactions"] # figures/Cool.circ.pairwise.trait.pdf

Cool.circ.trait.plotlist[[paste0("spain_circ")]] # figures/Cool.circ.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_circ")]] # figures/Cool.circ.trait.aus.pdf



#---- 3.6. PCA----
# PCA analysis
#https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/

country = "aus"
pca.plotlist <- list()
for( country in country.list){
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  Focal.group.df<- Cool.trait.df[[country]]$Focal.group.df
  Neigh.group.df<- Cool.trait.df[[country]]$Neigh.group.df
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] 
  library(Hmisc)
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  res2<-rcorr(as.matrix(trait.df))
  #view(flattenCorrMatrix(res2$r, res2$P))
  # for spain, SLA is significantly correclated with mean fecundity,Root mass density, SRA, water supp
  # & SRL with Canopy shape
  write.csv(flattenCorrMatrix(res2$r, res2$P),
            file=paste0("results/correlation.matrix.",country,".csv"))
  names(trait.df)
  if(country=="spain"){
    to.keep.trait <- c("Root mass density","Leaf C to N ratio",#"SRA",
                       "SLA","SRL","Canopy shape","Leaf area index"
    )
  }else{
    to.keep.trait<- c("Root mass density",#"Root length",
                      "SRL","C13 water use efficiency")
    
  }
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]  %>%
    dplyr::select(all_of(to.keep.trait))
  
  pca.trait.df <- MFA(trait.df, #%>%
                      #mutate(received =Focal.group.df$sigmoid.focal,
                      #given = Neigh.group.df$sigmoid.neigh),
                      group = c(rep(1, times= ncol(trait.df)+0)), 
                      type = c(rep("s", times= ncol(trait.df)+0)), # s = quantitative , n= factorial
                      name.group = c(colnames(trait.df)),#,"received","given"),
                      graph = T) 
  
  eig.val <- get_eigenvalue(pca.trait.df)
  head(eig.val)
  fviz_contrib(pca.trait.df, choice = "quanti.var", axes = 1, top = 20,palette = kelly(n=22))
  fviz_contrib(pca.trait.df, choice = "quanti.var", axes = 2, top = 20,palette = kelly(n=22))
  fviz_mfa_var(pca.trait.df, "quanti.var", palette = kelly(n=22), 
               col.var.sup = "violet", repel = TRUE)
  
  fviz_mfa_ind(pca.trait.df, col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE,visible = "quali.var")
  
  species.pca.df <- data.frame(xaxis=pca.trait.df$ind$coord[,1],
                               yaxis=pca.trait.df$ind$coord[,2]) %>%
    bind_cols(Focal.group.df) %>%
    bind_cols(Neigh.group.df) 
  trait.pca.df <- data.frame(xaxis=pca.trait.df$quanti.var$coord[,1],
                             yaxis=pca.trait.df$quanti.var$coord[,2]) %>%
    rownames_to_column("traitname")
  xperc <- eig.val[1,2]
  yperc <- eig.val[2,2]
  #if(country=="spain"){
  #species.pca.df <- data.frame(xaxis=-pca.trait.df$ind$coord[,2],
  #                             yaxis=pca.trait.df$ind$coord[,1]) %>%
  # bind_cols(Focal.group.df) %>%
  # bind_cols(Neigh.group.df) 
  #trait.pca.df <- data.frame(xaxis=-pca.trait.df$quanti.var$coord[,2],
  #                            yaxis=pca.trait.df$quanti.var$coord[,1]) %>%
  #  rownames_to_column("traitname")
  ## xperc <- eig.val[2,2]
  # yperc <- eig.val[1,2]
  #}
  
  dummy.df <- data.frame(xaxis=c(1.2,1.2),
                         yaxis=c(2.3,2),
                         focal=c("Received",
                                 "Given"))
  pca.plotlist[[country]] <- ggplot() +
    annotate(geom="point",x=dummy.df$xaxis,y=dummy.df$yaxis,
             shape=c(21,24),size=10,stroke=2) +
    annotate(geom="text",
             x=dummy.df$xaxis+0.5,y=dummy.df$yaxis,
             label=dummy.df$focal,size=6) +
    geom_point(data=species.pca.df,
               aes(x=xaxis-0.1,y=yaxis,fill=sigmoid.focal),
               shape=21,size=6) +
    geom_point(data=species.pca.df,
               aes(x=xaxis+0.1,y=yaxis,fill=sigmoid.neigh),
               shape=24,size=6) +
    geom_text(data=species.pca.df,
              aes(x=xaxis,y=yaxis+0.15,
                  label=addline_format(focal)),
              check_overlap = F,size=5)+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"),
                         limits=c(0,1),
                         breaks=c(0,1),
                         labels=c("100%\nFacilitatives",
                                  "100%\nCompetitives")) +
    geom_segment(data=trait.pca.df,
                 aes(x=0,y=0,xend=xaxis,yend=yaxis),
                 arrow=arrow(length = unit(0.5, "cm"), 
                             type = "closed"))   +
    geom_text(data=trait.pca.df,
              aes(x=xaxis+0.05,y=ifelse(yaxis>0, yaxis+0.05,yaxis-0.05),
                  label=traitname),
              check_overlap = F,size=5)+
    labs(x=paste0("PCA dimension ",round(xperc,digits=1) ,"%"),
         y=paste0("PCA dimension ",round(yperc,digits=1) ,"%"),
         fill="Ratio of Competitive to \nFacilitative interactions")+
    theme_bw()+
    theme(legend.key.size = unit(1, 'cm'),
          legend.position="bottom",
          legend.title.position = "top",
          panel.grid.minor = element_blank(),
          legend.title =element_text(size=20),
          legend.text =element_text(size=20),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          plot.margin = unit(c(1,0.5,0,0.5),"cm"))
  pca.plotlist[[country]]
}
pca.plotlist[["spain"]]
ggarrange(pca.plotlist[["spain"]],
          pca.plotlist[["aus"]],
          nrow=1,legend="bottom",
          common.legend = T, labels=c("a. Spain","b. Australia"),
          font.label = list(size = 20, color = "black", face = "bold", family = NULL)) # figures/PCA.interaction.pdf
