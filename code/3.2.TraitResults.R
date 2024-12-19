
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
#---- 0.1. Import results----
# parameter from models
load(paste0(home.dic,"results/Parameters_alpha.RData")) 
# Raw sigmoid
Param.sigm.df <- list()
Param.sigm.df.aus <- read.csv(paste0(project.dic,"results/Param.sigmoid.aus.csv.gz"))
Param.sigm.df.spain <- read.csv(paste0(project.dic,"results/Param.sigmoid.spain.csv.gz"))
Param.sigm.df$spain <- Param.sigm.df.spain
Param.sigm.df$aus<- Param.sigm.df.aus
save(Param.sigm.df,
     file=paste0(project.dic,"results/Param.sigm.df.RData"))
load(paste0(project.dic,"results/Param.sigm.df.RData"))
# realised interactions across time
load(paste0(project.dic,"results/Realised.Int.list.RData"))

# realised interactions for each year
load(paste0(project.dic,"results/Realised.Int.Year.list.RData")) 

load(paste0(project.dic,"results/Realised.Int.Obs.list.RData")) 
# PDSI
env_pdsi_aus <- read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) 
env_pdsi_spain <- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv")) 

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.Looking at Fac for answers---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.Traits Graph related- Looking at traits for answers---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1 With different grouping of functional group ----
Realised.Int.Group.plotlist <- list()
grouping.list <- list()
net.fct.group <- list()
country ="aus"
for( country in "aus"){
  trait.gr <- c("pol","above","below")
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  # i= 1
  for(i in 1:3){
    
    grouping.list[[i]] <- get(paste0("clean.data.",country))[[paste0(country,"_",trait.gr[i],"_grouping")]] %>%
      dplyr::select(final.code,group.cluster,name.cluster, ends_with(".cluster")) 
    names.clusters <- names(grouping.list[[i]])[!names(grouping.list[[i]]) %in% c("final.code","group.cluster","name.cluster")]
    
    # n = names.clusters
    for(n in c("name.cluster",names.clusters)){
      
      Realised.Int.list.df.with.group <-   Realised.Int.list[[country]]   %>%
        left_join(grouping.list[[i]] %>%
                    dplyr::select(final.code,n) %>%
                    rename("focal"="final.code",
                           "focal.name.cluster"=n)) %>%
        left_join(grouping.list[[i]] %>%
                    dplyr::select(final.code,n)%>%
                    rename("neigh"="final.code",
                           "neigh.name.cluster"=n)) %>%
        dplyr::filter(!is.na(focal.name.cluster)) %>%
        dplyr::filter(!is.na(neigh.name.cluster))
      
      Realised.Int.focal.neigh <- Realised.Int.list.df.with.group %>%
        group_by(neigh.name.cluster,focal.name.cluster)%>%
        summarise(RE.neg= length(realised.effect[realised.effect<1]),
                  RE.total = length(realised.effect),# high mean high competition
                  RE.mean = mean(realised.effect, na.rm = T),
                  RE.sd = sd(realised.effect,na.rm = T),
                  RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
                  RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
                  RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
        ungroup() %>%
        mutate(RE.ratio = RE.neg/RE.total) %>%
        as.data.frame()
      
      Realised.Int.focal <-Realised.Int.list.df.with.group %>%
        group_by(focal.name.cluster)%>%
        summarise(RE.neg= length(realised.effect[realised.effect<1]),
                  RE.total = length(realised.effect),# high mean high competition
                  RE.mean = mean(realised.effect, na.rm = T),
                  RE.sd = sd(realised.effect,na.rm = T),
                  RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
                  RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
                  RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
        ungroup() %>%
        mutate(RE.ratio = RE.neg/RE.total,
               neigh.name.cluster ="total") %>%
        as.data.frame() 
      
      Realised.Int.neigh <-Realised.Int.list.df.with.group %>%
        group_by(neigh.name.cluster)%>%
        summarise(RE.neg= length(realised.effect[realised.effect<1]),
                  RE.total = length(realised.effect),# high mean high competition
                  RE.mean = mean(realised.effect, na.rm = T),
                  RE.sd = sd(realised.effect,na.rm = T),
                  RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
                  RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
                  RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
        ungroup() %>%
        mutate(RE.ratio = RE.neg/RE.total,
               focal.name.cluster="total") %>%
        as.data.frame() 
      #heatmap
      Realised.Int.total <- Realised.Int.list.df.with.group %>%
        summarise(RE.neg= length(realised.effect[realised.effect<1]),
                  RE.total = length(realised.effect),# high mean high competition
                  RE.mean = mean(realised.effect, na.rm = T),
                  RE.sd = sd(realised.effect,na.rm = T),
                  RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
                  RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
                  RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
        ungroup() %>%
        mutate(RE.ratio = RE.neg/RE.total,
               focal.name.cluster ="total",
               neigh.name.cluster="total") %>%
        as.data.frame() 
      
      Realised.Int.Group.plotlist[[paste(country,"_",trait.gr[i],"_",n)]] <- Realised.Int.focal.neigh  %>%
        #bind_rows(Realised.Int.neigh,Realised.Int.focal, Realised.Int.total) %>%
        mutate(transparency = case_when(RE.mean < 0.05 ~ 0.2,
                                        (RE.mean >0.05 & RE.mean <=0.1) ~ 0.4,
                                        (RE.mean>0.1 &RE.mean <=0.25) ~ 0.6,
                                        (RE.mean >0.25 & RE.mean<=0.5) ~ 0.8,
                                        (RE.mean>0.5 ) ~ 1)) %>%
        mutate(neigh.name.cluster=factor(neigh.name.cluster, levels=c(levels(as.factor(grouping.list[[i]][,n])),
                                                                      "total")),
               focal.name.cluster= factor(focal.name.cluster,levels=c("total",
                                                                      rev(levels(as.factor(grouping.list[[i]][,n])))))) %>%
        ggplot(aes(x =neigh.name.cluster,y =focal.name.cluster)) + 
        geom_tile(aes(fill = RE.ratio,alpha=abs(RE.mean-1))) + 
        scale_x_discrete(position = "top") +
        scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                                   101, 
                                                   type = "continuous"),
                             limits=c(0,1))   +
        #annotate("rect",xmin = 0.5, xmax = 5.5, ymin = .5,
        #        ymax = 1.5, color="black",fill="transparent") + 
        #annotate("rect",ymin = 0.5, ymax = 5.5, xmin = 5- .5,
        #         xmax = 5 + .5, color="black",fill="transparent") + 
        theme_few() +
        scale_alpha_continuous(range = c(0.5, 1))+
        labs(alpha ="Mean effect",
             fill="Ratio of Comp/Fac",
             x="Giver",
             y="Receiver") +
        theme(legend.position ="right",
              axis.title = element_text(size=14),
              axis.text.x=element_text(angle=90,size=14),
              axis.text.y=element_text(size=14),
              axis.ticks =element_blank(),
              axis.line = element_blank())
      
      # Network
      
      ratio.mat <- Realised.Int.focal.neigh %>%
        dplyr::select(RE.ratio, neigh.name.cluster, focal.name.cluster) %>%
        spread(neigh.name.cluster,RE.ratio) %>%
        column_to_rownames("focal.name.cluster") %>%
        as.matrix()
      
      strength.mat <- Realised.Int.focal.neigh %>%
        dplyr::select(RE.Q5, neigh.name.cluster, focal.name.cluster) %>%
        mutate(RE.Q5 = abs(RE.Q5-1)) %>%
        spread(neigh.name.cluster,RE.Q5) %>%
        column_to_rownames("focal.name.cluster") %>%
        as.matrix()
      
      plot.network.gradient.int(ratio.mat,strength.mat,"",minimum.stength = 0.005)
      net.fct.group[[paste(country,"_",trait.gr[i],"_",n)]] <- recordPlot()
    }
  }
}
Realised.Int.Group.plotlist[[1]]#figures/Heatmap.trait.aus.poll.pdf
Realised.Int.Group.plotlist[[2]]#figures/Heatmap.trait.aus.above.pdf
Realised.Int.Group.plotlist[[3]] #figures/Heatmap.trait.aus.below.pdf
net.fct.group[[1]]#figures/Network.all.aus.poll.pdf
net.fct.group[[2]]#figures/Network.flower.size.aus.above.pdf
net.fct.group[[3]]#figures/Network.all.aus.above.pdf
net.fct.group[[4]]#figures/Network.SLA.aus.above.pdf
net.fct.group[[5]]#figures/Network.height.aus.above.pdf
net.fct.group[[6]]#figures/Network.area.aus.above.pdf
net.fct.group[[7]]#figures/Network.volume.aus.above.pdf
net.fct.group[[8]]#figures/Network.all.aus.below.pdf
net.fct.group[[9]]#figures/Network.srl.aus.below.pdf
net.fct.group[[10]]#figures/Network.height.root.aus.below.pdf

#---- 2.2 Along the fast to slow gradient ----
Realised.Int.FSspectrum.plotlist <- list()

country = "aus"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]  %>%
    rownames_to_column("species")
  if(country=="spain"){
    trait.df <- trait.df %>%
      bind_rows(data.frame(species="PAIN",
                           coord.slow.to.fast =-10))   
  }
  Realised.Int.sum <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  #hist(Realised.Int.sum$realised.effect,breaks=150)
  #Realised.Int.FSspectrum.plotlist[[country]] 
  #figures/Heatmap.N15_spain.pdf
  #figures/Heatmap.root.volume.less0.5_aus.pdf
  library(tibble)
  library(hrbrthemes)
  Realised.Int.FSspectrum.plotlist[[country]] <- Realised.Int.sum %>%
    #dplyr::filter(RE.Q5 < 1000) %>%
    mutate(neigh=factor(neigh, levels=c(trait.df$species[order(trait.df$coord.slow.to.fast)])),
           focal = factor(focal,levels=rev(trait.df$species[order(trait.df$coord.slow.to.fast)]))) %>%
    ggplot(aes(x =neigh,y =focal)) + 
    geom_tile(aes(fill = RE.ratio.sig,alpha=abs(RE.Q5.sig))) + 
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))   +
    theme_ipsum() +
    scale_alpha_continuous(range = c(0.4, 1))+
    labs(alpha ="Mean effect",
         fill="Ratio of Comp/Fac",
         x="Neighbour",
         y="Focal") +
    theme(legend.position ="bottom",
          axis.title = element_text(size=14),
          axis.text.x=element_text(angle=90,size=14),
          axis.text.y=element_text(size=14),
          axis.ticks =element_blank(),
          axis.line = element_blank())
  Realised.Int.FSspectrum.plotlist[[country]]
  
}
Realised.Int.FSspectrum.plotlist[["aus"]] #figures/Heatmap.FS.spectrum_aus.pdf

Realised.Int.FSspectrum.plotlist[["spain"]] #figures/Heatmap.FS.spectrum_spain.pdf

#---- 2.3 Along the dist between interactors----
Realised.Int.Dist.plotlist <- list()
country = "spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  
  Realised.Int.sum <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  head(trait.df)
  Dist.plot <- Realised.Int.sum %>%
    left_join(trait.df %>%
                dplyr::select(c("focal","neigh","dist"))) %>%
    #dplyr::filter(!dist==0) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.50 ~ "Comp",
                              RE.ratio.sig<=0.50 ~ "Fac",
                              T ~ "Neutral")) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(y =RE.Q5.sig,x =dist)) +
    geom_point(aes(fill=RE.ratio.sig, group=PosNeg,color=PosNeg),
               shape=21,
               alpha=0.5,
               size=4,
               stroke = 2) +
    geom_smooth(aes( group=PosNeg,color=PosNeg),
                method="glm",se = T,size=2) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual( values = c("#F21A00","#3B9AB2","#EBCC2A")) +
    theme_bw()
  library(ggExtra)
  Realised.Int.Dist.plotlist[[country]] <- ggMarginal(Dist.plot , type = "density",
                                                      groupColour = TRUE, groupFill = TRUE)
  Realised.Int.Dist.plotlist[[country]]
}
Realised.Int.Dist.plotlist[["aus"]] #figures/Dist.Q5_aus.pdf

Realised.Int.Dist.plotlist[["spain"]] #figures/Dist.Q5_spain.pdf


Realised.Int.Dist.Year.plotlist <- list()

country = "spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  
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
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  Realised.Int.Dist.Year.plotlist[[country]] <- Realised.Int.sum %>%
    left_join(trait.df %>%
                dplyr::select(c("focal","neigh","dist"))) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.50 ~ "Comp",
                              RE.ratio.sig<=0.50 ~ "Fac",
                              T ~ "Neutral")) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(y =RE.Q5.sig,x =dist)) +
    geom_point(aes(fill=RE.ratio.sig),position="jitter",
               shape=21,
               alpha=0.2,
               size=2,
               stroke = 2) +
    facet_wrap(.~ year) +
    geom_smooth(aes( group=as.factor(PosNeg),color=as.factor(PosNeg)),
                method="glm",se = T,size=2) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous")) +
    scale_color_manual(values=rev(c(wes_palette("Zissou1", 
                                                2, 
                                                type = "continuous")))) +
    #coord_cartesian(ylim=c(0.5,2),clip="on") +
    theme_bw()
  Realised.Int.Dist.Year.plotlist[[country]]
  
  # Grouped Scatter plot with marginal density plots
  #ggscatterhist(
  #test.df, x = "coord.slow.to.fast", y = "RE.ratio",
  #color = "PosNeg", size = 3, alpha = 0.6,
  #palette = rev(c(wes_palette("Zissou1", 
  #                            2, type = "continuous"))),
  #margin.params = list(fill = "PosNeg", color = "black", size = 0.2))
}
Realised.Int.Dist.plotlist[["aus"]] #figures/Ratio.Dist_aus.pdf
Realised.Int.Dist.Year.plotlist[["aus"]]
Realised.Int.Dist.plotlist[["spain"]] #figures/Ratio.Dist_spain.pdf
Realised.Int.Dist.Year.plotlist[["spain"]]

#---- 2.4 Along the trait dist between interactors----
Realised.Int.Dist.plotlist <- list()
country = "spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    if( i =="coord.slow.to.fast") next
    specific.trait.dist.n <-outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(Code.focal.list,each=length(Code.focal.list)),
             focal= rep(Code.focal.list,times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  trait.df <- trait.df %>%
    rownames_to_column("focal")
  Realised.Int.sum <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  
  Ratio.focal <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(focal) %>%
    summarise(RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  head(trait.df)
  test.df <-  Realised.Int.Obs.list[[country]] %>%
    full_join(specific.trait.dist,by=c("focal","neigh"),multiple = "all") 
  
  Dist.plot <-  Realised.Int.Obs.list[[country]] %>%
    full_join(specific.trait.dist,by=c("focal","neigh"),multiple = "all") %>%
    mutate(PosNeg = case_when(sigmoid<0 ~ "Comp",
                              sigmoid>0 ~ "Fac",
                              T ~ "Neutral"),
           grouping.fact = paste0(PosNeg,"_",trait)) %>%
    ggplot(aes(y =sigmoid,x =abs(trait.dist))) +
    geom_point(aes(fill=sigmoid),
               position="jitter",
               shape=21,
               alpha=0.3,
               size=2,
               stroke = 0) +
    geom_smooth(aes(group=trait, #grouping.fact,
                    linetype=trait),
                method="lm") +
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                   51, 
                                                   type = "continuous")),
                         limits=c(-0.35,0.35))+
    scale_color_manual(values = rev(wes_palette("Zissou1", 
                                                2, 
                                                type = "continuous")))+
    facet_wrap(.~focal,scale="free") +
    theme_bw() 
  Dist.plot
  # figures/Dist.trait.aus.pdf
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    #dplyr::filter(focal=="ARCA") %>%
    #dplyr::filter(!dist==0) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.51 ~ "Comp",
                              RE.ratio.sig<0.49 ~ "Fac",
                              T ~ "Neutral")) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(x =trait,y =trait.dist)) +
    geom_point(aes(fill=RE.ratio.sig),
               position="jitter",
               shape=21,
               alpha=0.3,
               size=2,
               stroke = 0) +
    stat_summary(aes(group=as.factor(PosNeg),color=as.factor(PosNeg)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual( values = c("#F21A00","#3B9AB2","#EBCC2A")) +
    facet_wrap(.~focal) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=66))
  Dist.plot
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    #mutate(focal = factor(focal,levels=rev(trait.df$focal[order(trait.df$coord.slow.to.fast)]))) %>%
    mutate(focal = factor(focal,levels=rev(Ratio.focal$focal[order(Ratio.focal$RE.ratio.sig)]))) %>%
    mutate(PosNeg = case_when(RE.Q5.sig<0 ~ "Comp",
                              RE.Q5.sig>0~ "Fac",
                              T ~ "Neutral")) %>%
    ggplot(aes(x =focal,y =trait.dist)) +
    geom_point(aes(fill=RE.Q5.sig),
               position="jitter",
               alpha=0.3,
               shape=21,
               stroke=0,
               size=2) +
    stat_summary(aes(group=as.factor(PosNeg),
                     color=as.factor(PosNeg)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    #scale_shape_manual(values=1:5) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual( values = c("#F21A00","#3B9AB2","#EBCC2A")) +
    theme_bw() +
    facet_wrap(.~trait)+
    theme(axis.text.x = element_text(angle=66))
  Dist.plot
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    #mutate(focal = factor(focal,levels=rev(trait.df$focal[order(trait.df$coord.slow.to.fast)]))) %>%
    mutate(focal = factor(focal,levels=rev(Ratio.focal$focal[order(Ratio.focal$RE.ratio.sig)]))) %>%
    mutate(PosNeg = case_when(RE.Q5.sig<0 ~ "Comp",
                              RE.Q5.sig>0 ~ "Fac",
                              T ~ "Neutral")) %>%
    dplyr::filter(PosNeg=="Fac") %>%
    ggplot(aes(x =trait,y =trait.dist)) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 geom = "line",size=0.5) +
    #scale_shape_manual(values=1:5) +
    scale_color_paletteer_d("MetBrewer::Benedictus") +
    scale_fill_paletteer_d("MetBrewer::Benedictus") +
    theme_bw() 
  Dist.plot
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    # mutate(focal = factor(focal,levels=rev(trait.df$focal[order(trait.df$coord.slow.to.fast)]))) %>%
    mutate(focal = factor(focal,levels=rev(Ratio.focal$focal[order(Ratio.focal$RE.ratio.sig)]))) %>%
    mutate(PosNeg =  case_when(RE.Q5.sig<0 ~ "Comp",
                               RE.Q5.sig>0 ~ "Fac",
                               T ~ "Neutral")) %>%
    dplyr::filter(PosNeg=="Comp") %>%
    ggplot(aes(x =trait,y =trait.dist)) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 geom = "line",size=0.5) +
    #scale_shape_manual(values=1:5) +
    scale_color_paletteer_d("MetBrewer::Benedictus") +
    scale_fill_paletteer_d("MetBrewer::Benedictus") +
    theme_bw() 
  Dist.plot
  
}

sem.trait.plotlist <- list()
country = "aus"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
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
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    if( i =="coord.slow.to.fast") next
    if(country=="spain"){Code.focal.list <- Code.focal.list[!Code.focal.list == "PAIN"] }
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(Code.focal.list,each=length(Code.focal.list)),
             focal= rep(Code.focal.list,times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  
  sem.df <- NULL
  for( sp in Code.focal.list){
    print(sp)
    test.df <-  Realised.Int.Obs.list[[country]] %>%
      full_join(specific.trait.dist,
                by=c("focal","neigh"), 
                multiple = "all") %>%
      # filter(focal ==sp) %>%
      filter(focal ==sp) %>%
      aggregate(sigmoid ~ density + focal + neigh + trait.dist + trait, median) %>%
      mutate(trait.dist =abs(trait.dist))%>%
      spread(trait,trait.dist) #%>%
    #left_join(env_pdsi %>%  mutate(year = as.numeric(year)),
    #        by=c("year"),multiple = "all")  
    if(country=="aus"){
      test.df <- test.df  %>%
        dplyr::select_if(~ sum(!is.na(.))>10) %>%
        drop_na()
      var.names <- names(test.df)[names(test.df) %in% c("sla","width.longest","srl","root.length",
                                                        "flower.size.numb","height","number.of.root.tips")]
      
    }else{  test.df <- test.df %>%select_if(~ sum(!is.na(.))>10)
    var.names <- names(test.df)[names(test.df) %in% c("root.volume.0.5","C13","CS","heigh","N15",
                                                      "SLA","TDMr","SRA")]
    }
    ggplot(test.df,aes(x=number.of.root.tips,y=sigmoid)) +
      geom_point() + geom_smooth(method="lm")
    head(test.df)
    sem.median  <- psem(
      lm(as.formula(paste("sigmoid ~ 1 +", paste(c(var.names,"density"), collapse= "+"))),
         test.df ),
      #MASS::glm.nb(density ~ aus_pdsi + pol.traits,
      #Realised.Int.Year.Sem.Ratio),
      #lm(density ~ 1 + Precip.extrem,
      #  test.df),
      test.df)
    
    #plot.psemhl(sem.median , correlation = T, layout = "tree") 
    #summary(sem.median)
    sem.df.n <- coefs(sem.median)[,1:8] %>%
      mutate(focal=sp)
    
    sem.df <- bind_rows(sem.df,sem.df.n)
  }
  
  Realised.Int.sum <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(focal) %>%
    summarise(RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  
  trait.df.dist <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  trait.df <- trait.df %>% rownames_to_column("focal")
  head(sem.df)
  sem.trait.plotlist[[country]] <- sem.df %>%
    dplyr::filter(Response=="sigmoid") %>%
    dplyr::filter(P.Value < 0.1) %>%
    #mutate(focal = factor(focal,levels=rev(unique(trait.df.dist$focal[order(trait.df.dist$dist)])))) %>%
    #mutate(focal = factor(focal,levels=rev(unique(trait.df$focal[order(trait.df$coord.slow.to.fast)])))) %>%
    mutate(focal = factor(focal,levels=rev(Realised.Int.sum$focal[order(Realised.Int.sum$RE.ratio.sig)]))) %>%
    #mutate(focal= factor(focal,levels=rev(Realised.Int.sum$neigh[order(Realised.Int.sum$RE.ratio.sig)]))) %>%
    ggplot(aes(x=as.factor(Predictor),
               y=Std.Estimate,
               group=as.factor(focal),
               color=as.factor(focal))) +
    geom_point() + 
    geom_line(size=2) +
    scale_color_paletteer_d("MetBrewer::Benedictus") +
    #scale_color_paletteer_d("MoMAColors::Avedon") +
    #scale_color_manual(values = c(wes_palette("Zissou1",
    #                                       12,type = "continuous")))+
    theme_bw() 
  
}
library(paletteer)
sem.trait.plotlist[["aus"]]
sem.trait.plotlist[["spain"]]
