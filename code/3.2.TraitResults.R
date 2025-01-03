
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
#---- 1.Looking at Interactions for answers---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.1. Make data df ----
country = "aus"
Cool.trait.df <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.Year.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  Focal.group.df <- Realised.Int.Obs.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.30 ~ "Receives Fac",
                                 sigmoid > 0.70 ~ "Receives Comp",
                               T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid")
  
  Neigh.group.df <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh, function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.30 ~ "Gives Fac",
                                 sigmoid > 0.70 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")

    
    
  if(country=="aus"){
    Neigh.group.df <-Neigh.group.df %>%
      mutate(neigh.org = case_when( sigmoid.neigh < 0.25 ~ "Gives Fac",
                                    sigmoid.neigh > 0.50 ~ "Gives Comp",
                                   T ~ "Gives both"))
    
    Focal.group.df <-   Focal.group.df %>%
      mutate(focal.org = case_when( sigmoid.focal< 0.25 ~ "Receives Fac",
                                    sigmoid.focal > 0.50 ~ "Receives Comp",
                                 T ~ "Receives both"))
  }
    
trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
library(vegan)
specific.trait.dist  <- NULL
for( i in names(trait.df)){
  #if( i =="coord.slow.to.fast") next
  specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
    as.data.frame() %>%
    gather(.,key="neigh",value="trait.dist") %>%
    mutate(trait.dist=scale(trait.dist),
           trait=i,
           neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
           focal= rep(rownames(trait.df),times=length(Code.focal.list)))
  # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
  #                                                       scale() %>% as.data.frame(),
  #                                                  na.rm = T,method="euclidean",diag=T))) %>%
  
  specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
}
Sym.group.df  <- NULL
for(i in 1:11){
  for(j in (i+1):10){
    df.focal.i.n <- Realised.Int.Obs.list[[country]] %>%
      dplyr::filter(focal==Code.focal.list[i] & neigh==Code.focal.list[j])
    prop.comp.on.i <- sum(df.focal.i.n$sigmoid < 0)/length(df.focal.i.n$sigmoid)
    df.focal.j.n <- Realised.Int.Obs.list[[country]] %>%
      dplyr::filter(focal==Code.focal.list[j] & neigh==Code.focal.list[i])
    prop.comp.on.j <- sum(df.focal.j.n$sigmoid < 0)/length(df.focal.j.n$sigmoid)
    Sym.group.df.n <- data.frame(focal=Code.focal.list[i] , neigh=Code.focal.list[j],
              prop.comp.on.i=prop.comp.on.i,prop.comp.on.j=prop.comp.on.j) %>%
      left_join(specific.trait.dist)
    
   Sym.group.df <- bind_rows( Sym.group.df, Sym.group.df.n)
   }
}
Sym.group.df <-  Sym.group.df %>%
  mutate(sym = case_when((prop.comp.on.i < 0.5 & prop.comp.on.j < 0.5) ~ "++",
                         (prop.comp.on.i > 0.5 & prop.comp.on.j > 0.5) ~ "--",
                         (prop.comp.on.i > 0.5 & prop.comp.on.j < 0.5) ~ "-+",
                         (prop.comp.on.i < 0.5 & prop.comp.on.j > 0.5) ~ "+-")) %>%
  mutate(trait.dist = case_when(sym == "-+" ~ -trait.dist,
                                T~ trait.dist))
  
  
if(country=="aus"){
specific.trait.dist <- specific.trait.dist %>%
  mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root biomass","Root volume","coord.slow.to.fast") ~ "1.BelowGround",
                                     trait %in% c("SLA","Stem height","Canopy Area") ~ "2.Aboveground",
                                     trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
}
if(country=="spain"){
  specific.trait.dist <- specific.trait.dist %>%
    mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","coord.slow.to.fast") ~ "1.BelowGround",
                                       trait %in% c("SLA","Water use efficiency","Canopy shape","Stem length","Ratio leafs","Leaf area",
                                                    "Leaf nitrogen cc") ~ "2.Aboveground",
                                       trait %in% c("Mean fecundity") ~ "3.Reproduction"))
}
trait.dist.df <- Realised.Int.Obs.list[[country]] %>%
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
            relationship ="many-to-many")

sum.trait.dist.median.df <-  trait.dist.df %>%
  mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                               Sigm.ratio< 0.5 ~"Fac",
                               Sigm.ratio== 0.5 ~"Neutre")) %>%
  group_by(compORfac,trait) %>%
  summarise(weighted.median = matrixStats::weightedMedian(trait.dist,Sigm.ratio,na.rm=T),
            weighted.mad = matrixStats::weightedMad(trait.dist,Sigm.ratio,na.rm=T)) %>%
  mutate(weightedQ1= weighted.median-weighted.mad,
         weightedQ9= weighted.median+weighted.mad)

sum.trait.dist.focal.df <- trait.dist.df %>%
  group_by(focal.org,trait) %>%
  summarise(median.trait.dist=median(trait.dist,na.rm=T),
            Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
            Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
            median.sigmoid=median(Sigm.Q5,na.rm=T),
            Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
            Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))%>%
  ungroup()
sum.trait.dist.neigh.df <- trait.dist.df %>%
  group_by(neigh.org,trait) %>%
  summarise(median.trait.dist=median(trait.dist,na.rm=T),
            Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
            Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
            median.sigmoid=median(Sigm.Q5,na.rm=T),
            Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
            Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))

Cool.trait.df[[country]] <- list(
  trait.dist.df=trait.dist.df,
  sum.trait.dist.focal.df=sum.trait.dist.focal.df,
  sum.trait.dist.neigh.df=sum.trait.dist.neigh.df,
  specific.trait.dist=specific.trait.dist,
  sum.trait.dist.median.df=sum.trait.dist.median.df,
  Sym.group.df=Sym.group.df)
}

#---- 1.2. Detailed plot of dist and median interaction ----
Cool.detailed.trait.plotlist <- list()
for( country in country.list){
  
Cool.detailed.trait.plotlist[[paste0(country,"_focal")]] <- ggplot() +
  geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
             aes(y=trait.dist,
                 x=Sigm.Q5,
                 group=as.factor(focal.org),
                 fill=Sigm.ratio,
                 color=as.factor(focal.org)),
             shape=21,
             size=1,
             alpha=0.2) + 
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df,
                  aes(y=median.trait.dist,
                      x=median.sigmoid,
                      xmin=Q10.sigmoid,
                      xmax=Q90.sigmoid,
                      color=as.factor(focal.org)),
                  shape=15,
                  size=1) + 
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df,
                  aes(y=median.trait.dist,
                      x=median.sigmoid,
                      ymin=Q10.trait.dist,
                      ymax=Q90.trait.dist,
                      color=as.factor(focal.org)),
                  shape=15,
                  size=1) +
  facet_wrap(.~trait,scale="free") + 
  labs(y="Focal trait value - Neigh trait value",
       x="Median sigmoid value given by Neigh")+
  scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                             101, 
                                             type = "continuous"))+
  scale_color_manual(values = c("#EBCC2A","#F21A00","#3B9AB2"))+
  theme_bw() 



Cool.detailed.trait.plotlist[[paste0(country,"_neigh")]]  <- ggplot() +
  geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
             aes(y=trait.dist,
                 x=Sigm.Q5,
                group=as.factor(neigh.org),
                fill=Sigm.ratio,
                color=as.factor(neigh.org)),
                shape=21,
                size=1,
                alpha=0.2) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        xmin=Q10.sigmoid,
                        xmax=Q90.sigmoid,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) + 
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df,
                  aes(y=median.trait.dist,
                      x=median.sigmoid,
                      ymin=Q10.trait.dist,
                      ymax=Q90.trait.dist,
                      color=as.factor(neigh.org)),
                  shape=15,
                  size=1) +
  facet_wrap(.~trait,scale="free") + 
  labs(y="Focal trait value - Neigh trait value",
       x="Median sigmoid value given by Neigh")+
  scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                             101, 
                                             type = "continuous"))+
  scale_color_manual(values = c("#EBCC2A","#F21A00","#3B9AB2"))+
  theme_bw() 
}

Cool.detailed.trait.plotlist[[paste0("spain_focal")]]
Cool.detailed.trait.plotlist[[paste0("spain_neigh")]]
Cool.detailed.trait.plotlist[[paste0("aus_focal")]]
Cool.detailed.trait.plotlist[[paste0("aus_neigh")]]

#---- 1.3. Stat Test ----
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
      dplyr::select(trait,compORfac,trait.dist)
    median.comp.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Comp" &
                                                                  sum.trait.dist.median.df$trait==trait.i]
    median.fac.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Fac" &
                                                                  sum.trait.dist.median.df$trait==trait.i]
    
    if(median.comp.n > median.fac.n){test.n ="greater"}else{test.n="less"}
    
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
#---- 1.4. Summary plot with dist - rect ----
Cool.rect.trait.plotlist <- list()
country="aus"
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
for( country in country.list){
  Test.focal.trait.df <-   Test.trait.list[[country]]$Test.focal.trait.df %>%
    dplyr::filter(p.value< 0.05) 
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  sum.trait.dist.focal.df<- Cool.trait.df[[country]]$sum.trait.dist.focal.df
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
Focal_sum <- ggplot()+
  geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                   mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
                   mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
             aes(x=trait.dist,  
                 y=trait,
                 fill= sigmoid.focal,
                 color=focal.org),
             position = position_dodge(width = 0.7),
             shape=21,size=2,alpha=0.2) +
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df %>%
                    mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))),
                  position = position_dodge(width = 0.7),
                  aes(x= median.trait.dist,
                      y=trait,
                      xmin=Q10.trait.dist,
                      xmax=Q90.trait.dist,
                      color=as.factor(focal.org)),
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
       color="Species grouping based on received interactions")+
  scale_color_manual(values = c("#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
  scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                          101, 
                                          type = "continuous"))+
  scale_alpha_continuous(range=c(0.1,1)) +
  theme_bw()  +
  guides(color = guide_legend(title.position = "top")) + 
  theme(legend.position="bottom",
        axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                             "bold","plain")))
Focal_sum 

Test.neigh.trait.df <-   Test.trait.list[[country]]$Test.neigh.trait.df %>%
  dplyr::filter(p.value< 0.05) 

Neigh_sum <- ggplot()+
  geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp")))%>%
                mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
             aes(x=-trait.dist, # reverse the distance to values on the right mean higher value of the trait  
                 y=trait,
                 fill= sigmoid.neigh,
                 color=as.factor(neigh.org)),
             position = position_dodge(width = 0.7),
             shape=21,size=2,alpha=0.2) +
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df %>%
                    mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))), # "Gives both"
                  position = position_dodge(width = 0.7),
                  aes(x=-median.trait.dist,
                      y=trait,
                      xmin=-Q10.trait.dist,
                      xmax=-Q90.trait.dist,
                      color=as.factor(neigh.org)),
                  shape=15,
                  size=1) +
  scale_y_discrete(expand=c(0.1,0.1))+
  annotate("text", x = 2, y=0.2, 
           label = "Neigh has higher \n trait value than focal",
           size=4) + 
  geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
               arrow = arrow(length = unit(0.2, "cm")))+
  annotate("text", x = -2, y=0.2, 
           label = "Neigh has lower \n trait value than focal",
           size=4) + 
  geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
               arrow = arrow(length = unit(0.2, "cm")))+
  labs(x="Focal trait value - Neigh trait value",
       y="trait",
       color="Species grouping based on given interactions")+
  scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                          101, 
                                          type = "continuous"))+
  scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
  theme_bw() +
  guides(color = guide_legend(title.position = "top")) + 
  theme(legend.position="bottom",
        axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.neigh.df$trait)) %in% levels(as.factor(Test.neigh.trait.df$trait)),
                                             "bold","plain")))



Cool.rect.trait.plotlist[[paste0(country,"_sum")]] <-
  ggarrange(Focal_sum,Neigh_sum,
            nrow=1,
            common.legend = F)

Test.trait.df <-   Test.trait.list[[country]]$Test.trait.df %>%
  dplyr::filter(p.value <= 0.1) 
Cool.rect.trait.plotlist[[paste0(country,"_Pairwise_interactions")]] <- ggplot()+
  geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
               mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                            Sigm.ratio< 0.5 ~"Fac",
                                            Sigm.ratio== 0.5 ~"Neutre")) %>%
               dplyr::filter(!compORfac =="Neutre") %>%
               mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
               mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
             aes(x=trait.dist,  
                 y=trait,
                 #alpha=Sigm.ratio,
                 #fill=compORfac, # focal.org,
                 color=compORfac), # focal.org,
             position = position_dodge(width = 0.7),
             shape=16,size=2,alpha=0.2) +
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df %>%
                    dplyr::filter(!compORfac =="Neutre"), 
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
  #scale_fill_gradientn(colours = wes_palette("Zissou1", 
  #                                       101, 
  #                                    type = "continuous"))+
  scale_alpha_continuous(range=c(0.1,1)) +
  theme_bw()  +
  guides(color = guide_legend(title.position = "top")) + 
  theme(legend.position="bottom",
        legend.title =element_text(size=16),
        legend.text =element_text(size=14),
        axis.title=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.trait.df$trait)),
                                             "bold","plain"),
                                 size=12))

}
Cool.rect.trait.plotlist[["aus_sum"]] 
Cool.rect.trait.plotlist[["spain_sum"]] 
#---- 1.5. Summary plot with dist - circ ----
Cool.circ.trait.plotlist <- list()
country="aus"
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
   Test.focal.trait.df <-   Test.trait.list[[country]]$Test.focal.trait.df %>%
    dplyr::filter(p.value< 0.05)
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  sum.trait.dist.focal.df<- Cool.trait.df[[country]]$sum.trait.dist.focal.df
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  
Focal_sum_circ <- ggplot()+
  geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
               mutate(focal = factor(focal, levels= Focal.group.df$focal[order(Focal.group.df$sigmoid)]))%>%
               mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
               mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
              aes(y=trait.dist,  
                 x=trait,
                 fill=sigmoid.focal,
                 color=focal.org),
             position = position_dodge(width = 0.7),
             shape=21,size=2,alpha=0.2) +
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df %>%
                    mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))), # "Receives both",
                  position = position_dodge(width = 0.7),
                  aes(y=median.trait.dist,
                      x=trait,
                      ymin=Q10.trait.dist,
                      ymax=Q90.trait.dist,
                      color=as.factor(focal.org)),
                  shape=16,
                  size=1) +
  geom_hline(yintercept=0,color="black") +
  annotate("text",x =3.4, y = 1, label = "Higher trait value",
           angle =12,color = "gray12",size = 3.5)+ # spain: 14, aus: 6
  geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 2),
                arrow = arrow(length = unit(0.5, "cm")))+
  annotate("text",x =3.6, y = -1, label = "Lower trait value",
           angle = 12,color = "gray12",size = 3.5)+
  geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = -2),
               arrow = arrow(length = unit(0.5, "cm")))+
  scale_y_continuous(breaks=c(0)) + 
  scale_x_discrete(labels=addline_format(unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))) +
  labs(fill="Percentage of Comp given",
       color="Species grouping based on given interactions")+
  scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                            101, 
                                            type = "continuous"))+
  scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
  coord_polar() +
  guides(color = guide_legend(title.position = "top"),
         fill = guide_legend(title.position = "top")) + 
  theme_bw() +
  theme(legend.position="bottom",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        # Use gray text for the region names
        axis.text.x = element_text(color = "gray12", size = 12,
                                   face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                                "bold","plain")),
        panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank()) 
Focal_sum_circ 

Test.neigh.trait.df  <-   Test.trait.list[[country]]$Test.neigh.trait.df %>%
  dplyr::filter(p.value <= 0.05)
  
Neigh_sum_circ <- ggplot()+
  geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
               #mutate(neigh= factor(neigh, levels= Neigh.group.df$neigh[order(Neigh.group.df$sigmoid)])) %>%
               mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))) %>%
               mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
              aes(y=-trait.dist,
                 x=trait,
                 fill= sigmoid.neigh,
                 color=as.factor(neigh.org)),
             position = position_dodge(width = 0.7),
             shape=21,size=2,alpha=0.2) +
  geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df %>%
                    mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))), # "Gives both"
                  position = position_dodge(width = 0.7),
                  aes(y=-median.trait.dist,
                      x=trait,
                      ymin=-Q10.trait.dist,
                      ymax=-Q90.trait.dist,
                      color=as.factor(neigh.org)),
                  shape=15,
                  size=1) +
  geom_hline(yintercept=0,color="black") +
  scale_y_continuous(breaks=c(0)) + 
  scale_x_discrete(labels=addline_format(unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))) +
   labs(fill="Percentage of Comp given",
        color="Species grouping based on given interactions")+
  scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                            101, 
                                            type = "continuous"))+
  scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
  coord_polar() +
  guides(color = guide_legend(title.position = "top"),
        fill = guide_legend(title.position = "top")) + 
  theme_bw() +
  theme(legend.position="bottom",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        # Use gray text for the region names
        panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x=element_text(color = "gray12", size = 12,
                                 face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.neigh.df$trait)) %in% levels(as.factor(Test.neigh.trait.df$trait)),
                                             "bold","plain")))
Neigh_sum_circ

Cool.circ.trait.plotlist[[paste0(country,"_circ")]] <-
  ggarrange(Focal_sum_circ,Neigh_sum_circ,
            nrow=1,
            common.legend = F, labels=c("a. Species as receiver","b. Species as impactor"))


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
               mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
               mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
             aes(y=trait.dist,  
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
  geom_text(data=Test.trait.df,aes(x=trait,y=3,label=labels.pvalue),size=12)+
  geom_hline(yintercept=0,color="black") +
  scale_y_continuous(breaks=c(0)) + 
  scale_x_discrete(labels=addline_format(unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))) +
  labs(color="Pairwise interaction mainly")+
  scale_color_manual(values = c("#F21A00","#3B9AB2"),
                     labels=c("Competitive","Facilitative"))+ # "#EBCC2A",
  coord_polar() +
  guides(color = guide_legend(title.position = "top"),
         fill = guide_legend(title.position = "top")) + 
  theme_bw() +
  theme(legend.position="bottom",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        # Use gray text for the region names
        panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title =element_text(size=16),
        legend.text =element_text(size=14),
        axis.text.x=element_text(color = "gray12", size = 12,
                                 face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.neigh.df$trait)) %in% levels(as.factor(Test.trait.df$trait)),
                                             "plain","plain")))
Cool.circ.trait.plotlist[[paste0(country,"_Pairwise_interactions")]]

}
Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]] #figures/Cool.circ.pairwise.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]] # figures/Cool.circ.pairwise.trait.aus.pdf

Cool.circ.trait.plotlist[["Pairwise_interactions"]] <-
  ggarrange(Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]],
            Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]],
            nrow=1,legend="bottom",
            common.legend = T, labels=c("a. Spain","b. Australia"))

Cool.circ.trait.plotlist["Pairwise_interactions"] # figures/Cool.circ.pairwise.trait.pdf

Cool.circ.trait.plotlist[[paste0("spain_circ")]] # figures/Cool.circ.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_circ")]] # figures/Cool.circ.trait.aus.pdf


#---- 1.6. Summary plot with symmetry----
Sym.trait.plotlist <- list()
country="aus"
for( country in country.list){
  Sym.group.df <- Cool.trait.df[[country]]$Sym.group.df %>%
    mutate(sym = case_when(sym == "-+" ~ "+-",
                           T ~ sym)) %>%
    dplyr::filter(!is.na(sym))
  Sym.trait.plotlist[[country]] <- ggplot() +
    geom_point(data=Sym.group.df,
               aes(x=abs(trait.dist),  
                   y=trait,
                   color=as.factor(sym)),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.8) +
    stat_summary(data=Sym.group.df ,
                    position = position_dodge(width = 0.7),
                    aes(x=abs(trait.dist),
                        y=trait,group=sym,color=as.factor(sym)),
                    fun=mean, 
                    fun.max = function(x) min(3,mean(x) + var(x)), #quantile(x,c(0.9)),
                    fun.min = function(x)  max(0,mean(x) - var(x)), #quantile(x,c(0.1)),
                    geom="pointrange",
                    shape=16,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1)) +
    scale_x_continuous(limits=c(0,3)) +
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Sym")+
    scale_color_manual(values = rev(c("#3B9AB2","#EBCC2A","#F21A00")))+ # "#EBCC2A",
    theme_bw()  +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                               "bold","plain")))
  
  
}
Sym.trait.plotlist[["aus"]]
Sym.trait.plotlist[["spain"]]
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.Looking at Interactions across time for answers---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1. Make data df ----
country = "spain"
Cool.trait.df <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.Year.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  Focal.group.df <- Realised.Int.Obs.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.30 ~ "Receives Fac",
                                 sigmoid > 0.70 ~ "Receives Comp",
                                 T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid") 

  
  Neigh.group.df <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh , function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.30 ~ "Gives Fac",
                                 sigmoid > 0.70 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")

  if(country=="aus"){
    Neigh.group.df <-Neigh.group.df %>%
      mutate(neigh.org = case_when( sigmoid.neigh < 0.25 ~ "Gives Fac",
                                    sigmoid.neigh > 0.50 ~ "Gives Comp",
                                    T ~ "Gives both"))
    
    Focal.group.df <-   Focal.group.df %>%
      mutate(focal.org = case_when( sigmoid.focal< 0.25 ~ "Receives Fac",
                                    sigmoid.focal > 0.50 ~ "Receives Comp",
                                    T ~ "Receives both"))
  }
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    #if( i =="coord.slow.to.fast") next
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
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root biomass","Root volume","coord.slow.to.fast") ~ "1.BelowGround",
                                         trait %in% c("SLA","Stem height","Canopy width","Canopy width 90deg") ~ "2.Aboveground",
                                         trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","coord.slow.to.fast") ~ "1.BelowGround",
                                         trait %in% c("SLA","Water use efficiency","Canopy shape","Stem length","Ratio leafs","Leaf area",
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
              relationship ="many-to-many")
  if(country=="spain"){
    trait.dist.df <- trait.dist.df%>% dplyr::filter(!focal=="PAIN") %>%
      dplyr::filter(!neigh=="PAIN") 
  }
  sum.trait.dist.focal.df <- trait.dist.df %>%
    group_by(focal.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))%>%
    ungroup()
  sum.trait.dist.neigh.df <- trait.dist.df %>%
    group_by(neigh.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))
  
  Cool.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    sum.trait.dist.focal.df=sum.trait.dist.focal.df,
    sum.trait.dist.neigh.df=sum.trait.dist.neigh.df,
    specific.trait.dist=specific.trait.dist)
}

#---- 2.2. Detailed plot of dist and median interaction ----
Cool.detailed.trait.plotlist <- list()
for( country in country.list){
  
  Cool.detailed.trait.plotlist[[paste0(country,"_focal")]] <- ggplot() +
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
               aes(y=trait.dist,
                   x=Sigm.Q5,
                   group=as.factor(focal.org),
                   fill=Sigm.ratio,
                   color=as.factor(focal.org)),
               shape=21,
               size=1,
               alpha=0.2) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        xmin=Q10.sigmoid,
                        xmax=Q90.sigmoid,
                        color=as.factor(focal.org)),
                    shape=15,
                    size=1) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        ymin=Q10.trait.dist,
                        ymax=Q90.trait.dist,
                        color=as.factor(focal.org)),
                    shape=15,
                    size=1) +
    facet_wrap(.~trait,scale="free") + 
    labs(y="Focal trait value - Neigh trait value",
         x="Median sigmoid value given by Neigh")+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual(values = c("#EBCC2A","#F21A00","#3B9AB2"))+
    theme_bw() 
  
  
  
  Cool.detailed.trait.plotlist[[paste0(country,"_neigh")]]  <- ggplot() +
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
               aes(y=trait.dist,
                   x=Sigm.Q5,
                   group=as.factor(neigh.org),
                   fill=Sigm.ratio,
                   color=as.factor(neigh.org)),
               shape=21,
               size=1,
               alpha=0.2) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        xmin=Q10.sigmoid,
                        xmax=Q90.sigmoid,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        ymin=Q10.trait.dist,
                        ymax=Q90.trait.dist,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) +
    facet_wrap(.~trait,scale="free") + 
    labs(y="Focal trait value - Neigh trait value",
         x="Median sigmoid value given by Neigh")+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual(values = c("#EBCC2A","#F21A00","#3B9AB2"))+
    theme_bw() 
}

Cool.detailed.trait.plotlist[[paste0("spain_focal")]]
Cool.detailed.trait.plotlist[[paste0("spain_neigh")]]
Cool.detailed.trait.plotlist[[paste0("aus_focal")]]
Cool.detailed.trait.plotlist[[paste0("aus_neigh")]]

#---- 2.3. Stat Test ----
Test.trait.df <- list()
for( country in country.list){
  trait.dist.df <- Cool.trait.df[[country]]$trait.dist.df
  Test.neigh.trait.df <- NULL
  Test.focal.trait.df <- NULL
  for( trait.i in levels(as.factor(trait.dist.df$trait))){
    if( trait.i =="coord.slow.to.fast") next
    for( i in levels(as.factor(trait.dist.df$focal.org))){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(!focal.org==i )  %>%
        dplyr::select(trait,focal.org,trait.dist)
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
                                    Test.neigh.trait.df=Test.neigh.trait.df)
}
#---- 2.4. Summary plot with dist - rect ----
Cool.rect.trait.plotlist <- list()
country="aus"
for( country in country.list){
  Test.focal.trait.df <-   Test.trait.list[[country]]$Test.focal.trait.df %>%
    dplyr::filter(p.value< 0.05) 
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  sum.trait.dist.focal.df<- Cool.trait.df[[country]]$sum.trait.dist.focal.df
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  Focal_sum <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(x=trait.dist,  
                   y=trait,
                   fill=sigmoid.focal,
                   color=focal.org),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df %>%
                      mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))),
                    position = position_dodge(width = 0.7),
                    aes(x=median.trait.dist,
                        y=trait,
                        xmin=Q10.trait.dist,
                        xmax=Q90.trait.dist,
                        color=as.factor(focal.org)),
                    shape=16,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1))+
    annotate("text", x = 2, y=0.2, 
             label = "Focal has higher \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    annotate("text", x = -2, y=0.2, 
             label = "Focal has lower \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Species grouping based on received interactions")+
    scale_color_manual(values = c("#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    theme_bw()  +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                               "bold","plain")))
  
  Test.neigh.trait.df <-   Test.trait.list[[country]]$Test.neigh.trait.df %>%
    dplyr::filter(p.value< 0.05) 
  
  Neigh_sum <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp")))%>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(x=-trait.dist, # reverse the distance to values on the right mean higher value of the trait  
                   y=trait,
                   fill= sigmoid.neigh,
                   color=as.factor(neigh.org)),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df %>%
                      mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))), # "Gives both"
                    position = position_dodge(width = 0.7),
                    aes(x=-median.trait.dist,
                        y=trait,
                        xmin=-Q10.trait.dist,
                        xmax=-Q90.trait.dist,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1))+
    annotate("text", x = 2, y=0.2, 
             label = "Neigh has higher \n trait value than focal",
             size=4) + 
    geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    annotate("text", x = -2, y=0.2, 
             label = "Neigh has lower \n trait value than focal",
             size=4) + 
    geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Species grouping based on given interactions")+
    scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                              101, 
                                              type = "continuous"))+
    scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    theme_bw() +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.neigh.df$trait)) %in% levels(as.factor(Test.neigh.trait.df$trait)),
                                               "bold","plain")))
  
  
  
  Cool.rect.trait.plotlist[[paste0(country,"_sum")]] <-
    ggarrange(Focal_sum,Neigh_sum,
              nrow=1,
              common.legend = F)
}
Cool.rect.trait.plotlist[["aus_sum"]] 
Cool.rect.trait.plotlist[["spain_sum"]] 
#---- 2.5. Summary plot with dist - circ ----
Cool.circ.trait.plotlist <- list()
country="aus"
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
  Test.focal.trait.df <-   Test.trait.list[[country]]$Test.focal.trait.df %>%
    dplyr::filter(p.value< 0.05)
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  sum.trait.dist.focal.df<- Cool.trait.df[[country]]$sum.trait.dist.focal.df
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  
  Focal_sum_circ <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(focal = factor(focal, levels= Focal.group.df$focal[order(Focal.group.df$sigmoid)]))%>%
                 mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(y=trait.dist,  
                   x=trait,
                   fill=sigmoid.focal,
                   color=focal.org),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df %>%
                      mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))), # "Receives both",
                    position = position_dodge(width = 0.7),
                    aes(y=median.trait.dist,
                        x=trait,
                        ymin=Q10.trait.dist,
                        ymax=Q90.trait.dist,
                        color=as.factor(focal.org)),
                    shape=16,
                    size=1) +
    geom_hline(yintercept=0,color="black") +
    annotate("text",x =3.4, y = 1, label = "Higher trait value",
             angle = 14,color = "gray12",size = 3.5)+ # spain: 14, aus: 6
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 2),
                 arrow = arrow(length = unit(0.5, "cm")))+
    annotate("text",x =3.6, y = -1, label = "Lower trait value",
             angle = 14,color = "gray12",size = 3.5)+
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = -2),
                 arrow = arrow(length = unit(0.5, "cm")))+
    scale_y_continuous(breaks=c(0)) + 
    scale_x_discrete(labels=addline_format(unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))) +
    labs(fill="Percentage of Comp given",
         color="Species grouping based on given interactions")+
    scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                              101, 
                                              type = "continuous"))+
    scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    coord_polar() +
    guides(color = guide_legend(title.position = "top"),
           fill = guide_legend(title.position = "top")) + 
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          axis.text.x = element_text(color = "gray12", size = 12,
                                     face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                                 "bold","plain")),
          panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank()) 
  Focal_sum_circ 
  
  Test.neigh.trait.df  <-   Test.trait.list[[country]]$Test.neigh.trait.df %>%
    dplyr::filter(p.value< 0.05)  
  
  Neigh_sum_circ <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 #mutate(neigh= factor(neigh, levels= Neigh.group.df$neigh[order(Neigh.group.df$sigmoid)])) %>%
                 mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))) %>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(y=-trait.dist,
                   x=trait,
                   fill= sigmoid.neigh,
                   color=as.factor(neigh.org)),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df %>%
                      mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))), # "Gives both"
                    position = position_dodge(width = 0.7),
                    aes(y=-median.trait.dist,
                        x=trait,
                        ymin=-Q10.trait.dist,
                        ymax=-Q90.trait.dist,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) +
    geom_hline(yintercept=0,color="black") +
    scale_y_continuous(breaks=c(0)) + 
    scale_x_discrete(labels=addline_format(unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))) +
    labs(fill="Percentage of Comp given",
         color="Species grouping based on given interactions")+
    scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                              101, 
                                              type = "continuous"))+
    scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    coord_polar() +
    guides(color = guide_legend(title.position = "top"),
           fill = guide_legend(title.position = "top")) + 
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x=element_text(color = "gray12", size = 12,
                                   face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.neigh.df$trait)) %in% levels(as.factor(Test.neigh.trait.df$trait)),
                                               "bold","plain")))
  Neigh_sum_circ
  
  Cool.circ.trait.plotlist[[paste0(country,"_circ")]] <-
    ggarrange(Focal_sum_circ,Neigh_sum_circ,
              nrow=1,
              common.legend = F, labels=c("a. Species as receiver","b. Species as impactor"))
}

Cool.circ.trait.plotlist[[paste0("spain_circ")]] # figures/Cool.circ.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_circ")]] # figures/Cool.circ.trait.aus.pdf


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.Traits Graph related- Looking at traits for answers---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.1 With different grouping of functional group ----
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

#---- 3.2 Along the fast to slow gradient ----
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

#---- 3.3 Along the dist between interactors----
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

#---- 3.4 Along the trait dist between interactors----
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
#---- 3.5 Sem -----
sem.trait.plotlist <- list()
country = "spain"
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
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  
  sem.df <- NULL
  for( sp in Code.focal.list){
    print(sp)
    test.df <-  Realised.Int.Obs.list[[country]] %>%
      left_join(specific.trait.dist,
                by=c("focal","neigh"), 
                relationship = "many-to-many") %>%
      #filter(focal ==sp) %>%
      filter(focal ==sp) %>%
      filter(!focal == neigh) %>%
      aggregate(sigmoid ~ year + density + focal + neigh + trait.dist + trait, median) %>%
      group_by(year, focal, neigh, trait.dist,trait) %>%
      summarise(sigmoid = median(sigmoid),
                density=median(density)) %>%
      ungroup() %>%
      mutate(trait.dist =abs(trait.dist))%>%
      spread(trait,trait.dist) %>%
      left_join(env_pdsi %>%  mutate(year = as.numeric(year)),
             by=c("year"),multiple = "all")  %>%
      as.data.frame()
    names(test.df) <- gsub(" ",".",colnames(test.df))
    
    if(country=="aus"){
      test.df <- test.df  %>%
        dplyr::select_if(~ sum(!is.na(.))>10) %>%
        drop_na()
      var.names <- names(test.df)[names(test.df) %in% c("sla","width.longest","srl","root.length",
                                                        "flower.size.numb","height","number.of.root.tips")]
      
    }else{  test.df <- test.df %>%select_if(~ sum(!is.na(.))>10)
    var.names <- names(test.df)[names(test.df) %in% c("Canopy.shape","coord.slow.to.fast","Leaf.area",
                                                      "Leaf.area.index",#"Mean.fecundity","Ratio.leafs","Root.diameter" ,"Water.use.efficiency","SRA","SRL","Stem.length"
                                                      "Leaf.nitrogen.cc" , "Root.mass.density")]
    }
    hist(test.df$sigmoid,breaks=150)
    sem.median  <- piecewiseSEM::psem(
      glm(as.formula(paste("sigmoid ~ 1 +", paste(c(var.names,"density","Precip.extrem"), collapse= "+"))),
         data=test.df, family= gaussian(link = "identity")),
      glm(density ~ Precip.extrem,
                   data=test.df,family= gaussian(link = "log")),
      test.df)
    
    #plot.psemhl(sem.median , correlation = T, layout = "tree") 
    summary(sem.median)
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


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.Traits Graph with time ---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1. Make data df ----
country = "spain"
Cool.trait.df <- list()
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
  
  Focal.group.df <- Realised.Int.Obs.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal + year, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.30 ~ "Receives Fac",
                                 sigmoid > 0.70 ~ "Receives Comp",
                                 T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid") %>%
    mutate(year=as.integer(year)) %>%
    left_join(env_pdsi%>% mutate(year=as.integer(year)) )
  hist(log(Focal.group.df$sigmoid.focal))
  test.glm <- glmer.nb(sigmoid.focal ~ Precip.extrem + (1|focal) ,data=Focal.group.df)
  summary(test.glm)
  ggplot(Focal.group.df,aes(y=sigmoid.focal, x=Precip.extrem,color=focal)) +
    geom_point()+
    geom_smooth(method="lm",se=F)
  
  Neigh.group.df <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh + year, function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.30 ~ "Gives Fac",
                                 sigmoid > 0.70 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")%>%
    mutate(year=as.integer(year)) %>%
    left_join(env_pdsi%>% mutate(year=as.integer(year)) )
  
  ggplot(Neigh.group.df,aes(y=sigmoid.neigh, x=Precip.extrem)) + geom_point()+
    geom_smooth(method="lm",se=F)
  
  if(country=="aus"){
    Neigh.group.df <-Neigh.group.df %>%
      mutate(neigh.org = case_when( sigmoid.neigh < 0.25 ~ "Gives Fac",
                                    sigmoid.neigh > 0.50 ~ "Gives Comp",
                                    T ~ "Gives both"))
    
    Focal.group.df <-   Focal.group.df %>%
      mutate(focal.org = case_when( sigmoid.focal< 0.25 ~ "Receives Fac",
                                    sigmoid.focal > 0.50 ~ "Receives Comp",
                                    T ~ "Receives both"))
  }
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    #if( i =="coord.slow.to.fast") next
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
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
      mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root biomass","Root volume","coord.slow.to.fast") ~ "1.BelowGround",
                                         trait %in% c("SLA","Stem height","Canopy width","Canopy width 90deg") ~ "2.Aboveground",
                                         trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","coord.slow.to.fast") ~ "1.BelowGround",
                                         trait %in% c("SLA","Water use efficiency","Canopy shape","Stem length","Ratio leafs","Leaf area",
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
              relationship ="many-to-many")
  if(country=="spain"){
    trait.dist.df <- trait.dist.df%>% dplyr::filter(!focal=="PAIN") %>%
      dplyr::filter(!neigh=="PAIN") 
  }
  sum.trait.dist.focal.df <- trait.dist.df %>%
    group_by(focal.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))%>%
    ungroup()
  sum.trait.dist.neigh.df <- trait.dist.df %>%
    group_by(neigh.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))
  
  Cool.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    sum.trait.dist.focal.df=sum.trait.dist.focal.df,
    sum.trait.dist.neigh.df=sum.trait.dist.neigh.df,
    specific.trait.dist=specific.trait.dist)
}
