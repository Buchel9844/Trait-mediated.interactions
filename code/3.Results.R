#---- Exploration of general results. 
# Graphs of abundances and raw sigmoid ----

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
#---- 1. Visualisation of Species Abundance ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.0 Clean data----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))
country.list <- c("aus","spain")
#---- 1.1 Abundances over time ----
widthplot = 25
abundance_plotlist <- NULL
for(country in country.list){
  #abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".preclean")]]
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  abundance_plotlist[[country]] <- abundance_summary %>%
    #filter( individuals > 0.001) %>%
    ggplot(aes(x=as.character(year), 
               y = individuals*widthplot,
               group=as.factor(species),
               color=as.factor(species))) +
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=2) +
    stat_summary(fun.y = mean,
                 geom = "line",size=1) +
    scale_y_log10(limits = c(1,100)) +
    scale_x_discrete("year",limits= year.levels) +
    scale_color_manual(values= col.df$color.name) +
    labs(color="species",y="Mean number of \nindividuals in 25x25cm plot",
         title=paste0("Density over time of annual plants in ",country)) +
    coord_cartesian( xlim = NULL, #ylim=c(0,750),
                     expand = TRUE, default = FALSE, clip = "on") +
    #scale_color_manual(values=safe_colorblind_palette) +
    theme_bw() +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=28),
           legend.text=element_text(size=20),
           legend.title=element_text(size=20),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_text(size=20, angle=66, hjust=1),
           axis.text.y= element_text(size=20),
           axis.title.x= element_text(size=24),
           axis.title.y= element_text(size=24),
           title=element_text(size=16))
}
#abundance_plotlist[["spain"]] # figures/Abundance_spain.pdf
#abundance_plotlist[["aus"]] #figures/Abundance_aus.pdf

#---- 1.2 Abundances over PDSI ----
widthplot = 25
abundance_pdsi_plotlist <- NULL
for(country in country.list){
  
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  env_pdsi <- read.csv(paste0(home.dic,"results/",country,"_env_pdsi.csv")) %>%
    group_by(year) %>% 
    mutate(PDSI.mean = mean(get(paste0(country,"_pdsi")), na.rm = T),
           PDSI.sd = sd(get(paste0(country,"_pdsi")),na.rm = T)) %>%
    dplyr::select(year,PDSI.mean,PDSI.sd) %>%
    filter(year %in%  year.levels) %>%
    unique() %>%
    ungroup() %>%
    mutate(PDSI.max = max(PDSI.mean)) %>%
    mutate_at("PDSI.mean",function(x) (x/ max(abs(x)))) %>% # scaling 
    as.data.frame()
  
  abundance_pdsi_plotlist[[country]] <- abundance_summary %>%
    mutate(year = as.numeric(year))%>%
    left_join(env_pdsi) %>%
    ggplot() +
    stat_summary(aes(x=exp(PDSI.mean), y = individuals*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange",size=2) +
    stat_summary(aes(x=exp(PDSI.mean), y = individuals*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 geom = "line",size=1) +
    scale_y_log10() +
    scale_color_manual(values= col.df$color.name) +
    labs(x="PDSI (exp)",
         color="species",y="Mean number of \nindividuals in 25x25cm plot",
         title=paste0("Density across PDSI of annual plants in ",country)) +
    coord_cartesian( xlim = NULL, #ylim = c(0,200),
                     expand = TRUE, default = FALSE, clip = "on") +
    #scale_color_manual(values=safe_colorblind_palette) +
    theme_bw() +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=28),
           legend.text=element_text(size=20),
           legend.title=element_text(size=20),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_text(size=20, angle=66, hjust=1),
           axis.text.y= element_text(size=20),
           axis.title.x= element_text(size=24),
           axis.title.y= element_text(size=24),
           title=element_text(size=16))
}
#abundance_pdsi_plotlist[["spain"]] #figures/Abundance_pdsi_spain.pdf
#abundance_pdsi_plotlist[["aus"]] #figures/Abundance_pdsi_aus.pdf

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2. Environmental Gradient  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
# PDSI
env_pdsi_aus <- read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) 
env_pdsi_spain <- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv")) 

#---- 2.1. Median ----
env_pdsi_plotlist <- NULL
env_pdsi_dflist <- NULL
for(country in country.list){
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  year.levels <-levels(as.factor( abundance_summary$year))
  
  
  env_pdsi <- read.csv(paste0(home.dic,"results/",country,"_env_pdsi.csv")) %>%
    group_by(year) %>% 
    summarise(PDSI.mean = mean(get(paste0(country,"_pdsi")), na.rm = T),
              PDSI.sd = sd(get(paste0(country,"_pdsi")),na.rm = T),
              PDSI.Q1 = quantile(get(paste0(country,"_pdsi")), c(0.1), na.rm = T),
              PDSI.Q5 = quantile(get(paste0(country,"_pdsi")), c(0.5), na.rm = T),
              PDSI.Q9 = quantile(get(paste0(country,"_pdsi")),c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    filter(year %in%  year.levels) %>%
    # mutate(PDSI.max = max(PDSI.mean)) %>%
    #mutate_at("PDSI.mean",function(x) (x/ max(abs(x)))) %>% # scaling 
    as.data.frame()
  env_pdsi_dflist[[country]] <-env_pdsi
  env_pdsi_plotlist[[country]] <-env_pdsi %>%
    gather(PDSI.mean,PDSI.sd,PDSI.Q1,PDSI.Q9,PDSI.Q5,
           key="PDSI",value="value") %>%
    ggplot(aes(x=year,y=value,group=PDSI,
               color=PDSI))+
    annotate("rect",xmin=min(env_pdsi$year),
             xmax=max(env_pdsi$year),
             ymin=0,ymax=max(env_pdsi$PDSI.Q9)+0.1,
             fill="lightblue",
             alpha=0.2) +
    geom_line(aes(linetype=PDSI)) +
    scale_color_manual(values=c("black","grey","grey","grey","red")) +
    scale_linetype_manual(values=c("solid","dashed","solid","dashed","solid"))+
    coord_cartesian(expand = FALSE) + 
    theme_bw()
  
}
env_pdsi_plotlist[["aus"]] #figures/env_pdsi_aus.pdf
env_pdsi_plotlist[["spain"]] #figures/env_pdsi_spain.pdf


env_pdsi %>%
  ggplot(aes(x=PDSI.sd,y=PDSI.sd))+
  geom_line() +
  geom_point()

env_pdsi_dflist[["aus"]] %>%
  ggplot(aes(y=PDSI.sd,x=PDSI.mean))+
  geom_line() +
  geom_point()
env_pdsi_dflist[["spain"]] %>%
  ggplot(aes(y=PDSI.sd,x=PDSI.mean))+
  geom_line() +
  geom_point()

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3. Visualisation Species interactions ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.0 Load data ----

RawData <- list()
Parameters <- list()
country.list <- c("aus","spain")

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  for(Code.focal in Code.focal.list){ #focal.levels
    
    load(file =paste0(project.dic,"results/Parameters_",Code.focal,"_",country,".Rdata"))
    #assign(paste0("parameter_",Code.focal),
    #       parameter)
    Parameters[[paste(country,"_",Code.focal)]] <- parameter
  }
}

save(Parameters,
     file=paste0(home.dic,"results/Parameters_alpha.RData"))
load(paste0(home.dic,"results/Parameters_alpha.RData"))


#---- 3.1 Lambda ----
color.year <- data.frame(year=c("2015","2016","2017","2018","2019","2020","2021"),
                         col.value=colorblind_pal()(8)[2:8])
plot.lambda <- list()
env_pdsi_med_list <- list()
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  env_pdsi_med_list[[country]] <- read.csv(paste0(home.dic,"results/",
                                                  country,"_env_pdsi.csv")) %>%
    group_by(year) %>% 
    mutate(PDSI.mean = mean(get(paste0(country,"_pdsi")), na.rm = T),
           PDSI.sd = sd(get(paste0(country,"_pdsi")),na.rm = T)) %>%
    dplyr::select(year,PDSI.mean,PDSI.sd) %>%
    unique() %>%
    ungroup() %>%
    mutate_at("PDSI.mean",function(x) (x/ max(abs(x)))) %>%
    as.data.frame()
  
  for(Code.focal in Code.focal.list){ #focal.levels
    
    year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
    
    plot.lambda[[Code.focal]] <- Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd %>%
      mutate_at(year.levels, ~ rowSums(cbind(., Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean))) %>%
      tidyr::gather(all_of(year.levels), key="year", value="lambda_sd") %>%
      mutate(year=as.numeric(year)) %>%
      left_join(env_pdsi_med_list[[country]], by="year") %>%
      #ggplot(aes(y=lambda_sd, x=as.numeric(year),
      #          group=year)) +
      ggplot(aes(y=lambda_sd*(25/15), x=as.numeric(PDSI.mean),
                 group=year)) +
      geom_boxplot(width=0.2) +
      geom_hline(yintercept=median(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean[,1])*(25/15)) + 
      labs(title = Code.focal, #title="Intrinsic growth rate of LEMA across years,\n as their mean PDSI",
           y= "intrinsic growth rate",
           x="PDSI",fill="year") +
      theme_bw() +
      scale_color_manual(values= color.year[which(color.year$year %in% year.levels),"col.value"]) +
      scale_fill_manual(values= color.year[which(color.year$year %in% year.levels),"col.value"]) +
      guides(fill=guide_legend(nrow = 1,
                               direction="horizontal",
                               byrow = TRUE,
                               title.hjust = 0.1),
             color="none") +
      #coord_cartesian(ylim=c(0,15)) +
      theme( legend.key.size = unit(1, 'cm'),
             legend.position = "bottom",
             strip.background = element_blank(),
             panel.background = element_rect(fill='transparent'), #transparent panel bg
             plot.background = element_rect(fill='transparent', color=NA),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_blank(),
             legend.text=element_text(size=8),
             legend.title=element_text(size=8),
             axis.text.x= element_text(size=8, angle=66, hjust=1),
             axis.text.y= element_text(size=8),
             axis.title.x= element_text(size=8),
             axis.title.y= element_text(size=8),
             title=element_text(size=10))
    
  }
}
#figures/PEAI_lambda.pdf
#plot.lambda[13:24]
plot.lambda.all <- ggarrange(plotlist=plot.lambda[1:12],
                             common.legend = T,
                             legend = "bottom")
ggsave(plot.lambda.all,
       file=paste0(home.dic,"figures/plot.lambda.aus.pdf"))

plot.lambda.all <- ggarrange(plotlist=plot.lambda[13:24],
                             common.legend = T,
                             legend = "bottom")

ggsave(plot.lambda.all,
       file=paste0(home.dic,"figures/plot.lambda.spain.pdf"))



#---- 3.3. Sigmoid representation ----
source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL

Param.sigm.df <- list()
country.list = c("aus","spain")
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  sigmoid.list.focal<- list()
  Param.sigm.country.df <- NULL
  for(Code.focal in Code.focal.list){ #focal.levels
    
    
    df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
    
    Sp.names = colnames(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
    year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
    print(paste0(country,Code.focal))
    
    test.sigmoid.all<- NULL
    test.sigmoid  <- NULL
    
    for( neigh in Sp.names){
      # print(country)
      alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                             neigh]
      
      alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                           neigh]
      
      alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                       neigh]
      
      param.neigh <- data.frame(neigh = neigh, 
                                country = country,
                                alpha_initial = alpha_initial,
                                alpha_slope = alpha_slope,
                                alpha_c=  alpha_c,
                                N_opt_mean = Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh],
                                focal=Code.focal)
      
      for (n in 1:nrow(param.neigh)){
        #if(n==1){print(n)}
        df_neigh_n <- data.frame(density=c(0:10),param.neigh[n,])
        
        
        df_neigh_n[,"sigmoid"] <- alpha_function4(df_neigh_n$alpha_initial,
                                                  df_neigh_n$alpha_slope,
                                                  df_neigh_n$alpha_c,
                                                  df_neigh_n$density,
                                                  df_neigh_n$N_opt_mean)
        
        test.sigmoid <- bind_rows(test.sigmoid,df_neigh_n)
        
      }
    }
    limits.y <- c(round(min(test.sigmoid$sigmoid),digits=1),
                  round(max(test.sigmoid$sigmoid),digits=1))
    
    sigmoid.list.focal[[Code.focal]] <- ggplot(test.sigmoid,
                                               aes(y=sigmoid,x=density)) + 
      stat_smooth(color = "black", size = 0.5, level = 0.999) +
      theme_bw() + 
      scale_x_continuous(breaks = c(0,5,10))+
      scale_y_continuous(limits=limits.y) +
      geom_hline(yintercept=0, color="black") +
      labs(y="",fill="",color="",
           #title=Code.neigh,
           x=paste0("")) +
      facet_wrap(.~neigh, nrow=1) +
      guides(color=guide_legend(nrow = 1,
                                direction="horizontal",
                                byrow = TRUE,
                                title.hjust = 0.1),
             fill=guide_legend(nrow = 1,
                               direction="horizontal",
                               byrow = TRUE,
                               title.hjust = 0.1)) + 
      theme( legend.key.size = unit(1, 'cm'),
             legend.position = "bottom",
             strip.background = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_blank(),
             strip.text = element_text(size=12),
             legend.text=element_text(size=12),
             legend.title=element_text(size=12),
             #axis.ticks.x=element_blank(),
             axis.text.x= element_blank(),#element_text(size=12, angle=66, hjust=1),
             axis.text.y=  element_text(size=12),
             axis.title.x= element_blank(),#element_text(size=12),
             axis.title.y= element_blank(),#element_text(size=12),
             title=element_text(size=12))
    
    write.csv(test.sigmoid,
              file=paste0("results/Param.sigmoid.",country,",",Code.focal,".csv"))
    
    Param.sigm.country.df <- bind_rows( Param.sigm.country.df,test.sigmoid)
  }
  
  #sigmoid.list[[country]] <- sigmoid.list.focal
  Param.sigm.df[[country]] <- Param.sigm.country.df
  write.csv(Param.sigm.country.df,
            file=paste0(project.dic,"results/Param.sigmoid.",country,".csv.gz"))
}


write.csv(Param.sigm.df$aus,
          file=paste0(project.dic,"results/Param.sigmoid.aus.csv.gz"))
write.csv(Param.sigm.df$spain,
          file=paste0(project.dic,"results/Param.sigmoid.spain.csv.gz"))


Param.sigm.df.aus <- read.csv(paste0(project.dic,"results/Param.sigmoid.aus.csv.gz"))
Param.sigm.df.spain <- read.csv(paste0(project.dic,"results/Param.sigmoid.spain.csv.gz"))
Param.sigm.df$spain <- Param.sigm.df.spain
Param.sigm.df$aus<- Param.sigm.df.aus

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


#---- 3.3.1 Raw sigmoid illustration ----

sigmoid.list <- list()
legend.plot.list <- list()
for(country in "aus"){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  sigmoid.list.focal<- list()
  Param.sigm.country.df <- NULL
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  legend.plot <- ggplot(col.df, aes(y=neigh, x=color.name,color=neigh,fill=neigh)) +
    geom_bar(stat="identity") +
    scale_color_manual(values =col.df$color.name)+
    scale_fill_manual(values =col.df$color.name)+
    theme_bw() + theme( legend.key.size = unit(1, 'cm'),
                        legend.position = "bottom") +
    guides(color=guide_legend(nrow = 2,
                              direction="horizontal",
                              byrow = TRUE,
                              title.hjust = 0.1),
           fill=guide_legend(nrow = 2,
                             direction="horizontal",
                             byrow = TRUE,
                             title.hjust = 0.1)) 
  legend.plot.list[[country]] <- get_legend(legend.plot)
  
  for(Code.focal in c("TRCY","WAAC","POLE")){ #Code.focal.list
    print(paste(Code.focal))
    test.sigmoid <- Param.sigm.df[[country]] %>%
      filter(focal ==Code.focal)
    
    limits.y <- c(max(round(min(test.sigmoid$sigmoid),digits=1),-1),
                  min(round(max(test.sigmoid$sigmoid),digits=1),1))
    
    sigmoid.list.focal[[Code.focal]] <-  ggplot(test.sigmoid,
                                                aes(y=sigmoid,x=density,
                                                    color=neigh,fill=neigh)) + 
      stat_smooth( size = 1.5, level = 0.999) +
      theme_bw() + 
      scale_x_continuous(breaks = c(0,5,10))+
      scale_y_continuous(breaks=c(seq(-0.2,0.2,0.2)),limits=c(-0.3,0.15)) +
      scale_color_manual(values =col.df$color.name[which(col.df$neigh %in% test.sigmoid$neigh)])+
      scale_fill_manual(values =col.df$color.name[which(col.df$neigh %in% test.sigmoid$neigh)])+
      geom_hline(yintercept=0, color="black") +
      labs(y="",fill="",color="",
           title=Code.focal,
           x=paste0("")) +
      guides(color=guide_legend(nrow = 1,
                                direction="horizontal",
                                byrow = TRUE,
                                title.hjust = 0.1),
             fill=guide_legend(nrow = 1,
                               direction="horizontal",
                               byrow = TRUE,
                               title.hjust = 0.1)) + 
      #coord_cartesian(ylim=limits.y) + 
      theme( legend.key.size = unit(1, 'cm'),
             legend.position = "bottom",
             strip.background = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_blank(),
             strip.text = element_text(size=12),
             legend.text=element_text(size=12),
             legend.title=element_text(size=12),
             #axis.ticks.x=element_blank(),
             axis.text.x= element_text(size=12),#element_text(size=12, angle=66, hjust=1),
             axis.text.y=  element_text(size=12),
             axis.title.x= element_blank(),#element_text(size=12),
             axis.title.y= element_blank(),#element_text(size=12),
             title=element_text(size=12,color=col.df$color.name[which(col.df$neigh ==Code.focal)]))
    
  }
  sigmoid.list[[country]] <-  sigmoid.list.focal
}

ggarrange(sigmoid.list[["aus"]]$POLE,sigmoid.list[["aus"]]$TRCY,sigmoid.list[["aus"]]$WAAC,
          #labels=species.spain,
          common.legend = F,
          legend="none",nrow=1) #preso_sigmoid.pdf


AUS.sigmoid<- ggarrange(ggarrange(plotlist =   sigmoid.list[["aus"]],
                                  #labels=species.spain,
                                  common.legend = F,
                                  legend="none"),
                        legend.plot.list[["aus"]],
                        nrow=2,
                        heights=c(2,0.15),
                        common.legend = F)

ggsave(AUS.sigmoid,
       heigh=25,
       width=30,
       units = "cm",
       file=paste0(home.dic,"figures/AUS.sigmoid.pdf"))

SPAIN.sigmoid<- ggarrange(ggarrange(plotlist =   sigmoid.list[["spain"]],
                                    #labels=species.spain,
                                    common.legend = F,
                                    legend="none"),
                          legend.plot.list[["spain"]],
                          nrow=2,
                          heights=c(2,0.15),
                          common.legend = F)

ggsave(SPAIN.sigmoid,
       heigh=25,
       width=30,
       units = "cm",
       file=paste0(home.dic,"figures/SPAIN.sigmoid.pdf"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
