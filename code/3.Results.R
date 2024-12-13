#---- Exploration of Results ----

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
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 4. Realised interactions  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL

Realised.Int.list <- list()
widthplot = 15
# for median of individuals 
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  
  if(country  =="aus"){
    abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".preclean")]] 
    abundance_df <-  abundance_df %>%
      aggregate(individuals ~ year + species + com_id, sum) %>% 
      mutate(individuals = individuals *widthplot) %>%
      spread(species, individuals)  %>%
      dplyr::select(all_of(c("year",Code.focal.list))) 
  }else{
    abundance_df <-  abundance_df %>%
      mutate(individuals = individuals*widthplot) %>%
      spread(species, individuals)
  }
  Realised.Int.focal<- list()
  Realised.Int.country.df <- NULL
  for(Code.focal in Code.focal.list){ #focal.levels
    
    df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
    
    year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
    print(paste0(country,Code.focal))
    
    abundance_short_focal_df <-  abundance_df %>%
      filter(Code.focal > 0 | !is.na(Code.focal)) %>%
      mutate_all(as.numeric) 
    
    SpNames <- names(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
    
    test.sigmoid.all <- NULL
    test.sigmoid  <- NULL
    
    for( neigh in  SpNames){
      print(neigh)
      neigh.abundance <- abundance_short_focal_df %>%
        dplyr::filter(!is.na(get(neigh))) %>%
        dplyr::select(neigh) %>%
        unlist() %>%
        as.vector()
      
      seq.abun.neigh <- seq(min(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)),
                            max(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)),
                            (max(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), 
                                          na.rm = T))-min(quantile(neigh.abundance,probs=c(0.25,0.5,0.75),
                                                                   na.rm = T)))/10)[1:10]
      if(sum(seq.abun.neigh==0| is.na(seq.abun.neigh))>2){
        
        seq.abun.neigh <- c(0, mean(neigh.abundance)-sd(neigh.abundance),mean(neigh.abundance),
                            mean(neigh.abundance)+sd(neigh.abundance))
        
        seq.abun.neigh <-    seq.abun.neigh[which(   seq.abun.neigh>=0)]
      }
      library(HDInterval)
      
      alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                             neigh]
      
      alpha_initial  <-  alpha_initial[which(alpha_initial >= quantile(alpha_initial,probs=c(0.10)) &
                                               alpha_initial <=  quantile(alpha_initial,probs=c(0.9)))][1:6400]
      
      
      alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                           neigh]
      
      alpha_slope  <-  alpha_slope[which(alpha_slope>= quantile(alpha_slope,probs=c(0.10)) &
                                           alpha_slope <=  quantile(alpha_slope,probs=c(0.9)))][1:6400]
      
      
      alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                       neigh]
      
      alpha_c  <-   alpha_c[which( alpha_c >= quantile( alpha_c,probs=c(0.10)) &
                                     alpha_c <=  quantile( alpha_c,probs=c(0.9)))][1:6400]
      
      N_opt = Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh]
      
      N_opt  <-   N_opt[which(N_opt>= quantile( N_opt,probs=c(0.10)) &
                                N_opt <=  quantile( N_opt,probs=c(0.9)))][1:6400]
      
      param.neigh <- data.frame(neigh = neigh, 
                                country = country,
                                alpha_initial = alpha_initial,
                                alpha_slope = alpha_slope,
                                alpha_c=  alpha_c,
                                N_opt_mean = N_opt ,
                                focal=Code.focal)
      
      for (n in 1:nrow(param.neigh)){
        #if(n==1){print(n)}
        df_neigh_n <- data.frame(density=seq.abun.neigh,param.neigh[n,])
        
        
        df_neigh_n[,"sigmoid"] <- alpha_function4(df_neigh_n$alpha_initial,
                                                  df_neigh_n$alpha_slope,
                                                  df_neigh_n$alpha_c,
                                                  df_neigh_n$density,
                                                  df_neigh_n$N_opt_mean)
        
        df_neigh_n[,"realised.effect"] <- exp(df_neigh_n$sigmoid*df_neigh_n$density)
        
        test.sigmoid <- bind_rows(test.sigmoid,df_neigh_n)
        
      }
    }
    
    write.csv(test.sigmoid,
              file=paste0(project.dic,"results/Realised.Int.",country,".",Code.focal,".csv"))
    
    Realised.Int.country.df <- bind_rows( Realised.Int.country.df,test.sigmoid)
    
    save(Realised.Int.country.df,
         file=paste0(project.dic,"results/Realised.Int.df_",country,".RData"))
    
  }
  
  #sigmoid.plot.list[[country]] <- sigmoid.plot.list.focal
  Realised.Int.list[[country]] <-  Realised.Int.country.df
}

save(Realised.Int.list,
     file=paste0(project.dic,"results/Realised.Int.list.RData"))


load(paste0(project.dic,"results/Realised.Int.list.RData"))

# for actually 
test.sigmoid.all  <- NULL

Realised.Int.Obs.list <- list()
widthplot = 15
country ="aus"
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  
  if(country  =="aus"){
    #abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".clean")]] 
    abundance_df <-  abundance_df %>%
      aggregate(individuals ~ year + species + com_id, sum) %>% 
      mutate(individuals = individuals *widthplot) %>%
      spread(species, individuals) 
  }else{
    abundance_df <-  abundance_df %>%
      mutate(individuals = individuals*widthplot) %>%
      spread(species, individuals)
  }
  Realised.Int.Obs.focal<- list()
  Realised.Int.Obs.country.df <- NULL
  for(Code.focal in Code.focal.list){ #focal.levels
    
    df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
    
    year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
    print(paste0(country,Code.focal))
    
    abundance_short_focal_df <-  abundance_df %>%
      dplyr::filter(get(Code.focal) > 0) %>%
      dplyr::filter( !is.na(get(Code.focal))) %>%
      mutate_at(Code.focal,as.numeric) 
    head(abundance_short_focal_df)
    SpNames <- names(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
    
    test.sigmoid.all <- NULL
    test.sigmoid  <- NULL
    
    for( neigh in  SpNames){
      print(neigh)
      neigh.abundance <- abundance_short_focal_df %>%
        dplyr::filter(!is.na(get(neigh))) %>%
        dplyr::select(year, com_id,neigh) 
      names(neigh.abundance)[3] <-"density"
      
      alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                             neigh]
      
      alpha_initial  <-  quantile(alpha_initial,probs=c(0.5)) 
      
      
      alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                           neigh]
      
      alpha_slope  <- quantile(alpha_slope,probs=c(0.5)) 
      
      alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                       neigh]
      
      alpha_c  <-    quantile(alpha_c,probs=c(0.5)) 
      
      N_opt = Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh]
      
      N_opt  <-  quantile(N_opt,probs=c(0.5))
      
      param.neigh <- neigh.abundance %>%
        mutate(neigh = neigh, 
               country = country,
               alpha_initial = alpha_initial,
               alpha_slope = alpha_slope,
               alpha_c=  alpha_c,
               N_opt_mean = N_opt ,
               focal=Code.focal) 
      
      for (n in 1:nrow(param.neigh)){
        #if(n==1){print(n)}
        df_neigh_n <- param.neigh[n,]
        
        param.neigh[n,"sigmoid"] <- alpha_function4(df_neigh_n$alpha_initial,
                                                    df_neigh_n$alpha_slope,
                                                    df_neigh_n$alpha_c,
                                                    df_neigh_n$density,
                                                    df_neigh_n$N_opt_mean)
        
        param.neigh[n,"realised.effect"] <- exp(param.neigh[n,"sigmoid"]*df_neigh_n$density)
        
        
      }
      test.sigmoid <- bind_rows(test.sigmoid, param.neigh)
      
    }
    
    write.csv(test.sigmoid,
              file=paste0(project.dic,"results/Realised.Int.Obs.",country,".",Code.focal,".csv"))
    
    Realised.Int.Obs.country.df <- bind_rows( Realised.Int.Obs.country.df,test.sigmoid)
    
    save(Realised.Int.Obs.country.df,
         file=paste0(project.dic,"results/Realised.Int.Obs.df_",country,".RData"))
    
  }
  
  #sigmoid.plot.list[[country]] <- sigmoid.plot.list.focal
  Realised.Int.Obs.list[[country]] <-  Realised.Int.Obs.country.df
}

save(Realised.Int.Obs.list,
     file=paste0(project.dic,"results/Realised.Int.Obs.list.RData"))


load(paste0(project.dic,"results/Realised.Int.Obs.list.RData"))



#---- 4.1. Visualisation for SUPP ----

for(country in country.list){
  Box.Plot.Realised.effect <- NULL
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  Box.Plot.Realised.effect <- Realised.Int.Obs.list[[country]] %>%
    mutate(intra.bi = case_when(focal==neigh ~ "INTRA",
                                T~"INTER")) %>%
    ggplot(aes(y=realised.effect,x=as.factor(neigh),
               color=as.factor(neigh),fill=intra.bi)) + 
    geom_hline(yintercept = 1, color="black",size=1) +
    geom_boxplot() +
    facet_wrap(.~focal,ncol=4,scale="free") +
    labs(y="Effect on intrinsic performance",
         x= "Interacting species",
         color="",
         fill="") +
    coord_cartesian(ylim=c(0,2),
                    expand = FALSE) + 
    theme_bw() +
    scale_color_manual(values =col.df$color.name)+
    scale_fill_manual(values =c("white","black"))+
    theme_bw() + theme( legend.key.size = unit(1, 'cm'),
                        legend.position = "bottom") +
    guides(color=guide_legend(nrow = 3,
                              direction="horizontal",
                              byrow = TRUE,
                              title.hjust = 0.1),
           fill=guide_legend(nrow = 2,
                             direction="horizontal",
                             byrow = TRUE,
                             title.hjust = 0.1)) +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           legend.title = element_text(size=12,color="black"),
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=12),
           legend.text=element_text(size=12),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_blank(),#element_text(size=12, angle=66, hjust=1),
           axis.text.y=  element_text(size=12),
           axis.title.x= element_blank(),
           axis.title.y= element_text(size=12))
  Box.Plot.Realised.effect
  ggsave(Box.Plot.Realised.effect,
         heigh=25,
         width=30,
         units = "cm",
         file=paste0(home.dic,"figures/","Boxplot.Obs.Realised.effect_",country,".pdf"))
}


#---- 4.3. Network visualasition----
colours  <- wes_palette("Zissou1", 101, type = "continuous")

plot.legend <- ggplot(data.frame(x=rep(c(1,2,3,4,5),
                                       times=50),
                                 y=1:250),
                      aes(as.factor(x),
                          y,fill=x)) +
  geom_point() +
  geom_line(aes(linewidth=as.factor(x), alpha=as.factor(x))) + 
  scale_linewidth_manual("Median effect on intrinsic performance",
                         values=c(0.5,1,1.5,2,3),
                         labels=c("< 5%","5%-10%","10%-25%","25%-50%",">50%")) +
  scale_alpha_manual("Median effect on intrinsic performance",
                     values=c(0.2,0.4,0.6,0.8,1),
                     labels=c("< 5%","5%-10%","10%-25%","25%-50%",">50%")) +
  scale_fill_gradientn("Ratio of facilitation and competitive effect",
                       colours = colours,
                       breaks=c(1,51,101),
                       labels=c("100% \nFacilitative",
                                "50/50",
                                "100% \nCompetitive"),
                       limits=c(1,101)) +
  guides(fill= guide_colourbar(title.position="top", title.hjust = 0.5),
         linewidth = guide_legend(title.position="top", 
                                  title.hjust = 0,
                                  nrow=2)) +
  theme_bw() +
  theme(legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        legend.text = element_text(size=16),
        legend.title = element_text(size=20))
plot.legend 

#load(paste0(project.dic,"results/Realised.Int.list.RData"))
net.country <- list()
# country = "spain"
#species.focal = "GORO"

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  ratio.mat <-   Realised.Int.Obs.list[[country]]  %>%
    aggregate(realised.effect ~ focal + neigh, function(x) length(which(x<1))/length(x)) %>% # percentage of negative interaction
    spread(neigh,realised.effect) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  
  strength.mat <- Realised.Int.Obs.list[[country]]  %>%
    aggregate(realised.effect ~ focal + neigh, function(x) abs(mean(x)-1)) %>%
    spread(neigh,realised.effect) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  # for the network to have row as receiver
  if(!dim(strength.mat)[1]==dim(strength.mat)[2]){
    df.all.comb <- expand.grid(Code.focal.list,Code.focal.list) %>%
      rename("focal"="Var1") %>%
      rename("neigh"="Var2")
    
    ratio.mat <-Realised.Int.Obs.list[[country]]  %>%
      aggregate(realised.effect ~ focal + neigh,
                function(x) length(which(x<1))/length(x)) %>% # percentage of negative interaction
      right_join(df.all.comb) %>%
      spread(neigh,realised.effect) %>%
      dplyr::select(-focal) %>%
      as.matrix()
    
    strength.mat <- Realised.Int.Obs.list[[country]]  %>%
      aggregate(realised.effect ~ focal + neigh, function(x) abs(mean(x)-1)) %>%
      right_join(df.all.comb) %>%
      spread(neigh,realised.effect) %>%
      dplyr::select(-focal) %>%
      as.matrix() %>%
      t()# for the network to have row has emitor and columns as receiver
  }
  plot.network.gradient.int(ratio.mat,strength.mat,"",0.01)
  net.country[[country]] <- recordPlot()
  
  
}
net.country$spain
# figures/Network.Realised.effect_spain.pdf
net.country$aus
# figures/Network.Realised.effect_aus.pdf
plot(get_legend(plot.legend))
#figures/Network.Realised.effect_legend.pdf
plot_grid(
  plot_grid(net.country$aus,net.country$spain,
            nrow = 1,labels=c("Australia","Spain")),
  get_legend(plot.legend),
  ncol = 1,
  rel_heights =c(1,0.2),
  labels = '')
# figures/Network.Realised.effect.pdf
ggsave( plot_grid(net.country[["aus"]],
                  get_legend(plot.legend),
                  ncol = 1,
                  rel_heights =c(1,0.24),
                  labels = "",
                  hjust = 0, vjust = 1),
        file=paste0(home.dic,"figures/Network.Obs.Realised.effect_aus.pdf")) 
ggsave( plot_grid(net.country[["spain"]],
                  get_legend(plot.legend),
                  ncol = 1,
                  rel_heights =c(1,0.24),
                  labels = "",
                  hjust = 0, vjust = 1),
        file=paste0(home.dic,"figures/Network.Obs.Realised.effect_spain.pdf")) 

# figures/Network.Realised.effect_aus.pdf")

library(qgraph)
plotwebfun<-function(web,metaweb=NULL,colmat,
                     layout = "circle", curveDefault = 0.4,
                     wtimes=200,vertex.size=8,edge.arrow.size=1,
                     minimum,theme,parallelAngleDefault,edge.width,
                     labels,curveAll){
  web<- t(web)
  if(!is.null(metaweb)){
    spcode<-sapply(colnames(web),
                   function(x)which(colnames(metaweb)==x))#could be useful
  }
  qgraph(web,color = "white", layout = layout,
         curveDefault = curveDefault, minimum=minimum,
         edge.width=edge.width,
         theme=theme,
         parallelAngleDefault=parallelAngleDefault,
         weighted=T,
         labels=labels,curveAll=curveAll)
}
plotwebfun( strength.mat, colmat=ratio.mat,
            minimum = 0,edge.width=2,
            parallelAngleDefault=pi/10,theme="Hollywood",
            labels= colnames(strength.mat),curveAll=F)

plot.network.gradient.int <- function(ratio.mat,strength.mat,
                                      title.y,minimum.stength){
  alphamat.pos <- unname(abs(round(strength.mat,3))) %>%
    replace(is.na(.),0)
  g <- igraph::graph_from_adjacency_matrix(t(alphamat.pos)  > 0,
                                           weighted=T, mode="directed")
  # Line width
  lwd.mat <-   as.numeric(round(strength.mat[alphamat.pos  > 0],3))
  E(g)$weight <- as.numeric(strength.mat[alphamat.pos  > 0])
  width.edge <- lwd.mat
  width.edge[lwd.mat <=0.05] <- 1 #1
  width.edge[lwd.mat >0.05 & lwd.mat <=0.1] <- 2 #3
  width.edge[lwd.mat >0.1 & lwd.mat <=0.25] <- 3#7
  width.edge[lwd.mat >0.25 & lwd.mat <=0.5] <- 5#7
  width.edge[lwd.mat >0.5 ] <- 7
  width.edge[lwd.mat <=0.01] <- 0.5 #1
  
  width.arrow <- lwd.mat
  width.arrow[lwd.mat <=0.05] <- 0.3#0.3
  width.arrow[lwd.mat >0.05 & lwd.mat <=0.1] <- 0.6#1
  width.arrow[lwd.mat >0.1 & lwd.mat <=0.25] <- 1
  width.arrow[lwd.mat >0.25 & lwd.mat <=0.5] <- 1.5
  width.arrow[lwd.mat >0.5] <- 2
  
  
  #width.arrow[lwd.mat <=summary(lwd.mat)["1st Qu."]] <- 0.2
  #width.arrow[lwd.mat >summary(lwd.mat)["1st Qu."] & lwd.mat <=summary(lwd.mat)["Median"]] <- 0.5
  #width.arrow[lwd.mat >summary(lwd.mat)["Median"] & lwd.mat <=summary(lwd.mat)["3rd Qu."]] <- 1
  #width.arrow[lwd.mat > summary(lwd.mat)["3rd Qu."] & lwd.mat <= max(lwd.mat)] <- 1.5
  
  # to oriented the edge of the loop 
  edgeloopAngles <- numeric(0)
  b <- 1
  M <- dim(strength.mat)[1]
  m <- 0
  for(col in 1:ncol(alphamat.pos)) {
    for(row in 1:nrow(alphamat.pos)) {
      if (row == col) {
        m <- m+1
      }
      if (alphamat.pos[row,col] > 0) {
        edgeloopAngles[[b]] <- 0
        
        if (row == col) {
          edgeloopAngles[[b]] <- (2.2 * pi * (M - m) / M)
        }
        b <- b+1
      }
    }
  }
  
  # Colour edge
  library(grDevices)
  col.vec <- round(as.numeric(ratio.mat[alphamat.pos  > 0]),3) * 1000 + 1 
  colours  <- wes_palette("Zissou1", 1001, type = "continuous")
  col.mat <-  colours[col.vec] 
  alpha.edge <-  lwd.mat
  alpha.edge[lwd.mat <=minimum.stength] <-0.2
  alpha.edge[lwd.mat >minimum.stength & lwd.mat <=2*minimum.stength] <- 0.4 #1
  alpha.edge[lwd.mat >2*minimum.stength & lwd.mat <=4*minimum.stength ] <- 0.6 #3
  alpha.edge[lwd.mat >4*minimum.stength& lwd.mat <=6*minimum.stength] <- 0.8 #3
  alpha.edge[lwd.mat >6*minimum.stength] <- 1
  
  col.mat[alpha.edge==0.2] <- adjustcolor( col.mat[alpha.edge==0.2],
                                           0.2)
  col.mat[which(alpha.edge==0.4)] <- adjustcolor( col.mat[which(alpha.edge==0.4)],
                                                  0.4)
  col.mat[alpha.edge==0.6] <- adjustcolor( col.mat[alpha.edge==0.6],
                                           0.6)
  col.mat[alpha.edge==0.8] <- adjustcolor( col.mat[alpha.edge==0.8],
                                           0.8)
  col.mat[alpha.edge==1] <- adjustcolor( col.mat[alpha.edge==1],
                                         1)
  E(g)$color <- col.mat
  
  # plot the network
  par(mar=c(0,3,0,0))
  plot(g,
       layout=layout_in_circle, #layout_nicely, 
       edge.curved=-.15,
       #margin=c(0.5,-0.8,0.1,-1),
       margin=c(0.2,0.2,0.2,0.2),
       vertex.label = colnames(strength.mat),
       vertex.label.family="Helvetica",   
       vertex.size = 30,
       vertex.label.color = "black",
       vertex.color = "transparent",
       vertex.frame.color = "black",
       edge.width = width.edge,
       edge.arrow.width =  1.5,#width.arrow,
       edge.loop.angle = edgeloopAngles)
  title(title.y, line = -2)
}

#---- 4.4. Table sum up ----
str(Realised.Int.list)
sum.up.df <- NULL
for(country in country.list){
  
  sum.up.df.n <- Realised.Int.list[[country]] %>%
    summarise(mean.effect = (mean(realised.effect)-1)*100,
              median.effect = (median(realised.effect)-1)*100,
              var.effect = (var(realised.effect))*100,
              max.positive.effect = (max(realised.effect)-1)*100,
              max.negative.effect = (min(realised.effect)-1)*100,
              count.positive = count(realised.effect >1),
              count.negative = count(realised.effect <1),
              count.total = count(realised.effect>0)) %>%
    mutate(proportion.positive =(count.positive/count.total) * 100,
           proportion.negative =(count.negative/count.total)* 100,
           proportion.neutre =100-(proportion.positive +proportion.negative),
           country = country,
           effect ="both",
           species = "All")
  
  
  sum.up.neigh.df.n <- Realised.Int.list[[country]] %>%
    group_by(neigh) %>%
    summarize(mean.effect = (mean(realised.effect)-1)*100,
              median.effect = (median(realised.effect)-1)*100,
              var.effect = (var(realised.effect))*100,
              max.positive.effect = (max(realised.effect)-1)*100,
              max.negative.effect = (min(realised.effect)-1)*100,
              count.neutre = count(realised.effect ==1),
              count.positive = count(realised.effect >1),
              count.negative = count(realised.effect <1),
              count.total = count(realised.effect>0)) %>%
    mutate(proportion.positive =(count.positive/count.total) * 100,
           proportion.negative =(count.negative/count.total)* 100,
           proportion.neutre =100-(proportion.positive +proportion.negative),
           effect ="given")%>%
    rename(species = "neigh")
  
  sum.up.focal.df.n <- Realised.Int.list[[country]] %>%
    group_by(focal) %>%
    summarize(mean.effect = (mean(realised.effect)-1)*100,
              median.effect = (median(realised.effect)-1)*100,
              var.effect = (var(realised.effect))*100,
              max.positive.effect = (max(realised.effect)-1)*100,
              max.negative.effect = (min(realised.effect)-1)*100,
              count.neutre = count(realised.effect ==1),
              count.positive = count(realised.effect >1),
              count.negative = count(realised.effect <1),
              count.total = count(realised.effect>0)) %>%
    mutate(proportion.positive =(count.positive/count.total) * 100,
           proportion.negative =(count.negative/count.total)* 100,
           proportion.neutre =100-(proportion.positive +proportion.negative),
           country = country,
           effect ="received") %>%
    rename(species = "focal")
  
  
  write.csv(bind_rows(sum.up.df.n,sum.up.focal.df.n,sum.up.neigh.df.n),
            file=paste0(home.dic,"results/Sum.up.species",country,".csv"))
  
}


sum.up.focal.df.n <- Realised.Int.list[[country]] %>%
  aggregate(realised.effect ~ focal, median)



#---- 4.5 With different grouping of functional group ----
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

#---- 4.6 Along the fast to slow gradient ----
Realised.Int.FSspectrum.plotlist <- list()

country = "spain"
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
    dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal) %>%
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
  #hist(Realised.Int.sum$realised.effect,breaks=150)
  Realised.Int.FSspectrum.plotlist[[country]] 
  #figures/Heatmap.N15_spain.pdf
  #figures/Heatmap.root.volume.less0.5_aus.pdf
  Realised.Int.sum %>%
    mutate(transparency = case_when(RE.mean < 0.05 ~ 0.2,
                                    (RE.mean >0.05 & RE.mean <=0.1) ~ 0.4,
                                    (RE.mean>0.1 &RE.mean <=0.25) ~ 0.6,
                                    (RE.mean >0.25 & RE.mean<=0.5) ~ 0.8,
                                    (RE.mean>0.5 ) ~ 1)) %>%
    mutate(neigh=factor(neigh, levels=c(trait.df$species[order(trait.df$TDMr)])),
           focal = factor(focal,levels=rev(trait.df$species[order(trait.df$TDMr)]))) %>%
    ggplot(aes(x =neigh,y =focal)) + 
    geom_tile(aes(fill = RE.ratio,alpha=abs(RE.mean-1))) + 
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))   +
    theme_few() +
    scale_alpha_continuous(range = c(0.5, 1))+
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
  
  
}
Realised.Int.FSspectrum.plotlist[["aus"]] #figures/Heatmap.FS.spectrum_aus.pdf

Realised.Int.FSspectrum.plotlist[["spain"]] #figures/Heatmap.FS.spectrum_spain.pdf

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 5.Network Metrics ---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 5.1 Heatmap Fac/Comp ----
Realised.Int.plotlist <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  Realised.Int.list.df <- Realised.Int.list[[country]]
  
  Realised.Int.sum <- Realised.Int.list.df %>%
    group_by(neigh,focal)%>%
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
  
  Realised.Int.focal <- Realised.Int.list.df %>%
    group_by(focal)%>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           neigh ="total") %>%
    as.data.frame() 
  
  Realised.Int.neigh <-Realised.Int.list.df%>%
    group_by(neigh)%>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           focal ="total") %>%
    as.data.frame() 
  
  Realised.Int.total <- Realised.Int.list.df %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           focal ="total",
           neigh="total") %>%
    as.data.frame() 
  
  
  Realised.Int.plotlist[[country]] <- Realised.Int.sum %>%
    bind_rows(Realised.Int.neigh,Realised.Int.focal, Realised.Int.total) %>%
    mutate(transparency = case_when(RE.mean < 0.05 ~ 0.2,
                                    (RE.mean >0.05 & RE.mean <=0.1) ~ 0.4,
                                    (RE.mean>0.1 &RE.mean <=0.25) ~ 0.6,
                                    (RE.mean >0.25 & RE.mean<=0.5) ~ 0.8,
                                    (RE.mean>0.5 ) ~ 1)) %>%
    mutate(neigh=factor(neigh, levels=c(Code.focal.list,"total")),
           focal = factor(focal,levels=c("total",rev(Code.focal.list)))) %>%
    ggplot(aes(x =neigh,y =focal)) + 
    geom_tile(aes(fill = RE.ratio,alpha=abs(RE.mean-1))) + 
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))   +
    annotate("rect",xmin = 0.5, xmax = 13.5, ymin = .5,
             ymax = 1.5, color="black",fill="transparent") + 
    annotate("rect",ymin = 0.5, ymax = 13.5, xmin = 13- .5,
             xmax = 13 + .5, color="black",fill="transparent") + 
    theme_few() +
    scale_alpha_continuous(range = c(0.5, 1))+
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
}
Realised.Int.plotlist[["aus"]] # figures/Heatmap.aus.pdf
Realised.Int.plotlist[["spain"]]  # figures/Heatmap.spain.pdf
#---- 5.2. Asymetrie of interactions ----
Realised.Int.plotlist <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  Realised.Int.list.df <- Realised.Int.list[[country]]
  
  Realised.Int.sum <- Realised.Int.list.df %>%
    group_by(neigh,focal)%>%
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
  
  strength.mat <- Realised.Int.sum %>%
    dplyr::select(neigh,focal,RE.Q5) %>%
    spread(neigh,RE.Q5) %>%
    column_to_rownames("focal")
  
  sym.mat <- symmetry.with.positive(strength.mat) %>%
    group_by(focal,symmetry) %>%
    #tally() %>%
    mutate(symmetry.names = case_when((symmetry =="++"|symmetry=="--")~"sym",
                                      (symmetry =="+-"|symmetry=="-+")~"asym"))
  
  ggplot( sym.mat, aes(x=symmetry)) + geom_bar()
  
}

symmetry.with.positive <- function(strength.mat){
  names.sp <- colnames(strength.mat)
  sym.mat <- data.frame(focal = rep(names.sp , each=length(names.sp )),
                        neigh = rep(names.sp , times=length(names.sp )))
  for(i in 1:ncol(strength.mat)){
    for(j in 1:ncol(strength.mat)){
      #print(paste0(i,j))
      if(i==j) next
      if(is.na(strength.mat[i,j])) next
      if(is.na(strength.mat[j,i])) next
      if(strength.mat[i,j] > 1 & strength.mat[j,i]>1){
        sym.mat$symmetry[which(sym.mat$focal==names.sp[i] &
                                 sym.mat$neigh==names.sp[i])] <- "++"
      }
      if(strength.mat[i,j] > 1 & strength.mat[j,i]<1){
        sym.mat$symmetry[which(sym.mat$focal==names.sp[i] &
                                 sym.mat$neigh==names.sp[i])] <- "+-"
      }
      if(strength.mat[i,j] < 1 & strength.mat[j,i]>1){
        sym.mat$symmetry[which(sym.mat$focal==names.sp[i] &
                                 sym.mat$neigh==names.sp[i])] <- "-+"
      }
      if(strength.mat[i,j] < 1 & strength.mat[j,i]<1){
        sym.mat$symmetry[which(sym.mat$focal==names.sp[i] &
                                 sym.mat$neigh==names.sp[i])] <- "--"
      }
    }
  }
  return(sym.mat)
}
