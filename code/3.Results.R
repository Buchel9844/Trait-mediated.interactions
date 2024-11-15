
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import packages----
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
#setwd("/home/lbuche/Eco_Bayesian/chapt3")
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Visualisation of Species Abundance ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2.0 Clean data----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))
country.list <- c("aus","spain")
#---- 2.1 Abundances over time ----
widthplot = 25
abundance_plotlist <- NULL
for(country in country.list){
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  if(country=="spain"){
    abundance_summary <- abundance_summary %>%
      gather( all_of(Code.focal.list), key="species",value="count")
  }
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  abundance_plotlist[[country]] <- ggplot() +
    stat_summary(data=abundance_summary,
                 aes(x=as.character(year), y = count*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=2) +
    stat_summary(data=abundance_summary,
                 aes(x=as.character(year), y = count*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 geom = "line",size=1) +
    scale_y_log10() +
    scale_x_discrete("year",limits= year.levels) +
    scale_color_manual(values= col.df$color.name) +
    labs(color="species",y="Mean number of \nindividuals in 25x25cm plot",
         title=paste0("Density over time of annual plants in ",country)) +
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
abundance_plotlist[["spain"]] # figures/Abundance_spain.pdf
abundance_plotlist[["aus"]] #figures/Abundance_aus.pdf

#---- 2.2 Abundances over PDSI ----
widthplot = 25
abundance_pdsi_plotlist <- NULL
for(country in country.list){
  
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  if(country=="spain"){
    abundance_summary <- abundance_summary %>%
      gather( all_of(Code.focal.list), key="species",value="count")
  }
  
  
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
    left_join(  env_pdsi ) %>%
    ggplot() +
    stat_summary(aes(x=exp(PDSI.mean), y = count*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange",size=2) +
    stat_summary(aes(x=exp(PDSI.mean), y = count*widthplot,
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
abundance_pdsi_plotlist[["spain"]] #figures/Abundance_pdsi_spain.pdf
abundance_pdsi_plotlist[["aus"]] #figures/Abundance_pdsi_aus.pdf

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. Visualisation Species interactions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
str(plot.lambda)
plot.lambda[["HOMA"]]
plot.lambda[["LEMA"]]
plot.lambda[["GORO"]]
plot.lambda[["WAAC"]]
plot.lambda[["PEAI"]]
plot.lambda[["ARCA"]]
plot.lambda[["TROR"]]
plot.lambda[["TRCY"]]
#figures/PEAI_lambda.pdf
length(Code.focal.list )
length(plot.lambda)
plot.lambda.all <- ggarrange(plotlist=plot.lambda[1:12],
                             common.legend = T,
                             legend = "bottom")
plot.lambda.all
ggsave(plot.lambda.all,
       file=paste0(home.dic,"figure/plot.lambda.pdf"))



#---- 3.3. Sigmoid representation ----
source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL

Param.sigm.df <- list()
country.list = c("aus","spain")
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  sigmoid.plot.list.focal<- list()
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
    
    sigmoid.plot.list.focal[[Code.focal]] <- ggplot(test.sigmoid,
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
  
  #sigmoid.plot.list[[country]] <- sigmoid.plot.list.focal
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

sigmoid.plot.list <- list()
legend.plot.list <- list()
country.list <- "aus"
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  sigmoid.plot.list.focal<- list()
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
  
  for(Code.focal in Code.focal.list){ #focal.levels
    print(paste(Code.focal))
    test.sigmoid <- Param.sigm.df[[country]] %>%
      filter(focal ==Code.focal)
    
    limits.y <- c(max(round(min(test.sigmoid$sigmoid),digits=1),-1),
                  min(round(max(test.sigmoid$sigmoid),digits=1),1))
    
    sigmoid.plot.list.focal[[Code.focal]] <-  ggplot(test.sigmoid,
                                                     aes(y=sigmoid,x=density,
                                                         color=neigh,fill=neigh)) + 
      stat_smooth( size = 1.5, level = 0.999) +
      theme_bw() + 
      scale_x_continuous(breaks = c(0,5,10))+
      #scale_y_continuous(breaks=c(seq(-1,1,0.4))) +
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
             axis.text.x= element_text(size=9),#element_text(size=12, angle=66, hjust=1),
             axis.text.y=  element_text(size=12),
             axis.title.x= element_blank(),#element_text(size=12),
             axis.title.y= element_blank(),#element_text(size=12),
             title=element_text(size=12,color=col.df$color.name[which(col.df$neigh ==Code.focal)]))
    
  }
  sigmoid.plot.list[[country]] <-  sigmoid.plot.list.focal
}


AUS.sigmoid.plot <- ggarrange(ggarrange(plotlist =   sigmoid.plot.list[["aus"]],
                                        #labels=species.spain,
                                        common.legend = F,
                                        legend="none"),
                              legend.plot.list[["aus"]],
                              nrow=2,
                              heights=c(2,0.15),
                              common.legend = F)

ggsave(AUS.sigmoid.plot,
       heigh=25,
       width=30,
       units = "cm",
       file=paste0(home.dic,"figures/AUS.sigmoid.plot.pdf"))

SPAIN.sigmoid.plot <- ggarrange(ggarrange(plotlist =   sigmoid.plot.list[["spain"]],
                                          #labels=species.spain,
                                          common.legend = F,
                                          legend="none"),
                                legend.plot.list[["spain"]],
                                nrow=2,
                                heights=c(2,0.15),
                                common.legend = F)

ggsave(SPAIN.sigmoid.plot,
       heigh=25,
       width=30,
       units = "cm",
       file=paste0(home.dic,"figures/SPAIN.sigmoid.plot.pdf"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. Realised interactions  ----

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL

Realised.Int.list <- list()

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  if(country  =="aus"){
    abundance_df <-  abundance_df %>%
      aggregate(count ~ year + species + id.plot +collector + scale.width, sum) %>%
      spread(species, count)  %>%
      dplyr::select(all_of(c("year",Code.focal.list))) 
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
      # print(country)
      neigh.abundance <- abundance_short_focal_df %>%
        dplyr::select(neigh) %>%
        dplyr::filter(!is.na(get(neigh))) 
      
      seq.abun.neigh <- seq(min(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)),
                            max(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)),
                            (max(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T))-min(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)))/10)
      
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
              file=paste0("results/Realised.Int.",country,",",Code.focal,".csv"))
    
    Realised.Int.country.df <- bind_rows( Realised.Int.country.df,test.sigmoid)
    
    save(Realised.Int.country.df,
         file=paste0(proejct.dic,"results/Realised.Int.df_",country,".RData"))
    
  }
  
  #sigmoid.plot.list[[country]] <- sigmoid.plot.list.focal
  Realised.Int.list[[country]] <-  Realised.Int.country.df
}

save(Realised.Int.list,
     file=paste0(project.dic,"results/Realised.Int.list.RData"))

load(paste0(project.dic,"results/Realised.Int.list.RData"))


#---- 3.1. Visualisation for SUPP ----

for(country in country.list){
  Box.Plot.Realised.effect <- NULL
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  Box.Plot.Realised.effect <- Realised.Int.list[[country]] %>%
    mutate(intra.bi = case_when(focal==neigh ~ "INTRA",
                                T~"INTER")) %>%
    ggplot(aes(y=realised.effect,x=as.factor(neigh),
               color=as.factor(neigh),fill=intra.bi)) + 
    geom_hline(yintercept = 1, color="black",size=1) +
    geom_boxplot() +
    facet_wrap(.~focal,ncol=3,scale="free") +
    labs(y="Effect on intrinsic performance",
         x= "Interacting species",
         color="",
         fill="") +
    #coord_cartesian(ylim=c(0.75,1.25),
    #               expand = FALSE) + 
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
         width=20,
         units = "cm",
         file=paste0(home.dic,"figures/","Boxplot.Realised.effect_",country,".pdf"))
}


#---- 3.3. Network visualasition----
library(igraph)
library(statnet)
library(intergraph)
library(ggraph)
colours  <- wes_palette("Zissou1", 101, type = "continuous")

plot.legend <- ggplot(data.frame(x=rep(c(1,2,3,4),
                                       times=50),
                                 y=1:100),
                      aes(as.factor(x),
                          y,fill=x)) +
  geom_point() +
  geom_line(aes(linewidth=as.factor(x))) + 
  scale_linewidth_discrete("Median effect on intrinsic performance",
                           labels=c("< 1%","1%-5%",
                                    "5%-10%",">10%")) +
  scale_fill_gradientn("Ratio of facilitation and competitive effect",
                       colours = colours,
                       breaks=c(1,51,101),
                       labels=c("100% \nFacilitative",
                                "50/50",
                                "100% \nCompetitive"),
                       limits=c(1,101)) +
  guides(fill= guide_colourbar(title.position="top", title.hjust = 0.5),
         linewidth = guide_legend(title.position="top", 
                                  title.hjust = 0.5,
                                  nrow=2)) +
  theme_bw() +
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.position = "bottom")
plot.legend 

net.country <- list()
# country = "aus"
species.focal = "GORO"
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  ratio.mat <-   Realised.Int.list[[country]]  %>%
    aggregate(realised.effect ~ focal + neigh, function(x) length(which(x<1))/length(x)) %>% # percentage of negative interaction
    #mutate(realised.effect = case_when(focal ==species.focal|neigh==species.focal ~ realised.effect,
    #                                 T~0)) %>%
    spread(neigh,realised.effect) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  
  strength.mat <- Realised.Int.list[[country]]  %>%
    aggregate(realised.effect ~ focal + neigh, function(x) abs(mean(x)-1)) %>%
    # mutate(realised.effect = case_when(focal ==species.focal|neigh==species.focal ~ realised.effect,
    #                                  T~0))%>%
    spread(neigh,realised.effect) %>%
    dplyr::select(-focal) %>%
    as.matrix() %>%
    t()# for the network to have row has emitor and columns as receiver
  
  alphamat.pos <- unname(abs(round(strength.mat,2))) %>%
    replace(is.na(.),0)
  g <- igraph::graph_from_adjacency_matrix(alphamat.pos  > 0)
  # Line width
  lwd.mat <-   as.numeric(round(strength.mat[alphamat.pos  > 0],2))
  E(g)$weight <- as.numeric(strength.mat[alphamat.pos  > 0])
  widths <- lwd.mat
  widths[lwd.mat <=0.01] <- 3 #1
  widths[lwd.mat >0.01 & lwd.mat <=0.05] <- 7 #3
  widths[lwd.mat >0.05 & lwd.mat <=0.1] <- 9#7
  widths[lwd.mat >0.1 ] <- 9
  
  width.arrow <- lwd.mat
  width.arrow[lwd.mat <=0.01] <- 2#0.3
  width.arrow[lwd.mat >0.01 & lwd.mat <=0.05] <- 2#1
  width.arrow[lwd.mat >0.05 & lwd.mat <=0.1] <- 4
  width.arrow[lwd.mat >0.1] <- 5
  
  #width.arrow[lwd.mat <=summary(lwd.mat)["1st Qu."]] <- 0.2
  #width.arrow[lwd.mat >summary(lwd.mat)["1st Qu."] & lwd.mat <=summary(lwd.mat)["Median"]] <- 0.5
  #width.arrow[lwd.mat >summary(lwd.mat)["Median"] & lwd.mat <=summary(lwd.mat)["3rd Qu."]] <- 1
  #width.arrow[lwd.mat > summary(lwd.mat)["3rd Qu."] & lwd.mat <= max(lwd.mat)] <- 1.5
  
  # to oriented the edge of the loop 
  edgeloopAngles <- numeric(0)
  b <- 1
  M <- dim(strength.mat)[1]
  m <- 0
  
  for(row in 1:nrow(alphamat.pos)) {
    for(col in 1:ncol(alphamat.pos)) {
      if (row == col) {
        m <- m+1
      }
      if (alphamat.pos[row,col] > 0) {
        edgeloopAngles[[b]] <- 0
        
        if (row == col) {
          edgeloopAngles[[b]] <- (3.2 * pi * (M - m) / M)
        }
        b <- b+1
      }
    }
  }
  
  # Colour edge
  col.vec <- round(as.numeric(ratio.mat[t(alphamat.pos)  > 0]),3) * 1000 + 1 
  colours  <- wes_palette("Zissou1", 1001, type = "continuous")
  #wes_palette("Zissou1")
  E(g)$color <- colours[col.vec] 
  #library(scales)
  #show_col(  colours[1] )
  #col.vec <- round(as.numeric(ratio.mat[t(alphamat.pos)  > 0]),3)
  #col.vec[col.vec>0.5] <- 1
  #col.vec[col.vec<=0.5] <- 2
  #colours  <- wes_palette("Zissou1", 1001, type = "continuous")
  #colours  <- c(colours[1000],colours[1])
  #E(g)$color <- colours[col.vec] 
  
  
  par(mar=c(0,0,0,0)+1)
  plot(g,
       layout=layout_in_circle, #layout_nicely, 
       edge.curved=-.3,
       margin=c(0.18,0,0.2,0),
       vertex.label = Code.focal.list,
       vertex.label.family="Helvetica",   
       vertex.size = 30,
       vertex.label.color = "black",
       vertex.color = "transparent",
       vertex.frame.color = "black",
       edge.width = widths,
       edge.arrow.width =  width.arrow,
       edge.loop.angle = edgeloopAngles)
  
  net.country[[country]] <- recordPlot()
  
  #ggsave( last_plot() ,
  #       heigh=40,
  #      units = "cm",
  #     file=paste0(home.dic,"figures/","Network.Realised.effect_",country,".pdf"))
}

plot_grid(net.country[[country]],
          get_legend(plot.legend),
          ncol = 1,
          rel_heights =c(1,0.2),
          labels = 'AUTO',
          hjust = 0, vjust = 1)

ggsave( plot_grid(net.country[[country]],
                  get_legend(plot.legend),
                  ncol = 1,
                  rel_heights =c(1,0.2),
                  labels = 'AUTO',
                  hjust = 0, vjust = 1),
        file=paste0(home.dic,"figures/","Network.Realised.effect_",country,".pdf")) 

# figures/Network.Realised.effect_aus.pdf")