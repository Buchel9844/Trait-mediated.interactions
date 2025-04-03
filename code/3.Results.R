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
library(Hmisc)

library(Polychrome)
library(knitr)
write_bib(x="stats") # creates citation bib
write_bib(x="PerformanceAnalytics")
write_bib(x="brms")
#setwd("/home/lbuche/Eco_Bayesian/chapt3")
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
home.dic <- ""
project.dic <- ""
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1. Visualisation of Species Abundance ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1.0 Clean data----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))
country.list <- c("aus","spain")

#-----1.1. Corrected abundance df ----
country="aus"
area.plot = pi*7.5^2
for(country in country.list){
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
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
    as.data.frame() %>%
    mutate(year=as.numeric(year))
  
  #abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".preclean")]]
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] %>%
    mutate(year=as.numeric(year)) %>%
    left_join(env_pdsi) %>%
    mutate(year=factor(year,levels=year.levels))  %>%
    dplyr::filter(!is.na(species)) %>%
    dplyr::filter(individuals>=0) %>%
    mutate(individuals.at.plot=individuals * area.plot)
  #group_by(year)
  ggplot( abundance_summary, aes(x=individuals.at.plot)) + geom_density()+
    facet_wrap(year~.,scale="free")
  
  # the coefficient if the sampling  effort differences between the years
  Correction.model.year <-  as.data.frame(confint(glmmTMB(individuals.at.plot ~ year + (1|species), 
                                                          data=abundance_summary,
                                                          family=nbinom2))) %>%
    as.data.frame() %>%
    mutate(year=c(year.levels,"random.factor")) %>%
    rownames_to_column("to.deleted") %>%
    dplyr::rename("Coef.year"="Estimate") %>%
    dplyr::filter(!year=="random.factor") %>%
    dplyr::select(year,Coef.year) %>%
    mutate(year=as.numeric(year)) 
  Correction.model.year$Coef.year[2:length(year.levels)] <- Correction.model.year$Coef.year[2:length(year.levels)] + Correction.model.year$Coef.year[1]
  
  
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] %>%
    mutate(year=as.numeric(year))
  
  abundance_summarycorrected <- abundance_summary %>%
    mutate(year=as.numeric(year)) %>%
    left_join(Correction.model.year) %>%
    mutate(individuals.at.plot =individuals * widthplot,
           corrected.density.at.plot = case_when(individuals > 0 ~ (individuals.at.plot-Coef.year),
                                         T ~ 0)) %>%
    mutate(corrected.density = corrected.density.at.plot/widthplot) %>%
    left_join(env_pdsi)
  write.csv(abundance_summarycorrected,
            file=paste0("results/Abundance.corrected.",country,".csv"))
  
  plot.abundance.density <- data.frame(
    abundance=c(abundance_summarycorrected$corrected.density,
                abundance_summary$individuals),
    year=c(abundance_summarycorrected$year,
           abundance_summary$year),
    nameab = c(rep("corrected", each=nrow(abundance_summarycorrected)),
               rep("observed", each=nrow(abundance_summary)))
  ) %>%
    ggplot( aes(x=abundance,fill=nameab,y=nameab)) +
    geom_density_ridges2(scale=1) +
    facet_wrap(year~.,scale="free")+
    theme_bw()
  plot.abundance.density
}
#---- 1.2 Abundances over time ----
widthplot = 25*25
abundance_plotlist <- NULL

country ="aus"    
richness_df <- NULL
for(country in country.list){
  #abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".preclean")]]
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  abundance_df  <- read.csv(paste0("results/Abundance.corrected.",country,".csv"))
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  df.split <- with(abundance_summary %>% 
               mutate(individuals.plot = individuals*widthplot), 
             split(abundance_summary %>% 
                     mutate(individuals.plot = individuals*widthplot),
                   list(year,species)))
  df.bb <- lapply( df.split, function(x) smean.cl.boot(x$individuals.plot,
                                                       conf.int=0.90, 
                                                       B=100,
                                                       na.rm=TRUE, reps=T))
  abundance_summary.for.ribbon <- do.call(rbind,df.bb) %>%
    unlist() %>%
    as.data.frame() %>%
    rownames_to_column("var.to.sep") %>%
    tidyr::separate(col="var.to.sep",into=c("year","species"),
                     extra = "merge") 
   

  richness_df[[country]] <- abundance_summary %>% 
    mutate(individuals.plot = individuals*widthplot) %>%
    filter(individuals.plot>0) %>%
    group_by(year, com_id) %>%
    summarise(n.richness = nlevels(as.factor(species))) %>%
    ungroup() %>%
    group_by(year) %>%
    summarise(richness.mean = mean(n.richness),
              richness.sd = sd(n.richness),
              richness.min = quantile(n.richness,0),
              richness.max = quantile(n.richness,1)) 
  
  #view(abundance_summary.for.ribbon)
  #view(abundance_summary %>%
   # aggregate(individuals ~ year + species, function(x) mean(x*widthplot)))
  abundance_plotlist[[country]] <- abundance_summary.for.ribbon %>%
    ggplot(aes(x=as.numeric(year), 
               y = Mean,
               group=as.factor(species),
               color=as.factor(species),
               fill=as.factor(species))) +
    geom_ribbon(aes(y=Mean, ymin=Lower,ymax=Upper),
                alpha=0.4) +
    #stat_summary(#fun.data = Hmisc::smean.cl.boot,
                 #fun.args=list(conf.int=0.95, na.rm=TRUE, reps=F),
                 #geom='ribbon',alpha=0.4) +
    scale_y_log10(breaks=c(.01,.1,1,10,100)) +
    scale_x_continuous() +
    scale_color_manual("species",values= col.df$color.name) +
    scale_fill_manual("species",values= col.df$color.name) +
    labs(y="Number of individuals \nin 25x25cm plot"
         ) +
    coord_cartesian( expand = F, default = FALSE, clip = "on") +
    theme_bw() +
    guides(fill = guide_legend(position="right",
                               ncol=1,override.aes = list(alpha=1)),
           color = guide_legend(position="right",ncol=1))+
    theme(legend.key.size = unit(0.5, 'cm'),
          #legend.position.inside =  c(1.08,0.5),
          legend.background = element_rect(fill="NA"),
          legend.box.background = element_rect(fill = alpha('grey98', .6)), 
          strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=28),
           legend.text=element_text(size=14),
           legend.title=element_blank(),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_blank(),
           axis.text.y= element_text(size=20),
           axis.title.x= element_blank(),
           axis.title.y= element_text(size=18),
           title=element_text(size=16),
           plot.margin=unit(c(1,1,0,0),"cm"))
abundance_plotlist[[country]]
}
abundance_plotlist[["spain"]] # figures/Abundance_spain.pdf
abundance_plotlist[["aus"]] #figures/Abundance_aus.pdf

write.csv(bind_rows(richness_df[["spain"]] %>%
          mutate(country="Spain"),
          richness_df[["aus"]]%>%
            mutate(country="Australia")),
          file="results/richness.df.csv")

#---- 1.3. Competition data ----
competition.plot.list <- list()
for(country in country.list){
Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
                     neigh = Code.focal.list)
if(country=="spain"){
competition_df <- get(paste0("clean.data.",country))[[paste0("competition_",country)]]  %>%
  gather(any_of(Code.focal.list),key="species",value="individuals") %>%
  mutate(individuals = case_when(focal==species ~ (individuals + 1),
                                 T~individuals)) %>%
  mutate(com_id = paste0(plot,"_", subplot)) %>%
  aggregate(individuals ~ year + com_id + species, sum) 
}else{
competition_df <- get(paste0("clean.data.",country))[[paste0("competition_",country)]]  %>%
  gather(any_of(Code.focal.list),key="species",value="individuals") %>%
  mutate(individuals = case_when(focal==species ~ (individuals + 1),
                                 T~individuals)) %>%
  mutate(individuals = (individuals/scale)*15) %>%
  dplyr::rename("com_id"="plot") %>%
  aggregate(individuals ~ year + com_id + species, sum) 
}

competition.plot  <- ggplot(competition_df,
                           aes(y=individuals, x= as.factor(year),
                               group=species, color=species)) +
  stat_summary( fun.y = median,
                fun.ymin = function(x) quantile(x,0.05), 
                fun.ymax = function(x) quantile(x,0.95), 
                geom = "pointrange",size=2,alpha=0.8,
                position=position_dodge2(width=0.5, reverse=T)) +
  scale_y_log10(breaks=c(.01,.1,1,10,100),
                limits=c(0.9,NA)) + 
  scale_color_manual("species",values= col.df$color.name) +
  scale_fill_manual("species",values= col.df$color.name) +
  labs(x="year",
       y="Median number of individuals \nin 15cm circle"
  ) +
  coord_cartesian( expand = F, default = FALSE, clip = "on") +
  theme_bw() +

  guides(fill = guide_legend(position="bottom",
                             nrow=1,override.aes = list(alpha=1)),
         color = guide_legend(position="bottom",nrow=1)) +
  theme(legend.key.size = unit(0.5, 'cm'),
        #legend.position.inside =  c(1.08,0.5),
        legend.background = element_rect(fill="NA"),
        legend.box.background = element_rect(fill = alpha('grey98', .6)), 
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=28),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.text.x= element_blank(),
        axis.text.y= element_text(size=20),
        #axis.title.x= element_blank(),
        axis.title.y= element_text(size=18),
        title=element_text(size=16),
        plot.margin=unit(c(1,1,1,1),"cm"))


fecundity.plot  <- get(paste0("clean.data.",country))[[paste0("competition_",country)]]%>%
  dplyr::filter(focal %in% Code.focal.list) %>%
  ggplot(aes(y=seed, x= as.factor(year),
                          group=focal, color=focal)) +
  stat_summary( fun.y = median,
                fun.ymin = function(x) quantile(x,0.05), 
                fun.ymax = function(x) quantile(x,0.95), 
                geom = "pointrange",size=2,alpha=0.8,
                position=position_dodge2(width=0.5, reverse=T)) +
  scale_y_log10(breaks=c(.01,.1,1,10,100,1000),
                limits=c(0.9,NA)) + 
  scale_color_manual("species",values= col.df$color.name) +
  scale_fill_manual("species",values= col.df$color.name) +
  labs(x="year",
       y="Median number of seeds per individuals"
  ) +
  coord_cartesian( expand = F, default = FALSE, clip = "on") +
  theme_bw() +
  guides(fill = guide_legend(position="bottom",
                             nrow=1,override.aes = list(alpha=1)),
         color = guide_legend(position="bottom",nrow=1)) +
  theme(legend.key.size = unit(0.5, 'cm'),
        #legend.position.inside =  c(1.08,0.5),
        legend.background = element_rect(fill="NA"),
        legend.box.background = element_rect(fill = alpha('grey98', .6)), 
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size=28),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.text.x= element_blank(),
        axis.text.y= element_text(size=20),
        #axis.title.x= element_blank(),
        axis.title.y= element_text(size=18),
        title=element_text(size=16),
        plot.margin=unit(c(1,1,1,1),"cm"))

competition.plot.list[[country]] <- ggarrange(competition.plot,
                                              fecundity.plot,ncol=2,
                                              labels=c("a. Abundances",
                                                       "b. Fecundity"),
                                              common.legend = T,
                                              legend = "bottom")
competition.plot.list[[country]]
}

competition.plot.list[["spain"]] # figures/supp/Spain_median_ind.pdf
competition.plot.list[["aus"]] #figures/supp/Aus_median_ind.pdf


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2. Visualisation Species interactions ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.0 Load data ----

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


#---- 2.1 Lambda ----
color.year <- data.frame(year=c("2015","2016","2017","2018","2019","2020","2021"),
                         col.value=colorblind_pal()(8)[2:8])
plot.lambda <- list()
env_pdsi_med_list <- list()
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  competition_df <- get(paste0("clean.data.",country))[[paste0("competition_",country)]] 
    
  
  for(Code.focal in Code.focal.list){ #focal.levels
    mean.fecundity <- competition_df %>% 
      dplyr::filter(focal ==Code.focal)
    year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
    
    plot.lambda[[Code.focal]] <- Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd %>%
      mutate_at(year.levels, ~ rowSums(cbind(., Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean))) %>%
      tidyr::gather(all_of(year.levels),
                    key="year", value="lambda_sd") %>%
      ggplot(aes(y=lambda_sd, 
                 x=as.factor(year))) +
      geom_boxplot(width=0.2) +
      geom_hline(yintercept=median(mean.fecundity$seed),
                 color="red") +
      geom_hline(yintercept=median(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean[,1])) + 
      labs(title = Code.focal, #title="Intrinsic growth rate of LEMA across years,\n as their mean PDSI",
           x= "year",
           y="Estimated intrinsic fecundity") +
      theme_bw() +
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
       file=paste0(home.dic,"figures/supp/plot.lambda.aus.pdf"))

plot.lambda.all <- ggarrange(plotlist=plot.lambda[13:24],
                             common.legend = T,
                             legend = "bottom")

ggsave(plot.lambda.all,
       file=paste0(home.dic,"figures/supp/plot.lambda.spain.pdf"))



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
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
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
