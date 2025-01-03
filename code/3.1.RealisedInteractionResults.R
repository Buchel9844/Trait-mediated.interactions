
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
#---- 1. Compute realised interactions  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL
#---- 1.1. Realised interactions across years  ----

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
      
      seq.abun.neigh <- seq(min(quantile(neigh.abundance,probs=c(0.10,0.5,0.80), na.rm = T)),
                            max(quantile(neigh.abundance,probs=c(0.10,0.5,0.80), na.rm = T)),
                            (max(quantile(neigh.abundance,probs=c(0.10,0.5,0.80), 
                                          na.rm = T))-min(quantile(neigh.abundance,probs=c(0.10,0.5,0.80),
                                                                   na.rm = T)))/10)[1:10]
      if(sum(seq.abun.neigh==0| is.na(seq.abun.neigh))>2){
        
        seq.abun.neigh <- c(0, mean(neigh.abundance)-sd(neigh.abundance),mean(neigh.abundance),
                            mean(neigh.abundance)+sd(neigh.abundance))
        
        seq.abun.neigh <-    seq.abun.neigh[which(   seq.abun.neigh>=0)]
      }
      library(HDInterval)
      
      alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                             neigh]
      
      #alpha_initial  <-  alpha_initial[which(alpha_initial >= quantile(alpha_initial,probs=c(0.10)) &
      #                                         alpha_initial <=  quantile(alpha_initial,probs=c(0.9)))][1:6400]
      alpha_initial  <-  quantile(alpha_initial,probs=c(0.50))
      
      alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                           neigh]
      
      #alpha_slope  <-  alpha_slope[which(alpha_slope>= quantile(alpha_slope,probs=c(0.10)) &
      #                                     alpha_slope <=  quantile(alpha_slope,probs=c(0.9)))][1:6400]
      
      alpha_slope  <- quantile(alpha_slope,probs=c(0.50))
      
      alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                       neigh]
      
      #alpha_c  <-   alpha_c[which( alpha_c >= quantile( alpha_c,probs=c(0.10)) &
       #                              alpha_c <=  quantile( alpha_c,probs=c(0.9)))][1:6400]
      alpha_c  <-quantile( alpha_c,probs=c(0.50)) 
      
      N_opt = Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh]
      
      #N_opt  <-   N_opt[which(N_opt>= quantile( N_opt,probs=c(0.10)) &
      #                          N_opt <=  quantile( N_opt,probs=c(0.9)))][1:6400]
      N_opt  <-quantile( N_opt,probs=c(0.50))
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
        
        df_neigh_n$sigmoid[which(df_neigh_n$density==0)] <- 0
        
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
#---- 1.2. Compute realised interactions based on time  ----

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL
Realised.Int.Year.list <- list()
widthplot = 15
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  if(country  =="aus"){
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
  Realised.Int.Year.focal<- list()
  Realised.Int.Year.country.df <- NULL
  for(Code.focal in Code.focal.list){ #focal.levels
    lambda = median(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean[,1], na.rm=T)
    
    df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
    
    print(paste0(country,Code.focal))
    
    abundance_short_focal_df <-  abundance_df %>%
      filter(Code.focal > 0 | !is.na(Code.focal)) %>%
      mutate_at(Code.focal.list,as.numeric) 
    year.levels <- levels(as.factor(abundance_short_focal_df$year))
    
    SpNames <- names(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
    
    test.sigmoid.all <- NULL
    test.sigmoid  <- NULL
    for( y in year.levels){
      for( neigh.sp in  SpNames){
        print(paste(neigh.sp,y))
        
        neigh.abundance <- abundance_short_focal_df %>%
          filter(year ==y) %>%
          dplyr::select(neigh.sp) %>%
          filter(!is.na(get(neigh.sp)))
        if(nrow(neigh.abundance)==0) next
        
        seq.abun.neigh <- seq(min(quantile(neigh.abundance,probs=c(0.10,0.5,0.80), na.rm = T)),
                              max(quantile(neigh.abundance,probs=c(0.10,0.5,0.80), na.rm = T)),
                              (max(quantile(neigh.abundance,probs=c(0.10,0.5,0.80), na.rm = T))-min(quantile(neigh.abundance,probs=c(0.10,0.5,0.80), na.rm = T)))/100)
        
        param.neigh <- data.frame(neigh = neigh.sp, 
                                  country = country,
                                  alpha_initial = median(df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                                                                neigh.sp]),
                                  alpha_slope =  median(df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                                                               neigh.sp]),
                                  alpha_c=  median(df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                                                          neigh.sp]),
                                  N_opt_mean =median(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh.sp]),
                                  focal=Code.focal)
        
        for (n in 1:nrow(param.neigh)){
          #if(n==1){print(n)}
          df_neigh_n <- data.frame(density=seq.abun.neigh,param.neigh[n,])
          
          
          df_neigh_n[,"sigmoid"] <- alpha_function4(df_neigh_n$alpha_initial,
                                                    df_neigh_n$alpha_slope,
                                                    df_neigh_n$alpha_c,
                                                    df_neigh_n$density,
                                                    df_neigh_n$N_opt_mean)
          df_neigh_n$sigmoid[which(df_neigh_n$density==0)] <- 0
          
          df_neigh_n[,"sigmoid.std"] <- (df_neigh_n$sigmoid)/log(lambda)
          
          df_neigh_n[,"realised.effect"] <- exp(df_neigh_n$sigmoid*df_neigh_n$density)
          df_neigh_n[,"realised.effect.std"] <- (df_neigh_n$sigmoid*df_neigh_n$density)/log(lambda)
          
          df_neigh_n$year <- as.numeric(y)
          
          test.sigmoid <- bind_rows(test.sigmoid,df_neigh_n)
          
        }
      }
    }
    
    write.csv(test.sigmoid,
              file=paste0("results/Realised.Int.Year",country,",",Code.focal,".csv"))
    
    Realised.Int.Year.country.df <- bind_rows(Realised.Int.Year.country.df,test.sigmoid)
    
    save(Realised.Int.Year.country.df,
         file=paste0(project.dic,"results/Realised.Int.Year.df_",country,".RData"))
    
  }
  
  #sigmoid.list[[country]] <- sigmoid.list.focal
  Realised.Int.Year.list[[country]] <-  Realised.Int.Year.country.df
}

save(Realised.Int.Year.list,
     file=paste0(project.dic,
                 "results/Realised.Int.Year.list.RData"))
load(paste0(project.dic,
            "results/Realised.Int.Year.list.RData"))


#---- 1.3. Realised interactions for obs  ----

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
  # Code.focal ="ARCA"
  for(Code.focal in Code.focal.list){ #focal.levels
    lambda = median(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean[,1], na.rm=T)
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
        if(df_neigh_n$density==0){
          param.neigh[n,"sigmoid"] <- 0
        }
        
        param.neigh[n,"realised.effect"] <- exp(param.neigh[n,"sigmoid"]*df_neigh_n$density)
        param.neigh[n,"realised.effect.std.lambda"] <- (param.neigh[n,"sigmoid"]*df_neigh_n$density)/log(lambda)
        
        
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




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2. Display Interactions  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1. Boxplot of RI  ----

for(country in country.list){
  Box.Plot.Realised.effect <- NULL
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  Box.Plot.Realised.effect <- Realised.Int.Obs.list[[country]] %>%
    mutate(intra.bi = case_when(focal==neigh ~ "INTRA",
                                T~"INTER")) %>%
    ggplot(aes(y=realised.effect.std.lambda,x=as.factor(neigh),
               color=as.factor(neigh),fill=intra.bi)) + 
    #geom_hline(yintercept = 1, color="black",size=1) +
    geom_boxplot() +
    facet_wrap(.~focal,ncol=4,scale="free") +
    labs(y="Effect on intrinsic performance",
         x= "Interacting species",
         color="",
         fill="") +
    #coord_cartesian(ylim=c(0,2),
    #              expand = FALSE) + 
    scale_color_manual(values =col.df$color.name)+
    scale_fill_manual(values =c("white","black"))+
    theme_bw() +
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


#---- 2.2. Network visualasition ----
colours  <- wes_palette("Zissou1", 101, type = "continuous")

plot.legend <- ggplot(data.frame(x=rep(c(1,2,3,4,5),
                                       times=50),
                                 y=1:250),
                      aes(as.factor(x),
                          y,fill=x)) +
  geom_point() +
  geom_line(aes(linewidth=as.factor(x), alpha=as.factor(x))) + 
  scale_linewidth_manual("Standardized median effect",
                         values=c(0.5,1,1.5,2,3),
                         labels=c("< 5%","5%-10%","10%-25%","25%-50%",">50%")) +
  scale_alpha_manual("Standardized median effect",
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
# country = "aus"
#species.focal = "GORO"

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  ratio.mat <-   Realised.Int.Obs.list[[country]]  %>%
    aggregate(sigmoid ~ focal + neigh, function(x) length(which(x<0))/length(x)) %>% # percentage of negative interaction
    spread(neigh,sigmoid) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  
  strength.mat <- Realised.Int.Obs.list[[country]]  %>%
    aggregate(sigmoid ~ focal + neigh, function(x) abs(mean(x))) %>%
    spread(neigh,sigmoid) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  # for the network to have row as receiver
  if(!dim(strength.mat)[1]==dim(strength.mat)[2]){
    df.all.comb <- expand.grid(Code.focal.list,Code.focal.list) %>%
      rename("focal"="Var1") %>%
      rename("neigh"="Var2")
    
    ratio.mat <-Realised.Int.Obs.list[[country]]  %>%
      aggregate(sigmoid ~ focal + neigh,
                function(x) length(which(x<0))/length(x)) %>% # percentage of negative interaction
      right_join(df.all.comb) %>%
      spread(neigh,sigmoid ) %>%
      dplyr::select(-focal) %>%
      as.matrix()
    
    strength.mat <- Realised.Int.Obs.list[[country]]  %>%
      aggregate(sigmoid ~ focal + neigh, function(x) abs(mean(x))) %>%
      right_join(df.all.comb) %>%
      spread(neigh,sigmoid) %>%
      dplyr::select(-focal) %>%
      as.matrix() %>%
      t()# for the network to have row has emitor and columns as receiver
  }
  plot.network.gradient.int(ratio.mat,strength.mat,"",0.01)
  par(mar = rep(0, 4))
  net.country[[country]] <- recordPlot()
  
  
}
net.country$spain
# figures/Network.Realised.effect_spain.pdf
net.country$aus
layout(1)
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
ggsave( plot_grid(net.country$aus,
                  get_legend(plot.legend),
                  ncol = 1,
                  rel_heights =c(1,0.24),
                  labels = "",
                  hjust = 0, vjust = 1),
        file=paste0(home.dic,"figures/Network.Obs.Sigmoid.effect_aus.pdf")) 
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

#---- 2.3. Table sum up ----
str(Realised.Int.list)
sum.up.df <- NULL
for(country in country.list){
  
  sum.up.df.n <- Realised.Int.Obs.list[[country]] %>%
    summarise(mean.effect = (mean(sigmoid)),
              median.effect = (median(sigmoid)),
              var.effect = (var(sigmoid)),
              max.positive.effect = (max(sigmoid)),
              max.negative.effect = (min(sigmoid)),
              count.positive = length(sigmoid[sigmoid >0]),
              count.negative = length(sigmoid[sigmoid  <0]),
              count.total = length(sigmoid)) %>%
    mutate(proportion.positive =(count.positive/count.total),
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           country = country,
           effect ="both",
           species = "All")
  
  
  sum.up.neigh.df.n <- Realised.Int.list[[country]] %>%
    group_by(neigh) %>%
    summarise(mean.effect = (mean(sigmoid)),
              median.effect = (median(sigmoid)),
              var.effect = (var(sigmoid)),
              max.positive.effect = (max(sigmoid)),
              max.negative.effect = (min(sigmoid)),
              count.positive = length(sigmoid[sigmoid >0]),
              count.negative = length(sigmoid[sigmoid  <0]),
              count.total = length(sigmoid)) %>%
    mutate(proportion.positive =(count.positive/count.total),
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           country = country,
           effect ="given")%>%
    rename(species = "neigh")
  
  sum.up.focal.df.n <- Realised.Int.list[[country]] %>%
    group_by(focal) %>%
    summarise(mean.effect = (mean(sigmoid)),
              median.effect = (median(sigmoid)),
              var.effect = (var(sigmoid)),
              max.positive.effect = (max(sigmoid)),
              max.negative.effect = (min(sigmoid)),
              count.positive = length(sigmoid[sigmoid >0]),
              count.negative = length(sigmoid[sigmoid  <0]),
              count.total = length(sigmoid)) %>%
    mutate(proportion.positive =(count.positive/count.total),
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           country = country,
           effect ="received") %>%
    rename(species = "focal")
  
  
  write.csv(bind_rows(sum.up.df.n,sum.up.focal.df.n,sum.up.neigh.df.n),
            file=paste0(home.dic,"results/Sum.up.species",country,".csv"))
  
}


sum.up.focal.df.n <- Realised.Int.list[[country]] %>%
  aggregate(realised.effect ~ focal, median)
view(bind_rows(sum.up.df.n,sum.up.focal.df.n,sum.up.neigh.df.n))


#---- 2.4 Heatmap Fac/Comp ----
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
#---- 2.5. Asymetrie of interactions ----
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

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3. Display Interactions for each year ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.1. Df for summary ----

sum.up.time.df <- NULL
Ratio.Year.list <- NULL
Ratio.Year.Plot <- NULL
country ="spain"
for(country in country.list){
  
  sum.up.df.n <- Realised.Int.Obs.list[[country]] %>%
    group_by(year) %>%
    summarize(mean.effect = (mean(realised.effect)-1),
              median.effect = (median(realised.effect)-1),
              var.effect = (var(realised.effect)),
              max.positive.effect = (max(realised.effect)-1),
              max.negative.effect = (min(realised.effect)-1),
              count.neutre = count(realised.effect ==1),
              count.positive = count(realised.effect >1),
              count.negative = count(realised.effect <1),
              count.total = count(realised.effect>0)) %>%
    mutate(proportion.positive =(count.positive/count.total) ,
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           country = country,
           effect ="both",
           species = "All")
  
  
  sum.up.neigh.df.n <- Realised.Int.Obs.list[[country]] %>%
    group_by(year,neigh) %>%
    summarize(mean.effect = (mean(realised.effect)-1),
              median.effect = (median(realised.effect)-1),
              var.effect = (var(realised.effect)),
              max.positive.effect = (max(realised.effect)-1),
              max.negative.effect = (min(realised.effect)-1),
              count.neutre = count(realised.effect ==1),
              count.positive = count(realised.effect >1),
              count.negative = count(realised.effect <1),
              count.total = count(realised.effect>0)) %>%
    mutate(proportion.positive =(count.positive/count.total) ,
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           effect ="given")%>%
    rename(species = "neigh")
  
  sum.up.focal.df.n <- Realised.Int.Obs.list[[country]] %>%
    group_by(year,focal) %>%
    summarize(mean.effect = (mean(realised.effect)-1),
              median.effect = (median(realised.effect)-1),
              var.effect = (var(realised.effect)),
              max.positive.effect = (max(realised.effect)-1),
              max.negative.effect = (min(realised.effect)-1),
              count.neutre = count(realised.effect ==1),
              count.positive = count(realised.effect >1),
              count.negative = count(realised.effect <1),
              count.total = count(realised.effect>0)) %>%
    mutate(proportion.positive =(count.positive/count.total),
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           country = country,
           effect ="received") %>%
    rename(species = "focal")
  
  year.levels <-levels(as.factor(Realised.Int.Obs.list[[country]]$year))
  
  
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
  
  Ratio.Year.list[[country]] <- bind_rows(sum.up.df.n,sum.up.focal.df.n,sum.up.neigh.df.n) %>%
    left_join(env_pdsi)
  
  write.csv( Ratio.Year.list[[country]] ,
             file=paste0(home.dic,"results/Sum.up.species",country,".csv"))
  # Asymmetry  plot 
  strength.df <- Realised.Int.Obs.list[[country]] %>%
    group_by(year,com_id,neigh,focal) %>%
    summarize(mean.effect = median(sigmoid))
  # plot 
  sym.df <- symmetry.with.positive.for.df(strength.df) 
  #group_by(focal,symmetry) %>%
  #tally() %>%
  save( sym.df,
        file=paste0(home.dic,"results/Sym.df.",country,".Rdata"))
  load(paste0(home.dic,"results/Sym.df.",country,".Rdata"))
  yuplim <- data.frame(country=c("aus","spain"),
                       up =c(1,1))
  sym.sum.df <- sym.df %>%
    mutate(symmetry.2=case_when(is.na(symmetry)~"00",
                                T~symmetry)) %>%
    filter(!symmetry.2 =="00") %>%
    group_by(year,com_id) %>%
    tally() %>%
    rename("total"="n")
  
  asym.plot <- sym.df %>%
    group_by(year,com_id,symmetry) %>%
    tally() %>%
    ungroup() %>%
    left_join(sym.sum.df) %>%
    mutate(percentage= n/total) %>% # 12*12/2- 6
    dplyr::filter(!is.na(symmetry)) %>%
    mutate(symmetry.names = case_when((symmetry =="++"|symmetry=="--"|symmetry=="00")~"sym",
                                      T~"asym")) %>%
    mutate(symmetry= factor(symmetry, levels=c("+-","++","-+","--","00",
                                               "-0","0-","+0","0+"))) %>%
    ggplot(aes(x=as.numeric(year),y=percentage)) +
    stat_summary(fun.y = mean,aes(color=symmetry),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=1) +
    stat_summary(fun.y = mean,aes(color=symmetry),
                 geom = "line",size=1) +
    stat_summary(fun.y = mean,aes(color=symmetry.names),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=1) +
    stat_summary(fun.y = mean,aes(color=symmetry.names),
                 geom = "line",size=1) +
    scale_color_manual(values=c("#CEDB9C","#9C9EDE","#8CA252","#5254A3",
                                #"grey","grey","grey","grey",
                                "#637939","#393B79")) + # stepped2() colors 
    labs(y="Percentage of assymmetry\nin the network",
         x="",
         title="",
         color="") +
    scale_x_continuous(breaks= as.numeric(year.levels)) +
    scale_y_continuous(expand = c(0, 0),
                       limits=c(0,yuplim$up[which(yuplim$country==country)]),
                       labels = scales::percent_format(accuracy = 1)) +
    theme_bw()+
    theme(legend.position = c(0.6, 0.9),
          legend.direction = "horizontal",
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black") ,
          axis.title = element_text(size=14),
          axis.text = element_text(size=14),
          panel.grid.minor = element_blank())
  asym.plot
  
  #ggplot( sym.df , aes(x=symmetry)) + geom_bar()
  # PDSI plot 
  pdsi.plot <- env_pdsi %>%
    gather(PDSI.mean,PDSI.sd,PDSI.Q1,PDSI.Q9,PDSI.Q5,
           key="PDSI",value="value") %>%
    ggplot()+
    annotate("rect",xmin=min(env_pdsi$year),
             xmax=max(env_pdsi$year),
             ymin=0,ymax=max(env_pdsi$PDSI.Q9)+0.1,
             fill="lightblue",
             alpha=0.2) +
    geom_line(aes(x=year,
                  y=value,group=PDSI,linetype=PDSI,
                  color=PDSI),size=1.5) +
    scale_color_manual(values=c("black","grey","grey","grey","brown")) +
    scale_linetype_manual(values=c("solid","dashed","solid","dashed","solid"))+
    scale_x_continuous(breaks= as.numeric(year.levels)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y="PDSI\n") +
    theme_bw() +
    theme(legend.position = c(0.5, 0.9),
          legend.direction = "horizontal",
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black") ,
          axis.title = element_text(size=14),
          axis.text = element_text(size=14),
          panel.grid.minor = element_blank())
  # strength
  strength.plot <- ggplot()+ 
    geom_line(data=Ratio.Year.list[[country]] %>%
                filter(effect =="both" & species == "All"),
              aes(x=year,y=median.effect,
                  color=proportion.negative),
              size=2) +
    geom_line(data=Ratio.Year.list[[country]] %>%
                filter(effect =="received"),
              aes(x=year,y=median.effect, 
                  group=species,
                  color=proportion.negative),
              alpha=0.2) +
    scale_color_gradientn(colours = wes_palette("Zissou1", 
                                                101, 
                                                type = "continuous"),
                          limits = c(0, 1)) +
    labs(y="Percentage of Intrisic fitness changed by\n
         the effect of i on j",
         x="",
         title="",
         color="") +
    scale_x_continuous(breaks= as.numeric(year.levels)) +
    scale_y_continuous(expand = c(0, 0),
                       labels = scales::percent_format(accuracy = 1)) +
    theme_bw()+
    theme(legend.position = c(0.8, 0.2),
          legend.direction = "horizontal",
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black") ,
          axis.title = element_text(size=14),
          axis.text = element_text(size=14),
          panel.grid.minor = element_blank())
  strength.plot
  # ratio
  ratio.plot <- ggplot()+ 
    geom_line(data=Ratio.Year.list[[country]] %>%
                filter(effect =="both" & species == "All"),
              aes(x=year,y=proportion.positive), 
              color="#3B9AB2", size=2) +
    geom_line(data=Ratio.Year.list[[country]] %>%
                filter(effect =="received"),
              aes(x=year,y=proportion.positive, 
                  group=species,color=species),
              color="#3B9AB2",alpha=0.2) +
    geom_line(data=Ratio.Year.list[[country]] %>%
                filter(effect =="both" & species == "All"),
              aes(x=year,y=proportion.negative), 
              color="#F21A00", size=2) +
    geom_line(data=Ratio.Year.list[[country]] %>%
                filter(effect =="received"),
              aes(x=year,y=proportion.negative, 
                  group=species,color=species),
              color="#F21A00",alpha=0.2) +
    labs(y="Percentage of +/- sign\nin the effect of i on j",
         x="",
         title=country) +
    scale_x_continuous(breaks= as.numeric(year.levels)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw()+
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=14),
          panel.grid.minor = element_blank())
  
  Ratio.Year.Plot[[country]]  <- ggarrange(ratio.plot,strength.plot,asym.plot,pdsi.plot,
                                           ncol=2,nrow=2)
  Ratio.Year.Plot[[country]] 
  ggsave(ggarrange(ratio.plot,strength.plot,asym.plot,pdsi.plot,
                   ncol=2,nrow=2),
         width=35,
         height = 25,
         units = "cm",
         file=paste0(home.dic,"figures/RatioOverTime.",country,".pdf"))
  
}
Ratio.Year.Plot[["spain"]] 
Ratio.Year.Plot[["aus"]] 

ratio.plot
#---- 3.2 Network based on time  ----
# country = "spain"
#species.focal = "GORO"
load(paste0(project.dic,"results/Realised.Int.Year.list.RData"))
net.country.year <- list()
for(country in country.list){
  net.country.year.country <- list()
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  year.levels <- levels(as.factor(abundance_df$year))
  
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  for(y in year.levels){
    ratio.mat <-Realised.Int.Obs.list[[country]]  %>%
      filter(year ==y) %>%
      aggregate(realised.effect ~ focal + neigh,
                function(x) length(which(x<1))/length(x)) %>% # percentage of negative interaction
      spread(neigh,realised.effect) %>%
      dplyr::select(-focal) %>%
      as.matrix()
    
    strength.mat <- Realised.Int.Obs.list[[country]]  %>%
      filter(year ==y ) %>%
      aggregate(realised.effect ~ focal + neigh, function(x) abs(mean(x)-1)) %>%
      spread(neigh,realised.effect) %>%
      dplyr::select(-focal) %>%
      as.matrix() # for the network to have row has emitor and columns as receiver
    if(!dim(strength.mat)[1]==dim(strength.mat)[2]){
      df.all.comb <- expand.grid(Code.focal.list,Code.focal.list) %>%
        rename("focal"="Var1") %>%
        rename("neigh"="Var2")
      
      ratio.mat <-Realised.Int.Obs.list[[country]]  %>%
        filter(year ==y) %>%
        aggregate(realised.effect ~ focal + neigh,
                  function(x) length(which(x<1))/length(x)) %>% # percentage of negative interaction
        right_join(df.all.comb) %>%
        spread(neigh,realised.effect) %>%
        dplyr::select(-focal) %>%
        as.matrix()
      
      strength.mat <- Realised.Int.Obs.list[[country]]  %>%
        filter(year ==y ) %>%
        aggregate(realised.effect ~ focal + neigh, function(x) abs(median(x)-1)) %>%
        right_join(df.all.comb) %>%
        spread(neigh,realised.effect) %>%
        dplyr::select(-focal) %>%
        as.matrix() %>%
        t()# for the network to have row has emitor and columns as receiver
    }
    plot.network.gradient.int(ratio.mat,strength.mat,y,0.01)
    net.country.year.country[[y]] <- recordPlot()
    ggsave( plot_grid(net.country.year.country[[y]]),
            file=paste0(home.dic,"figures/Network.",country,"_",y,".pdf"))
  }
  net.country.year[[country]] <-  net.country.year.country
}

net.country.year$aus[1] #figures/Network.aus_2010.pdf
net.country.year$aus[2]#figures/Network.aus_2011.pdf
net.country.year$aus[3]#figures/Network.aus_2014.pdf
net.country.year$aus[4]#figures/Network.aus_2015.pdf
net.country.year$aus[5]#figures/Network.aus_2016.pdf
net.country.year$aus[6]#figures/Network.aus_2017.pdf
net.country.year$aus[7] #figures/Network.aus_2018.pdf
net.country.year$aus[8]#figures/Network.aus_2020.pdf
net.country.year$aus[9]#figures/Network.aus_2022.pdf
net.country.year$aus[10]#figures/Network.aus_2023.pdf
plot_grid(net.country.year.country[[1]],
          net.country.year.country[[2]],
          labels=c("2010","2011" ,"2014", "2015", "2016" ,
                   "2017", "2018","2020","2022" ,"2023"),
          ncol= 2)
net.country.year$spain[1]#figures/Network.spain_2015.pdf
net.country.year$spain[2]#figures/Network.spain_2016.pdf
net.country.year$spain[3]#figures/Network.spain_2017.pdf
net.country.year$spain[4]#figures/Network.spain_2018.pdf
net.country.year$spain[5]#figures/Network.spain_2019.pdf
net.country.year$spain[6]#figures/Network.spain_2020.pdf
net.country.year$spain[7]#figures/Network.spain_2021.pdf
net.country.year$spain[8]#figures/Network.spain_2022.pdf
net.country.year$spain[9]#figures/Network.spain_2023.pdf



rm(ratio.mat)
rm(strength.mat)

#---- 3.3. Realised heatmap across time heat map ----
Realised.Int.Year.plotlist <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  Realised.Int.Year.df <- Realised.Int.Year.list[[country]] %>%
    filter(realised.effect<1000)
  Realised.Int.Year.sum <- Realised.Int.Year.df  %>%
    group_by(neigh,focal,year)%>%
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
  
  
  Realised.Int.Year.neigh <- Realised.Int.Year.df %>%
    group_by(neigh,year)%>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           focal = "total") %>%
    as.data.frame() 
  
  Realised.Int.Year.focal <- Realised.Int.Year.df %>%
    group_by(focal,year)%>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           neigh="total") %>%
    as.data.frame() 
  
  Realised.Int.Year.total <- Realised.Int.Year.df %>%
    group_by(year) %>%
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
  
  head(Realised.Int.Year.sum)
  Realised.Int.Year.plotlist[[country]] <- Realised.Int.Year.sum %>%
    bind_rows(Realised.Int.Year.neigh,
              Realised.Int.Year.focal,
              Realised.Int.Year.total) %>%
    mutate(neigh=factor(neigh, levels=c(Code.focal.list,"total")),
           focal = factor(focal,levels=c("total",rev(Code.focal.list)))) %>%
    mutate(year=as.numeric(year)) %>%
    left_join(env_pdsi_dflist[[country]]) %>%
    ggplot(aes(x =neigh,y =focal)) + 
    geom_tile(aes(fill = RE.ratio,
                  alpha=abs(RE.mean-1))) + 
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))   +
    theme_few() +
    annotate("rect",xmin = 0.5, xmax = 13.5, ymin = 1- .5,
             ymax = 1 + .5, color="black",fill="transparent") + 
    annotate("rect",ymin = 0.5, ymax = 13.5, xmin = 13- .5,
             xmax = 13 + .5, color="black",fill="transparent") + 
    labs(alpha ="Mean effect",
         fill="Ratio of Comp/Fac",
         x="neihg",y="focal") +
    scale_alpha(range=c(0.5,1)) +
    facet_wrap(.~ year,nrow=2 ) +
    #facet_wrap(.~ PDSI.mean,nrow=2 ) +
    theme(legend.position ="bottom",
          axis.title = element_blank(),
          axis.text.x=element_text(angle=90,size=14),
          axis.text.y=element_text(size=14),
          axis.ticks =element_blank(),
          axis.line = element_blank())
}
Realised.Int.Year.plotlist[["aus"]] # figures/Heatmap.PDSI.aus.pdf
Realised.Int.Year.plotlist[["spain"]] # figures/Heatmap.PDSI.spain.pdf

