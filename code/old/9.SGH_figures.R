#---- show results as gradient of time/PDSI----

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1. SET UP: Import packages----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
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
#---- 1.2. Import clean data----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))
country.list <- c("aus","spain")

#---- 1.3. Import results----
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
#---- 4. Realised interactions based on tim
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 4.0. Compute realised interactions based on time  ----

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
        
        seq.abun.neigh <- seq(min(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)),
                              max(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)),
                              (max(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T))-min(quantile(neigh.abundance,probs=c(0.25,0.5,0.75), na.rm = T)))/100)
        
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
#---- 4.1. Df for summary ----

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
#---- 4.2 Network based on time  ----
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

#---- 4.3. Realised heatmap across time heat map ----
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


#---- 4.3. Realised heatmap across time with F-S spectrum ----
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
    dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal,year) %>%
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
  #figures/Heatmap.Time.FS.spectrum_spain.pdf
  #figures/Heatmap.Time.FS.spectrum_aus.pdf
  Realised.Int.sum %>%
    mutate(transparency = case_when(RE.mean < 0.05 ~ 0.2,
                                    (RE.mean >0.05 & RE.mean <=0.1) ~ 0.4,
                                    (RE.mean>0.1 &RE.mean <=0.25) ~ 0.6,
                                    (RE.mean >0.25 & RE.mean<=0.5) ~ 0.8,
                                    (RE.mean>0.5 ) ~ 1)) %>%
    mutate(neigh=factor(neigh, levels=c(trait.df$species[order(trait.df$coord.slow.to.fast)])),
           focal = factor(focal,levels=rev(trait.df$species[order(trait.df$coord.slow.to.fast)]))) %>%
    ggplot(aes(x =neigh,y =focal)) + 
    geom_tile(aes(fill = RE.ratio,alpha=abs(RE.mean-1))) + 
    scale_x_discrete(position = "top") +
    facet_wrap(~year) + 
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




# ---- Functions----

symmetry.with.positive.for.df <- function(strength.df){
  names.sp <- levels(as.factor(strength.df$focal))
  yearlevels <- levels(as.factor(strength.df$year))
  sym.mat <- NULL
  df.all.comb <- expand.grid(names.sp,names.sp) %>%
    rename("focal"="Var1") %>%
    rename("neigh"="Var2") 
  for(y in yearlevels){
    #print(y)
    strength.df.n <- strength.df %>%
      dplyr::filter(year == as.numeric(y)) 
    comlevels <- levels(as.factor(strength.df.n$com_id))
    for(comi in comlevels){
      #print(comi)
    strength.df.n <- strength.df %>%
      dplyr::filter(year == as.numeric(y)) %>%
      dplyr::filter(com_id == comi) %>%
      right_join(df.all.comb %>% mutate(year=as.numeric(y),
                                        com_id = comi ),
                 multiple = "all") %>%
      spread(neigh,mean.effect) %>%
      column_to_rownames("focal") %>%
      dplyr::select(-year)  %>%
      mutate_all( ~replace(., is.na(.), 0))
    sym.mat.n <- data.frame(focal = rep(names.sp , each=length(names.sp )),
                            neigh = rep(names.sp , times=length(names.sp )),
                            year= as.numeric(y),
                            com_id = comi)
  for(i in 1:nrow(strength.df.n)){
    for(j in (i+1):nrow(strength.df.n)){ # just to have the upper levels
      if(j>nrow(strength.df.n))next
     # print(paste0(i,j))
      sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                 sym.mat.n$neigh==names.sp[j])] <-NA
      if(i==j) next
      if(is.na(strength.df.n[i,j])) next
      if(is.na(strength.df.n[j,i])) next
      if(strength.df.n[i,j] > 1 & strength.df.n[j,i]>1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                 sym.mat.n$neigh==names.sp[j])] <- "++"
      }
      if(strength.df.n[i,j] > 1 & strength.df.n[j,i]<1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                 sym.mat.n$neigh==names.sp[j])] <- "+-"
      }
      if(strength.df.n[i,j] < 1 & strength.df.n[j,i]>1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                 sym.mat.n$neigh==names.sp[j])] <- "-+"
      }
      if(strength.df.n[i,j] > 1 & strength.df.n[j,i]==1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                   sym.mat.n$neigh==names.sp[j])] <- "+0"
      }
      if(strength.df.n[i,j] < 1 & strength.df.n[j,i]==1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                   sym.mat.n$neigh==names.sp[j])] <- "-0"
      }
      if(strength.df.n[i,j] == 1 & strength.df.n[j,i]>1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                   sym.mat.n$neigh==names.sp[j])] <- "0+"
      }
      if(strength.df.n[i,j] == 1 & strength.df.n[j,i]<1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                   sym.mat.n$neigh==names.sp[j])] <- "0-"
      }
      if(strength.df.n[i,j] < 1 & strength.df.n[j,i]<1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                 sym.mat.n$neigh==names.sp[j])] <- "--"
      }
      if(strength.df.n[i,j] == 1 & strength.df.n[j,i]==1){
        sym.mat.n$symmetry[which(sym.mat.n$focal==names.sp[i] &
                                   sym.mat.n$neigh==names.sp[j])] <- "00"
      }
      }
    }
    sym.mat <- bind_rows(sym.mat,sym.mat.n)
    }
  }
  return(sym.mat)
}
