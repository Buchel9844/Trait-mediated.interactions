
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import packages----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages("HDInterval")
library("HDInterval")
#install.packages("tidyverse")
library("tidyverse")
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

setwd("/home/lbuche/Eco_Bayesian/chapt3")
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Visualisation Species interactions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2.0 Load data ----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))

RawData <- list()
Parameters <- list()
country.list <- as.character(args[1]) #c("aus","spain")

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  for(Code.focal in Code.focal.list){ #focal.levels
    
  load(file =paste0(home.dic,"results/Parameters_",Code.focal,"_",country,".Rdata"))
  #assign(paste0("parameter_",Code.focal),
  #       parameter)
  Parameters[[paste(country,"_",Code.focal)]] <- parameter
  }
}

# environemtnal data
spain_env_pdsi<- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv"))

spain_env_pdsi_med <- spain_env_pdsi %>%
  dplyr::filter(month >3 & month < 8) %>%
  aggregate(spain_pdsi ~ year, median) 

year.levels <- names(parameter$df_lambda_sd)


#---- 2.1 Lambda ----
color.year <- data.frame(year=c("2015","2016","2017","2018","2019","2020","2021"),
                         col.value=colorblind_pal()(8)[2:8])
plot.lambda <- list()

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  for(Code.focal in Code.focal.list){ #focal.levels
    
  year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
  
  plot.lambda[[Code.focal]] <- Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd %>%
    mutate_at(year.levels, ~ rowSums(cbind(., Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean))) %>%
    gather(key="year", value="lambda_sd") %>%
    mutate(year=as.numeric(year)) %>%
    left_join(spain_env_pdsi_med, by="year") %>%
    ggplot(aes(y=lambda_sd, x=as.numeric(spain_pdsi),
               group=spain_pdsi)) +
    geom_boxplot(aes(fill=as.factor(year),color=as.factor(year))) +
    geom_hline(yintercept=median(Parameters[[Code.focal]]$df_lambda_mean[,1])) + 
    labs(title = Code.focal, #title="Intrinsic growth rate of LEMA across years,\n as their mean PDSI",
         y= "intrinsic growth rate (log)",
         x="PDSI",fill="year") +
    theme_bw() +
    scale_color_manual(values= color.year[which(color.year$year %in% year.levels),"col.value"]) +
    scale_fill_manual(values= color.year[which(color.year$year %in% year.levels),"col.value"]) +
    guides(fill=guide_legend(nrow = 1,
                             direction="horizontal",
                             byrow = TRUE,
                             title.hjust = 0.1),
           color="none") +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           strip.background = element_blank(),
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
plot.lambda.all <- ggarrange(plotlist=plot.lambda,
                             common.legend = T,
                             legend = "bottom")
ggsave(plot.lambda.all,
       file=paste0(home.dic,"figures/plot.lambda.pdf"))


  
#---- 2.2. Sigmoid representation ----
source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL
sigmoid.plot.list <- list()
Param.sigm.df <- list()

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  sigmoid.plot.list.focal<- list()
  Param.sigm.country.df <- NULL
  for(Code.focal in Code.focal.list){ #focal.levels
    
  
  df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
  
  Sp.names = colnames(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
  year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
  print(country,Code.focal)
  
  test.sigmoid.all<- NULL
  test.sigmoid  <- NULL
  
  for( neigh in Sp.names){
      print(country)
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
        if(n==1){print(n)}
        df_neigh_n <- data.frame(density=c(0:10),param.neigh[n,])
        
        
        df_neigh_n[,"sigmoid"] <- alpha_function4(df_neigh_n$alpha_initial,
                                                  df_neigh_n$alpha_slope,
                                                  df_neigh_n$alpha_c,
                                                  df_neigh_n$density,
                                                  df_neigh_n$N_opt_mean)
        
        test.sigmoid <- bind_rows(test.sigmoid,df_neigh_n)
        
      }
    }
  limits.y <- c(round(min(test.sigmoid$sigmoid),digits=1)-0.1,
                round(max(test.sigmoid$sigmoid),digits=1)+0.1)
  
  sigmoid.plot.list.focal[[Code.focal]] <- ggplot(test.sigmoid,
                                            aes(y=sigmoid,x=density)) + 
    stat_smooth(stat_smooth(color = "black", size = 0.5, level = 0.999)) +
    theme_bw() + 
    scale_x_continuous(breaks = c(0,5,10))+
    scale_y_continuous(limits=limits.y) +
    geom_hline(yintercept=0, color="black") +
    labs(y="",fill="",color="",
         #title=Code.neigh,
         x=paste0("")) +
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
           axis.text.y=  element_blank(), #element_text(size=12),
           axis.title.x= element_blank(),#element_text(size=12),
           axis.title.y= element_blank(),#element_text(size=12),
           title=element_text(size=12))
  
  Param.sigm.country.df <- bind_rows( Param.sigm.country.df,test.sigmoid)
  }
  sigmoid.plot.list[[country]] <- sigmoid.plot.list.focal
  Param.sigm.df[[country]] <- Param.sigm.country.df
}

save( Param.sigm.df,
     file="results/Param.sigmoid.rData")
save( sigmoid.plot.list,
      file="results/Plot.sigmoid.rData")

load(file=paste0(home.dic,"results/Param.sigmoid.rData"))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


#---- 2.2.1 Raw sigmoid illustration ----

ggsave(ggarrange(plotlist = sigmoid.plot.list[["aus"]],
                 #labels=species.spain,
                 common.legend = T, ncol=1,
                 legend="bottom"),
       width=21,
       heigh=40,
       units = "cm",
       file=paste0(home.dic,"figures/AUS.sigmoid.plot.pdf"))

ggsave(ggarrange(plotlist = sigmoid.plot.list[["spain"]],
                 #labels=species.spain,
                 common.legend = T, ncol=1,
                 legend="bottom"),
       width=21,
        heigh=40,
        units = "cm",
        file=paste0(home.dic,"figures/Spain.sigmoid.plot.pdf"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. Quantile and mean effect per year and family ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(Code.focal in c("LEMA")){ #focal.levels
  for(country in c("All")){ #year.levels
    
    # load(paste0(home.dic,"results/chapt3/Inclusion",Code.focal,"_",country,".Rdata"))
  }
}



range.spmatrix <- function(SpMatrix,year.veclevels){
  SpMatrix.out <- NULL
  SpMatrix.year <- as.data.frame(SpMatrix) %>%
    mutate(year = as.numeric(year.veclevels)) 
  for(n in levels(as.factor(year.veclevels))){
    SpMatrix.n <- SpMatrix.year %>%
      dplyr::filter(year == n) %>%
      summary() %>%
      as.data.frame() %>%
      select(!"Var1") %>%
      rename("neigh"=Var2) %>%
      separate(Freq, sep=":", c("param","value")) %>%
      mutate_if(is.character, str_trim) %>%
      mutate_if(is.factor, str_trim) %>%
      spread(param, value) %>%
      mutate(year=n)
    
    SpMatrix.out <- bind_rows( SpMatrix.out, SpMatrix.n)
  }
  return(SpMatrix.out)
}


SpMatrix.out <- range.spmatrix(SpMatrix,year.veclevels) 


range.effect <- NULL
for( n in levels(as.factor(test.sigmoid$neigh))){
  for( y in levels(as.factor(SpMatrix.out$year))){
    print(paste(n," ",y))
    test.sigmoid.n <- test.sigmoid %>% 
      dplyr::filter( neigh == n  & year == y & density==0) %>%
      select(!"density") %>%
      unique()
    
    SpMatrix.out.n <- SpMatrix.out %>% 
      mutate(neigh = as.character(neigh)) %>%
      dplyr::filter( neigh == n  & year == y)  %>%
      select("Min.","1st Qu.", "Mean","Median","3rd Qu.","Max.") %>%
      gather("range", "density") %>%
      mutate(range = c("Min","1stQu","Mean","Median","3rdQu","Max"))
    
    range.name <- SpMatrix.out.n$range[1:length(unique(SpMatrix.out.n$density))]
    
    range.out <- alpha_function4(rep(test.sigmoid.n$alpha_initial,each=nrow(SpMatrix.out.n)),
                                 rep(test.sigmoid.n$alpha_slope, each=nrow(SpMatrix.out.n)),
                                 rep(test.sigmoid.n$alpha_c, each=nrow(SpMatrix.out.n)),
                                 rep(as.numeric(SpMatrix.out.n$density),
                                     times=nrow(test.sigmoid.n)),
                                 rep(test.sigmoid.n$N_opt_mean, each=nrow(SpMatrix.out.n)))
    
    range.n.y <-  data.frame(effect.raw = range.out,
                             density = rep(as.numeric(SpMatrix.out.n$density),
                                           times=nrow(test.sigmoid.n)),
                             range = rep(SpMatrix.out.n$range,
                                         times=nrow(test.sigmoid.n))) %>%
      aggregate(effect.raw ~ density + range,  function(x) c(median = median(x), 
                                                             sd = sd(x))) %>%
      mutate(effect.raw.median = .[[3]][,1],
             effect.raw.sd= .[[3]][,2]) %>%
      select(density,range ,effect.raw.median,effect.raw.sd) %>%
      mutate(effect.on.lambda.median = effect.raw.median * density ,
             effect.on.lambda.sd = effect.raw.sd * density,
             neigh=n,
             year=y)
    
    range.effect <-  bind_rows( range.effect, range.n.y )
  }
}

range.effect.wider <- pivot_wider(data = range.effect, 
                                  id_cols = c(neigh,year), 
                                  names_from = range, 
                                  values_from = c("effect.raw.median", "effect.raw.sd",
                                                  "effect.on.lambda.median",
                                                  "effect.on.lambda.sd")) %>%
  as.data.frame()  %>%
  mutate(year=as.numeric(year)) %>%
  left_join(spain_env_pdsi_med, by="year")

# add envi data 
spain_env_pdsi_med <- spain_env_pdsi %>%
  dplyr::filter(month >3 & month < 8) %>%
  aggregate(spain_pdsi ~ year, median) 

range.boxplot <- range.effect %>%
  mutate(year=as.numeric(year)) %>%
  left_join(spain_env_pdsi_med, by="year") %>%
  ggplot(aes( x=spain_pdsi,
              y=effect.on.lambda.median,
              group=as.factor(year),
              color=as.factor(year))) +
  geom_boxplot() +
  labs(y="Effect on LEMA intrinsic growth rate",
       color="year",
       x="PDSI")  +
  facet_wrap(.~neigh,scale="free") +
  scale_color_manual(values=safe_colorblind_palette) +
  theme_bw() 
range.boxplot
ggplotly(range.boxplot)   
ggsave(range.boxplot,
       file="figures/range.boxplot.pdf")


range.pointplot <- range.effect.wider %>%
  ggplot(aes( x=spain_pdsi,group=as.factor(year),
              color=as.factor(year))) +
  #geom_point(aes(y=effect.on.lambda.median_3rdQu),size=4,
  #           position=position_dodge(width=0.5)) +
  # geom_point(aes(y=effect.on.lambda.median_Median),size=4,
  #           shape=17,
  #          position=position_dodge(width=0.5)) +
  geom_pointrange(aes(y=exp(effect.on.lambda.median_Median),
                      ymax=exp(effect.on.lambda.median_3rdQu),
                      ymin=exp(effect.on.lambda.median_1stQu)),
                  alpha=1, 
                  size=1,
                  position=position_dodge(width=0.5)) +
  geom_pointrange(aes(y=exp(effect.on.lambda.median_Mean),
                      ymax=exp(effect.on.lambda.median_Max),
                      ymin=exp(effect.on.lambda.median_Min)),
                  alpha=0.5, 
                  size=1,
                  shape=1,
                  position=position_dodge(width=0.5)) +
  scale_color_manual(values=safe_colorblind_palette) +
  theme_bw() + 
  facet_wrap(.~neigh,scale="free") +
  geom_hline(yintercept=1, color="black") +
  labs(y="Effect on LEMA intrinsic growth rate",
       color="year",
       x="PDSI") +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=20),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=20, angle=66, hjust=1),
         axis.text.y= element_text(size=20),
         axis.title.x= element_text(size=22),
         axis.title.y= element_text(size=22),
         title=element_text(size=16))
range.pointplot
ggplotly(range.pointplot)   
ggsave(range.pointplot,
       file="figures/range.pointplot.pdf")


range.plot <- range.effect.wider %>%
  ggplot(aes( x=spain_pdsi, color=neigh)) +
  #geom_point(aes(y=effect.on.lambda.median_3rdQu),size=4,
  #           position=position_dodge(width=0.5)) +
  # geom_point(aes(y=effect.on.lambda.median_Median),size=4,
  #           shape=17,
  #          position=position_dodge(width=0.5)) +
  geom_pointrange(aes(y=exp(effect.on.lambda.median_Median),
                      ymax=exp(effect.on.lambda.median_3rdQu),
                      ymin=exp(effect.on.lambda.median_1stQu)),
                  alpha=0.9, 
                  size=1,
                  position=position_dodge(width=0.5)) +
  scale_color_manual(values=safe_colorblind_palette) +
  theme_bw() + 
  geom_hline(yintercept=1, color="black") +
  labs(y="Effect on LEMA intrinsic growth rate\n 
       based on mean and max abundances observed",
       x="PDSI",
       color="neighbours'\n identity") +
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
         axis.title.x= element_text(size=22),
         axis.title.y= element_text(size=22),
         title=element_text(size=16))
range.plot
ggplotly(range.plot) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 4. Model of interaction according to PDSI and neighbourhood ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SpMatrix.year <- as.data.frame(SpMatrix) %>%
  mutate(year = as.numeric(year.veclevels)) 


interaction.effect <- NULL
for( n in levels(as.factor(test.sigmoid$neigh))){
  for( y in levels(as.factor(SpMatrix.out$year))){
    print(paste(n," ",y))
    test.sigmoid.n <- test.sigmoid %>% 
      dplyr::filter( neigh == n  & year == y & density==0) %>%
      select(!"density") %>%
      unique()
    
    SpMatrix.year.n <- SpMatrix.year %>% 
      dplyr::select(all_of(n),year) %>%
      dplyr::filter( year == as.numeric(y)) 
    
    
    interaction.effect.out <- alpha_function4(rep(test.sigmoid.n$alpha_initial,each=nrow(SpMatrix.year.n)),
                                              rep(test.sigmoid.n$alpha_slope, each=nrow(SpMatrix.year.n)),
                                              rep(test.sigmoid.n$alpha_c, each=nrow(SpMatrix.year.n)),
                                              rep(as.numeric(SpMatrix.year.n[,n]),
                                                  times=nrow(test.sigmoid.n)),
                                              rep(test.sigmoid.n$N_opt_mean, each=nrow(SpMatrix.year.n)))
    
    interaction.effect.n.y <-  data.frame(effect.raw =  interaction.effect.out,
                                          density = rep(as.numeric(SpMatrix.year.n[,n]),
                                                        times=nrow(test.sigmoid.n)))  %>%
      mutate(effect.on.lambda = effect.raw * density,
             neigh=n,
             year=y)
    
    interaction.effect <-  bind_rows(interaction.effect , interaction.effect.n.y )
  }
}

ggplot(interaction.effect, aes(x = exp(effect.on.lambda))) +
  geom_histogram(colour = "#8B5A00", fill = "#CD8500") +
  theme_bw() 

spain_env_pdsi<- read.csv("results/spain_env_pdsi.csv")

spain_env_pdsi_med <- spain_env_pdsi %>%
  dplyr::filter(month >3 & month < 8) %>%
  aggregate(spain_pdsi ~ year, median) 

interaction.effect.env <- interaction.effect %>%
  mutate(year=as.numeric(year)) %>%
  left_join(spain_env_pdsi_med, by="year") %>%
  mutate(density=scale(density))

head(interaction.effect.env)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. Abundance across time ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


abundance_spain <- read.csv(paste0("data/abundance_spain.csv"),
                            header = T,stringsAsFactors = F, sep=",",
                            na.strings=c("","NA"))

plant_code_spain <- read.csv(paste0( "data/plant_code_spain.csv"),
                             header = T, stringsAsFactors = F, sep=",",
                             na.strings = c("","NA"))

abundance_spain_short <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") 


colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", 
                             "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#D55E00",
                             "#882255", "#661100", "#6699CC", 
                             "#888888","#009E73","#0072B2","#E69F00")
scales::show_col(safe_colorblind_palette)

abundance_spain_plot <- abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  ggplot(aes(x=as.character(year), y = individuals,
             group=as.factor(code.analysis),
             color=as.factor(code.analysis))) + 
  stat_summary(fun = mean,
               fun.min = function(x) quantile(x,0.025), 
               fun.max = function(x) quantile(x,0.975), 
               geom = "pointrange",size=2) +
  stat_summary(fun = mean,
               geom = "line",size=1) +
  labs(color="Groups",
       y="Averaged number of individuals \n in 1meter squarred plot (log10)",
       x="year",
       title="Density over time of annual plants in Caracoles") +
  #coord_cartesian( xlim = NULL, ylim = c(0,500),expand = TRUE, default = FALSE, clip = "on") +
  scale_color_manual(values=safe_colorblind_palette) +
  scale_y_log10() +
  #scale_y_continuous(limits=c(0,10)) +
  theme_bw() +
  guides(fill=guide_legend(nrow = 1,
                           direction="horizontal",
                           byrow = TRUE,
                           title.hjust = 0.1)) +
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
abundance_spain_plot
plotly::ggplotly(abundance_spain_plot)
library(plotly)
#figures/abundance_spain_plot.pdf


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
