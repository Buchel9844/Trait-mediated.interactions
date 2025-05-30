
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
#install.packages("ggraph")
library(Polychrome)
#setwd("/home/lbuche/Eco_Bayesian/chapt3")
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
project.dic <- ""
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"

load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))
country.list <- c("aus","spain")
load(paste0(home.dic,"results/Parameters_alpha.RData"))

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1. Compute theoretical interactions  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL
#---- 2.1. Theoretical interactions  ----

Theoretical.Int.list <- list()
areaplot = pi*7.5^2
areaplot = 15 

# for median of individuals 
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  competition_df <- get(paste0("clean.data.",country))[[paste0("competition_",country)]] 
    if(country=="aus"){
      competition_df <- competition_df %>%
        gather(any_of(Code.focal.list),
               key="species",value="individuals") %>%
        mutate(individuals = (individuals/scale)*15) %>%
        dplyr::rename("com_id"="plot") %>%
        aggregate(individuals ~ year + com_id + species + focal, sum)  %>%
        spread(species, individuals) 
    }
  Theoretical.Int.country.df<- NULL
  for(Code.focal in Code.focal.list){ #focal.levels
    lambda = median(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean[,1], na.rm=T)
    
    df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
  
    
    competition_short_focal_df <-   competition_df%>%
      filter(focal ==Code.focal) %>%
      mutate_at(Code.focal.list,as.numeric)
    
    
    SpNames <- names(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
    
    test.sigmoid  <- NULL
      for( neigh.sp in  SpNames){
        print(paste(Code.focal,neigh.sp))

        
        
        neigh.abundance <- competition_short_focal_df %>%
          dplyr::select(neigh.sp) %>%
          dplyr::filter(!is.na(get(neigh.sp)))%>%
          dplyr::filter(get(neigh.sp) > 0)
        
        if(nrow(neigh.abundance)==0) next
        
        density.neigh <- c(0,quantile(neigh.abundance[,1],c(0.1,0.5,0.9))) 
  
    df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
    
    SpNames <- names(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
    

      alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                             neigh.sp]
      
      #alpha_initial  <-  alpha_initial[which(alpha_initial >= quantile(alpha_initial,probs=c(0.10)) &
      #                                         alpha_initial <=  quantile(alpha_initial,probs=c(0.9)))][1:6400]
      alpha_initial  <-  quantile(alpha_initial,probs=c(0.50))
      
      alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                           neigh.sp]
      
      #alpha_slope  <-  alpha_slope[which(alpha_slope>= quantile(alpha_slope,probs=c(0.10)) &
      #                                     alpha_slope <=  quantile(alpha_slope,probs=c(0.9)))][1:6400]
      
      alpha_slope  <- quantile(alpha_slope,probs=c(0.50))
      
      alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                       neigh.sp]
      
      #alpha_c  <-   alpha_c[which( alpha_c >= quantile( alpha_c,probs=c(0.10)) &
      #                              alpha_c <=  quantile( alpha_c,probs=c(0.9)))][1:6400]
      alpha_c  <-quantile( alpha_c,probs=c(0.50)) 
      
      N_opt = Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh.sp]
      
      #N_opt  <-   N_opt[which(N_opt>= quantile( N_opt,probs=c(0.10)) &
      #                          N_opt <=  quantile( N_opt,probs=c(0.9)))][1:6400]
      N_opt  <-quantile( N_opt,probs=c(0.50))
      param.neigh <- data.frame(neigh = neigh.sp, 
                                country = country,
                                alpha_initial = alpha_initial,
                                alpha_slope = alpha_slope,
                                alpha_c=  alpha_c,
                                N_opt_mean = N_opt ,
                                focal=Code.focal,
                                lambda=lambda,
                                density=density.neigh,
                                density.quantile=c("intercept","low","medium","high"))
      
      for (n in 1:4){
        #if(n==1){print(n)}
        df_neigh_n <-  param.neigh[n,]
        
        
        df_neigh_n[,"theoretical.effect"] <- alpha_function4(df_neigh_n$alpha_initial,
                                                  df_neigh_n$alpha_slope,
                                                  df_neigh_n$alpha_c,
                                                  df_neigh_n$density,
                                                  df_neigh_n$N_opt_mean)

        test.sigmoid <- bind_rows(test.sigmoid,df_neigh_n)
        
      }
    }
    
    
    Theoretical.Int.country.df <- bind_rows( Theoretical.Int.country.df,test.sigmoid)
    
    save(Theoretical.Int.country.df,
         file=paste0(project.dic,"results/Theoretical.Int.df_",country,".RData"))
    
  }
  
  #sigmoid.plot.list[[country]] <- sigmoid.plot.list.focal
  Theoretical.Int.list[[country]] <-  Theoretical.Int.country.df
}

save(Theoretical.Int.list,
     file=paste0(project.dic,"results/Theoretical.Int.list.RData"))


load(paste0(project.dic,"results/Theoretical.Int.list.RData"))
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3. Display Interactions  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#---- 3.1. Network function ----
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
  col.vec <- round(as.numeric(ratio.mat[alphamat.pos  > 0]),1) * 10 + 1 
  colours  <- wes_palette("Zissou1", 11, type = "continuous")
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
colours  <- wes_palette("Zissou1", 101, type = "continuous")

plot.legend <- ggplot(data.frame(x=rep(c(1,2,3),
                                       times=50),
                                 y=1:150),
                      aes(as.factor(x),
                          y,fill=x)) +
  geom_point() +
  geom_line(aes(linewidth=as.factor(x), alpha=as.factor(x))) + 

  scale_fill_gradientn("Sign of interactions",#"Ratio of facilitation and competitive effect",
                       colours = colours,
                       breaks=c(1,51,101),
                       labels=c("100% \nFacilitation",
                                "50/50",
                                "100% \nCompetition"),
                       limits=c(1,101)) +
  scale_linewidth_manual("Strength of interactions",#"Standardized median effect",
                         values=c(0.5,1,1.5),#,2,3),
                         labels=c("< 5%","5%-10%","10%-25%"))+#,"25%-50%",">50%")) +
  scale_alpha_manual("Strength of interactions",#"Standardized median effect",
                     values=c(0.2,0.4,0.6),#,0.8,1),
                     labels=c("< 5%","5%-10%","10%-25%"))+#,"25%-50%",">50%")) +
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

#---- 3.2. Network per density ----
net.country <- list()
density.quantile.name <- c("intercept","low","medium","high")
for(country in country.list){
  for(dq in density.quantile.name){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  ratio.mat <-   Theoretical.Int.list[[country]]  %>%
    dplyr::filter(density.quantile==dq) %>%
    dplyr::rename("sigmoid"="theoretical.effect") %>%
    group_by(focal, neigh) %>%
    summarise(count.positive = length(sigmoid[sigmoid > 0]),
              count.negative = length(sigmoid[sigmoid < 0]),
              count.neutral = length(sigmoid[sigmoid  ==0]),
              count.total = length(sigmoid)) %>%
    mutate(sigmoid = count.negative/(count.total-count.neutral)) %>%
    ungroup() %>%
    dplyr::select(-c('count.positive','count.negative','count.neutral','count.total')) %>%
    spread(neigh,sigmoid) %>%
    dplyr::select(-focal) %>% # percentage of negative interaction
    as.matrix()
  
  
  strength.mat <- Theoretical.Int.list[[country]]  %>%
    dplyr::filter(density.quantile==dq) %>%
    dplyr::rename("sigmoid"="theoretical.effect") %>%
    aggregate(sigmoid ~ focal + neigh, function(x) abs(median(x))) %>%
    spread(neigh,sigmoid) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  # for the network to have row as receiver
  if(!dim(strength.mat)[1]==dim(strength.mat)[2]){
    df.all.comb <- expand.grid(Code.focal.list,Code.focal.list) %>%
      rename("focal"="Var1") %>%
      rename("neigh"="Var2")
    
    ratio.mat <-Theoretical.Int.list[[country]]  %>%
      dplyr::filter(density.quantile==dq) %>%
      dplyr::rename("sigmoid"="theoretical.effect") %>%
      group_by(focal, neigh) %>%
      summarise(count.positive = length(sigmoid[sigmoid > 0]),
                count.negative = length(sigmoid[sigmoid < 0]),
                count.neutral = length(sigmoid[sigmoid  ==0]),
                count.total = length(sigmoid)) %>%
      mutate(sigmoid = count.negative/(count.total-count.neutral)) %>%
      ungroup() %>%
      dplyr::select(-c('count.positive','count.negative','count.neutral','count.total')) %>%
      right_join(df.all.comb) %>%
      spread(neigh,sigmoid) %>%
      dplyr::select(-focal) %>% # percentage of negative interaction
      as.matrix()
    
    strength.mat <- Theoretical.Int.list[[country]]  %>%
      dplyr::filter(density.quantile==dq) %>%
      dplyr::rename("sigmoid"="theoretical.effect") %>%
      aggregate(sigmoid ~ focal + neigh, function(x) abs(mean(x))) %>%
      right_join(df.all.comb) %>%
      spread(neigh,sigmoid) %>%
      dplyr::select(-focal) %>%
      as.matrix() %>%
      t()# for the network to have row has emitor and columns as receiver
  }
  plot.network.gradient.int(ratio.mat,strength.mat,"",0.01)
  par(mar = rep(0, 4))
  net.country[[paste0(country,"_",dq)]] <- recordPlot()
  
    }
}
net.country[[paste0("aus_",density.quantile.name[1])]] #figures/networks/Aus_intercept.pdf
net.country[[paste0("aus_",density.quantile.name[2])]]#figures/networks/Aus_low.pdf
net.country[[paste0("aus_",density.quantile.name[3])]]#figures/networks/Aus_medium.pdf
net.country[[paste0("aus_",density.quantile.name[4])]]#figures/networks/Aus_high.pdf

net.country[[paste0("spain_",density.quantile.name[1])]]#figures/networks/Spain_intercept.pdf
net.country[[paste0("spain_",density.quantile.name[2])]]#figures/networks/Spain_low.pdf
net.country[[paste0("spain_",density.quantile.name[3])]]#figures/networks/Spain_medium.pdf
net.country[[paste0("spain_",density.quantile.name[4])]]#figures/networks/Spain_high.pdf


#---- 3.3. Table sum up for intra vs inter ----
density.quantile.name <- c("intercept","low","medium","high") 

str(Realised.Int.list)
sum.up.df <- NULL
std <- function(x) sd(x)/sqrt(length(x))
for(country in country.list){
  
  sum.up.intra.df.n <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(neigh==focal)%>%
    dplyr::rename("sigmoid"="theoretical.effect")%>%
    group_by(density.quantile) %>%
    summarise(mean.effect = (mean(sigmoid)),
              median.effect = (median(sigmoid)),
              std.effect = (std(sigmoid)),
              max.positive.effect = (max(sigmoid)),
              max.negative.effect = (min(sigmoid)),
              count.positive = length(sigmoid[sigmoid >0.001]),
              count.negative = length(sigmoid[sigmoid  <0.001]),
              count.total = length(sigmoid)) %>%
    ungroup()%>%
    mutate(proportion.positive =(count.positive/count.total),
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           country = country,
           effect ="both",
           species = "All",
           interaction="intra") %>%
    mutate(density.quantile.number=as.numeric(factor(density.quantile,
                                                       level=density.quantile.name)))
  
  #view(sum.up.intra.df.n )
  
  sum.up.inter.df.n <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh==focal)%>%
    dplyr::rename("sigmoid"="theoretical.effect")%>%
    group_by(density.quantile) %>%
    summarise(mean.effect = (mean(sigmoid)),
              median.effect = (median(sigmoid)),
              std.effect = (std(sigmoid)),
              max.positive.effect = (max(sigmoid)),
              max.negative.effect = (min(sigmoid)),
              count.positive = length(sigmoid[sigmoid >0]),
              count.negative = length(sigmoid[sigmoid  <0]),
              count.total = length(sigmoid)) %>%
    ungroup()%>%
    mutate(proportion.positive =(count.positive/count.total),
           proportion.negative =(count.negative/count.total),
           proportion.neutre =1-(proportion.positive +proportion.negative),
           country = country,
           effect ="both",
           species = "All",
           interaction="inter") %>%
    mutate(density.quantile.number=as.numeric(factor(density.quantile,
                                   level=density.quantile.name))+4)
  
  #view(sum.up.inter.df.n )
  write.csv(bind_rows(sum.up.intra.df.n,sum.up.inter.df.n) %>%
              mutate(proportion.positive = round(proportion.positive*100,digits=3),
                     proportion.negative= round(proportion.negative*100,digits=3)),
            file=paste0(home.dic,"results/Sum.up.Intra.Inter.",country,".csv"))
  
}

write.csv(read.csv(file=paste0(home.dic,"results/Sum.up.Intra.Inter.aus.csv")) %>%
  mutate(country="Australia") %>%
  bind_rows(read.csv(file=paste0(home.dic,"results/Sum.up.Intra.Inter.spain.csv"))%>%
              mutate(country="Spain"))%>%
    dplyr::select(density.quantile.number,country,density.quantile,interaction,proportion.positive,proportion.negative,
                  mean.effect,std.effect) %>%
    mutate(proportion.positive = round(proportion.positive, digits = 2),
           proportion.negative= round(proportion.negative, digits = 2),
           mean.effect= round(mean.effect, digits = 3),
           std.effect= round(std.effect, digits = 3),
           interaction = case_when(interaction=="intra"~"con-",
                                   interaction=="inter"~"hetero-"),
           density.quantile.number= case_when(country=="Spain"~ density.quantile.number+8,
                                              T ~density.quantile.number)),
  file=paste0(home.dic,"results/SUPP.SUM.UP.csv"))


# side figure for FIG3
color.vec <- c(proportion.positive =wes_palette("Zissou1", 2, type = "continuous")[1],
               proportion.negative =wes_palette("Zissou1", 2, type = "continuous")[2])
country="aus"
read.csv(file=paste0(home.dic,"results/Sum.up.Intra.Inter.",country,".csv")) %>%
  dplyr::filter(density.quantile %in% c("low","high"))%>%
  gather(proportion.positive,proportion.negative,
         key="proportion",value="prop") %>%
  mutate(label.position = case_when(proportion=="proportion.positive" ~ prop/2,
                                    T ~ 0)) %>%
  mutate(label.text = case_when((proportion=="proportion.positive" & prop>0) ~ paste(as.character(round(prop,digits = 2)),"%"), 
                                    T ~ "")) %>%
  filter(density.quantile=="low")%>%
  ggplot(aes(group=proportion,fill=proportion, 
             y= prop, x=interaction)) +
  geom_bar(stat="identity",position="stack") +
  geom_text(aes(label=label.text ,
                y=label.position,
                x=interaction),
            size=20) +
  scale_x_discrete(labels=c("Hetero-","Con-")) +
  scale_fill_manual(values=color.vec) +
  theme_void() +
  theme(axis.text.x=element_text(size=65),
        legend.position = "none",
        strip.text=element_blank())

#figures/networks/bar_aus_high.png
#figures/networks/bar_aus_low.png
#figures/networks/bar_spain_high.png
#figures/networks/bar_spain_low.png  
  


