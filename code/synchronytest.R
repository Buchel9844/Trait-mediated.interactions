# Synchrony project: team Jeli
# Authors: Jeremy Collings and Lisa Buche
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. Setup ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import packages
library(forecast) # for ARIMA function
library(lmtest) # for ARIMA eigen value
library(ppcor) # for variable correlation
library(tsvr) # for synchrony 
library(tidyverse) # for synchrony 
library(ggplot2) 
library(ggthemes)
library(ggpattern)
library(wesanderson)
library(colorspace)
library(broom)
library(wsyn)
library(readxl)
library(corrplot)


abundance_spain_synch <-abundance_spain %>%
  rename("code.plant"=species) %>%
  left_join(plant_code_spain, by="code.plant") %>%
  dplyr::filter( code.analysis %in% final.species.list.spain ) %>%
  aggregate(individuals ~ code.analysis + year + plot + subplot , sum) %>%
  mutate(individuals =individuals,
         com_id = paste(plot,subplot,sep="_")) %>%
  aggregate(individuals ~ year + code.analysis, mean)  %>%
  spread(code.analysis, individuals) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Function------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df.abundance =abundance_spain_synch
nsurrogs=100
#### 2. Make function to have correlation pvalue with surrogate method ####
cor.sur.mtest <- function(df.abundance,nsurrogs, ...) {
  mat <- as.matrix(df.abundance)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      vec.i <- mat[, i]
      vec.j <- mat[, j]
      #adjust to mean=0 (required for surrog function)	
      vec.i <- vec.i-mean(vec.i); vec.j <- vec.j-mean(vec.j) 
      
      #make surrogates
      #choice between surrtype=aaft or fft dependent on normality
      #aaft recommended for non-normal (visually checked)
      #hist(JR);hist(KC)
      sur.vec.i <- surrog(vec.i, nsurrogs, "aaft", syncpres=T) #list
      sur.vec.j <- surrog(vec.j, nsurrogs, "aaft", syncpres=T)
      
      #correlation for each pair of surrogates (one for each site)
      for (S in 1:nsurrogs){
        rhoSurrogs[S] <- cor(sur.vec.i[[S]], sur.vec.j[[S]])
        synchSurrogs[S] <- tsvreq_classic(t(as.matrix(data.frame(Ni=exp(sur.vec.i[[S]]), 
                                                                 Nj=exp(sur.vec.j[[S]])))))
      }
      #real data correlation
      rhoTrue <- cor(vec.i, vec.j, use="complete.obs")
      
      #fraction of surrogate correlations that are greater than true
      pval <- sum(rhoSurrogs > rhoTrue) / nsurrogs
      p.mat[i,j] <- pval
      p.mat[j,i] <- pval
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Synchrony test ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create all combinaison of pair possible to test for synchrony between them
specieslist.sp <- as.data.frame(expand.grid(final.species.list.spain,
                                                 final.species.list.spain)) %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2))


# Find synchrony long and short for all species pair in JR

df.synchrony <- NULL
for ( i in 1:nrow(specieslist.sp)){
  df.synchrony.i <- t(as.matrix(data.frame(Ni= abundance_spain_synch[,specieslist.sp$Var1[i]],
                                           Nj =abundance_spain_synch[,specieslist.sp$Var2[i]]))) 
  # you can have more than two species at a time, but it is a start I guess
  
  
  tryCatch( { vr.trial <- tsvreq_classic(df.synchrony.i); print(res) }
            , error = function(e) {an.error.occured <<- TRUE})
  if(exists("vr.trial")){
    aggresShort <- aggts(vr.trial, vr.trial$ts[vr.trial$ts<4])[[3]]
    aggresLong <- aggts(vr.trial, vr.trial$ts[vr.trial$ts>=4])[[3]]
    if(aggresShort>1|aggresLong> 1){
      synchrony.significance <- "Long and short synchrony"
      if(aggresShort>1 & aggresLong < 1){
        synchrony.significance <-"Short synchrony"
      }
      if(aggresShort<1 & aggresLong > 1){
        synchrony.significance <-"Long synchrony"
      }
    }else{synchrony.significance <- "No"}
  }else{
    aggresLong <- NA
    aggresShort <- NA
    synchrony.significance <- NA
  }
  
  df.synchrony.i <- data.frame( couple = paste(specieslist.sp$Var1[i],
                                               specieslist.sp$Var2[i],sep="_"),
                                species1=specieslist.sp$Var1[i],
                                species2=specieslist.sp$Var2[i],
                                aggresShort= aggresShort,
                                aggresLong= aggresLong,
                                synchrony.significance=synchrony.significance)
  df.synchrony <- bind_rows(df.synchrony,df.synchrony.i)
}
head(df.synchrony)

write_csv(df.synchrony,
          "results/df.synchrony.csv")
#~~~~~~~~~~~
#---- 3. Visualizing Synchrony ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long.df.synchrony <- pivot_longer(df.synchrony, 
                                     cols = c("aggresShort", "aggresLong"), 
                                     names_to = c("timescale"), 
                                     values_to = c("varRatio"))


VarRatio_hist <-ggplot(data = long.df.synchrony[which(long.df.synchrony$varRatio != 2), ], 
                          aes(x = varRatio, fill = timescale)) +
  geom_histogram(position = "identity", alpha = .5) + 
  xlab("Variance Ratio") + 
  theme_classic(base_size = 15) + ylab("Frequency") + 
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) + 
  scale_fill_manual(values = c("#FFBE0B", "#7340A0"), 
                    name = "Timescale", labels = c("Long", "Short")) 
VarRatio_hist
ggsave("figures/VarRatio_hist.pdf",
       plot = VarRatio_hist )


table(long.df.synchrony$varRatio[which(long.df.synchrony != 2)] > 1, 
      long.df.synchrony$timescale[which(long.df.synchrony != 2)])

# generally less synchrony, but maybe marginally more at short timescales??


VarRatio <- ggplot(data = long.df.synchrony[which(long.df.synchrony$varRatio != 2), ], 
                      aes(x = timescale, y = varRatio)) +
  geom_violin() + geom_jitter() + 
  xlab("Timescale") + 
  theme_classic(base_size = 15) + ylab("Variance Ratio") +
  scale_x_discrete(labels = c("Long", "Short")) + 
  geom_hline(yintercept = 1, linetype = "dashed", size = 1)
VarRatio
ggsave("figures/VarRatio.pdf",
       plot = VarRatio)

# generally much more synchrony, especially at shorter timescales

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 4. Networks ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(igraph)

####### Short -----

# Create data
df.synchrony.matrix.short <- df.synchrony %>%
  dplyr::filter(species1 != species2) %>%
  dplyr::select(species1,species2,aggresShort) %>%
  spread(species2,aggresShort) %>%
  column_to_rownames("species1")%>%
  as.matrix()

network.short.pos <- unname(abs(df.synchrony.matrix.short))
diag(network.short.pos)=0 
network.short <- graph_from_adjacency_matrix(network.short.pos>0 & (network.short.pos>1.2| network.short.pos<0.8),
                                             mode     = "undirected",
                                             weighted = TRUE)
widths <- E(network.short)$weight*2
color.mat <-  as.numeric(network.short.pos[which(network.short.pos>0 & (network.short.pos>1.2| network.short.pos<0.8))])
color.mat <- unique(color.mat)
color.mat[which(color.mat > 1)] <- 1
color.mat[which(color.mat < 1)] <- 3
E(network.short)$lty <- color.mat

pdf(file = "figures/short.network.pdf", width = 6, height = 6)

plot(network.short , layout=layout.circle,
     main="Short synchrony",
     vertex.label = colnames(df.synchrony.matrix.short),
     vertex.frame.color = "transparent",
     vertex.label.family="Helvetica",
     #vertex.label.cex = 1,
     #vertex.label.dist= 2,
     #vertex.label.degree = c(pi/2,pi/2,-pi/2),
     vertex.label.color = "black",
     vertex.color = "grey80",
     vertex.size = 30,
     edge.width = widths,
     edge.color = "black",
     edge.arrow.size = 0.5,
     edge.curved = TRUE)

dev.off()
####### Long -----

# Create data
df.synchrony.matrix.long <- df.synchrony %>%
  dplyr::filter(species1 != species2) %>%
  dplyr::select(species1,species2,aggresLong) %>%
  spread(species2,aggresLong) %>%
  column_to_rownames("species1")%>%
  as.matrix()

network.long.pos <- unname(abs(df.synchrony.matrix.long))
diag(network.long.pos)=0 
network.long <- graph_from_adjacency_matrix(network.long.pos>0 & (network.long.pos>1.2| network.long.pos<0.8),
                                             mode     = "undirected",
                                             weighted = TRUE)
widths <- E(network.long)$weight*2
color.mat <-  as.numeric(network.long.pos[which(network.long.pos>0 & (network.long.pos>1.2| network.long.pos<0.8))])
color.mat <- unique(color.mat)
color.mat[which(color.mat > 1)] <- 1
color.mat[which(color.mat < 1)] <- 3
E(network.long)$lty <- color.mat

pdf(file = "figures/long.network.pdf", width = 6, height = 6)

plot(network.long , layout=layout.circle,
     main="Long synchrony",
     vertex.label = colnames(df.synchrony.matrix.long),
     vertex.frame.color = "transparent",
     vertex.label.family="Helvetica",
     #vertex.label.cex = 1,
     #vertex.label.dist= 2,
     #vertex.label.degree = c(pi/2,pi/2,-pi/2),
     vertex.label.color = "black",
     vertex.color = "grey80",
     vertex.size = 30,
     edge.width = widths,
     edge.color = "black",
     edge.arrow.size = 0.5,
     edge.curved = TRUE)

dev.off()

