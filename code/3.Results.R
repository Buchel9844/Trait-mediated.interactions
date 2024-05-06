
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

home.dic <- ""
year.int = "All"    
Code.focal = "LEMA"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Visualisation Species interactions ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load(file= paste0(home.dic,"results/Parameters_",Code.focal,"_All.Rdata"))
load(file= paste0(home.dic,"results/Inclusion",Code.focal,"_",year.int,".Rdata"))


year.levels <- names(parameter$df_lambda_sd)

# Let's look at lambdas
parameter$df_lambda_sd %>%
  mutate_at(year.levels, ~ rowSums(cbind(., parameter$df_lambda_mean))) %>%
  gather(key="year", value="lambda_sd") %>%
  ggplot() +
  geom_boxplot(aes(y=lambda_sd, x=year)) +
  geom_hline(yintercept=median(parameter$df_lambda_mean[,1])) + 
  labs(main="Intrinsic growth rate of LEMA across years",
       y= "intrinsic growth rate: lambda") +
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

df_alpha_generic_param = parameter$df_alpha_generic_param

Sp.names = names(df_alpha_generic_param[1:11])

source("code/PopProjection_toolbox.R")

test.sigmoid  <- NULL
for( neigh in Sp.names){
  for(year.int in year.levels){
    print(neigh)
    print(year.int)

alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                           neigh]
if(Inclusion$Inclusion_alpha_initial[year.int,neigh]>0){
alpha_initial = alpha_initial + parameter$df_alpha_initial[,paste(year.int,neigh,sep="_")]
}

alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                       neigh]
if(Inclusion$Inclusion_alpha_slope[year.int,neigh]>0){
  alpha_slope = alpha_slope + parameter$df_alpha_slope[,paste(year.int,neigh,sep="_")]
}

alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                     neigh]
if(Inclusion$Inclusion_c[year.int,neigh]>0){
  alpha_c = alpha_c + parameter$df_alpha_c[,paste(year.int,neigh,sep="_")]
}


param.neigh <- data.frame(neigh = neigh, 
                          year = year.int,
                          alpha_initial = alpha_initial,
                          alpha_slope = alpha_slope,
                          alpha_c=  alpha_c,
                          N_opt_mean = parameter$df_N_opt_mean[,neigh])

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
}

save(test.sigmoid,
     file="results/test.sigmoid.rData")
load("results/test.sigmoid.rData")
family.to.keep.spain <- levels(as.factor(test.sigmoid$neigh))
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

safe_colorblind_palette[family.to.keep.spain=="conspecific"]<- "black"

#---- 2.1 Raw sigmoid functions per year and family ----
sigmoid.plot <- ggplot(test.sigmoid, aes(y=sigmoid, x=density,group=neigh)) + 
  geom_smooth(aes(color=neigh,fill=neigh), size=2) +
  scale_color_manual(values=safe_colorblind_palette) +
  scale_fill_manual(values=safe_colorblind_palette) +
  theme_bw() + 
  geom_hline(yintercept=0, color="black") +
  facet_wrap(year~.,nrow=1) +
  labs(y="per capita effect on LEMA",
       x="Neighbours' density",
       color="neighbours'\n identity",
       fill="neighbours'\n identity") +
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
sigmoid.plot 
  
ggsave(sigmoid.plot ,
       "figures/sigmoid.plot.pdf")

#---- 2.1 Min/Max effect per year and family ----
for(Code.focal in c("LEMA")){ #focal.levels
  for(year.int in c("All")){ #year.levels
    
# load(paste0(home.dic,"results/chapt3/Inclusion",Code.focal,"_",year.int,".Rdata"))
  }
}



range.spmatrix <- function(SpMatrix,year.veclevels){
  SpMatrix.out <- NULL
  SpMatrix.year <- as.data.frame(SpMatrix) %>%
    mutate(year = as.numeric(year.veclevels)) 
  for(n in levels(as.factor(year.veclevels))){
  SpMatrix.n <- SpMatrix.year %>%
    filter(year == n) %>%
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
      select("Min.", "Max.", "Mean") %>%
      gather("range", "density") %>%
      mutate(range = c("min","max","mean"))
      
    range.name <-  SpMatrix.out.n$range[1:length(unique(SpMatrix.out.n$density))]
    
    range.out <- alpha_function4(rep(test.sigmoid.n$alpha_initial,each=3),
                    rep(test.sigmoid.n$alpha_slope,each=3),
                    rep(test.sigmoid.n$alpha_c,each=3),
                    rep(as.numeric(SpMatrix.out.n$density),
                        times=nrow(test.sigmoid.n)),
                    rep(test.sigmoid.n$N_opt_mean,each=3))
    
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
  as.data.frame()

range.plot <- range.effect.wider %>%
  ggplot(aes( x=year,group=neigh,color=neigh)) +
  geom_point(aes(y=exp(effect.on.lambda.median_mean)),size=4,
             position=position_dodge(width=0.5)) +
  geom_point(aes(y=exp(effect.on.lambda.median_max)),size=4,
             shape=17,
             position=position_dodge(width=0.5)) +
  geom_errorbar(aes(y=exp(effect.on.lambda.median_max),
                    ymax=exp(effect.on.lambda.median_mean),
                    ymin=exp(effect.on.lambda.median_max)),
                width = 1,
                size=5, 
                position=position_dodge(width=0.5)) +
  scale_color_manual(values=safe_colorblind_palette) +
  theme_bw() + 
  geom_hline(yintercept=1, color="black") +
  labs(y="Effect on LEMA intrinsic growth rate\n based on mean and max abundances observed",
       x="year",
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
ggsave(range.plot,
       "figures/range.plot.pdf")

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
  left_join(plant_code_spain) %>%
    mutate( family = tolower(family)) %>%
  
    dplyr::filter(!is.na(family)) %>%
    dplyr::filter(family %in% family.to.keep.spain) %>%
    dplyr::filter(!is.na(individuals )) %>%
    aggregate(individuals ~ code.plant + year + plot+ subplot+
              species + family, max)
  
  

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77",  "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

abundance_spain_plot <- ggplot() + 
    stat_summary(data=abundance_spain_short ,
                 aes(x=as.character(year), y = individuals,
                     group=as.factor(family),color=as.factor(family)),
                 fun = mean,
                 fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x), 
                 geom = "pointrange",size=2) +
    stat_summary(data=abundance_spain_short,
                 aes(x=as.character(year), y = individuals,
                     group=as.factor(family),color=as.factor(family)),
                 fun = mean,
                 geom = "line",size=1) +
    stat_summary(data=abundance_spain_short[which(abundance_spain_short$code.plant=="LEMA"),],
                 aes(x=as.character(year), y = individuals),
                 fun = mean, group="LEMA",
                 fun.min = function(x) mean(x) - sd(x), 
                 fun.max = function(x) mean(x) + sd(x),
                 color="black",geom = "pointrange",size=2) +
    stat_summary(data=abundance_spain_short[which(abundance_spain_short$code.plant=="LEMA"),],
                 aes(x=as.character(year), y = individuals),
                 fun = mean, group="LEMA",
                 color="black", geom = "line", size=1) +
    labs(color="family",y="Mean number of \nindividuals in 1meter squarred plot",
         x="year",
         title="Density over time of annual plants in Caracoles") +
    #coord_cartesian( xlim = NULL, ylim = c(0,500),expand = TRUE, default = FALSE, clip = "on") +
    scale_color_manual(values=safe_colorblind_palette) +
  scale_y_log10() +
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
abundance_spain_plot
plotly::ggplotly(community_plot)
  
#figures/abundance_spain_plot.pdf


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

