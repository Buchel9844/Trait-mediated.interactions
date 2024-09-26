
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
