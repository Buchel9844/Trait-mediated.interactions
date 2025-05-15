#old stuff 




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
# 1. Data prep
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#abundance_aus_var <-abundance_aus.preclean   %>%
# aggregate(individuals ~ species + year, function(x) sd(x[!x==0 & !is.na(x)])) 

# for Australia, sample 10 community for each year based on the observation mean and var of the species
#n.year = levels(as.factor(abundance_aus.preclean$year))
#abundance_aus.clean <- data.frame(year=rep(n.year,each=10),
#                            com_id = rep(1:10,times=length(n.year)))

#for(n in final.species.list.aus){
# abundance_aus.clean[,n]<- NA
# for(y in  n.year){
#  abundance.vec.realistic.value <- rlnorm(100,meanlog=as.numeric((abundance_aus_med %>% filter( year==y & species ==n))$individuals),
#                                          sdlog=as.numeric((abundance_aus_var %>% filter( year==y & species ==n))$individuals))
#  if(sum(is.na(abundance.vec.realistic.value))>1){
#    abundance.vec.realistic.value <- rlnorm(100,meanlog=as.numeric(mean((abundance_aus_med %>% filter(species ==n))$individuals,na.rm=T)),
#                                             sdlog=as.numeric(mean((abundance_aus_var %>% filter(species ==n))$individuals,na.rm=T)))
#    
#  }
#  abundance_aus.clean[which(abundance_aus.clean$year==y),n] <- abundance.vec.realistic.value[abundance.vec.realistic.value< 100 & 
#                                                                                                abundance.vec.realistic.value> 0.001][1:10]
#   
# }
#}

#----4.2.Category based on reproduction strategy----
plant_pol.traits_aus <- aus_traits_df %>%
  mutate(mean.seed.mass.mg = case_when(is.na(mean.seed.mass.mg)~seed.dry.mass,
                                       T~mean.seed.mass.mg)) %>%
  dplyr::select(final.code,
                flower.size.numb,
                Fecundity,
                surv,mean.seed.mass.mg )  %>%
  column_to_rownames("final.code")
str(plant_pol.traits_aus )
view(plant_pol.traits_aus)

# PCA analysis
#https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/
res.pol.traits <- MFA(plant_pol.traits_aus, 
                      group = c(1,1,1,1),#(1,1,1,1,1,1,1,1,1,1), 
                      type = c("s","s","s","s"),#c("n","n","n","s","s","s","s","s","n","n"),
                      name.group = c("flower","fecundity","survival.seed","seed.mass"),
                      #"flower.size","flower.colour","flowers.per.plant","Dispersal",
                      #"Fecundity","germ",
                      #"surv",
                      # "mean.seed.mass.mg","seed.dry.mass"),
                      graph = T)

eig.val <- get_eigenvalue(res.pol.traits)
head(eig.val)
fviz_contrib(res.pol.traits, choice = "quanti.var", axes = 1, top = 20,palette = "jco")
fviz_contrib(res.pol.traits, choice = "quanti.var", axes = 2, top = 20,palette = "jco")
fviz_mfa_var(res.pol.traits, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)

fviz_mfa_ind(res.pol.traits, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
# cluster analysis
library(vegan)
plant_pol.traits_aus.dist <- vegdist(plant_pol.traits_aus %>% scale() %>% 
                                       as.data.frame(),
                                     na.rm = T,method="euclidean")

hc <- hclust(plant_pol.traits_aus.dist ,method = "average", members = NULL) 
plot(hc)

#flower.size.numb. range 0-1;1-2; > 2


aus_pol_grouping <-  plant_pol.traits_aus %>%
  rownames_to_column("final.code") %>%
  mutate(group.cluster =case_when(final.code %in% c("ARCA","POAR","WAAC","POLE") ~ 1,
                                  final.code %in% c("PEAI") ~ 2,
                                  final.code %in% c("GORO","GOBE","LARO") ~ 3,
                                  final.code %in% c("PLDE","HYGL","MIMY","TRCY") ~ 4)) %>%
  mutate(name.cluster =case_when(group.cluster==1 ~ "Big flower",
                                 group.cluster ==2~ "Grass-like",
                                 group.cluster ==3~ "Medium flower with big seeds",
                                 group.cluster ==4~ "Multiple small flowers with medium seeds")) %>%
  mutate(flower.size.cluster =case_when( final.code == "PEAI" ~ NA,
                                         flower.size.numb>2 ~ "Big flower",
                                         flower.size.numb <= 2 & flower.size.numb> 1 ~ "Medium flower",
                                         flower.size.numb <= 1 ~ "Small flower")) 



#----4.3.Category based on above ground strategy----

plant_abovetraits_aus <- aus_traits_df %>%
  dplyr::select(final.code,
                #aboverground.biomass.mg,
                height.mm,sla.mm2.mg,
                width.longest.mm, width.90.from.longest.mm,
                #delta.C13.discrimination.permill
  )  %>%
  column_to_rownames("final.code")
str(plant_abovetraits_aus)
view(plant_abovetraits_aus)

res.abovetraits.aus <- MFA(plant_abovetraits_aus, 
                           group = rep(1, times= ncol(plant_abovetraits_aus)), 
                           type = rep("s", times= ncol(plant_abovetraits_aus)), # s = quantitative , n= factorial
                           name.group = colnames(plant_abovetraits_aus),
                           graph = T)

eig.val <- get_eigenvalue(res.abovetraits.aus)
head(eig.val)
fviz_contrib(res.abovetraits.aus, choice = "quanti.var", axes = 1, top = 20,palette = "jco")
fviz_contrib(res.abovetraits.aus, choice = "quanti.var", axes = 2, top = 20,palette = "jco")
fviz_mfa_var(res.abovetraits.aus, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)

fviz_mfa_ind(res.abovetraits.aus, partial = "all",
             #col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
fviz_mfa_ind(res.abovetraits.aus,col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
# cluster analysis
plant_abovetraits_aus.dist <- vegdist(plant_abovetraits_aus %>% scale() %>% as.data.frame(),
                                      na.rm = T,method="euclidean")
hc <- hclust(plant_abovetraits_aus.dist ,method = "average", members = NULL) 
plot(hc)

hist(plant_abovetraits_aus$sla.mm2.mg)
hist(aus_above_grouping$area)
hist(aus_above_grouping$volume)
aus_above_grouping <-  plant_abovetraits_aus  %>%
  rownames_to_column("final.code") %>%
  mutate(area = 0.5*width.longest.mm*height.mm) %>%
  mutate(volume = (1/3)*width.longest.mm*height.mm*width.90.from.longest.mm) %>%
  mutate(group.cluster = case_when(final.code %in% c("ARCA","GORO","TRCY") ~ 1, 
                                   final.code %in% c("POAR","WAAC","LARO","HYGL") ~ 2, 
                                   final.code %in% c("PLDE","GOBE","POLE","MIMY","PEAI") ~ 3)) %>% 
  mutate(name.cluster = case_when(group.cluster==1 ~ "Creeping",
                                  group.cluster ==2~ "Tall",
                                  group.cluster ==3~ "Medium"))%>%
  mutate(SLA.cluster = case_when(sla.mm2.mg > 30 ~ "Thick leaves",
                                 sla.mm2.mg < 25~ "Fin leaves",
                                 T ~ "Medium leaves")) %>%
  mutate(heigh.cluster = case_when(height.mm > 100 ~ "Tall stem",
                                   height.mm < 50 ~ "Small stem",
                                   T ~ "Medium stem")) %>%
  mutate(Area.cluster = case_when(area > 3000 ~ "Large area",
                                  area < 1000 ~ "Small area",
                                  T ~ "Medium area")) %>%
  mutate(Volume.cluster = case_when(volume < 10000 ~ "Small Volume",
                                    volume > 1000 & volume < 100000 ~ "Medium Volume",
                                    T ~ "Large Volume"))


aus_above_grouping <-  plant_abovetraits_aus  %>%
  rownames_to_column("final.code") %>%
  mutate(group.cluster = case_when(final.code %in% c("ARCA","GORO","TRCY") ~ 1, 
                                   final.code %in% c("POAR","WAAC","LARO","HYGL") ~ 2, 
                                   final.code %in% c("PLDE","GOBE","POLE","MIMY","PEAI") ~ 3)) %>% 
  mutate(name.cluster = case_when(group.cluster==1 ~ "Creeping",
                                  group.cluster ==2~ "Tall",
                                  group.cluster ==3~ "Medium")) 

#----4.4.Category based on below ground strategy----

plant_belowtraits_aus <- aus_traits_df %>%
  #mutate(biomass.ratio = root.biomass/aboverground.biomass.mg) %>% # if big, invest more below groung
  dplyr::select(final.code,
                total.root.length.cm,
                number.of.root.tips,
                root.volume.less.than.0.5mm.diameter.mm3,
                root.biomass,
                srl)  %>%
  column_to_rownames("final.code")
str(plant_belowtraits_aus)
view(plant_belowtraits_aus)

res.belowtraits.aus <- MFA(plant_belowtraits_aus, 
                           group = rep(1, times= ncol(plant_belowtraits_aus)), 
                           type = rep("s", times= ncol(plant_belowtraits_aus)), # s = quantitative , n= factorial
                           name.group = names(plant_belowtraits_aus),
                           graph = T)

eig.val <- get_eigenvalue(res.belowtraits.aus)
head(eig.val)
fviz_contrib(res.belowtraits.aus, choice = "quanti.var", axes = 1, top = 20,palette = "jco")
fviz_contrib(res.belowtraits.aus, choice = "quanti.var", axes = 2, top = 20,palette = "jco")
fviz_mfa_var(res.belowtraits.aus, "quanti.var", palette = "jco", 
             col.var.sup = "violet", repel = TRUE)

fviz_mfa_ind(res.belowtraits.aus, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
# cluster analysis
library(vegan)
plant_belowtraits_aus.dist <- vegdist(plant_belowtraits_aus %>% scale() %>% as.data.frame(),
                                      na.rm = T,method="euclidean")
hc <- hclust(plant_belowtraits_aus.dist ,method = "average", members = NULL) 
plot(hc)

hist(plant_belowtraits_aus$root.biomass)
aus_below_grouping <-  plant_belowtraits_aus  %>%
  rownames_to_column("final.code") %>%
  mutate(group.cluster = case_when(final.code %in% c("POAR","WAAC") ~ 1,
                                   final.code %in% c("PEAI") ~ 2,
                                   final.code %in% c("ARCA") ~ 3,
                                   final.code %in% c("LARO","GOBE","MIMY","HYGL",
                                                     "GORO","POLE","TRCY","PLDE") ~ 4)) %>%
  mutate(name.cluster = case_when(group.cluster==1 ~ "Heavy Taproot",
                                  group.cluster ==2~ "Grass-like/Thin and spread",
                                  group.cluster ==3~ "Long and thin taproot",
                                  group.cluster ==4~ "Small and thin Taproot")) %>%
  mutate(srl.cluster = case_when( srl < 5 ~ "small srl",
                                  srl >10 ~ "high srl",
                                  T ~ "medium srl"))  %>%
  mutate(root.length.cluster = case_when( total.root.length.cm > 50 ~ "deep root",
                                          total.root.length.cm < 20 ~ "small root",
                                          T ~ "medium root")) 


#---- 4.4. Gather trait mat ----
trait.dist_aus.df <- as.data.frame(as.matrix(plant_pol.traits_aus.dist)) %>%
  rownames_to_column(var="focal") %>%
  gather(all_of(final.species.list.aus),
         key="neigh",value="pol.traits")  %>%
  left_join(as.data.frame(as.matrix(plant_belowtraits_aus.dist)) %>%
              rownames_to_column(var="focal") %>%
              gather(all_of(final.species.list.aus),
                     key="neigh",value="below.traits")) %>%
  left_join(as.data.frame(as.matrix(plant_abovetraits_aus.dist)) %>%
              rownames_to_column(var="focal") %>%
              gather(all_of(final.species.list.aus),
                     key="neigh",value="above.traits")) %>%
  left_join(as.data.frame(as.matrix(aus.traits.dist)) %>%
              rownames_to_column(var="focal") %>%
              gather(all_of(final.species.list.aus),
                     key="neigh",value="dist"))



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
# 2.Results analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 1.1. Make data df ----
country = "aus"
Cool.trait.df <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  Focal.group.df <- Realised.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.50 ~ "Receives Fac",
                                 sigmoid > 0.50 ~ "Receives Comp",
                                 T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid")
  
  Neigh.group.df <- Realised.Int.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh, function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.50 ~ "Gives Fac",
                                 sigmoid > 0.50 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")
  
  
  
  if(country=="aus"){
    Neigh.group.df <-Neigh.group.df %>%
      mutate(neigh.org = case_when( sigmoid.neigh < 0.50 ~ "Gives Fac",
                                    sigmoid.neigh > 0.50 ~ "Gives Comp",
                                    T ~ "Gives both"))
    
    Focal.group.df <-   Focal.group.df %>%
      mutate(focal.org = case_when( sigmoid.focal< 0.5 ~ "Receives Fac",
                                    sigmoid.focal > 0.50 ~ "Receives Comp",
                                    T ~ "Receives both"))
  }
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    #if( i =="Fast to slow spectrum") next
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  Sym.group.df  <- NULL
  for(i in 1:11){
    for(j in (i+1):10){
      df.focal.i.n <- Realised.Int.list[[country]] %>%
        dplyr::filter(focal==Code.focal.list[i] & neigh==Code.focal.list[j])
      prop.comp.on.i <- sum(df.focal.i.n$sigmoid < 0)/length(df.focal.i.n$sigmoid)
      df.focal.j.n <- Realised.Int.list[[country]] %>%
        dplyr::filter(focal==Code.focal.list[j] & neigh==Code.focal.list[i])
      prop.comp.on.j <- sum(df.focal.j.n$sigmoid < 0)/length(df.focal.j.n$sigmoid)
      Sym.group.df.n <- data.frame(focal=Code.focal.list[i] , neigh=Code.focal.list[j],
                                   prop.comp.on.i=prop.comp.on.i,prop.comp.on.j=prop.comp.on.j) %>%
        left_join(specific.trait.dist)
      
      Sym.group.df <- bind_rows( Sym.group.df, Sym.group.df.n)
    }
  }
  Sym.group.df <-  Sym.group.df %>%
    mutate(sym = case_when((prop.comp.on.i < 0.5 & prop.comp.on.j < 0.5) ~ "++",
                           (prop.comp.on.i > 0.5 & prop.comp.on.j > 0.5) ~ "--",
                           (prop.comp.on.i > 0.5 & prop.comp.on.j < 0.5) ~ "-+",
                           (prop.comp.on.i < 0.5 & prop.comp.on.j > 0.5) ~ "+-")) %>%
    mutate(trait.dist = case_when(sym == "-+" ~ -trait.dist,
                                  T~ trait.dist))
  
  
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Fast to slow spectrum","Root length","Root tips","Root biomass","Root volume",
                                            "Mean fecundity","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root biomass","Root volume","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Stem height","Canopy shape") ~ "2.Aboveground",
                                         trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Fast to slow spectrum","Root diameter","Root mass density","SRA",
                                            "Water use efficiency", "Mean fecundity","Leaf nitrogen cc","Ratio leafs",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Water use efficiency","Canopy shape","Stem height","Ratio leafs","Leaf area",
                                                      "Leaf nitrogen cc") ~ "2.Aboveground",
                                         trait %in% c("Mean fecundity") ~ "3.Reproduction"))
  }
  trait.dist.df <- Realised.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    group_by(neigh,focal) %>%
    summarise(Sigm.Q5= median(sigmoid),
              Sigm.neg= length(sigmoid[sigmoid<0]),
              Sigm.total = length(sigmoid)#,# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(Sigm.ratio= Sigm.neg/Sigm.total) %>% 
    left_join(specific.trait.dist %>% 
                left_join(Focal.group.df%>% dplyr::select(focal, focal.org,sigmoid.focal), multiple="all") %>%
                left_join(Neigh.group.df %>% dplyr::select(neigh, neigh.org,sigmoid.neigh), multiple="all"),
              relationship ="many-to-many")
  
  sum.trait.dist.median.df <-  trait.dist.df %>%
    mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                 Sigm.ratio< 0.5 ~"Fac",
                                 Sigm.ratio== 0.5 ~"Neutre")) %>%
    group_by(compORfac,trait) %>%
    mutate(Sigm.ratio = abs(Sigm.ratio -0.5)) %>%
    summarise(weighted.median = matrixStats::weightedMedian(trait.dist,Sigm.ratio,na.rm=T),
              weighted.mad = matrixStats::weightedMad(trait.dist,Sigm.ratio,na.rm=T)) %>%
    mutate(weightedQ1= weighted.median-weighted.mad,
           weightedQ9= weighted.median+weighted.mad)
  
  sum.trait.dist.focal.df <- trait.dist.df %>%
    group_by(focal.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))%>%
    ungroup()
  sum.trait.dist.neigh.df <- trait.dist.df %>%
    group_by(neigh.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))
  
  Cool.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    sum.trait.dist.focal.df=sum.trait.dist.focal.df,
    sum.trait.dist.neigh.df=sum.trait.dist.neigh.df,
    specific.trait.dist=specific.trait.dist,
    sum.trait.dist.median.df=sum.trait.dist.median.df,
    Sym.group.df=Sym.group.df)
}

#---- 1.2. Detailed plot of dist and median interaction ----
Cool.detailed.trait.plotlist <- list()
for( country in country.list){
  
  Cool.detailed.trait.plotlist[[paste0(country,"_focal")]] <- ggplot() +
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
               aes(y=trait.dist,
                   x=Sigm.Q5,
                   group=as.factor(focal.org),
                   fill=Sigm.ratio,
                   color=as.factor(focal.org)),
               shape=21,
               size=1,
               alpha=0.2) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        xmin=Q10.sigmoid,
                        xmax=Q90.sigmoid,
                        color=as.factor(focal.org)),
                    shape=15,
                    size=1) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        ymin=Q10.trait.dist,
                        ymax=Q90.trait.dist,
                        color=as.factor(focal.org)),
                    shape=15,
                    size=1) +
    facet_wrap(.~trait,scale="free") + 
    labs(y="Focal trait value - Neigh trait value",
         x="Median sigmoid value given by Neigh")+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual(values = c("#EBCC2A","#F21A00","#3B9AB2"))+
    theme_bw() 
  
  
  
  Cool.detailed.trait.plotlist[[paste0(country,"_neigh")]]  <- ggplot() +
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
               aes(y=trait.dist,
                   x=Sigm.Q5,
                   group=as.factor(neigh.org),
                   fill=Sigm.ratio,
                   color=as.factor(neigh.org)),
               shape=21,
               size=1,
               alpha=0.2) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        xmin=Q10.sigmoid,
                        xmax=Q90.sigmoid,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df,
                    aes(y=median.trait.dist,
                        x=median.sigmoid,
                        ymin=Q10.trait.dist,
                        ymax=Q90.trait.dist,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) +
    facet_wrap(.~trait,scale="free") + 
    labs(y="Focal trait value - Neigh trait value",
         x="Median sigmoid value given by Neigh")+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual(values = c("#EBCC2A","#F21A00","#3B9AB2"))+
    theme_bw() 
}

Cool.detailed.trait.plotlist[[paste0("spain_focal")]]
Cool.detailed.trait.plotlist[[paste0("spain_neigh")]]
Cool.detailed.trait.plotlist[[paste0("aus_focal")]]
Cool.detailed.trait.plotlist[[paste0("aus_neigh")]]

#---- 1.3. Stat Test ----
Test.trait.list <- list()
for( country in country.list){
  trait.dist.df <- Cool.trait.df[[country]]$trait.dist.df
  sum.trait.dist.median.df <- Cool.trait.df[[country]]$sum.trait.dist.median.df
  Test.neigh.trait.df <- NULL
  Test.focal.trait.df <- NULL
  Test.trait.df <- NULL
  for( trait.i in levels(as.factor(trait.dist.df$trait))){
    trait.dist.df.i <-  trait.dist.df %>%
      dplyr::filter(trait==trait.i) %>%
      mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                   Sigm.ratio< 0.5 ~"Fac",
                                   Sigm.ratio== 0.5 ~"Neutre")) %>%
      dplyr::filter(!compORfac =="Neutre") %>%
      dplyr::select(trait,compORfac,trait.dist)
    median.comp.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Comp" &
                                                                sum.trait.dist.median.df$trait==trait.i]
    median.fac.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Fac" &
                                                               sum.trait.dist.median.df$trait==trait.i]
    
    if(median.comp.n > median.fac.n){test.n ="greater"}else{test.n="less"}
    #test.n ="two.sided"
    test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$compORfac,
                        paired = F,
                        alternative = test.n)
    Test.trait.df.n <- data.frame(trait=trait.i,
                                  p.value = test$p.value,
                                  stat= test$statistic,
                                  alternative=test.n
    )
    Test.trait.df <- bind_rows(Test.trait.df,Test.trait.df.n)
    for( i in levels(as.factor(trait.dist.df$focal.org))){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(!focal.org==i )  %>%
        dplyr::select(trait,focal.org,trait.dist)
      if(nlevels(as.factor(trait.dist.df.i$focal.org))<2)next
      #https://www.r-bloggers.com/2020/06/wilcoxon-test-in-r-how-to-compare-2-groups-under-the-non-normality-assumption-2/
      test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$focal.org,
                          paired = F) # I am not sure if they should be paired or not 
      Test.focal.trait.df.n <- data.frame(trait=trait.i,
                                          focal.org.1 = levels(as.factor(trait.dist.df.i$focal.org))[1],
                                          focal.org.2= levels(as.factor(trait.dist.df.i$focal.org))[2],
                                          p.value = test$p.value,
                                          stat= test$statistic
      )
      Test.focal.trait.df <- bind_rows(Test.focal.trait.df,Test.focal.trait.df.n)
    }
    
    for( i in levels(as.factor(trait.dist.df$neigh.org))){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(!neigh.org==i )  %>%
        dplyr::select(trait,neigh.org,trait.dist)
      if(nlevels(as.factor(trait.dist.df.i$neigh.org))<2)next
      test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$neigh.org,
                          paired = F) # I am not sure if they should be paired or not 
      Test.neigh.trait.df.n <- data.frame(trait=trait.i,
                                          neigh.org.1 = levels(as.factor(trait.dist.df.i$neigh.org))[1],
                                          neigh.org.2= levels(as.factor(trait.dist.df.i$neigh.org))[2],
                                          p.value = test$p.value,
                                          stat= test$statistic)
      Test.neigh.trait.df <- bind_rows(Test.neigh.trait.df,Test.neigh.trait.df.n)
    }
  }
  Test.trait.list[[country]] <- list( Test.focal.trait.df= Test.focal.trait.df,
                                      Test.neigh.trait.df=Test.neigh.trait.df,
                                      Test.trait.df=Test.trait.df)
}
Test.trait.list[["spain"]]$Test.trait.df
Test.trait.list[["aus"]]$Test.trait.df
#---- 1.4. Summary plot with dist - rect ----
Cool.rect.trait.plotlist <- list()
country="aus"
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
for( country in country.list){
  Test.focal.trait.df <-   Test.trait.list[[country]]$Test.focal.trait.df %>%
    dplyr::filter(p.value< 0.05) 
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  sum.trait.dist.focal.df<- Cool.trait.df[[country]]$sum.trait.dist.focal.df
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  Focal_sum <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(x=trait.dist,  
                   y=trait,
                   fill= sigmoid.focal,
                   color=focal.org),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df %>%
                      mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))),
                    position = position_dodge(width = 0.7),
                    aes(x= median.trait.dist,
                        y=trait,
                        xmin=Q10.trait.dist,
                        xmax=Q90.trait.dist,
                        color=as.factor(focal.org)),
                    shape=16,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1)) +
    annotate("text", x = 2, y=0.2, 
             label = "Focal has higher \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    annotate("text", x = -2, y=0.2, 
             label = "Focal has lower \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Species grouping based on received interactions")+
    scale_color_manual(values = c("#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_alpha_continuous(range=c(0.1,1)) +
    theme_bw()  +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                               "bold","plain")))
  Focal_sum 
  
  Test.neigh.trait.df <-   Test.trait.list[[country]]$Test.neigh.trait.df %>%
    dplyr::filter(p.value< 0.05) 
  
  Neigh_sum <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp")))%>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(x=-trait.dist, # reverse the distance to values on the right mean higher value of the trait  
                   y=trait,
                   fill= sigmoid.neigh,
                   color=as.factor(neigh.org)),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df %>%
                      mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))), # "Gives both"
                    position = position_dodge(width = 0.7),
                    aes(x=-median.trait.dist,
                        y=trait,
                        xmin=-Q10.trait.dist,
                        xmax=-Q90.trait.dist,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1))+
    annotate("text", x = 2, y=0.2, 
             label = "Neigh has higher \n trait value than focal",
             size=4) + 
    geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    annotate("text", x = -2, y=0.2, 
             label = "Neigh has lower \n trait value than focal",
             size=4) + 
    geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Species grouping based on given interactions")+
    scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                              101, 
                                              type = "continuous"))+
    scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    theme_bw() +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.neigh.df$trait)) %in% levels(as.factor(Test.neigh.trait.df$trait)),
                                               "bold","plain")))
  
  
  
  Cool.rect.trait.plotlist[[paste0(country,"_sum")]] <-
    ggarrange(Focal_sum,Neigh_sum,
              nrow=1,
              common.legend = F)
  
  Test.trait.df <-   Test.trait.list[[country]]$Test.trait.df
  Cool.rect.trait.plotlist[[paste0(country,"_Pairwise_interactions")]] <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                              Sigm.ratio< 0.5 ~"Fac",
                                              Sigm.ratio== 0.5 ~"Neutre")) %>%
                 dplyr::filter(!compORfac =="Neutre") %>%
                 mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
                 mutate(trait = factor(trait, levels = levels(specific.trait.dist$trait))),
               aes(x=trait.dist,  
                   y=trait,
                   #alpha=Sigm.ratio,
                   #fill=compORfac, # focal.org,
                   color=compORfac), # focal.org,
               position = position_dodge(width = 0.7),
               shape=16,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df %>%
                      dplyr::filter(!compORfac =="Neutre") %>%
                      mutate(trait = factor(trait, levels = levels(specific.trait.dist$trait))),
                    position = position_dodge(width = 0.7),
                    aes(x=weighted.median, # median.trait.dist,
                        y=trait,
                        xmin=weightedQ1,#Q10.trait.dist,
                        xmax=weightedQ9,#Q90.trait.dist,
                        color=as.factor(compORfac)),
                    shape=16,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1)) +
    annotate("text", x = 2, y=0.2, 
             label = "Focal has higher \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    annotate("text", x = -2, y=0.2, 
             label = "Focal has lower \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Pairwise interaction mainly")+#"Species grouping based on received interactions")+
    scale_color_manual(values = c("#F21A00","#3B9AB2"),
                       labels=c("Competitive","Facilitative"))+ # "#EBCC2A",
    scale_alpha_continuous(range=c(0.1,1)) +
    theme_bw()  +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          legend.title =element_text(size=16),
          legend.text =element_text(size=14),
          axis.title=element_text(size=14),
          axis.text.x=element_text(color = "gray12", size = 16,
                                   face=ifelse(Test.trait.df$trait %in% Test.trait.df$trait[Test.trait.df$p.value<0.10],
                                               "bold","plain")))
  
}
Cool.rect.trait.plotlist[["aus_sum"]] 
Cool.rect.trait.plotlist[["spain_sum"]] 
Cool.rect.trait.plotlist[["aus_Pairwise_interactions"]]
Cool.rect.trait.plotlist[["spain_Pairwise_interactions"]]
#---- 1.5. Summary plot with dist - circ ----
Cool.circ.trait.plotlist <- list()
country="aus"
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
  Test.focal.trait.df <-   Test.trait.list[[country]]$Test.focal.trait.df
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  sum.trait.dist.focal.df<- Cool.trait.df[[country]]$sum.trait.dist.focal.df
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  
  Focal_sum_circ <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(focal = factor(focal, levels= Focal.group.df$focal[order(Focal.group.df$sigmoid)]))%>%
                 mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(y=trait.dist,  
                   x=trait,
                   fill=sigmoid.focal,
                   color=focal.org),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.focal.df %>%
                      mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))), # "Receives both",
                    position = position_dodge(width = 0.7),
                    aes(y=median.trait.dist,
                        x=trait,
                        ymin=Q10.trait.dist,
                        ymax=Q90.trait.dist,
                        color=as.factor(focal.org)),
                    shape=16,
                    size=1) +
    geom_hline(yintercept=0,color="black") +
    annotate("text",x =3.4, y = 1, label = "Higher trait value",
             angle =12,color = "gray12",size = 3.5)+ # spain: 14, aus: 6
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 2),
                 arrow = arrow(length = unit(0.5, "cm")))+
    annotate("text",x =3.6, y = -1, label = "Lower trait value",
             angle = 12,color = "gray12",size = 3.5)+
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = -2),
                 arrow = arrow(length = unit(0.5, "cm")))+
    scale_y_continuous(breaks=c(0)) + 
    scale_x_discrete(labels=addline_format(unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))) +
    labs(fill="Percentage of Comp given",
         color="Species grouping based on given interactions")+
    scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                              101, 
                                              type = "continuous"))+
    scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    coord_polar() +
    guides(color = guide_legend(title.position = "top"),
           fill = guide_legend(title.position = "top")) + 
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          axis.text.x = element_text(color = "gray12", size = 12,
                                     face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                                 "bold","plain")),
          panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank()) 
  Focal_sum_circ 
  
  Test.neigh.trait.df  <-   Test.trait.list[[country]]$Test.neigh.trait.df 
  
  Neigh_sum_circ <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 #mutate(neigh= factor(neigh, levels= Neigh.group.df$neigh[order(Neigh.group.df$sigmoid)])) %>%
                 mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))) %>%
                 mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(y=-trait.dist,
                   x=trait,
                   fill= sigmoid.neigh,
                   color=as.factor(neigh.org)),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.neigh.df %>%
                      mutate(neigh.org = factor(neigh.org, levels=c("Gives Fac","Gives both","Gives Comp"))), # "Gives both"
                    position = position_dodge(width = 0.7),
                    aes(y=-median.trait.dist,
                        x=trait,
                        ymin=-Q10.trait.dist,
                        ymax=-Q90.trait.dist,
                        color=as.factor(neigh.org)),
                    shape=15,
                    size=1) +
    geom_hline(yintercept=0,color="black") +
    scale_y_continuous(breaks=c(0)) + 
    scale_x_discrete(labels=addline_format(unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))) +
    labs(fill="Percentage of Comp given",
         color="Species grouping based on given interactions")+
    scale_fill_gradientn(colours= wes_palette("Zissou1", 
                                              101, 
                                              type = "continuous"))+
    scale_color_manual(values = c( "#3B9AB2","#EBCC2A","#F21A00"))+ # "#EBCC2A",
    coord_polar() +
    guides(color = guide_legend(title.position = "top"),
           fill = guide_legend(title.position = "top")) + 
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x=element_text(color = "gray12", size = 12,
                                   face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.neigh.df$trait)) %in% levels(as.factor(Test.neigh.trait.df$trait)),
                                               "bold","plain")))
  Neigh_sum_circ
  
  Cool.circ.trait.plotlist[[paste0(country,"_circ")]] <-
    ggarrange(Focal_sum_circ,Neigh_sum_circ,
              nrow=1,
              common.legend = F, labels=c("a. Species as receiver","b. Species as impactor"))
  
  
  Test.trait.df <-   Test.trait.list[[country]]$Test.trait.df %>%
    mutate(labels.pvalue = case_when((p.value <= 0.1 & p.value > 0.05 )~ "*",
                                     (p.value <= 0.05 & p.value > 0.01 )~ "**",
                                     p.value <= 0.01 ~ "***",
                                     T~ ""))
  
  Cool.circ.trait.plotlist[[paste0(country,"_Pairwise_interactions")]] <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                              Sigm.ratio< 0.5 ~"Fac",
                                              Sigm.ratio== 0.5 ~"Neutre")) %>%
                 dplyr::filter(!compORfac =="Neutre") %>%
                 mutate(focal.org = factor(focal.org, levels=c("Receives Fac","Receives both","Receives Comp"))) %>%
                 mutate(trait = factor(trait, levels = levels(specific.trait.dist$trait))),
               #mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(y=trait.dist,  
                   x=trait,
                   #alpha=Sigm.ratio,
                   #fill=compORfac, # focal.org,
                   color=compORfac), # focal.org,
               position = position_dodge(width = 0.7),
               shape=16,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df %>%
                      dplyr::filter(!compORfac =="Neutre"), 
                    position = position_dodge(width = 0.7),
                    aes(y=weighted.median, # median.trait.dist,
                        x=trait,
                        ymin=weightedQ1,#Q10.trait.dist,
                        ymax=weightedQ9,#Q90.trait.dist,
                        color=as.factor(compORfac)),
                    shape=16,
                    size=1) +
    #geom_text(data=Test.trait.df,aes(x=trait,y=3,label=labels.pvalue),size=12)+
    geom_hline(yintercept=0,color="black") +
    scale_y_continuous(breaks=c(0),expand = c(0, 0)) + 
    scale_x_discrete(labels=addline_format(paste(Test.trait.df$trait,Test.trait.df$labels.pvalue))) +
    labs(color="Pairwise interaction mainly")+
    scale_color_manual(values = c("#F21A00","#3B9AB2"),
                       labels=c("Competitive","Facilitative"))+ # "#EBCC2A",
    coord_polar() +
    annotate("text",x =3.4, y = 1.3, label = "Focal>Neigh",
             angle =0,color = "gray12",size = 4.5)+ # spain: 14, aus: 6
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 3),
                 arrow = arrow(length = unit(0.5, "cm")))+
    annotate("text",x =3.64, y = -1.3, label = "Focal<Neigh",#"trait value",
             angle = 0,color = "gray12",size = 4.5)+
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = -3),
                 arrow = arrow(length = unit(0.5, "cm")))+
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           fill = guide_legend(title.position = "top",
                               nrow=2)) + 
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title =element_text(size=20),
          legend.text =element_text(size=20),
          axis.text.x=element_text(color = "gray12", size = 16,
                                   face=ifelse(Test.trait.df$trait %in% Test.trait.df$trait[Test.trait.df$p.value<0.10],
                                               "bold","plain")),
          plot.margin = unit(c(0,0,0,0),"cm"))
  Cool.circ.trait.plotlist[[paste0(country,"_Pairwise_interactions")]]
  
}
Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]] #figures/Cool.circ.pairwise.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]] # figures/Cool.circ.pairwise.trait.aus.pdf

Cool.circ.trait.plotlist[["Pairwise_interactions"]] <-
  ggarrange(Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]],
            Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]],
            nrow=1,legend="bottom",
            common.legend = T, labels=c("a. Spain","b. Australia"),
            font.label = list(size = 20, color = "black", face = "bold", family = NULL))

Cool.circ.trait.plotlist["Pairwise_interactions"] # figures/Cool.circ.pairwise.trait.pdf

Cool.circ.trait.plotlist[[paste0("spain_circ")]] # figures/Cool.circ.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_circ")]] # figures/Cool.circ.trait.aus.pdf


#---- 1.6. Summary plot with symmetry----
Sym.trait.plotlist <- list()
country="aus"
for( country in country.list){
  Sym.group.df <- Cool.trait.df[[country]]$Sym.group.df %>%
    mutate(sym = case_when(sym == "-+" ~ "+-",
                           T ~ sym)) %>%
    dplyr::filter(!is.na(sym))
  Sym.trait.plotlist[[country]] <- ggplot() +
    geom_point(data=Sym.group.df,
               aes(x=abs(trait.dist),  
                   y=trait,
                   color=as.factor(sym)),
               position = position_dodge(width = 0.7),
               shape=21,size=2,alpha=0.8) +
    stat_summary(data=Sym.group.df ,
                 position = position_dodge(width = 0.7),
                 aes(x=abs(trait.dist),
                     y=trait,group=sym,color=as.factor(sym)),
                 fun=mean, 
                 fun.max = function(x) min(3,mean(x) + var(x)), #quantile(x,c(0.9)),
                 fun.min = function(x)  max(0,mean(x) - var(x)), #quantile(x,c(0.1)),
                 geom="pointrange",
                 shape=16,
                 size=1) +
    scale_y_discrete(expand=c(0.1,0.1)) +
    scale_x_continuous(limits=c(0,3)) +
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Sym")+
    scale_color_manual(values = rev(c("#3B9AB2","#EBCC2A","#F21A00")))+ # "#EBCC2A",
    theme_bw()  +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          axis.text.y=element_text(face=ifelse(levels(as.factor(Cool.trait.df[[country]]$sum.trait.dist.focal.df$trait)) %in% levels(as.factor(Test.focal.trait.df$trait)),
                                               "bold","plain")))
  
  
}
Sym.trait.plotlist[["aus"]]
Sym.trait.plotlist[["spain"]]


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.Traits Graph related- Looking at traits for answers---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.1 With different grouping of functional group ----
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

#---- 3.2 Along the fast to slow gradient ----
Realised.Int.FSspectrum.plotlist <- list()

country = "aus"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]  %>%
    rownames_to_column("species")
  if(country=="spain"){
    trait.df <- trait.df %>%
      bind_rows(data.frame(species="PAIN",
                           Fast to slow spectrum =-10))   
  }
  Realised.Int.sum <- Realised.Int.Year.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  #hist(Realised.Int.sum$realised.effect,breaks=150)
  #Realised.Int.FSspectrum.plotlist[[country]] 
  #figures/Heatmap.N15_spain.pdf
  #figures/Heatmap.root.volume.less0.5_aus.pdf
  library(tibble)
  library(hrbrthemes)
  Realised.Int.FSspectrum.plotlist[[country]] <- Realised.Int.sum %>%
    #dplyr::filter(RE.Q5 < 1000) %>%
    mutate(neigh=factor(neigh, levels=c(trait.df$species[order(trait.df$Fast to slow spectrum)])),
           focal = factor(focal,levels=rev(trait.df$species[order(trait.df$Fast to slow spectrum)]))) %>%
    ggplot(aes(x =neigh,y =focal)) + 
    geom_tile(aes(fill = RE.ratio.sig,alpha=abs(RE.Q5.sig))) + 
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))   +
    theme_ipsum() +
    scale_alpha_continuous(range = c(0.4, 1))+
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
  Realised.Int.FSspectrum.plotlist[[country]]
  
}
Realised.Int.FSspectrum.plotlist[["aus"]] #figures/Heatmap.FS.spectrum_aus.pdf

Realised.Int.FSspectrum.plotlist[["spain"]] #figures/Heatmap.FS.spectrum_spain.pdf

#---- 3.3 Along the dist between interactors----
Realised.Int.Dist.plotlist <- list()
country = "spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  
  Realised.Int.sum <- Realised.Int.Year.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  head(trait.df)
  Dist.plot <- Realised.Int.sum %>%
    left_join(trait.df %>%
                dplyr::select(c("focal","neigh","dist"))) %>%
    #dplyr::filter(!dist==0) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.50 ~ "Comp",
                              RE.ratio.sig<=0.50 ~ "Fac",
                              T ~ "Neutral")) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(y =RE.Q5.sig,x =dist)) +
    geom_point(aes(fill=RE.ratio.sig, group=PosNeg,color=PosNeg),
               shape=21,
               alpha=0.5,
               size=4,
               stroke = 2) +
    geom_smooth(aes( group=PosNeg,color=PosNeg),
                method="glm",se = T,size=2) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual( values = c("#F21A00","#3B9AB2","#EBCC2A")) +
    theme_bw()
  library(ggExtra)
  Realised.Int.Dist.plotlist[[country]] <- ggMarginal(Dist.plot , type = "density",
                                                      groupColour = TRUE, groupFill = TRUE)
  Realised.Int.Dist.plotlist[[country]]
}
Realised.Int.Dist.plotlist[["aus"]] #figures/Dist.Q5_aus.pdf

Realised.Int.Dist.plotlist[["spain"]] #figures/Dist.Q5_spain.pdf


Realised.Int.Dist.Year.plotlist <- list()

country = "spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  
  Realised.Int.sum <- Realised.Int.Year.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal,year) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  Realised.Int.Dist.Year.plotlist[[country]] <- Realised.Int.sum %>%
    left_join(trait.df %>%
                dplyr::select(c("focal","neigh","dist"))) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.50 ~ "Comp",
                              RE.ratio.sig<=0.50 ~ "Fac",
                              T ~ "Neutral")) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(y =RE.Q5.sig,x =dist)) +
    geom_point(aes(fill=RE.ratio.sig),position="jitter",
               shape=21,
               alpha=0.2,
               size=2,
               stroke = 2) +
    facet_wrap(.~ year) +
    geom_smooth(aes( group=as.factor(PosNeg),color=as.factor(PosNeg)),
                method="glm",se = T,size=2) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous")) +
    scale_color_manual(values=rev(c(wes_palette("Zissou1", 
                                                2, 
                                                type = "continuous")))) +
    #coord_cartesian(ylim=c(0.5,2),clip="on") +
    theme_bw()
  Realised.Int.Dist.Year.plotlist[[country]]
  
  # Grouped Scatter plot with marginal density plots
  #ggscatterhist(
  #test.df, x = "Fast to slow spectrum", y = "RE.ratio",
  #color = "PosNeg", size = 3, alpha = 0.6,
  #palette = rev(c(wes_palette("Zissou1", 
  #                            2, type = "continuous"))),
  #margin.params = list(fill = "PosNeg", color = "black", size = 0.2))
}
Realised.Int.Dist.plotlist[["aus"]] #figures/Ratio.Dist_aus.pdf
Realised.Int.Dist.Year.plotlist[["aus"]]
Realised.Int.Dist.plotlist[["spain"]] #figures/Ratio.Dist_spain.pdf
Realised.Int.Dist.Year.plotlist[["spain"]]

#---- 3.4 Along the trait dist between interactors----
Realised.Int.Dist.plotlist <- list()
country = "spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  trait.df <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    if( i =="Fast to slow spectrum") next
    specific.trait.dist.n <-outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(Code.focal.list,each=length(Code.focal.list)),
             focal= rep(Code.focal.list,times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  trait.df <- trait.df %>%
    rownames_to_column("focal")
  Realised.Int.sum <- Realised.Int.Year.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(neigh,focal) %>%
    summarise(RE.neg= length(realised.effect[realised.effect<1]),
              RE.total = length(realised.effect),# high mean high competition
              RE.mean = mean(realised.effect, na.rm = T),
              RE.sd = sd(realised.effect,na.rm = T),
              RE.Q1 = quantile(realised.effect, c(0.1), na.rm = T),
              RE.Q5 = quantile(realised.effect, c(0.5), na.rm = T),
              RE.Q9 = quantile(realised.effect,c(0.9),  na.rm = T),
              RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
              RE.mean.sig = mean(sigmoid, na.rm = T),
              RE.sd.sig = sd(sigmoid,na.rm = T),
              RE.Q1.sig = quantile(sigmoid, c(0.1), na.rm = T),
              RE.Q5.sig = quantile(sigmoid, c(0.5), na.rm = T),
              RE.Q9.sig = quantile(sigmoid,c(0.9),  na.rm = T)) %>%
    ungroup() %>%
    mutate(RE.ratio = RE.neg/RE.total,
           RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  
  Ratio.focal <- Realised.Int.Year.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(focal) %>%
    summarise(RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  head(trait.df)
  test.df <-  Realised.Int.Year.list[[country]] %>%
    full_join(specific.trait.dist,by=c("focal","neigh"),multiple = "all") 
  
  Dist.plot <-  Realised.Int.Year.list[[country]] %>%
    full_join(specific.trait.dist,by=c("focal","neigh"),multiple = "all") %>%
    mutate(PosNeg = case_when(sigmoid<0 ~ "Comp",
                              sigmoid>0 ~ "Fac",
                              T ~ "Neutral"),
           grouping.fact = paste0(PosNeg,"_",trait)) %>%
    ggplot(aes(y =sigmoid,x =abs(trait.dist))) +
    geom_point(aes(fill=sigmoid),
               position="jitter",
               shape=21,
               alpha=0.3,
               size=2,
               stroke = 0) +
    geom_smooth(aes(group=trait, #grouping.fact,
                    linetype=trait),
                method="lm") +
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                   51, 
                                                   type = "continuous")),
                         limits=c(-0.35,0.35))+
    scale_color_manual(values = rev(wes_palette("Zissou1", 
                                                2, 
                                                type = "continuous")))+
    facet_wrap(.~focal,scale="free") +
    theme_bw() 
  Dist.plot
  # figures/Dist.trait.aus.pdf
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    #dplyr::filter(focal=="ARCA") %>%
    #dplyr::filter(!dist==0) %>%
    mutate(PosNeg = case_when(RE.ratio.sig>0.51 ~ "Comp",
                              RE.ratio.sig<0.49 ~ "Fac",
                              T ~ "Neutral")) %>%
    #dplyr::filter(!dist==0) %>%
    ggplot(aes(x =trait,y =trait.dist)) +
    geom_point(aes(fill=RE.ratio.sig),
               position="jitter",
               shape=21,
               alpha=0.3,
               size=2,
               stroke = 0) +
    stat_summary(aes(group=as.factor(PosNeg),color=as.factor(PosNeg)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual( values = c("#F21A00","#3B9AB2","#EBCC2A")) +
    facet_wrap(.~focal) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=66))
  Dist.plot
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    #mutate(focal = factor(focal,levels=rev(trait.df$focal[order(trait.df$Fast to slow spectrum)]))) %>%
    mutate(focal = factor(focal,levels=rev(Ratio.focal$focal[order(Ratio.focal$RE.ratio.sig)]))) %>%
    mutate(PosNeg = case_when(RE.Q5.sig<0 ~ "Comp",
                              RE.Q5.sig>0~ "Fac",
                              T ~ "Neutral")) %>%
    ggplot(aes(x =focal,y =trait.dist)) +
    geom_point(aes(fill=RE.Q5.sig),
               position="jitter",
               alpha=0.3,
               shape=21,
               stroke=0,
               size=2) +
    stat_summary(aes(group=as.factor(PosNeg),
                     color=as.factor(PosNeg)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    #scale_shape_manual(values=1:5) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual( values = c("#F21A00","#3B9AB2","#EBCC2A")) +
    theme_bw() +
    facet_wrap(.~trait)+
    theme(axis.text.x = element_text(angle=66))
  Dist.plot
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    #mutate(focal = factor(focal,levels=rev(trait.df$focal[order(trait.df$Fast to slow spectrum)]))) %>%
    mutate(focal = factor(focal,levels=rev(Ratio.focal$focal[order(Ratio.focal$RE.ratio.sig)]))) %>%
    mutate(PosNeg = case_when(RE.Q5.sig<0 ~ "Comp",
                              RE.Q5.sig>0 ~ "Fac",
                              T ~ "Neutral")) %>%
    dplyr::filter(PosNeg=="Fac") %>%
    ggplot(aes(x =trait,y =trait.dist)) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 geom = "line",size=0.5) +
    #scale_shape_manual(values=1:5) +
    scale_color_paletteer_d("MetBrewer::Benedictus") +
    scale_fill_paletteer_d("MetBrewer::Benedictus") +
    theme_bw() 
  Dist.plot
  
  Dist.plot <- Realised.Int.sum %>%
    right_join(specific.trait.dist,by=c("focal","neigh"),
               multiple = "all") %>%
    # mutate(focal = factor(focal,levels=rev(trait.df$focal[order(trait.df$Fast to slow spectrum)]))) %>%
    mutate(focal = factor(focal,levels=rev(Ratio.focal$focal[order(Ratio.focal$RE.ratio.sig)]))) %>%
    mutate(PosNeg =  case_when(RE.Q5.sig<0 ~ "Comp",
                               RE.Q5.sig>0 ~ "Fac",
                               T ~ "Neutral")) %>%
    dplyr::filter(PosNeg=="Comp") %>%
    ggplot(aes(x =trait,y =trait.dist)) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=0.5) +
    stat_summary(aes(group=as.factor(focal),
                     color=as.factor(focal)),
                 fun.y  = function(x) quantile(x,0.5),
                 geom = "line",size=0.5) +
    #scale_shape_manual(values=1:5) +
    scale_color_paletteer_d("MetBrewer::Benedictus") +
    scale_fill_paletteer_d("MetBrewer::Benedictus") +
    theme_bw() 
  Dist.plot
  
}
#---- 3.5 Sem -----
sem.trait.plotlist <- list()
country = "spain"
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.Year.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
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
    as.data.frame()
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  
  sem.df <- NULL
  for( sp in Code.focal.list){
    print(sp)
    test.df <-  Realised.Int.Year.list[[country]] %>%
      left_join(specific.trait.dist,
                by=c("focal","neigh"), 
                relationship = "many-to-many") %>%
      #filter(focal ==sp) %>%
      filter(focal ==sp) %>%
      filter(!focal == neigh) %>%
      aggregate(sigmoid ~ year + density + focal + neigh + trait.dist + trait, median) %>%
      group_by(year, focal, neigh, trait.dist,trait) %>%
      summarise(sigmoid = median(sigmoid),
                density=median(density)) %>%
      ungroup() %>%
      mutate(trait.dist =abs(trait.dist))%>%
      spread(trait,trait.dist) %>%
      left_join(env_pdsi %>%  mutate(year = as.numeric(year)),
                by=c("year"),multiple = "all")  %>%
      as.data.frame()
    names(test.df) <- gsub(" ",".",colnames(test.df))
    
    if(country=="aus"){
      test.df <- test.df  %>%
        dplyr::select_if(~ sum(!is.na(.))>10) %>%
        drop_na()
      var.names <- names(test.df)[names(test.df) %in% c("sla","width.longest","srl","root.length",
                                                        "flower.size.numb","height","number.of.root.tips")]
      
    }else{  test.df <- test.df %>%select_if(~ sum(!is.na(.))>10)
    var.names <- names(test.df)[names(test.df) %in% c("Canopy.shape","Fast to slow spectrum","Leaf.area",
                                                      "Leaf.area.index",#"Mean.fecundity","Ratio.leafs","Root.diameter" ,"Water.use.efficiency","SRA","SRL","Stem.length"
                                                      "Leaf.nitrogen.cc" , "Root.mass.density")]
    }
    hist(test.df$sigmoid,breaks=150)
    sem.median  <- piecewiseSEM::psem(
      glm(as.formula(paste("sigmoid ~ 1 +", paste(c(var.names,"density","Precip.extrem"), collapse= "+"))),
          data=test.df, family= gaussian(link = "identity")),
      glm(density ~ Precip.extrem,
          data=test.df,family= gaussian(link = "log")),
      test.df)
    
    #plot.psemhl(sem.median , correlation = T, layout = "tree") 
    summary(sem.median)
    sem.df.n <- coefs(sem.median)[,1:8] %>%
      mutate(focal=sp)
    
    sem.df <- bind_rows(sem.df,sem.df.n)
  }
  
  Realised.Int.sum <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    group_by(focal) %>%
    summarise(RE.neg.sig= length(sigmoid[sigmoid<0]),
              RE.total.sig = length(sigmoid),# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(RE.ratio.sig = RE.neg.sig/RE.total.sig) %>%
    as.data.frame()
  
  
  trait.df.dist <- get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]] 
  trait.df <- trait.df %>% rownames_to_column("focal")
  head(sem.df)
  sem.trait.plotlist[[country]] <- sem.df %>%
    dplyr::filter(Response=="sigmoid") %>%
    dplyr::filter(P.Value < 0.1) %>%
    #mutate(focal = factor(focal,levels=rev(unique(trait.df.dist$focal[order(trait.df.dist$dist)])))) %>%
    #mutate(focal = factor(focal,levels=rev(unique(trait.df$focal[order(trait.df$Fast to slow spectrum)])))) %>%
    mutate(focal = factor(focal,levels=rev(Realised.Int.sum$focal[order(Realised.Int.sum$RE.ratio.sig)]))) %>%
    #mutate(focal= factor(focal,levels=rev(Realised.Int.sum$neigh[order(Realised.Int.sum$RE.ratio.sig)]))) %>%
    ggplot(aes(x=as.factor(Predictor),
               y=Std.Estimate,
               group=as.factor(focal),
               color=as.factor(focal))) +
    geom_point() + 
    geom_line(size=2) +
    scale_color_paletteer_d("MetBrewer::Benedictus") +
    #scale_color_paletteer_d("MoMAColors::Avedon") +
    #scale_color_manual(values = c(wes_palette("Zissou1",
    #                                       12,type = "continuous")))+
    theme_bw() 
  
}
library(paletteer)
sem.trait.plotlist[["aus"]]
sem.trait.plotlist[["spain"]]


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.Traits Graph with time ---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1. Make data df ----
country = "spain"
Cool.trait.df <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.Year.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
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
    as.data.frame()
  
  Focal.group.df <- Realised.Int.Obs.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal + year, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.30 ~ "Receives Fac",
                                 sigmoid > 0.70 ~ "Receives Comp",
                                 T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid") %>%
    mutate(year=as.integer(year)) %>%
    left_join(env_pdsi%>% mutate(year=as.integer(year)) )
  hist(log(Focal.group.df$sigmoid.focal))
  test.glm <- glmer.nb(sigmoid.focal ~ Precip.extrem + (1|focal) ,data=Focal.group.df)
  summary(test.glm)
  ggplot(Focal.group.df,aes(y=sigmoid.focal, x=Precip.extrem,color=focal)) +
    geom_point()+
    geom_smooth(method="lm",se=F)
  
  Neigh.group.df <- Realised.Int.Obs.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh + year, function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.30 ~ "Gives Fac",
                                 sigmoid > 0.70 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")%>%
    mutate(year=as.integer(year)) %>%
    left_join(env_pdsi%>% mutate(year=as.integer(year)) )
  
  ggplot(Neigh.group.df,aes(y=sigmoid.neigh, x=Precip.extrem)) + geom_point()+
    geom_smooth(method="lm",se=F)
  
  if(country=="aus"){
    Neigh.group.df <-Neigh.group.df %>%
      mutate(neigh.org = case_when( sigmoid.neigh < 0.25 ~ "Gives Fac",
                                    sigmoid.neigh > 0.50 ~ "Gives Comp",
                                    T ~ "Gives both"))
    
    Focal.group.df <-   Focal.group.df %>%
      mutate(focal.org = case_when( sigmoid.focal< 0.25 ~ "Receives Fac",
                                    sigmoid.focal > 0.50 ~ "Receives Comp",
                                    T ~ "Receives both"))
  }
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    #if( i =="Fast to slow spectrum") next
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root mass density","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Stem height","Canopy width","Canopy width 90deg") ~ "2.Aboveground",
                                         trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Water use efficiency","Canopy shape","Stem height","Leaf C to N ratio","Leaf area",
                                                      "Leaf nitrogen cc") ~ "2.Aboveground",
                                         trait %in% c("Mean fecundity") ~ "3.Reproduction"))
  }
  trait.dist.df <- Realised.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    group_by(neigh,focal) %>%
    summarise(Sigm.Q5= median(sigmoid),
              Sigm.neg= length(sigmoid[sigmoid<0]),
              Sigm.total = length(sigmoid)#,# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(Sigm.ratio= Sigm.neg/Sigm.total) %>% 
    left_join(specific.trait.dist %>% 
                left_join(Focal.group.df%>% dplyr::select(focal, focal.org,sigmoid.focal), multiple="all") %>%
                left_join(Neigh.group.df %>% dplyr::select(neigh, neigh.org,sigmoid.neigh), multiple="all"),
              relationship ="many-to-many")
  if(country=="spain"){
    trait.dist.df <- trait.dist.df%>% dplyr::filter(!focal=="PAIN") %>%
      dplyr::filter(!neigh=="PAIN") 
  }
  sum.trait.dist.focal.df <- trait.dist.df %>%
    group_by(focal.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))%>%
    ungroup()
  sum.trait.dist.neigh.df <- trait.dist.df %>%
    group_by(neigh.org,trait) %>%
    summarise(median.trait.dist=median(trait.dist,na.rm=T),
              Q10.trait.dist=quantile(trait.dist,c(0.10),na.rm=T),
              Q90.trait.dist=quantile(trait.dist,c(0.90),na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T))
  
  Cool.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    sum.trait.dist.focal.df=sum.trait.dist.focal.df,
    sum.trait.dist.neigh.df=sum.trait.dist.neigh.df,
    specific.trait.dist=specific.trait.dist)
}
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- PCA results ---
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#---- 1.4.PCA axes----
PCA.trait.df <- list()
density.quantile.name <- c("intercept","low","medium","high")
data_normalized <- list()
res.traits <- list()
country="aus"
library(FactoMineR)
library(factoextra)
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  # make trait data frame
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] %>%
    dplyr::select(-c("Mean fecundity"))
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=trait.dist,
             scaled.trait.dist=as.vector(scale(trait.dist)),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list))) %>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>% 
                  rownames_to_column("focal") %>% dplyr::rename("receiver.trait"=i) %>%
                  mutate(receiver.trait = as.vector(scale(receiver.trait))))%>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>%
                  rownames_to_column("neigh") %>% dplyr::rename("emitter.trait"=i)%>%
                  mutate(emitter.trait = as.vector(scale(emitter.trait))))
    
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
    
    
  }
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","SLA","Stem height"))) 
    toremove <- c("Root diameter",
                  "Flower width","Seed mass",
                  "Leaf area index",
                  "Canopy shape","SRA","Root tips")
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","SLA","Stem height"))) 
    toremove <- c("Root diameter",
                  "Flower width","Seed mass",
                  "Leaf area index",
                  "Canopy shape","Root tips")
  }
  #country="aus"
  #country="spain"
  
  
  data_normalized[[country]] <- specific.trait.dist %>% 
    dplyr::select(trait,receiver.trait,focal) %>%
    dplyr::filter(!trait %in% toremove ) %>%
    unique()%>%
    spread(trait,receiver.trait) %>%
    column_to_rownames("focal")
  data_normalized[[country]] 
  
  #data_normalized[["spain"]] 
  #data_normalized[["aus"]] 
  #https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/
  ifelse(country=="aus",c(2,3,1),c(3,2,1))
  res.traits[[country]] <- MFA(data_normalized[[country]], 
                               group = if(country=="aus"){c(3,2,1)}else{c(3,3,1)},#(1,1,1,1,1,1,1,1,1,1), 
                               type = c("s","s","s"),#c("n","n","n","s","s","s","s","s","n","n"),
                               name.group = c("Root cooperation spectrum",
                                              "Drought sensitivity",
                                              "Plant size"),
                               graph = T)
  
  #FactoMineR::PCA(data_normalized[[country]])
  
  #country="aus"
  #country="spain"
  eig.val <- get_eigenvalue(res.traits[[country]])
  head(eig.val)
  # pdf(paste0("figures/pca/",country,"PCAaxes.pdf"))
  chart.Correlation(  data_normalized[[country]], histogram=TRUE, pch=19)
  
  fviz_contrib(res.traits[[country]],
               choice = "quanti.var", axes = 1, top = 20,palette = "jco")
  fviz_contrib(res.traits[[country]], choice = "quanti.var",
               axes = 2, top = 20,palette = "jco")
  fviz_contrib(res.traits[[country]], choice = "quanti.var", axes = 3, top = 20,palette = "jco")
  fviz_mfa_var(res.traits[[country]], "quanti.var", palette = "jco", 
               col.var.sup = "violet", repel = TRUE)
  
  fviz_mfa_ind(res.traits[[country]], col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE)
  #dev.off()
  # cluster analysis
  library(vegan)
  if(country=="spain"){
    PCA.trait.df[[country]]  <-  data.frame("BelowGroundStrategy"=res.traits[[country]]$ind$coord[,1],
                                            "PlantSize"=res.traits[[country]]$ind$coord[,2],
                                            "AboveGroundStrategy"=res.traits[[country]]$ind$coord[,3],
                                            "BelowGroundTrait"=data_normalized[[country]]$SRA,
                                            "PlantSizeTrait"=data_normalized[[country]]$`Stem height`,
                                            "AboveGroundTrait"=-1*data_normalized[[country]]$SLA)# negative sign to be in parrallel with SPAIN
  }
  
  if(country=="aus"){
    PCA.trait.df[[country]]  <-  data.frame("PlantSize"=-1*res.traits[[country]]$ind$coord[,1],
                                            "AboveGroundStrategy"=-1*res.traits[[country]]$ind$coord[,2],
                                            "BelowGroundStrategy"=res.traits[[country]]$ind$coord[,3],
                                            "BelowGroundTrait"=data_normalized[[country]]$`Root length`,
                                            "PlantSizeTrait"=data_normalized[[country]]$`Stem height`,
                                            "AboveGroundTrait"=data_normalized[[country]]$`C13 water use efficiency`) 
  }
}

#---- 1.5.PCA GLM- FOCAL----
glm.trait.pca.list <- list()
glm.trait.list <- list()
glm.new.data.list <- list()
glm.new.data.trait.list<- list()
for( country in country.list){
  
  trait.pca.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    left_join(as.data.frame(PCA.trait.df[[country]] %>%
                              rownames_to_column("focal")),
              relationship ="many-to-many")
  
  glm.trait.pca.df <- NULL
  glm.trait.df <- NULL
  glm.new.data.df <- NULL
  glm.new.data.trait.df <- NULL
  #glm.new.data.df.3D <- NULL
  #glm.new.data.df.3D.i <- data.frame(AboveGroundStrategy = c(0,seq(min(trait.pca.df$AboveGroundStrategy),max(trait.pca.df$AboveGroundStrategy),by=(max(trait.pca.df$AboveGroundStrategy)-min(trait.pca.df$AboveGroundStrategy))/100)),
  #                                   BelowGroundStrategy = c(0,seq(min(trait.pca.df$BelowGroundStrategy),max(trait.pca.df$BelowGroundStrategy),by=(max(trait.pca.df$BelowGroundStrategy)-min(trait.pca.df$BelowGroundStrategy))/100)))
  #glm.new.data.df.3D.i <- expand.grid(glm.new.data.df.3D.i)
  #glm.new.data.df.3D.i$PlantSize <- glm.new.data.df.3D.i$AboveGroundStrategy + glm.new.data.df.3D.i$BelowGroundStrategy
  #glm.new.data.df.3D.i <- glm.new.data.df.3D.i %>%
  #  dplyr::filter(PlantSize>min(trait.pca.df$PlantSize) &
  #          PlantSize < max(trait.pca.df$PlantSize))
  #glm.new.data.df.i <- data.frame(PlantSize=c(0,seq(min(trait.pca.df$PlantSize),max(trait.pca.df$PlantSize),by=(max(trait.pca.df$PlantSize)-min(trait.pca.df$PlantSize))/100)),# interaction Size-Above
  #                                AboveGroundStrategy = c(0,seq(min(trait.pca.df$AboveGroundStrategy),max(trait.pca.df$AboveGroundStrategy),by=(max(trait.pca.df$AboveGroundStrategy)-min(trait.pca.df$AboveGroundStrategy))/100)),
  #                                BelowGroundStrategy = c(0,seq(min(trait.pca.df$BelowGroundStrategy),max(trait.pca.df$BelowGroundStrategy),by=(max(trait.pca.df$BelowGroundStrategy)-min(trait.pca.df$BelowGroundStrategy))/100)))
  
  #glm.new.data.df.i <- expand.grid(glm.new.data.df.i)
  
  glm.new.data.df.i <- data.frame(PlantSize=c(runif(50000,min = min(trait.pca.df$PlantSize), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                              #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                              runif(50000,min = 0, max = max(trait.pca.df$PlantSize))),# interaction Size-Above # quantile(trait.pca.df$PlantSize,0.66)
                                  AboveGroundStrategy = runif(100000, min =min(trait.pca.df$AboveGroundStrategy),
                                                              max=max(trait.pca.df$AboveGroundStrategy)),
                                  BelowGroundStrategy = runif(100000, min =min(trait.pca.df$BelowGroundStrategy),
                                                              max=max(trait.pca.df$BelowGroundStrategy)),
                                  PlantSizeCat =rep(c("BelowAverage","AboveAverage"),each=50000),
                                  BelowGroundTrait=runif(100000, min =min(trait.pca.df$BelowGroundTrait,na.rm=T),
                                                         max=max(trait.pca.df$BelowGroundTrait,na.rm=T)),
                                  PlantSizeTrait=c(runif(50000,min = min(trait.pca.df$PlantSizeTrait,na.rm=T), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                                   #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                                   runif(50000,min = 0, max = max(trait.pca.df$PlantSizeTrait))),
                                  AboveGroundTrait=runif(100000, min =min(trait.pca.df$AboveGroundTrait,na.rm=T),
                                                         max=max(trait.pca.df$AboveGroundTrait,na.rm=T))) #"AroundAverage",
  glm.new.data.df.i <- glm.new.data.df.i %>%
    dplyr::filter(PlantSize>min(trait.pca.df$PlantSize) &
                    PlantSize < max(trait.pca.df$PlantSize))
  
  
  for(n in density.quantile.name){
    trait.pca.df.i <-  trait.pca.df %>%
      dplyr::filter(density.quantile==n) %>%
      dplyr::select(focal,theoretical.effect,density.quantile,
                    PlantSize,AboveGroundStrategy,BelowGroundStrategy,
                    BelowGroundTrait,PlantSizeTrait,AboveGroundTrait) %>%
      mutate(focal=as.factor(focal))
    
    glm.trait.pca.i <- glmmTMB(theoretical.effect ~  PlantSize*AboveGroundStrategy*BelowGroundStrategy, 
                               trait.pca.df.i,
                               family="gaussian")
    summary(glm.trait.pca.i)
    glm.trait.pca.df.i <- as.data.frame(confint(glm.trait.pca.i)) %>%
      mutate(density.quantile=n) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
                    "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.pca.df <- bind_rows(glm.trait.pca.df,glm.trait.pca.df.i)
    
    
    glm.new.data.df.i$theoretical.effect.pca <- predict(glm.trait.pca.i,
                                                        glm.new.data.df.i)
    
    glm.trait.i <- glmmTMB(theoretical.effect ~  PlantSizeTrait*AboveGroundTrait*BelowGroundTrait, 
                           trait.pca.df.i,
                           family="gaussian")
    summary(glm.trait.i)
    glm.trait.df.i <- as.data.frame(confint(glm.trait.i)) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
                    "Q97.5"="97.5 %") %>%
      mutate(density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.df <- bind_rows(glm.trait.df,glm.trait.df.i)
    
    
    glm.new.data.df.i$theoretical.effect <- predict(glm.trait.i,
                                                    glm.new.data.df.i)
    glm.new.data.df <- bind_rows( glm.new.data.df, glm.new.data.df.i%>%
                                    mutate(density.quantile=n))
    
    
  }
  glm.trait.pca.list[[country]] <-  glm.trait.pca.df
  glm.trait.list[[country]] <-  glm.trait.df
  glm.new.data.list[[country]] <-glm.new.data.df
  glm.new.data.trait.list[[country]] <- glm.new.data.trait.df
}
view(glm.trait.pca.list[["aus"]])  
view(glm.trait.pca.list[["spain"]])
view(glm.trait.list[["aus"]])  
view(glm.trait.list[["spain"]])
#---- 1.6.PCA plot glm-FOCAL----
# Above ground and below ground
glm.plot.list <- list()
# for focal species
for( country in country.list){
  if(country=="aus"){
    limits.vec = c(-0.4,0.4)
    ylabname = "WUE (13C)"
    xlabname = "Root length"
  }else{limits.vec = c(-0.4,0.4)
  ylabname = "WUE (13C, inversed SLA)"
  xlabname = "SRA"}
  focal.pca.df.i <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    aggregate(theoretical.effect~ focal, median) %>%
    left_join(as.data.frame(PCA.trait.df[[country]] %>%
                              rownames_to_column("focal")))%>%
    mutate(PlantSizeCat = case_when(PlantSizeTrait < 0 ~ "BelowAverage", #quantile(trait.pca.df$PlantSize,0.33)
                                    PlantSizeTrait > 0 ~ "AboveAverage", #quantile(trait.pca.df$PlantSize,0.66)
                                    T~"AroundAverage"))
  
  glm.new.data.df.i <- glm.new.data.list[[country]]%>%
    filter(theoretical.effect  > min(Theoretical.Int.list[[country]]$theoretical.effect) &
             theoretical.effect < max(Theoretical.Int.list[[country]]$theoretical.effect))
  glm.plot.list[[country]] <-glm.new.data.df.i  %>%
    dplyr::filter(density.quantile %in%c("intercept"))%>%
    ggplot(aes(y=AboveGroundTrait, #Strategy,
               x=BelowGroundTrait, #Strategy,
               fill=theoretical.effect)) +
    geom_tile(width = 0.04,
              height = 0.04,
              alpha=0.8)+
    geom_point(data=focal.pca.df.i,
               aes(y=AboveGroundTrait, #Strategy,
                   x=BelowGroundTrait, #Strategy,
                   fill=theoretical.effect),
               shape=21,size=3) +
    geom_text(data=focal.pca.df.i,
              aes(y=AboveGroundTrait, #Strategy, 
                  x=BelowGroundTrait, #Strategy,
                  label=focal),
              size=4,hjust=0.5,vjust=-0.5) +
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                   101, 
                                                   type = "continuous")),
                         limits=limits.vec)+
    facet_wrap(density.quantile~PlantSizeCat,nrow=1)+
    theme_clean()+
    labs(title=country,
         y =   ylabname, #"Above ground strategy, \nfrom less to more drought sensitive",
         x =   xlabname, #"Below ground strategy, \nfrom more to less cooperative",
         subtitle = "Plant height", #"Plant size, from high/fast to small/slow growing",
         fill="Predicted interaction effect on \nof focal species with that trait combinaison")+
    theme(legend.position="bottom",
          legend.key.size = unit(1, 'cm'),
          legend.background = element_blank(),
          plot.background = element_blank(),
          axis.text =element_text(size=20),
          axis.title =element_text(size=20),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title.position = "top",
          legend.title =element_text(size=20),
          legend.text =element_text(size=20))
  
  glm.plot.list[[country]]
  
  
}

glm.plot.list[["spain"]] # figures/pca/Flat3D.focal.spain.pdf
glm.plot.list[["aus"]] # figures/pca/Flat3D.focal.aus.pdf

ggarrange(glm.plot.list[["spain"]],
          glm.plot.list[["aus"]],
          common.legend = T,legend="bottom",
          nrow=2, 
          labels=c("a. Spain","b. Australia"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
# figures/pca/Flat3D.focal.intercept.pdf
levels(as.factor(Chapt1_Parameters_values$parameter))
names(Chapt1_Parameters_values)
ParameterValuesForJurg <- Chapt1_Parameters_values %>%
  filter(parameter %in% c("Intraspecific","Plant - plant") &
           year== "2021" &
           complexity.plant =="family") %>%
  dplyr::rename("Iterations"="X" ) %>%
  dplyr::select(Iterations,focal,lambdas,parameter,estimate) %>%
  spread(parameter,estimate)
view(ParameterValuesForJurg)
write.csv(ParameterValuesForJurg,
          file="/Users/lisabuche/Downloads/ParameterValuesForJurg")
#---- 1.7.PCA  GLM-NEIGH----
glm.trait.pca.neigh.list <- list()
glm.trait.neigh.list <- list()
glm.new.data.neigh.list <- list()
glm.new.data.neigh.3D.list<- list()
for( country in country.list){
  
  trait.pca.df <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    left_join(as.data.frame(PCA.trait.df[[country]] %>%
                              rownames_to_column("neigh")),
              relationship ="many-to-many")
  
  
  glm.new.data.df.i <- data.frame(PlantSize=c(runif(50000,min = min(trait.pca.df$PlantSize), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                              #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                              runif(50000,min = 0, max = max(trait.pca.df$PlantSize))),# interaction Size-Above # quantile(trait.pca.df$PlantSize,0.66)
                                  AboveGroundStrategy = runif(100000, min =min(trait.pca.df$AboveGroundStrategy),
                                                              max=max(trait.pca.df$AboveGroundStrategy)),
                                  BelowGroundStrategy = runif(100000, min =min(trait.pca.df$BelowGroundStrategy),
                                                              max=max(trait.pca.df$BelowGroundStrategy)),
                                  PlantSizeCat =rep(c("BelowAverage","AboveAverage"),each=50000),
                                  BelowGroundTrait=runif(100000, min =min(trait.pca.df$BelowGroundTrait,na.rm=T),
                                                         max=max(trait.pca.df$BelowGroundTrait,na.rm=T)),
                                  PlantSizeTrait=c(runif(50000,min = min(trait.pca.df$PlantSizeTrait,na.rm=T), max=0), #max = quantile(trait.pca.df$PlantSize,0.33)),
                                                   #runif(50000,min = quantile(trait.pca.df$PlantSize,0.33), max = quantile(trait.pca.df$PlantSize,0.66)),
                                                   runif(50000,min = 0, max = max(trait.pca.df$PlantSizeTrait))),
                                  AboveGroundTrait=runif(100000, min =min(trait.pca.df$AboveGroundTrait,na.rm=T),
                                                         max=max(trait.pca.df$AboveGroundTrait,na.rm=T))) #"AroundAverage",
  glm.new.data.df.i <- glm.new.data.df.i %>%
    dplyr::filter(PlantSize>min(trait.pca.df$PlantSize) &
                    PlantSize < max(trait.pca.df$PlantSize))
  
  glm.trait.df <- NULL
  glm.trait.pca.df<- NULL
  glm.new.data.df <- NULL
  for(n in density.quantile.name){
    trait.pca.df.i <-  trait.pca.df %>%
      dplyr::filter(density.quantile==n) %>%
      dplyr::select(neigh,theoretical.effect,density.quantile,
                    PlantSize,AboveGroundStrategy,BelowGroundStrategy,
                    PlantSizeTrait,AboveGroundTrait,BelowGroundTrait) %>%
      mutate(neigh=as.factor(neigh))
    
    glm.trait.pca.i <- glmmTMB(theoretical.effect ~  PlantSize*AboveGroundStrategy*BelowGroundStrategy, 
                               trait.pca.df.i,
                               family="gaussian")
    summary(glm.trait.pca.i)
    glm.trait.pca.df.i <- as.data.frame(confint(glm.trait.pca.i)) %>%
      mutate(density.quantile=n) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
                    "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.pca.df <- bind_rows(glm.trait.pca.df,glm.trait.pca.df.i)
    
    
    glm.new.data.df.i$theoretical.effect.pca <- predict(glm.trait.pca.i,
                                                        glm.new.data.df.i)
    glm.trait.i <- glmmTMB(theoretical.effect ~  PlantSizeTrait*AboveGroundTrait*BelowGroundTrait, 
                           trait.pca.df.i,
                           family="gaussian")
    summary(glm.trait.i)
    glm.trait.df.i <- as.data.frame(confint(glm.trait.i)) %>%
      mutate(density.quantile=n) %>%
      rownames_to_column("parameters") %>%
      dplyr::rename("Q2.5"="2.5 %",
                    "Q97.5"="97.5 %") %>%
      mutate(trait=trait.i,
             density.quantile=n,
             signif=case_when((Q2.5<0 & Q97.5<0 )~"*",
                              (Q2.5>0 & Q97.5>0 )~"*",
                              T~""))
    
    glm.trait.df <- bind_rows(glm.trait.df,glm.trait.df.i)
    
    
    glm.new.data.df.i$theoretical.effect <- predict(glm.trait.i,
                                                    glm.new.data.df.i)
    glm.new.data.df <- bind_rows( glm.new.data.df, glm.new.data.df.i%>%
                                    mutate(density.quantile=n))
  }
  glm.trait.pca.neigh.list[[country]] <-  glm.trait.pca.df
  glm.trait.neigh.list[[country]] <-  glm.trait.df
  glm.new.data.neigh.list[[country]] <-glm.new.data.df
}
#---- 1.8.PCA plot glm-NEIGH----
glm.plot.neigh.list <- list()
for( country in country.list){
  if(country=="aus"){
    limits.vec = c(-0.4,0.4)
    ylabname = "13C, WUE"
    xlabname = "Root length"
  }else{limits.vec = c(-0.4,0.4)
  ylabname = "SLA - Inversed 13C, WUE"
  xlabname = "SRA"}
  neigh.pca.df.i <- Theoretical.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>% 
    aggregate(theoretical.effect~ neigh , median) %>%
    left_join(as.data.frame(PCA.trait.df[[country]] %>%
                              rownames_to_column("neigh")))%>%
    mutate(PlantSizeCat = case_when(PlantSizeTrait < 0 ~ "BelowAverage",
                                    PlantSizeTrait > 0 ~ "AboveAverage",
                                    T~"AroundAverage"))
  
  glm.new.data.df.i <- glm.new.data.neigh.list[[country]]%>%
    filter(theoretical.effect  > min(Theoretical.Int.list[[country]]$theoretical.effect) &
             theoretical.effect < max(Theoretical.Int.list[[country]]$theoretical.effect))
  
  glm.plot.neigh.list[[country]] <-glm.new.data.df.i  %>%
    dplyr::filter(density.quantile %in%c("intercept"))%>%
    ggplot(aes(y=AboveGroundTrait, #Strategy, 
               x=BelowGroundTrait, #Strategy, 
               fill=theoretical.effect)) +
    geom_tile(width = 0.04,
              height = 0.04,
              alpha=0.8)+
    geom_point(data= neigh.pca.df.i,
               aes(y=AboveGroundTrait, #Strategy, 
                   x=BelowGroundTrait, #Strategy, 
                   fill=theoretical.effect),
               shape=21,size=3) +
    geom_text(data= neigh.pca.df.i,
              aes(y=AboveGroundTrait, #Strategy,  
                  x=BelowGroundTrait, #Strategy, 
                  label=neigh),
              size=4,hjust=0.5,vjust=-0.5) +
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                   101, 
                                                   type = "continuous")),
                         limits=limits.vec)+
    facet_wrap(density.quantile~PlantSizeCat,nrow=1)+
    theme_clean()+
    labs(title=country,
         y = ylabname, #"Above ground strategy, \nfrom less to more drought sensitive",
         x =  xlabname, #"Below ground strategy, \nfrom more to less cooperative",
         subtitle = "Plant height",#"Plant size, from high/fast to small/slow growing",
         fill="Predicted interaction effect on \nof focal species with that trait combinaison")+
    theme(legend.position="bottom",
          legend.key.size = unit(2, 'cm'),
          legend.background = element_blank(),
          plot.background = element_blank(),
          axis.text =element_text(size=20),
          axis.title =element_text(size=20),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title.position = "top",
          legend.title =element_text(size=20),
          legend.text =element_text(size=20))
  
  glm.plot.neigh.list[[country]]
  
  
}
glm.plot.neigh.list[["aus"]]# figures/pca/Flat3D.neigh.aus.pdf
glm.plot.neigh.list[["spain"]]# figures/pca/Flat3D.neigh.spain.pdf

ggarrange(glm.plot.neigh.list[["spain"]],
          glm.plot.neigh.list[["aus"]],
          common.legend = T,legend="bottom",
          nrow=2, 
          #labels=c("a. Spain","b. Australia"),
          font.label = list(size = 20, color = "black", 
                            face = "bold", family = NULL))
# figures/pca/Flat3D.neigh.intercept.pdf

#----Not use----
data.3d <- glm.new.data.3D.list[[country]]%>%
  dplyr::filter(density.quantile=="intercept") %>%
  dplyr::select(-density.quantile) %>%
  as.matrix()%>%
  as.data.frame()%>%
  mutate_all(function(x) as.numeric(x))

str(data.3d)
ThreeDfig <- plot_ly(x=~PlantSize, y= ~AboveGroundStrategy, 
                     z=~BelowGroundStrategy, data=data.3d,
                     #type="scatter3d", mode="markers", 
                     #color=~theoretical.effect)
                     type='mesh3d',
                     flatshading=T,
                     intensity =~theoretical.effect,
                     colors=colorRamp(rev(wes_palette("Zissou1", 
                                                      round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*200)+1, 
                                                      type = "continuous")))) %>%
  add_trace(x=~PlantSize, y= ~AboveGroundStrategy, 
            z=~BelowGroundStrategy, data=focal.pca.df.i,
            color= ~theoretical.effect,mode = "markers",
            type = "scatter3d",
            marker = list(size = 5),
            colors=colorRamp(rev(wes_palette("Zissou1", 
                                             round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*200)+1, 
                                             type = "continuous"))))
ThreeDfig
library(scatterplot3d)
focal.pca.df.i <-  PCA.trait.df[[country]] %>%
  rownames_to_column("focal")

colors.gradient <- rev(wes_palette("Zissou1", 
                                   round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*200)+1, 
                                   type = "continuous"))
focal.colors.gradient <- colors.gradient[(round(focal.pca.df.i$theoretical.effect*100)+round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*100)+1)]
colors.gradient <- colors.gradient[(round(data.3d[,4]*100)+round(max(abs(min(data.3d[,4])),max(data.3d[,4]))*100)+1)]


s3d <- scatterplot3d(data.3d[,1:3], angle = 60,
                     pch = 16, color=colors.gradient,
                     grid=TRUE,
                     main=country,
                     xlab = "Above ground strategy, from less to more drought sensitive",
                     ylab = "Below ground strategy, from more to less cooperative",
                     zlab = "Plant size, from small/slow to high/fast growing")

s3d$points3d(focal.pca.df.i$AboveGroundStrategy,
             focal.pca.df.i$BelowGroundStrategy,
             focal.pca.df.i$PlantSize, 
             col = focal.colors.gradient , type = "h", pch = 8)

text(s3d$xyz.convert(focal.pca.df.i$AboveGroundStrategy,
                     focal.pca.df.i$BelowGroundStrategy,
                     focal.pca.df.i$PlantSize),
     labels=focal.pca.df.i$focal,
     col='tomato', cex=2)
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1. Compute realised interactions  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL
#---- 1.1. Realised interactions across years  ----

Realised.Int.list <- list()
areaplot = 15
# for median of individuals 
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  #abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  abundance_df  <- read.csv(paste0("results/Abundance.corrected.",country,".csv")) %>%
    rename("obs.individuals"="individuals")%>%
    rename("individuals"="corrected.density")%>%
    aggregate(individuals ~ year + species + com_id, sum) %>% 
    mutate(individuals = individuals *areaplot) %>%
    spread(species, individuals) 
  
  
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
        dplyr::filter(get(neigh)>0) %>%
        dplyr::select(neigh) %>%
        unlist() %>%
        as.vector()
      
      seq.abun.neigh <- seq(min(quantile(neigh.abundance,probs=c(0.10,0.5,0.90), na.rm = T)),
                            max(quantile(neigh.abundance,probs=c(0.10,0.5,0.90), na.rm = T)),
                            (max(quantile(neigh.abundance,probs=c(0.10,0.5,0.90), 
                                          na.rm = T))-min(quantile(neigh.abundance,probs=c(0.10,0.5,0.90),
                                                                   na.rm = T)))/10)[1:10]
      seq.abun.neigh <-    seq.abun.neigh[which(   seq.abun.neigh>=0)]
      seq.abun.neigh <-    seq.abun.neigh[!is.na(seq.abun.neigh)]
      
      if(length( seq.abun.neigh )<10) next()
      
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
        
        #df_neigh_n$sigmoid[which(df_neigh_n$density==0)] <- 0
        
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
areaplot = 15
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  #abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  abundance_df  <- read.csv(paste0("results/Abundance.corrected.",country,".csv")) %>%
    rename("obs.individuals"="individuals")%>%
    rename("individuals"="corrected.density")%>%
    aggregate(individuals ~ year + species + com_id, sum) %>% 
    mutate(individuals = individuals *areaplot) %>%
    spread(species, individuals) 
  
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
areaplot = 15
country ="aus"
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] %>%
    #abundance_df  <- read.csv(paste0("results/Abundance.corrected.",country,".csv")) %>%
    #rename("obs.individuals"="individuals")%>%
    #rename("individuals"="corrected.density")%>%
    aggregate(individuals ~ year + species + com_id, sum) %>% 
    mutate(individuals = individuals *areaplot) %>%
    spread(species, individuals) 
  
  Realised.Int.Obs.focal<- list()
  Realised.Int.Obs.country.df <- NULL
  # Code.focal ="ARCA"
  for(Code.focal in Code.focal.list){ #focal.levels
    lambda = median(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_mean[,1], na.rm=T)
    df_alpha_generic_param = Parameters[[paste(country,"_",Code.focal)]]$df_alpha_generic_param
    
    year.levels <- colnames(Parameters[[paste(country,"_",Code.focal)]]$df_lambda_sd)
    print(paste0(country,Code.focal))
    
    abundance_short_focal_df <-  abundance_df %>%
      dplyr::filter(year ==y) %>%
      dplyr::filter(get(Code.focal) > 0) %>%
      dplyr::filter( !is.na(get(Code.focal))) %>%
      mutate_at(Code.focal,as.numeric) 
    head(abundance_short_focal_df)
    SpNames <- names(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt)
    
    test.sigmoid.all <- NULL
    test.sigmoid  <- NULL
    
    for( y in year.levels){
      for( neigh.sp in  SpNames){
        print(paste(neigh.sp,y))
        
        neigh.abundance <- abundance_short_focal_df %>%
          dplyr::filter(!is.na(get(neigh.sp))) %>%
          dplyr::select(year, com_id,neigh.sp) 
        
        if(nrow(neigh.abundance)==0) next
        
        names(neigh.abundance)[3] <-"density"
        
        alpha_initial = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                               neigh.sp]
        
        alpha_initial  <-  quantile(alpha_initial,probs=c(0.5)) 
        
        
        alpha_slope = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                             neigh.sp]
        
        alpha_slope  <- quantile(alpha_slope,probs=c(0.5)) 
        
        alpha_c = df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                         neigh.sp]
        
        alpha_c  <-    quantile(alpha_c,probs=c(0.5)) 
        
        N_opt = Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh.sp]
        
        N_opt  <-  quantile(N_opt,probs=c(0.5))
        
        param.neigh <- neigh.abundance %>%
          mutate(neigh = neigh.sp, 
                 country = country,
                 alpha_initial = alpha_initial,
                 alpha_slope = alpha_slope,
                 alpha_c=  alpha_c,
                 N_opt_mean = N_opt ,
                 focal=Code.focal) 
        
        for (n in 1:nrow(param.neigh)){
          #if(n==1){print(n)}
          df_neigh_n <- param.neigh[n,]
          
          df_neigh_n$sigmoid <- alpha_function4(df_neigh_n$alpha_initial,
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


load(paste0(project.dic,"results/Realised.Int.list.RData"))

#---- 1.4. Compute realised interactions for obs based on time  ----

source(paste0(home.dic,"code/PopProjection_toolbox.R"))
test.sigmoid.all  <- NULL
Realised.Obs.Year.list <- list()
areaplot = 15
for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  #abundance_df <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]] 
  abundance_df  <- read.csv(paste0("results/Abundance.corrected.",country,".csv")) %>%
    rename("obs.individuals"="individuals")%>%
    rename("individuals"="corrected.density")%>%
    aggregate(individuals ~ year + species + com_id, sum) %>% 
    mutate(individuals = individuals *areaplot) %>%
    spread(species, individuals) 
  
  Realised.Obs.Year.focal<- list()
  Realised.Obs.Year.country.df <- NULL
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
          dplyr::filter(year ==y) %>%
          dplyr::select(neigh.sp) %>%
          filter(!is.na(get(neigh.sp))) %>%
          dplyr::filter(get(neigh.sp)>0) %>%
          unique()
        
        minmax.vec <-quantile(neigh.abundance,probs=c(0.05,0.95), na.rm = T)
        
        neigh.abundance <- neigh.abundance %>%
          dplyr::filter(get(neigh.sp) > minmax.vec[1] & get(neigh.sp) < minmax.vec[2])
        N_opt_mean = as.numeric(unlist(median(Parameters[[paste(country,"_",Code.focal)]]$df_N_opt[,neigh.sp])))
        
        param.neigh <- neigh.abundance %>%
          rename("density"= as.character(neigh.sp)) %>%
          mutate(year=as.numeric(y),
                 neigh = neigh.sp, 
                 country = country,
                 alpha_initial = median(df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_initial"),
                                                               neigh.sp]),
                 alpha_slope =  median(df_alpha_generic_param[which(df_alpha_generic_param$parameter =="alpha_slope"),
                                                              neigh.sp]),
                 alpha_c =  median(df_alpha_generic_param[which(df_alpha_generic_param$parameter =="c"),
                                                          neigh.sp]),
                 N_opt_mean = N_opt_mean,
                 focal=Code.focal)
        
        for (n in 1:nrow(param.neigh)){
          #if(n==1){print(n)}
          param.neigh[n,"sigmoid"] <- alpha_function4(param.neigh$alpha_initial[n],
                                                      param.neigh$alpha_slope[n],
                                                      param.neigh$alpha_c[n],
                                                      param.neigh$density[n],
                                                      param.neigh$N_opt_mean[n])
          
          
        }
        test.sigmoid <- bind_rows(test.sigmoid,param.neigh)
      }
    }
    
    Realised.Obs.Year.country.df <- bind_rows(Realised.Obs.Year.country.df,test.sigmoid)
    
    save(Realised.Obs.Year.country.df,
         file=paste0(project.dic,"results/Realised.Obs.Year.df_",country,".RData"))
    
  }
  
  #sigmoid.list[[country]] <- sigmoid.list.focal
  Realised.Obs.Year.list[[country]] <-  Realised.Obs.Year.country.df
}

save(Realised.Obs.Year.list,
     file=paste0(project.dic,
                 "results/Realised.Obs.Year.list.RData"))
load(paste0(project.dic,
            "results/Realised.Obs.Year.list.RData"))
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3. Display Interactions  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#---- 3.1. Network fucnction ----
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
#figures/networks/Legends
#load(paste0(project.dic,"results/Realised.Int.list.RData"))
net.country <- list()
# country = "aus"
#species.focal = "GORO"

for(country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  ratio.mat <-   Realised.Int.list[[country]]  %>%
    group_by(focal, neigh) %>%
    summarise(count.positive = length(sigmoid[sigmoid > 0]),
              count.negative = length(sigmoid[sigmoid < 0]),
              count.neutral = length(sigmoid[sigmoid  ==0]),
              count.total = length(sigmoid)) %>%
    mutate(sigmoid = count.negative/(count.total-count.neutral)) %>%
    ungroup() %>%
    dplyr::select(-c('count.positive','count.negative','count.neutral','count.total')) %>%
    spread(neigh,sigmoid) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  
  strength.mat <- Realised.Int.list[[country]]  %>%
    aggregate(sigmoid ~ focal + neigh, function(x) abs(median(x))) %>%
    spread(neigh,sigmoid) %>%
    dplyr::select(-focal) %>%
    as.matrix()
  # for the network to have row as receiver
  if(!dim(strength.mat)[1]==dim(strength.mat)[2]){
    df.all.comb <- expand.grid(Code.focal.list,Code.focal.list) %>%
      rename("focal"="Var1") %>%
      rename("neigh"="Var2")
    
    ratio.mat <- Realised.Int.list[[country]]  %>%
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
    
    strength.mat <- Realised.Int.list[[country]]  %>%
      aggregate(sigmoid ~ focal + neigh, function(x) abs(median(x))) %>%
      right_join(df.all.comb) %>%
      spread(neigh,sigmoid) %>%
      column_to_rownames(var="focal") %>%
      as.matrix() %>%
      t()# for the network to have row has emitor and columns as receiver
  }
  plot.network.gradient.int(ratio.mat,strength.mat,"",0.0001)
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
# figures/metworks/Network.Realised.effect.pdf
ggsave( plot_grid(net.country$aus,
                  get_legend(plot.legend),
                  ncol = 1,
                  rel_heights =c(1,0.24),
                  labels = "",
                  hjust = 0, vjust = 1),
        file=paste0(home.dic,"figures/networks/Network.Obs.Sigmoid.effect_aus.pdf")) 
ggsave( plot_grid(net.country[["spain"]],
                  get_legend(plot.legend),
                  ncol = 1,
                  rel_heights =c(1,0.24),
                  labels = "",
                  hjust = 0, vjust = 1),
        file=paste0(home.dic,"figures/networks/Network.Obs.Realised.effect_spain.pdf")) 

#---- 3.3. Table sum up ----
str(Realised.Int.list)
sum.up.df <- NULL
for(country in country.list){
  
  sum.up.df.n <- Realised.Int.list[[country]] %>%
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
    dplyr::summarise(mean.effect = (mean(sigmoid)),
                     median.effect = (median(sigmoid)),
                     var.effect = (var(sigmoid)),
                     max.positive.effect = (max(sigmoid)),
                     max.negative.effect = (min(sigmoid)),
                     count.positive = length(sigmoid[sigmoid >0]),
                     count.negative = length(sigmoid[sigmoid  <0]),
                     count.total = length(sigmoid)) %>%
    dplyr::mutate(proportion.positive =(count.positive/count.total),
                  proportion.negative =(count.negative/count.total),
                  proportion.neutre =1-(proportion.positive +proportion.negative),
                  country = country,
                  effect ="given") %>%
    dplyr::rename("species" = "neigh")
  
  sum.up.focal.df.n <- Realised.Int.list[[country]] %>%
    group_by(focal) %>%
    dplyr::summarise(mean.effect = (mean(sigmoid)),
                     median.effect = (median(sigmoid)),
                     var.effect = (var(sigmoid)),
                     max.positive.effect = (max(sigmoid)),
                     max.negative.effect = (min(sigmoid)),
                     count.positive = length(sigmoid[sigmoid >0]),
                     count.negative = length(sigmoid[sigmoid  <0]),
                     count.total = length(sigmoid)) %>%
    dplyr::mutate(proportion.positive =(count.positive/count.total),
                  proportion.negative =(count.negative/count.total),
                  proportion.neutre =1-(proportion.positive +proportion.negative),
                  country = country,
                  effect ="received") %>%
    dplyr::rename("species" = "focal")
  
  
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
  
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.Looking at Interactions across time for answers----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2.1. Make data df ----
country = "spain"
Cool.trait.year.df <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.Year.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  Focal.group.df <- Realised.Obs.Year.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal +year, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.50 ~ "Receives Fac",
                                 sigmoid > 0.50 ~ "Receives Comp",
                                 T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid") 
  
  
  Neigh.group.df <- Realised.Obs.Year.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh +year , function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.50 ~ "Gives Fac",
                                 sigmoid > 0.50 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=trait.dist,
             scaled.trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list))) %>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>% 
                  rownames_to_column("focal") %>% rename("receiver.trait"=i) %>%
                  mutate(receiver.trait = scale(receiver.trait)))%>%
      left_join(trait.df %>% dplyr::select(all_of(i)) %>%
                  rownames_to_column("neigh") %>% rename("emitter.trait"=i)%>%
                  mutate(emitter.trait = scale(emitter.trait)))
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root mass density","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Stem height","Canopy shape") ~ "2.Aboveground",
                                         trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Water use efficiency","Canopy shape","Stem height","Leaf C to N ratio","Leaf area",
                                                      "Leaf nitrogen cc") ~ "2.Aboveground",
                                         trait %in% c("Mean fecundity") ~ "3.Reproduction"))
  }
  
  trait.detail.dist.df <- Realised.Obs.Year.list[[country]] %>%
    dplyr::filter(!is.na(density)) %>%
    left_join(specific.trait.dist,
              relationship ="many-to-many") 
  
  trait.dist.df <- Realised.Obs.Year.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    group_by(neigh,focal,year) %>%
    summarise(Sigm.Q5= median(sigmoid),
              Sigm.neg= length(sigmoid[sigmoid<0]),
              Sigm.total = length(sigmoid)#,# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(Sigm.ratio= Sigm.neg/Sigm.total) %>% 
    left_join(specific.trait.dist,
              relationship ="many-to-many") %>%
    mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                 Sigm.ratio< 0.5 ~"Fac",
                                 Sigm.ratio== 0.5 ~"Neutre"))
  
  sum.trait.dist.median.df <-  trait.dist.df %>%
    group_by(compORfac,trait,year) %>%
    mutate(Sigm.ratio = abs(Sigm.ratio -0.5)) %>%
    summarise(weighted.median = matrixStats::weightedMedian(scaled.trait.dist,Sigm.ratio,na.rm=T),
              weighted.mad = matrixStats::weightedMad(scaled.trait.dist,Sigm.ratio,na.rm=T)) %>%
    mutate(weightedQ1= weighted.median-weighted.mad,
           weightedQ9= weighted.median+weighted.mad)
  
  
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
    as.data.frame()
  
  
  Cool.trait.year.df[[country]] <- list(
    trait.detail.dist.df=trait.detail.dist.df,
    trait.dist.df=trait.dist.df,
    specific.trait.dist=specific.trait.dist,
    sum.trait.dist.median.df=sum.trait.dist.median.df,
    env_pdsi= env_pdsi)
  
}

#---- 2.2. Detailed plot of dist and year interaction ----
Cool.detailed.trait.year.plotlist <- list()
country="aus"
for( country in country.list){
  year.level <- nlevels(as.factor(env_pdsi$year))
  
  Precipitation.trait.test.df <- Cool.trait.year.df[[country]]$trait.detail.dist.df %>%
    mutate(year= as.numeric(year)) %>%
    left_join(Cool.trait.year.df[[country]]$env_pdsi %>%mutate(year= as.numeric(year))%>%
                dplyr::select(year,Precip.extrem))
  
  Test.trait.year.df <- NULL
  for( trait.i in levels(as.factor(Precipitation.trait.test.df$trait))){
    Precipitation.trait.i.df <- Precipitation.trait.test.df %>%
      dplyr::filter(trait ==trait.i)
    
    glm.precipitation.trait.i <-glm(sigmoid~ emitter.trait*Precip.extrem + receiver.trait*Precip.extrem + abs(scaled.trait.dist)*Precip.extrem, 
                                    family="gaussian",
                                    data=Precipitation.trait.i.df,
                                    na.action=na.omit)
    summary(glm.precipitation.trait.i)
    Test.trait.year.df.i <- as.data.frame(summary(glm.precipitation.trait.i)$coef) %>%
      rownames_to_column("parameters") %>%
      rename("estimate"="Estimate",
             "std.error"="Std. Error",
             "t.value"="t value",
             "pvalue"="Pr(>|t|)") %>%
      mutate(trait=trait.i,
             signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                              (pvalue<0.05 & pvalue >0.01 )~"**",
                              (pvalue<0.01)~"***",
                              T~""),
             RMSE = sqrt(mean(glm.precipitation.trait.i$residuals^2)))
    
    Test.trait.year.df <- bind_rows(Test.trait.year.df,Test.trait.year.df.i)
  }
  Test.trait.year.list[[country]] <- Test.trait.year.df
  levels(as.factor(Test.trait.year.df$parameters))
  sign.trait.df <- Test.trait.year.df %>% dplyr::filter(!signif =="") %>%
    dplyr::filter(parameters %in% c("Precip.extrem",
                                    "Precip.extrem:abs(scaled.trait.dist)",
                                    "Precip.extrem:receiver.trait",
                                    "emitter.trait:Precip.extrem")) %>%
    mutate(parameters=factor(parameters,
                             levels=c("Precip.extrem:receiver.trait",
                                      "emitter.trait:Precip.extrem",
                                      "Precip.extrem",
                                      "Precip.extrem:abs(scaled.trait.dist)"))) %>%
    mutate(n = as.numeric(as.factor(parameters))) %>%
    mutate(traitnum  = as.numeric(as.factor(trait))) %>%
    # Add scaling factor by n group 
    mutate(y_num = traitnum + case_when(n == 1 ~  0,
                                        n == 2 ~ -0.2,
                                        n == 3 ~   0.2,
                                        n == 4 ~   0.4)) 
  
  Cool.detailed.trait.year.plotlist[[paste0(country)]]<- sign.trait.df %>%
    ggplot(aes(y=y_num,
               x=estimate,
               shape=as.factor(parameters),
               color=as.factor(parameters),
               fill=as.factor(parameters))) +
    geom_pointrange(aes(xmin=estimate - std.error,
                        xmax=estimate + std.error),
                    size=1.5,alpha=0.7,
                    position=position_dodge(width=0.2)) +
    scale_y_continuous(breaks = unique(sign.trait.df$traitnum), 
                       labels = unique(sign.trait.df$trait),
                       limits=c(0,11.4)) +
    scale_color_colorblind("parameter")+
    scale_fill_colorblind("parameter")+
    scale_shape_manual("parameter",values=c(15,16,17,18,1,2,3))+
    theme_bw() +
    geom_vline(xintercept=0) + 
    labs(x="estimate",
         y=paste0("trait")) +
    guides(color = guide_legend(title.position = "top",
                                nrow=2)) +
    #annotate("text", x = -0.012, y=0.2, 
    #        label = "Increase in competition",
    #       size=4) + 
    #geom_segment(aes(x = 0, y = 0, xend = 0.015, yend = 0),
    #            arrow = arrow(length = unit(0.4, "cm")), color="black")+
    #annotate("text", x = 0.012, y=0.2, 
    #         label = "Increase in facilitation",
    #         size=4) + 
    ##geom_segment(aes(x = 0, y = 0, xend = -0.015, yend = 0),
    #             arrow = arrow(length = unit(0.4, "cm")),color="black") +
    theme(legend.position="bottom",
          legend.title =element_text(size=18),
          legend.text =element_text(size=16),
          axis.title=element_text(size=16),
          axis.text = element_text(size=14),
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.major.y = element_line(colour = 'black', linetype = 'dashed'),
          panel.grid.minor = element_blank())
  Cool.detailed.trait.year.plotlist[[paste0(country)]] 
  
}

Cool.detailed.trait.year.plotlist[[paste0("spain")]]
Cool.detailed.trait.year.plotlist[[paste0("aus")]]
ggarrange(Cool.detailed.trait.year.plotlist[[paste0("spain")]],
          Cool.detailed.trait.year.plotlist[[paste0("aus")]],
          nrow=1,legend="bottom",label.x = 0.2,
          align=c("h"),
          common.legend = T, labels=c("b.Spain","a.Autralia"),
          font.label = list(size = 20, color = "black", face = "bold", family = NULL))
#figures/GLM.Precip.traits.pdf
#---- 2.2. Q5 plot of dist and year interaction ----
Cool.detailed.trait.year.plotlist <- list()
country="aus"
for( country in country.list){
  year.level <- nlevels(as.factor(env_pdsi$year))
  
  Precipitation.trait.test.df <- Cool.trait.year.df[[country]]$trait.dist.df %>%
    dplyr::filter(!compORfac =="Neutre") %>%
    left_join(Cool.trait.year.df[[country]]$env_pdsi %>%mutate(year= as.numeric(year)))
  
  trait.dist.df <- Cool.trait.year.df[[country]]$trait.dist.df
  Test.trait.year.df <- NULL
  for( trait.i in levels(as.factor(trait.dist.df$trait))){
    Precipitation.trait.i.df <- Precipitation.trait.test.df %>%
      dplyr::filter(trait ==trait.i) %>%
      dplyr::filter(!is.na(trait.dist))
    
    glm.precipitation.trait.i <- glm(Sigm.Q5 ~ scaled.trait.dist +  scaled.trait.dist*Precip.extrem,
                                     data=Precipitation.trait.i.df)
    Test.trait.year.df.i <- as.data.frame(summary(glm.precipitation.trait.i)$coef) %>%
      rownames_to_column("parameters") %>%
      rename("estimate"="Estimate",
             "std.error"="Std. Error",
             "t.value"="t value",
             "pvalue"="Pr(>|t|)") %>%
      mutate(trait=trait.i,
             signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                              (pvalue<0.05 & pvalue >0.01 )~"**",
                              (pvalue<0.01)~"***",
                              T~""))
    Test.trait.year.df <- bind_rows(Test.trait.year.df,Test.trait.year.df.i)
  }
  Test.trait.year.list[[country]] <- Test.trait.year.df
  sign.trait.df <- Test.trait.year.df %>% dplyr::filter(!signif =="") %>%
    dplyr::filter(parameters %in% c("scaled.trait.dist","scaled.trait.dist:Precip.extrem","Precip.extrem"))
  
  Cool.detailed.trait.year.plotlist[[paste0(country)]] <- ggplot(data=Cool.trait.year.df[[country]]$trait.dist.df %>%
                                                                   dplyr::filter(!compORfac =="Neutre") %>%
                                                                   left_join(Cool.trait.year.df[[country]]$env_pdsi %>%mutate(year= as.numeric(year))),
                                                                 aes(x=emitter.trait,
                                                                     y=Sigm.Q5, #as.factor(year),
                                                                     group=as.factor(Precip.extrem),
                                                                     fill=Sigm.Q5)) +
    geom_point(shape=21, color="black",size=1, alpha=0.5,position = position_dodge(width = 0.2)) + 
    geom_smooth(aes(color=Precip.extrem),se=F,method="lm",alpha=0.5) +
    facet_wrap(.~trait , scale="free") + 
    labs(x="emitter trait",#x="Focal trait value - Neigh trait value",
         y="Species interactions",
         fill="Ratio of Competitive to \nFacilitative interactions")+
    scale_fill_gradientn(colours = rev(wes_palette("Zissou1", 
                                                   101, 
                                                   type = "continuous")))+
    scale_color_gradientn(colours = rev(viridis::viridis(year.level)))+
    theme_bw() 
  
  Cool.detailed.trait.year.plotlist[[paste0(country)]] 
  
}

Test.trait.year.list[["spain"]]
view(Test.trait.year.list[["aus"]])
Cool.detailed.trait.year.plotlist[[paste0("spain")]]
Cool.detailed.trait.year.plotlist[[paste0("aus")]]


#---- 2.3 Detailed plot of dist and year interaction ----
Cool.detailed.trait.year.plotlist <- list()
country="spain"
Test.trait.year.list <- list()
Cool.detailed.trait.precip.plotlist<-list()
for( country in country.list){
  
  Precipitation.trait.test.df <- Cool.trait.year.df[[country]]$trait.dist.df %>%
    dplyr::filter(!compORfac =="Neutre") %>%
    left_join(Cool.trait.year.df[[country]]$env_pdsi %>%mutate(year= as.numeric(year)))
  
  trait.dist.df <- Cool.trait.year.df[[country]]$trait.dist.df
  Test.trait.year.df <- NULL
  for( trait.i in levels(as.factor(trait.dist.df$trait))){
    Precipitation.trait.i.df <- Precipitation.trait.test.df %>%
      dplyr::filter(trait ==trait.i) %>%
      dplyr::filter(!is.na(trait.dist))
    
    library(plotly)
    plot_ly(Precipitation.trait.i.df,
            x=~scaled.trait.dist[,1],
            y=~Precip.extrem, 
            z=~Sigm.ratio,  
            type="scatter3d", mode="markers",color=~Sigm.ratio)%>% 
      layout(scene = list(
        xaxis = list(title = 'Focal trait value - neigh trait value'), 
        yaxis = list(title = 'Precipitation extrem (positive=wet)'),
        zaxis = list(title = 'Ratio of facilitation to comp')))
    
    glm.precipitation.trait.i <- glm(scaled.trait.dist~ compORfac + compORfac*Precip.extrem,
                                     data=Precipitation.trait.i.df)
    Test.trait.year.df.i <- as.data.frame(summary(glm.precipitation.trait.i)$coef) %>%
      rownames_to_column("parameters") %>%
      rename("estimate"="Estimate",
             "std.error"="Std. Error",
             "t.value"="t value",
             "pvalue"="Pr(>|t|)") %>%
      mutate(trait=trait.i,
             signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                              (pvalue<0.05 & pvalue >0.01 )~"**",
                              (pvalue<0.01)~"***",
                              T~""))
    Test.trait.year.df <- bind_rows(Test.trait.year.df,Test.trait.year.df.i)
  }
  Test.trait.year.list[[country]] <- Test.trait.year.df
  sign.trait.df <- Test.trait.year.df %>% dplyr::filter(!signif =="") %>%
    dplyr::filter(parameters %in% c("compORfacFac:Precip.extrem","Precip.extrem")) 
  
  Cool.detailed.trait.precip.plotlist[[paste0(country)]] <- Cool.trait.year.df[[country]]$trait.dist.df %>%
    dplyr::filter(!compORfac =="Neutre") %>%
    dplyr::filter(trait %in% sign.trait.df$trait) %>%
    left_join(Cool.trait.year.df[[country]]$env_pdsi %>%mutate(year= as.numeric(year))) %>%
    ggplot() +
    geom_point(aes(y=scaled.trait.dist,
                   x=Precip.extrem, #as.factor(year),
                   group=as.factor(compORfac ),
                   fill=Sigm.ratio),
               color="white",
               shape=21,
               size=3,
               alpha=0.5,
               position = position_dodge(width = 10)) + 
    geom_smooth(aes(y=scaled.trait.dist,
                    x=Precip.extrem, #as.factor(year),
                    group=as.factor(compORfac ),
                    color=compORfac),
                method="lm")+
    geom_label(data=sign.trait.df %>%
                 mutate(compORfac = case_when(parameters=="Precip.extrem"~"Comp",
                                              parameters=="compORfacFac:Precip.extrem"~"Fac")),
               aes(y= ifelse(estimate>0,1,-1),
                   label=paste0(format(estimate, scientific = TRUE)," ",signif),
                   #label=paste0(fancy_scientific(estimate)," ",signif),
                   color=compORfac),
               x=50,size=6) + 
    facet_wrap(.~trait , scale="free") + 
    labs(y="Focal trait value - Neigh trait value",
         x="Precipitation extrem",
         color="Median for each type of interaction",
         fill="Ratio of Competitive to Facilitative interactions")+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"),
                         limits=c(0,1),
                         breaks=c(0,1),
                         labels=c("100%\nFacilitatives",
                                  "100%\nCompetitives"))+
    scale_color_manual(values = c("#F21A00","#3B9AB2"))+
    guides(color="none") +
    theme_bw() +
    theme(legend.key.size = unit(1, 'cm'),
          legend.position="bottom",
          legend.title.position = "top",
          panel.grid.minor = element_blank(),
          legend.title =element_text(size=20),
          legend.text =element_text(size=20),
          axis.text=element_text(size=20),
          strip.text=element_text(size=20),
          axis.title=element_text(size=20),
          plot.margin = unit(c(1,0.5,0,0.5),"cm"))
  Cool.detailed.trait.precip.plotlist[[paste0(country)]] 
  
}
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  #l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
#view(Test.trait.year.list[["aus"]])
#view(Test.trait.year.list[["spain"]])

ggarrange(Cool.detailed.trait.precip.plotlist[[paste0("spain")]],
          Cool.detailed.trait.precip.plotlist[[paste0("aus")]],
          nrow=2,legend="bottom",heights=c(0.5,1),
          common.legend = T, labels=c("a. Spain","b. Australia"),
          font.label = list(size = 20, color = "black", face = "bold", family = NULL)) # figures/PCA.interaction.pdf

# figures/Traits.with.precip.pdf

# MATRIX MANTEL TEST 

summary(glm.inter.trait.i)
X <- trait.dist.df.i %>%
  dplyr::select(focal,neigh,emitter.trait) %>%
  #filter(!neigh=="PAIN" & !focal=="PAIN") %>% 
  spread(neigh,emitter.trait) %>%
  dplyr::select(-focal) %>% 
  as.matrix() 
diag(X) <- 0
X[is.na(X)] <- 0
Y <- trait.dist.df.i %>%
  dplyr::select(focal,neigh,theoretical.effect) %>%
  #filter(!neigh=="PAIN" & !focal=="PAIN") %>% 
  spread(neigh,theoretical.effect) %>%
  dplyr::select(-focal) %>% 
  as.matrix() 
diag(Y) <- 0
view(Y)
vegan::mantel(Y, X)
install.packages("ecodist")
library(ecodist)
ecodist::mantel(Y~X)
dist(trait.dist.df.i$theoretical.effect)
trait.dist.df.i <- trait.dist.df.i[complete.cases(trait.dist.df.i),]
mantel(theoretical.effect ~ receiver.trait, 
       data=trait.dist.df.i, nperm=10) 

dist.matrix.i <- MRM(dist(theoretical.effect) ~ dist(receiver.trait) + dist(emitter.trait) + dist(scaled.trait.dist), 
                     data=trait.dist.df.i, nperm=10,mrank=F)

dist.matrix.i
#summary(glm.inter.trait.i)
# GLM COMP
variable.to.keep <- Inter.trait.df.i$parameters[which(!Inter.trait.df.i$signif=="")]
variable.to.keep <- variable.to.keep[!variable.to.keep == "(Intercept)"]
if(!length(variable.to.keep) == 0){
  trait.dist.df.i <- trait.dist.df.i[complete.cases(trait.dist.df.i),]
  glm.inter.trait.zero.i <- glm(theoretical.effect ~ 1, trait.dist.df.i,
                                na.action = na.omit,
                                family="gaussian")
  
  glm.inter.trait.line.i <- glm(formula(paste("theoretical.effect ~ ",paste0(variable.to.keep, collapse="+"))), 
                                trait.dist.df.i,
                                na.action = na.omit,
                                family="gaussian")
  
  glm.inter.trait.quad.i <- glm(formula(paste0("theoretical.effect ~ ",paste0(variable.to.keep, collapse="+"),"+",
                                               paste0(paste0("I(",variable.to.keep,"^2)"), collapse="+"))),
                                trait.dist.df.i,
                                family="gaussian")
  glm.inter.trait.summary.i <- bind_rows(
    as.data.frame(lrtest( glm.inter.trait.zero.i,glm.inter.trait.line.i)) %>%
      mutate(model.comp =c("null","null vs linear"),
             model=c("glm.inter.trait.zero.i","glm.inter.trait.line.i")),
    as.data.frame(lrtest(  glm.inter.trait.zero.i,glm.inter.trait.quad.i))%>%
      mutate(model.comp =c("H0","null vs quadratic"),
             model=c("glm.inter.trait.zero.i","glm.inter.trait.quad.i")),
    as.data.frame(lrtest(  glm.inter.trait.line.i,glm.inter.trait.quad.i))%>%
      mutate(model.comp =c("H0","linear vs quadratic"),
             model=c("glm.inter.trait.line","glm.inter.trait.quad.i"))) %>%
    dplyr::filter(!model.comp=="H0") %>%
    rename("pvalue"="Pr(>Chisq)") %>%
    mutate(signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                            (pvalue<0.05 & pvalue >0.01 )~"**",
                            (pvalue<0.01)~"***",
                            T~" "),
           trait=trait.i,
           density.quantile=n)  
  glm.inter.trait.summary <- bind_rows( glm.inter.trait.summary,
                                        glm.inter.trait.summary.i)
  
  if(glm.inter.trait.summary.i$signif[which(glm.inter.trait.summary.i$model.comp=="null vs linear")]==" "){
    model.to.keep <- "glm.inter.trait.zero.i" 
  }else{model.to.keep <- "glm.inter.trait.i"
  }
  best.model = model.to.keep }else{
    model.to.keep <- "glm.inter.trait.i"
    best.model <- "no significant variables"
  }
Inter.trait.df.i <- as.data.frame(summary(get(model.to.keep))$coef) %>%
  rownames_to_column("parameters") %>%
  rename("estimate"="Estimate",
         "std.error"="Std. Error",
         "t.value"="t value",
         "pvalue"="Pr(>|t|)") %>%
  mutate(trait=trait.i,
         best.model =best.model,
         density.quantile=n,
         signif=case_when((pvalue<0.1 & pvalue >0.05 )~"*",
                          (pvalue<0.05 & pvalue >0.01 )~"**",
                          (pvalue<0.01)~"***",
                          T~""),
         RMSE = sqrt(mean(get(model.to.keep)$residuals^2)))
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.Looking at INTER Interactions----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 3.1. Make data df ----
country = "aus" 
Cool.trait.df <- list()
for( country in country.list){
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  year.levels <-levels(as.factor(Realised.Int.list[[country]]$year))
  col.df <- data.frame(color.name = unname(kelly())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  Focal.group.df <- Realised.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ focal, function(x) length(which(x<0))/length(x)) %>%
    mutate(focal.org = case_when(sigmoid < 0.50 ~ "Receives Fac",
                                 sigmoid > 0.50 ~ "Receives Comp",
                                 T ~ "Receives both"))%>%
    rename("sigmoid.focal" = "sigmoid")
  
  Neigh.group.df <- Realised.Int.list[[country]] %>%
    #dplyr::filter(realised.effect<10) %>%
    dplyr::filter(!neigh ==focal) %>%
    aggregate(sigmoid ~ neigh, function(x) length(which(x<0))/length(x)) %>%
    mutate(neigh.org = case_when(sigmoid < 0.50 ~ "Gives Fac",
                                 sigmoid > 0.50 ~ "Gives Comp",
                                 T ~ "Gives both")) %>%
    rename("sigmoid.neigh" = "sigmoid")
  
  
  
  if(country=="aus"){
    Neigh.group.df <-Neigh.group.df %>%
      mutate(neigh.org = case_when( sigmoid.neigh < 0.50 ~ "Gives Fac",
                                    sigmoid.neigh > 0.50 ~ "Gives Comp",
                                    T ~ "Gives both"))
    
    Focal.group.df <-   Focal.group.df %>%
      mutate(focal.org = case_when( sigmoid.focal< 0.5 ~ "Receives Fac",
                                    sigmoid.focal > 0.50 ~ "Receives Comp",
                                    T ~ "Receives both"))
  }
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]
  library(vegan)
  specific.trait.dist  <- NULL
  for( i in names(trait.df)){
    specific.trait.dist.n <- outer(trait.df[,i], trait.df[,i], '-') %>%
      as.data.frame() %>%
      gather(.,key="neigh",value="trait.dist") %>%
      mutate(trait.dist=trait.dist,
             scaled.trait.dist=scale(trait.dist),
             trait=i,
             neigh= rep(rownames(trait.df),each=length(Code.focal.list)),
             focal= rep(rownames(trait.df),times=length(Code.focal.list)))
    # specific.trait.dist.n <- as.data.frame(as.matrix(vegdist(trait.df %>% dplyr::select(i) %>% 
    #                                                       scale() %>% as.data.frame(),
    #                                                  na.rm = T,method="euclidean",diag=T))) %>%
    
    specific.trait.dist <- bind_rows(specific.trait.dist,specific.trait.dist.n)
  }
  
  if(country=="aus"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root tips","Root mass density","Root length",
                                            "Mean fecundity","C13 water use efficiency","Flower width","Seed mass",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","Root length","Root tips","Root mass density","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Stem height","Canopy shape") ~ "2.Aboveground",
                                         trait %in% c("Flower width","Mean fecundity","Seed mass")~ "3.Reproduction"))
  }
  if(country=="spain"){
    specific.trait.dist <- specific.trait.dist %>%
      mutate(trait = factor(trait, levels=c("SRL","Root diameter","Root mass density","SRA",
                                            "Mean fecundity","C13 water use efficiency","Leaf C to N ratio","Leaf area index",
                                            "Canopy shape","Stem height","SLA"))) %>%
      mutate(category.traits = case_when(trait %in% c("SRL","SRA","Root mass density","Root diameter","Leaf area index","Fast to slow spectrum") ~ "1.BelowGround",
                                         trait %in% c("SLA","Water use efficiency","Canopy shape","Stem height","Leaf C to N ratio","Leaf area",
                                                      "Leaf nitrogen cc") ~ "2.Aboveground",
                                         trait %in% c("Mean fecundity") ~ "3.Reproduction"))
  }
  trait.dist.df <- Realised.Int.list[[country]] %>%
    dplyr::filter(!neigh ==focal) %>%
    group_by(neigh,focal) %>%
    summarise(Sigm.Q5= median(sigmoid),
              Sigm.neg= length(sigmoid[sigmoid<0]),
              Sigm.total = length(sigmoid)#,# high mean high competition
    ) %>%
    ungroup() %>%
    mutate(Sigm.ratio= Sigm.neg/Sigm.total) %>% 
    left_join(specific.trait.dist %>% 
                left_join(Focal.group.df%>% dplyr::select(focal, focal.org,sigmoid.focal), multiple="all") %>%
                left_join(Neigh.group.df %>% dplyr::select(neigh, neigh.org,sigmoid.neigh), multiple="all"),
              relationship ="many-to-many") %>%
    mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                 Sigm.ratio< 0.5 ~"Fac",
                                 Sigm.ratio== 0.5 ~"Neutre"))
  
  sum.trait.dist.median.df <-  trait.dist.df %>%
    group_by(compORfac,trait) %>%
    mutate(Sigm.ratio = abs(Sigm.ratio -0.5)) %>%
    summarise(weighted.median = matrixStats::weightedMedian(scaled.trait.dist,Sigm.ratio,na.rm=T),
              weighted.mad = matrixStats::weightedMad(scaled.trait.dist,Sigm.ratio,na.rm=T),
              median.sigmoid=median(Sigm.Q5,na.rm=T),
              Q10.sigmoid=quantile(Sigm.Q5,c(0.10),na.rm=T),
              Q90.sigmoid=quantile(Sigm.Q5,c(0.90),na.rm=T)) %>%
    mutate(weightedQ1= weighted.median-weighted.mad,
           weightedQ9= weighted.median+weighted.mad)
  
  
  Cool.trait.df[[country]] <- list(
    trait.dist.df=trait.dist.df,
    specific.trait.dist=specific.trait.dist,
    sum.trait.dist.median.df=sum.trait.dist.median.df,
    Focal.group.df=Focal.group.df,
    Neigh.group.df=Neigh.group.df,
    Sym.group.df=Sym.group.df)
}

#---- 3.2. Detailed plot of dist and median interaction ----
Cool.detailed.trait.plotlist <- list()
for( country in country.list){
  
  Cool.detailed.trait.plotlist[[paste0(country)]] <- ggplot() +
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df,
               aes(y=scaled.trait.dist,
                   x=Sigm.Q5,
                   group=as.factor(compORfac),
                   fill=Sigm.ratio,
                   color=as.factor(compORfac)),
               shape=21,
               size=1,
               alpha=0.2) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df,
                    aes(y=weighted.median,
                        x=median.sigmoid,
                        xmin=Q10.sigmoid,
                        xmax=Q90.sigmoid,
                        color=as.factor(compORfac)),
                    shape=15,
                    size=1) + 
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df,
                    aes(y=weighted.median,
                        x=median.sigmoid,
                        ymin=weightedQ1,
                        ymax=weightedQ9,
                        color=as.factor(compORfac)),
                    shape=15,
                    size=1) +
    facet_wrap(.~trait,scale="free") + 
    labs(y="Focal trait value - Neigh trait value",
         x="Median sigmoid value given by Neigh")+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"))+
    scale_color_manual(values = c("#F21A00","#3B9AB2","#EBCC2A"))+
    theme_bw() 
  
}
Cool.detailed.trait.plotlist[[paste0("spain")]]
Cool.detailed.trait.plotlist[[paste0("aus")]]

#---- 3.3. Stat Test ----
Test.trait.list <- list()
for( country in country.list){
  trait.dist.df <- Cool.trait.df[[country]]$trait.dist.df
  sum.trait.dist.median.df <- Cool.trait.df[[country]]$sum.trait.dist.median.df
  Test.neigh.trait.df <- NULL
  Test.focal.trait.df <- NULL
  Test.trait.df <- NULL
  for( trait.i in levels(as.factor(trait.dist.df$trait))){
    trait.dist.df.i <-  trait.dist.df %>%
      dplyr::filter(trait==trait.i) %>%
      mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                   Sigm.ratio< 0.5 ~"Fac",
                                   Sigm.ratio== 0.5 ~"Neutre")) %>%
      dplyr::filter(!compORfac =="Neutre") %>%
      dplyr::select(trait,compORfac,trait.dist,Sigm.ratio)
    
    median.comp.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Comp" &
                                                                sum.trait.dist.median.df$trait==trait.i]
    median.fac.n <- sum.trait.dist.median.df$weighted.median[sum.trait.dist.median.df$compORfac=="Fac" &
                                                               sum.trait.dist.median.df$trait==trait.i]
    
    if(median.comp.n > median.fac.n){test.n ="greater"}else{test.n="less"}
    #test.n ="two.sided"
    test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$compORfac,
                        paired = F,
                        alternative = test.n)
    Test.trait.df.n <- data.frame(trait=trait.i,
                                  p.value = test$p.value,
                                  stat= test$statistic,
                                  alternative=test.n
    )
    Test.trait.df <- bind_rows(Test.trait.df,Test.trait.df.n)
    for( i in levels(as.factor(trait.dist.df$focal.org))){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(!focal.org==i )  %>%
        dplyr::select(trait,focal.org,trait.dist)
      if(nlevels(as.factor(trait.dist.df.i$focal.org))<2)next
      #https://www.r-bloggers.com/2020/06/wilcoxon-test-in-r-how-to-compare-2-groups-under-the-non-normality-assumption-2/
      test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$focal.org,
                          paired = F) # I am not sure if they should be paired or not 
      Test.focal.trait.df.n <- data.frame(trait=trait.i,
                                          focal.org.1 = levels(as.factor(trait.dist.df.i$focal.org))[1],
                                          focal.org.2= levels(as.factor(trait.dist.df.i$focal.org))[2],
                                          p.value = test$p.value,
                                          stat= test$statistic
      )
      Test.focal.trait.df <- bind_rows(Test.focal.trait.df,Test.focal.trait.df.n)
    }
    
    for( i in levels(as.factor(trait.dist.df$neigh.org))){
      trait.dist.df.i <-  trait.dist.df %>%
        dplyr::filter(trait==trait.i) %>%
        dplyr::filter(!neigh.org==i )  %>%
        dplyr::select(trait,neigh.org,trait.dist)
      if(nlevels(as.factor(trait.dist.df.i$neigh.org))<2)next
      test <- wilcox.test(trait.dist.df.i$trait.dist ~ trait.dist.df.i$neigh.org,
                          paired = F) # I am not sure if they should be paired or not 
      Test.neigh.trait.df.n <- data.frame(trait=trait.i,
                                          neigh.org.1 = levels(as.factor(trait.dist.df.i$neigh.org))[1],
                                          neigh.org.2= levels(as.factor(trait.dist.df.i$neigh.org))[2],
                                          p.value = test$p.value,
                                          stat= test$statistic)
      Test.neigh.trait.df <- bind_rows(Test.neigh.trait.df,Test.neigh.trait.df.n)
    }
  }
  Test.trait.list[[country]] <- list( Test.focal.trait.df= Test.focal.trait.df,
                                      Test.neigh.trait.df=Test.neigh.trait.df,
                                      Test.trait.df=Test.trait.df)
}
Test.trait.list[["spain"]]$Test.trait.df
Test.trait.list[["aus"]]$Test.trait.df
#---- 3.4. Summary plot with dist - rect ----
Cool.rect.trait.plotlist <- list()
country="aus"
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
for( country in country.list){
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  
  Test.trait.df <-   Test.trait.list[[country]]$Test.trait.df
  Cool.rect.trait.plotlist[[paste0(country,"_Pairwise_interactions")]] <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                              Sigm.ratio< 0.5 ~"Fac",
                                              Sigm.ratio== 0.5 ~"Neutre")) %>%
                 dplyr::filter(!compORfac =="Neutre") %>%
                 mutate(trait = factor(trait, levels = levels(specific.trait.dist$trait))),
               aes(x=scaled.trait.dist,  
                   y=trait,
                   #alpha=Sigm.ratio,
                   #fill=compORfac, # focal.org,
                   color=compORfac), # focal.org,
               position = position_dodge(width = 0.7),
               shape=16,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df %>%
                      dplyr::filter(!compORfac =="Neutre") %>%
                      mutate(trait = factor(trait, 
                                            levels = levels(specific.trait.dist$trait))),
                    position = position_dodge(width = 0.7),
                    aes(x=weighted.median, # median.trait.dist,
                        y=trait,
                        xmin=weightedQ1,#Q10.trait.dist,
                        xmax=weightedQ9,#Q90.trait.dist,
                        color=as.factor(compORfac)),
                    shape=16,
                    size=1) +
    scale_y_discrete(expand=c(0.1,0.1)) +
    annotate("text", x = 2, y=0.2, 
             label = "Focal has higher \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = 1.5, y = 0.2, xend = 2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm")))+
    annotate("text", x = -2, y=0.2, 
             label = "Focal has lower \n trait value than neigh",
             size=4) + 
    geom_segment(aes(x = -1.5, y = 0.2, xend = -2.5, yend = 0.2),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    labs(x="Focal trait value - Neigh trait value",
         y="trait",
         color="Pairwise interaction mainly")+#"Species grouping based on received interactions")+
    scale_color_manual(values = c("#F21A00","#3B9AB2"),
                       labels=c("Competitive","Facilitative"))+ # "#EBCC2A",
    scale_alpha_continuous(range=c(0.1,1)) +
    theme_bw()  +
    guides(color = guide_legend(title.position = "top")) + 
    theme(legend.position="bottom",
          legend.title =element_text(size=16),
          legend.text =element_text(size=14),
          axis.title=element_text(size=14),
          axis.text.x=element_text(color = "gray12", size = 16,
                                   face=ifelse(Test.trait.df$trait %in% Test.trait.df$trait[Test.trait.df$p.value<0.10],
                                               "bold","plain")))
  
}
Cool.rect.trait.plotlist[["aus_Pairwise_interactions"]]
Cool.rect.trait.plotlist[["spain_Pairwise_interactions"]]
#---- 3.5. Summary plot with dist - circ ----
Cool.circ.trait.plotlist <- list()
country="aus"
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
for( country in country.list){
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  
  Test.trait.df <-   Test.trait.list[[country]]$Test.trait.df %>%
    mutate(labels.pvalue = case_when((p.value <= 0.1 & p.value > 0.05 )~ "*",
                                     (p.value <= 0.05 & p.value > 0.01 )~ "**",
                                     p.value <= 0.01 ~ "***",
                                     T~ ""))
  
  Cool.circ.trait.plotlist[[paste0(country,"_Pairwise_interactions")]] <- ggplot()+
    geom_point(data=Cool.trait.df[[country]]$trait.dist.df %>%
                 mutate(compORfac = case_when(Sigm.ratio> 0.5 ~"Comp",
                                              Sigm.ratio< 0.5 ~"Fac",
                                              Sigm.ratio== 0.5 ~"Neutre")) %>%
                 dplyr::filter(!compORfac =="Neutre") %>%
                 mutate(trait = factor(trait, levels = levels(specific.trait.dist$trait))),
               #mutate(trait = factor(trait, levels = unique(specific.trait.dist$trait[order(specific.trait.dist$category.traits)]))),
               aes(y=scaled.trait.dist,  
                   x=trait,
                   #alpha=Sigm.ratio,
                   #fill=compORfac, # focal.org,
                   color=compORfac), # focal.org,
               position = position_dodge(width = 0.7),
               shape=16,size=2,alpha=0.2) +
    geom_pointrange(data=Cool.trait.df[[country]]$sum.trait.dist.median.df %>%
                      dplyr::filter(!compORfac =="Neutre"), 
                    position = position_dodge(width = 0.7),
                    aes(y=weighted.median, # median.trait.dist,
                        x=trait,
                        ymin=weightedQ1,#Q10.trait.dist,
                        ymax=weightedQ9,#Q90.trait.dist,
                        color=as.factor(compORfac)),
                    shape=16,
                    size=1) +
    #geom_text(data=Test.trait.df,aes(x=trait,y=3,label=labels.pvalue),size=12)+
    geom_hline(yintercept=0,color="black") +
    scale_y_continuous(breaks=c(0),expand = c(0, 0)) + 
    scale_x_discrete(labels=addline_format(paste(Test.trait.df$trait,Test.trait.df$labels.pvalue))) +
    labs(color="Pairwise interaction mainly")+
    scale_color_manual(values = c("#F21A00","#3B9AB2"),
                       labels=c("Competitive","Facilitative"))+ # "#EBCC2A",
    coord_polar() +
    annotate("text",x =3.4, y = 1.3, label = "Focal>Neigh",
             angle =-8,color = "gray12",size = 4.5)+ # spain: 14, aus: 6
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 3),
                 arrow = arrow(length = unit(0.5, "cm")))+
    annotate("text",x =3.64, y = -1.3, label = "Focal<Neigh",#"trait value",
             angle = -8,color = "gray12",size = 4.5)+
    geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = -3),
                 arrow = arrow(length = unit(0.5, "cm")))+
    guides(color = guide_legend(title.position = "top",
                                nrow=2),
           fill = guide_legend(title.position = "top",
                               nrow=2)) + 
    theme_bw() +
    theme(legend.position="bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          # Use gray text for the region names
          panel.background = element_blank(), #element_rect(fill = "white", color = "white"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title =element_text(size=20),
          legend.text =element_text(size=20),
          axis.text.x=element_text(color = "gray12", size = 16,
                                   face=ifelse(Test.trait.df$trait %in% Test.trait.df$trait[Test.trait.df$p.value<0.10],
                                               "bold","plain")),
          plot.margin = unit(c(0,0,0,0),"cm"))
  Cool.circ.trait.plotlist[[paste0(country,"_Pairwise_interactions")]]
  
}
Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]] #figures/Cool.circ.pairwise.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]] # figures/Cool.circ.pairwise.trait.aus.pdf

Cool.circ.trait.plotlist[["Pairwise_interactions"]] <-
  ggarrange(Cool.circ.trait.plotlist[[paste0("spain_Pairwise_interactions")]],
            Cool.circ.trait.plotlist[[paste0("aus_Pairwise_interactions")]],
            nrow=1,legend="bottom",
            common.legend = T, labels=c("a. Spain","b. Australia"),
            font.label = list(size = 20, color = "black", face = "bold", family = NULL))

Cool.circ.trait.plotlist["Pairwise_interactions"] # figures/Cool.circ.pairwise.trait.pdf

Cool.circ.trait.plotlist[[paste0("spain_circ")]] # figures/Cool.circ.trait.spain.pdf
Cool.circ.trait.plotlist[[paste0("aus_circ")]] # figures/Cool.circ.trait.aus.pdf



#---- 3.6. PCA----
# PCA analysis
#https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/116-mfa-multiple-factor-analysis-in-r-essentials/

country = "aus"
pca.plotlist <- list()
for( country in country.list){
  trait.dist.df <-  Cool.trait.df[[country]]$trait.dist.df 
  specific.trait.dist <- Cool.trait.df[[country]]$specific.trait.dist 
  Focal.group.df<- Cool.trait.df[[country]]$Focal.group.df
  Neigh.group.df<- Cool.trait.df[[country]]$Neigh.group.df
  
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]] 
  library(Hmisc)
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  res2<-rcorr(as.matrix(trait.df))
  #view(flattenCorrMatrix(res2$r, res2$P))
  # for spain, SLA is significantly correclated with mean fecundity,Root mass density, SRA, water supp
  # & SRL with Canopy shape
  write.csv(flattenCorrMatrix(res2$r, res2$P),
            file=paste0("results/correlation.matrix.",country,".csv"))
  names(trait.df)
  if(country=="spain"){
    to.keep.trait <- c("Root mass density","Leaf C to N ratio",#"SRA",
                       "SLA","SRL","Canopy shape","Leaf area index"
    )
  }else{
    to.keep.trait<- c("Root mass density",#"Root length",
                      "SRL","C13 water use efficiency")
    
  }
  trait.df <- get(paste0("clean.data.",country))[["plant_traits"]]  %>%
    dplyr::select(all_of(to.keep.trait))
  
  pca.trait.df <- MFA(trait.df, #%>%
                      #mutate(received =Focal.group.df$sigmoid.focal,
                      #given = Neigh.group.df$sigmoid.neigh),
                      group = c(rep(1, times= ncol(trait.df)+0)), 
                      type = c(rep("s", times= ncol(trait.df)+0)), # s = quantitative , n= factorial
                      name.group = c(colnames(trait.df)),#,"received","given"),
                      graph = T) 
  
  eig.val <- get_eigenvalue(pca.trait.df)
  head(eig.val)
  fviz_contrib(pca.trait.df, choice = "quanti.var", axes = 1, top = 20,palette = kelly(n=22))
  fviz_contrib(pca.trait.df, choice = "quanti.var", axes = 2, top = 20,palette = kelly(n=22))
  fviz_mfa_var(pca.trait.df, "quanti.var", palette = kelly(n=22), 
               col.var.sup = "violet", repel = TRUE)
  
  fviz_mfa_ind(pca.trait.df, col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE,visible = "quali.var")
  
  species.pca.df <- data.frame(xaxis=pca.trait.df$ind$coord[,1],
                               yaxis=pca.trait.df$ind$coord[,2]) %>%
    bind_cols(Focal.group.df) %>%
    bind_cols(Neigh.group.df) 
  trait.pca.df <- data.frame(xaxis=pca.trait.df$quanti.var$coord[,1],
                             yaxis=pca.trait.df$quanti.var$coord[,2]) %>%
    rownames_to_column("traitname")
  xperc <- eig.val[1,2]
  yperc <- eig.val[2,2]
  #if(country=="spain"){
  #species.pca.df <- data.frame(xaxis=-pca.trait.df$ind$coord[,2],
  #                             yaxis=pca.trait.df$ind$coord[,1]) %>%
  # bind_cols(Focal.group.df) %>%
  # bind_cols(Neigh.group.df) 
  #trait.pca.df <- data.frame(xaxis=-pca.trait.df$quanti.var$coord[,2],
  #                            yaxis=pca.trait.df$quanti.var$coord[,1]) %>%
  #  rownames_to_column("traitname")
  ## xperc <- eig.val[2,2]
  # yperc <- eig.val[1,2]
  #}
  
  dummy.df <- data.frame(xaxis=c(1.2,1.2),
                         yaxis=c(2.3,2),
                         focal=c("Received",
                                 "Given"))
  pca.plotlist[[country]] <- ggplot() +
    annotate(geom="point",x=dummy.df$xaxis,y=dummy.df$yaxis,
             shape=c(21,24),size=10,stroke=2) +
    annotate(geom="text",
             x=dummy.df$xaxis+0.5,y=dummy.df$yaxis,
             label=dummy.df$focal,size=6) +
    geom_point(data=species.pca.df,
               aes(x=xaxis-0.1,y=yaxis,fill=sigmoid.focal),
               shape=21,size=6) +
    geom_point(data=species.pca.df,
               aes(x=xaxis+0.1,y=yaxis,fill=sigmoid.neigh),
               shape=24,size=6) +
    geom_text(data=species.pca.df,
              aes(x=xaxis,y=yaxis+0.15,
                  label=addline_format(focal)),
              check_overlap = F,size=5)+
    scale_fill_gradientn(colours = wes_palette("Zissou1", 
                                               101, 
                                               type = "continuous"),
                         limits=c(0,1),
                         breaks=c(0,1),
                         labels=c("100%\nFacilitatives",
                                  "100%\nCompetitives")) +
    geom_segment(data=trait.pca.df,
                 aes(x=0,y=0,xend=xaxis,yend=yaxis),
                 arrow=arrow(length = unit(0.5, "cm"), 
                             type = "closed"))   +
    geom_text(data=trait.pca.df,
              aes(x=xaxis+0.05,y=ifelse(yaxis>0, yaxis+0.05,yaxis-0.05),
                  label=traitname),
              check_overlap = F,size=5)+
    labs(x=paste0("PCA dimension ",round(xperc,digits=1) ,"%"),
         y=paste0("PCA dimension ",round(yperc,digits=1) ,"%"),
         fill="Ratio of Competitive to \nFacilitative interactions")+
    theme_bw()+
    theme(legend.key.size = unit(1, 'cm'),
          legend.position="bottom",
          legend.title.position = "top",
          panel.grid.minor = element_blank(),
          legend.title =element_text(size=20),
          legend.text =element_text(size=20),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          plot.margin = unit(c(1,0.5,0,0.5),"cm"))
  pca.plotlist[[country]]
}
pca.plotlist[["spain"]]
ggarrange(pca.plotlist[["spain"]],
          pca.plotlist[["aus"]],
          nrow=1,legend="bottom",
          common.legend = T, labels=c("a. Spain","b. Australia"),
          font.label = list(size = 20, color = "black", face = "bold", family = NULL)) # figures/PCA.interaction.pdf


#---- 3. PLOTS----
#---- 3.1. Boxplot of RI  ----

for(country in country.list){
  Box.Plot.Realised.effect <- NULL
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
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



