#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 1. SET UP: Import packages----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#install.packages("lavaan", dependencies=TRUE)
library(lavaan)
#install.packages("emmeans")
library(emmeans)
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/piecewiseSEM/piecewiseSEM_2.1.0.tar.gz"
#install.packages(packageurl, repos=NULL, type="source",dependencies = F)
library(piecewiseSEM)

home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/chapt3/"
home.dic <- "/home/lbuche/Eco_Bayesian/chapt3/"
#---- 1.2. Import clean data----
load(file=paste0(home.dic,"data/clean.data.aus.RData"))
load(file=paste0(home.dic,"data/clean.data.spain.RData"))
country.list <- c("aus","spain")
#---- 1.3. Import results---
# parameter from models
load(paste0(home.dic,"results/Parameters_alpha.RData")) 
# Raw sigmoid
Param.sigm.df <- list()
Param.sigm.df.aus <- read.csv(paste0(project.dic,"results/Param.sigmoid.aus.csv.gz"))
Param.sigm.df.spain <- read.csv(paste0(project.dic,"results/Param.sigmoid.spain.csv.gz"))
Param.sigm.df$spain <- Param.sigm.df.spain
Param.sigm.df$aus<- Param.sigm.df.aus
# realised interactions across time
load(paste0(project.dic,"results/Realised.Int.list.RData"))

# realised interactions for each year
load(paste0(project.dic,"results/Realised.Int.Year.list.RData")) 

# PDSI
env_pdsi_aus <- read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) 
env_pdsi_spain <- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv")) 
for( country in "aus"){
}
trait.gr <- c("pol","above","below")
Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
for(i in 1:3){
  grouping.list[[i]] <- get(paste0("clean.data.",country))[[paste0(country,"_",trait.gr[i],"_grouping")]] %>%
    dplyr::select(final.code,group.cluster,name.cluster) 
}
Realised.Int.Year.sem <- Realised.Int.Year.list[[country]] %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, median),
            multiple = "all") %>%
  left_join(get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]])
str( Realised.Int.Year.sem)

m1b <-   '
  # regressions
    realised.effect ~ 1 + density + aus_pdsi + pol.traits + below.traits + above.traits
    density ~ aus_pdsi + pol.traits + below.traits + above.traits
  # variance (optional)
    realised.effect ~~ realised.effect
'
library(lavaan)
fit1b <- sem(m1b, data=Realised.Int.Year.sem)
summary(fit1b)


#SEM
# Hanlun CODE ====
plot.psemhl <- function(x, return=FALSE,
                        node_attrs = data.frame(shape = "circle", color = "black",
                                                fillcolor = "white"),
                        edge_attrs = data.frame(style = "solid", color="black"),
                        ns_dashed = T, alpha=0.05,
                        show = "std", digits = 2, layout="tree",
                        add_edge_label_spaces = T, correlation=F,weight=4,nodewidth=0.65,nodeheight=0.5,...
){
  library(DiagrammeR)
  #get the coefficients table
  ctab <- coefs(x)
  ctab$Response <- as.character(ctab$Response)
  ctab$Predictor <- as.character(ctab$Predictor)
  if(!correlation){
    er.na <- which(ctab$Std.Error=="-")
    char <- sapply(gregexpr("~~",ctab$Response),function(x)x[1])
    char.na<- which(char==1)
    ctab <- ctab[-unique(er.na,char.na),]
  }
  #make a nodes DF
  unique_nodes <- unique(c(ctab$Response, ctab$Predictor))
  nodes <- create_node_df(n = length(unique_nodes),
                          nodes = unique_nodes,
                          type = "lower",
                          label = unique_nodes,
                          width = nodewidth,
                          height = nodeheight)
  nodes <- cbind(nodes, node_attrs)
  nodes[] <- lapply(nodes, as.character)
  nodes$id <- as.numeric(nodes$id)
  #make an edges DF
  edges <- create_edge_df(
    from = match(ctab$Predictor, unique_nodes),
    to = match(ctab$Response, unique_nodes))
  edges <- data.frame(edges, edge_attrs)
  edges[] <- lapply(edges, as.character)
  edges$id <- as.numeric(edges$id)
  edges$from <- as.numeric(edges$from)
  edges$to <- as.numeric(edges$to)
  if(ns_dashed) edges$style[which(ctab$P.Value>alpha)] <- "dashed"
  if(show == "std") edges$label = round(ctab$`Std.Estimate`, digits)
  if(show == "unstd") edges$label = round(ctab$Estimate, digits)
  if(add_edge_label_spaces) edges$label = paste0(" ", edges$label, " ")
  edges$color[as.numeric(edges$label)<0]<-"#CC0000"
  edges$color[as.numeric(edges$label)>0]<-"navy"
  edges$fontcolor<- edges$color
  #turn into a graph
  sem_graph <- create_graph(nodes, edges, directed=TRUE)%>%
    select_edges() %>%
    set_edge_attrs_ws(
      edge_attr = value,
      value = weight*(abs(as.numeric(edges$label))))%>%
    copy_edge_attrs(
      edge_attr_from = value,
      edge_attr_to = penwidth)
  if(return) return(sem_graph)else
    render_graph(sem_graph, layout=layout,...)
}

#All variables should be numeric
rownames(Realised.Int.Year.sem) <- 1:nrow(Realised.Int.Year.sem)

# realised effect
Realised.Int.Year.sem <- Realised.Int.Year.list[[country]] %>%
  #filter(realised.effect < 100) %>%
  mutate(realised.effect = log(realised.effect)) %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, mean)%>%
              mutate(year.before= year-1) %>%
              left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
                          aggregate(aus_pdsi ~ year, mean)%>%
                          rename("year.before"= "year",
                                 "aus_pdsi.tminus1" = "aus_pdsi"))%>%
              mutate(aus_pdsi.dif = aus_pdsi - aus_pdsi.tminus1),
            multiple = "all") %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, mean)%>%
              rename("aus_pdsi.mean" ="aus_pdsi"),
            multiple = "all") %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, var)%>%
              rename("aus_pdsi.var" ="aus_pdsi"),
            multiple = "all") %>%
  left_join(get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]])


sem.realised.effect <- psem(
  glm(realised.effect ~ 1 + density   + pol.traits + below.traits + above.traits,
      "gaussian",Realised.Int.Year.sem),
  glm(density ~ aus_pdsi.mean  + below.traits + pol.traits + above.traits, "gaussian",
      Realised.Int.Year.sem),
  #glm(aus_pdsi.var ~ aus_pdsi.mean , "gaussian",
  #     Realised.Int.Year.sem),
  Realised.Int.Year.sem
)
AIC(sem.realised.effect)
summary(sem.realised.effect)
rsquared(sem.realised.effect, method = "mcfadden")
plot(sem.realised.effect)
plot.psemhl(sem.realised.effect, correlation = T, layout = "tree")

# Ratio 
Realised.Int.Year.Sem.Ratio <- Realised.Int.Year.list[[country]] %>%
  filter(realised.effect < 100) %>%
  group_by(neigh,focal,year)%>%
  summarise(RE.neg= length(realised.effect[realised.effect<1]),
            RE.total = length(realised.effect)) %>%
  ungroup() %>%
  mutate(RE.ratio = RE.neg/RE.total) %>%
  as.data.frame() %>%
  left_join(Realised.Int.Year.list[[country]] %>%
              aggregate(density ~ year + focal + neigh, mean)%>%
              rename("density.mean" ="density")) %>%
  left_join(Realised.Int.Year.list[[country]] %>%
              aggregate(density ~ year + focal + neigh, var)%>%
              rename("density.var" ="density")) %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, mean) %>%
              rename("aus_pdsi.mean" ="aus_pdsi"),
            multiple = "all") %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, var)%>%
              rename("aus_pdsi.var" ="aus_pdsi"),
            multiple = "all") %>%
  left_join(get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]])
str(Realised.Int.Year.Sem.Ratio)
sem.ratio <- psem(
  glm(RE.ratio ~ 1 + density.mean +  density.var + below.traits +
        above.traits ,
      "gaussian",Realised.Int.Year.Sem.Ratio),
  glm(density.mean ~  aus_pdsi.mean  , "gaussian",
      Realised.Int.Year.Sem.Ratio),
  glm(density.var ~ density.mean +aus_pdsi.mean, "gaussian",
      Realised.Int.Year.Sem.Ratio),
  #glm(aus_pdsi.var ~ aus_pdsi.mean, "gaussian",
  #    Realised.Int.Year.Sem.Ratio),
  Realised.Int.Year.Sem.Ratio
)
AIC(sem.ratio )
plot.psemhl(sem.ratio , correlation = T, layout = "tree")
summary(sem.ratio )
rsquared(sem.ratio , method = "mcfadden")
plot(sem.ratio )


# Median
Realised.Int.Year.Sem.Median <- Realised.Int.Year.list[[country]] %>%
  filter(realised.effect < 100) %>%
  group_by(neigh,focal,year)%>%
  summarise(RE.0.5= abs(median(realised.effect,na.rm=T)-1)) %>%
  ungroup() %>%
  as.data.frame() %>%
  left_join(Realised.Int.Year.list[[country]] %>%
              aggregate(density ~ year + focal + neigh, mean)%>%
              rename("density.mean" ="density")) %>%
  left_join(Realised.Int.Year.list[[country]] %>%
              aggregate(density ~ year + focal + neigh, var)%>%
              rename("density.var" ="density")) %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, mean) %>%
              rename("aus_pdsi.mean" ="aus_pdsi"),
            multiple = "all") %>%
  left_join(read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) %>%
              aggregate(aus_pdsi ~ year, var)%>%
              rename("aus_pdsi.var" ="aus_pdsi"),
            multiple = "all") %>%
  left_join(get(paste0("clean.data.",country))[[paste0("trait.dist_",country,".df")]])
str(Realised.Int.Year.Sem.Ratio)
sem.median  <- psem(
  glm(RE.0.5 ~ 1 + density.mean + density.var,
      "gaussian",Realised.Int.Year.Sem.Median ),
  #MASS::glm.nb(density ~ aus_pdsi + pol.traits,
  #Realised.Int.Year.Sem.Ratio),
  glm(density.mean ~  aus_pdsi.mean  , "gaussian",
      Realised.Int.Year.Sem.Ratio),
  glm(density.var ~ density.mean +aus_pdsi.mean, "gaussian",
      Realised.Int.Year.Sem.Ratio),
  Realised.Int.Year.Sem.Median 
)

AIC(sem.median )
plot.psemhl(sem.median , correlation = T, layout = "tree")
summary(sem.median )
rsquared(sem.median , method = "mcfadden")
plot(sem.median )

