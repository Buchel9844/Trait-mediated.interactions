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
library(scPDSI)
library(SPEI)

setwd("/home/lbuche/Eco_Bayesian/chapt3")
home.dic <- "" #"/Users/lisabuche/Documents/Projects/Facilitation_gradient/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Import envi data  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Good reference to understand drought index: https://nhess.copernicus.org/articles/23/3543/2023/
#---- SPAIN ----
# https://www.juntadeandalucia.es/agriculturaypesca/ifapa/riaweb/web/estacion/41/20#his
# SE20ETo mm/dia
# SE20ETo mm

env_spain <-  read.csv("data/spain_rawdata/env_spain.csv",sep=",") %>%
  separate("FECHA", sep="/", into=c("day","month","year")) %>%
  filter(year >= 2005 & year <= 2023 ) %>%
  select(year,month,Se20Precip,Se20ETo) %>%
  mutate_all(as.numeric) %>%
  aggregate(.~ year + month, sum) %>%
  mutate(year.month = factor(paste0(year,sep="_",month)))
    

view(env_spain )
line_spain_env <- ggplot(env_spain) +
  geom_line(aes(y=Se20ETo,x=factor(year.month),group=1),
            color="black",size=2) +
  geom_line(aes(y=Se20Precip,x=factor(year.month),group=1),
            color="grey",size=2)+
  scale_x_discrete(breaks=c("2001_1","2002_1","2003_1","2004_1",
                            "2005_1","2006_1","2007_1","2008_1",
                            "2009_1","2010_1","2011_1",
                            "2012_1","2013_1","2014_1","2015_1","2016_1",
                            "2017_1","2018_1","2019_1", "2020_1","2021_1",
                            "2022_1", "2023_1","2024_1"),
                   labels=c("2001","2002","2003","2004",
                            "2005","2006","2007","2008",
                            "2009","2010","2011",
                            "2012","2013","2014",
                            "2015","2016","2017","2018",
                            "2019","2020","2021",
                            "2022","2023","2024")) +
  labs(title="Montly Precipitation (grey) and ETo (black)",x="",
       y="mm/month",
       fill="conditions") +
  theme_bw()+
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
         axis.title.x= element_text(size=20),
         axis.title.y= element_text(size=20),
         title=element_text(size=16))
line_spain_env

# compute monthly conventional Palmer Drought Severity Index (PDSI)
# https://www.rdocumentation.org/packages/scPDSI/versions/0.1.3/topics/pdsi
#install.packages('/Users/lisabuche/Downloads/scPDSI_0.1.2.tar.gz', 
 #                repo=NULL,
 #                dependencies=TRUE)

library(scPDSI)

P <- aggregate(Se20Precip ~ month +year, env_spain,sum)[,"Se20Precip"]
PE <- aggregate(Se20ETo  ~ month + year, env_spain,sum)[,"Se20ETo"]
sc_pdsi <- pdsi(P, PE, start = 2005,end=2023,sc=F,AWC = 50)
plot(sc_pdsi)
length(sc_pdsi$X)

spain_env_pdsi <- data.frame(spain_pdsi = as.numeric(sc_pdsi$X),
           year=rep(c(2005:2023),each=12),
           month= rep(c(1:12),times=19)) %>%
  mutate(year.month = factor(paste0(year,sep="_",month),
                             levels=paste0(year,sep="_",month))) %>%
  left_join(env_spain) %>%
  rename("prec"="Se20Precip")


view(spain_env_pdsi)
write_csv(spain_env_pdsi,
          "results/spain_env_pdsi.csv")

# VISUALISATION 

boxplot_spain_env_pdsi <- spain_env_pdsi %>%
  filter(month >2 & month < 9) %>%
  filter(!is.na( spain_pdsi)) %>%
  mutate(fill.pdsi = ifelse(0 >= spain_pdsi, "dry","wet")) %>%
  ggplot(aes(y=spain_pdsi)) +
  geom_vline(xintercept = c("2005","2006","2007","2008",
                            "2009","2010","2011",
                            "2012","2013","2014",
                            "2015","2016","2017","2018",
                            "2019","2020","2021",
                            "2022","2023","2024"),
             color="grey",alpha=0.7) +
  geom_boxplot(aes(x=as.factor(year),fill=fill.pdsi)) +
  geom_hline(yintercept=0, color="black") +
  labs(y="PDSI",title="Average from March to August",x="",
       fill="conditions") +
  theme_bw()+
  scale_fill_manual(values=c("wet" = "blue", 
                             "dry" = "orange")) +
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
         axis.title.x= element_text(size=20),
         axis.title.y= element_text(size=20),
         title=element_text(size=16))
boxplot_spain_env_pdsi

hist_spain_env_pdsi <- spain_env_pdsi %>%
  filter(!is.na( spain_pdsi)) %>%
  mutate(fill.pdsi = ifelse(0 >= spain_pdsi, "dry","wet")) %>%
  ggplot() +
  geom_vline(xintercept = c("2005_1","2006_1","2007_1","2008_1",
                            "2009_1","2010_1","2011_1",
                            "2012_1","2013_1","2014_1","2015_1","2016_1",
                            "2017_1","2018_1","2019_1", "2020_1","2021_1",
                            "2022_1", "2023_1","2024_1"),
             color="grey",alpha=0.7) +
  geom_bar(stat="identity",
           aes(y=spain_pdsi,
               x=as.factor(year.month),
               fill=fill.pdsi)) +
  scale_fill_manual(values=c("wet" = "blue", 
                             "dry" = "orange")) +
  scale_x_discrete(breaks=c("2005_1","2006_1","2007_1","2008_1",
                            "2009_1","2010_1","2011_1",
                            "2012_1","2013_1","2014_1","2015_1","2016_1",
                            "2017_1","2018_1","2019_1", "2020_1","2021_1",
                            "2022_1", "2023_1","2024_1"),
                   labels=c("2005","2006","2007","2008",
                            "2009","2010","2011",
                            "2012","2013","2014",
                            "2015","2016","2017","2018",
                            "2019","2020","2021",
                            "2022","2023","2024")) +
   labs(y="PDSI",x="",title="Monthly value per year",
       fill="conditions") +
  geom_hline(yintercept=0, color="black") +
  theme_bw() +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_blank(),
         strip.text = element_text(size=28),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=20, angle=66, hjust=1),
         axis.text.y= element_text(size=20),
         axis.title.x= element_text(size=20),
         axis.title.y= element_text(size=20),
         title=element_text(size=16))
hist_spain_env_pdsi

spain_env_plot <- list(line_spain_env=line_spain_env,
                     boxplot_spain_env_pdsi=boxplot_spain_env_pdsi,
                     hist_spain_env_pdsi=hist_spain_env_pdsi)
plot_spain_env_pdsi <- ggarrange( plotlist =spain_env_plot,
                                labels=c("a.","b.","c."),
                                align = c( "v"),
          label.x = 0,
          label.y = 1,
          font.label = list(size = 24,
                            color = "black",
                            face = "bold", 
                            family = NULL),
          nrow=3,
          common.legend = T,
          legend="bottom")
plot_spain_env_pdsi
ggsave(plot_spain_env_pdsi,
       file ="figures/plot_spain_env_pdsi.pdf")

#---- AUS ----
#http://www.bom.gov.au/watl/eto/tables/wa/daily.shtml
# Monthly Precipitation and ETO were averaged over the following weather stations:
station_aus <- c("Bowgada", "Five.Gums", "High.Fields", "Latham", "Carnamah", "Wubin", "Dalwallinu.North", "Perenjori")
# Average across all station for daily rain - to reduce uncertainty
# sum of that average for every month
env_prec_aus <- read.csv("data/aus_rawdata/monthly_rainfall.csv",sep=",")

head(env_prec_aus)
env_prec_aus <- read.csv("data/aus_rawdata/Perenjori_monthly_prec.csv",sep=",")  %>%
  mutate(station="Perenjori") %>%
  bind_rows(read.csv("data/aus_rawdata/Latham_monthly_prec.csv",sep=",")%>%
              mutate(station="Latham")) %>%
  gather(Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec,
         key="month",value="prec") %>%
  mutate(month = case_when(month=="Jan"~1,
                           month=="Feb"~2,
                           month=="Mar"~3,
                           month=="Apr"~4,
                           month=="May"~5,
                           month=="Jun"~6,
                           month=="Jul"~7,
                           month=="Aug"~8,
                           month=="Sep"~9,
                           month=="Oct"~10,
                           month=="Nov"~11,
                           month=="Dec"~12)) %>%
  mutate(prec=as.numeric(prec)) %>%
  aggregate(prec ~ month + year, function(x){mean(x,na.rm=T)})
view(env_prec_aus)
# and ETo via TERN :
# https://tern-landscapes.earthengine.app/view/cmrset-landsat-v22 for all stations
env_ETo_aus <- read.csv("data/aus_rawdata/ETo_aus.csv",sep=",") %>%
  mutate(month = case_when(month=="Jan"~1,
                           month=="Feb"~2,
                           month=="Mar"~3,
                           month=="Apr"~4,
                           month=="May"~5,
                           month=="Jun"~6,
                           month=="Jul"~7,
                           month=="Aug"~8,
                           month=="Sep"~9,
                           month=="Oct"~10,
                           month=="Nov"~11,
                           month=="Dec"~12))
env_ETo_aus$sumETo <- rowMeans(env_ETo_aus[,station_aus])
head(env_ETo_aus)
# merge prec and ETo data set
env_aus <- left_join(env_ETo_aus,env_prec_aus, by=c("month", "year")) %>%
  filter(year >= 2004 & year < 2024 ) %>%
  mutate(year.month = factor(paste0(year,sep="_",month),
                             levels=paste0(year,sep="_",month)))

head(env_aus)
line_aus_env <- ggplot(env_aus) +
  geom_line(aes(y=sumETo,x=factor(year.month),group=1),
           color="black",size=2) +
  geom_line(aes(y=prec,x=factor(year.month),group=1),
           color="grey",size=2)+
  scale_x_discrete(breaks=c("2001_1","2002_1","2003_1","2004_1",
                            "2005_1","2006_1","2007_1","2008_1",
                            "2009_1","2010_1","2011_1",
                            "2012_1","2013_1","2014_1","2015_1","2016_1",
                            "2017_1","2018_1","2019_1", "2020_1","2021_1",
                            "2022_1", "2023_1","2024_1"),
                   labels=c("2001","2002","2003","2004",
                            "2005","2006","2007","2008",
                            "2009","2010","2011",
                            "2012","2013","2014",
                            "2015","2016","2017","2018",
                            "2019","2020","2021",
                            "2022","2023","2024")) +
  labs(title="Montly Precipitation (grey) and ETo (black)",x="year",
       y="mm/month",x="",fill="conditions") +
  theme_bw()+
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
         axis.title.x= element_text(size=20),
         axis.title.y= element_text(size=20),
         title=element_text(size=16))
line_aus_env
# compute SDI
#https://cran.rstudio.com/web/packages/drought/drought.pdf
install.packages("drought")
library(drought)
X= env_aus[,"prec"] # 10-year monthly data
Yc <- ACCU(X,ts=6) # Compute the 6 month accumulated series
fit1 <- SDI(X,ts=6) # Get the standardized drought index (or SPI)
z = matrix(t(fit1$SDI),ncol=1)
Res <- RunDS(z, -1)# Get drought duration and severity based on threshold SPI=-1
Y = env_aus[,"sumETo"] # 10-year monthly data
fit2<-MSDI(X,Y,ts=6) # Compute the 6 month Multivariate Standardized Drought Index (MSDI)
fit2$MSDI #Get the empirical MSDI

#compute PDSI
P_aus <- env_aus[,"prec"]
PE_aus <-  env_aus[,"sumETo"]
# if sc= F then  p = 0.897 and q = 1/3 according to Palmer (1965) 
# we don't have enough data for the pdsi to calibrate correctly 
sc_pdsi_aus <- pdsi(P=P_aus, PE=PE_aus, AWC = 50,
                    start = 2004,sc = F)
plot(sc_pdsi_aus, index = "PHDI")
plot(sc_pdsi_aus, index = "WPLM")

length(sc_pdsi_aus$X)

aus_env_pdsi <- data.frame(aus_pdsi = as.numeric(sc_pdsi_aus$X),
                             year=rep(c(2004:2023),each=12),
                             month= rep(c(1:12),times=20)) %>%
  mutate(year.month = factor(paste0(year,sep="_",month),
                             levels=paste0(year,sep="_",month))) %>%
  left_join(env_aus %>%
              dplyr::select(month,year,sumETo,prec)%>%
              aggregate(.~ year + month,sum))

write_csv(aus_env_pdsi,
          "results/aus_env_pdsi.csv")

# VISUALISATION 

boxplot_aus_env_pdsi <- aus_env_pdsi %>%
  filter(month > 6 & month < 11) %>%
  mutate(fill.pdsi = ifelse(0 >= aus_pdsi, "dry","wet")) %>%
  ggplot(aes(y=aus_pdsi)) +
  geom_vline(xintercept = c("2001","2002","2003","2004",
                            "2005","2006","2007","2008",
                            "2009","2010","2011",
                            "2012","2013","2014",
                            "2015","2016","2017","2018",
                            "2019","2020","2021",
                            "2022","2023","2024"),
             color="grey",alpha=0.7) +
  geom_boxplot(aes(x=as.factor(year),fill=fill.pdsi)) +
  geom_hline(yintercept=0, color="black") +
  labs(y="PDSI",title="Average from June to October",x="",
       fill="conditions") +
  theme_bw()+
  scale_fill_manual(values=c("wet" = "blue", 
                             "dry" = "orange")) +
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
         axis.title.x= element_text(size=20),
         axis.title.y= element_text(size=20),
         title=element_text(size=16))
boxplot_aus_env_pdsi

hist_aus_env_pdsi <- aus_env_pdsi %>%
  left_join(env_aus) %>%
  mutate(fill.pdsi = ifelse(0 >= aus_pdsi, "dry","wet")) %>%
  ggplot() +
  geom_vline(xintercept = c("2001_1","2002_1","2003_1","2004_1",
                            "2005_1","2006_1","2007_1","2008_1",
                            "2009_1","2010_1","2011_1",
                            "2012_1","2013_1","2014_1","2015_1","2016_1",
                            "2017_1","2018_1","2019_1", "2020_1","2021_1",
                            "2022_1", "2023_1","2024_1"),
             color="grey",alpha=0.7) +
  geom_bar(stat="identity",
           aes(y=aus_pdsi,
               x=as.factor(year.month),
               fill=fill.pdsi)) +
  #geom_point(aes(y=prec/10,
   #             x=as.factor(year.month))) +
  scale_fill_manual(values=c("wet" = "blue", 
                             "dry" = "orange")) +
  scale_x_discrete(breaks=c("2001_1","2002_1","2003_1","2004_1",
                            "2005_1","2006_1","2007_1","2008_1",
                            "2009_1","2010_1","2011_1",
                            "2012_1","2013_1","2014_1","2015_1","2016_1",
                            "2017_1","2018_1","2019_1", "2020_1","2021_1",
                            "2022_1", "2023_1","2024_1"),
                   labels=c("2001","2002","2003","2004",
                            "2005","2006","2007","2008",
                            "2009","2010","2011",
                            "2012","2013","2014",
                            "2015","2016","2017","2018",
                            "2019","2020","2021",
                            "2022","2023","2024")) +
  #scale_y_continuous(
 #   name = "PSDI",
  #  sec.axis = sec_axis(~.*10, name="Precipitation"),
  #  limits = c(-10,10)) +
  labs(y="PDSI",title="Monthly value per year",x="",
       fill="conditions") +
  geom_hline(yintercept=0, color="black") +
  theme_bw() +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_blank(),
         strip.text = element_text(size=28),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=20, angle=66, hjust=1),
         axis.text.y= element_text(size=20),
         axis.title.x= element_text(size=20),
         axis.title.y= element_text(size=20),
         title=element_text(size=16))
hist_aus_env_pdsi

aus_env_plot <- list(line_aus_env=line_aus_env,
                     boxplot_aus_env_pdsi=boxplot_aus_env_pdsi,
                     hist_aus_env_pdsi=hist_aus_env_pdsi)
plot_aus_env_pdsi <- ggarrange( plotlist =aus_env_plot,
                                 labels=c("a.","b.","c."),
                                 align = c( "v"),
                                 label.x = 0,
                                 label.y = 1,
                                 font.label = list(size = 24,
                                                   color = "black",
                                                   face = "bold", 
                                                   family = NULL),
                                 nrow=3,
                                 common.legend = T,
                                 legend="bottom")
plot_aus_env_pdsi
ggsave(plot_aus_env_pdsi,
       file="figures/plot_aus_env_pdsi.pdf")

#----old----
#---- 1.3 Precipitation extrem over time ----
Precipitation.plot.list <- list()
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
    mutate(Precip.extrem = Precip - mean(Precip),
           Precip.extrem.Q1 = Precip - (mean(Precip) + sd(Precip)),
           Precip.extrem.Q9 = Precip - (mean(Precip) + sd(Precip))) %>%
    rename("year"="year.thermique") %>%
    mutate(PDSI.max = max(PDSI.mean)) %>%
    mutate_at("PDSI.mean",function(x) (x/ max(abs(x)))) %>% # scaling 
    as.data.frame() %>%
    mutate(year=as.numeric(year))
  
  
  Precipitation.plot.list[[country]] <- env_pdsi %>%
    ggplot(aes(y=Precip.extrem,x=year))+
    geom_line(size=2) +
    theme_bw() +
    labs(y="Precipitation \n extreme (mm)")+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_y_continuous(breaks=c(-100,0,100)) +
    scale_x_continuous(breaks=c(2011,2013,2015,2017,2019,2021,2023))+ 
    coord_cartesian( expand = F, default = FALSE, clip = "on") +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=28),
           legend.text=element_text(size=20),
           legend.title=element_text(size=20),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_text(size=18, angle=66, hjust=1),
           axis.text.y= element_text(size=18),
           axis.title.x= element_blank(),
           axis.title.y= element_text(size=18),
           title=element_text(size=18),
           plot.margin=unit(c(1,1,0,0),"cm"))
}

ggarrange(ggarrange(abundance_plotlist[["aus"]],
                    Precipitation.plot.list[["aus"]],
                    nrow=2,align="v",
                    heights = c(1.5,1)),
          ggarrange(abundance_plotlist[["spain"]],
                    Precipitation.plot.list[["spain"]],
                    nrow=2,align="v",
                    heights = c(1.5,1)),
          nrow=2,align="v",
          heights = c(1,1),
          label.x = c(-0.06,-.04),
          label.y = 1.01,
          font.label = list(size = 26, color = "black", 
                            face = "bold", family = NULL),
          ncol=1,labels=c("a. Australia","b. Spain"))

# figures/main/Abundance.over.time.pdf
#height =862
#width =953
#---- 1.4 Corrected abundances over time with precipitation extrem in it ----
widthplot = 25*25
area.plot = pi*7.5^2
abundance_corrected_plotlist <- NULL

for(country in country.list){
  abundance_summarycorrected  <- read.csv(paste0("results/Abundance.corrected.",country,".csv"))
  
  upper <- log10(max(abundance_summarycorrected$corrected.density.at.plot)) /max(abundance_summarycorrected$Precip.extrem)
  upper.precip <- max(max(abundance_summarycorrected$Precip.extrem),
                      abs(min(abundance_summarycorrected$Precip.extrem)))
  if(country=="spain"){
    legend.position.vec <- c(0.8, 0.85)
  }else{legend.position.vec <- c(0.8, 0.15)}
  abundance_corrected_plotlist[[country]] <- abundance_summarycorrected %>%
    ggplot(aes(x=year)) +
    stat_summary(aes(x=year,
                     y = Precip.extrem),
                 fun.y = median,
                 color="black",size=3,
                 geom = "line") + 
    stat_summary(aes(y = log10(corrected.density.at.plot)/ upper,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = median,
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=2,alpha=0.8,
                 position=position_dodge2(width=0.5, reverse=T)) +
    stat_summary(aes(y = log10(corrected.density.at.plot) / upper,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = median,alpha=0.5,
                 geom = "line",size=1,
                 position=position_dodge2(width=0.5, reverse=T)) +
    scale_y_continuous("Precipitation extrem",
                       position = "right", 
                       #limits = c(-upper.precip,upper.precip),
                       sec.axis = sec_axis(~10^(.* upper),
                                           name="Median number of \nindividuals in 25x25cm plot",
                                           breaks= c(0.1,1,10,100))) +
    scale_x_continuous("") +
    scale_color_manual(values= col.df$color.name) +
    labs(color="",#title=paste0("Density over time of annual plants in ",country)
    ) +
    coord_cartesian( xlim = NULL, #ylim=c(0,750),
                     expand = TRUE, default = FALSE, clip = "on") +
    #scale_color_manual(values=safe_colorblind_palette) +
    theme_bw() +
    guides(
      colour = guide_legend(position = "inside",ncol=3))+
    theme( legend.key.size = unit(0.8, 'cm'),
           legend.position.inside =  legend.position.vec,
           legend.background = element_rect(fill="NA"),
           legend.box.background = element_rect(
             fill = alpha('grey98', .6), 
             size = 0.1, linetype = "solid", color = "#333333"
           ),
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=28),
           legend.text=element_text(size=16),
           legend.title=element_blank(),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_text(size=20, angle=66, hjust=1),
           axis.text.y= element_text(size=20),
           axis.title.x= element_blank(),
           axis.title.y= element_text(size=24),
           title=element_text(size=16),
           plot.margin=unit(c(1,0,0,0),"cm"))
  abundance_corrected_plotlist[[country]] 
  breakfun <- function(x) {
    10^scales::extended_breaks()(log10(x))
  }
  
}
abundance_corrected_plotlist[["spain"]] # figures/Abundance_spain.pdf
abundance_corrected_plotlist[["aus"]] #figures/Abundance_aus.pdf

ggarrange(abundance_corrected_plotlist[["aus"]],
          abundance_corrected_plotlist[["spain"]],
          nrow=2, common.legend = F,
          font.label = list(size = 30, color = "black",
                            face = "bold", family = NULL),
          label.x = c(-0.06,-0.05),
          label.y = 1.02,
          labels=c("a. Australia", 
                   "b. Spain"))
#figures/Abundances.pdsi.pdf
#---- 1.5 Overall abundance over time ----
widthplot = 25*25
ALLabundance_corrected_plotlist <- NULL
country ="aus"
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
    mutate(individuals=individuals * widthplot)
  #group_by(year)
  ggplot( abundance_summary, aes(x=individuals)) + geom_density()+
    facet_wrap(year~.,scale="free")
  
  # the coefficient if the sampling  effort differences between the years
  Correction.model.year <-  as.data.frame(confint(glmmTMB(individuals ~ year + (1|species), 
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
    mutate(indiviuals.at.plot = individuals*widthplot,
           corrected.density = case_when(individuals > 0 ~ ((individuals*widthplot)-Coef.year),
                                         T ~ 0)) %>%
    left_join(env_pdsi)
  
  write.csv(abundance_summarycorrected,
            file=paste0("results/Abundance.corrected.",country,".csv"))
  upper <- log10(max(abundance_summarycorrected$corrected.density)) /max(abundance_summarycorrected$Precip.extrem)
  upper.precip <- max(max(abundance_summarycorrected$Precip.extrem),
                      abs(min(abundance_summarycorrected$Precip.extrem)))
  if(country=="spain"){
    legend.position.vec <- c(0.8, 0.85)
  }else{legend.position.vec <- c(0.8, 0.15)}
  ALLabundance_corrected_plotlist[[country]] <- abundance_summarycorrected %>%
    ggplot(aes(x=year)) +
    stat_summary(aes(x=year,
                     y = Precip.extrem),
                 fun.y = median,
                 color="black",size=3,
                 geom = "line") + 
    stat_summary(aes(y = log10(corrected.density)/ upper),
                 color="red",
                 fun.y = median,
                 fun.ymin = function(x) quantile(x,0.05), 
                 fun.ymax = function(x) quantile(x,0.95), 
                 geom = "pointrange",size=2,alpha=0.8,
                 position=position_dodge2(width=0.5, reverse=T)) +
    stat_summary(aes(y = log10(corrected.density) / upper),
                 color="red",
                 fun.y = median,alpha=0.5,
                 geom = "line",size=1,
                 position=position_dodge2(width=0.5, reverse=T)) +
    scale_y_continuous("Precipitation extrem",
                       position = "right", 
                       #limits = c(-upper.precip,upper.precip),
                       sec.axis = sec_axis(~10^(.* upper),
                                           name="Median number of \nindividuals in 25x25cm plot",
                                           breaks= c(0.1,1,10,100))) +
    scale_x_continuous("") +
    #scale_color_manual(values=safe_colorblind_palette) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(size=28),
          legend.text=element_text(size=16),
          legend.title=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.text.x= element_text(size=20, angle=66, hjust=1),
          axis.text.y= element_text(size=20),
          axis.title.x= element_blank(),
          axis.title.y= element_text(size=24),
          title=element_text(size=16),
          plot.margin=unit(c(1,0,0,0),"cm"))
  ALLabundance_corrected_plotlist[[country]] 
  breakfun <- function(x) {
    10^scales::extended_breaks()(log10(x))
  }
  
}
ALLabundance_corrected_plotlist[["spain"]] # figures/Abundance_spain.pdf
ALLabundance_corrected_plotlist[["aus"]] #figures/Abundance_aus.pdf

ggarrange(ALLabundance_corrected_plotlist[["aus"]],
          ALLabundance_corrected_plotlist[["spain"]],
          nrow=2, common.legend = F,
          font.label = list(size = 30, color = "black",
                            face = "bold", family = NULL),
          label.x = c(-0.06,-0.05),
          label.y = 1.02,
          labels=c("a. Australia", 
                   "b. Spain"))
#figures/Abundances.pdsi.pdf
#---- 1.6 Corrected abundances of each species over time ----
widthplot = 25*25
Spabundance_plotlist <- NULL

for(country in country.list){
  abundance_summary <- read.csv(paste0("results/Abundance.corrected.",country,".csv"))
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
  Spabundance_plotlist[[country]] <- abundance_summary %>%
    gather(individuals,corrected.density, key="corrected.or.not.density",value="density") %>%
    ggplot() +
    geom_density_ridges2(aes(x=density*widthplot,
                             y=corrected.or.not.density,
                             fill=species),alpha=0.5, scale=1) + 
    facet_wrap(.~species,scales="free_x",nrow=3) + 
    scale_color_manual(values= col.df$color.name) +
    scale_fill_manual(values= col.df$color.name) + 
    scale_y_discrete(labels=c("Corrected","Observed")) + 
    theme_bw() +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "none",
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
           title=element_text(size=16),
           plot.margin=unit(c(1,0,0,0),"cm"))
}
Spabundance_plotlist[["spain"]] #figures/SpAbundance_spain.pdf
Spabundance_plotlist[["aus"]]  #figures/SpAbundance_aus.pdf

#---- 1.7 Abundances over PDSI ----
widthplot = 25
abundance_pdsi_plotlist <- NULL
for(country in country.list){
  
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
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
  
  abundance_pdsi_plotlist[[country]] <- abundance_summary %>%
    mutate(year = as.numeric(year))%>%
    left_join(env_pdsi) %>%
    ggplot() +
    stat_summary(aes(x=Precip.extrem , y = individuals*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = median,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange",size=2) +
    stat_summary(aes(x=Precip.extrem , y = individuals*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 geom = "line",size=1) +
    scale_y_log10() +
    scale_color_manual(values= col.df$color.name) +
    labs(x="Precipitation extrem (positive=wet)",
         color="species",y="Median number of \nindividuals in 25x25cm plot"#,
         #title=paste0("Density across PDSI of annual plants in ",country)
    ) +
    coord_cartesian( xlim = NULL, #ylim = c(0,200),
                     expand = TRUE, default = FALSE, clip = "on") +
    #scale_color_manual(values=safe_colorblind_palette) +
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
}
abundance_pdsi_plotlist[["spain"]] #figures/Abundance_pdsi_spain.pdf
abundance_pdsi_plotlist[["aus"]] #figures/Abundance_pdsi_aus.pdf

#---- 1.8 Corrected abundances over PDSI ----
widthplot = 25
abundance_corrected_pdsi_plotlist <- NULL
for(country in country.list){
  
  abundance_summary <- read.csv(paste0("results/Abundance.corrected.",country,".csv"))
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  
  year.levels <-levels(as.factor( abundance_summary$year))
  col.df <- data.frame(color.name = unname(kelly.colors())[3:(length(Code.focal.list)+2)],
                       neigh = Code.focal.list)
  
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
  
  abundance_corrected_pdsi_plotlist[[country]] <- abundance_summary %>%
    mutate(year = as.numeric(year))%>%
    left_join(env_pdsi) %>%
    ggplot() +
    stat_summary(aes(x=Precip.extrem , y = corrected.density*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange",size=2) +
    stat_summary(aes(x=Precip.extrem , y = corrected.density*widthplot,
                     group=as.factor(species),
                     color=as.factor(species)),
                 fun.y = mean,
                 geom = "line",size=1) +
    scale_y_log10() +
    scale_color_manual(values= col.df$color.name) +
    labs(x="Precipitation extrem (positive=wet)",
         color="species",y="Median number of \nindividuals in 25x25cm plot"#,
         #title=paste0("Density across PDSI of annual plants in ",country)
    ) +
    coord_cartesian( xlim = NULL, #ylim = c(0,200),
                     expand = TRUE, default = FALSE, clip = "on") +
    #scale_color_manual(values=safe_colorblind_palette) +
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
}
abundance_corrected_pdsi_plotlist[["spain"]] #figures/Abundance_pdsi_spain.pdf
abundance_corrected_pdsi_plotlist[["aus"]] #figures/Abundance_pdsi_aus.pdf
ggarrange(abundance_pdsi_plotlist[["spain"]],
          abundance_corrected_pdsi_plotlist[["spain"]],
          nrow=2, common.legend = T,legend="bottom",
          font.label = list(size = 20, color = "black",
                            face = "bold", family = NULL),
          label.x = 0,
          label.y = 1,
          labels=c("a. Observed raw density", 
                   "b. Observed desnity for year/sampling effort effect"))


ggarrange(abundance_pdsi_plotlist[["aus"]],
          abundance_corrected_pdsi_plotlist[["aus"]],
          nrow=2, common.legend = T,legend="bottom",
          font.label = list(size = 20, color = "black",
                            face = "bold", family = NULL),
          label.x = 0,
          label.y = 1,
          labels=c("a. Observed raw density", 
                   "b. Observed desnity for year/sampling effort effect"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#---- 2. Environmental Gradient  ----
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
# PDSI
env_pdsi_aus <- read.csv(paste0(home.dic,"results/aus_env_pdsi.csv")) 
env_pdsi_spain <- read.csv(paste0(home.dic,"results/spain_env_pdsi.csv")) 

#---- 2.1. Median ----
env_pdsi_plotlist <- NULL
env_pdsi_dflist <- NULL
for(country in country.list){
  abundance_summary <- get(paste0("clean.data.",country))[[paste0("abundance_",country,".summary")]]
  Code.focal.list <- get(paste0("clean.data.",country))[[paste0("species_",country)]]
  
  year.levels <-levels(as.factor( abundance_summary$year))
  
  
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
  env_pdsi_dflist[[country]] <-env_pdsi
  env_pdsi_plotlist[[country]] <-env_pdsi %>%
    gather(PDSI.mean,PDSI.sd,PDSI.Q1,PDSI.Q9,PDSI.Q5,
           key="PDSI",value="value") %>%
    ggplot(aes(x=year,y=value,group=PDSI,
               color=PDSI))+
    annotate("rect",xmin=min(env_pdsi$year),
             xmax=max(env_pdsi$year),
             ymin=0,ymax=max(env_pdsi$PDSI.Q9)+0.1,
             fill="lightblue",
             alpha=0.2) +
    geom_line(aes(linetype=PDSI)) +
    scale_color_manual(values=c("black","grey","grey","grey","red")) +
    scale_linetype_manual(values=c("solid","dashed","solid","dashed","solid"))+
    coord_cartesian(expand = FALSE) + 
    theme_bw()
  
}
env_pdsi_plotlist[["aus"]] #figures/env_pdsi_aus.pdf
env_pdsi_plotlist[["spain"]] #figures/env_pdsi_spain.pdf


env_pdsi %>%
  ggplot(aes(x=PDSI.sd,y=PDSI.sd))+
  geom_line() +
  geom_point()

env_pdsi_dflist[["aus"]] %>%
  ggplot(aes(y=PDSI.sd,x=PDSI.mean))+
  geom_line() +
  geom_point()
env_pdsi_dflist[["spain"]] %>%
  ggplot(aes(y=PDSI.sd,x=PDSI.mean))+
  geom_line() +
  geom_point()


