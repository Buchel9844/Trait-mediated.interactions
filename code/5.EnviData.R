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



