install.packages("scPDSI")
install.packages("SPEI")

# Import envi data

env_spain <- read.csv("data/env_spain.csv",sep=",")

env_spain <- env_spain %>%
  separate("FECHA", sep="/", into=c("day","month","year")) %>%
  mutate_all(as.numeric)  %>%
  filter(year > 2014 & year < 2022)
str(env_spain) 

levels(as.factor(env_spain$year))

# compute monthly conventional Palmer Drought Severity Index (PDSI)
# https://www.rdocumentation.org/packages/scPDSI/versions/0.1.3/topics/pdsi
#install.packages('/Users/lisabuche/Downloads/scPDSI_0.1.2.tar.gz', 
 #                repo=NULL,
 #                dependencies=TRUE)

library(scPDSI)

P <- aggregate(Se20Precip ~ month +year, env_spain,sum)[,"Se20Precip"]
PE <- aggregate(Se20ETo  ~ month+year, env_spain,sum)[,"Se20ETo"]
sc_pdsi <- pdsi(P, PE, start = 2015)


spain_env_pdsi <- data.frame(spain_pdsi = as.numeric(sc_pdsi$X),
           year=rep(c(2015:2021),each=12),
           month= rep(c(1:12),times=7)) %>%
  mutate(year.month = factor(paste0(year,sep="_",month),
                             levels=paste0(year,sep="_",month)))

write_csv(spain_env_pdsi,
          "data/spain_env_pdsi.pdf")

# VISUALISATION
boxplot_spain_env_pdsi <- spain_env_pdsi %>%
  ggplot(aes(y=spain_pdsi)) +
  geom_boxplot(aes(x=as.factor(year))) +
  geom_hline(yintercept=0, color="black") +
  labs(y="PDSI",x="") +
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
         axis.title.x= element_text(size=24),
         axis.title.y= element_text(size=24),
         title=element_text(size=16))
boxplot_spain_env_pdsi

hist_spain_env_pdsi <- spain_env_pdsi %>%
  mutate(fill.pdsi = ifelse(0 > spain_pdsi, "dry","wet")) %>%
  ggplot() +
  geom_bar(stat="identity",
           aes(y=spain_pdsi,
               x=as.factor(year.month),
               fill=fill.pdsi)) +
  scale_fill_manual(values=c("wet" = "blue", "dry" = "orange")) +
  scale_x_discrete(breaks=c("2015_6","2016_6","2017_6","2018_6",
                              "2019_6","2020_6","2021_6"),
                   labels=c("2015","2016","2017","2018",
                            "2019","2020","2021")) +
  geom_vline(xintercept = c("2015_1","2016_1",
             "2017_1","2018_1","2019_1", "2020_1","2021_1"),
             color="grey",alpha=0.7) +
  labs(y="PDSI",x="Monthly value per year",fill="conditions") +
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
         axis.title.x= element_text(size=24),
         axis.title.y= element_text(size=24),
         title=element_text(size=16))
hist_spain_env_pdsi


plot_spain_env_pdsi <- ggarrange(boxplot_spain_env_pdsi,
          hist_spain_env_pdsi,
          labels=c("a.","b."),
          align = c( "v"),
          label.x = 0,
          label.y = 1,
          font.label = list(size = 24,
                            color = "black",
                            face = "bold", 
                            family = NULL),
          nrow=2)
plot_spain_env_pdsi
ggsave(plot_spain_env_pdsi,
       "figures/plot_spain_env_pdsi.pdf")


