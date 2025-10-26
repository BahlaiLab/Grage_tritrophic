## Kathryn Grage
## January 30th
## Tritrophic analysis figures and seeing what exists

## first figures. need to display the slopes of tritrophic level one and two but over what?
## years? window lengths?

## establish a gradient for amount of positive and amount of negative

## lets just use hubbard brook for now, keep it simple

## read in our csv file
getwd()
hubbard_data <- read.csv("model_output/model_hubbard_brook.csv")
library(tidyverse)
hubbard_3 <- hubbard_data %>%
  filter(N_years == "3") %>% 
  ggplot(aes(x=start_year, y=slope))+
  geom_point(na.rm=TRUE) +
  ggtitle("Hubbard Three Years")+
  scale_x_continuous(limits=c(1969,2016))+
  xlab("starting year")+
  facet_wrap(~trophic_level) 
hubbard_3

hubbard_4 <- hubbard_data %>%
  filter(N_years == "4") %>% 
  ggplot(aes(x=start_year, y=slope))+
  geom_point(na.rm=TRUE) +
  ggtitle("Hubbard Four Years")+
  scale_x_continuous(limits=c(1969,2016))+
  xlab("starting year")+
  facet_wrap(~trophic_level) 
hubbard_4

hubbard_5 <- hubbard_data %>%
  filter(N_years == "5") %>% 
  ggplot(aes(x=start_year, y=slope))+
  geom_point(na.rm=TRUE) +
  ggtitle("Hubbard Five Years")+
  scale_x_continuous(limits=c(1969,2016))+
  xlab("starting year")+
  facet_wrap(~trophic_level) 
hubbard_5

hubbard_6 <- hubbard_data %>%
  filter(N_years == "6") %>% 
  ggplot(aes(x=start_year, y=slope))+
  geom_point(na.rm=TRUE) +
  ggtitle("Hubbard Six Years")+
  scale_x_continuous(limits=c(1969,2016))+
  xlab("starting year")+
  facet_wrap(~trophic_level) 
hubbard_6

hubbard_7 <- hubbard_data %>%
  filter(N_years == "7") %>% 
  ggplot(aes(x=start_year, y=slope))+
  geom_point(na.rm=TRUE) +
  ggtitle("Hubbard Seven Years")+
  scale_x_continuous(limits=c(1969,2016))+
  xlab("starting year")+
  facet_wrap(~trophic_level) 
hubbard_7

hubbard_8 <- hubbard_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=slope))+
  geom_point(na.rm=TRUE) +
  ggtitle("Hubbard Eight Years")+
  scale_x_continuous(limits=c(1969,2016))+
  xlab("starting year")+
  facet_wrap(~trophic_level) 
hubbard_8

## establish a heat map gradient for slope values
## psuedocode:
## average values to per year using a pipe
## include a standard deviation to that!
## question: how do you include a standard deviation to a heat map?

library(tidyverse)
library(viridisLite)
library(ggpubr)
hubbard_3_heat <- hubbard_data %>%
  filter(N_years == "3") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_gradient(magma) +  
  guides(fill=guide_legend(title="Slope")) +
  labs(title = "Three Year Interval - Hubbard Brook",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
hubbard_3_heat

hubbard_4_heat <- hubbard_data %>%
  filter(N_years == "4") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_gradient(magma) +  
  guides(fill=guide_legend(title="Slope")) +
  labs(title = "Four Year Interval - Hubbard Brook",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal() 
hubbard_4_heat

hubbard_5_heat <- hubbard_data %>%
  filter(N_years == "5") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_gradient(magma) +  
  guides(fill=guide_legend(title="Slope")) +
  labs(title = "Five Year Interval - Hubbard Brook",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal() 
hubbard_5_heat

hubbard_6_heat <- hubbard_data %>%
  filter(N_years == "6") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_gradient(magma) +  
  guides(fill=guide_legend(title="Slope")) +
  labs(title = "Six Year Interval - Hubbard Brook",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal() 
hubbard_6_heat

hubbard_7_heat <- hubbard_data %>%
  filter(N_years == "7") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_gradient(magma) +  
  guides(fill=guide_legend(title="Slope")) +
  labs(title = "Seven Year Interval - Hubbard Brook",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal() 
hubbard_7_heat

hubbard_8_heat <- hubbard_data %>%
  filter(N_years == "8") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_gradient(magma) +  
  guides(fill=guide_legend(title="Slope")) +
  labs(title = "Eight Year Interval - Hubbard Brook",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal() 
hubbard_8_heat

library(viridis)
## to only include dates where there is data for every trophic level
hubbard_8_heat <- hubbard_data %>%
  filter(N_years == "8") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") +   guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1991,2009))+
  labs(title = "Hubbard Brook ", 
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
hubbard_8_heat

library(tidyverse)
library(RColorBrewer)
#combine all locations into a data frame
all_sites <- read.csv("model_output/model_all_LTER_sites.csv")
all_heat <- all_sites %>%
  filter(N_years == "3"|N_years == "4"|N_years == "5"|N_years == "6"|
           N_years == "7"|N_years == "8") %>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu") + 
  #scale_fill_viridis(option="mako") + 
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1991,2009))+
  labs(title = "Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid(rows=vars(site), cols=vars(N_years))
all_heat

most_heat <- all_sites %>%
  filter(N_years == "5"|N_years == "6"|
           N_years == "7"|N_years == "8") %>%
  filter(site=="sbcoastal"|site=="hbrook"|site=="konza")%>%
  #group_by(start_year) %>% 
  #summarise(Mean=mean(slope)) %>%
  ## above would make fill = mean instead
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu") + 
  #scale_fill_viridis(option="mako") + 
  scale_fill_gradient2(low = "red", mid="white", high="blue")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1991,2009))+
  labs(title = "Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid(rows=vars(site), cols=vars(N_years))
most_heat





hubbard_heat <- hubbard_data %>%
  filter(N_years == "3"|N_years == "4"|
           N_years == "5"|N_years == "6"|N_years == "7"|N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1991,2009))+
  labs(title = "Hubbard Brook Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid(cols=vars(N_years))
hubbard_heat

konza_data <- read.csv("model_output/model_konza.csv")

konza_heat <- konza_data %>%
  filter(N_years == "3"|N_years == "4"|
           N_years == "5"|N_years == "6"|N_years == "7"|N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = ("red"), mid = "white", high = ("darkgreen"),
    midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1995,2010))+
  labs(title = "Konza Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid(cols=vars(N_years))
konza_heat

ntl_data <- read.csv("model_output/model_north_temperate_lakes.csv")

ntl_heat <- ntl_data %>%
  filter(N_years == "1"|N_years == "2"|N_years == "3"|N_years == "4"|
           N_years == "5"|N_years == "6"|N_years == "7"|N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1991,2009))+
  labs(title = "North Temperate Lakes Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid(cols=vars(N_years))
ntl_heat

sbc_data <- read.csv("model_output/model_santa_barbara_coastal.csv")

sbc_heat <- sbc_data %>%
  filter(N_years == "1"|N_years == "2"|N_years == "3"|N_years == "4"|
           N_years == "5"|N_years == "6"|N_years == "7"|N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1995,2015))+
  labs(title = "Santa Barbara Coastal Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid(cols=vars(N_years))
sbc_heat

sbc_heat_8 <- sbc_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1999,2014))+
  labs(title = "Santa Barbara Coastal",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid()
sbc_heat_8

kp_heat_8 <- sbc_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1999,2014))+
  labs(title = "Santa Barbara Coastal",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_grid()
sbc_heat_8
## facet by LTER site going DOWN
## facet by 3, 4, 5 ACROSS
## make sure using correct producer/consumer combos
## reach out to SITE directors to inquire
## kanza prairie could be a good start


#PLANS
## regression on regressions
## correlation between correlations
## Could use a Morans I or a smoothing stats test for getting degrees of spatial autocorrelation in data. 

## Basically testing to see if things are more similar/correlated than chance itself
## Could give us a means of calculating how correlated different groups are
library(tidyverse)
#install.packages("forecast", "gridExtra")
library(forecast)
library(gridExtra)

hubbard_data <- read.csv("model_output/model_hubbard_brook.csv")
#str(hubbard_data)
#hubbard_data <- as.numeric(hubbard_data$start_year)

############## HUB - 8 ###############  


hub_8 <- hubbard_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_8 <- hub_8[hub_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_8 <- hub_8[hub_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_8 <- hub_8[hub_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(hub_prim_8$avg_slope, hub_secon_8$avg_slope)
## between sec and tert
ccf(hub_secon_8$avg_slope, hub_tert_8$avg_slope)
## between sec and tert
ccf(hub_prim_8$avg_slope, hub_tert_8$avg_slope)

## maybe use grid extra to get them on the same image
library(gridExtra)
library(forecast)
## cross correlation plot in ggplot between prim and sec
ps_hb_8 <- ggCcf(hub_prim_8$avg_slope, hub_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between sec and tert
st_hb_8 <- ggCcf(hub_secon_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between prim and tert
pt_hb_8 <- ggCcf(hub_prim_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="Years (Lag)", y = "Correlation") +
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
grid.arrange(ps_hb_8, st_hb_8, pt_hb_8, ncol=3, 
             top=("Cross correlations of population trajectories between trophic levels at Hubbard Brook"))
hubb_8_cross <-grid.arrange(ps_hb_8, st_hb_8, pt_hb_8, ncol=3, 
                            top="Hubbard Brook")
############## HUB - 7 ###############    
hub_7 <- hubbard_data %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_7 <- hub_7[hub_7$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_7 <- hub_7[hub_7$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_7 <- hub_7[hub_7$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(hub_prim_7$avg_slope, hub_secon_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_secon_7$avg_slope, hub_tert_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_prim_7$avg_slope, hub_tert_7$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_hb_7 <- ggCcf(hub_prim_7$avg_slope, hub_secon_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Hubbard Brook 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_hb_7 <- ggCcf(hub_secon_7$avg_slope, hub_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Hubbard Brook 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_hb_7 <- ggCcf(hub_prim_7$avg_slope, hub_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Hubbard Brook 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_hb_7, st_hb_7, pt_hb_7, ncol=3)
############## HUB - 6 ###############  
hub_6 <- hubbard_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_6 <- hub_6[hub_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_6 <- hub_6[hub_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_6 <- hub_6[hub_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(hub_prim_6$avg_slope, hub_secon_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_secon_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_prim_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_hb_6 <- ggCcf(hub_prim_6$avg_slope, hub_secon_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Hubbard Brook 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_hb_6 <- ggCcf(hub_secon_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Hubbard Brook 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_hb_6 <- ggCcf(hub_prim_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Hubbard Brook 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_hb_6, st_hb_6, pt_hb_6, ncol=3)
############## HUB - 5 ###############  
hub_5 <- hubbard_data %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_5 <- hub_5[hub_5$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_5 <- hub_5[hub_5$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_5 <- hub_5[hub_5$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(hub_prim_5$avg_slope, hub_secon_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_secon_5$avg_slope, hub_tert_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_prim_5$avg_slope, hub_tert_5$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_hb_5 <- ggCcf(hub_prim_5$avg_slope, hub_secon_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Hubbard Brook 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_hb_5 <- ggCcf(hub_secon_5$avg_slope, hub_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Hubbard Brook 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_hb_5 <- ggCcf(hub_prim_5$avg_slope, hub_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Hubbard Brook 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_hb_5, st_hb_5, pt_hb_5, ncol=3)


############## HUB - 4 ###############  


hub_4 <- hubbard_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_4 <- hub_4[hub_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_4 <- hub_4[hub_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_4 <- hub_4[hub_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(hub_prim_4$avg_slope, hub_secon_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_secon_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_prim_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_hb_4 <- ggCcf(hub_prim_4$avg_slope, hub_secon_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Hubbard Brook 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_hb_4 <- ggCcf(hub_secon_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Hubbard Brook 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_hb_4 <- ggCcf(hub_prim_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Hubbard Brook 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_hb_4, st_hb_4, pt_hb_4, ncol=3)


############## HUB - 3 ############### 


hub_3 <- hubbard_data %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_3 <- hub_3[hub_3$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_3 <- hub_3[hub_3$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_3 <- hub_3[hub_3$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(hub_prim_3$avg_slope, hub_secon_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_secon_3$avg_slope, hub_tert_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(hub_prim_3$avg_slope, hub_tert_3$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_hb_3 <- ggCcf(hub_prim_3$avg_slope, hub_secon_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Hubbard Brook 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_hb_3 <- ggCcf(hub_secon_3$avg_slope, hub_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Hubbard Brook 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_hb_3 <- ggCcf(hub_prim_3$avg_slope, hub_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Hubbard Brook 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_hb_3, st_hb_3, pt_hb_3, ncol=3)

############## ALL ##############
grid.arrange(ps_hb_8, ps_hb_7, ps_hb_6, ps_hb_5, ps_hb_4, ps_hb_3,  
             st_hb_8, st_hb_7, st_hb_6, st_hb_5, st_hb_4, st_hb_3,
             pt_hb_8, pt_hb_7, pt_hb_6, pt_hb_5, pt_hb_4, pt_hb_3,
             nrow=3)

konza_data <- read.csv("model_output/model_konza.csv")


############## KONZ - 8 ###############  


konz_8 <- konza_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_8 <- konz_8[konz_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_8 <- konz_8[konz_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_8 <- konz_8[konz_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

library(forecast)
library(gridExtra)

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_8$avg_slope, konz_secon_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kz_8 <- ggCcf(konz_prim_8$avg_slope, konz_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between sec and tert
st_kz_8 <- ggCcf(konz_secon_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between prim and tert
pt_kz_8 <- ggCcf(konz_prim_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="Years (Lag)", y = "Correlation") +
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
grid.arrange(ps_kz_8, st_kz_8, pt_kz_8, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Konza Prairie")
konz_8_cross <-grid.arrange(ps_kz_8, st_kz_8, pt_kz_8, ncol=3, 
             top="Konza Prairie")

## cross correlation plot in ggplot between prim and sec
ps_kz_8_alone <- ggCcf(konz_prim_8$avg_slope, konz_secon_8$avg_slope, lag.max = 4, 
                          type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = "Primary and secondary ", x= " ", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x=element_text(size=8))
## between sec and tert
st_kz_8_alone <- ggCcf(konz_secon_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4, 
                          type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = "Secondary and tertiary", x ="Years (Lag) ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.x=element_text(size=8))
## between prim and tert
pt_kz_8_alone <- ggCcf(konz_prim_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4, 
                          type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = "Primary and tertiary ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank(), axis.title.x=element_text(size=8))
## figure code: 005
grid.arrange(ps_kz_8_alone, st_kz_8_alone, pt_kz_8_alone, ncol=3, 
             top="Konza Prairie")

############## KONZ - 7 ############### 


konz_7 <- konza_data %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_7 <- konz_7[konz_7$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_7 <- konz_7[konz_7$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_7 <- konz_7[konz_7$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_7$avg_slope, konz_secon_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kz_7 <- ggCcf(konz_prim_7$avg_slope, konz_secon_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kz_7 <- ggCcf(konz_secon_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kz_7 <- ggCcf(konz_prim_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kz_7, st_kz_7, pt_kz_7, ncol=3)


############## KONZ - 6 ###############  


konz_6 <- konza_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_6 <- konz_6[konz_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_6 <- konz_6[konz_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_6 <- konz_6[konz_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_6$avg_slope, konz_secon_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kz_6 <- ggCcf(konz_prim_6$avg_slope, konz_secon_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kz_6 <- ggCcf(konz_secon_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kz_6 <- ggCcf(konz_prim_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kz_6, st_kz_6, pt_kz_6, ncol=3)


############## KONZ - 5 ############### 


konz_5 <- konza_data %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_5 <- konz_5[konz_5$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_5 <- konz_5[konz_5$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_5 <- konz_5[konz_5$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_5$avg_slope, konz_secon_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kz_5 <- ggCcf(konz_prim_5$avg_slope, konz_secon_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kz_5 <- ggCcf(konz_secon_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kz_5 <- ggCcf(konz_prim_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kz_5, st_kz_5, pt_kz_5, ncol=3)


############## KONZ - 4 ###############  


konz_4 <- konza_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_4 <- konz_4[konz_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_4 <- konz_4[konz_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_4 <- konz_4[konz_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_4$avg_slope, konz_secon_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kz_4 <- ggCcf(konz_prim_4$avg_slope, konz_secon_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kz_4 <- ggCcf(konz_secon_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kz_4 <- ggCcf(konz_prim_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kz_4, st_kz_4, pt_kz_4, ncol=3)


############## KONZ - 3 ############### 


konz_3 <- konza_data %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_3 <- konz_3[konz_3$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_3 <- konz_3[konz_3$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_3 <- konz_3[konz_3$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_3$avg_slope, konz_secon_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kz_3 <- ggCcf(konz_prim_3$avg_slope, konz_secon_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kz_3 <- ggCcf(konz_secon_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kz_3 <- ggCcf(konz_prim_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kz_3, st_kz_3, pt_kz_3, ncol=3)

############## ALL ##############
grid.arrange(ps_kz_8, ps_kz_7, ps_kz_6, ps_kz_5, ps_kz_4, ps_kz_3,  
             st_kz_8, st_kz_7, st_kz_6, st_kz_5, st_kz_4, st_kz_3,
             pt_kz_8, pt_kz_7, pt_kz_6, pt_kz_5, pt_kz_4, pt_kz_3,
             nrow=3)


santa_data <- read.csv("model_output/model_santa_barbara_coastal.csv")


############## SANT - 8 ###############  


sant_8 <- santa_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_8 <- sant_8[sant_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_8 <- sant_8[sant_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_8 <- sant_8[sant_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(sant_prim_8$avg_slope, sant_secon_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_secon_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_prim_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_sb_8 <- ggCcf(sant_prim_8$avg_slope, sant_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between sec and tert
st_sb_8 <- ggCcf(sant_secon_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between prim and tert
pt_sb_8 <- ggCcf(sant_prim_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="Years (Lag)", y = "Correlation") +
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
grid.arrange(ps_sb_8, st_sb_8, pt_sb_8, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Santa Barbara Coastal" 
                          )
sant_8_cross <-grid.arrange(ps_sb_8, st_sb_8, pt_sb_8, ncol=3, 
                            top="Santa Barbara Coastal")

############## SANT - 7 ############### 


sant_7 <- santa_data %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_7 <- sant_7[sant_7$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_7 <- sant_7[sant_7$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_7 <- sant_7[sant_7$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(sant_prim_7$avg_slope, sant_secon_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_secon_7$avg_slope, sant_tert_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_prim_7$avg_slope, sant_tert_7$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_sb_7 <- ggCcf(sant_prim_7$avg_slope, sant_secon_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Santa Barbara 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_sb_7 <- ggCcf(sant_secon_7$avg_slope, sant_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Santa Barbara 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_sb_7 <- ggCcf(sant_prim_7$avg_slope, sant_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Santa Barbara 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_sb_7, st_sb_7, pt_sb_7, ncol=3)


############## SANT - 6 ###############  


sant_6 <- santa_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_6 <- sant_6[sant_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_6 <- sant_6[sant_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_6 <- sant_6[sant_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(sant_prim_6$avg_slope, sant_secon_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_secon_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_prim_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_sb_6 <- ggCcf(sant_prim_6$avg_slope, sant_secon_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Santa Barbara 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_sb_6 <- ggCcf(sant_secon_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Santa Barbara 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_sb_6 <- ggCcf(sant_prim_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Santa Barbara 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_sb_6, st_sb_6, pt_sb_6, ncol=3)


############## SANT - 5 ############### 


sant_5 <- santa_data %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_5 <- sant_5[sant_5$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_5 <- sant_5[sant_5$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_5 <- sant_5[sant_5$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(sant_prim_5$avg_slope, sant_secon_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_secon_5$avg_slope, sant_tert_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_prim_5$avg_slope, sant_tert_5$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_sb_5 <- ggCcf(sant_prim_5$avg_slope, sant_secon_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Santa Barbara 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_sb_5 <- ggCcf(sant_secon_5$avg_slope, sant_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Santa Barbara 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_sb_5 <- ggCcf(sant_prim_5$avg_slope, sant_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Santa Barbara 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_sb_5, st_sb_5, pt_sb_5, ncol=3)


############## SANT - 4 ###############  


sant_4 <- santa_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_4 <- sant_4[sant_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_4 <- sant_4[sant_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_4 <- sant_4[sant_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(sant_prim_4$avg_slope, sant_secon_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_secon_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_prim_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_sb_4 <- ggCcf(sant_prim_4$avg_slope, sant_secon_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Santa Barbara 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_sb_4 <- ggCcf(sant_secon_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Santa Barbara 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_sb_4 <- ggCcf(sant_prim_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Santa Barbara 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_sb_4, st_sb_4, pt_sb_4, ncol=3)


############## SANT - 3 ############### 


sant_3 <- santa_data %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_3 <- sant_3[sant_3$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_3 <- sant_3[sant_3$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_3 <- sant_3[sant_3$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(sant_prim_3$avg_slope, sant_secon_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_secon_3$avg_slope, sant_tert_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(sant_prim_3$avg_slope, sant_tert_3$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_sb_3 <- ggCcf(sant_prim_3$avg_slope, sant_secon_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Santa Barbara 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_sb_3 <- ggCcf(sant_secon_3$avg_slope, sant_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Santa Barbara 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_sb_3 <- ggCcf(sant_prim_3$avg_slope, sant_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Santa Barbara 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_sb_3, st_sb_3, pt_sb_3, ncol=3)

############## ALL ##############
grid.arrange(ps_sb_8, ps_sb_7, ps_sb_6, ps_sb_5, ps_sb_4, ps_sb_3,  
             st_sb_8, st_sb_7, st_sb_6, st_sb_5, st_sb_4, st_sb_3,
             pt_sb_8, pt_sb_7, pt_sb_6, pt_sb_5, pt_sb_4, pt_sb_3,
             nrow=3)

temp_data <- read.csv("model_output/model_north_temperate_lakes.csv")


############## TEMP - 8 ###############  


temp_8 <- temp_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
temp_prim_8 <- temp_8[temp_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_secon_8 <- temp_8[temp_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_tert_8 <- temp_8[temp_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(temp_prim_8$avg_slope, temp_secon_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_secon_8$avg_slope, temp_tert_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_prim_8$avg_slope, temp_tert_8$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_nt_8 <- ggCcf(temp_prim_8$avg_slope, temp_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary",
       subtitle ="North Temp Lakes 8 year", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_nt_8 <- ggCcf(temp_secon_8$avg_slope, temp_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="North Temp Lakes 8 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_nt_8 <- ggCcf(temp_prim_8$avg_slope, temp_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="North Temp Lakes 8 year  ", x ="Years (Lag)", y = "Correlation") +
  theme_minimal()
grid.arrange(ps_nt_8, st_nt_8, pt_nt_8, ncol=3)


############## TEMP - 7 ############### 


temp_7 <- temp_data %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
temp_prim_7 <- temp_7[temp_7$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_secon_7 <- temp_7[temp_7$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_tert_7 <- temp_7[temp_7$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(temp_prim_7$avg_slope, temp_secon_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_secon_7$avg_slope, temp_tert_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_prim_7$avg_slope, temp_tert_7$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_nt_7 <- ggCcf(temp_prim_7$avg_slope, temp_secon_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="North Temp Lakes 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_nt_7 <- ggCcf(temp_secon_7$avg_slope, temp_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="North Temp Lakes 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_nt_7 <- ggCcf(temp_prim_7$avg_slope, temp_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="North Temp Lakes 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_nt_7, st_nt_7, pt_nt_7, ncol=3)


############## TEMP - 6 ###############  


temp_6 <- temp_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
temp_prim_6 <- temp_6[temp_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_secon_6 <- temp_6[temp_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_tert_6 <- temp_6[temp_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(temp_prim_6$avg_slope, temp_secon_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_secon_6$avg_slope, temp_tert_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_prim_6$avg_slope, temp_tert_6$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_nt_6 <- ggCcf(temp_prim_6$avg_slope, temp_secon_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="North Temp Lakes 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_nt_6 <- ggCcf(temp_secon_6$avg_slope, temp_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="North Temp Lakes 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_nt_6 <- ggCcf(temp_prim_6$avg_slope, temp_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="North Temp Lakes 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_nt_6, st_nt_6, pt_nt_6, ncol=3)


############## TEMP - 5 ############### 


temp_5 <- temp_data %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
temp_prim_5 <- temp_5[temp_5$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_secon_5 <- temp_5[temp_5$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_tert_5 <- temp_5[temp_5$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(temp_prim_5$avg_slope, temp_secon_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_secon_5$avg_slope, temp_tert_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_prim_5$avg_slope, temp_tert_5$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_nt_5 <- ggCcf(temp_prim_5$avg_slope, temp_secon_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="North Temp Lakes 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_nt_5 <- ggCcf(temp_secon_5$avg_slope, temp_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="North Temp Lakes 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_nt_5 <- ggCcf(temp_prim_5$avg_slope, temp_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="North Temp Lakes 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_nt_5, st_nt_5, pt_nt_5, ncol=3)


############## TEMP - 4 ###############  


temp_4 <- temp_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
temp_prim_4 <- temp_4[temp_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_secon_4 <- temp_4[temp_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_tert_4 <- temp_4[temp_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(temp_prim_4$avg_slope, temp_secon_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_secon_4$avg_slope, temp_tert_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_prim_4$avg_slope, temp_tert_4$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_nt_4 <- ggCcf(temp_prim_4$avg_slope, temp_secon_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="North Temp Lakes 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_nt_4 <- ggCcf(temp_secon_4$avg_slope, temp_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="North Temp Lakes 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_nt_4 <- ggCcf(temp_prim_4$avg_slope, temp_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="North Temp Lakes 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_nt_4, st_nt_4, pt_nt_4, ncol=3)


############## TEMP - 3 ############### 


temp_3 <- temp_data %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
temp_prim_3 <- temp_3[temp_3$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_secon_3 <- temp_3[temp_3$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
temp_tert_3 <- temp_3[temp_3$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(temp_prim_3$avg_slope, temp_secon_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_secon_3$avg_slope, temp_tert_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(temp_prim_3$avg_slope, temp_tert_3$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_nt_3 <- ggCcf(temp_prim_3$avg_slope, temp_secon_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="North Temp Lakes 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_nt_3 <- ggCcf(temp_secon_3$avg_slope, temp_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="North Temp Lakes 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_nt_3 <- ggCcf(temp_prim_3$avg_slope, temp_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="North Temp Lakes 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_nt_3, st_nt_3, pt_nt_3, ncol=3)

############## ALL ##############
grid.arrange(ps_nt_8, ps_nt_7, ps_nt_6, ps_nt_5, ps_nt_4, ps_nt_3,  
             st_nt_8, st_nt_7, st_nt_6, st_nt_5, st_nt_4, st_nt_3,
             pt_nt_8, pt_nt_7, pt_nt_6, pt_nt_5, pt_nt_4, pt_nt_3,
             nrow=3)


## poster graphics
library(gridExtra)
library(tidyverse)
library(forecast)
library(grid)
hubbard_data <- read.csv("model_output/model_hubbard_brook.csv")
hub_8 <- hubbard_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_8 <- hub_8[hub_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_8 <- hub_8[hub_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_8 <- hub_8[hub_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_hb_8 <- ggCcf(hub_prim_8$avg_slope, hub_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x=element_blank())
## between sec and tert
st_hb_8 <- ggCcf(hub_secon_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank())

## between prim and tert
pt_hb_8 <- ggCcf(hub_prim_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="", y = " ") +
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank())
grid.arrange(ps_hb_8, st_hb_8, pt_hb_8, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Hubbard Brook" 
                          )
hubb_8_cross <-grid.arrange(ps_hb_8, st_hb_8, pt_hb_8, ncol=3, 
                            top="Hubbard Brook")


konzakat_data <- read.csv("model_output/model_konzakat.csv")

konzkat_8 <- konzakat_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konzkat_prim_8 <- konzkat_8[konzkat_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_8 <- konzkat_8[konzkat_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_8 <- konzkat_8[konzkat_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_kzkat_8 <- ggCcf(konzkat_prim_8$avg_slope, konzkat_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x= " ", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title.y=element_text(size=8), axis.title.x=element_blank())
## between sec and tert
st_kzkat_8 <- ggCcf(konzkat_secon_8$avg_slope, konzkat_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title=element_blank())
## between prim and tert
pt_kzkat_8 <- ggCcf(konzkat_prim_8$avg_slope, konzkat_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title=element_blank())
grid.arrange(ps_kzkat_8, st_kzkat_8, pt_kzkat_8, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Konza Prairie"
                          )
konzkat_8_cross <-grid.arrange(ps_kzkat_8, st_kzkat_8, pt_kzkat_8, ncol=3, 
                            top="Konza Prairie")

santa_data <- read.csv("model_output/model_santa_barbara_coastal.csv")

sant_8 <- santa_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_8 <- sant_8[sant_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_8 <- sant_8[sant_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_8 <- sant_8[sant_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_sb_8 <- ggCcf(sant_prim_8$avg_slope, sant_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", x ="", y = "Correlation")+
  ggtitle("Primary and secondary")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_blank(), axis.title.y=element_text(size=8), axis.title.x = element_text(size=8))
## between sec and tert
st_sb_8 <- ggCcf(sant_secon_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  ylim (-.75,.75)+
  xlim(-8,8)+
  labs(title="Secondary and tertiary", x ="Years (Lag)", y = " ")+
  ggtitle("Secondary and tertiary")+
  theme(plot.title=element_blank(), axis.title.x=element_text(size=8), axis.title.y=element_blank())
## between prim and tert
pt_sb_8 <- ggCcf(sant_prim_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary", x ="", y = " ")+
  ggtitle("Primary and tertiary")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_blank(), axis.title.x=element_text(size=8), axis.title.y=element_blank())
grid.arrange(ps_sb_8, st_sb_8, pt_sb_8, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Santa Barbara Coastal")
sant_8_cross <-grid.arrange(ps_sb_8, st_sb_8, pt_sb_8, ncol=3, 
                            top="Santa Barbara Coastal")

sant_8_cross_sing <-grid.arrange(ps_sb_8, st_sb_8, pt_sb_8, ncol=3)

## figure code: 000 ###
grid.arrange(hubb_8_cross, konzkat_8_cross, sant_8_cross,  nrow=3, top="Eight Year Window")
## figure code: 001 ##
grid.arrange(hubb_8_cross, konzkat_8_cross, sant_8_cross,  nrow=3, 
             top="Cross correlations of population trajectories between trophic levels at different LTER sites")
grid.arrange(sant_8_cross_sing,  
             top="Cross correlations of population trajectories between trophic levels ")
sant_8_cross

############################
ps_sb_8 <- ggCcf(sant_prim_8$avg_slope, sant_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="Lag (years)", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x=element_text(size=8))
## between sec and tert
st_sb_8 <- ggCcf(sant_secon_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="", y = "Lag (years) ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10))

## between prim and tert
pt_sb_8 <- ggCcf(sant_prim_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="", y = "Lag (years) ") +
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10))
grid.arrange(ps_hb_8, st_hb_8, pt_hb_8, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels")
hubb_8_cross <-grid.arrange(ps_hb_8, st_hb_8, pt_hb_8, ncol=3, 
                            top="Hubbard Brook")

### lets make some complete plots so they can exist separatley

## cross correlation plot in ggplot between prim and sec
ps_hb_8_alone <- ggCcf(hub_prim_8$avg_slope, hub_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x=element_text(size=8))
## between sec and tert
st_hb_8_alone <- ggCcf(hub_secon_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="Years (Lag)", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.x=element_text(size=8))

## between prim and tert
pt_hb_8_alone <- ggCcf(hub_prim_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="", y = " ") +
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.x=element_text(size=8))

## figure code: 002 ##
grid.arrange(ps_hb_8_alone, st_hb_8_alone, pt_hb_8_alone, ncol=3, 
             top="Hubbard Brook")


## cross correlation plot in ggplot between prim and sec
ps_kzkat_8_alone <- ggCcf(konzkat_prim_8$avg_slope, konzkat_secon_8$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = "Primary and secondary ", x= " ", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x=element_text(size=8))
## between sec and tert
st_kzkat_8_alone <- ggCcf(konzkat_secon_8$avg_slope, konzkat_tert_8$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = "Secondary and tertiary", x ="Years (Lag) ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.x=element_text(size=8))
## between prim and tert
pt_kzkat_8_alone <- ggCcf(konzkat_prim_8$avg_slope, konzkat_tert_8$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = "Primary and tertiary ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank(), axis.title.x=element_text(size=8))
## figure code: 003
grid.arrange(ps_kzkat_8_alone, st_kzkat_8_alone, pt_kzkat_8_alone, ncol=3, 
             top="Konza Prairie")


## cross correlation plot in ggplot between prim and sec
ps_sb_8_alone <- ggCcf(sant_prim_8$avg_slope, sant_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", x ="", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x = element_text(size=8))
## between sec and tert
st_sb_8_alone <- ggCcf(sant_secon_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  ylim (-.75,.75)+
  xlim(-8,8)+
  labs(title="Secondary and tertiary", x ="Years (Lag)", y = " ")+
  theme(plot.title=element_text(size=10), axis.title.x=element_text(size=8), axis.title.y=element_blank())
## between prim and tert
pt_sb_8_alone <- ggCcf(sant_prim_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary", x ="", y = " ")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_text(size=10), axis.title.x=element_text(size=8), axis.title.y=element_blank())

## figure code: 004
grid.arrange(ps_sb_8_alone, st_sb_8_alone, pt_sb_8_alone, ncol=3, 
             top="Santa Barbara Coastal")


######### trying out new adjustments

konzakat_data <- read.csv("model_output/model_konzakat.csv")


############## KONZ - 8 ###############  


konz_8 <- konzakat_data %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_8 <- konz_8[konz_8$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_8 <- konz_8[konz_8$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_8 <- konz_8[konz_8$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

library(forecast)
library(gridExtra)
library(ggplot2)
library(grid)

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_8$avg_slope, konz_secon_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kzkat_8 <- ggCcf(konz_prim_8$avg_slope, konz_secon_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between sec and tert
st_kzkat_8 <- ggCcf(konz_secon_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="Years (Lag)", y = "Correlation")+
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
## between prim and tert
pt_kzkat_8 <- ggCcf(konz_prim_8$avg_slope, konz_tert_8$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="Years (Lag)", y = "Correlation") +
  ylim (-.75,.75)+
  theme(plot.title=element_text(size=10))
grid.arrange(ps_kzkat_8, st_kzkat_8, pt_kzkat_8, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Konza Prairie")
konz_8_cross <-grid.arrange(ps_kzkat_8, st_kzkat_8, pt_kzkat_8, ncol=3, 
                            top="Konza Prairie")

############## KONZ - 7 ############### 


konz_7 <- konzakat_data %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_7 <- konz_7[konz_7$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_7 <- konz_7[konz_7$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_7 <- konz_7[konz_7$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_7$avg_slope, konz_secon_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kzkat_7 <- ggCcf(konz_prim_7$avg_slope, konz_secon_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kzkat_7 <- ggCcf(konz_secon_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kzkat_7 <- ggCcf(konz_prim_7$avg_slope, konz_tert_7$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 7 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kzkat_7, st_kzkat_7, pt_kzkat_7, ncol=3)


############## KONZ - 6 ###############  


konz_6 <- konzakat_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_6 <- konz_6[konz_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_6 <- konz_6[konz_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_6 <- konz_6[konz_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_6$avg_slope, konz_secon_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kzkat_6 <- ggCcf(konz_prim_6$avg_slope, konz_secon_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kzkat_6 <- ggCcf(konz_secon_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kzkat_6 <- ggCcf(konz_prim_6$avg_slope, konz_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 6 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kzkat_6, st_kzkat_6, pt_kzkat_6, ncol=3)


############## KONZ - 5 ############### 


konz_5 <- konzakat_data %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_5 <- konz_5[konz_5$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_5 <- konz_5[konz_5$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_5 <- konz_5[konz_5$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_5$avg_slope, konz_secon_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kzkat_5 <- ggCcf(konz_prim_5$avg_slope, konz_secon_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kzkat_5 <- ggCcf(konz_secon_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kzkat_5 <- ggCcf(konz_prim_5$avg_slope, konz_tert_5$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 5 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kzkat_5, st_kzkat_5, pt_kzkat_5, ncol=3)


############## KONZ - 4 ###############  


konz_4 <- konzakat_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_4 <- konz_4[konz_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_4 <- konz_4[konz_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_4 <- konz_4[konz_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_4$avg_slope, konz_secon_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kzkat_4 <- ggCcf(konz_prim_4$avg_slope, konz_secon_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kzkat_4 <- ggCcf(konz_secon_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kzkat_4 <- ggCcf(konz_prim_4$avg_slope, konz_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 4 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kzkat_4, st_kzkat_4, pt_kzkat_4, ncol=3)


############## KONZ - 3 ############### 


konz_3 <- konzakat_data %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konz_prim_3 <- konz_3[konz_3$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_secon_3 <- konz_3[konz_3$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konz_tert_3 <- konz_3[konz_3$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

par(mfrow=c(3,1))
## cross correlation plot in base R between prim and sec
ccf(konz_prim_3$avg_slope, konz_secon_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_secon_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4)
## between sec and tert
ccf(konz_prim_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4)

## maybe use grid extra to get them on the same image
## cross correlation plot in ggplot between prim and sec
ps_kzkat_3 <- ggCcf(konz_prim_3$avg_slope, konz_secon_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary  ",
       subtitle ="Konza 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between sec and tert
st_kzkat_3 <- ggCcf(konz_secon_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ",
       subtitle ="Konza 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
## between prim and tert
pt_kzkat_3 <- ggCcf(konz_prim_3$avg_slope, konz_tert_3$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ",
       subtitle ="Konza 3 year  ", x ="Years (Lag)", y = "Correlation")+
  theme_minimal()
grid.arrange(ps_kzkat_3, st_kzkat_3, pt_kzkat_3, ncol=3)

############## ALL ##############
grid.arrange(ps_kzkat_8, ps_kzkat_7, ps_kzkat_6, ps_kzkat_5, ps_kzkat_4, ps_kzkat_3,  
             st_kzkat_8, st_kzkat_7, st_kzkat_6, st_kzkat_5, st_kzkat_4, st_kzkat_3,
             pt_kzkat_8, pt_kzkat_7, pt_kzkat_6, pt_kzkat_5, pt_kzkat_4, pt_kzkat_3,
             nrow=3)



konzakat_data <- read.csv("model_output/model_konzakat.csv")

## figure code: 006
#konzakat_heat

konza_data <- read.csv("model_output/model_konza.csv")

konza_heat <- konza_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  #scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("blue"),
  #                     midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1995,2015))+
  labs(title = "Konza Prairie Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()



konza_heat <- konza_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("darkgreen"),
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1995,2015))+
  labs(title = "Konza Prairie Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
## figure code: 006
konza_heat

hubbard_data <- read.csv("model_output/model_hubbard_brook.csv")

hubbard_heat <- hubbard_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  #scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("darkgreen"),
  #                     midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1990,2013))+
  labs(title = "Hubbard Brook Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
## figure code: 007
hubbard_heat

sbc_heat <- sbc_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  scale_fill_distiller(palette="RdBu", trans="reverse") + 
  #scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("darkgreen"),
  #                     midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1998,2015))+
  labs(title = "Santa Barbara Coastal Tritrophic Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
## figure code: 008
sbc_heat

konza_heat

## figure code 010
sbc_heat_5 <- sbc_data %>%
  filter(N_years == "5") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = "red", mid="white", high="blue")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1998,2015))+
  labs(title = "Santa Barbara Coastal 5 Year Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
sbc_heat_5
## figure code 011
sbc_heat_3 <- sbc_data %>%
  filter(N_years == "3") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("darkgreen"),
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1998,2015))+
  labs(title = "Santa Barbara Coastal 3 Year Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
sbc_heat_3

## figure code 013
hubbard_heat_5 <- hubbard_data %>%
  filter(N_years == "5") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = "red", mid="white", high="blue")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1990,2013))+
  labs(title = "Hubbard Brook 5 Year Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
hubbard_heat_5

## figure code 014
hubbard_heat_3 <- hubbard_data %>%
  filter(N_years == "3") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("darkgreen"),
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1990,2013))+
  labs(title = "Hubbard Brook 3 Year Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
hubbard_heat_3

## figure code 015
konz_heat_8 <- konza_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = "red", mid="white", high="blue")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1995,2015))+
  labs(title = "Konza Prairie",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
konz_heat_8
## figure code 012
hubbard_heat_8 <- hubbard_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = "red", mid="white", high="blue")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1990,2013))+
  labs(title = "Hubbard Brook",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
hubbard_heat_8
## figure code 009 
sbc_heat_8 <- sbc_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = "red", mid="white", high="blue")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1998,2015))+
  labs(title = "Santa Barbara Coastal",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
sbc_heat_8
## December 9 new plot
konza_data$site <- "Konza Prairie (B)"
hubbard_data$site <- "Hubbard Brook (A)"
sbc_data$site <- "Santa Barbara Coastal (C)"
total_data<- rbind(konza_data, hubbard_data, sbc_data)
gg_all <- total_data %>%
  filter(N_years == "8") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = "#377EB8", mid="white", high="darkred")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1990,2015))+
  labs(x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()+
  facet_wrap("site", nrow=3)
gg_all
ggarrange(konz_heat_8, hubbard_heat_8, sbc_heat_8,
          ncol = 1, nrow = 3)
## figure code 016
konz_heat_5 <- konza_data %>%
  filter(N_years == "5") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("darkgreen"),
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1995,2015))+
  labs(title = "Konza Prairie 5 Year Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
konz_heat_5


## figure code 017
konz_heat_3 <- konza_data %>%
  filter(N_years == "3") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = slope),colour = "white") +
  #scale_fill_distiller(palette="RdBu", trans="reverse") + 
  scale_fill_gradient2(low = ("darkred"), mid = "white", high = ("darkgreen"),
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  guides(fill=guide_legend(title="Slope")) +
  scale_x_continuous(limits=c(1995,2015))+
  labs(title = "Konza Prairie 3 Year Trends",
       x = "Start Year", y = "Trophic Level") +
  theme_bw() + theme_minimal()
konz_heat_3

## figure code #018
hubbard_data <- read.csv("model_output/model_hubbard_brook.csv")
hub_6 <- hubbard_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_6 <- hub_6[hub_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_6 <- hub_6[hub_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_6 <- hub_6[hub_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_hb_6 <- ggCcf(hub_prim_6$avg_slope, hub_secon_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x=element_blank())
## between sec and tert
st_hb_6 <- ggCcf(hub_secon_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank())

## between prim and tert
pt_hb_6 <- ggCcf(hub_prim_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="", y = " ") +
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank())
grid.arrange(ps_hb_6, st_hb_6, pt_hb_6, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Hubbard Brook (6)")
hubb_6_cross <-grid.arrange(ps_hb_6, st_hb_6, pt_hb_6, ncol=3, 
                            top="Hubbard Brook")


konzakat_data <- read.csv("model_output/model_konzakat.csv")

konzkat_6 <- konzakat_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konzkat_prim_6 <- konzkat_6[konzkat_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_6 <- konzkat_6[konzkat_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_6 <- konzkat_6[konzkat_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_kzkat_6 <- ggCcf(konzkat_prim_6$avg_slope, konzkat_secon_6$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x= " ", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title.y=element_text(size=8), axis.title.x=element_blank())
## between sec and tert
st_kzkat_6 <- ggCcf(konzkat_secon_6$avg_slope, konzkat_tert_6$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title=element_blank())
## between prim and tert
pt_kzkat_6 <- ggCcf(konzkat_prim_6$avg_slope, konzkat_tert_6$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title=element_blank())
grid.arrange(ps_kzkat_6, st_kzkat_6, pt_kzkat_6, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Konza Prairie (6)")
konzkat_6_cross <-grid.arrange(ps_kzkat_6, st_kzkat_6, pt_kzkat_6, ncol=3, 
                               top="Konza Prairie")

santa_data <- read.csv("model_output/model_santa_barbara_coastal.csv")

sant_6 <- santa_data %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_6 <- sant_6[sant_6$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_6 <- sant_6[sant_6$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_6 <- sant_6[sant_6$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_sb_6 <- ggCcf(sant_prim_6$avg_slope, sant_secon_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", x ="", y = "Correlation")+
  ggtitle("Primary and secondary")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_blank(), axis.title.y=element_text(size=8), axis.title.x = element_text(size=8))
## between sec and tert
st_sb_6 <- ggCcf(sant_secon_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  ylim (-.75,.75)+
  xlim(-8,8)+
  labs(title="Secondary and tertiary", x ="Years (Lag)", y = " ")+
  ggtitle("Secondary and tertiary")+
  theme(plot.title=element_blank(), axis.title.x=element_text(size=8), axis.title.y=element_blank())
## between prim and tert
pt_sb_6 <- ggCcf(sant_prim_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary", x ="", y = " ")+
  ggtitle("Primary and tertiary")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_blank(), axis.title.x=element_text(size=8), axis.title.y=element_blank())
grid.arrange(ps_sb_6, st_sb_6, pt_sb_6, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Santa Barbara Coastal (5)")
sant_6_cross <-grid.arrange(ps_sb_6, st_sb_6, pt_sb_6, ncol=3, 
                            top="Santa Barbara Coastal")

sant_6_cross_sing <-grid.arrange(ps_sb_6, st_sb_6, pt_sb_6, ncol=3)

## figure code: 018 ###
grid.arrange(hubb_6_cross, konzkat_6_cross, sant_6_cross,  nrow=3, top="Six Year Window")



## figure code #019
hubbard_data <- read.csv("model_output/model_hubbard_brook.csv")
hub_4 <- hubbard_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
hub_prim_4 <- hub_4[hub_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_4 <- hub_4[hub_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_4 <- hub_4[hub_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_hb_4 <- ggCcf(hub_prim_4$avg_slope, hub_secon_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", fontsize = 10, x ="", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title.y=element_text(size=8), axis.title.x=element_blank())
## between sec and tert
st_hb_4 <- ggCcf(hub_secon_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Secondary and tertiary  ", x ="", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank())

## between prim and tert
pt_hb_4 <- ggCcf(hub_prim_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary  ", x ="", y = " ") +
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_text(size=10), axis.title=element_blank())
grid.arrange(ps_hb_4, st_hb_4, pt_hb_4, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Hubbard Brook (3)")
hubb_4_cross <-grid.arrange(ps_hb_4, st_hb_4, pt_hb_4, ncol=3, 
                            top="Hubbard Brook")


konzakat_data <- read.csv("model_output/model_konzakat.csv")

konzkat_4 <- konzakat_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
konzkat_prim_4 <- konzkat_4[konzkat_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_4 <- konzkat_4[konzkat_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_4 <- konzkat_4[konzkat_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_kzkat_4 <- ggCcf(konzkat_prim_4$avg_slope, konzkat_secon_4$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x= " ", y = "Correlation")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title.y=element_text(size=8), axis.title.x=element_blank())
## between sec and tert
st_kzkat_4 <- ggCcf(konzkat_secon_4$avg_slope, konzkat_tert_4$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title=element_blank())
## between prim and tert
pt_kzkat_4 <- ggCcf(konzkat_prim_4$avg_slope, konzkat_tert_4$avg_slope, lag.max = 4, 
                    type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title = " ", x =" ", y = " ")+
  ylim (-.75,.75)+
  xlim(-10,10)+
  theme(plot.title=element_blank(), axis.title=element_blank())
grid.arrange(ps_kzkat_4, st_kzkat_4, pt_kzkat_4, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Konza Prairie (3)")
konzkat_4_cross <-grid.arrange(ps_kzkat_4, st_kzkat_4, pt_kzkat_4, ncol=3, 
                               top="Konza Prairie")

santa_data <- read.csv("model_output/model_santa_barbara_coastal.csv")

sant_4 <- santa_data %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) 
# need to manually check the start and end dates in common so arrays match up
sant_prim_4 <- sant_4[sant_4$trophic_level == 'Primary',] %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_4 <- sant_4[sant_4$trophic_level == 'Secondary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_4 <- sant_4[sant_4$trophic_level == 'Tertiary',]%>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

## cross correlation plot in ggplot between prim and sec
ps_sb_4 <- ggCcf(sant_prim_4$avg_slope, sant_secon_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and secondary", x ="", y = "Correlation")+
  ggtitle("Primary and secondary")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_blank(), axis.title.y=element_text(size=8), axis.title.x = element_text(size=8))
## between sec and tert
st_sb_4 <- ggCcf(sant_secon_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  ylim (-.75,.75)+
  xlim(-8,8)+
  labs(title="Secondary and tertiary", x ="Years (Lag)", y = " ")+
  ggtitle("Secondary and tertiary")+
  theme(plot.title=element_blank(), axis.title.x=element_text(size=8), axis.title.y=element_blank())
## between prim and tert
pt_sb_4 <- ggCcf(sant_prim_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4, 
                 type = c("correlation", "covariance"), plot = TRUE, na.action = na.contiguous) +
  labs(title="Primary and tertiary", x ="", y = " ")+
  ggtitle("Primary and tertiary")+
  ylim (-.75,.75)+
  xlim(-8,8)+
  theme(plot.title=element_blank(), axis.title.x=element_text(size=8), axis.title.y=element_blank())
grid.arrange(ps_sb_4, st_sb_4, pt_sb_4, ncol=3, 
             top="Cross correlations of population trajectories between trophic levels at Santa Barbara Coastal (4)")
sant_4_cross <-grid.arrange(ps_sb_4, st_sb_4, pt_sb_4, ncol=3, 
                            top="Santa Barbara Coastal")

sant_4_cross_sing <-grid.arrange(ps_sb_4, st_sb_4, pt_sb_4, ncol=3)

## figure code: 019 ###
grid.arrange(hubb_4_cross, konzkat_4_cross, sant_4_cross,  nrow=3, top="Four Year Window")


################################################
## fitting curves
################################################

## okay pseudocode
## create a dataframe to store each set of outputs
## 6 years x 3 sites x 3 troph combos = 54 data frames
## something like, 

# make all these damn subgroups

sant_prim_8 <- santa_data[santa_data$trophic_level == 'Primary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_prim_7 <- santa_data[santa_data$trophic_level == 'Primary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_prim_6 <- santa_data[santa_data$trophic_level == 'Primary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_prim_5 <- santa_data[santa_data$trophic_level == 'Primary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_prim_4 <- santa_data[santa_data$trophic_level == 'Primary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_prim_3 <- santa_data[santa_data$trophic_level == 'Primary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_prim_8 <- hubbard_data[hubbard_data$trophic_level == 'Primary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_prim_7 <- hubbard_data[hubbard_data$trophic_level == 'Primary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_prim_6 <- hubbard_data[hubbard_data$trophic_level == 'Primary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_prim_5 <- hubbard_data[hubbard_data$trophic_level == 'Primary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_prim_4 <- hubbard_data[hubbard_data$trophic_level == 'Primary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_prim_3 <- hubbard_data[hubbard_data$trophic_level == 'Primary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_prim_8 <- konzakat_data[konzakat_data$trophic_level == 'Primary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_prim_7 <- konzakat_data[konzakat_data$trophic_level == 'Primary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_prim_6 <- konzakat_data[konzakat_data$trophic_level == 'Primary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_prim_5 <- konzakat_data[konzakat_data$trophic_level == 'Primary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_prim_4 <- konzakat_data[konzakat_data$trophic_level == 'Primary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_prim_3 <- konzakat_data[konzakat_data$trophic_level == 'Primary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
## secon
sant_secon_8 <- santa_data[santa_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_7 <- santa_data[santa_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_6 <- santa_data[santa_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_5 <- santa_data[santa_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_4 <- santa_data[santa_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_secon_3 <- santa_data[santa_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_8 <- hubbard_data[hubbard_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_7 <- hubbard_data[hubbard_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_6 <- hubbard_data[hubbard_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_5 <- hubbard_data[hubbard_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_4 <- hubbard_data[hubbard_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_secon_3 <- hubbard_data[hubbard_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_8 <- konzakat_data[konzakat_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_7 <- konzakat_data[konzakat_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_6 <- konzakat_data[konzakat_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_5 <- konzakat_data[konzakat_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_4 <- konzakat_data[konzakat_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_secon_3 <- konzakat_data[konzakat_data$trophic_level == 'Secondary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
##tert
sant_tert_8 <- santa_data[santa_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_7 <- santa_data[santa_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_6 <- santa_data[santa_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_5 <- santa_data[santa_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_4 <- santa_data[santa_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sant_tert_3 <- santa_data[santa_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_8 <- hubbard_data[hubbard_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_7 <- hubbard_data[hubbard_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_6 <- hubbard_data[hubbard_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_5 <- hubbard_data[hubbard_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_4 <- hubbard_data[hubbard_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hub_tert_3 <- hubbard_data[hubbard_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_8 <- konzakat_data[konzakat_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "8" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_7 <- konzakat_data[konzakat_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "7" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_6 <- konzakat_data[konzakat_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "6" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_5 <- konzakat_data[konzakat_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "5" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_4 <- konzakat_data[konzakat_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "4" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konzkat_tert_3 <- konzakat_data[konzakat_data$trophic_level == 'Tertiary',] %>%
  filter(N_years == "3" & start_year>=1982 & start_year<=2020) %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))


## sb ps
sb8ps <- ccf(sant_prim_8$avg_slope, sant_secon_8$avg_slope, lag.max = 4)
sb7ps <- ccf(sant_prim_7$avg_slope, sant_secon_7$avg_slope, lag.max = 4)
sb6ps <- ccf(sant_prim_6$avg_slope, sant_secon_6$avg_slope, lag.max = 4)
sb5ps <- ccf(sant_prim_5$avg_slope, sant_secon_5$avg_slope, lag.max = 4)
sb4ps <- ccf(sant_prim_4$avg_slope, sant_secon_4$avg_slope, lag.max = 4)
sb3ps <- ccf(sant_prim_3$avg_slope, sant_secon_3$avg_slope, lag.max = 4)
##sb st
sb8st <- ccf(sant_secon_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4)
sb7st <- ccf(sant_secon_7$avg_slope, sant_tert_7$avg_slope, lag.max = 4)
sb6st <- ccf(sant_secon_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4)
sb5st <- ccf(sant_secon_5$avg_slope, sant_tert_5$avg_slope, lag.max = 4)
sb4st <- ccf(sant_secon_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4)
sb3st <- ccf(sant_secon_3$avg_slope, sant_tert_3$avg_slope, lag.max = 4)
## sb tp
sb8tp <- ccf(sant_prim_8$avg_slope, sant_tert_8$avg_slope, lag.max = 4)
sb7tp <- ccf(sant_prim_7$avg_slope, sant_tert_7$avg_slope, lag.max = 4)
sb6tp <- ccf(sant_prim_6$avg_slope, sant_tert_6$avg_slope, lag.max = 4)
sb5tp <- ccf(sant_prim_5$avg_slope, sant_tert_5$avg_slope, lag.max = 4)
sb4tp <- ccf(sant_prim_4$avg_slope, sant_tert_4$avg_slope, lag.max = 4)
sb3tp <- ccf(sant_prim_3$avg_slope, sant_tert_3$avg_slope, lag.max = 4)
## hb ps
hb8ps <- ccf(hub_prim_8$avg_slope, hub_secon_8$avg_slope, lag.max = 4)
hb7ps <- ccf(hub_prim_7$avg_slope, hub_secon_7$avg_slope, lag.max = 4)
hb6ps <- ccf(hub_prim_6$avg_slope, hub_secon_6$avg_slope, lag.max = 4)
hb5ps <- ccf(hub_prim_5$avg_slope, hub_secon_5$avg_slope, lag.max = 4)
hb4ps <- ccf(hub_prim_4$avg_slope, hub_secon_4$avg_slope, lag.max = 4)
hb3ps <- ccf(hub_prim_3$avg_slope, hub_secon_3$avg_slope, lag.max = 4)
## hb st
hb8st <- ccf(hub_secon_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4)
hb7st <- ccf(hub_secon_7$avg_slope, hub_tert_7$avg_slope, lag.max = 4)
hb6st <- ccf(hub_secon_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4)
hb5st <- ccf(hub_secon_5$avg_slope, hub_tert_5$avg_slope, lag.max = 4)
hb4st <- ccf(hub_secon_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4)
hb3st <- ccf(hub_secon_3$avg_slope, hub_tert_3$avg_slope, lag.max = 4)
## hb tp
hb8tp <- ccf(hub_prim_8$avg_slope, hub_tert_8$avg_slope, lag.max = 4)
hb7tp <- ccf(hub_prim_7$avg_slope, hub_tert_7$avg_slope, lag.max = 4)
hb6tp <- ccf(hub_prim_6$avg_slope, hub_tert_6$avg_slope, lag.max = 4)
hb5tp <- ccf(hub_prim_5$avg_slope, hub_tert_5$avg_slope, lag.max = 4)
hb4tp <- ccf(hub_prim_4$avg_slope, hub_tert_4$avg_slope, lag.max = 4)
hb3tp <- ccf(hub_prim_3$avg_slope, hub_tert_3$avg_slope, lag.max = 4)
## kz ps
kp8ps <- ccf(konzkat_prim_8$avg_slope, konzkat_secon_8$avg_slope, lag.max = 4)
kp7ps <- ccf(konzkat_prim_7$avg_slope, konzkat_secon_7$avg_slope, lag.max = 4)
kp6ps <- ccf(konzkat_prim_6$avg_slope, konzkat_secon_6$avg_slope, lag.max = 4)
kp5ps <- ccf(konzkat_prim_5$avg_slope, konzkat_secon_5$avg_slope, lag.max = 4)
kp4ps <- ccf(konzkat_prim_4$avg_slope, konzkat_secon_4$avg_slope, lag.max = 4)
kp3ps <- ccf(konzkat_prim_3$avg_slope, konzkat_secon_3$avg_slope, lag.max = 4)
#kz st
kp8st <- ccf(konzkat_secon_8$avg_slope, konzkat_tert_8$avg_slope, lag.max = 4)
kp7st <- ccf(konzkat_secon_7$avg_slope, konzkat_tert_7$avg_slope, lag.max = 4)
kp6st <- ccf(konzkat_secon_6$avg_slope, konzkat_tert_6$avg_slope, lag.max = 4)
kp5st <- ccf(konzkat_secon_5$avg_slope, konzkat_tert_5$avg_slope, lag.max = 4)
kp4st <- ccf(konzkat_secon_4$avg_slope, konzkat_tert_4$avg_slope, lag.max = 4)
kp3st <- ccf(konzkat_secon_3$avg_slope, konzkat_tert_3$avg_slope, lag.max = 4)
## kz tp
kp8tp <- ccf(konzkat_prim_8$avg_slope, konzkat_tert_8$avg_slope, lag.max = 4)
kp7tp <- ccf(konzkat_prim_7$avg_slope, konzkat_tert_7$avg_slope, lag.max = 4)
kp6tp <- ccf(konzkat_prim_6$avg_slope, konzkat_tert_6$avg_slope, lag.max = 4)
kp5tp <- ccf(konzkat_prim_5$avg_slope, konzkat_tert_5$avg_slope, lag.max = 4)
kp4tp <- ccf(konzkat_prim_4$avg_slope, konzkat_tert_4$avg_slope, lag.max = 4)
kp3tp <- ccf(konzkat_prim_3$avg_slope, konzkat_tert_3$avg_slope, lag.max = 4)

## cool. looks the same. Exported great.

SB8yrPS <- data.frame(lag=sb8ps$lag,CCF=sb8ps$acf)
SB8yrST <- data.frame(lag=sb8st$lag,CCF=sb8st$acf)
SB8yrTP <- data.frame(lag=sb8tp$lag,CCF=sb8tp$acf)
HB8yrPS <- data.frame(lag=hb8ps$lag,CCF=hb8ps$acf)
HB8yrST <- data.frame(lag=hb8st$lag,CCF=hb8st$acf)
HB8yrTP <- data.frame(lag=hb8tp$lag,CCF=hb8tp$acf)
KP8yrPS <- data.frame(lag=kp8ps$lag,CCF=kp8ps$acf)
KP8yrST <- data.frame(lag=kp8st$lag,CCF=kp8st$acf)
KP8yrTP <- data.frame(lag=kp8tp$lag,CCF=kp8tp$acf)

SB7yrPS <- data.frame(lag=sb7ps$lag,CCF=sb7ps$acf)
SB7yrST <- data.frame(lag=sb7st$lag,CCF=sb7st$acf)
SB7yrTP <- data.frame(lag=sb7tp$lag,CCF=sb7tp$acf)
HB7yrPS <- data.frame(lag=hb7ps$lag,CCF=hb7ps$acf)
HB7yrST <- data.frame(lag=hb7st$lag,CCF=hb7st$acf)
HB7yrTP <- data.frame(lag=hb7tp$lag,CCF=hb7tp$acf)
KP7yrPS <- data.frame(lag=kp7ps$lag,CCF=kp7ps$acf)
KP7yrST <- data.frame(lag=kp7st$lag,CCF=kp7st$acf)
KP7yrTP <- data.frame(lag=kp7tp$lag,CCF=kp7tp$acf)

SB6yrPS <- data.frame(lag=sb6ps$lag,CCF=sb6ps$acf)
SB6yrST <- data.frame(lag=sb6st$lag,CCF=sb6st$acf)
SB6yrTP <- data.frame(lag=sb6tp$lag,CCF=sb6tp$acf)
HB6yrPS <- data.frame(lag=hb6ps$lag,CCF=hb6ps$acf)
HB6yrST <- data.frame(lag=hb6st$lag,CCF=hb6st$acf)
HB6yrTP <- data.frame(lag=hb6tp$lag,CCF=hb6tp$acf)
KP6yrPS <- data.frame(lag=kp6ps$lag,CCF=kp6ps$acf)
KP6yrST <- data.frame(lag=kp6st$lag,CCF=kp6st$acf)
KP6yrTP <- data.frame(lag=kp6tp$lag,CCF=kp6tp$acf)

SB5yrPS <- data.frame(lag=sb5ps$lag,CCF=sb5ps$acf)
SB5yrST <- data.frame(lag=sb5st$lag,CCF=sb5st$acf)
SB5yrTP <- data.frame(lag=sb5tp$lag,CCF=sb5tp$acf)
HB5yrPS <- data.frame(lag=hb5ps$lag,CCF=hb5ps$acf)
HB5yrST <- data.frame(lag=hb5st$lag,CCF=hb5st$acf)
HB5yrTP <- data.frame(lag=hb5tp$lag,CCF=hb5tp$acf)
KP5yrPS <- data.frame(lag=kp5ps$lag,CCF=kp5ps$acf)
KP5yrST <- data.frame(lag=kp5st$lag,CCF=kp5st$acf)
KP5yrTP <- data.frame(lag=kp5tp$lag,CCF=kp5tp$acf)

SB4yrPS <- data.frame(lag=sb4ps$lag,CCF=sb4ps$acf)
SB4yrST <- data.frame(lag=sb4st$lag,CCF=sb4st$acf)
SB4yrTP <- data.frame(lag=sb4tp$lag,CCF=sb4tp$acf)
HB4yrPS <- data.frame(lag=hb4ps$lag,CCF=hb4ps$acf)
HB4yrST <- data.frame(lag=hb4st$lag,CCF=hb4st$acf)
HB4yrTP <- data.frame(lag=hb4tp$lag,CCF=hb4tp$acf)
KP4yrPS <- data.frame(lag=kp4ps$lag,CCF=kp4ps$acf)
KP4yrST <- data.frame(lag=kp4st$lag,CCF=kp4st$acf)
KP4yrTP <- data.frame(lag=kp4tp$lag,CCF=kp4tp$acf)

SB3yrPS <- data.frame(lag=sb3ps$lag,CCF=sb3ps$acf)
#, site="SB", trophCombo="PS", year=3)
SB3yrST <- data.frame(lag=sb3st$lag,CCF=sb3st$acf)
SB3yrTP <- data.frame(lag=sb3tp$lag,CCF=sb3tp$acf)
HB3yrPS <- data.frame(lag=hb3ps$lag,CCF=hb3ps$acf)
HB3yrST <- data.frame(lag=hb3st$lag,CCF=hb3st$acf)
HB3yrTP <- data.frame(lag=hb3tp$lag,CCF=hb3tp$acf)
KP3yrPS <- data.frame(lag=kp3ps$lag,CCF=kp3ps$acf)
KP3yrST <- data.frame(lag=kp3st$lag,CCF=kp3st$acf)
KP3yrTP <- data.frame(lag=kp3tp$lag,CCF=kp3tp$acf)

merge(SB8yrPS,SB8yrST)


sine_model <- function(x, A, B, C, D) {
  A * sin(B * x + C) + D
}
library(nls2)
# Fit the model to data

library(ggplot2)
library(dplyr)

fit_and_plot_model <- function(data, formula, start_list, plot_title, dataset_name) {
  fit <- nls(formula, data = data, start = start_list, control= nls.control(maxiter = 500))
  data$predicted <- predict(fit, newdata = data)
  maxcorr <- max(abs(data$CCF))
  mincorr <- min(abs(data$CCF))
  meancorr <-mean(abs(data$CCF))
  bestlag <- abs(data$lag[which(abs(data$CCF)==maxcorr)])
  rmse <- sqrt(mean((data$CCF - data$predicted)^2))
  ss_total <- sum((data$CCF - mean(data$CCF))^2)
  ss_res <- sum((data$CCF - data$predicted)^2)
  r_squared <- 1 - (ss_res / ss_total)
  metrics_df <- data.frame(Name = dataset_name, RMSE = rmse, R_squared = r_squared, 
                           bestlag, maxcorr=maxcorr, mincorr=mincorr, meancorr=meancorr)
  metrics_df$Year <- str_sub(metrics_df$Name, 3, 3)
  metrics_df$Site <- str_sub(metrics_df$Name, 1, 2)
  metrics_df$Combo <- str_sub(metrics_df$Name, 6, 7)
  plot <- ggplot(data, aes(x = lag, y = CCF)) +
    geom_point(color = "blue") +       # Original data points
    geom_line(aes(y = predicted), color = "red") +   # Fitted sine curve
    labs(x = "Lag", y = "CCF", title = plot_title)
  list(plot = plot, metrics = metrics_df)
}

### SB PS
SB8yrPSplot <- fit_and_plot_model(SB8yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                   start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                   plot_title = "Santa Barbara, Prim/secon, 8 years",
                                  dataset_name = "SB8yrPS")
print(SB8yrPSplot$plot)
print(SB8yrPSplot$metrics)
SB8yrPSplot$metrics
SB7yrPSplot <- fit_and_plot_model(SB7yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                   start_list = list(A = 0.5, B = 1, C = 1.3, D = 0.05),
                                   plot_title = "Santa Barbara, Prim/secon, 7 years",
                                  dataset_name = "SB7yrPS")
print(SB7yrPSplot$plot)
SB6yrPSplot <- fit_and_plot_model(SB6yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                   start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                   plot_title = "Santa Barbara, Prim/secon, 6 years",
                                  dataset_name = "SB6yrPS")
print(SB6yrPSplot$plot)
SB5yrPSplot <- fit_and_plot_model(SB5yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                   start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                   plot_title = "Santa Barbara, Prim/secon, 5 years",
                                  dataset_name = "SB5yrPS")
print(SB5yrPSplot$plot)
SB4yrPSplot <- fit_and_plot_model(SB4yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                   start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                   plot_title = "Santa Barbara, Prim/secon, 4 years",
                                  dataset_name = "SB4yrPS")
print(SB4yrPSplot$plot)

## HB PS
HB8yrPSplot <- fit_and_plot_model(HB8yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/secon, 8 years",
                                  dataset_name = "HB8yrPS")
print(HB8yrPSplot$plot)
HB7yrPSplot <- fit_and_plot_model(HB7yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/secon, 7 years",
                                  dataset_name = "HB7yrPS")
print(HB7yrPSplot$plot)
HB6yrPSplot <- fit_and_plot_model(HB6yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/secon, 6 years",
                                  dataset_name = "HB6yrPS")
print(HB6yrPSplot$plot)
HB5yrPSplot <- fit_and_plot_model(HB5yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/secon, 5 years",
                                  dataset_name = "HB5yrPS")
print(HB5yrPSplot$plot)
HB4yrPSplot <- fit_and_plot_model(HB4yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/secon, 4 years",
                                  dataset_name = "HB4yrPS")

print(HB4yrPSplot$plot)
## KP PS
KP8yrPSplot <- fit_and_plot_model(KP8yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                  plot_title = "Konza Prairie, Prim/secon, 8 years",
                                  dataset_name = "KP8yrPS")
print(KP8yrPSplot$plot)
KP7yrPSplot <- fit_and_plot_model(KP7yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.6, B = 100, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/secon, 7 years",
                                  dataset_name = "KP7yrPS")
print(KP7yrPSplot$plot)
KP6yrPSplot <- fit_and_plot_model(KP6yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 100, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/secon, 6 years",
                                  dataset_name = "KP6yrPS")
print(KP6yrPSplot$plot)
KP5yrPSplot <- fit_and_plot_model(KP5yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 100, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/secon, 5 years",
                                  dataset_name = "KP5yrPS")
print(KP5yrPSplot$plot)
KP4yrPSplot <- fit_and_plot_model(KP4yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1.0, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/secon, 4 years",
                                  dataset_name = "KP4yrPS")
print(KP4yrPSplot$plot)

### SB ST
SB8yrSTplot <- fit_and_plot_model(SB8yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Secon/Tert, 8 years",
                                  dataset_name = "SB8yrST")
print(SB8yrSTplot$plot)
print(SB8yrSTplot$metrics)
SB7yrSTplot <- fit_and_plot_model(SB7yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.1, B = 0.5, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Secon/Tert, 7 years",
                                  dataset_name = "SB7yrST")
print(SB7yrSTplot$plot)
SB6yrSTplot <- fit_and_plot_model(SB6yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Secon/Tert, 6 years",
                                  dataset_name = "SB6yrST")
print(SB6yrSTplot$plot)
SB5yrSTplot <- fit_and_plot_model(SB5yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Secon/Tert, 5 years",
                                  dataset_name = "SB5yrST")
print(SB5yrSTplot$plot)
SB4yrSTplot <- fit_and_plot_model(SB4yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Secon/Tert, 4 years",
                                  dataset_name = "SB4yrST")
print(SB4yrSTplot$plot)
## HB ST
HB8yrSTplot <- fit_and_plot_model(HB8yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 100, C = 2, D = 0),
                                  plot_title = "Hubbard Brook, Secon/Tert, 8 years",
                                  dataset_name = "HB8yrST")
print(HB8yrSTplot$plot)
HB7yrSTplot <- fit_and_plot_model(HB7yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 100, C = 2, D = 0),
                                  plot_title = "Hubbard Brook, Secon/Tert, 7 years",
                                  dataset_name = "HB7yrST")
print(HB7yrSTplot$plot)
HB6yrSTplot <- fit_and_plot_model(HB6yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.6, B = .75, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Secon/Tert, 6 years",
                                  dataset_name = "HB6yrST")
print(HB6yrSTplot$plot)
HB5yrSTplot <- fit_and_plot_model(HB5yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.6, B = 0.75, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Secon/Tert, 5 years",
                                  dataset_name = "HB5yrST")
print(HB5yrSTplot$plot)
HB4yrSTplot <- fit_and_plot_model(HB4yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  # alternatively, and more frequent, 0.25, B = 1.5,
                                  plot_title = "Hubbard Brook, Secon/Tert, 4 years",
                                  dataset_name = "HB4yrST")

print(HB4yrSTplot$plot)
## KP ST
KP8yrSTplot <- fit_and_plot_model(KP8yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Secon/Tert, 8 years",
                                  dataset_name = "KP8yrST")
print(KP8yrSTplot$plot)
KP7yrSTplot <- fit_and_plot_model(KP7yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.1, B = 0.5, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Secon/Tert, 7 years",
                                  dataset_name = "KP7yrST")
print(KP7yrSTplot$plot)
KP6yrSTplot <- fit_and_plot_model(KP6yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.4, B = 0.5, C = 1, D = 0.05),
                                  plot_title = "Konza Prairie, Secon/Tert, 6 years",
                                  dataset_name = "KP6yrST")
print(KP6yrSTplot$plot)
KP5yrSTplot <- fit_and_plot_model(KP5yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Secon/Tert, 5 years",
                                  dataset_name = "KP5yrST")
print(KP5yrSTplot$plot)
KP4yrSTplot <- fit_and_plot_model(KP4yrST, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Secon/Tert, 4 years",
                                  dataset_name = "KP4yrST")
print(KP4yrSTplot$plot)
###### TP
##HB TP
SB8yrTPplot <- fit_and_plot_model(SB8yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Prim/tert, 8 years",
                                  dataset_name = "SB8yrTP")
print(SB8yrTPplot$plot)
SB7yrTPplot <- fit_and_plot_model(SB7yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Prim/tert, 7 years",
                                  dataset_name = "SB7yrTP")
print(SB7yrTPplot$plot)
SB6yrTPplot <- fit_and_plot_model(SB6yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Prim/tert, 6 years",
                                  dataset_name = "SB6yrTP")
print(SB6yrTPplot$plot)
SB5yrTPplot <- fit_and_plot_model(SB5yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Prim/tert, 5 years",
                                  dataset_name = "SB5yrTP")
print(SB5yrTPplot$plot)
SB4yrTPplot <- fit_and_plot_model(SB4yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1.2, C = 0, D = 0),
                                  plot_title = "Santa Barbara, Prim/tert, 4 years",
                                  dataset_name = "SB4yrTP")
print(SB4yrTPplot$plot)

## HB TP
HB8yrTPplot <- fit_and_plot_model(HB8yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/tert, 8 years",
                                  dataset_name = "HB8yrTP")
print(HB8yrTPplot$plot)
HB7yrTPplot <- fit_and_plot_model(HB7yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.1, B = 0.5, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/tert, 7 years",
                                  dataset_name = "HB7yrTP")
print(HB7yrTPplot$plot)
HB6yrTPplot <- fit_and_plot_model(HB6yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.75, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/tert, 6 years",
                                  dataset_name = "HB6yrTP")
print(HB6yrTPplot$plot)
HB5yrTPplot <- fit_and_plot_model(HB5yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.75, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/tert, 5 years",
                                  dataset_name = "HB5yrTP")
print(HB5yrTPplot$plot)
HB4yrTPplot <- fit_and_plot_model(HB4yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.75, C = 0, D = 0),
                                  plot_title = "Hubbard Brook, Prim/tert, 4 years",
                                  dataset_name = "HB4yrTP")

print(HB4yrTPplot$plot)
## KP TP
KP8yrTPplot <- fit_and_plot_model(KP8yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/tert, 8 years",
                                  dataset_name = "KP8yrTP")
print(KP8yrTPplot$plot)
KP7yrTPplot <- fit_and_plot_model(KP7yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.1, B = 0.5, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/tert, 7 years",
                                  dataset_name = "KP7yrTP")
print(KP7yrTPplot$plot)
KP6yrTPplot <- fit_and_plot_model(KP6yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/tert, 6 years",
                                  dataset_name = "KP6yrTP")
print(KP6yrTPplot$plot)
print(KP6yrTPplot$metrics)
KP5yrTPplot <- fit_and_plot_model(KP5yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 0.5, C = 0, D = 0),
                                  plot_title = "Konza Prairie, Prim/tert, 5 years",
                                  dataset_name = "KP5yrTP")
print(KP5yrTPplot$plot)
KP4yrTPplot <- fit_and_plot_model(KP4yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.45, B = 1, C = -4, D = 0.05),
                                  plot_title = "Konza Prairie, Prim/tert, 4 years",
                                  dataset_name = "KP4yrTP")
print(KP4yrTPplot$plot)



library(ggpubr)
mfrow=c(3,3)
ggarrange(SB8yrPSplot$plot,SB8yrSTplot$plot,SB8yrTPplot$plot,HB8yrPSplot$plot,HB8yrSTplot$plot,
          HB8yrTPplot$plot,KP8yrPSplot$plot,KP8yrSTplot$plot, KP8yrTPplot$plot,  ncol = 3, nrow =3)

ggarrange(SB7yrPSplot$plot,SB7yrSTplot$plot,SB7yrTPplot$plot,HB7yrPSplot$plot,HB7yrSTplot$plot,
          HB7yrTPplot$plot,KP7yrPSplot$plot,KP7yrSTplot$plot, KP7yrTPplot$plot,  ncol = 3, nrow =3)

ggarrange(SB6yrPSplot$plot,SB6yrSTplot$plot,SB6yrTPplot$plot,HB6yrPSplot$plot,HB6yrSTplot$plot,
          HB6yrTPplot$plot,KP6yrPSplot$plot,KP6yrSTplot$plot, KP6yrTPplot$plot,  ncol = 3, nrow =3)

ggarrange(SB5yrPSplot$plot,SB5yrSTplot$plot,SB5yrTPplot$plot,HB5yrPSplot$plot,HB5yrSTplot$plot,
          HB5yrTPplot$plot,KP5yrPSplot$plot,KP5yrSTplot$plot, KP5yrTPplot$plot,  ncol = 3, nrow =3)

ggarrange(SB4yrPSplot$plot,SB4yrSTplot$plot,SB4yrTPplot$plot,HB4yrPSplot$plot,HB4yrSTplot$plot,
          HB4yrTPplot$plot,KP4yrPSplot$plot,KP4yrSTplot$plot, KP4yrTPplot$plot,  ncol = 3, nrow =3)


#not run, using minimum overlap of 4 years
# ggarrange(SB3yrPSplot$plot,SB3yrSTplot$plot,SB3yrTPplot,HB3yrPSplot,HB3yrSTplot,
#           HB3yrTPplot$plot,KP3yrPSplot$plot,KP3yrSTplot, KP3yrTPplot,  ncol = 3, nrow =3)

df <- bind_rows(SB4yrPSplot$metrics,SB4yrSTplot$metrics,SB4yrTPplot$metrics,HB4yrPSplot$metrics,HB4yrSTplot$metrics,
                HB4yrTPplot$metrics,KP4yrPSplot$metrics,KP4yrSTplot$metrics, KP4yrTPplot$metrics,SB5yrPSplot$metrics,
                SB5yrSTplot$metrics,SB5yrTPplot$metrics,HB5yrPSplot$metrics,HB5yrSTplot$metrics,
                HB5yrTPplot$metrics,KP5yrPSplot$metrics,KP5yrSTplot$metrics, KP5yrTPplot$metrics,SB6yrPSplot$metrics,
                SB6yrSTplot$metrics,SB6yrTPplot$metrics,HB6yrPSplot$metrics,HB6yrSTplot$metrics,
                HB6yrTPplot$metrics,KP6yrPSplot$metrics,KP6yrSTplot$metrics, KP6yrTPplot$metrics, SB7yrPSplot$metrics,
                SB7yrSTplot$metrics,SB7yrTPplot$metrics,HB7yrPSplot$metrics,HB7yrSTplot$metrics,
                HB7yrTPplot$metrics,KP7yrPSplot$metrics,KP7yrSTplot$metrics, KP7yrTPplot$metrics, SB8yrPSplot$metrics,
                SB8yrSTplot$metrics,SB8yrTPplot$metrics,HB8yrPSplot$metrics,HB8yrSTplot$metrics,
                HB8yrTPplot$metrics,KP8yrPSplot$metrics,KP8yrSTplot$metrics, KP8yrTPplot$metrics)
df
library(stringr)

write.csv(df, file="SineFit.csv")

ggplot(df, aes(as.numeric(Year), R_squared, color=as.factor(Site)))+
  geom_point()+geom_smooth(method="lm", se=F)+facet_wrap(~bestlag)


PrimSecon <- df %>%
  filter(Combo=="PS") %>%
  pivot_longer(cols = c(R_squared, RMSE), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x=as.numeric(Year), y=Value, color = Metric))+
  geom_point()+geom_smooth(method="lm", se=F)+
  facet_grid(~Site)
PrimSecon       

SeconTert <- df %>%
  filter(Combo=="ST") %>%
  pivot_longer(cols = c(R_squared, RMSE), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x=as.numeric(Year), y=Value, color = Metric))+
  geom_point()+geom_smooth(method="lm", se=F)+
  facet_grid(~Site)
SeconTert  

PrimTert <- df %>%
  filter(Combo=="TP") %>%
  pivot_longer(cols = c(R_squared, RMSE), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x=as.numeric(Year), y=Value, color = Metric))+
  geom_point()+geom_smooth(method="lm", se=F)+
  facet_grid(~Site)
PrimTert  


library(ggplot2)
library(ggpmisc)
library(tidyr)
df$Year <- as.numeric(df$Year)
All <- df %>%
  pivot_longer(cols = c(R_squared, RMSE, maxcorr, mincorr), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x=Year, y=Value, color = Metric))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE) + # Add this line for the trend lines
  facet_grid(Site ~ Combo, labeller = labeller(
    Site = c(SB = "Santa Barbara", HB = "Hubbard Brook", KP = "Konza Prairie"), # Example renaming for SITE
    Combo = c(PS = "Primary & Secondary", TP = "Primary & Tertiary", ST = "Secondary and Tertiary") # Example renaming for COMBO
  )) # Facet by SITE across and COMBO down with custom labels
All 
library(gridExtra)
library(purrr)
library(broom)

models <- df %>%
  pivot_longer(cols = c(R_squared, RMSE), names_to = "Metric", values_to = "Value") %>%
  group_by(Site, Combo, Metric) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(Value ~ Year, data = .x)))

# Now to extract the R-squared values and tidy the models dataframe
model_summaries <- models %>%
  mutate(glance = map(model, glance)) %>%
  unnest(glance) %>%
  select(Site, Combo, Metric, r.squared)

print(model_summaries)

library(gt)
gt_table <- gt(model_summaries) %>%
  tab_header(
    title = "Model Summaries",
    subtitle = "R-squared values for various combinations of Site, Combo, and Metric"
  ) %>%
  cols_label(
    Site = "Site",
    Combo = "Combination",
    Metric = "Metric",
    r.squared = "R-squared"
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_options(
    heading.background.color = "#D3D3D3",
    table.font.size = "small"
  )

# Print the gt table
print(gt_table)



