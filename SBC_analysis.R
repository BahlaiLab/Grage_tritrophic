# analysis for sbc site
#do analysis and create figures for each site then bring the panels for each analysis piece together in another file
#let's make our code like russian dolls


#bring in data
#for sbc we will select algae in carp lot as as primary, mobile inverts in carp for secondary,
#and predatory fish in carp for tertiary


sbcdata.primary<-read.csv(file="cleaned_data/SBC_producer_algae_carp.csv")
sbcdata.secondary<-read.csv(file="cleaned_data/SBC_consumer_minvert_carp.csv")
sbcdata.tertiary<-read.csv(file="cleaned_data/SBC_predator_fish_carp.csv")

#issue with column names
names(sbcdata.primary)[names(sbcdata.primary) == "biomass"] <- "Abundance"
names(sbcdata.secondary)[names(sbcdata.secondary) == "biomass"] <- "Abundance"
names(sbcdata.tertiary)[names(sbcdata.tertiary) == "biomass"] <- "Abundance"
names(sbcdata.primary)[names(sbcdata.primary) == "YEAR"] <- "Year"
names(sbcdata.secondary)[names(sbcdata.secondary) == "YEAR"] <- "Year"
names(sbcdata.tertiary)[names(sbcdata.tertiary) == "YEAR"] <- "Year"

#merge three trophic levels into a single frame

library(dplyr)


sbcdata <- bind_rows(
  sbcdata.primary %>%
    mutate(Trophic.level = "Primary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE)),
  
  sbcdata.secondary %>%
    mutate(Trophic.level = "Secondary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE)),
  
  sbcdata.tertiary %>%
    mutate(Trophic.level = "Tertiary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE))
)
#create a timeseries plot and annotate it

library(ggplot2)
library(grid)

#summarize data so we don't have so many points on the plot
sbc_summary <- sbcdata %>%
  group_by(Year, Trophic.level) %>%
  summarize(meanAbundance = mean(normAbundance, na.rm = TRUE), .groups = "drop")

timeseriesplot.sbcdata <- ggplot(sbc_summary, aes(
  x = Year, 
  y = meanAbundance, 
  fill = Trophic.level, 
  color = Trophic.level, 
  shape = Trophic.level, 
  linetype = Trophic.level
)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 3) +
  scale_color_manual(values = c("tan4", "orange3", "firebrick4"), name = "Trophic level") +
  scale_fill_manual(values = c("tan", "orange", "firebrick1"), name = "Trophic level") +
  scale_shape_manual(values = c(23, 22, 21), name = "Trophic level") +
  scale_linetype_manual(values = c("solid", "solid", "solid"), name = "Trophic level") +
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.85)
  ) +
  ylab("Normalized abundance")

timeseriesplot.sbcdata



#run bad_breakup on each trophic level
#####################################
#load bad breakup script
source_github <- function(u) {
  # load package
  require(RCurl)
  
  # read script lines from website
  script <- getURL(u, ssl.verifypeer = FALSE)
  
  # parse lines and evaluate in the global environment
  eval(parse(text = script))
}

source("https://raw.githubusercontent.com/BahlaiLab/bad_breakup_2/master/R_model/bad_breakup_script.R")

options(scipen=10)

model.sbcdata.primary <- multiple_breakups(sbcdata.primary)
model.sbcdata.secondary <- multiple_breakups(sbcdata.secondary)
model.sbcdata.tertiary <- multiple_breakups(sbcdata.tertiary)

#and for the record, what's the stability time for each trophic level? 
stability_time(sbcdata.primary) #8
stability_time(sbcdata.secondary)#10
stability_time(sbcdata.tertiary) #6

#and for the record, what's the proportion wih a misleading trajectory before stability for each trophic level? 
proportion_wrong_before_stability(sbcdata.primary, significance=0.05) #85%
proportion_wrong_before_stability(sbcdata.secondary, significance=0.05)#67%
proportion_wrong_before_stability(sbcdata.tertiary, significance=0.05) #7%

linefit(standardize(sbcdata.primary)) #startyear Ndata Nyears slope slopeSE, slopeP etc
linefit(standardize(sbcdata.secondary))
linefit(standardize(sbcdata.tertiary))

#combine data back into a single frame

model.sbcdata.primary$site <- rep(c("sbc"),each = 190)
model.sbcdata.primary$trophic_level <- rep(c("Primary"),each = 190)

model.sbcdata.secondary$site <- rep(c("sbc"),each = 190)
model.sbcdata.secondary$trophic_level <- rep(c("Secondary"),each = 190)

model.sbcdata.tertiary$site <- rep(c("sbc"),each = 190)
model.sbcdata.tertiary$trophic_level <- rep(c("Tertiary"),each = 190)

#now merge all dataframes together
model.sbcdata <- rbind(model.sbcdata.primary, model.sbcdata.secondary, model.sbcdata.tertiary)

library(forcats)
library(stringr)

#create heatmap for each study duration

#use only 4, 7, 10y study durations because 10 year max stability time 
all_heat.sbcdata <- model.sbcdata %>%
  filter(N_years == "4"|N_years == "7"|N_years == "10") %>%
  ggplot(aes(x=start_year, y=trophic_level)) + geom_tile(aes(fill = as.numeric(slope)),colour = "white") +
  scale_fill_distiller(palette="RdBu") + 
  #scale_fill_viridis(option="mako") + 
  guides(fill = guide_colorbar(title = "Slope")) +
  labs(y = "Trophic level") +
  scale_x_continuous( name = "Start year",
                      sec.axis = sec_axis(~ ., name = "Study duration"))+
  theme_bw() + theme_minimal()+
  theme(text = element_text(), axis.text.x.top =  element_blank(), 
        panel.spacing = unit(2, "lines"))+ 
  facet_grid(cols=vars(N_years))
all_heat.sbcdata

library(svglite)


#perform cross-correlation analysis
#subset data
##prim

sbc_prim_10 <- model.sbcdata[model.sbcdata$trophic_level == 'Primary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_prim_9 <- model.sbcdata[model.sbcdata$trophic_level == 'Primary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_prim_8 <- model.sbcdata[model.sbcdata$trophic_level == 'Primary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_prim_7 <- model.sbcdata[model.sbcdata$trophic_level == 'Primary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_prim_6 <- model.sbcdata[model.sbcdata$trophic_level == 'Primary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_prim_5 <- model.sbcdata[model.sbcdata$trophic_level == 'Primary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_prim_4 <- model.sbcdata[model.sbcdata$trophic_level == 'Primary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

##sec

sbc_sec_10 <- model.sbcdata[model.sbcdata$trophic_level == 'Secondary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_sec_9 <- model.sbcdata[model.sbcdata$trophic_level == 'Secondary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_sec_8 <- model.sbcdata[model.sbcdata$trophic_level == 'Secondary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_sec_7 <- model.sbcdata[model.sbcdata$trophic_level == 'Secondary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_sec_6 <- model.sbcdata[model.sbcdata$trophic_level == 'Secondary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_sec_5 <- model.sbcdata[model.sbcdata$trophic_level == 'Secondary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_sec_4 <- model.sbcdata[model.sbcdata$trophic_level == 'Secondary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))


##tert

sbc_tert_10 <- model.sbcdata[model.sbcdata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_tert_9 <- model.sbcdata[model.sbcdata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_tert_8 <- model.sbcdata[model.sbcdata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_tert_7 <- model.sbcdata[model.sbcdata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_tert_6 <- model.sbcdata[model.sbcdata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_tert_5 <- model.sbcdata[model.sbcdata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
sbc_tert_4 <- model.sbcdata[model.sbcdata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

#Now run the cross correlation

sbc10ps <- ccf(sbc_prim_10$avg_slope, sbc_sec_10$avg_slope, lag.max = 4)
sbc9ps <- ccf(sbc_prim_9$avg_slope, sbc_sec_9$avg_slope, lag.max = 4)
sbc8ps <- ccf(sbc_prim_8$avg_slope, sbc_sec_8$avg_slope, lag.max = 4)
sbc7ps <- ccf(sbc_prim_7$avg_slope, sbc_sec_7$avg_slope, lag.max = 4)
sbc6ps <- ccf(sbc_prim_6$avg_slope, sbc_sec_6$avg_slope, lag.max = 4)
sbc5ps <- ccf(sbc_prim_5$avg_slope, sbc_sec_5$avg_slope, lag.max = 4)
sbc4ps <- ccf(sbc_prim_4$avg_slope, sbc_sec_4$avg_slope, lag.max = 4)

##st

sbc10st <- ccf(sbc_sec_10$avg_slope, sbc_tert_10$avg_slope, lag.max = 4)
sbc9st <- ccf(sbc_sec_9$avg_slope, sbc_tert_9$avg_slope, lag.max = 4)
sbc8st <- ccf(sbc_sec_8$avg_slope, sbc_tert_8$avg_slope, lag.max = 4)
sbc7st <- ccf(sbc_sec_7$avg_slope, sbc_tert_7$avg_slope, lag.max = 4)
sbc6st <- ccf(sbc_sec_6$avg_slope, sbc_tert_6$avg_slope, lag.max = 4)
sbc5st <- ccf(sbc_sec_5$avg_slope, sbc_tert_5$avg_slope, lag.max = 4)
sbc4st <- ccf(sbc_sec_4$avg_slope, sbc_tert_4$avg_slope, lag.max = 4)

## tp

sbc10tp <- ccf(sbc_prim_10$avg_slope, sbc_tert_10$avg_slope, lag.max = 4)
sbc9tp <- ccf(sbc_prim_9$avg_slope, sbc_tert_9$avg_slope, lag.max = 4)
sbc8tp <- ccf(sbc_prim_8$avg_slope, sbc_tert_8$avg_slope, lag.max = 4)
sbc7tp <- ccf(sbc_prim_7$avg_slope, sbc_tert_7$avg_slope, lag.max = 4)
sbc6tp <- ccf(sbc_prim_6$avg_slope, sbc_tert_6$avg_slope, lag.max = 4)
sbc5tp <- ccf(sbc_prim_5$avg_slope, sbc_tert_5$avg_slope, lag.max = 4)
sbc4tp <- ccf(sbc_prim_4$avg_slope, sbc_tert_4$avg_slope, lag.max = 4)

#prepare data for plots

sbc10yrPS <- data.frame(lag=sbc10ps$lag,CCF=sbc10ps$acf)
sbc10yrST <- data.frame(lag=sbc10st$lag,CCF=sbc10st$acf)
sbc10yrTP <- data.frame(lag=sbc10tp$lag,CCF=sbc10tp$acf)

sbc9yrPS <- data.frame(lag=sbc9ps$lag,CCF=sbc9ps$acf)
sbc9yrST <- data.frame(lag=sbc9st$lag,CCF=sbc9st$acf)
sbc9yrTP <- data.frame(lag=sbc9tp$lag,CCF=sbc9tp$acf)

sbc8yrPS <- data.frame(lag=sbc8ps$lag,CCF=sbc8ps$acf)
sbc8yrST <- data.frame(lag=sbc8st$lag,CCF=sbc8st$acf)
sbc8yrTP <- data.frame(lag=sbc8tp$lag,CCF=sbc8tp$acf)

sbc7yrPS <- data.frame(lag=sbc7ps$lag,CCF=sbc7ps$acf)
sbc7yrST <- data.frame(lag=sbc7st$lag,CCF=sbc7st$acf)
sbc7yrTP <- data.frame(lag=sbc7tp$lag,CCF=sbc7tp$acf)

sbc6yrPS <- data.frame(lag=sbc6ps$lag,CCF=sbc6ps$acf)
sbc6yrST <- data.frame(lag=sbc6st$lag,CCF=sbc6st$acf)
sbc6yrTP <- data.frame(lag=sbc6tp$lag,CCF=sbc6tp$acf)

sbc5yrPS <- data.frame(lag=sbc5ps$lag,CCF=sbc5ps$acf)
sbc5yrST <- data.frame(lag=sbc5st$lag,CCF=sbc5st$acf)
sbc5yrTP <- data.frame(lag=sbc5tp$lag,CCF=sbc5tp$acf)

sbc4yrPS <- data.frame(lag=sbc4ps$lag,CCF=sbc4ps$acf)
sbc4yrST <- data.frame(lag=sbc4st$lag,CCF=sbc4st$acf)
sbc4yrTP <- data.frame(lag=sbc4tp$lag,CCF=sbc4tp$acf)

#now analyse and plot the cross correlation analysis
sine_model <- function(x, A, B, C, D) {
  A * sin(B * x + C) + D
}
library(minpack.lm)
# Fit the model to data
library(readr)


fit_and_plot_model <- function(data, formula, start_list, plot_title, dataset_name) {
  fit <- nlsLM(formula, data = data, start = start_list, control= nls.control(maxiter = 5000))
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
  metrics_df$Year <- parse_number(metrics_df$Name)
  metrics_df$Site <- str_sub(metrics_df$Name, 1, 2)
  metrics_df$Combo <- str_sub(metrics_df$Name, -2,-1)
  plot <- ggplot(data, aes(x = lag, y = CCF)) +
    geom_point(color = "blue") +       # Original data points
    geom_line(aes(y = predicted), color = "red") +   # Fitted sine curve
    labs(x = NULL, y =  NULL, title = NULL)+
    theme_bw()
  list(plot = plot, metrics = metrics_df)
}

### sbc PS



sbc10yrPSplot <- fit_and_plot_model(sbc10yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "sbc10yrPS")
#print(sbc10yrPSplot$plot)
print(sbc10yrPSplot$metrics)

sbc9yrPSplot <- fit_and_plot_model(sbc9yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc9yrPS")
#print(sbc9yrPSplot$plot)
print(sbc9yrPSplot$metrics)

sbc8yrPSplot <- fit_and_plot_model(sbc8yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc8yrPS")
#print(sbc8yrPSplot$plot)
print(sbc8yrPSplot$metrics)

sbc7yrPSplot <- fit_and_plot_model(sbc7yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc7yrPS")
#print(sbc7yrPSplot$plot)
print(sbc7yrPSplot$metrics)

sbc6yrPSplot <- fit_and_plot_model(sbc6yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc6yrPS")
#print(sbc6yrPSplot$plot)
print(sbc6yrPSplot$metrics)

sbc5yrPSplot <- fit_and_plot_model(sbc5yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc5yrPS")
#print(sbc5yrPSplot$plot)
print(sbc5yrPSplot$metrics)

sbc4yrPSplot <- fit_and_plot_model(sbc4yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc4yrPS")
#print(sbc4yrPSplot$plot)
print(sbc4yrPSplot$metrics)

### sbc ST


sbc10yrSTplot <- fit_and_plot_model(sbc10yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "sbc10yrST")
#print(sbc10yrSTplot$plot)
print(sbc10yrSTplot$metrics)

sbc9yrSTplot <- fit_and_plot_model(sbc9yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 2, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc9yrST")
#print(sbc9yrSTplot$plot)
print(sbc9yrSTplot$metrics)

sbc8yrSTplot <- fit_and_plot_model(sbc8yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc8yrST")
#print(sbc8yrSTplot$plot)
print(sbc8yrSTplot$metrics)

sbc7yrSTplot <- fit_and_plot_model(sbc7yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                     plot_title = "Sec/tert, 7 years",
                                     dataset_name = "sbc7yrST")
#print(sbc7yrSTplot$plot)
print(sbc7yrSTplot$metrics)

sbc6yrSTplot <- fit_and_plot_model(sbc6yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc6yrST")
#print(sbc6yrSTplot$plot)
print(sbc6yrSTplot$metrics)

sbc5yrSTplot <- fit_and_plot_model(sbc5yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc5yrST")
#print(sbc5yrSTplot$plot)
print(sbc5yrSTplot$metrics)

sbc4yrSTplot <- fit_and_plot_model(sbc4yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc4yrST")
#print(sbc4yrSTplot$plot)
print(sbc4yrSTplot$metrics)

### sbc TP


sbc10yrTPplot <- fit_and_plot_model(sbc10yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "sbc10yrTP")
#print(sbc10yrTPplot$plot)
print(sbc10yrTPplot$metrics)

sbc9yrTPplot <- fit_and_plot_model(sbc9yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc9yrTP")
#print(sbc9yrTPplot$plot)
print(sbc9yrTPplot$metrics)


sbc8yrTPplot <- fit_and_plot_model(sbc8yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc8yrTP")
#print(sbc8yrTPplot$plot)
print(sbc8yrTPplot$metrics)

sbc7yrTPplot <- fit_and_plot_model(sbc7yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc7yrTP")
#print(sbc7yrTPplot$plot)
print(sbc7yrTPplot$metrics)

sbc6yrTPplot <- fit_and_plot_model(sbc6yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc6yrTP")
#print(sbc6yrTPplot$plot)
print(sbc6yrTPplot$metrics)

sbc5yrTPplot <- fit_and_plot_model(sbc5yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc5yrTP")
#print(sbc5yrTPplot$plot)
print(sbc5yrTPplot$metrics)

sbc4yrTPplot <- fit_and_plot_model(sbc4yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 1, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "sbc4yrTP")
#print(sbc4yrTPplot$plot)
print(sbc4yrTPplot$metrics)

library(ggpubr)
library(cowplot)
# ggarrange(sbc8yrPSplot$plot,sbc8yrSTplot$plot,sbc8yrTPplot$plot, 
#           #sbc7yrPSplot$plot,sbc7yrSTplot$plot,sbc7yrTPplot$plot,
#           sbc6yrPSplot$plot,sbc6yrSTplot$plot,sbc6yrTPplot$plot,
#           #sbc5yrPSplot$plot,sbc5yrSTplot$plot,sbc5yrTPplot$plot,
#           sbc4yrPSplot$plot,sbc4yrSTplot$plot,sbc4yrTPplot$plot,
#           ncol = 3, nrow=3)

prim.sec<-text_grob("Primary\n and\n Secondary", size = 9)
sec.tert<-text_grob("Secondary\n and\n Tertiary", size = 9)
prim.tert<-text_grob("Primary\n and\n Tertiary", size = 9)

pssine<-plot_grid(sbc4yrPSplot$plot, sbc7yrPSplot$plot, sbc10yrPSplot$plot, prim.sec, 
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
stsine<-plot_grid(sbc4yrSTplot$plot, sbc7yrSTplot$plot, sbc10yrSTplot$plot, sec.tert,
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
ptsine<-plot_grid( sbc4yrTPplot$plot, sbc7yrTPplot$plot, sbc10yrTPplot$plot, prim.tert,
                   ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))

studyduration<-text_grob("Study duration", size=11)
duration4<-text_grob("   4", size=10)
durationmid<-text_grob("   7", size=10)
durationtop<-text_grob("   10", size=10)
blankspace<-text_grob("", size=10)

durationbar<-plot_grid(blankspace, duration4, durationmid, durationtop, blankspace,
                       ncol = 5, nrow=1, rel_widths=c(0.07,0.3,0.3,0.3,0.1))
durationlabel<-plot_grid(blankspace, studyduration, blankspace,ncol = 3, nrow=1, rel_widths=c(0.07,0.98,0.1))

plots_together<-plot_grid(pssine, stsine, ptsine, ncol = 1)

corlabel<-text_grob("Cross-correlation", size=12, rot=90)
laglabel<-text_grob("Years of lag", size=12)

plots_together_lab<-plot_grid(corlabel, plots_together, blankspace, laglabel,  ncol=2,
                              rel_heights=c(0.95, 0.07), rel_widths=c(0.05, 0.95) )

plots_together_lab


allccf.sbcdata<-plot_grid(durationlabel, durationbar, plots_together_lab,   ncol=1, rel_heights=c(0.1, 0.08, 1.9))

allccf.sbcdata

#compile fit metrics and make plots
sbcdf <- bind_rows(sbc10yrPSplot$metrics,sbc10yrSTplot$metrics,sbc10yrTPplot$metrics,
                     sbc9yrPSplot$metrics,sbc9yrSTplot$metrics,sbc9yrTPplot$metrics,
                     sbc8yrPSplot$metrics,sbc8yrSTplot$metrics,sbc8yrTPplot$metrics, 
                     sbc7yrPSplot$metrics,sbc7yrSTplot$metrics,sbc7yrTPplot$metrics,
                     sbc6yrPSplot$metrics,sbc6yrSTplot$metrics,sbc6yrTPplot$metrics,
                     sbc5yrPSplot$metrics,sbc5yrSTplot$metrics,sbc5yrTPplot$metrics,
                     sbc4yrPSplot$metrics,sbc4yrSTplot$metrics,sbc4yrTPplot$metrics)
sbcdf


library(ggpmisc)
library(tidyr)
sbcdf$Year <- as.numeric(sbcdf$Year)
allderived.sbcdata <- sbcdf %>%
  pivot_longer(cols = c(R_squared, RMSE, maxcorr), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x=Year, y=Value, color = Metric, shape=Metric, linetype=Metric))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE) + # Add this line for the trend lines
  facet_wrap(~ Combo, labeller = labeller(
    Combo = c(PS = "Primary and Secondary", TP = "Primary and Tertiary", ST = "Secondary and Tertiary") # Example renaming for COMBO
  ))+ # Facet by  COMBO with custom labels
  theme_classic()+
  scale_color_manual(values=c("#1b9e77", "#d95f02","#7570b3"), 
                     labels=c("Maximum\n correlation", expression(paste("R"^2)), "RMSE"))+
  scale_shape_manual(values=c(5, 6, 7),
                     labels=c("Maximum\n correlation", expression(paste("R"^2)), "RMSE"))+
  scale_linetype_manual(values=c("dashed", "solid", "dotted"),
                        labels=c("Maximum\n correlation", expression(paste("R"^2)), "RMSE"))+
  labs(x="Study duration")+
  theme(plot.margin = unit(c(1, 0, 0.3, 1.1), "cm"))

allderived.sbcdata 

#ok let's try to bring this together into a single figure

timeseriesplot.sbcdata
all_heat.sbcdata
allccf.sbcdata
allderived.sbcdata

complete_analysis.sbcdata<-plot_grid(timeseriesplot.sbcdata, all_heat.sbcdata, allccf.sbcdata, allderived.sbcdata, 
                                       ncol=2,  rel_heights=c(1,1,1,1), labels="AUTO", axis="b", align="v")

complete_analysis.sbcdata

pdf("figures/sbcfig1.pdf", width = 14.5, height = 9)
complete_analysis.sbcdata
dev.off()
