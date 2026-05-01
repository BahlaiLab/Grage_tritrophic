# analysis for hubbard site
#do analysis and create figures for each site then bring the panels for each analysis piece together in another file
#let's make our code like russian dolls
#doing hubbard second


#bring in data
#for hubbard we will select litter mass as primary, caterpillar biomass in maple for secondary,
#and total birds for tertiary


hubbarddata.primary<-read.csv(file="cleaned_data/hubbard_producer_litter_mass.csv")
hubbarddata.secondary<-read.csv(file="cleaned_data/hubbard_herbivore_maple_biomass.csv")
hubbarddata.tertiary<-read.csv(file="cleaned_data/hubbard_omnivore_bird_total.csv")

#issue with column names
names(hubbarddata.primary)[names(hubbarddata.primary) == "Avlitter"] <- "Abundance"
names(hubbarddata.secondary)[names(hubbarddata.secondary) == "biomass"] <- "Abundance"
names(hubbarddata.tertiary)[names(hubbarddata.tertiary) == "total"] <- "Abundance"
names(hubbarddata.primary)[names(hubbarddata.primary) == "YEAR"] <- "Year"
names(hubbarddata.secondary)[names(hubbarddata.secondary) == "Year"] <- "Year"
names(hubbarddata.tertiary)[names(hubbarddata.tertiary) == "Year"] <- "Year"

#merge three tropic levels into a single frame

library(dplyr)


hubbarddata <- bind_rows(
  hubbarddata.primary %>%
    mutate(Trophic.level = "Primary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE)),
  
  hubbarddata.secondary %>%
    mutate(Trophic.level = "Secondary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE)),
  
  hubbarddata.tertiary %>%
    mutate(Trophic.level = "Tertiary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE))
)
#create a timeseries plot and annotate it

library(ggplot2)
library(grid)


hubbard_summary <- hubbarddata %>%
  group_by(Year, Trophic.level) %>%
  summarize(meanAbundance = mean(normAbundance, na.rm = TRUE), .groups = "drop")

timeseriesplot.hubbarddata <- ggplot(hubbard_summary, aes(
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

timeseriesplot.hubbarddata




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

model.hubbarddata.primary <- multiple_breakups(hubbarddata.primary)
model.hubbarddata.secondary <- multiple_breakups(hubbarddata.secondary)
model.hubbarddata.tertiary <- multiple_breakups(hubbarddata.tertiary)

#and for the record, what's the stability time for each trophic level? 
stability_time(hubbarddata.primary) #11
stability_time(hubbarddata.secondary)#13
stability_time(hubbarddata.tertiary) #33

#and for the record, what's the proportion wih a misleading trajectory before stability for each trophic level? 
proportion_wrong_before_stability(hubbarddata.primary, significance=0.05) #77%
proportion_wrong_before_stability(hubbarddata.secondary, significance=0.05)#64%
proportion_wrong_before_stability(hubbarddata.tertiary, significance=0.05) #58%

linefit(standardize(hubbarddata.primary)) #startyear Ndata Nyears slope slopeSE, slopeP etc
linefit(standardize(hubbarddata.secondary))
linefit(standardize(hubbarddata.tertiary))

#combine data back into a single frame

model.hubbarddata.primary$site <- rep(c("hubbard"),each = 351)
model.hubbarddata.primary$trophic_level <- rep(c("Primary"),each = 351)

model.hubbarddata.secondary$site <- rep(c("hubbard"),each = 496)
model.hubbarddata.secondary$trophic_level <- rep(c("Secondary"),each = 496)

model.hubbarddata.tertiary$site <- rep(c("hubbard"),each = 1035)
model.hubbarddata.tertiary$trophic_level <- rep(c("Tertiary"),each = 1035)

#now merge all dataframes together
model.hubbarddata <- rbind(model.hubbarddata.primary, model.hubbarddata.secondary, model.hubbarddata.tertiary)

library(forcats)
library(stringr)

#create heatmap for each study duration
#because one stability time>> overlap of time series, use second highest stability time: 13

#use only 4, 8, 13y study duration 
all_heat.hubbarddata <- model.hubbarddata %>%
  filter(N_years == "4"|N_years == "8"|N_years == "13") %>%
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
all_heat.hubbarddata

library(svglite)


#perform cross-correlation analysis
#subset data
##prim

hubbard_prim_13 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_12 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_11 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_10 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_9 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_8 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_7 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_6 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_5 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_prim_4 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Primary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

##sec

hubbard_sec_13 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_12 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_11 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_10 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_9 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_8 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_7 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_6 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_5 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_sec_4 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Secondary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))


##tert

hubbard_tert_13 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_12 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_11 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_10 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_9 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_8 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_7 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_6 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_5 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
hubbard_tert_4 <- model.hubbarddata[model.hubbarddata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

#Now run the cross correlation

hubbard13ps <- ccf(hubbard_prim_13$avg_slope, hubbard_sec_13$avg_slope, lag.max = 4)
hubbard12ps <- ccf(hubbard_prim_12$avg_slope, hubbard_sec_12$avg_slope, lag.max = 4)
hubbard11ps <- ccf(hubbard_prim_11$avg_slope, hubbard_sec_11$avg_slope, lag.max = 4)
hubbard10ps <- ccf(hubbard_prim_10$avg_slope, hubbard_sec_10$avg_slope, lag.max = 4)
hubbard9ps <- ccf(hubbard_prim_9$avg_slope, hubbard_sec_9$avg_slope, lag.max = 4)
hubbard8ps <- ccf(hubbard_prim_8$avg_slope, hubbard_sec_8$avg_slope, lag.max = 4)
hubbard7ps <- ccf(hubbard_prim_7$avg_slope, hubbard_sec_7$avg_slope, lag.max = 4)
hubbard6ps <- ccf(hubbard_prim_6$avg_slope, hubbard_sec_6$avg_slope, lag.max = 4)
hubbard5ps <- ccf(hubbard_prim_5$avg_slope, hubbard_sec_5$avg_slope, lag.max = 4)
hubbard4ps <- ccf(hubbard_prim_4$avg_slope, hubbard_sec_4$avg_slope, lag.max = 4)

##st

hubbard13st <- ccf(hubbard_sec_13$avg_slope, hubbard_tert_13$avg_slope, lag.max = 4)
hubbard12st <- ccf(hubbard_sec_12$avg_slope, hubbard_tert_12$avg_slope, lag.max = 4)
hubbard11st <- ccf(hubbard_sec_11$avg_slope, hubbard_tert_11$avg_slope, lag.max = 4)
hubbard10st <- ccf(hubbard_sec_10$avg_slope, hubbard_tert_10$avg_slope, lag.max = 4)
hubbard9st <- ccf(hubbard_sec_9$avg_slope, hubbard_tert_9$avg_slope, lag.max = 4)
hubbard8st <- ccf(hubbard_sec_8$avg_slope, hubbard_tert_8$avg_slope, lag.max = 4)
hubbard7st <- ccf(hubbard_sec_7$avg_slope, hubbard_tert_7$avg_slope, lag.max = 4)
hubbard6st <- ccf(hubbard_sec_6$avg_slope, hubbard_tert_6$avg_slope, lag.max = 4)
hubbard5st <- ccf(hubbard_sec_5$avg_slope, hubbard_tert_5$avg_slope, lag.max = 4)
hubbard4st <- ccf(hubbard_sec_4$avg_slope, hubbard_tert_4$avg_slope, lag.max = 4)

## tp

hubbard13tp <- ccf(hubbard_prim_13$avg_slope, hubbard_tert_13$avg_slope, lag.max = 4)
hubbard12tp <- ccf(hubbard_prim_12$avg_slope, hubbard_tert_12$avg_slope, lag.max = 4)
hubbard11tp <- ccf(hubbard_prim_11$avg_slope, hubbard_tert_11$avg_slope, lag.max = 4)
hubbard10tp <- ccf(hubbard_prim_10$avg_slope, hubbard_tert_10$avg_slope, lag.max = 4)
hubbard9tp <- ccf(hubbard_prim_9$avg_slope, hubbard_tert_9$avg_slope, lag.max = 4)
hubbard8tp <- ccf(hubbard_prim_8$avg_slope, hubbard_tert_8$avg_slope, lag.max = 4)
hubbard7tp <- ccf(hubbard_prim_7$avg_slope, hubbard_tert_7$avg_slope, lag.max = 4)
hubbard6tp <- ccf(hubbard_prim_6$avg_slope, hubbard_tert_6$avg_slope, lag.max = 4)
hubbard5tp <- ccf(hubbard_prim_5$avg_slope, hubbard_tert_5$avg_slope, lag.max = 4)
hubbard4tp <- ccf(hubbard_prim_4$avg_slope, hubbard_tert_4$avg_slope, lag.max = 4)

#prepare data for plots

hubbard13yrPS <- data.frame(lag=hubbard13ps$lag,CCF=hubbard13ps$acf)
hubbard13yrST <- data.frame(lag=hubbard13st$lag,CCF=hubbard13st$acf)
hubbard13yrTP <- data.frame(lag=hubbard13tp$lag,CCF=hubbard13tp$acf)

hubbard12yrPS <- data.frame(lag=hubbard12ps$lag,CCF=hubbard12ps$acf)
hubbard12yrST <- data.frame(lag=hubbard12st$lag,CCF=hubbard12st$acf)
hubbard12yrTP <- data.frame(lag=hubbard12tp$lag,CCF=hubbard12tp$acf)

hubbard11yrPS <- data.frame(lag=hubbard11ps$lag,CCF=hubbard11ps$acf)
hubbard11yrST <- data.frame(lag=hubbard11st$lag,CCF=hubbard11st$acf)
hubbard11yrTP <- data.frame(lag=hubbard11tp$lag,CCF=hubbard11tp$acf)

hubbard10yrPS <- data.frame(lag=hubbard10ps$lag,CCF=hubbard10ps$acf)
hubbard10yrST <- data.frame(lag=hubbard10st$lag,CCF=hubbard10st$acf)
hubbard10yrTP <- data.frame(lag=hubbard10tp$lag,CCF=hubbard10tp$acf)

hubbard9yrPS <- data.frame(lag=hubbard9ps$lag,CCF=hubbard9ps$acf)
hubbard9yrST <- data.frame(lag=hubbard9st$lag,CCF=hubbard9st$acf)
hubbard9yrTP <- data.frame(lag=hubbard9tp$lag,CCF=hubbard9tp$acf)

hubbard8yrPS <- data.frame(lag=hubbard8ps$lag,CCF=hubbard8ps$acf)
hubbard8yrST <- data.frame(lag=hubbard8st$lag,CCF=hubbard8st$acf)
hubbard8yrTP <- data.frame(lag=hubbard8tp$lag,CCF=hubbard8tp$acf)

hubbard7yrPS <- data.frame(lag=hubbard7ps$lag,CCF=hubbard7ps$acf)
hubbard7yrST <- data.frame(lag=hubbard7st$lag,CCF=hubbard7st$acf)
hubbard7yrTP <- data.frame(lag=hubbard7tp$lag,CCF=hubbard7tp$acf)

hubbard6yrPS <- data.frame(lag=hubbard6ps$lag,CCF=hubbard6ps$acf)
hubbard6yrST <- data.frame(lag=hubbard6st$lag,CCF=hubbard6st$acf)
hubbard6yrTP <- data.frame(lag=hubbard6tp$lag,CCF=hubbard6tp$acf)

hubbard5yrPS <- data.frame(lag=hubbard5ps$lag,CCF=hubbard5ps$acf)
hubbard5yrST <- data.frame(lag=hubbard5st$lag,CCF=hubbard5st$acf)
hubbard5yrTP <- data.frame(lag=hubbard5tp$lag,CCF=hubbard5tp$acf)

hubbard4yrPS <- data.frame(lag=hubbard4ps$lag,CCF=hubbard4ps$acf)
hubbard4yrST <- data.frame(lag=hubbard4st$lag,CCF=hubbard4st$acf)
hubbard4yrTP <- data.frame(lag=hubbard4tp$lag,CCF=hubbard4tp$acf)

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
    ylim(-1,1)+
    theme_bw()
  list(plot = plot, metrics = metrics_df)
}

### hubbard PS

hubbard13yrPSplot <- fit_and_plot_model(hubbard13yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard13yrPS")
#print(hubbard13yrPSplot$plot)
print(hubbard13yrPSplot$metrics)

hubbard12yrPSplot <- fit_and_plot_model(hubbard12yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard12yrPS")
#print(hubbard12yrPSplot$plot)
print(hubbard12yrPSplot$metrics)

hubbard11yrPSplot <- fit_and_plot_model(hubbard11yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard11yrPS")
#print(hubbard11yrPSplot$plot)
print(hubbard11yrPSplot$metrics)

hubbard10yrPSplot <- fit_and_plot_model(hubbard10yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard10yrPS")
#print(hubbard10yrPSplot$plot)
print(hubbard10yrPSplot$metrics)

hubbard9yrPSplot <- fit_and_plot_model(hubbard9yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 1, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard9yrPS")
#print(hubbard9yrPSplot$plot)
print(hubbard9yrPSplot$metrics)

hubbard8yrPSplot <- fit_and_plot_model(hubbard8yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 1, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard8yrPS")
#print(hubbard8yrPSplot$plot)
print(hubbard8yrPSplot$metrics)

hubbard7yrPSplot <- fit_and_plot_model(hubbard7yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 1, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard7yrPS")
#print(hubbard7yrPSplot$plot)
print(hubbard7yrPSplot$metrics)

hubbard6yrPSplot <- fit_and_plot_model(hubbard6yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard6yrPS")
#print(hubbard6yrPSplot$plot)
print(hubbard6yrPSplot$metrics)

hubbard5yrPSplot <- fit_and_plot_model(hubbard5yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard5yrPS")
#print(hubbard5yrPSplot$plot)
print(hubbard5yrPSplot$metrics)

hubbard4yrPSplot <- fit_and_plot_model(hubbard4yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard4yrPS")
#print(hubbard4yrPSplot$plot)
print(hubbard4yrPSplot$metrics)

### hubbard ST



hubbard13yrSTplot <- fit_and_plot_model(hubbard13yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard13yrST")
#print(hubbard13yrSTplot$plot)
print(hubbard13yrSTplot$metrics)

hubbard12yrSTplot <- fit_and_plot_model(hubbard12yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 2, B = 20, C = -1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard12yrST")
#print(hubbard12yrSTplot$plot)
print(hubbard12yrSTplot$metrics)

hubbard11yrSTplot <- fit_and_plot_model(hubbard11yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 2, B = 1, C = -1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard11yrST")
#print(hubbard11yrSTplot$plot)
print(hubbard11yrSTplot$metrics)

hubbard10yrSTplot <- fit_and_plot_model(hubbard10yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard10yrST")
#print(hubbard10yrSTplot$plot)
print(hubbard10yrSTplot$metrics)

hubbard9yrSTplot <- fit_and_plot_model(hubbard9yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 2, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard9yrST")
#print(hubbard9yrSTplot$plot)
print(hubbard9yrSTplot$metrics)

hubbard8yrSTplot <- fit_and_plot_model(hubbard8yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard8yrST")
#print(hubbard8yrSTplot$plot)
print(hubbard8yrSTplot$metrics)

hubbard7yrSTplot <- fit_and_plot_model(hubbard7yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                     plot_title = "Sec/tert, 7 years",
                                     dataset_name = "hubbard7yrST")
#print(hubbard7yrSTplot$plot)
print(hubbard7yrSTplot$metrics)

hubbard6yrSTplot <- fit_and_plot_model(hubbard6yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 1, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard6yrST")
#print(hubbard6yrSTplot$plot)
print(hubbard6yrSTplot$metrics)

hubbard5yrSTplot <- fit_and_plot_model(hubbard5yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard5yrST")
#print(hubbard5yrSTplot$plot)
print(hubbard5yrSTplot$metrics)

hubbard4yrSTplot <- fit_and_plot_model(hubbard4yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard4yrST")
#print(hubbard4yrSTplot$plot)
print(hubbard4yrSTplot$metrics)

### hubbard TP

hubbard13yrTPplot <- fit_and_plot_model(hubbard13yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard13yrTP")
#print(hubbard13yrTPplot$plot)
print(hubbard13yrTPplot$metrics)

hubbard12yrTPplot <- fit_and_plot_model(hubbard12yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard12yrTP")
#print(hubbard12yrTPplot$plot)
print(hubbard12yrTPplot$metrics)

hubbard11yrTPplot <- fit_and_plot_model(hubbard11yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard11yrTP")
#print(hubbard11yrTPplot$plot)
print(hubbard11yrTPplot$metrics)

hubbard10yrTPplot <- fit_and_plot_model(hubbard10yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "hubbard10yrTP")
#print(hubbard10yrTPplot$plot)
print(hubbard10yrTPplot$metrics)

hubbard9yrTPplot <- fit_and_plot_model(hubbard9yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard9yrTP")
#print(hubbard9yrTPplot$plot)
print(hubbard9yrTPplot$metrics)


hubbard8yrTPplot <- fit_and_plot_model(hubbard8yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard8yrTP")
#print(hubbard8yrTPplot$plot)
print(hubbard8yrTPplot$metrics)

hubbard7yrTPplot <- fit_and_plot_model(hubbard7yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard7yrTP")
#print(hubbard7yrTPplot$plot)
print(hubbard7yrTPplot$metrics)

hubbard6yrTPplot <- fit_and_plot_model(hubbard6yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard6yrTP")
#print(hubbard6yrTPplot$plot)
print(hubbard6yrTPplot$metrics)

hubbard5yrTPplot <- fit_and_plot_model(hubbard5yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard5yrTP")
#print(hubbard5yrTPplot$plot)
print(hubbard5yrTPplot$metrics)

hubbard4yrTPplot <- fit_and_plot_model(hubbard4yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "hubbard4yrTP")
#print(hubbard4yrTPplot$plot)
print(hubbard4yrTPplot$metrics)

library(ggpubr)
library(cowplot)
# ggarrange(hubbard8yrPSplot$plot,hubbard8yrSTplot$plot,hubbard8yrTPplot$plot, 
#           #hubbard7yrPSplot$plot,hubbard7yrSTplot$plot,hubbard7yrTPplot$plot,
#           hubbard6yrPSplot$plot,hubbard6yrSTplot$plot,hubbard6yrTPplot$plot,
#           #hubbard5yrPSplot$plot,hubbard5yrSTplot$plot,hubbard5yrTPplot$plot,
#           hubbard4yrPSplot$plot,hubbard4yrSTplot$plot,hubbard4yrTPplot$plot,
#           ncol = 3, nrow=3)

prim.sec<-text_grob("Primary\n and\n Secondary", size = 9)
sec.tert<-text_grob("Secondary\n and\n Tertiary", size = 9)
prim.tert<-text_grob("Primary\n and\n Tertiary", size = 9)

pssine<-plot_grid(hubbard4yrPSplot$plot, hubbard8yrPSplot$plot, hubbard13yrPSplot$plot, prim.sec, 
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
stsine<-plot_grid(hubbard4yrSTplot$plot, hubbard8yrSTplot$plot, hubbard13yrSTplot$plot, sec.tert,
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
ptsine<-plot_grid( hubbard4yrTPplot$plot, hubbard8yrTPplot$plot, hubbard13yrTPplot$plot, prim.tert,
                   ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))

studyduration<-text_grob("Study duration", size=11)
duration4<-text_grob("   4", size=10)
durationmid<-text_grob("   8", size=10)
durationtop<-text_grob("   13", size=10)
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


allccf.hubbarddata<-plot_grid(durationlabel, durationbar, plots_together_lab,   ncol=1, rel_heights=c(0.1, 0.08, 1.9))

allccf.hubbarddata

#compile fit metrics and make plots
hubbarddf <- bind_rows(hubbard13yrPSplot$metrics,hubbard13yrSTplot$metrics,hubbard13yrTPplot$metrics,
                     hubbard12yrPSplot$metrics,hubbard12yrSTplot$metrics,hubbard12yrTPplot$metrics,
                     hubbard11yrPSplot$metrics,hubbard11yrSTplot$metrics,hubbard11yrTPplot$metrics,
                     hubbard10yrPSplot$metrics,hubbard10yrSTplot$metrics,hubbard10yrTPplot$metrics,
                     hubbard9yrPSplot$metrics,hubbard9yrSTplot$metrics,hubbard9yrTPplot$metrics,
                     hubbard8yrPSplot$metrics,hubbard8yrSTplot$metrics,hubbard8yrTPplot$metrics, 
                     hubbard7yrPSplot$metrics,hubbard7yrSTplot$metrics,hubbard7yrTPplot$metrics,
                     hubbard6yrPSplot$metrics,hubbard6yrSTplot$metrics,hubbard6yrTPplot$metrics,
                     hubbard5yrPSplot$metrics,hubbard5yrSTplot$metrics,hubbard5yrTPplot$metrics,
                     hubbard4yrPSplot$metrics,hubbard4yrSTplot$metrics,hubbard4yrTPplot$metrics)
hubbarddf


library(ggpmisc)
library(tidyr)
hubbarddf$Year <- as.numeric(hubbarddf$Year)
allderived.hubbarddata <- hubbarddf %>%
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

allderived.hubbarddata 

#ok let's try to bring this together into a single figure

timeseriesplot.hubbarddata
all_heat.hubbarddata
allccf.hubbarddata
allderived.hubbarddata

complete_analysis.hubbarddata<-plot_grid(timeseriesplot.hubbarddata, all_heat.hubbarddata, allccf.hubbarddata, allderived.hubbarddata, 
                                       ncol=2,  rel_heights=c(1,1,1,1), labels="AUTO", axis="b", align="v")

complete_analysis.hubbarddata

pdf("figures/hubbardfig1.pdf", width = 14.5, height = 9)
complete_analysis.hubbarddata
dev.off()
