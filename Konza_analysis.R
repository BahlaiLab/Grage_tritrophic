# analysis for Konza site
#do analysis and create figures for each site then bring the panels for each analysis piece together in another file
#let's make our code like russian dolls


#bring in data
#for Konza we will select control grass as primary, total grasshoppers in ungrazed plots for secondary,
#and mammals in ungrazed for tertiary


konzadata.primary<-read.csv(file="cleaned_data/Konza_producer_control_grass.csv")
konzadata.secondary<-read.csv(file="cleaned_data/Konza_herbivore_ungrazed_grasshopper_total.csv")
konzadata.tertiary<-read.csv(file="cleaned_data/Konza_omnivore_ungrazed_mammal_total.csv")

#issue with column names
names(konzadata.primary)[names(konzadata.primary) == "avg.LIVEGRASS"] <- "Abundance"
names(konzadata.secondary)[names(konzadata.secondary) == "TOTAL"] <- "Abundance"
names(konzadata.tertiary)[names(konzadata.tertiary) == "TOTAL"] <- "Abundance"
names(konzadata.primary)[names(konzadata.primary) == "RECYEAR"] <- "Year"
names(konzadata.secondary)[names(konzadata.secondary) == "RECYEAR"] <- "Year"
names(konzadata.tertiary)[names(konzadata.tertiary) == "RECYEAR"] <- "Year"

#merge three tropic levels into a single frame

library(dplyr)


konzadata <- bind_rows(
  konzadata.primary %>%
    mutate(Trophic.level = "Primary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE)),
  
  konzadata.secondary %>%
    mutate(Trophic.level = "Secondary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE)),
  
  konzadata.tertiary %>%
    mutate(Trophic.level = "Tertiary",
           normAbundance = Abundance / mean(Abundance, na.rm = TRUE))
)
#create a timeseries plot and annotate it

library(ggplot2)
library(grid)


konza_summary <- konzadata %>%
  group_by(Year, Trophic.level) %>%
  summarize(meanAbundance = mean(normAbundance, na.rm = TRUE), .groups = "drop")

timeseriesplot.konzadata <- ggplot(konza_summary, aes(
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

timeseriesplot.konzadata



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

model.konzadata.primary <- multiple_breakups(konzadata.primary)
model.konzadata.secondary <- multiple_breakups(konzadata.secondary)
model.konzadata.tertiary <- multiple_breakups(konzadata.tertiary)

#and for the record, what's the stability time for each trophic level? 
stability_time(konzadata.primary) #15
stability_time(konzadata.secondary)#12
stability_time(konzadata.tertiary) #17

#and for the record, what's the proportion wih a misleading trajectory before stability for each trophic level? 
proportion_wrong_before_stability(konzadata.primary, significance=0.05) #54%
proportion_wrong_before_stability(konzadata.secondary, significance=0.05)#46%
proportion_wrong_before_stability(konzadata.tertiary, significance=0.05) #9%

linefit(standardize(konzadata.primary)) #startyear Ndata Nyears slope slopeSE, slopeP etc
linefit(standardize(konzadata.secondary))
linefit(standardize(konzadata.tertiary))

#combine data back into a single frame

model.konzadata.primary$site <- rep(c("konza"),each = 406)
model.konzadata.primary$trophic_level <- rep(c("Primary"),each = 406)

model.konzadata.secondary$site <- rep(c("konza"),each = 276)
model.konzadata.secondary$trophic_level <- rep(c("Secondary"),each = 276)

model.konzadata.tertiary$site <- rep(c("konza"),each = 496)
model.konzadata.tertiary$trophic_level <- rep(c("Tertiary"),each = 496)

#now merge all dataframes together
model.konzadata <- rbind(model.konzadata.primary, model.konzadata.secondary, model.konzadata.tertiary)

library(forcats)
library(stringr)

#create heatmap for each study duration

#use only 4, 10, 17y study durations because 17 year max stability time 
all_heat.konzadata <- model.konzadata %>%
  filter(N_years == "4"|N_years == "10"|N_years == "17") %>%
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
all_heat.konzadata

library(svglite)


#perform cross-correlation analysis
#subset data
##prim
konza_prim_17 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "17") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_16 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "16") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_15 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "15") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_14 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "14") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_13 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_12 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_11 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_10 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_9 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_8 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_7 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_6 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_5 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_prim_4 <- model.konzadata[model.konzadata$trophic_level == 'Primary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

##sec
konza_sec_17 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "17") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_16 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "16") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_15 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "15") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_14 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "14") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_13 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_12 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_11 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_10 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_9 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_8 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_7 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_6 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_5 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_sec_4 <- model.konzadata[model.konzadata$trophic_level == 'Secondary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))


##tert
konza_tert_17 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "17") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_16 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "16") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_15 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "15") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_14 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "14") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_13 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_12 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_11 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_10 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_9 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_8 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_7 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_6 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_5 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
konza_tert_4 <- model.konzadata[model.konzadata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

#Now run the cross correlation
konza17ps <- ccf(konza_prim_17$avg_slope, konza_sec_17$avg_slope, lag.max = 4)
konza16ps <- ccf(konza_prim_16$avg_slope, konza_sec_16$avg_slope, lag.max = 4)
konza15ps <- ccf(konza_prim_15$avg_slope, konza_sec_15$avg_slope, lag.max = 4)
konza14ps <- ccf(konza_prim_14$avg_slope, konza_sec_14$avg_slope, lag.max = 4)
konza13ps <- ccf(konza_prim_13$avg_slope, konza_sec_13$avg_slope, lag.max = 4)
konza12ps <- ccf(konza_prim_12$avg_slope, konza_sec_12$avg_slope, lag.max = 4)
konza11ps <- ccf(konza_prim_11$avg_slope, konza_sec_11$avg_slope, lag.max = 4)
konza10ps <- ccf(konza_prim_10$avg_slope, konza_sec_10$avg_slope, lag.max = 4)
konza9ps <- ccf(konza_prim_9$avg_slope, konza_sec_9$avg_slope, lag.max = 4)
konza8ps <- ccf(konza_prim_8$avg_slope, konza_sec_8$avg_slope, lag.max = 4)
konza7ps <- ccf(konza_prim_7$avg_slope, konza_sec_7$avg_slope, lag.max = 4)
konza6ps <- ccf(konza_prim_6$avg_slope, konza_sec_6$avg_slope, lag.max = 4)
konza5ps <- ccf(konza_prim_5$avg_slope, konza_sec_5$avg_slope, lag.max = 4)
konza4ps <- ccf(konza_prim_4$avg_slope, konza_sec_4$avg_slope, lag.max = 4)

##st
konza17st <- ccf(konza_sec_17$avg_slope, konza_tert_17$avg_slope, lag.max = 4)
konza16st <- ccf(konza_sec_16$avg_slope, konza_tert_16$avg_slope, lag.max = 4)
konza15st <- ccf(konza_sec_15$avg_slope, konza_tert_15$avg_slope, lag.max = 4)
konza14st <- ccf(konza_sec_14$avg_slope, konza_tert_14$avg_slope, lag.max = 4)
konza13st <- ccf(konza_sec_13$avg_slope, konza_tert_13$avg_slope, lag.max = 4)
konza12st <- ccf(konza_sec_12$avg_slope, konza_tert_12$avg_slope, lag.max = 4)
konza11st <- ccf(konza_sec_11$avg_slope, konza_tert_11$avg_slope, lag.max = 4)
konza10st <- ccf(konza_sec_10$avg_slope, konza_tert_10$avg_slope, lag.max = 4)
konza9st <- ccf(konza_sec_9$avg_slope, konza_tert_9$avg_slope, lag.max = 4)
konza8st <- ccf(konza_sec_8$avg_slope, konza_tert_8$avg_slope, lag.max = 4)
konza7st <- ccf(konza_sec_7$avg_slope, konza_tert_7$avg_slope, lag.max = 4)
konza6st <- ccf(konza_sec_6$avg_slope, konza_tert_6$avg_slope, lag.max = 4)
konza5st <- ccf(konza_sec_5$avg_slope, konza_tert_5$avg_slope, lag.max = 4)
konza4st <- ccf(konza_sec_4$avg_slope, konza_tert_4$avg_slope, lag.max = 4)

## tp
konza17tp <- ccf(konza_prim_17$avg_slope, konza_tert_17$avg_slope, lag.max = 4)
konza16tp <- ccf(konza_prim_16$avg_slope, konza_tert_16$avg_slope, lag.max = 4)
konza15tp <- ccf(konza_prim_15$avg_slope, konza_tert_15$avg_slope, lag.max = 4)
konza14tp <- ccf(konza_prim_14$avg_slope, konza_tert_14$avg_slope, lag.max = 4)
konza13tp <- ccf(konza_prim_13$avg_slope, konza_tert_13$avg_slope, lag.max = 4)
konza12tp <- ccf(konza_prim_12$avg_slope, konza_tert_12$avg_slope, lag.max = 4)
konza11tp <- ccf(konza_prim_11$avg_slope, konza_tert_11$avg_slope, lag.max = 4)
konza10tp <- ccf(konza_prim_10$avg_slope, konza_tert_10$avg_slope, lag.max = 4)
konza9tp <- ccf(konza_prim_9$avg_slope, konza_tert_9$avg_slope, lag.max = 4)
konza8tp <- ccf(konza_prim_8$avg_slope, konza_tert_8$avg_slope, lag.max = 4)
konza7tp <- ccf(konza_prim_7$avg_slope, konza_tert_7$avg_slope, lag.max = 4)
konza6tp <- ccf(konza_prim_6$avg_slope, konza_tert_6$avg_slope, lag.max = 4)
konza5tp <- ccf(konza_prim_5$avg_slope, konza_tert_5$avg_slope, lag.max = 4)
konza4tp <- ccf(konza_prim_4$avg_slope, konza_tert_4$avg_slope, lag.max = 4)

#prepare data for plots
konza17yrPS <- data.frame(lag=konza17ps$lag,CCF=konza17ps$acf)
konza17yrST <- data.frame(lag=konza17st$lag,CCF=konza17st$acf)
konza17yrTP <- data.frame(lag=konza17tp$lag,CCF=konza17tp$acf)

konza16yrPS <- data.frame(lag=konza16ps$lag,CCF=konza16ps$acf)
konza16yrST <- data.frame(lag=konza16st$lag,CCF=konza16st$acf)
konza16yrTP <- data.frame(lag=konza16tp$lag,CCF=konza16tp$acf)

konza15yrPS <- data.frame(lag=konza15ps$lag,CCF=konza15ps$acf)
konza15yrST <- data.frame(lag=konza15st$lag,CCF=konza15st$acf)
konza15yrTP <- data.frame(lag=konza15tp$lag,CCF=konza15tp$acf)

konza14yrPS <- data.frame(lag=konza14ps$lag,CCF=konza14ps$acf)
konza14yrST <- data.frame(lag=konza14st$lag,CCF=konza14st$acf)
konza14yrTP <- data.frame(lag=konza14tp$lag,CCF=konza14tp$acf)

konza13yrPS <- data.frame(lag=konza13ps$lag,CCF=konza13ps$acf)
konza13yrST <- data.frame(lag=konza13st$lag,CCF=konza13st$acf)
konza13yrTP <- data.frame(lag=konza13tp$lag,CCF=konza13tp$acf)

konza12yrPS <- data.frame(lag=konza12ps$lag,CCF=konza12ps$acf)
konza12yrST <- data.frame(lag=konza12st$lag,CCF=konza12st$acf)
konza12yrTP <- data.frame(lag=konza12tp$lag,CCF=konza12tp$acf)

konza11yrPS <- data.frame(lag=konza11ps$lag,CCF=konza11ps$acf)
konza11yrST <- data.frame(lag=konza11st$lag,CCF=konza11st$acf)
konza11yrTP <- data.frame(lag=konza11tp$lag,CCF=konza11tp$acf)

konza10yrPS <- data.frame(lag=konza10ps$lag,CCF=konza10ps$acf)
konza10yrST <- data.frame(lag=konza10st$lag,CCF=konza10st$acf)
konza10yrTP <- data.frame(lag=konza10tp$lag,CCF=konza10tp$acf)

konza9yrPS <- data.frame(lag=konza9ps$lag,CCF=konza9ps$acf)
konza9yrST <- data.frame(lag=konza9st$lag,CCF=konza9st$acf)
konza9yrTP <- data.frame(lag=konza9tp$lag,CCF=konza9tp$acf)

konza8yrPS <- data.frame(lag=konza8ps$lag,CCF=konza8ps$acf)
konza8yrST <- data.frame(lag=konza8st$lag,CCF=konza8st$acf)
konza8yrTP <- data.frame(lag=konza8tp$lag,CCF=konza8tp$acf)

konza7yrPS <- data.frame(lag=konza7ps$lag,CCF=konza7ps$acf)
konza7yrST <- data.frame(lag=konza7st$lag,CCF=konza7st$acf)
konza7yrTP <- data.frame(lag=konza7tp$lag,CCF=konza7tp$acf)

konza6yrPS <- data.frame(lag=konza6ps$lag,CCF=konza6ps$acf)
konza6yrST <- data.frame(lag=konza6st$lag,CCF=konza6st$acf)
konza6yrTP <- data.frame(lag=konza6tp$lag,CCF=konza6tp$acf)

konza5yrPS <- data.frame(lag=konza5ps$lag,CCF=konza5ps$acf)
konza5yrST <- data.frame(lag=konza5st$lag,CCF=konza5st$acf)
konza5yrTP <- data.frame(lag=konza5tp$lag,CCF=konza5tp$acf)

konza4yrPS <- data.frame(lag=konza4ps$lag,CCF=konza4ps$acf)
konza4yrST <- data.frame(lag=konza4st$lag,CCF=konza4st$acf)
konza4yrTP <- data.frame(lag=konza4tp$lag,CCF=konza4tp$acf)

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

### konza PS

konza17yrPSplot <- fit_and_plot_model(konza17yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza17yrPS")
#print(konza17yrPSplot$plot)
print(konza17yrPSplot$metrics)

konza16yrPSplot <- fit_and_plot_model(konza16yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza16yrPS")
#print(konza16yrPSplot$plot)
print(konza16yrPSplot$metrics)

konza15yrPSplot <- fit_and_plot_model(konza15yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza15yrPS")
#print(konza15yrPSplot$plot)
print(konza15yrPSplot$metrics)



konza14yrPSplot <- fit_and_plot_model(konza14yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza14yrPS")
#print(konza14yrPSplot$plot)
print(konza14yrPSplot$metrics)

konza13yrPSplot <- fit_and_plot_model(konza13yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza13yrPS")
#print(konza13yrPSplot$plot)
print(konza13yrPSplot$metrics)

konza12yrPSplot <- fit_and_plot_model(konza12yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza12yrPS")
#print(konza12yrPSplot$plot)
print(konza12yrPSplot$metrics)

konza11yrPSplot <- fit_and_plot_model(konza11yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza11yrPS")
#print(konza11yrPSplot$plot)
print(konza11yrPSplot$metrics)

konza10yrPSplot <- fit_and_plot_model(konza10yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza10yrPS")
#print(konza10yrPSplot$plot)
print(konza10yrPSplot$metrics)

konza9yrPSplot <- fit_and_plot_model(konza9yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza9yrPS")
#print(konza9yrPSplot$plot)
print(konza9yrPSplot$metrics)

konza8yrPSplot <- fit_and_plot_model(konza8yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza8yrPS")
#print(konza8yrPSplot$plot)
print(konza8yrPSplot$metrics)

konza7yrPSplot <- fit_and_plot_model(konza7yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza7yrPS")
#print(konza7yrPSplot$plot)
print(konza7yrPSplot$metrics)

konza6yrPSplot <- fit_and_plot_model(konza6yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza6yrPS")
#print(konza6yrPSplot$plot)
print(konza6yrPSplot$metrics)

konza5yrPSplot <- fit_and_plot_model(konza5yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza5yrPS")
#print(konza5yrPSplot$plot)
print(konza5yrPSplot$metrics)

konza4yrPSplot <- fit_and_plot_model(konza4yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza4yrPS")
#print(konza4yrPSplot$plot)
print(konza4yrPSplot$metrics)

### konza ST

konza17yrSTplot <- fit_and_plot_model(konza17yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza17yrST")
#print(konza17yrSTplot$plot)
print(konza17yrSTplot$metrics)

konza16yrSTplot <- fit_and_plot_model(konza16yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 2, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza16yrST")
#print(konza16yrSTplot$plot)
print(konza16yrSTplot$metrics)

konza15yrSTplot <- fit_and_plot_model(konza15yrST, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza15yrST")
#print(konza15yrSTplot$plot)
print(konza15yrSTplot$metrics)

konza14yrSTplot <- fit_and_plot_model(konza14yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza14yrST")
#print(konza14yrSTplot$plot)
print(konza14yrSTplot$metrics)

konza13yrSTplot <- fit_and_plot_model(konza13yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza13yrST")
#print(konza13yrSTplot$plot)
print(konza13yrSTplot$metrics)

konza12yrSTplot <- fit_and_plot_model(konza12yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 2, B = 20, C = -1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza12yrST")
#print(konza12yrSTplot$plot)
print(konza12yrSTplot$metrics)

konza11yrSTplot <- fit_and_plot_model(konza11yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 2, B = 1, C = -1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza11yrST")
#print(konza11yrSTplot$plot)
print(konza11yrSTplot$metrics)

konza10yrSTplot <- fit_and_plot_model(konza10yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza10yrST")
#print(konza10yrSTplot$plot)
print(konza10yrSTplot$metrics)

konza9yrSTplot <- fit_and_plot_model(konza9yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 2, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza9yrST")
#print(konza9yrSTplot$plot)
print(konza9yrSTplot$metrics)

konza8yrSTplot <- fit_and_plot_model(konza8yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza8yrST")
#print(konza8yrSTplot$plot)
print(konza8yrSTplot$metrics)

konza7yrSTplot <- fit_and_plot_model(konza7yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                    plot_title = "Sec/tert, 7 years",
                                    dataset_name = "konza7yrST")
#print(konza7yrSTplot$plot)
print(konza7yrSTplot$metrics)

konza6yrSTplot <- fit_and_plot_model(konza6yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza6yrST")
#print(konza6yrSTplot$plot)
print(konza6yrSTplot$metrics)

konza5yrSTplot <- fit_and_plot_model(konza5yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza5yrST")
#print(konza5yrSTplot$plot)
print(konza5yrSTplot$metrics)

konza4yrSTplot <- fit_and_plot_model(konza4yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza4yrST")
#print(konza4yrSTplot$plot)
print(konza4yrSTplot$metrics)

### konza TP

konza17yrTPplot <- fit_and_plot_model(konza17yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza17yrTP")
#print(konza17yrTPplot$plot)
print(konza17yrTPplot$metrics)

konza16yrTPplot <- fit_and_plot_model(konza16yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza16yrTP")
#print(konza16yrTPplot$plot)
print(konza16yrTPplot$metrics)

konza15yrTPplot <- fit_and_plot_model(konza15yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                      start_list = list(A = 0.5, B = 1, C = -1, D = 0),
                                      plot_title = "",
                                      dataset_name = "konza15yrTP")
#print(konza15yrTPplot$plot)
print(konza15yrTPplot$metrics)


konza14yrTPplot <- fit_and_plot_model(konza14yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.1, B = 0.6, C = 1.5, D = -0.3),
                                     plot_title = "",
                                     dataset_name = "konza14yrTP")
#print(konza14yrTPplot$plot)
print(konza14yrTPplot$metrics)

konza13yrTPplot <- fit_and_plot_model(konza13yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza13yrTP")
#print(konza13yrTPplot$plot)
print(konza13yrTPplot$metrics)

konza12yrTPplot <- fit_and_plot_model(konza12yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza12yrTP")
#print(konza12yrTPplot$plot)
print(konza12yrTPplot$metrics)

konza11yrTPplot <- fit_and_plot_model(konza11yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza11yrTP")
#print(konza11yrTPplot$plot)
print(konza11yrTPplot$metrics)

konza10yrTPplot <- fit_and_plot_model(konza10yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "konza10yrTP")
#print(konza10yrTPplot$plot)
print(konza10yrTPplot$metrics)

konza9yrTPplot <- fit_and_plot_model(konza9yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza9yrTP")
#print(konza9yrTPplot$plot)
print(konza9yrTPplot$metrics)


konza8yrTPplot <- fit_and_plot_model(konza8yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza8yrTP")
#print(konza8yrTPplot$plot)
print(konza8yrTPplot$metrics)

konza7yrTPplot <- fit_and_plot_model(konza7yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza7yrTP")
#print(konza7yrTPplot$plot)
print(konza7yrTPplot$metrics)

konza6yrTPplot <- fit_and_plot_model(konza6yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza6yrTP")
#print(konza6yrTPplot$plot)
print(konza6yrTPplot$metrics)

konza5yrTPplot <- fit_and_plot_model(konza5yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza5yrTP")
#print(konza5yrTPplot$plot)
print(konza5yrTPplot$metrics)

konza4yrTPplot <- fit_and_plot_model(konza4yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "konza4yrTP")
#print(konza4yrTPplot$plot)
print(konza4yrTPplot$metrics)

library(ggpubr)
library(cowplot)
# ggarrange(konza8yrPSplot$plot,konza8yrSTplot$plot,konza8yrTPplot$plot, 
#           #konza7yrPSplot$plot,konza7yrSTplot$plot,konza7yrTPplot$plot,
#           konza6yrPSplot$plot,konza6yrSTplot$plot,konza6yrTPplot$plot,
#           #konza5yrPSplot$plot,konza5yrSTplot$plot,konza5yrTPplot$plot,
#           konza4yrPSplot$plot,konza4yrSTplot$plot,konza4yrTPplot$plot,
#           ncol = 3, nrow=3)

prim.sec<-text_grob("Primary\n and\n Secondary", size = 9)
sec.tert<-text_grob("Secondary\n and\n Tertiary", size = 9)
prim.tert<-text_grob("Primary\n and\n Tertiary", size = 9)

pssine<-plot_grid(konza4yrPSplot$plot, konza10yrPSplot$plot, konza17yrPSplot$plot, prim.sec, 
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
stsine<-plot_grid(konza4yrSTplot$plot, konza10yrSTplot$plot, konza17yrSTplot$plot, sec.tert,
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
ptsine<-plot_grid( konza4yrTPplot$plot, konza10yrTPplot$plot, konza17yrTPplot$plot, prim.tert,
                   ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))

studyduration<-text_grob("Study duration", size=11)
duration4<-text_grob("   4", size=10)
durationmid<-text_grob("   10", size=10)
durationtop<-text_grob("   17", size=10)
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


allccf.konzadata<-plot_grid(durationlabel, durationbar, plots_together_lab,   ncol=1, rel_heights=c(0.1, 0.08, 1.9))

allccf.konzadata

#compile fit metrics and make plots
konzadf <- bind_rows(konza17yrPSplot$metrics,konza17yrSTplot$metrics,konza17yrTPplot$metrics,
                    konza16yrPSplot$metrics,konza16yrSTplot$metrics,konza16yrTPplot$metrics,
                    konza15yrPSplot$metrics,konza15yrSTplot$metrics,konza15yrTPplot$metrics,
                    konza14yrPSplot$metrics,konza14yrSTplot$metrics,konza14yrTPplot$metrics,
                    konza13yrPSplot$metrics,konza13yrSTplot$metrics,konza13yrTPplot$metrics,
                    konza12yrPSplot$metrics,konza12yrSTplot$metrics,konza12yrTPplot$metrics,
                    konza11yrPSplot$metrics,konza11yrSTplot$metrics,konza11yrTPplot$metrics,
                    konza10yrPSplot$metrics,konza10yrSTplot$metrics,konza10yrTPplot$metrics,
                    konza9yrPSplot$metrics,konza9yrSTplot$metrics,konza9yrTPplot$metrics,
                    konza8yrPSplot$metrics,konza8yrSTplot$metrics,konza8yrTPplot$metrics, 
                    konza7yrPSplot$metrics,konza7yrSTplot$metrics,konza7yrTPplot$metrics,
                    konza6yrPSplot$metrics,konza6yrSTplot$metrics,konza6yrTPplot$metrics,
                    konza5yrPSplot$metrics,konza5yrSTplot$metrics,konza5yrTPplot$metrics,
                    konza4yrPSplot$metrics,konza4yrSTplot$metrics,konza4yrTPplot$metrics)
konzadf


library(ggpmisc)
library(tidyr)
konzadf$Year <- as.numeric(konzadf$Year)
allderived.konzadata <- konzadf %>%
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

allderived.konzadata 

#ok let's try to bring this together into a single figure

timeseriesplot.konzadata
all_heat.konzadata
allccf.konzadata
allderived.konzadata

complete_analysis.konzadata<-plot_grid(timeseriesplot.konzadata, all_heat.konzadata, allccf.konzadata, allderived.konzadata, 
                                      ncol=2,  rel_heights=c(1,1,1,1), labels="AUTO", axis="b", align="v")

complete_analysis.konzadata

pdf("figures/konzafig1.pdf", width = 14.5, height = 9)
complete_analysis.konzadata
dev.off()
