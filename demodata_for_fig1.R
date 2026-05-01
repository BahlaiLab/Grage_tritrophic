#demo data and analysis for figure 1
#data is a 30 year hypothetical timeseries with three trophic levels
#primary is year zero, secondary has a year lag of 1 and is transformed by dividing by 1.5
#tertiary is lagged from secondary by two years and is transformed by dividing by 2

#bring in data
demodata<-read.csv(file="cleaned_data/demodata.csv")

#prep it for analysis
demodata1<-demodata
demodata1$Trophic.level<-NULL

demodata.primary<-demodata1[which(demodata$Trophic.level=="Primary"),]
demodata.secondary<-demodata1[which(demodata$Trophic.level=="Secondary"),]
demodata.tertiary<-demodata1[which(demodata$Trophic.level=="Tertiary"),]

#create a timeseries plot and annotate it

library(ggplot2)
library(grid)
library(dplyr)

timeseriesplot.demodata<-ggplot(demodata, aes(x=Year, y=Abundance, fill=Trophic.level, 
                                     color=Trophic.level, shape=Trophic.level, linetype=Trophic.level)) +
  geom_line(linewidth=0.6)+
  geom_point(size=3)+
  scale_color_manual(values=c("tan4", "orange3","firebrick4"), name="Trophic level")+
  scale_fill_manual(values=c("tan","orange","firebrick1"), name="Trophic level")+
  scale_shape_manual(values=c(23,22,21), name="Trophic level")+
  scale_linetype_manual(values=c("solid","solid", "solid"), name="Trophic level")+
  theme_classic()+theme(legend.position = "inside",legend.position.inside=c(0.1, 0.85))+
  annotate('segment', x = 7, y = 98, xend = 7, yend = 160, color="black", lty="dashed") +
  annotate('segment', x = 8, y = 66, xend = 8, yend = 128, color="black", lty="dashed") +
  annotate('segment', x = 10, y = 43, xend = 10, yend = 160, color="black", lty="dashed") +
  annotate('text', x = 7.9, y = 135, 
           label =  'italic(L[P~"&"~S])',
           size = 3, 
           color = "black",
           parse=T)  +
  annotate('segment', x = 7, xend = 8, y = 128, yend = 128,
    arrow = arrow(ends = "both", length = unit(0.05, "inches")))+
  annotate('text', x = 9, y = 121, 
           label =  'italic(L[S~"&"~T])',
           size = 3, 
           color = "black",
           parse=T)  +
  annotate('segment', x = 8, xend = 10, y = 116, yend = 116,
           arrow = arrow(ends = "both", length = unit(0.05, "inches")))+
  annotate('text', x = 8.5, y = 150, 
           label =  'italic(L[P~"&"~T])',
           size = 3, 
           color = "black",
           parse=T)  +
  annotate('segment', x = 7, xend = 10, y = 145, yend = 145,
           arrow = arrow(ends = "both", length = unit(0.05, "inches")))




timeseriesplot.demodata




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

model.demodata.primary <- multiple_breakups(demodata.primary)
model.demodata.secondary <- multiple_breakups(demodata.secondary)
model.demodata.tertiary <- multiple_breakups(demodata.tertiary)

#and for the record, what's the stability time for each trophic level? 
stability_time(demodata.primary) #14
stability_time(demodata.secondary)#14
stability_time(demodata.tertiary) #11

#and for the record, what's the proportion wih a misleading trajectory before stability for each trophic level? 
proportion_wrong_before_stability(demodata.primary, significance=0.05) #99%
proportion_wrong_before_stability(demodata.secondary, significance=0.05)#98%
proportion_wrong_before_stability(demodata.tertiary, significance=0.05) #97%

linefit(standardize(demodata.primary)) #startyear Ndata Nyears slope slopeSE, slopeP etc
linefit(standardize(demodata.secondary))
linefit(standardize(demodata.tertiary))

#combine data back into a single frame

model.demodata.primary$site <- rep(c("demo"),each = 406)
model.demodata.primary$trophic_level <- rep(c("Primary"),each = 406)

model.demodata.secondary$site <- rep(c("demo"),each = 406)
model.demodata.secondary$trophic_level <- rep(c("Secondary"),each = 406)

model.demodata.tertiary$site <- rep(c("demo"),each = 406)
model.demodata.tertiary$trophic_level <- rep(c("Tertiary"),each = 406)

#now merge all dataframes together
model.demodata <- rbind(model.demodata.primary, model.demodata.secondary, model.demodata.tertiary)

library(forcats)
library(stringr)

#create heatmap for each study duration

#use only 4, 9, 14y study durations
all_heat.demodata <- model.demodata %>%
  filter(N_years == "4"|N_years == "9"|N_years == "14") %>%
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
all_heat.demodata

library(svglite)


#perform cross-correlation analysis
#subset data
##prim
demo_prim_14 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "14") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_13 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_12 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_11 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_10 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_9 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_8 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_7 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_6 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_5 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_prim_4 <- model.demodata[model.demodata$trophic_level == 'Primary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

##sec
demo_sec_14 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "14") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_13 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_12 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_11 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_10 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_9 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_8 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_7 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_6 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_5 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_sec_4 <- model.demodata[model.demodata$trophic_level == 'Secondary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))


##tert
demo_tert_14 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "14") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_13 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "13") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_12 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "12") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_11 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "11") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_10 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "10") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_9 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "9") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_8 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "8") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_7 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "7") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_6 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "6") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_5 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "5") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))
demo_tert_4 <- model.demodata[model.demodata$trophic_level == 'Tertiary',] %>%
  filter(N_years == "4") %>%
  group_by(start_year) %>%
  summarise(avg_slope = mean(slope))

#Now run the cross correlation
demo14ps <- ccf(demo_prim_14$avg_slope, demo_sec_14$avg_slope, lag.max = 4)
demo13ps <- ccf(demo_prim_13$avg_slope, demo_sec_13$avg_slope, lag.max = 4)
demo12ps <- ccf(demo_prim_12$avg_slope, demo_sec_12$avg_slope, lag.max = 4)
demo11ps <- ccf(demo_prim_11$avg_slope, demo_sec_11$avg_slope, lag.max = 4)
demo10ps <- ccf(demo_prim_10$avg_slope, demo_sec_10$avg_slope, lag.max = 4)
demo9ps <- ccf(demo_prim_9$avg_slope, demo_sec_9$avg_slope, lag.max = 4)
demo8ps <- ccf(demo_prim_8$avg_slope, demo_sec_8$avg_slope, lag.max = 4)
demo7ps <- ccf(demo_prim_7$avg_slope, demo_sec_7$avg_slope, lag.max = 4)
demo6ps <- ccf(demo_prim_6$avg_slope, demo_sec_6$avg_slope, lag.max = 4)
demo5ps <- ccf(demo_prim_5$avg_slope, demo_sec_5$avg_slope, lag.max = 4)
demo4ps <- ccf(demo_prim_4$avg_slope, demo_sec_4$avg_slope, lag.max = 4)

##st
demo14st <- ccf(demo_sec_14$avg_slope, demo_tert_14$avg_slope, lag.max = 4)
demo13st <- ccf(demo_sec_13$avg_slope, demo_tert_13$avg_slope, lag.max = 4)
demo12st <- ccf(demo_sec_12$avg_slope, demo_tert_12$avg_slope, lag.max = 4)
demo11st <- ccf(demo_sec_11$avg_slope, demo_tert_11$avg_slope, lag.max = 4)
demo10st <- ccf(demo_sec_10$avg_slope, demo_tert_10$avg_slope, lag.max = 4)
demo9st <- ccf(demo_sec_9$avg_slope, demo_tert_9$avg_slope, lag.max = 4)
demo8st <- ccf(demo_sec_8$avg_slope, demo_tert_8$avg_slope, lag.max = 4)
demo7st <- ccf(demo_sec_7$avg_slope, demo_tert_7$avg_slope, lag.max = 4)
demo6st <- ccf(demo_sec_6$avg_slope, demo_tert_6$avg_slope, lag.max = 4)
demo5st <- ccf(demo_sec_5$avg_slope, demo_tert_5$avg_slope, lag.max = 4)
demo4st <- ccf(demo_sec_4$avg_slope, demo_tert_4$avg_slope, lag.max = 4)

## tp
demo14tp <- ccf(demo_prim_14$avg_slope, demo_tert_14$avg_slope, lag.max = 4)
demo13tp <- ccf(demo_prim_13$avg_slope, demo_tert_13$avg_slope, lag.max = 4)
demo12tp <- ccf(demo_prim_12$avg_slope, demo_tert_12$avg_slope, lag.max = 4)
demo11tp <- ccf(demo_prim_11$avg_slope, demo_tert_11$avg_slope, lag.max = 4)
demo10tp <- ccf(demo_prim_10$avg_slope, demo_tert_10$avg_slope, lag.max = 4)
demo9tp <- ccf(demo_prim_9$avg_slope, demo_tert_9$avg_slope, lag.max = 4)
demo8tp <- ccf(demo_prim_8$avg_slope, demo_tert_8$avg_slope, lag.max = 4)
demo7tp <- ccf(demo_prim_7$avg_slope, demo_tert_7$avg_slope, lag.max = 4)
demo6tp <- ccf(demo_prim_6$avg_slope, demo_tert_6$avg_slope, lag.max = 4)
demo5tp <- ccf(demo_prim_5$avg_slope, demo_tert_5$avg_slope, lag.max = 4)
demo4tp <- ccf(demo_prim_4$avg_slope, demo_tert_4$avg_slope, lag.max = 4)

#prepare data for plots
demo14yrPS <- data.frame(lag=demo14ps$lag,CCF=demo14ps$acf)
demo14yrST <- data.frame(lag=demo14st$lag,CCF=demo14st$acf)
demo14yrTP <- data.frame(lag=demo14tp$lag,CCF=demo14tp$acf)

demo13yrPS <- data.frame(lag=demo13ps$lag,CCF=demo13ps$acf)
demo13yrST <- data.frame(lag=demo13st$lag,CCF=demo13st$acf)
demo13yrTP <- data.frame(lag=demo13tp$lag,CCF=demo13tp$acf)

demo12yrPS <- data.frame(lag=demo12ps$lag,CCF=demo12ps$acf)
demo12yrST <- data.frame(lag=demo12st$lag,CCF=demo12st$acf)
demo12yrTP <- data.frame(lag=demo12tp$lag,CCF=demo12tp$acf)

demo11yrPS <- data.frame(lag=demo11ps$lag,CCF=demo11ps$acf)
demo11yrST <- data.frame(lag=demo11st$lag,CCF=demo11st$acf)
demo11yrTP <- data.frame(lag=demo11tp$lag,CCF=demo11tp$acf)

demo10yrPS <- data.frame(lag=demo10ps$lag,CCF=demo10ps$acf)
demo10yrST <- data.frame(lag=demo10st$lag,CCF=demo10st$acf)
demo10yrTP <- data.frame(lag=demo10tp$lag,CCF=demo10tp$acf)

demo9yrPS <- data.frame(lag=demo9ps$lag,CCF=demo9ps$acf)
demo9yrST <- data.frame(lag=demo9st$lag,CCF=demo9st$acf)
demo9yrTP <- data.frame(lag=demo9tp$lag,CCF=demo9tp$acf)

demo8yrPS <- data.frame(lag=demo8ps$lag,CCF=demo8ps$acf)
demo8yrST <- data.frame(lag=demo8st$lag,CCF=demo8st$acf)
demo8yrTP <- data.frame(lag=demo8tp$lag,CCF=demo8tp$acf)

demo7yrPS <- data.frame(lag=demo7ps$lag,CCF=demo7ps$acf)
demo7yrST <- data.frame(lag=demo7st$lag,CCF=demo7st$acf)
demo7yrTP <- data.frame(lag=demo7tp$lag,CCF=demo7tp$acf)

demo6yrPS <- data.frame(lag=demo6ps$lag,CCF=demo6ps$acf)
demo6yrST <- data.frame(lag=demo6st$lag,CCF=demo6st$acf)
demo6yrTP <- data.frame(lag=demo6tp$lag,CCF=demo6tp$acf)

demo5yrPS <- data.frame(lag=demo5ps$lag,CCF=demo5ps$acf)
demo5yrST <- data.frame(lag=demo5st$lag,CCF=demo5st$acf)
demo5yrTP <- data.frame(lag=demo5tp$lag,CCF=demo5tp$acf)

demo4yrPS <- data.frame(lag=demo4ps$lag,CCF=demo4ps$acf)
demo4yrST <- data.frame(lag=demo4st$lag,CCF=demo4st$acf)
demo4yrTP <- data.frame(lag=demo4tp$lag,CCF=demo4tp$acf)

#now analyse and plot the cross correlation analysis
sine_model <- function(x, A, B, C, D) {
  A * sin(B * x + C) + D
}
library(minpack.lm)
# Fit the model to data
library(readr)


fit_and_plot_model <- function(data, formula, start_list, plot_title, dataset_name) {
  fit <- nlsLM(formula, data = data, start = start_list, control= nls.control(maxiter = 500))
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

### demo PS
demo14yrPSplot <- fit_and_plot_model(demo14yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo14yrPS")
#print(demo14yrPSplot$plot)
print(demo14yrPSplot$metrics)

demo13yrPSplot <- fit_and_plot_model(demo13yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo13yrPS")
#print(demo13yrPSplot$plot)
print(demo13yrPSplot$metrics)

demo12yrPSplot <- fit_and_plot_model(demo12yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo12yrPS")
#print(demo12yrPSplot$plot)
print(demo12yrPSplot$metrics)

demo11yrPSplot <- fit_and_plot_model(demo11yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo11yrPS")
#print(demo11yrPSplot$plot)
print(demo11yrPSplot$metrics)

demo10yrPSplot <- fit_and_plot_model(demo10yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo10yrPS")
#print(demo10yrPSplot$plot)
print(demo10yrPSplot$metrics)

demo9yrPSplot <- fit_and_plot_model(demo9yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo9yrPS")
#print(demo9yrPSplot$plot)
print(demo9yrPSplot$metrics)

demo8yrPSplot <- fit_and_plot_model(demo8yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                  start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                  plot_title = "",
                                  dataset_name = "demo8yrPS")
#print(demo8yrPSplot$plot)
print(demo8yrPSplot$metrics)

demo7yrPSplot <- fit_and_plot_model(demo7yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo7yrPS")
#print(demo7yrPSplot$plot)
print(demo7yrPSplot$metrics)

demo6yrPSplot <- fit_and_plot_model(demo6yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo6yrPS")
#print(demo6yrPSplot$plot)
print(demo6yrPSplot$metrics)

demo5yrPSplot <- fit_and_plot_model(demo5yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo5yrPS")
#print(demo5yrPSplot$plot)
print(demo5yrPSplot$metrics)

demo4yrPSplot <- fit_and_plot_model(demo4yrPS, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo4yrPS")
#print(demo4yrPSplot$plot)
print(demo4yrPSplot$metrics)

### demo ST
demo14yrSTplot <- fit_and_plot_model(demo14yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo14yrST")
#print(demo14yrSTplot$plot)
print(demo14yrSTplot$metrics)

demo13yrSTplot <- fit_and_plot_model(demo13yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo13yrST")
#print(demo13yrSTplot$plot)
print(demo13yrSTplot$metrics)

demo12yrSTplot <- fit_and_plot_model(demo12yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo12yrST")
#print(demo12yrSTplot$plot)
print(demo12yrSTplot$metrics)

demo11yrSTplot <- fit_and_plot_model(demo11yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo11yrST")
#print(demo11yrSTplot$plot)
print(demo11yrSTplot$metrics)

demo10yrSTplot <- fit_and_plot_model(demo10yrST, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo10yrST")
#print(demo10yrSTplot$plot)
print(demo10yrSTplot$metrics)

demo9yrSTplot <- fit_and_plot_model(demo9yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo9yrST")
#print(demo9yrSTplot$plot)
print(demo9yrSTplot$metrics)

demo8yrSTplot <- fit_and_plot_model(demo8yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo8yrST")
#print(demo8yrSTplot$plot)
print(demo8yrSTplot$metrics)

demo7yrSTplot <- fit_and_plot_model(demo7yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "Sec/tert, 7 years",
                                    dataset_name = "demo7yrST")
#print(demo7yrSTplot$plot)
print(demo7yrSTplot$metrics)

demo6yrSTplot <- fit_and_plot_model(demo6yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo6yrST")
#print(demo6yrSTplot$plot)
print(demo6yrSTplot$metrics)

demo5yrSTplot <- fit_and_plot_model(demo5yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo5yrST")
#print(demo5yrSTplot$plot)
print(demo5yrSTplot$metrics)

demo4yrSTplot <- fit_and_plot_model(demo4yrST, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo4yrST")
#print(demo4yrSTplot$plot)
print(demo4yrSTplot$metrics)

### demo TP
demo14yrTPplot <- fit_and_plot_model(demo14yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo14yrTP")
#print(demo14yrTPplot$plot)
print(demo14yrTPplot$metrics)

demo13yrTPplot <- fit_and_plot_model(demo13yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo13yrTP")
#print(demo13yrTPplot$plot)
print(demo13yrTPplot$metrics)

demo12yrTPplot <- fit_and_plot_model(demo12yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo12yrTP")
#print(demo12yrTPplot$plot)
print(demo12yrTPplot$metrics)

demo11yrTPplot <- fit_and_plot_model(demo11yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 0.5, C = -1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo11yrTP")
#print(demo11yrTPplot$plot)
print(demo11yrTPplot$metrics)

demo10yrTPplot <- fit_and_plot_model(demo10yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                     start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                     plot_title = "",
                                     dataset_name = "demo10yrTP")
#print(demo10yrTPplot$plot)
print(demo10yrTPplot$metrics)

demo9yrTPplot <- fit_and_plot_model(demo9yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.1, B = 0.5, C = -1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo9yrTP")
#print(demo9yrTPplot$plot)
print(demo9yrTPplot$metrics)


demo8yrTPplot <- fit_and_plot_model(demo8yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo8yrTP")
#print(demo8yrTPplot$plot)
print(demo8yrTPplot$metrics)

demo7yrTPplot <- fit_and_plot_model(demo7yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.1, B = 0.5, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo7yrTP")
#print(demo7yrTPplot$plot)
print(demo7yrTPplot$metrics)

demo6yrTPplot <- fit_and_plot_model(demo6yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo6yrTP")
#print(demo6yrTPplot$plot)
print(demo6yrTPplot$metrics)

demo5yrTPplot <- fit_and_plot_model(demo5yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo5yrTP")
#print(demo5yrTPplot$plot)
print(demo5yrTPplot$metrics)

demo4yrTPplot <- fit_and_plot_model(demo4yrTP, CCF ~ sine_model(lag, A, B, C, D),
                                    start_list = list(A = 0.5, B = 1, C = 1, D = 0),
                                    plot_title = "",
                                    dataset_name = "demo4yrTP")
#print(demo4yrTPplot$plot)
print(demo4yrTPplot$metrics)

library(ggpubr)
library(cowplot)
# ggarrange(demo8yrPSplot$plot,demo8yrSTplot$plot,demo8yrTPplot$plot, 
#           #demo7yrPSplot$plot,demo7yrSTplot$plot,demo7yrTPplot$plot,
#           demo6yrPSplot$plot,demo6yrSTplot$plot,demo6yrTPplot$plot,
#           #demo5yrPSplot$plot,demo5yrSTplot$plot,demo5yrTPplot$plot,
#           demo4yrPSplot$plot,demo4yrSTplot$plot,demo4yrTPplot$plot,
#           ncol = 3, nrow=3)

prim.sec<-text_grob("Primary\n and\n Secondary", size = 9)
sec.tert<-text_grob("Secondary\n and\n Tertiary", size = 9)
prim.tert<-text_grob("Primary\n and\n Tertiary", size = 9)

pssine<-plot_grid(demo4yrPSplot$plot, demo9yrPSplot$plot, demo14yrPSplot$plot, prim.sec, 
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
stsine<-plot_grid(demo4yrSTplot$plot, demo9yrSTplot$plot, demo14yrSTplot$plot, sec.tert,
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))
ptsine<-plot_grid( demo4yrTPplot$plot, demo9yrTPplot$plot, demo14yrTPplot$plot, prim.tert,
                  ncol = 4, nrow=1, rel_widths=c(0.3,0.3,0.3,0.1))

studyduration<-text_grob("Study duration", size=11)
duration4<-text_grob("   4", size=10)
durationmid<-text_grob("   9", size=10)
durationtop<-text_grob("   14", size=10)
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


allccf.demodata<-plot_grid(durationlabel, durationbar, plots_together_lab,   ncol=1, rel_heights=c(0.1, 0.08, 1.9))

allccf.demodata

#compile fit metrics and make plots
demodf <- bind_rows(demo14yrPSplot$metrics,demo14yrSTplot$metrics,demo14yrTPplot$metrics,
                    demo13yrPSplot$metrics,demo13yrSTplot$metrics,demo13yrTPplot$metrics,
                    demo12yrPSplot$metrics,demo12yrSTplot$metrics,demo12yrTPplot$metrics,
                    demo11yrPSplot$metrics,demo11yrSTplot$metrics,demo11yrTPplot$metrics,
                    demo10yrPSplot$metrics,demo10yrSTplot$metrics,demo10yrTPplot$metrics,
                    demo9yrPSplot$metrics,demo9yrSTplot$metrics,demo9yrTPplot$metrics,
                    demo8yrPSplot$metrics,demo8yrSTplot$metrics,demo8yrTPplot$metrics, 
                    demo7yrPSplot$metrics,demo7yrSTplot$metrics,demo7yrTPplot$metrics,
                    demo6yrPSplot$metrics,demo6yrSTplot$metrics,demo6yrTPplot$metrics,
                    demo5yrPSplot$metrics,demo5yrSTplot$metrics,demo5yrTPplot$metrics,
                    demo4yrPSplot$metrics,demo4yrSTplot$metrics,demo4yrTPplot$metrics)
demodf


library(ggpmisc)
library(tidyr)
demodf$Year <- as.numeric(demodf$Year)
allderived.demodata <- demodf %>%
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

allderived.demodata 

#ok let's try to bring this together into a single figure

timeseriesplot.demodata
all_heat.demodata
allccf.demodata
allderived.demodata

complete_analysis.demodata<-plot_grid(timeseriesplot.demodata, all_heat.demodata, allccf.demodata, allderived.demodata, 
                             ncol=2,  rel_heights=c(1,1,1,1), labels="AUTO", axis="b", align="v")

complete_analysis.demodata

pdf("figures/demofig1.pdf", width = 14.5, height = 9)
complete_analysis.demodata
dev.off()

