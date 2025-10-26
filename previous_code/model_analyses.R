#####################################
# Exploratory analysis of model outputs for LTER sites

library(ggplot2)
library(lattice)

# import the data
mod.lter <- read.csv(file="model_output/model_all_LTER_sites.csv")

str(mod.lter)
levels(as.factor(mod.lter$trophic_level))

boxplot(slope ~ trophic_level, data = mod.lter)
boxplot(p_value ~ trophic_level, data = mod.lter)
boxplot(slope_se ~ trophic_level, data = mod.lter)
boxplot(r_square ~ trophic_level, data = mod.lter)
boxplot(adj_r_square ~ trophic_level, data = mod.lter)

ggplot(mod.lter, aes(start_year, slope, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(start_year, slope_se, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(start_year, p_value, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(start_year, r_square, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(start_year, adj_r_square, colour = trophic_level)) + geom_point() + theme_bw()

ggplot(mod.lter, aes(N_years, slope, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(N_years, slope_se, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(N_years, p_value, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(N_years, r_square, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.lter, aes(N_years, adj_r_square, colour = trophic_level)) + geom_point() + theme_bw()

options(scipen = 999)
mod.lter$p_value <- round(mod.lter$p_value, digits = 10)
mod.lter$intercept_p_value <- round(mod.lter$intercept_p_value, digits = 10)

str(mod.lter)
mod.lter$trophic_level <- as.factor(mod.lter$trophic_level)
levels(mod.lter$trophic_level)
#mod.lter$trophic_level <- factor(mod.lter$trophic_level,
#                                 levels = c("producer", "herbivore", "consumer", "predator"))
mod.lter$N_years <- as.numeric(mod.lter$N_years)
mod.lter$start_year <- as.numeric(mod.lter$start_year)

stripParams <- list(cex=2, lines=1.5)

png("figures/LTER_SITES.png", width = 2000, height = 600, pointsize = 20)
xyplot(mod.lter$N_years ~ mod.lter$slope | mod.lter$trophic_level,
       par.strip.text = stripParams,
       par.settings = list(strip.background = list(col = "lightgray")),
       col = ifelse(mod.lter$p_value <= 0.05,'red','black'),
       pch = 19, cex = 1.2,
       xlab = list(label = "Slope", cex = 2.2),
       ylab = list(label = "Number of years", cex = 2.2),
       scales=list(cex = 1.5),
       abline=c(v=0, lwd = 3))
dev.off()


model.grassmass

## Konza Prairie
mod.konza <- read.csv(file="model_output/model_konza.csv")
options(scipen = 999)
mod.konza$p_value <- round(mod.konza$p_value, digits = 10)
mod.konza$intercept_p_value <- round(mod.konza$intercept_p_value, digits = 10)

str(mod.konza)
mod.konza$trophic_level <- as.factor(mod.konza$trophic_level)
levels(mod.konza$trophic_level)
mod.konza$N_years <- as.numeric(mod.konza$N_years)
mod.konza$start_year <- as.numeric(mod.konza$start_year)

png("figures/Konza_Prairie.png", width = 1600, height = 600, pointsize = 20)
xyplot(mod.konza$N_years ~ mod.konza$slope | mod.konza$trophic_level,
       par.strip.text = stripParams,
       par.settings = list(strip.background = list(col = "lightgray")),
       col = ifelse(mod.konza$p_value <= 0.05,'red','black'),
       pch = 19, cex = 1.2,
       xlab = list(label = "Slope", cex = 2.2),
       ylab = list(label = "Number of years", cex = 2.2),
       scales=list(cex = 1.8),
       abline=c(v=0, lwd = 3))
dev.off()

ggplot(mod.konza, aes(N_years, slope, colour = trophic_level)) + geom_point() + theme_classic()
ggplot(mod.konza, aes(N_years, slope_se, colour = trophic_level)) + geom_point() + theme_classic()
ggplot(mod.konza, aes(N_years, p_value, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.konza, aes(N_years, r_square, colour = trophic_level)) + geom_point() + theme_bw()
ggplot(mod.konza, aes(N_years, adj_r_square, colour = trophic_level)) + geom_point() + theme_bw()


## Hubbard Brook
mod.hbrook <- read.csv(file="model_output/model_hubbard_brook.csv")
options(scipen = 999)
mod.hbrook$p_value <- round(mod.hbrook$p_value, digits = 10)
mod.hbrook$intercept_p_value <- round(mod.hbrook$intercept_p_value, digits = 10)

str(mod.hbrook)
mod.hbrook$trophic_level <- as.factor(mod.hbrook$trophic_level)
levels(mod.hbrook$trophic_level)
mod.hbrook$N_years <- as.numeric(mod.hbrook$N_years)
mod.hbrook$start_year <- as.numeric(mod.hbrook$start_year)

png("figures/Hubbard_Brook.png", width = 1600, height = 600, pointsize = 20)
xyplot(mod.hbrook$N_years ~ mod.hbrook$slope | mod.hbrook$trophic_level,
       par.strip.text = stripParams,
       par.settings = list(strip.background = list(col = "lightgray")),
       col = ifelse(mod.hbrook$p_value <= 0.05,'red','black'),
       pch = 19, cex = 1.2,
       xlab = list(label = "Slope", cex = 2.2),
       ylab = list(label = "Number of years", cex = 2.2),
       scales=list(cex = 1.5),
       abline=c(v=0, lwd = 3))
dev.off()

## North Temperate Lakes
mod.ntlakes <- read.csv(file="model_output/model_north_temperate_lakes.csv")
options(scipen = 999)
mod.ntlakes$p_value <- round(mod.ntlakes$p_value, digits = 10)
mod.ntlakes$intercept_p_value <- round(mod.ntlakes$intercept_p_value, digits = 10)

str(mod.ntlakes)
mod.ntlakes$trophic_level <- as.factor(mod.ntlakes$trophic_level)
levels(mod.ntlakes$trophic_level)
mod.ntlakes$N_years <- as.numeric(mod.ntlakes$N_years)
mod.ntlakes$start_year <- as.numeric(mod.ntlakes$start_year)

png("figures/North_Temperate_Lakes.png", width = 1600, height = 600, pointsize = 20)
xyplot(mod.ntlakes$N_years ~ mod.ntlakes$slope | mod.ntlakes$trophic_level,
       par.strip.text = stripParams,
       par.settings = list(strip.background = list(col = "lightgray")),
       col = ifelse(mod.ntlakes$p_value <= 0.05,'red','black'),
       pch = 19, cex = 1.2,
       xlab = list(label = "Slope", cex = 2.2),
       ylab = list(label = "Number of years", cex = 2.2),
       scales=list(cex = 1.5),
       abline=c(v=0, lwd = 3))
dev.off()


## Santa Barbara Coastal
mod.sbcoastal <- read.csv(file="model_output/model_santa_barbara_coastal.csv")
options(scipen = 999)
mod.sbcoastal$p_value <- round(mod.sbcoastal$p_value, digits = 10)
mod.sbcoastal$intercept_p_value <- round(mod.sbcoastal$intercept_p_value, digits = 10)

str(mod.sbcoastal)
mod.sbcoastal$trophic_level <- as.factor(mod.sbcoastal$trophic_level)
levels(mod.sbcoastal$trophic_level)
mod.sbcoastal$N_years <- as.numeric(mod.sbcoastal$N_years)
mod.sbcoastal$start_year <- as.numeric(mod.sbcoastal$start_year)

png("figures/Santa_Barbara_Coastal.png", width = 1600, height = 600, pointsize = 20)
xyplot(mod.sbcoastal$N_years ~ mod.sbcoastal$slope | mod.sbcoastal$trophic_level,
       par.strip.text = stripParams,
       par.settings = list(strip.background = list(col = "lightgray")),
       col = ifelse(mod.sbcoastal$p_value <= 0.05,'red','black'),
       pch = 19, cex = 1.2,
       xlab = list(label = "Slope", cex = 2.2),
       ylab = list(label = "Number of years", cex = 2.2),
       scales=list(cex = 1.5),
       abline=c(v=0, lwd = 3))
dev.off()


## kathryn 04/04/2023
## remove other columns
# max true and min true are the funky ones
## must be the subset of data im using
## ask what form of data needs to be imputed

out <- multiple_breakups(prim.sbcoastal)
count<-nrow(out)
true_error<-(out[count,5])*(sqrt(out[count, 2]))
true_error

out_2 <- multiple_breakups(lampyrid_es)
count<-nrow(out_2)
true_error<-(out_2[count,5])*(sqrt(out_2[count, 2]))
true_error

prim.sbcoastal <- subset(mod.sbcoastal, trophic_level=="Primary" )
secon.sbcoastal <- subset(mod.sbcoastal, trophic_level=="Secondary")
tert.sbcoastal <- subset(mod.sbcoastal, trophic_level=="Tertiary")

pyramid.sb.prim <- pyramid_plot(prim.sbcoastal, rsq_points=TRUE)
plot(pyramid.sb.prim)
pyramid.sb.secon <- pyramid_plot(secon.sbcoastal, rsq_points=TRUE)
plot(pyramid.sb.secon)
pyramid.sb.tert <- pyramid_plot(tert.sbcoastal, rsq_points=TRUE)
plot(pyramid.sb.tert)

## aaaand may 28th update. all combos. lets add titles too
## nevermind. remember this doesnt work? go back to the tritrophic_analysis script


stability_time(prim.sbcoastal, min_percent = 95, error_multiplyer = 1)
stability_time(secon.sbcoastal, min_percent = 95, error_multiplyer = 1)
stability_time(tert.sbcoastal, min_percent = 95, error_multiplyer = 1)

max_true(prim.sbcoastal)
wrongness <- wrongness_plot(prim.sbcoastal)
plot(wrongness)





