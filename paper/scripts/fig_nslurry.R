rm(list = ls())

library("lubridate")
library("readxl")
library("ggplot2")
library("svglite")
library("dplyr")
library('tidyr')

path_dat <- '../data/dat_comp.xlsx'
dat_analysis <- read_excel(path_dat, sheet = "analysis")

dat_long <- mutate(group_by(dat_analysis, period), time = as.numeric(date - min(date)) / (60 * 60 * 24) + 7)
dat_long <- pivot_longer(dat_long, c(TN, TAN, VS, DM, pH), names_to = 'compound', values_to = 'value')

dat_long$treatment <- gsub("control", "Control (C)", gsub("frequentflushing", "Weekly flushing (WF)", gsub("slurryfunnels", "Slurry funnels (SF)", gsub("slurrytrays", "Slurry trays (ST)", dat_long$treatment)))) 

plot_dat = dat_long %>% filter(compound %in% c("TN", "TAN"))

plot_dat$period <- as.factor(plot_dat$period)

p <- ggplot(plot_dat, aes(x = time, y = value)) + 
  geom_point(aes(col = period), size = 0.5) + 
  facet_grid(cols = vars(treatment), rows = vars(compound), scales = "free_y") + ylab(bquote(Concentration~(gN~kg^-1))) + xlab("Days") + theme_bw()

svglite("../figures/Nslurry.svg", width = 19/2.54, height = 8/2.54)
  p
dev.off()

png('../figures/Nslurry.png', width = 19/2.54, height = 8/2.54, units = 'in', res = 600)
grid::grid.draw(p)
dev.off()

pdf('../figures/Nslurry.pdf', width = 19/2.54, height = 8/2.54)
grid::grid.draw(p)
dev.off()

## stats ##

TN_dat <- plot_dat %>% filter(compound == 'TN')

p_slope_TN <- summary(lm(value~time, TN_dat))$coefficients[2,4]
r_TN_time <- cor(TN_dat$time[!is.na(TN_dat$value)], TN_dat$value[!is.na(TN_dat$value)])

TAN_dat <- plot_dat %>% filter(compound == 'TAN')

p_slope_TAN <- summary(lm(value~time, TAN_dat))$coefficients[2,4]
r_TAN_time <- cor(TAN_dat$time[!is.na(TAN_dat$value)], TAN_dat$value[!is.na(TAN_dat$value)])

