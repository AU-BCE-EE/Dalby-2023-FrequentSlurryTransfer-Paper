# Parameter estimation for ABM using control section of SB33 data

rm(list = ls())

path_fun <- '../ABM_functions/'
ff <- list.files(path = path_fun) 
ff <- ff[ff != 'x.R']
for (i in ff) source(paste0(path_fun,i))

library('readxl')
library('gridExtra')
library('svglite')
library('ggplot2')
library('lubridate')
library('dplyr')
library('tidyr')

meas_C <- as.data.frame(read_excel("../data/dat_simple.xlsx", sheet = "C"))
meas_FF <- as.data.frame(read_excel("../data/dat_simple.xlsx", sheet = "FF"))
meas_SF <- as.data.frame(read_excel("../data/dat_simple.xlsx", sheet = "SF"))
meas_ST <- as.data.frame(read_excel("../data/dat_simple.xlsx", sheet = "ST"))

per <- read.csv('../data/periods.csv')

opt_pars <- as.data.frame(read_excel("../data/opt_pars.xlsx"))

# Some cleaning
per$date.time.in <- parse_date_time(paste(per$date.in, '12:00'), orders = 'dmyHM')
per$date.time.out <- parse_date_time(paste(per$date.out, '12:00'), orders = 'dmyHM')
per$date.in <- dmy(per$date.in)
per$date.out <- dmy(per$date.out)

# C slurry_mass data frame
meas_C$time <- meas_C$days
meas_C$slurry_mass <- meas_C$mass

slurry_mass_dat_C <- meas_C[!is.na(meas_C$mass), c('date', 'time', 'slurry_mass')]
wash_dates_C <-ymd_hms(c("2020-08-13 07:00:00 UTC", "2020-11-05 02:00:00 UTC", "2021-01-28 02:00:00 UTC"))
slurry_mass_dat_C$wash_water <- 0
slurry_mass_dat_C$wash_water[slurry_mass_dat_C$date %in% wash_dates_C] <- 3500
slurry_mass_dat_C$date <- NULL
slurry_mass_dat_C <- slurry_mass_dat_C[!duplicated(slurry_mass_dat_C),]

temp_dat_C <- meas_C[!is.na(meas_C$temp), c('time', 'temp')]

# FF slurry mass data frame
meas_FF$time <- meas_FF$days
meas_FF$slurry_mass <- meas_FF$mass

slurry_mass_dat_FF <- meas_FF[!is.na(meas_FF$mass), c('date', 'time', 'slurry_mass')]
wash_dates_FF <-ymd_hms(c("2020-08-13 05:00:00 UTC", "2020-11-05 04:00:00 UTC", "2021-01-28 04:00:00 UTC"))
slurry_mass_dat_FF$wash_water <- 0
slurry_mass_dat_FF$wash_water[slurry_mass_dat_FF$date %in% wash_dates_FF] <- 3500
slurry_mass_dat_FF$date <- NULL
slurry_mass_dat_FF <- slurry_mass_dat_FF[!duplicated(slurry_mass_dat_FF),]

temp_dat_FF <- meas_FF[!is.na(meas_FF$temp), c('time', 'temp')]

# SF slurry mass data frame
meas_SF$time <- meas_SF$days
meas_SF$slurry_mass <- meas_SF$mass

slurry_mass_dat_SF <- meas_SF[!is.na(meas_SF$mass), c('date', 'time', 'slurry_mass')]
wash_dates_SF <-ymd_hms(c("2020-08-13 07:00:00 UTC", "2020-11-04 07:00:00 UTC", "2021-01-27 07:00:00 UTC"))
slurry_mass_dat_SF$wash_water <- 0
slurry_mass_dat_SF$wash_water[slurry_mass_dat_SF$date %in% wash_dates_SF] <- 3500
slurry_mass_dat_SF$date <- NULL
slurry_mass_dat_SF <- slurry_mass_dat_SF[!duplicated(slurry_mass_dat_SF),]

temp_dat_SF <- meas_FF[!is.na(meas_FF$temp), c('time', 'temp')]

# FF slurry mass data frame

meas_ST$time <- meas_ST$days
meas_ST$slurry_mass <- meas_ST$mass

slurry_mass_dat_ST <- meas_ST[!is.na(meas_ST$mass), c('date', 'time', 'slurry_mass')]
wash_dates_ST <-ymd_hms(c("2020-08-13 02:00:00 UTC", "2020-11-05 02:00:00 UTC", "2021-01-28 02:00:00 UTC"))
slurry_mass_dat_ST$wash_water <- 0
slurry_mass_dat_ST$wash_water[slurry_mass_dat_ST$date %in% wash_dates_ST] <- 3500
slurry_mass_dat_ST$date <- NULL
slurry_mass_dat_ST <- slurry_mass_dat_ST[!duplicated(slurry_mass_dat_ST),]

temp_dat_ST <- meas_FF[!is.na(meas_FF$temp), c('time', 'temp')]


########################## Setup model parameters ##############################

wthr_pars = list(temp_air_C = 20, RH = 90, rain = 0, pres_kpa = 101, rs = 10)
evap_pars = list(evap = 0.5 * et(temp_C = wthr_pars$temp_air_C, pres_kpa = wthr_pars$pres_kpa, rs = wthr_pars$rs))
grz_pars = list(graze_start = "May",
                graze_days = 0,
                graze_hours = 0)
grp_pars = list(grps = c('m0','m1','m2', 'sr1'),
                yield = c(default = 0.05, sr1 = 0.065),
                xa_fresh = c(default = 0.02),
                xa_init = c(all = 0.0001),
                decay_rate = c(all = 0.02),
                ks_coefficient = c(default = 1, sr1 = 0.4),
                qhat_opt = c(m0 = 1.5, m1 = 3.6, m2 = 5.6 , m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8),
                T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                T_min = c(m0 = 0, m1 = 10, m2 = 10, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                ki_NH3_min = c(all = 0.015),
                ki_NH3_max = c(all = 0.13),
                ki_NH4_min = c(all = 2.7),
                ki_NH4_max = c(all = 4.8),
                pH_upr = c(all = 8.0),
                pH_lwr = c(default = 6.0))
mic_pars = list(ks_SO4 = 0.0067,
                ki_H2S_meth = 0.23,
                ki_H2S_sr = 0.25,
                alpha_opt = c(xa_dead= 0.02, urea = 70, VSd = 0.02),
                alpha_T_min = c(xa_dead= 0, urea = 0, VSd = 0),
                alpha_T_opt = c(xa_dead= 50, urea = 50, VSd = 50),
                alpha_T_max = c(xa_dead= 60, urea = 60, VSd = 60))
conc_fresh = list(S2 = 0.01, urea = 2.4, SO4 = 0.2, TAN = 0.63,
                  VFA = 2.83, xa_dead = 0, VSd = 73.8)
# C parms
mng_pars_C = list(slurry_prod_rate = 0,  
                  slurry_mass = slurry_mass_dat_C,          
                  storage_depth = 0.6,        
                  resid_depth = 0.035,         
                  floor_area = 22.09,
                  area = 22.09, 
                  empty_int = 35,
                  temp_C = temp_dat_C,
                  wash_water = 0,
                  wash_int = NA,
                  rest_d = 7,
                  RH = 90, 
                  cover = NA,
                  resid_enrich = 0.9,
                  slopes = c(urea = NA, slurry_prod_rate = NA),
                  scale = c(ks_coefficient = opt_pars$scale.ks_coefficient, qhat_opt = opt_pars$scale.qhat_opt, xa_fresh = opt_pars$scale.xa_fresh, yield = 1, alpha_opt = opt_pars$scale.alpha_opt))

man_pars_C = list(conc_fresh = conc_fresh, pH = 6.88, dens = 1000)

# FF parms

mng_pars_FF = list(slurry_prod_rate = 0,  
                  slurry_mass = slurry_mass_dat_FF,          
                  storage_depth = 0.6,        
                  resid_depth = 0.035,         
                  floor_area = 22.09,
                  area = 22.09, 
                  empty_int = 35,
                  temp_C = temp_dat_FF,
                  wash_water = 0,
                  wash_int = NA,
                  rest_d = 7,
                  RH = 90, 
                  cover = NA,
                  resid_enrich = 0.9,
                  slopes = c(urea = NA, slurry_prod_rate = NA),
                  scale = c(ks_coefficient = opt_pars$scale.ks_coefficient, qhat_opt = opt_pars$scale.qhat_opt, xa_fresh = opt_pars$scale.xa_fresh, yield = 1, alpha_opt = opt_pars$scale.alpha_opt))

man_pars_FF = list(conc_fresh = conc_fresh, pH = 6.82, dens = 1000)

# SF parms
mng_pars_SF = list(slurry_prod_rate = 0,  
                  slurry_mass = slurry_mass_dat_SF,          
                  storage_depth = 0.6,        
                  resid_depth = 0.035,         
                  floor_area = 22.09,
                  area = 22.09, 
                  empty_int = 35,
                  temp_C = temp_dat_FF,
                  wash_water = 0,
                  wash_int = NA,
                  rest_d = 7,
                  RH = 90, 
                  cover = NA,
                  resid_enrich = 0.0, # set to 0, because there is no containier and no residuals in the pipes. 
                  slopes = c(urea = NA, slurry_prod_rate = NA),
                  scale = c(ks_coefficient = opt_pars$scale.ks_coefficient, qhat_opt = opt_pars$scale.qhat_opt, xa_fresh = opt_pars$scale.xa_fresh, yield = 1, alpha_opt = opt_pars$scale.alpha_opt))

man_pars_SF = list(conc_fresh = conc_fresh, pH = 6.75, dens = 1000)
# ST parms

mng_pars_ST = list(slurry_prod_rate = 0,  
                   slurry_mass = slurry_mass_dat_ST,          
                   storage_depth = 0.6,        
                   resid_depth = 0.035,         
                   floor_area = 22.09,
                   area = 22.09, 
                   empty_int = 35,
                   temp_C = temp_dat_FF,
                   wash_water = 0,
                   wash_int = NA,
                   rest_d = 7,
                   RH = 90, 
                   cover = NA,
                   resid_enrich = 0.9,
                   slopes = c(urea = NA, slurry_prod_rate = NA),
                   scale = c(ks_coefficient = opt_pars$scale.ks_coefficient, qhat_opt = opt_pars$scale.qhat_opt, xa_fresh = opt_pars$scale.xa_fresh, yield = 1, alpha_opt = opt_pars$scale.alpha_opt))

man_pars_ST = list(conc_fresh = conc_fresh, pH = 6.83, dens = 1000)



################### adjust for optimizer and exclude period 4 ###################
# C plot
days_sim <- max(meas_C$days[meas_C$period == 3], na.rm = T)
meas_C$day <- ceiling(meas_C$days)
dat_C <- as.data.frame(summarise(group_by(meas_C, day), CH4_emis_rate = mean(CH4_emis_rate, na.rm = TRUE), VFA_conc = mean(vfa / 1000 * 0.93, na.rm = TRUE)))
dat_C <- dat_C[dat_C$day <= ceiling(days_sim), ] 
times_C <- dat_C$day
pred_C <- abm(days = max(times_C), times = times_C, man_pars = man_pars_C, mng_pars = mng_pars_C, grp_pars = grp_pars, evap_pars = evap_pars, wthr_pars = wthr_pars)[, c('time', 'CH4_emis_rate', 'VFA_conc', 'slurry_mass')]

# FF plot 
days_sim <- max(meas_FF$days[meas_FF$period == 3], na.rm = T)
meas_FF$day <- ceiling(meas_FF$days)
dat_FF <- as.data.frame(summarise(group_by(meas_FF, day), CH4_emis_rate = mean(CH4_emis_rate, na.rm = TRUE), VFA_conc = mean(vfa / 1000 * 0.93, na.rm = TRUE)))
dat_FF <- dat_FF[dat_FF$day <= ceiling(days_sim), ] 
times_FF <- dat_FF$day
pred_FF <- abm(days = max(times_FF), times = times_FF, man_pars = man_pars_FF, mng_pars = mng_pars_FF, grp_pars = grp_pars, evap_pars = evap_pars, wthr_pars = wthr_pars)[, c('time', 'CH4_emis_rate', 'VFA_conc', 'slurry_mass')]

# SF plot 
days_sim <- max(meas_SF$days[meas_SF$period == 3], na.rm = T)
meas_SF$day <- ceiling(meas_SF$days)
dat_SF <- as.data.frame(summarise(group_by(meas_SF, day), CH4_emis_rate = mean(CH4_emis_rate, na.rm = TRUE), VFA_conc = mean(vfa / 1000 * 0.93, na.rm = TRUE)))
dat_SF <- dat_SF[dat_SF$day <= ceiling(days_sim), ] 
times_SF <- dat_SF$day
pred_SF <- abm(days = max(times_SF), times = times_SF, man_pars = man_pars_SF, mng_pars = mng_pars_SF, grp_pars = grp_pars, evap_pars = evap_pars, wthr_pars = wthr_pars)[, c('time', 'CH4_emis_rate', 'VFA_conc', 'slurry_mass')]

# ST 
days_sim <- max(meas_ST$days[meas_ST$period == 3], na.rm = T)
meas_ST$day <- ceiling(meas_ST$days)
dat_ST <- as.data.frame(summarise(group_by(meas_ST, day), CH4_emis_rate = mean(CH4_emis_rate, na.rm = TRUE), VFA_conc = mean(vfa / 1000 * 0.93, na.rm = TRUE)))
dat_ST <- dat_ST[dat_ST$day <= ceiling(days_sim), ] 
times_ST <- dat_ST$day
pred_ST <- abm(days = max(times_ST), times = times_ST, man_pars = man_pars_ST, mng_pars = mng_pars_ST, grp_pars = grp_pars, evap_pars = evap_pars, wthr_pars = wthr_pars)[, c('time', 'CH4_emis_rate', 'VFA_conc', 'slurry_mass')]

pred_C$treat <- "Control (C)"
pred_FF$treat <- "Weekly flushing (WF)"
pred_SF$treat <- "Slurry funnels (SF)"
pred_ST$treat <- "Slurry trays (ST)"

pred_all <- rbind(pred_C, pred_FF, pred_SF, pred_ST)
pred_all$type <- "pred"
pred_all <- pred_all[-which(pred_all$time %in% 162:165),]

dat_C$treat <- "Control (C)"
dat_FF$treat <- "Weekly flushing (WF)"
dat_SF$treat <- "Slurry funnels (SF)"
dat_ST$treat <- "Slurry trays (ST)" 


dat_all <- rbind(dat_C, dat_FF, dat_SF, dat_ST)
dat_all <- rename(dat_all, time = day)
dat_all$type <- "dat"

dat_pred <- rbind(pred_all, dat_all)
dat_pred$CH4_emis_rate_norm <- dat$CH4_emis_rate / dat$slurry_mass
dat_pred_long <- pivot_longer(dat_pred, c(CH4_emis_rate, VFA_conc, CH4_emis_rate_norm), names_to = 'compound', values_to = 'value')
dat_pred_long$type <- as.factor(dat_pred_long$type)

vline <- as.numeric(difftime(c(per$date.in, per$date.out), min(per$date.in), units = 'days'))
per$day.in <- as.numeric(difftime(per$date.in, min(per$date.in), units = 'days'))
per$day.out <- as.numeric(difftime(per$date.out, min(per$date.in), units = 'days'))
vline <- vline[vline <= per$day.out[3]]

# Sort out NAs in between periods
# NTS: In dat_stacked.csv there is overlap between periods!
# NTS: These days are arbitrary
# NTS: We need definitive timing on periods applied everywhere!
btwn.per <- (dat_pred_long$time > per$day.out[1] &  dat_pred_long$time < per$day.in[2]) |
            (dat_pred_long$time > per$day.out[2] & dat_pred_long$time < per$day.in[3]) |
            (dat_pred_long$time > per$day.out[3])
dat_pred_long$value[btwn.per] <- NA

dat_pred_long <- subset(dat_pred_long, time <= per$day.out[3])

# Drop NAs
dat_pred_long <- subset(dat_pred_long, !is.na(value))

# Duplicate NAs in between periods for VFA data and some CH4 in control that apparently don't have any rows
ipna <- expand.grid(time = vline, 
                    treat = unique(dat_pred_long$treat), 
                    type = unique(dat_pred_long$type), 
                    compound = unique(dat_pred_long$compound), 
                    value = NA)
dat_pred_long <- rbind(dat_pred_long, ipna)

sub_dat_CH4 <- filter(dat_pred_long, type == 'dat', compound == 'CH4_emis_rate')
sub_pred_CH4 <- filter(dat_pred_long, type == 'pred', compound == 'CH4_emis_rate')

sub_dat_VFA <- filter(dat_pred_long, type == 'dat', compound == 'VFA_conc')
sub_pred_VFA <- filter(dat_pred_long, type == 'pred', compound == 'VFA_conc')

sub_dat_CH4$treat <- factor(sub_dat_CH4$treat, levels = c('Control (C)' = 'Control (C)', 'Weekly flushing (WF)' = 'Weekly flushing (WF)', 
                                                          'Slurry funnels (SF)' = 'Slurry funnels (SF)', 'Slurry trays (ST)' = 'Slurry trays (ST)'))
sub_pred_CH4$treat <- factor(sub_pred_CH4$treat, levels = c('Control (C)' = 'Control (C)', 'Weekly flushing (WF)' = 'Weekly flushing (WF)', 
                                                          'Slurry funnels (SF)' = 'Slurry funnels (SF)', 'Slurry trays (ST)' = 'Slurry trays (ST)'))
sub_dat_VFA$treat <- factor(sub_dat_VFA$treat, levels = c('Control (C)' = 'Control (C)', 'Weekly flushing (WF)' = 'Weekly flushing (WF)', 
                                                           'Slurry funnels (SF)' = 'Slurry funnels (SF)', 'Slurry trays (ST)' = 'Slurry trays (ST)'))
sub_pred_VFA$treat <- factor(sub_pred_VFA$treat, levels = c('Control (C)' = 'Control (C)', 'Weekly flushing (WF)' = 'Weekly flushing (WF)', 
                                                           'Slurry funnels (SF)' = 'Slurry funnels (SF)', 'Slurry trays (ST)' = 'Slurry trays (ST)'))

p1 <- ggplot() + 
  geom_point(data = sub_dat_CH4, aes(x = time, y = value), size = 0.4) + 
  geom_line(data = sub_pred_CH4, aes(x = time, y = value), size = 0.4, color = 'red') +
  facet_grid(treat~., scales = "free_y") + theme_bw() + theme(legend.position="none") + 
  labs(x = 'Days', y = expression('Methane emission rate'~(g~d^'-1')), colour = '') + 
  geom_vline(xintercept = vline, lty = 2, col = 'gray65')

p2 <- ggplot() + 
  geom_point(data = sub_dat_VFA, aes(x = time, y = value), size = 0.4) + 
  geom_line(data = sub_pred_VFA, aes(x = time, y = value), size = 0.4, color = 'red') +
  facet_grid(treat~.) + theme_bw() + theme(legend.position="none") + 
  labs(x = 'Days', y = expression('VFA concentration'~(gCOD~kg^'-1')), colour = '') + scale_color_manual(name = "type", values=c("black", "red")) + 
  geom_vline(xintercept = vline, lty = 2, col = 'gray65')

svglite("../figures/fig_valid.svg", width = 17/2.54, height = 19/2.54)
  fig_valid <- (grid.arrange(p1,p2, ncol=2))
dev.off()

png('../figures/fig_valid.png',  width = 17/2.54, height = 19/2.54, units = 'in', res = 600)
grid::grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2)))
dev.off()

pdf('../figures/fig_valid.pdf',   width = 17/2.54, height = 19/2.54)
grid::grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2)))
dev.off()

  
reduction_FF_dat <- (1 - mean(dat_FF$CH4_emis_rate[!is.na(dat_FF$CH4_emis_rate)]) / mean(dat_C$CH4_emis_rate[!is.na(dat_C$CH4_emis_rate)])) * 100
reduction_FF_pred <- (1 - mean(pred_FF$CH4_emis_rate[!is.na(dat_FF$CH4_emis_rate)], na.rm = T) / mean(pred_C$CH4_emis_rate[!is.na(dat_C$CH4_emis_rate)])) * 100

reduction_SF_dat <- (1 - mean(dat_SF$CH4_emis_rate[!is.na(dat_SF$CH4_emis_rate)]) / mean(dat_C$CH4_emis_rate[!is.na(dat_C$CH4_emis_rate)])) * 100
reduction_SF_pred <- (1 - mean(pred_SF$CH4_emis_rate[!is.na(dat_SF$CH4_emis_rate)], na.rm = T) / mean(pred_C$CH4_emis_rate[!is.na(dat_C$CH4_emis_rate)])) * 100

reduction_ST_dat <- (1 - mean(dat_ST$CH4_emis_rate[!is.na(dat_ST$CH4_emis_rate)]) / mean(dat_C$CH4_emis_rate[!is.na(dat_C$CH4_emis_rate)])) * 100
reduction_ST_pred <- (1 - mean(pred_ST$CH4_emis_rate[!is.na(dat_ST$CH4_emis_rate)], na.rm = T) / mean(pred_C$CH4_emis_rate[!is.na(dat_C$CH4_emis_rate)])) * 100


mean(pred_C$CH4_emis_rate[!is.na(dat_C$CH4_emis_rate)], na.rm = T)
mean(pred_FF$CH4_emis_rate[!is.na(dat_FF$CH4_emis_rate)], na.rm = T)
mean(pred_SF$CH4_emis_rate[!is.na(dat_SF$CH4_emis_rate)], na.rm = T)
mean(pred_ST$CH4_emis_rate[!is.na(dat_ST$CH4_emis_rate)], na.rm = T)

mean(dat_C$VFA_conc[!is.na(dat_C$VFA_conc)])
mean(dat_FF$VFA_conc[!is.na(dat_FF$VFA_conc)])
mean(dat_SF$VFA_conc[!is.na(dat_SF$VFA_conc)])
mean(dat_ST$VFA_conc[!is.na(dat_ST$VFA_conc)])

mean(pred_C$VFA_conc[which(!is.na(dat_C$VFA_conc))-1])
mean(pred_FF$VFA_conc[which(!is.na(dat_FF$VFA_conc))-1])
mean(pred_SF$VFA_conc[which(!is.na(dat_SF$VFA_conc))-1])
mean(pred_ST$VFA_conc[which(!is.na(dat_ST$VFA_conc))-1])

