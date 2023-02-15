# Parameter estimation for ABM using control section of SB33 data

rm(list = ls())

path_fun <- '../../ABM_functions/'
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
library('openxlsx')
library('Rcpp')

sourceCpp('../../C++/CTM_cpp.cpp')

meas_C <- as.data.frame(read_excel("../data/dat_simple.xlsx", sheet = "C"))

meas_C$time <- meas_C$days
meas_C$slurry_mass <- meas_C$mass

slurry_mass_dat_C <- meas_C[!is.na(meas_C$mass), c('date', 'time', 'slurry_mass')]
wash_dates_C <-ymd_hms(c("2020-08-13 07:00:00 UTC", "2020-11-05 02:00:00 UTC", "2021-01-28 02:00:00 UTC"))
slurry_mass_dat_C$wash_water <- 0
slurry_mass_dat_C$wash_water[slurry_mass_dat_C$date %in% wash_dates_C] <- 3500
slurry_mass_dat_C$date <- NULL
slurry_mass_dat_C <- slurry_mass_dat_C[!duplicated(slurry_mass_dat_C),]

temp_dat_C <- meas_C[!is.na(meas_C$temp), c('time', 'temp')]
temp_dat_C <- data.frame(data.table::data.table(temp_dat_C)[, .(time = mean(time), temp = mean(temp)), by = .(timebin = floor(time / 1))])[,-c(1)]

wthr_pars = list(temp_air_C = 20, RH = 90, rain = 0, pres_kpa = 101, rs = 10)
evap_pars = list(evap = 0.5 * et(temp_C = wthr_pars$temp_air_C, pres_kpa = wthr_pars$pres_kpa, rs = wthr_pars$rs))

grz_pars = list(graze_start = "may",
                graze_days = 0,
                graze_hours = 0)
mic_pars = list(ks_SO4 = 0.0067,
                ki_H2S_meth = 0.23,
                ki_H2S_sr = 0.25,
                alpha_opt = c(xa_dead= 0.02, starch = 0.2, CF = 0.02, CP = 0.02, urea = 70, NDF = 0.02, iNDF = 0, VSd = 0.02),
                alpha_T_min = c(xa_dead= 0, starch = 0, CF = 0, CP = 0, urea = 0, NDF = 0, iNDF = 0, VSd = 0),
                alpha_T_opt = c(xa_dead= 50, starch = 50, CF = 50, CP = 50, urea = 50, NDF = 50, iNDF = 50, VSd = 50),
                alpha_T_max = c(xa_dead= 60, starch = 60, CF = 60, CP = 60, urea = 60, NDF = 60, iNDF = 60, VSd = 60))

chem_pars = list(COD_conv = c(CH4 = 0.2507, xa_dead = 0.73, NDF = 0.84, iNDF = 0.65, starch = 0.85, 
                              CF = 0.35, CP = 0.65, VFA = 0.93, S = 0.5015, VS = 0.69, CO2_anaer = 0.53, CO2_aer = 1.1, CO2_sr = 1.2, CO2_ureo = 1.57,
                              N_CP = 0.1014, C_xa_dead = 0.358, C_NDF = 0.376, C_iNDF = 0.358
                              , C_starch = 0.377, C_CF = 0.265, C_CP = 0.359 , C_VFA = 0.374, C_VSd = 0.344, C_N_urea = 0.429), 
                 kl = c(NH3 = 52, NH3_floor = 22, H2S = 0.02)) 
# C parms
# make parms data frame for looping through in sensitivity analysis
norm <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
fraction <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.3, 1.5, 1.7, 1.9)

ks_coefficient <- c(fraction, rep(norm, 6))
xa_fresh <- c(norm, fraction, rep(norm, 5))
yield <- c(rep(norm, 2), fraction, rep(norm, 4))
decay_rate <- c(rep(norm, 3), fraction, rep(norm, 3))
resid_enrich <- c(rep(norm, 4), fraction, rep(norm, 2))
alpha_opt <- c(rep(norm, 5), fraction, norm)
qhat_opt <- c(rep(norm, 6), fraction)

pars_dat <- data.frame(ks_coefficient = ks_coefficient, xa_fresh = xa_fresh, 
                       yield = yield, decay_rate = decay_rate, resid_enrich = resid_enrich,
                       alpha_opt = alpha_opt, qhat_opt = qhat_opt) 

l <- nrow(pars_dat)

norm_alpha <- rep(c(1, 1, 5, 5), each = l)
norm_qhat <- rep(c(1, 0.4, 1, 0.4), each = l)

pars_dat_all <- data.frame(rbind(pars_dat, pars_dat, pars_dat, pars_dat))
pars_dat_all$alpha_opt <- pars_dat_all$alpha_opt * norm_alpha
pars_dat_all$qhat_opt <- pars_dat_all$qhat_opt * norm_qhat

pred_all <- NULL

for (i in 1:nrow(pars_dat_all)){
  
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
                    resid_enrich = 0.9 * pars_dat_all$resid_enrich[i],
                    slopes = c(urea = NA, slurry_prod_rate = NA),
                    scale = c(ks_coefficient = pars_dat_all[i, 'ks_coefficient'] * 1, 
                              qhat_opt = pars_dat_all[i, 'qhat_opt'] * 1, 
                              xa_fresh = pars_dat_all[i, 'xa_fresh'] * 1, 
                              yield = pars_dat_all[i, 'yield'] * 1, 
                              alpha_opt = pars_dat_all[i, 'alpha_opt'] * 1))
  
  grp_pars = list(grps = c('m0','m1','m2', 'sr1'),
                  yield = c(default = 0.05, sr1 = 0.065),
                  xa_fresh = c(default =  pars_dat_all[i, 'decay_rate'] * 0.02),
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
  
  man_pars_C = list(conc_fresh = list(S2 = 0.01, urea = 2.4, SO4 = 0.2, TAN = 0.63, starch = 0, 
                                      VFA = 2.83, xa_dead = 0, CF = 0, CP = 0, NDF = 0, iNDF = 15, VSd = 73.8, VSd_A = 44.4), pH = 6.88, dens = 1000)
  
 
  ################### adjust for optimizer and exclude period 4 ###################

  days_sim <- max(meas_C$days[meas_C$period == 3], na.rm = T)
  meas_C$day <- ceiling(meas_C$days)
  dat_C <- as.data.frame(summarise(group_by(meas_C, day), CH4_emis_rate = mean(CH4_emis_rate, na.rm = TRUE), VFA_conc = mean(vfa / 1000 * 0.93, na.rm = TRUE)))
  dat_C <- dat_C[dat_C$day <= ceiling(days_sim), ] 
  times_C <- dat_C$day
  pred_C <- abm(days = max(times_C), times = times_C, man_pars = man_pars_C, mng_pars = mng_pars_C, grp_pars = grp_pars, evap_pars = evap_pars, wthr_pars = wthr_pars)
  pred_C <- mean(pred_C$CH4_emis_rate, na.rm = T)
  pred_all <- c(pred_all, pred_C)
}

pred_dat <- cbind(pars_dat_all, CH4 = pred_all)

pred_dat_mod <- pred_dat
l <- nrow(pred_dat_mod)
pred_dat_mod$mode <- c(rep(1,l/4), rep(2, l/4), rep(3, l/4), rep(4, l/4)) 
pred_dat_mod[pred_dat_mod == 1] <- NA
pred_dat_mod[pred_dat_mod == 0.4] <- NA
pred_dat_mod[pred_dat_mod == 5] <- NA

pred_dat_mod$mode[is.na(pred_dat_mod$mode)] <- 1

pred_dat_long <- pred_dat_mod %>% pivot_longer(cols = ks_coefficient:qhat_opt, names_to = 'parm', values_to = 'value') %>%
  filter(!is.na(value))

norm_alpha <- rep(c(1, 1, 5, 5), each = l)
norm_qhat <- rep(c(1, 0.4, 1, 0.4), each = l)

pred_dat_long


divide <- rep(1, nrow(pred_dat_long))
divide[which((pred_dat_long$mode == 3 | pred_dat_long$mode == 4) & pred_dat_long$parm == 'alpha_opt')] <- 5
divide[which((pred_dat_long$mode == 2 | pred_dat_long$mode == 4) & pred_dat_long$parm == 'qhat_opt')] <- 0.4

pred_dat_long$divide <- NULL
pred_dat_long$divide <- divide

pred_dat_long <- mutate(pred_dat_long, perc_value = (value/divide -1) * 100)

pars_long <- pred_dat[,-c(8)]
rows_norm_CH4 <- c(which(rowSums(pars_long) %in% c(7, 11, 6.2, 10.2)))
norm_CH4 <- rep(unique(pred_dat[rows_norm_CH4, 'CH4']), each = nrow(pred_dat_long)/4)
pred_dat_long$norm_CH4 <- norm_CH4
pred_dat_long <- mutate(pred_dat_long, perc_CH4 = (CH4 / norm_CH4 -1) *100) 

write.xlsx(pred_dat_long, '../data/sensitivity.xlsx')
