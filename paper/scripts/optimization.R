# Parameter estimation for ABM using control section of SB33 data

rm(list = ls())
path_fun <- '../ABM_functions/' 
ff <- list.files(path = path_fun) 
ff <- ff[ff != 'x.R']
for (i in ff) source(paste0(path_fun,i))

abm_packages()

meas <- as.data.frame(read_excel("../data/dat_simple.xlsx", sheet = "C"))

# handle variable slurry and temperature

# C slurry_mass data frame

meas$time <- meas$days
meas$slurry_mass <- meas$mass

slurry_mass_dat <- meas[!is.na(meas$mass), c('date', 'time', 'slurry_mass')]
wash_dates <-ymd_hms(c("2020-08-13 07:00:00 UTC", "2020-11-05 02:00:00 UTC", "2021-01-28 02:00:00 UTC"))
slurry_mass_dat$wash_water <- 0
slurry_mass_dat$wash_water[slurry_mass_dat$date %in% wash_dates] <- 3500
slurry_mass_dat$date <- NULL
slurry_mass_dat <- slurry_mass_dat[!duplicated(slurry_mass_dat),]

temp_dat <- meas[!is.na(meas$temp), c('time', 'temp')]
# Averge temperature by 10 d periods to increase speed.
temp_dat2 <- data.frame(data.table::data.table(temp_dat)[, .(time = mean(time), temp = mean(temp)), by = .(timebin = floor(time / 1))])

# C parms

wthr_pars = list(temp_air_C = 20, RH = 90, rain = 0, pres_kpa = 101, rs = 10)
evap_pars = list(evap = 0.5 * et(temp_C = wthr_pars$temp_air_C, pres_kpa = wthr_pars$pres_kpa, rs = wthr_pars$rs))
mng_pars = list(slurry_prod_rate = 0,  
                slurry_mass = slurry_mass_dat,          
                storage_depth = 0.6,        
                resid_depth = 0.035,         
                floor_area = 22.09,
                area = 22.09, 
                empty_int = 35,
                temp_C = temp_dat,
                wash_water = 0,
                wash_int = NA,
                rest_d = 7,
                RH = 90, 
                cover = NA,
                resid_enrich = 0.9,
                slopes = c(urea = NA, slurry_prod_rate = NA),
                scale = c(ks_coefficient = 1, qhat_opt = 0.4, xa_fresh = 4, yield = 1, alpha_opt = 1))
grz_pars = list(graze_start = "may",
                graze_days = 0,
                graze_hours = 0)
man_pars = list(conc_fresh = list(S2 = 0.01, urea = 2.4, SO4 = 0.2, TAN = 0.63,
                                  VFA = 2.83, xa_dead = 0, VSd = 73.8), pH = 6.88, dens = 1000)

grp_pars = list(grps = c('m0','m1','m2','sr1'),
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

chem_pars = list(COD_conv = c(CH4 = 0.2507, xa_dead = 0.73, VFA = 0.93, S = 0.5015, VS = 0.69, CO2_anaer = 0.53, CO2_aer = 1.1, CO2_sr = 1.2, CO2_ureo = 1.57,
                              C_xa_dead = 0.358,
                              C_VFA = 0.374, C_VSd = 0.344, C_N_urea = 0.429), 
                 kl = c(NH3 = 52, NH3_floor = 22, H2S = 0.02)) 

################### adjust for optimizer and exclude period 4 ###################

# Get emission data frame with average daily emission CH4 rate
days_sim <- max(meas$days[meas$period == 3], na.rm = T)
meas$day <- ceiling(meas$days)
dat <- as.data.frame(summarise(group_by(meas, day), CH4_emis_rate = mean(CH4_emis_rate, na.rm = TRUE), VFA_conc = mean(vfa / 1000 * 0.93, na.rm = TRUE)))
dat <- dat[dat$day <= ceiling(days_sim), ] 

################### run optimization calculation ###################

# Parameter estimation with L-BFGS-B 
new_pars <- data.frame(scale.qhat_opt = 0.313, scale.xa_fresh = 3.13, scale.alpha_opt = 2.47, scale.ks_coefficient = 1.17)
res <- 1 
loops <- 1 
pars.cal <- log10(new_pars)

for (i in 1:loops){
  
  # Weight so *total* weight for two variables (correcting for standard deviation in resCalc) are equal
  # Otherwise, change `1 *` to e.g., `3 *`.
  weights <- data.frame(CH4_emis_rate = 1 * (nc <- !is.na(dat$CH4_emis_rate)) / sum(nc), 
                        VFA_conc = 1 * (nv <- !is.na(dat$VFA_conc)) / sum(nv))
  
  # run calibration function
  cal <- optim(par = pars.cal, 
               fn = resCalc,
               dat = dat,
               weights = weights,
               to = c('CH4_emis_rate', 'VFA_conc'), 
               mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars, wthr_pars = wthr_pars, evap_pars = evap_pars,
               plot = F,
               method = 'L-BFGS-B',
               lower = log10(c(0.2, 1, 0.5, 0.5)), upper = log10(c(1, 50, 10, 20)),
               control = list(reltol = 0.01),
               hessian = TRUE
        )

  res1 <- cal$value
  res <- rbind(res, res1)
  new_pars1 <- 10^cal$par
  new_pars1 <- new_pars1[c("scale.qhat_opt", "scale.xa_fresh", "scale.alpha_opt", "scale.ks_coefficient")] # set to original order
  new_pars <- rbind(new_pars, new_pars1)
  pars.cal <- log10(new_pars[i + 1, sample(1:ncol(new_pars))]) # resample in preparation for next call. 
}

new_pars <- cbind(new_pars, res)
rmin <- which(new_pars$res == min(new_pars$res[-c(1)]))[1]
opt_pars <- new_pars[rmin,]
write.xlsx(opt_pars, '../data/opt_pars.xlsx') 
