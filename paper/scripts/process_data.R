rm(list = ls())

library("readxl")
library("lubridate")
library("openxlsx")

path_dat <- '../data/dat_comp.xlsx'
source('../ABM_functions/interpm.R')

# import data
dat <- read_excel(path_dat, sheet = "emis")
dat_ent <- read_excel(path_dat, sheet = "enteric")
dat_temp <- read_excel(path_dat, sheet = "temp")
dat_mass <- read_excel(path_dat, sheet = "mass")
dat_NH3 <- read_excel(path_dat, sheet = "emis NH3")
dat_analysis <- read_excel(path_dat, sheet = "analysis")
dat_weights <- read_excel(path_dat, sheet = "weights")
dat_vfa <- read_excel(path_dat, sheet = "vfa")
dat_odor <- read_excel(path_dat, sheet = "odor")

per <- read.csv('../data/periods.csv')

# add derived output to dat
dat$CH4_rate <- dat$CH4E * dat$pigs
dat$doy <- yday(dat$date)
dat$airflow <- dat$airflow * 24 # from m3/h to m3/day

C <- dat[dat$treatment == "control",]
FF <- dat[dat$treatment == "frequentflushing",]
ST <- dat[dat$treatment == "slurrytrays",]
SF <- dat[dat$treatment == "slurryfunnels",]

# interpolate enteric methane prod, temperature, slurry mass to periods and treatment data
dat_ent_C <- dat_ent[which(dat_ent$treatment == "control"),]
C$enteric <- as.data.frame(approx(y = dat_ent_C$CH4_enterisk * dat_ent_C$No_pigs, x = dat_ent_C$date, xout = C$date))[,2]
dat_temp_C <- dat_temp[which(dat_temp$treatment == "control"), c("date", "temp", "period", "treatment")]
C <- merge(C, dat_temp_C, all = T)
dat_mass_C <- dat_mass[which(dat_mass$treatment == "control"), c("date", "mass", "period", "treatment")]
C <- merge(C, dat_mass_C, all = T)
dat_NH3_C <- dat_NH3[which(dat_NH3$treatment == "control"),]
C$NH3 <- as.data.frame(approx(y = dat_NH3_C$NH3, x = dat_NH3_C$date, xout = C$date))[,2]
dat_weights_C <- dat_weights[which(dat_weights$treatment == "control"),]
C$weights <- as.data.frame(approx(y = dat_weights_C$weights, x = dat_weights_C$date, xout = C$date))[,2]
dat_vfa_C <- as.data.frame(dat_vfa[which(dat_vfa$treatment == 'control'),])
C <- merge(C, dat_vfa_C, all = T)
dat_H2S_C <- dat_odor[which(dat_odor$treatment == 'Control'),]
C$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_C$M35E, x = dat_H2S_C$date, xout = C$date))[,2]# g H2S /day

dat_ent_FF <- dat_ent[which(dat_ent$treatment == "frequentflushing"),]
FF$enteric <- as.data.frame(approx(y = dat_ent_FF$CH4_enterisk * dat_ent_FF$No_pigs, x = dat_ent_FF$date, xout = FF$date))[,2]
dat_temp_FF <- dat_temp[which(dat_temp$treatment == "frequentflushing"), c("date", "temp", "period", "treatment")]
FF <- merge(FF, dat_temp_FF, all = T)
dat_mass_FF <- dat_mass[which(dat_mass$treatment == "frequentflushing"), c("date", "mass", "period", "treatment")]
FF <- merge(FF, dat_mass_FF, all = T)
dat_NH3_FF <- dat_NH3[which(dat_NH3$treatment == "frequentflushing"),]
FF$NH3 <- as.data.frame(approx(y = dat_NH3_FF$NH3, x = dat_NH3_FF$date, xout = FF$date))[,2]
dat_weights_FF <- dat_weights[which(dat_weights$treatment == "frequentflushing"),]
FF$weights <- as.data.frame(approx(y = dat_weights_FF$weights, x = dat_weights_FF$date, xout = FF$date))[,2]
dat_vfa_FF <- as.data.frame(dat_vfa[which(dat_vfa$treatment == 'frequentflushing'),])
FF <- merge(FF, dat_vfa_FF, all = T)
dat_H2S_FF <- dat_odor[which(dat_odor$treatment == 'Frequentflushing'),]
FF$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_FF$M35E, x = dat_H2S_FF$date, xout = FF$date))[,2] # g H2S /day

dat_ent_SF <- dat_ent[which(dat_ent$treatment == "slurryfunnels"),]
SF$enteric <- as.data.frame(approx(y = dat_ent_SF$CH4_enterisk * dat_ent_SF$No_pigs, x = dat_ent_SF$date, xout = SF$date))[,2]
dat_temp_SF <- dat_temp[which(dat_temp$treatment == "frequentflushing"), c("date", "temp", "period", "treatment")]
SF <- merge(SF, dat_temp_SF, all = T)
dat_mass_SF <- dat_mass[which(dat_mass$treatment == "slurryfunnels"), c("date", "mass", "period", "treatment")]
SF <- merge(SF, dat_mass_SF, all = T)
dat_NH3_SF <- dat_NH3[which(dat_NH3$treatment == "slurryfunnels"),]
SF$NH3 <- as.data.frame(approx(y = dat_NH3_SF$NH3, x = dat_NH3_SF$date, xout = SF$date))[,2]
dat_weights_SF <- dat_weights[which(dat_weights$treatment == "slurryfunnels"),]
SF$weights <- as.data.frame(approx(y = dat_weights_SF$weights, x = dat_weights_SF$date, xout = SF$date))[,2]
dat_vfa_SF <- as.data.frame(dat_vfa[which(dat_vfa$treatment == 'slurryfunnels'),])
SF <- merge(SF, dat_vfa_SF, all = T)
dat_H2S_SF <- dat_odor[which(dat_odor$treatment == 'Slurryfunnels'),]
SF$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_SF$M35E, x = dat_H2S_SF$date, xout = SF$date))[,2]# g H2S /day

dat_ent_ST <- dat_ent[which(dat_ent$treatment == "slurrytrays"),]
ST$enteric <- as.data.frame(approx(y = dat_ent_ST$CH4_enterisk * dat_ent_ST$No_pigs, x = dat_ent_ST$date, xout = ST$date))[,2]
dat_temp_ST <- dat_temp[which(dat_temp$treatment == "frequentflushing"), c("date", "temp", "period", "treatment")]
ST <- merge(ST, dat_temp_ST, all = T)
dat_mass_ST <- dat_mass[which(dat_mass$treatment == "slurrytrays"), c("date", "mass", "period", "treatment")]
ST <- merge(ST, dat_mass_ST, all = T)
dat_NH3_ST <- dat_NH3[which(dat_NH3$treatment == "slurrytrays"),]
ST$NH3 <- as.data.frame(approx(y = dat_NH3_ST$NH3, x = dat_NH3_ST$date, xout = ST$date))[,2]
dat_weights_ST <- dat_weights[which(dat_weights$treatment == "slurrytrays"),]
ST$weights <- as.data.frame(approx(y = dat_weights_ST$weights, x = dat_weights_ST$date, xout = ST$date))[,2]
dat_vfa_ST <- as.data.frame(dat_vfa[which(dat_vfa$treatment == 'slurrytrays'),])
ST <- merge(ST, dat_vfa_ST, all = T)
dat_H2S_ST <- dat_odor[which(dat_odor$treatment == 'Slurrytrays'),]
ST$H2S_emis_rate <- as.data.frame(approx(y = dat_H2S_ST$M35E, x = dat_H2S_ST$date, xout = ST$date))[,2] # g H2S /day


C$CO2_enteric <- (0.136 * C$weights^0.573 ) * 1000
FF$CO2_enteric <- (0.136 * FF$weights^0.573 ) * 1000
ST$CO2_enteric <- (0.136 * ST$weights^0.573 ) * 1000
SF$CO2_enteric <- (0.136 * SF$weights^0.573 ) * 1000

# correct for enteric CH4, final units in g slurry-CH4 per section per day 
C$CH4_emis_rate <- C$CH4_rate - C$enteric
FF$CH4_emis_rate <- FF$CH4_rate - FF$enteric
SF$CH4_emis_rate <- SF$CH4_rate - SF$enteric
ST$CH4_emis_rate <- ST$CH4_rate - ST$enteric

#rename NH3
C$NH3_emis_rate <- C$NH3
FF$NH3_emis_rate <- FF$NH3
SF$NH3_emis_rate <- SF$NH3
ST$NH3_emis_rate <- ST$NH3

# correct for enteric CO2, units of CO2 per pig per day
C$CO2_emis_rate <- C$CO2E - C$CO2_enteric
FF$CO2_emis_rate <- FF$CO2E - FF$CO2_enteric
SF$CO2_emis_rate <- SF$CO2E - SF$CO2_enteric
ST$CO2_emis_rate <- ST$CO2E - ST$CO2_enteric

# Add *definitive* period info (overwrites existing period info)
per$date.time.in <- parse_date_time(paste(per$date.in, '12:00'), orders = 'dmyHM')
per$date.time.out <- parse_date_time(paste(per$date.out, '12:00'), orders = 'dmyHM')

# Pull out "day 0" for *everything*
day0 <- min(per$date.time.in)

C$days <- as.numeric(difftime(C$date, day0, units = "days"))
FF$days <- as.numeric(difftime(FF$date, day0, units = "days"))
SF$days <- as.numeric(difftime(SF$date, day0, units = "days"))
ST$days <- as.numeric(difftime(ST$date, day0, units = "days"))

# add start concentration of VFA, mg/kg slurry
C$vfa[1] <- 2830
FF$vfa[1] <- 2830
SF$vfa[1] <- 2830
ST$vfa[1] <- 2830

dat_simple <- list(C = C, FF = FF, SF = SF, ST = ST)

dat_simple$period <- NA

for (i in 1:nrow(per)) {
  dat_simple$period[dat_simple$date > per$date.time.in[i] & dat_simple$date < per$date.time.out[i]] <- per$period[i]
}

# Export data
write.xlsx(dat_simple, "../data/dat_simple.xlsx")

# Stack data for plots
dat_stacked <- rbind(C, FF, SF, ST)

dat_stacked$period <- NA

for (i in 1:nrow(per)) {
  dat_stacked$period[dat_stacked$date > per$date.time.in[i] & dat_stacked$date < per$date.time.out[i]] <- per$period[i]
}

dat_stacked$days <- as.numeric(difftime(dat_stacked$date, day0, units = 'days'))

write.csv(dat_stacked, '../data/dat_stacked.csv', row.names = FALSE)

