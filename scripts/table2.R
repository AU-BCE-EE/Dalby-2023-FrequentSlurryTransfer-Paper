rm(list = ls())

source('../scripts/process_data.R')

library('dplyr')
library('tidyr')

comp_summary <- summarise(group_by(dat_analysis, period, treatment), DM = mean(DM, na.rm = TRUE), VS = mean(VS, na.rm = TRUE), 
                      pH = mean(pH, na.rm = TRUE), TN = mean(TN, na.rm = TRUE), TAN = mean(TAN, na.rm = TRUE))

comp_summary_year <- summarise(group_by(comp_summary, treatment), mean_DM = mean(DM, na.rm = TRUE), mean_VS = mean(VS, na.rm = TRUE), 
                          mean_pH = mean(pH, na.rm = TRUE), mean_TN = mean(TN, na.rm = TRUE), mean_TAN = mean(TAN, na.rm = TRUE),
                          sd_DM = sd(DM, na.rm = TRUE), sd_VS = sd(VS, na.rm = TRUE), 
                          sd_pH = sd(pH, na.rm = TRUE), sd_TN = sd(TN, na.rm = TRUE), sd_TAN = sd(TAN, na.rm = TRUE))

vfa_dat_sum <- summarise(group_by(dat_vfa, period, treatment), vfa = mean(vfa/1000, na.rm = TRUE))
vfa_summary <- summarise(group_by(vfa_dat_sum, treatment), mean = mean(vfa, na.rm = TRUE), std = sd(vfa, na.rm = TRUE))


#temperatures

temp_summary <- summarise(group_by(dat_stacked, period, treatment), temp_air = mean(temp_air, na.rm = TRUE), temp = mean(temp, na.rm = TRUE))

temp_summary_year <- summarise(group_by(temp_summary, treatment), mean_temp_air = mean(temp_air, na.rm = TRUE), mean_temp = mean(temp, na.rm = TRUE),
                          std_temp_air = sd(temp_air, na.rm = TRUE), std_temp = sd(temp, na.rm = TRUE))

write.csv(comp_summary, '../output/comp_summary.csv', row.names = FALSE)
write.csv(comp_summary_year, '../output/comp_summary_year.csv', row.names = FALSE)
write.csv(temp_summary, '../output/temp_summary.csv', row.names = FALSE)
write.csv(temp_summary_year, '../output/temp_summary_year.csv', row.names = FALSE)


# calculation of slurry retention times

slurry_retention_dat = NULL

for (w in 1:4){

treat_dat = NULL

  for (o in c('control', 'frequentflushing', 'slurryfunnels', 'slurrytrays')){

test_dat <- dat_mass[which(dat_mass$period == w & dat_mass$treatment == o),]

time = NULL
mass_age = NULL
ave_time = NULL

for (i in 2:nrow(test_dat)){
  if(test_dat$mass[i] > test_dat$mass[i-1]){
    mass_age1 <- test_dat$mass[i] - test_dat$mass[i-1]
    time1 <- (test_dat$date[i] - test_dat$date[i-1])
    if(i == 2 | i == 3){
      ave_time1 <- time1/2
    } else {
    ave_time1 <- (time1 * mass_age1/(mass_age[i-2] + mass_age1) + (time[i-2] + time1) * mass_age[i-2] / (mass_age[i-2] + mass_age1))/2 
    }
  } else if (test_dat$mass[i] < test_dat$mass[i-1]) {
    mass_age1 <- test_dat$mass[i]
    time1 <- test_dat$date[i] - test_dat$date[i-2]
    if(i == 2 | i == 3){
      ave_time1 <- time1/2
    } else {
    ave_time1 <- time1 - time[i-2] + ave_time[i-2]
    }
  }
  time <- rbind(time, time1)
  mass_age <- rbind(mass_age, mass_age1)
  ave_time <- rbind(ave_time, ave_time1)
  mass_weighed <- mass_age * ave_time
  ret_time <- rbind(sum(mass_weighed)/sum(mass_age))
}

treat_dat <- rbind(treat_dat, ret_time)

  }

slurry_retention_dat <- cbind(slurry_retention_dat, treat_dat)

}

colnames(slurry_retention_dat) <- c(1,2,3,4)
rownames(slurry_retention_dat) <- c('control', 'frequentflushing', 'slurryfunnels', 'slurrytrays')

write.csv(slurry_retention_dat, '../output/slurry_retention.csv', row.names = FALSE)

C.mean <- mean(slurry_retention_dat[1,])
WF.mean <- mean(slurry_retention_dat[2,])
SF.mean <- mean(slurry_retention_dat[3,])
ST.mean <- mean(slurry_retention_dat[4,])

C.sd <- sd(slurry_retention_dat[1,])
WF.sd <- sd(slurry_retention_dat[2,])
SF.sd <- sd(slurry_retention_dat[3,])
ST.sd <- sd(slurry_retention_dat[4,])


# Slurry mass average in pits

mass_summary <- dat_stacked %>% select(mass, period, treatment) %>% filter(!is.na(mass), !is.na(period)) %>% 
  group_by(treatment, period) %>% summarise(mass.mean = mean(mass, na.rm =T)) %>% group_by(treatment) %>% 
  summarise(mass_sd = sd(mass.mean, na.rm = TRUE), mass.mean = mean(mass.mean, na.rm = T))

write.csv(mass_summary, '../output/slurry_mass_summary.csv', row.names = FALSE)
