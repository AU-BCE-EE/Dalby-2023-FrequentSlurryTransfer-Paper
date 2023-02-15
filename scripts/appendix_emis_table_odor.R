rm(list = ls())

# Load packages
library('dplyr')
library('tidyr')
library('readxl')
library('reshape2')

dat <- read.csv('../data/dat_stacked.csv')
dat_odor <- read_excel('../data/dat_comp.xlsx', sheet = "odor")


# Some cleaning
dat_odor$treatment <- tolower(dat_odor$treatment)
dat_odor <- rename(dat_odor, H2S = M35E, MT = M49E, TMA = M60E, HAC = M61E, 
                   HPA = M75E, HBA = M89E, HPENA = M103E, MP = M109E, SKAT = M132E)
dat_odor$period[dat_odor$period == 2] <- 4

# Drop obs with no periods
dat <- subset(dat, !is.na(period))

# Convert date/times to dates
dat$date.time <- dat$date
dat$date <- as.Date(dat$date)

dat_odor$date.time <- dat_odor$date
dat_odor$date <- as.Date(dat_odor$date)

# Add pig numbers to odor data by *period* (no change over time, apparently some pig count data missing)
pigs <- summarise(group_by(dat, period, treatment), pigs = mean(pigs, na.rm = TRUE))
dat_odor <- merge(dat_odor, pigs, by = c('period', 'treatment'))

dol <- melt(dat_odor, id.vars = c('date.time', 'date', 'treatment', 'period', 'pigs'),
            measure.vars = c('H2S','MT','TMA','HAC','HPA','HBA','HPENA','MP','SKAT', 'OAVE'))

# Calculate per pig per d only for H2S (not odor)
dol$value[dol$variable != 'OAVE'] <- (dol$value / dol$pigs)[dol$variable != 'OAVE']

odor_summ <- summarise(group_by(dol, variable, treatment, period), mn = mean(value), n = length(value))

# Summarise over full year (4 periods)
odor_summ_yr <- summarise(group_by(odor_summ, variable, treatment), amn = mean(mn), s = sd(mn), n = length(mn))
# Get full year emission reductions
odor_summ_yr <- mutate(group_by(odor_summ_yr, variable), amn.control = amn[treatment == 'control'], 
                       rel.red = round(100 * (amn.control - amn) / amn.control, 1))


write.csv(odor_summ, '../output/odor_summ.csv', row.names = FALSE)
write.csv(odor_summ_yr, '../output/odor_summ_year.csv', row.names = FALSE)
