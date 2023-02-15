# Creates means etc. for table 1 in the paper

rm(list = ls())

# Load packages
library('dplyr')
library('tidyr')
library('readxl')
library('reshape2')

# Load data
dat <- read.csv('../data/dat_stacked.csv')
dat_odor <- read_excel('../data/dat_comp.xlsx', sheet = "odor")

# Some cleaning
dat_odor$treatment <- tolower(dat_odor$treatment)
dat_odor <- rename(dat_odor, H2S = M35E)

# Period two in raw dat of odor is actually period 4 in other datasets. Correct here
dat_odor$period[dat_odor$period == 2] <- 4

# Drop obs with no periods
dat <- subset(dat, !is.na(period))

# Convert date/times to dates
dat$date.time <- dat$date
dat$date <- as.Date(dat$date)

dat_odor$date.time <- dat_odor$date
dat_odor$date <- as.Date(dat_odor$date)

# Add pig numbers to odor data by date now
# Note that some odor obs are dropped--dates must be missing in pigs data (would be better to interpolate pig numbers and not lose odor measurements)
pigs <- summarise(group_by(dat, period, treatment, date), pigs = mean(pigs, na.rm = TRUE))
dat_odor <- merge(dat_odor, pigs, by = c('period', 'treatment', 'date'))

# Reshape for summary
dl <- melt(dat, id.vars = c('date.time', 'date', 'treatment', 'period', 'pigs'),
           measure.vars = c('CH4_rate', 'CH4_emis_rate', 'NH3_emis_rate', 'CO2E', 'CO2_emis_rate'))

# Normalize to "per pig"
# Apparently CO2 is already per pig
dl$value[!grepl('CO2', dl$variable)] <- (dl$value / dl$pigs)[!grepl('CO2', dl$variable)]

# Reshape odor
dol <- melt(dat_odor, id.vars = c('date.time', 'date', 'treatment', 'period', 'pigs'),
            measure.vars = c('H2S', 'OAVE'))

# Calculate per pig per d only for H2S (not odor)
dol$value[dol$variable == 'H2S'] <- (dol$value / dol$pigs)[dol$variable == 'H2S']

# Combine odor data with other emission data and drop missing values
dl <- rbind(dl, dol)
dl <- subset(dl, !is.na(value))

# Summarise by period
emis_summ <- summarise(group_by(dl, variable, treatment, period), mn = mean(value), n = length(value))

# Summarise over full year (4 periods)
emis_summ_yr <- summarise(group_by(emis_summ, variable, treatment), amn = mean(mn), s = sd(mn), n = length(mn))
# Get full year emission reductions
emis_summ_yr <- mutate(group_by(emis_summ_yr, variable), amn.control = amn[treatment == 'control'], 
                       rel.red = round(100 * (amn.control - amn) / amn.control, 1))
                                      
# Export
write.csv(emis_summ, '../output/emis_summary.csv', row.names = FALSE)
write.csv(emis_summ_yr, '../output/emis_summary_year.csv', row.names = FALSE)
