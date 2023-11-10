# Create emission rate figure

rm(list = ls())

# Packages
library(ggplot2)
library(data.table)
library(lubridate)
library(reshape2)

# Other functions
source('../functions/rbindf.R')
source('../functions/dfcombos.R')
source('../functions/ggsave2x.R')

# Language for date labels
Sys.setlocale("LC_ALL", "English")

# Load data
dat <- fread('../data/dat_stacked.csv')
per <- read.csv('../data/periods.csv')

# Some cleaning
dat$date.time <- ymd_hms(dat$date) 
dat$doy <- as.integer(as.character(dat$date.time, format = '%j'))
dat$date <- as.Date(dat$date.time)

per$date.time.in <- dmy_hm(paste(per$date.in, '12:00')) 
per$date.time.out <- dmy_hm(paste(per$date.out, '12:00')) 
per$date.in <- dmy(per$date.in)
per$date.out <- dmy(per$date.out)

# Pull out separate temeprature data
tdat <- subset(dat, !is.na(temp) & (treatment == 'control' | treatment == 'frequentflushing'))

dat <- subset(dat, !is.na(period))

# Normalize emission to 1 pig (already done for CO2)
dat$NH3_emis_rate <- dat$NH3_emis_rate  / dat$pigs
dat$H2S_emis_rate <- dat$H2S_emis_rate  / dat$pigs
dat$CH4_emis_rate <- dat$CH4_emis_rate  / dat$pigs

# Interpolate missing slurry mass
source('interpm.R')
dat <- interpm(dat, 'date', 'mass', by = 'treatment')
dat$CH4_emis_rate_norm <- dat$CH4_emis_rate * dat$pigs / dat$mass
# Reshape
dl <- as.data.table(melt(dat, id.vars = c('treatment', 'date', 'date.time', 'doy', 'period', 'pigs'),
           measure.vars = c('CH4_emis_rate', 'CH4_emis_rate_norm', 'CO2_emis_rate', 'NH3_emis_rate', 'H2S_emis_rate')))

dl$variable <- factor(dl$variable, levels = c('CH4_emis_rate', 'CH4_emis_rate_norm', 'CO2_emis_rate', 'NH3_emis_rate', 'H2S_emis_rate'),
                      labels = c('Methane', 'Methane normalized', 'Carbon dioxide', 'Ammonia', 'Hydrogen sulfide'))
dl$treatment <- factor(dl$treatment, levels = c('control', 'frequentflushing', 'slurryfunnels', 'slurrytrays'),
                       labels = c('Control (C)', 'Weekly\nflushing (WF)\n', 'Slurry\nfunnels (SF)\n', 'Slurry\ntrays (ST)\n'))

# Daily average emission rate
dl$date.group <- as.integer(dl$date - min(dl$date)) %/% 5
# Daily for all except . . .
dl1 <- dl[variable != 'Hydrogen sulfide', .(emis.ave = mean(value)), by = .(treatment, period, date, pigs, variable)]
# . . . use multi-day average for H2S
dl2 <- dl[variable == 'Hydrogen sulfide', .(emis.ave = mean(value), date = mean(date)), 
          by = .(treatment, period, date.group, pigs, variable)]
dl <- rbindf(dl1, dl2)

# Emission rate per pig
dl$emis.pig.d <- dl$emis.ave / dl$pigs

dl$value[is.na(dl$period)] <- NA

# Drop H2S in missing periods
dl <- subset(dl, variable != 'Hydrogen sulfide' | !period %in% 2:3)

# Drop NAs
dl <- subset(dl, !is.na(emis.ave) & !is.na(treatment))
# Add NAs for start and end of periods to break lines in plots
dl <- rbindf(dl, expand.grid(date = c(per$date.in, per$date.out), 
                             treatment = na.omit(unique(dl$treatment)), 
                             variable = unique(dl$variable), 
                             value = NA))

dl <- dl[order(dl$variable, dl$date), ]
x <- subset(dl, variable == 'Methane normalized' & 
            date < as.POSIXct('2020 11 01', format = '%Y %m %d') &
            date > as.POSIXct('2020 08 13', format = '%Y %m %d'))
y <- subset(x, date == max(date))

ggplot(x, aes(date, emis.ave, colour = treatment)) +
  geom_line() +
  geom_label(data = y, aes(label = treatment)) +
  labs(x = 'Date', y = expression('Normalized'~CH[4]~'emission rate'~(g~kg^'-1'~d^'-1'))) +
  theme(legend.position = 'none')
 
ep <- ggplot(dl, aes(date, emis.ave, colour = treatment)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.7, size = 0.3) +
  facet_grid(variable ~ ., scales = 'free') +
  labs(x = 'Date', y = expression('Emission rate'~(g~pig^'-1'~d^'-1')), colour = '') +
  xlim(as.Date(c('2020-05-30', '2021-04-30'))) +
  geom_vline(xintercept = c(per$date.in, per$date.out), lty = 2, col = 'gray65') +
  theme_bw() +
  theme(legend.position = 'top',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

# Sort out temperature data
ta <- tdat[, .(temp.ave = mean(temp, na.rm = TRUE)), by = .(date, period, treatment)]
ta <- rbindf(ta, data.frame(date = c(per$date.in, per$date.out), period = NA, temp.ave = NA))
ta$temp.ave[is.na(ta$period)] <- NA
ta$variable <- 'Temperature'
ta <- subset(ta, !is.na(temp.ave) | is.na(period))
ta$treatment <- factor(ta$treatment, levels = c('control', 'frequentflushing', 'slurryfunnels', 'slurrytrays'),
                       labels = c('Control (C)', 'Weekly\nflushing (WF)\n', 'Slurry\nfunnels (SF)\n', 'Slurry\ntrays (ST)\n'))
ta <- ta[!is.na(ta$treatment), ]


# Make some date labels
dd <- expand.grid(2020:2021, 1:12, 1)
ddl <- as.Date(paste(dd[, 1], dd[, 2], dd[, 3], sep = '-'))
cc <- scales::hue_pal()(4)[1:2]  
tp <- ggplot(ta, aes(date, temp.ave, colour = treatment)) +
  geom_line() + 
  geom_point(alpha = 0.7, size = 0.3) +
  labs(x = 'Date', y = expression('Slurry temp.'~('\u00b0'*C)), colour = '') +
  facet_grid(variable ~ ., scales = 'free') +
  geom_vline(xintercept = c(per$date.in, per$date.out), lty = 2, col = 'gray65') +
  scale_x_date(breaks = ddl, labels = format(ddl, '%b'), limits = as.Date(c('2020-05-30', '2021-04-30'))) +
  scale_colour_manual(values = cc) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = 'none') 

# NTS: #@! stupid ggplotGrob() makes some kind of table and the number of columns need to match to use
# NTS: rbind() as below
# NTS: the theme(legend.position =) is needed in tp to  get 12 x 10 instead of the 12 x 12 we get without it
# NTS: I cannot find any description of the 12 x 10 etc. and str() gives a massive object with 
# NTS: without clear info on what the 12 x 10 or 12 x 12 is 

png('../figures/fig_emis_temp.png', height = 7, width = 6.5, units = 'in', res = 600)
  grid::grid.draw(rbind(ggplotGrob(ep), ggplotGrob(tp)))
dev.off()

pdf('../figures/fig_emis_temp.pdf', height = 7, width = 6.5)
  grid::grid.draw(rbind(ggplotGrob(ep), ggplotGrob(tp)))
dev.off()

