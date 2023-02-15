# Creates model sensitivity figure

rm(list = ls())

library('readxl')
library('gridExtra')
library('svglite')
library('ggplot2')
library('dplyr')
library('tidyr')
library('data.table')

source('../functions/ggsave2x.R')
source('../functions/dfcombos.R')

# Load data
ref_emis <- fread('../data/ref_emis.csv')
dat_stacked <- read.csv('../data/dat_stacked.csv')

#slurry produced in barn
sum_barn_slurry <- dat_stacked %>% select(mass, period, treatment) %>% filter(!is.na(mass), period != 4, !is.na(period)) %>% 
  mutate(diff = c(0,diff(mass))) %>% group_by(treatment, period) %>% filter(diff > 0) %>%  summarise(cum_slurry = sum(diff))

#slurry flushed to storage
summary_barn_parms <- data.frame(read.csv('../data/summary_barn_parms.csv'))

# extra parms for unit conversion in storage
summary_storage_parms <- data.frame(read.csv('../data/summary_storageVar_parms.csv'))
mean_storage_slurry <- summary_storage_parms[summary_storage_parms$parm == 'slurry_mass_mean', ]

#emissions from barn and storage
summary_barn <- data.frame(read.csv('../data/summary_barn.csv'))
summary_barn <- summary_barn[summary_barn$alpha_opt == 1 & summary_barn$qhat_opt == 1,][1,]
names <- colnames(summary_barn)
summary_barn <- data.frame(matrix(summary_barn, 36, 6, byrow = T))
colnames(summary_barn) <- names
summary_barn <- as.data.frame(apply(summary_barn, 2, as.numeric))

summary_storage <- data.frame(read.csv('../data/summary_storageVar.csv'))

meas_C <- as.data.frame(read_excel('../data/dat_simple.xlsx', sheet = "C"))

# convert emission from averaged g CH4/day to kg CH4/year/m3 slurry excreted from animal
sum_barn_slurry <- data.frame(summarise(group_by(sum_barn_slurry, treatment), slurry_cum = sum(cum_slurry/1000, na.rm = T)))
summary_barn$C <- summary_barn$C/1000 * 365/ (sum_barn_slurry$slurry_cum[sum_barn_slurry$treatment == 'control'] * 365/ summary_barn_parms$scores[summary_barn_parms$treat == 'C' & summary_barn_parms$parm == 'time'])
summary_barn$FF <- summary_barn$FF/1000 * 365/ (sum_barn_slurry$slurry_cum[sum_barn_slurry$treatment == 'frequentflushing'] * 365 / summary_barn_parms$scores[summary_barn_parms$treat == 'FF' & summary_barn_parms$parm == 'time'])
summary_barn$SF <- summary_barn$SF/1000 * 365/ (sum_barn_slurry$slurry_cum[sum_barn_slurry$treatment == 'slurryfunnels'] * 365 / summary_barn_parms$scores[summary_barn_parms$treat == 'SF' & summary_barn_parms$parm == 'time'])
summary_barn$ST <- summary_barn$ST/1000 * 365/ (sum_barn_slurry$slurry_cum[sum_barn_slurry$treatment == 'slurrytrays'] * 365 / summary_barn_parms$scores[summary_barn_parms$treat == 'ST' & summary_barn_parms$parm == 'time'])

# kg CH4 / year / ton slurry ab barn
# convert emission from averaged g CH4/day to kg CH4/year/m3 slurry excreted ab barn
sum_barn_eff_slurry <- summary_barn_parms[summary_barn_parms$parm == 'slurry_mass_eff',] %>% mutate(scores = scores / 1000)
summary_storage$C <- summary_storage$C/1000 * 365 / (sum_barn_eff_slurry$scores[sum_barn_eff_slurry$treat == 'C'] * 365/ summary_barn_parms$scores[summary_barn_parms$treat == 'C' & summary_barn_parms$parm == 'time']) # kg CH4 / year / ton slurry ab barn
summary_storage$FF <- summary_storage$FF/1000 * 365 / (sum_barn_eff_slurry$scores[sum_barn_eff_slurry$treat == 'FF'] * 365 / summary_barn_parms$scores[summary_barn_parms$treat == 'FF' & summary_barn_parms$parm == 'time']) # kg CH4 / year / ton slurry ab barn
summary_storage$SF <- summary_storage$SF/1000 * 365 / (sum_barn_eff_slurry$scores[sum_barn_eff_slurry$treat == 'SF'] * 365 / summary_barn_parms$scores[summary_barn_parms$treat == 'SF' & summary_barn_parms$parm == 'time']) # kg CH4 / year / ton slurry ab barn
summary_storage$ST <- summary_storage$ST/1000 * 365 / (sum_barn_eff_slurry$scores[sum_barn_eff_slurry$treat == 'ST'] * 365 / summary_barn_parms$scores[summary_barn_parms$treat == 'ST' & summary_barn_parms$parm == 'time']) # kg CH4 / year / ton slurry ab barn

# kg CH4 / year / ton slurry present in storage
#mean_storage_slurry <- mean_storage_slurry %>% mutate(scores = scores/1000)
#summary_storage$C <- summary_storage$C/1000 * 365 / mean_storage_slurry$scores[mean_storage_slurry$treat == 'C']
#summary_storage$FF <- summary_storage$FF/1000 * 365 / mean_storage_slurry$scores[mean_storage_slurry$treat == 'FF']
#summary_storage$SF <- summary_storage$SF/1000 * 365 / mean_storage_slurry$scores[mean_storage_slurry$treat == 'SF']
#summary_storage$ST <- summary_storage$ST/1000 * 365 / mean_storage_slurry$scores[mean_storage_slurry$treat == 'ST']


summary_net <- cbind(summary_barn[, !names(summary_barn) %in% c('alpha_opt','qhat_opt')] + summary_storage[, !names(summary_storage) %in% c('alpha_opt','qhat_opt')], summary_barn[, c('alpha_opt', 'qhat_opt')])


summary_dat <- rbind(cbind(summary_barn, source = 'Barn'), cbind(summary_storage, source = 'Storage'), cbind(summary_net, source = 'Net'))
summary_dat$reduction_FF <- (1 - summary_dat$FF/summary_dat$C) * 100
summary_dat$reduction_SF <- (1 - summary_dat$SF/summary_dat$C) * 100
summary_dat$reduction_ST <- (1 - summary_dat$ST/summary_dat$C) * 100

# Clean and reshape

ref_emis$animal <- factor(ref_emis$animal, levels = c('Finisher', 'Gestating sow', 'Farrowing sow', 'Weaner'))

summary_dat_long <- summary_dat %>% pivot_longer(c(reduction_FF, reduction_SF, reduction_ST),
                                         names_to = 'treat', values_to = 'reduction') 

summary_dat_long$source <- factor(summary_dat_long$source, levels = c("Barn", "Storage", "Net")) 

new.lab <- as_labeller(c(reduction_FF = "W.~flush.~(WF)", reduction_SF = "S.~funnels.~(SF)", reduction_ST= "Slurry~trays~(ST)", 
                          Barn = "Barn", Storage = "Storage", Net = "Net"), label_parsed)


# More data processing for plot
summary_dat_long <- data.table(summary_dat_long)
summary_dat_long$varied <- NA
summary_dat_long$varied[summary_dat_long$alpha_opt == 1] <- 'a qhat'
summary_dat_long$varied[summary_dat_long$qhat_opt == 1] <- 'b alpha'
summary_dat_long$varied[summary_dat_long$qhat_opt != 1 & summary_dat_long$alpha_opt != 1] <- 'c both'
summary_dat_long$x <- summary_dat_long$qhat_opt
summary_dat_long$x[summary_dat_long$qhat_opt == 1] <- summary_dat_long$alpha_opt[summary_dat_long$qhat_opt == 1]
summary_dat_long$x <- log10(summary_dat_long$x)
ref1 <- summary_dat_long[summary_dat_long$qhat_opt == 1 & summary_dat_long$alpha_opt == 1, ] 
ref1$varied <- NULL
ref1 <- dfcombos(ref1, data.frame(varied = unique(summary_dat_long$varied)))
summary_dat_long <- rbind(summary_dat_long, ref1)
summary_dat_long <- summary_dat_long[order(summary_dat_long$varied, summary_dat_long$qhat_opt * summary_dat_long$alpha_opt), ]

ref <- unique(summary_dat_long[, c('C', 'alpha_opt', 'qhat_opt', 'varied', 'x', 'source')])
ref$treat <- 'Control (C)'


# Plots
cols <- scales::hue_pal()(3)
dd <- subset(ref, alpha_opt == 1 & qhat_opt == 1)
refplot <- ggplot(ref, aes(x, C, colour = varied)) + 
             geom_point(data = dd, colour = 'gray55', size = 2) + 
             geom_path() + 
             theme_bw() +
             facet_grid(source~factor(treat)) +
             labs(x = expression('Log'[10]~'par. adj.'), y = expression('Emission rate'~(kg~CH[4]~m^'3'~yr^'-1')), tag = 'b.') +
             scale_colour_manual(name = 'Parameter', values = cols, 
                                 labels = expression(q['max, opt'], alpha['opt'], 'Both')) +
             theme(legend.position = 'none')
refplot
dd <- subset(summary_dat_long, alpha_opt == 1 & qhat_opt == 1)
redplot <- ggplot(summary_dat_long, aes(x, reduction, colour = varied)) + 
             geom_point(data = dd, colour = 'gray55', size = 2) + 
             geom_path() + 
             theme_bw() +
             facet_grid(source~treat, labeller = new.lab, scales = 'free_y') +
             labs(x = expression('Log'[10]~'parameter adjustment'), y = 'Relative emission reduction (%)', tag = 'c.') +
             theme(legend.position = 'right') +
             scale_colour_manual(name = 'Parameter', values = cols, 
                                 labels = expression(q['max, opt'], alpha['opt'], 'Both'))
redplot
refemisplot <- ggplot(ref_emis, aes(x = animal, y = CH4)) + geom_boxplot(varwidth = TRUE) + geom_jitter(width = 0.2) + theme_bw() + 
  labs(x ="", y = expression('Methane emission rate'~(g~pig^'-1'~d^'-1')), colour = "", tag = 'a.') + 
  stat_summary(fun.y ="mean", color="red")

p_all <- grid.arrange(refemisplot, refplot, redplot, widths = c(1.5,4), heights = c(2.5,4), layout_matrix = rbind(c(1), c(2,3)))


svglite("../figures/fig_predict.svg", width = 18/2.54, height = 18/2.54)
grid::grid.draw(p_all)
dev.off()

png('../figures/fig_predict.png',  width = 18/2.54, height = 18/2.54, units = 'in', res = 600)
grid::grid.draw(p_all)
dev.off()

pdf('../figures/fig_predict.pdf',  width = 18/2.54, height = 18/2.54)
grid::grid.draw(p_all)
dev.off()

# summarise results
sum_reference <- filter(summary_dat_long, alpha_opt == 1, qhat_opt == 1)
summarise(group_by(summary_dat_long, treat, source), mean(C, na.rm = T), min(reduction, na.rm = T), max(reduction, na.rm = T), mean(reduction, na.rm = T), median(reduction, na.rm = T))

# get results from reference situation and emission ratio
barn_reference <- data.table(summary_barn)[alpha_opt == 1.0 & qhat_opt == 1.0, C]   
storage_reference <- data.table(summary_storage)[alpha_opt == 1.0 & qhat_opt == 1.0, C]   
barn_reference/(storage_reference + barn_reference)

# check reasons for reductions in storage with frequent flushing techniques
summary_storage_parms <- data.frame(read.csv('../data/summary_storage_parms.csv'))
ggplot(summary_storage_parms, aes(x = treat, y = scores)) + geom_point() + facet_wrap(~parm, scales = "free")





