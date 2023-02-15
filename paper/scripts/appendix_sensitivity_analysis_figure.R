rm(list = ls())

library('readxl')
library('gridExtra')
library('svglite')
library('ggplot2')
library('lubridate')
library('dplyr')
library('tidyr')

source('../functions/ggsave2x.R')

dat <- read_excel('../data/sensitivity.xlsx')

dat$mode <- factor(dat$mode)
new.lab <- as_labeller(c("1"= "alpha[opt]~0.02~and~q[max.opt]~reference", "2"= "alpha[opt]~0.1~and~q[max.opt]~reference", "3"= "alpha[opt]~0.02~and~q[max.opt]~0.4~x~reference", "4"= "alpha[opt]~0.1~and~q[max.opt]~0.4~x~reference"), label_parsed)

p <- ggplot(dat = dat, aes(x = perc_value, y = perc_CH4, col = as.factor(parm))) + geom_line() + geom_point() + 
  scale_color_discrete(labels=c(expression(alpha[opt]), expression(k[di]), expression(K[S.coef]), expression(a[enrich]), expression(q[max.opt]), expression(C[Xi.in]), expression(Y[i]))) +
  xlab("Parameter change (%)") + 
  ylab(bquote(CH[4]~emission~change~("%"))) + 
  facet_wrap(~mode, nrow = 2, ncol = 2, scales = "free", labeller = new.lab) + 
  theme(legend.title = element_blank())

svglite('../figures/fig_sensitivity.svg', width = 16/2.54, height = 10/2.54)
 p
dev.off() 

png('../figures/fig_sensitivity.png',  width = 16/2.54, height = 8/2.54, units = 'in', res = 600)
 p
dev.off()

pdf('../figures/fig_sensitivity.pdf',  width = 16/2.54, height = 8/2.54)
 p 
dev.off()
