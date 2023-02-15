# Runs stats, saves pdf

library(rmarkdown)
library(FSA)
render('../scripts/stats.Rmd', output_dir = '../output')
render('../scripts/stats_select.Rmd', output_dir = '../output')

