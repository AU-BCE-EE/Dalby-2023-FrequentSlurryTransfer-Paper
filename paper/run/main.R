# Main script for running all the code related to the paper
 
# Get raw data from excel file and preprocess
source('../scripts/process_data.R')
# Figure 2
source('../scripts/fig_emis.R')
# Figure 3
source('../scripts/fig_nslurry.R')
# Table 1
source('../scripts/table1.R')
# statistics Table 1
source('../scripts/stats.R')
# Table 2
source('../scripts/table2.R')
# Run optimization of code (Takes several hours to execute)
source('../scripts/optimization.R')
# Figure 4 and Table 4
source('../scripts/fig_valid.R')
# Run prediction scenarios (Takes several hours to execute)
source()