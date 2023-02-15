# Main script for running all the code related to the paper:
# Dalby et al., Simple management changes drastically reduce pig house methane emission 
# in combined experimental and modelling study. Environ. Sci. Technol. 2023, XXXX, XXX, 
# DOI: 10.1021/acs.est.2c08891

### Code for manuscript ###

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
source('../scripts/predictions.R')
# Figure 5
source('../scripts/fig_predict.R')

### Code for Supporting information ###

# Figure S1
source('../scripts/appendix_qhat_dynamics.R')
# Table S5
source('../scripts/appendix_emis_table_odor.R')
# Run sensitivity analysis
source('../scripts/appendix_sensitivity_analysis.R')
# Figure S2
source('../scripts/appendix_sensitivity_analysis_figure.R')

