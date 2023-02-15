abm_packages <- function(){
  x<-c('Rcpp', 'readxl', 'dplyr', 'rioja',
       'zoo', 'lubridate', 'openxlsx', 'deSolve', 
       'tidyr', 'textclean', 'data.table')
  lapply(x, require, character.only = TRUE)
}

