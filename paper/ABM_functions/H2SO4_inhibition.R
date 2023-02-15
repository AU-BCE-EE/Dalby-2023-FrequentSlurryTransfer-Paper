H2SO4_inhibition <- function(SO4_conc){
  td <- data.frame(SO4_conc = c(0, 0.392251223, 0.686439641, 1.046003263, 1.470942088, 1.961256117, 4), 
                   H2SO4_inhib = c(1, 0.54, 0.32, 0.19, 0.11, 0.04, 0.001))
  H2SO4_inhib <- approx(td$SO4_conc, td$H2SO4_inhib, xout = SO4_conc, rule = 2)$y   
}