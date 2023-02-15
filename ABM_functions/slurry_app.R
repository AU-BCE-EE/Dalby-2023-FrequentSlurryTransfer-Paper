slurry_app <- function(days, begin, specs, from, slurry){ 
  
  years <- floor(days/365)
  
    if(from == 'storage'){
      a <- 's'
    }
    
    if(from == 'digestate'){
      a <- 'd'
    }
  
  mon <- specs[, grepl(paste('^app_t[0-9]', a, sep=""), names(specs))]
  apps <- doy(month = mon, begin = begin)$day
  rems <- as.numeric(specs[, grepl(paste('^app[0-9]', a, sep=""), names(specs))])
  
  t_add <- rep(seq(0, days - 365, 365), each = length(apps))
  appsl <- c(rep(apps, times = years) + t_add, tail(t_add,1) + 365)
  remsl <- rep(rems, length(appsl))
  
  if(from == 'storage'){
    for (i in appsl){
      slurry$slurry_mass[slurry$time > i] <- 
      slurry$slurry_mass[slurry$time > i] - 
      slurry$slurry_mass[slurry$time == i] * remsl[[which(appsl == i)]]
    }
    return(slurry_mass_dat = slurry)
  }
  
  if(from == 'digestate'){
    for (i in appsl){
      slurry$slurry_mass[slurry$time > i] <- 
        slurry$slurry_mass[slurry$time > i] - 
        slurry$slurry_mass[slurry$time == i] * remsl[[which(appsl == i)]]
    }
    return(slurry_digestate_dat = slurry)
  }
}