interpmult <- function(days, dat){

dat_approx <- data.frame(matrix(ncol = ncol(dat), nrow = 0))
colnames(dat_approx) <- colnames(dat) 
ys = colnames(dat[, unlist(lapply(dat, is.numeric))])
ys <- ys[-which(ys == "time")]

xout <- c(1:days)

for (i in ys){
dat_approx[xout, i] <- approx(dat[, "time"], dat[, i], xout = xout, rule = 2)$y
}

dat_approx[, "time"] <- xout
dat_approx[, colnames(select_if(dat, is.character))] <- select_if(dat, is.character)

if(ncol(dat_approx) != ncol(dat)){
  stop("some columns could not be interpolated")
}
  
return(dat_approx)
}
