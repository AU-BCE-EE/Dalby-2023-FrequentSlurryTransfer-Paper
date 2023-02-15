# Functions for parameter estimation

resCalc <- function(p, dat, weights = 1, to, plot = FALSE, ...){
  # Use ... to pass fixed parameters

  # Extract observations, keeping missing values (fixed interval needed to compare to measurements)
  obs <- obs.orig <- dat[, to]
  times <- dat$day

  # Get days from times (no real effect in abm())
  days <- max(times)

  if (length(unique(diff(times))) > 1) {
    warning('Measurement data frame dat does not have fixed interval.')
  }

  if(!all(is.numeric(weights) | length(weights) > 0 | !all(is.na(weights)))) stop('Problem with weights argument.')
  if(any(is.na(weights)) || any(is.null(weights))) stop('weights are NA or NULL')
  if(length(weights) == 1) weights <- rep(weights, nrow(dat))

  obs[weights == 0] <- NA # To prevent warning on NaNs

  if(any(is.na(obs[weights > 0]))) stop('NA values in observations obs, not just where weights = 0.')

  # Back-transform parameters
  p <- 10^p
  out <- abm(days = days, times = times, delta_t = delta_t, add_pars = p, ...)
  pred <- pred.orig <- out[, to]

  # Center and scale measurements and predicted values by *measured* mean sd
  ss <- apply(obs, 2, sd, na.rm = TRUE)
  sm <- apply(obs, 2, mean, na.rm = TRUE)
  obs <- scale(obs, center = sm, scale = ss)
  pred <- scale(pred, center = sm, scale = ss)
  
  res <- weights*(pred - obs) 
  res[weights == 0] <- 0 

  if(any(is.na(res))) {
    warning('NA in residuals, see rows ', paste(which(is.na(res)), collapse = ', '), '.\nReplacing with 0.')
    res[is.na(res)] <- 0
  }

  obj <- sum(abs(res))

  cat('Parameter estimates: ', signif(p, 3), '\n')
  cat('Objective: ', signif(obj, 3), '\n')

  if (plot) {
    par(mfrow = c(length(to), 1))
    for (v in to) {
      plot(times, obs.orig[, v], ylab = v)
      lines(times, pred.orig[, v], ylab = v)
    }
  }

  return(obj)

}

