# pH-dependent chemical speciation 

eqSpec <- function(tot = tot, temp.c = 20, of = 'a', ll = -14, ul = 0, pH, adjpH = 'H2CO3') {

  # Check arguments
  if(length(tot) != 7 || any(!names(tot)%in%c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc'))) stop('Check tot argument. Should be a numeric vector with elements ', "'H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc'")

  # Define other functions
  pHSpec <- function(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z) {
  # Iteratively solve, updating ionic strength each time
    i <- sum(0.5*tot[tot>0]) # Very crude guess
    b <- 999
    k <- 10^l.k
    l.a <- 0*l.k # Just for length and names
    l.a['H.'] <- l.a.h
    a <- 10^l.a
    j <- 0  
    di <- 999
    while (di/i>1E-4){ #abs(log10(i/b))>log10(1.001)) {
      j <- j + 1
      b <- i
      l.g<- -a.dh*z^2*sqrt(i)/(1+b.dh*a.par[1, ]*sqrt(i)) + a.par[2, ]*i
      g <- 10^l.g
  
      l.a['K.']   <-log10(tot['K.']) + l.g['K.']
      l.a['Na.']  <-log10(tot['Na.']) + l.g['Na.']
      l.a['Cl.']  <-log10(tot['Cl.']) + l.g['Cl.']
      l.a['NH4.'] <-log10(tot['NH3']*k['NH4.']*g['NH3']*a['H.']*g['NH4.']/(g['NH4.'] + k['NH4.']*g['NH3']*a['H.']) )
      l.a['NH3']  <-l.a['NH4.'] - l.a['H.'] - l.k['NH4.']
      l.a['H2CO3'] <- log10(tot['H2CO3']*a['H.']*g['H2CO3']/(k['HCO3.']*g['H2CO3']/g['HCO3.'] + a['H.'] + k['CO3.2']*g['H2CO3']/(g['CO3.2']*a['H.']) + a['H.']*k['CO2']*g['H2CO3']/g['CO2'] ) )
      l.a['HCO3.'] <- l.k['HCO3.'] + l.a['H2CO3'] - l.a['H.']
      l.a['CO3.2'] <- l.k['CO3.2'] + l.a['H2CO3'] - 2*l.a['H.']
      l.a['CO2']  <-l.k['CO2'] + l.a['H2CO3']
      l.a['OH.']  <- 0 -l.a['H.'] + l.k['OH.']
      l.a['Ac.']  <- log10(k['Ac.']*g['HAc']*tot['HAc']*g['Ac.']/(a['H.']*g['Ac.'] + k['Ac.']*g['HAc']))
      l.a['HAc']  <- l.a['Ac.'] + l.a['H.'] - l.k['Ac.']
  
      l.m <- l.a - l.g
      m <- 10^l.m
      i <- sum(0.5*m*z^2)
      di <- abs(i - b)
    }
    a <- 10^l.a
    cb <- sum(z*m)
	totk <- c(tot[1:2], m[3:4], tot[4:7])
	totk[3] <- totk[3] + m['HCO3.'] + m['CO3.2']
    list(m = m, a = a, g = g, i = i, l.m = l.m, l.a = l.a, l.g = l.g, tot = tot, totk = totk, cb = cb, i.its = j)
  }

  # Calculates error in total H for a given pH guess
  HBalErr <- function(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z) {
    m <- pHSpec(l.a.h = l.a.h, a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z)$m
    abs(tot['H.'] - sum(s[, 1]*m)) # All calculated totals given by t(s)%*%m. Alternative would be charge balance.
  }

  # Now code for speciation
  if(!exists('n.calls.spec')) n.calls.spec <<- 0
  n.calls.spec <<- n.calls.spec + 1
  # Temperature
  temp.k <- 273.15+temp.c

  # Henry's law constants
  kh.CO2<- 10^(108.38578 +0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
  kh.NH3<- 10^(-3.51645 -0.0013637*temp.k + 1701.35/temp.k)
 
  # Equilibrium constants
  temp.k <- 273.15+temp.c
  l.k <- c('H.' = 0, 'NH3' = 0, 'H2CO3' = 0, 'CO2' = 2.778151, 'K.' = 0, 'Na.' = 0, 'Cl.' = 0, 'HAc' = 0, 
      'OH.'= -4.2195 -2915.16/temp.k, 
      'NH4.'= 0.0905 + 2729.31/temp.k, 
      'HCO3.'= -353.5305 -0.06092*temp.k + 21834.37/temp.k + 126.8339*log10(temp.k) -1684915/temp.k^2, 
      'CO3.2'= -461.4176 -0.093448*temp.k + 26986.16/temp.k + 165.7595*log10(temp.k) -2248629/temp.k^2, 
      'Ac.' = -4.8288 + 21.42/temp.k)
  # Matrix of stoichiometric coefficients. Not used for much, but is good for keeping track of reactions
  s <- matrix(c(1, 0,0, 0,0, 0,0, 0,-1, 1,-1, -2, -1,  0, 1,0, 0,0, 0,0, 0,0, 1,0, 0,0,  0, 0,1, 1,0, 0,0, 0,0, 0,1, 1,0, 
              0, 0,0, 0,1, 0,0, 0,0, 0,0, 0,0,      0, 0,0, 0,0, 1,0, 0,0, 0,0, 0,0,  0, 0,0, 0,0, 0,1, 0,0, 0,0, 0,0,  
	      0, 0,0, 0,0, 0,0, 1,0, 0,0, 0,1), 
            nrow = 13, dimnames = list(c("H.", "NH3", 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', "OH.", "NH4.", 'HCO3.', 'CO3.2', 'Ac.'), 
                                  c("H.", 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc'))
  )
  # Species charge
  z <- c('H.' = 1, 'NH3' = 0, 'H2CO3' = 0, 'CO2' = 0, 'K.' = 1, 'Na.' = 1, 'Cl.' = -1, 'HAc' = 0, 'OH.' = -1, 'NH4.' = 1, 'HCO3.' = -1, 'CO3.2' = -2, 'Ac.' = -1)
  # Parameters for calculation of activity coefficients, using extended Debye-Huckel equation (see PHREEQC manual)
  a.par <- matrix(c(9, 0,0, 0.1, 0,0.1, 0,0.1, 3,0, 4.25, 0,3, 0,0, 0.1, 3.5, 0,2.5, 0,5.4, 0,4.5, 0,4.5, 0), nrow = 2, dimnames = list(c('a', 'b'), c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))

  # Calculate a.dh and b.dh parameters for Debye-Huckel equation. Dielectric constant, density rearranged from PHREEQC code
  # Dielectric constant
  de <- 2727.586 + 0.6224107*temp.k - 1075.112*log10(temp.k) - 52000.87/temp.k
  # Water density
  c.d <- 647.26 - temp.k
  d.H2O <- (1 + .1342489*c.d^(1/3) - 3.946263E-3*c.d)/(3.1975 - 0.3151548*c.d^(1/3) - 1.203374E-3*c.d + 7.48908E-13*c.d^4)
  # a.dh and b.dh are A and B in Debye-Huckel equation. Equations are from Tuesdell & Jones (1974)
  a.dh <- 1.82483E6*d.H2O^0.5/(de*temp.k)^(3/2)
  b.dh <- 50.2916*d.H2O^0.5/(de*temp.k)^0.5

  # If pH not specified (assumed to be typical case) it is calculated, if specified, KOH and HAc is added to reach specified value, based on proton balance
  if(missing(pH)) {
    sol <- optimize(f = HBalErr, interval = c(ll, ul), a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z, tol = 1E-15)
    if(sol$objective>5E-7) {
      print(sol)
      stop('Around line 160, optimize didn\'t converge: complete results above, specified limits: ', ll, ' ', ul)
    }
    l.a.h <- sol$minimum
  } else { # pH is specified, acid or base adjusted to match it
    l.a.h<- -pH
    dhb <- 999
    hb <- tot['H.'] - sum(s[, 1]*pHSpec(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z)$m)
    while(dhb>1E-10) { # Must be solved iteratively because i will change
      if(adjpH == 'HAc') {
        tot['HAc'] <- tot['HAc'] - hb  
      } else if(adjpH == 'HCl') {
        tot['Cl.'] <- tot['Cl.'] - hb
        tot['H.'] <- tot['H.'] - hb
      } else if(adjpH == 'H2CO3') {
        tot['H2CO3'] <- tot['H2CO3'] - hb
      } else if(adjpH == 'KOH') {
          tot['K.'] <- tot['K.'] + hb 
          tot['H.'] <- tot['H.'] - hb
      } else if(adjpH == 'NaOH') {
          tot['Na.'] <- tot['Na.'] + hb 
          tot['H.'] <- tot['H.'] - hb
      } else if(adjpH == 'NH3') {
        tot['NH3'] <- tot['NH3'] + hb
       } else stop('adjpH argument must be "HAc", "H2CO3", "KOH", or "HCl" but is ', adjpH)
      hb2 <- tot['H.'] - sum(s[, 1]*pHSpec(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z)$m)
      dhb <- abs(hb2-hb)
      hb <- hb2
    }
  }

  out <- pHSpec(l.a.h, a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z)
  if(abs(out$cb)>5E-8) warning('Charge balance off in eqSpec. Check results. cb =', out$cb)
  tot.g <- t(s)%*%out$m
  p.CO2 <- out$a['CO2']/kh.CO2
  p.NH3 <- out$a['NH3']/kh.NH3
  names(p.CO2) <- NULL
  if(of == 'a') return(out$a)
  if(of == 'm') return(out$m)
  if(of == 'k') return(l.k)
  if(of == 'all') return(c(out, p.CO2 = p.CO2, p.NH3 = p.NH3))
}

