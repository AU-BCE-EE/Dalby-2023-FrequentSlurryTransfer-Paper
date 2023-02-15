
abm <- function(
  days = 365,                                # Number of days to run
  delta_t = 1,                               # Time step for output
  times = NA,
  wthr_pars = list(temp_air_C = 20, RH = 90, rain = 1.9, pres_kpa = 101, rs = 20),
  evap_pars = list(evap = 0.5 * et(temp_C = wthr_pars$temp_air_C, pres_kpa = wthr_pars$pres_kpa, rs = wthr_pars$rs)),                # mm/d
  mng_pars = list(slurry_prod_rate = 1000,   # kg/d
                  slurry_mass = 1000,           # Initial slurry mass (kg) NTS: convert to depth??
                  storage_depth = 4,         # Storge structure depth, assued to be maximum slurry depth (m)
                  resid_depth = 0.2,         # Residual slurry depth after emptying (m)
                  floor_area = 11,           # Currently for NH3 loss from barn floor (nothing to do with pit/tank floor)
                  area = 11,                 # Area (assume vertical sides) (m2)
                  empty_int = 35,            # (d)
                  temp_C = 20,
                  wash_water = 0,            
                  wash_int = NA,
                  rest_d = 0,
                  cover = NA,
                  resid_enrich = 0.9,
                  slopes = c(urea = NA, slurry_prod_rate = NA),
                  scale = c(ks_coefficient = 1.170719751, qhat_opt = 0.316190792, xa_fresh = 3.142802667, yield = 1, alpha_opt = 2.477011426)),
  grz_pars = list(graze_start = "may",
                    graze_days = 0,
                    graze_hours = 0),
  man_pars = list(conc_fresh = list(S2 = 0.01, urea = 2.4, SO4 = 0.2, TAN = 0.6,
                                    VFA = 2.83, xa_dead = 0, VSd = 75), pH = 7, dens = 1000), 
  grp_pars = list(grps = c('m0','m1','m2', 'sr1'),
                  yield = c(default = 0.05, sr1 = 0.065),
                  xa_fresh = c(default = 0.02),
                  xa_init = c(all = 0.0001),
                  decay_rate = c(all = 0.02),
                  ks_coefficient = c(default = 1, sr1 = 0.4),
                  qhat_opt = c(m0 = 1.5, m1 = 3.6, m2 = 5.6 , m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8),
                  T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                  T_min = c(m0 = 0, m1 = 10, m2 = 10, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                  T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                  ki_NH3_min = c(all = 0.015),
                  ki_NH3_max = c(all = 0.13),
                  ki_NH4_min = c(all = 2.7),
                  ki_NH4_max = c(all = 4.8),
                  pH_upr = c(all = 8.0),
                  pH_lwr = c(default = 6.0)),
  mic_pars = list(ks_SO4 = 0.0067,
                  ki_H2S_meth = 0.23,
                  ki_H2S_sr = 0.25,
                  alpha_opt = c(xa_dead= 0.02, urea = 70, VSd = 0.02),
                  alpha_T_min = c(xa_dead= 0, urea = 0, VSd = 0),
                  alpha_T_opt = c(xa_dead= 50, urea = 50, VSd = 50),
                  alpha_T_max = c(xa_dead= 60, urea = 60, VSd = 60)),
  chem_pars = list(COD_conv = c(CH4 = 0.2507, xa_dead = 0.73, VFA = 0.93, S = 0.5015, VS = 0.69, CO2_anaer = 0.53, CO2_aer = 1.1, CO2_sr = 1.2, CO2_ureo = 1.57,
                                C_xa_dead = 0.358, C_VFA = 0.374, C_VSd = 0.344, C_N_urea = 0.429), 
                   kl = c(NH3 = 52, NH3_floor = 22, H2S = 0.02)), 
  
  
  # kl = mass transfer coefficient (liquid phase units) in m/d
  add_pars = NULL,
  pars = NULL,
  startup = -Inf,
  starting = NULL,
  approx_method_temp = 'linear', # NTS: combine these in 1 vector?
  approx_method_pH = 'linear',
  approx_method_vent_air = 'linear',
  par_key = '\\.',
  value = 'ts',                              # Type of output
  warn = TRUE) {
  
  
  # If starting conditions are provided from a previous simulation, move them to pars
  # Note that additional state variables are extracted from starting in abm_*.R
  if (!is.null(starting) & is.data.frame(starting)) {
    message('Using starting conditions from `starting` argument')
    grp_pars[['xa_init']] <- as.numeric(starting[nrow(starting), paste0(names(grp_pars[['qhat_opt']]), '_conc')])
    names(grp_pars[['xa_init']]) <- names(grp_pars[['qhat_opt']])
    mng_pars['slurry_mass'] <- starting[nrow(starting), 'slurry_mass']
  }
  
  # Combine pars to make extraction and pass to rates() easier
  if (is.null(pars)) {
    pars <- c(wthr_pars, evap_pars, mng_pars, man_pars, grp_pars, mic_pars, chem_pars, grz_pars, list(days = days))
  }

  # Create error if batch time is not determined and slurry rate should increase over a batch
  if (is.na(pars$wash_int) && !is.na(pars$slopes['slurry_prod_rate'])){
    stop('wash interval cannot be "NA" when slurry production rate increase over a batch of animals, since wash interval defines a batch time')
  }
  # Sort out parameter inputs
  # Note: pe.pars = add_pars that use par.element approach, these are converted to normal (simple) add_par elements here
  # Note: sa.pars = normal (simple) add_pars that do not need to be converted
  # Note: Use of [] vs. [[]] affect how code works--needs to work for both lists and vector pars
  # Note: par.element approach is only designed to work for vector elements
  if (!is.null(add_pars) && length(add_pars) > 0 && any(ii <- grepl(par_key, names(add_pars)))) {
    pe.pars <- add_pars[ii]
    sa.pars <- add_pars[!ii]
    apnames <- names(pe.pars)
    pe.pars <- as.numeric(pe.pars)
    split.pars <- strsplit(apnames, par_key)
    pnames <- sapply(split.pars, '[[', 1)
    enames <- sapply(split.pars, '[[', 2)
    names(pe.pars) <- enames
    pe.pars <- split(pe.pars, pnames)
    add_pars <- c(sa.pars, pe.pars)
  }

  # If any additional parameters were added (or modified) using add_pars, update them in pars list here
  # Needs to work in a case where default is all but e.g., m1 is given in add_pars (see def stuff below)
  grp_par_nms <- names(grp_pars)
  grp_par_nms <- grp_par_nms[grp_par_nms != 'grps']
  if (!is.null(add_pars) && length(add_pars) > 0) {
    if (any(bad.names <- !names(add_pars) %in% names(pars))) {
      stop ('Some `add_pars` names not recognized as valid parameters: ', names(add_pars)[bad.names]) 
    }
    # Add in pars (or replace existing elements unless it is time series data added)
    for (i in names(add_pars)) {
      if (!is.data.frame(add_pars[[i]]) && length(pars[[i]]) > 1) {
        pars[[i]][names(add_pars[[i]])] <- add_pars[[i]]
      } else {
        def <- pars[[i]]['all']
        pars[[i]] <- add_pars[[i]]
        if (i %in% grp_par_nms) {
          pars[[i]]['default'] <- def
        }
      }
    }
  }
  
  # Unlike others, grps in add_pars *will* override default vector (i.e., can be used to remove groups)
  if ('grps' %in% names(add_pars)) {
    pars$grps <- add_pars$grps
  }
  
  # NTS: Add comment on this (FD?)
  
  # Below code is used only when we have variable concentration of methanogens in the fresh slurry.
  # In that case xa_fresh['time'] will need to be defined in a data.frame, but the block starting below this one 
  # will remove xa_fresh["time"], so we save it as xa_fresh_time before this happens. xa_fresh_time is used in L189 
  
  if(is.data.frame(pars$xa_fresh)) xa_fresh_time <- pars$xa_fresh["time"]
  
  # Fill in default values for grp_pars if keyword name `default` or `all` is used
  # Note: Microbial groups are defined by grps element
  # Note: `default` does *not* work with add_pars argument because grps are already defined in defaults
  # Note: But `all` *does* work
  grp_nms <- pars$grps
  for (i in grp_par_nms) {
    ppo <- pars[[i]]
    p_nms <- names(pars[[i]])
    if (any(p_nms == 'default')) {
      pars[[i]][grp_nms] <- pars[[i]]['default']
      if (any(p_nms != 'default')) {
        pars[[i]][p_nms[p_nms != 'default']] <- ppo[p_nms[p_nms != 'default']]
      }
    }
    if (any(p_nms == 'all')) {
      pars[[i]][grp_nms] <- pars[[i]]['all']
    }
    # Fix order, drop default element if present, drop unused names
    pars[[i]] <- pars[[i]][grp_nms]
    # Check for missing values
    if (any(is.na(pars[[i]]))) stop('Missing grp_pars elements in ', i, '.')
  }
  
  # Check grp arguments, including order of element names in some pars
  # After above block, this should be redundant
  checkGrpNames(pars)
  
  # Convert temperature constants to K if needed
  pars <- tempsC2K(pars, ll = 200)
  
  # Create time-variable functions
  # Note that pars$x must be numeric constant or df with time (col 1) and var (2)
  temp_C_fun <- makeTimeFunc(pars$temp_C, approx_method = approx_method_temp)
  temp_air_C_fun <- makeTimeFunc(pars$temp_air_C, approx_method = approx_method_temp)
  pH_fun <- makeTimeFunc(pars$pH, approx_method = approx_method_pH)
  
  # Inhibition function
  td <- data.frame(SO4_conc = c(0, 0.392251223, 0.686439641, 1.046003263, 1.470942088, 1.961256117, 4), 
                   H2SO4_inhib = c(1, 0.5407, 0.3372, 0.2326, 0.1092, 0.0442, 0.001))
  H2SO4_inhibition_fun <- approxfun(td$SO4_conc, td$H2SO4_inhib, rule = 2)
  
  # when variable conc fresh is used, conc_fresh_start needs to be defined to pass for y (initial condition) vector in abm reg or abm variable 
  if (is.data.frame(pars$conc_fresh)){
    conc_fresh_starting <- pars$conc_fresh[1, !grepl("time", colnames(pars$conc_fresh))]
    names(conc_fresh_starting) <- gsub(pattern = "conc_fresh.", replacement = "", x = names(conc_fresh_starting))
    pars$conc_fresh_starting <- as.list(conc_fresh_starting)
  } else {
    pars$conc_fresh_start <- pars$conc_fresh
  }
  
  # when variable xa_fresh is used xa_init should be the same as starting xa_fresh
  if (is.data.frame(pars$xa_fresh)) {
    pars$xa_init <- as.numeric(pars$xa_fresh[1, grepl("m[0-9]|sr[0-9]", names(pars$xa_fresh))])
    names(pars$xa_init) <- names(pars$xa_fresh)[!grepl("time", names(pars$xa_fresh))]
    pars$xa_fresh["time"] <- xa_fresh_time
  }

  # Convert some supplied parameters
  # Maximum slurry mass in kg
  pars$max_slurry_mass <- pars$storage_depth * pars$area * pars$dens
  pars$resid_mass <- pars$resid_depth / pars$storage_depth * pars$max_slurry_mass

  # Cover effect on NH3 emission rate
  # reduction from cover 
  pars$EF_NH3 <- coverfun(pars$cover)
  
  # calculate grazing interval of year if needed
  if(pars$graze_days > 0){
    pars$graze_int <- c(doy(pars$graze_start)$day, doy(pars$graze_start)$day + pars$graze_days)
  }
  
  if(is.data.frame(pars$slurry_mass)){
    # If missing, set initial slurry mass to 0
    if (pars$slurry_mass[1, 'time'] > 0) {
      pars$slurry_mass <- rbind(c(0, 0), pars$slurry_mass)
    }
    slurry_mass_init <- pars$slurry_mass[1, 'slurry_mass']
  } else {
    slurry_mass_init <- pars$slurry_mass
  }
  
  if (slurry_mass_init == 0) {
    slurry_mass_init <- 1E-10
  }
  
  y <- c(xa = pars$xa_init * slurry_mass_init,
         slurry_mass = slurry_mass_init, 
         xa_dead = pars$conc_fresh_start[['xa_dead']] * slurry_mass_init, 
         VSd = pars$conc_fresh_start[['VSd']] * slurry_mass_init,
         starch = pars$conc_fresh_start[['starch']] * slurry_mass_init,
         VFA = pars$conc_fresh_start[['VFA']] * slurry_mass_init, 
         urea = pars$conc_fresh_start[['urea']] * slurry_mass_init, 
         TAN = pars$conc_fresh_start[['TAN']] * slurry_mass_init, 
         sulfate = pars$conc_fresh_start[['SO4']] * slurry_mass_init, 
         sulfide = pars$conc_fresh_start[['S2']] * slurry_mass_init, 
         NH3_emis_cum = 0, CH4_emis_cum = 0, CO2_emis_cum = 0, 
         COD_conv_cum = 0, COD_conv_cum_meth = 0, 
         COD_conv_cum_respir = 0, COD_conv_cum_sr = 0)
  
  if (!is.null(starting) & is.data.frame(starting)) {
    start.vars <- c('slurry_mass', 'xa_dead', 'VSd', 'VFA', 'urea', 'TAN', 'sulfate', 'sulfide')
    y[start.vars]  <- starting[nrow(starting), start.vars]
  }  
  
  if (is.numeric(pars$slurry_mass)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- abm_regular(days = days, delta_t = delta_t, y = y, pars = pars, starting = starting, temp_C_fun = temp_C_fun, temp_air_C_fun = temp_air_C_fun, pH_fun = pH_fun,  H2SO4_inhibition_fun = H2SO4_inhibition_fun)
  } else if (is.data.frame(pars$slurry_mass)) {
    # Option 2: Everything based on given slurry mass vs. time
    dat <- abm_variable(days = days, delta_t = delta_t, times = times, y = y, pars = pars, warn = warn, temp_C_fun = temp_C_fun, temp_air_C_fun = temp_air_C_fun, pH_fun = pH_fun, H2SO4_inhibition_fun = H2SO4_inhibition_fun)
  } 
  
  colnames(dat) <- gsub("conc_fresh.","conc_fresh_", colnames(dat))
  
  # Calculate concentrations where relevant
  #conc.names <- names(dat)[!grepl('conc|time|slurry_mass|inhib|qhat|CH4_emis_cum', names(dat))]
  mic_names <- pars$grps
  eff_names <- names(dat[grepl("_eff$", names(dat))])
  eff_conc_names <- eff_names[eff_names != "slurry_mass_eff"]
  conc_names <-  c('TAN', 'xa_dead', 'urea', 'VSd', 'VFA', 'sulfide', 'sulfate', mic_names)
  dat_conc <- dat[, conc_names]/(dat$slurry_mass)
  dat_eff_conc <- dat[, eff_conc_names]/(dat$slurry_mass_eff)
  names(dat_conc) <- paste0(names(dat_conc), '_conc')
  names(dat_eff_conc) <- paste0(names(dat_eff_conc), '_conc')
  dat <- cbind(dat, dat_conc, dat_eff_conc)
  
  # Add temperature and pH
  dat$temp_C <- temp_C_fun(dat$time)
  dat$temp_air_C <- temp_air_C_fun(dat$time)
  if (is.numeric(pars$pH) | is.data.frame(pars$pH)) {
    dat$pH <- pH_fun(dat$time)
  } else if (pars$pH == 'calc'){
    dat$pH <- H2SO4_titrat(dat$sulfate_conc, class_anim = "pig")$pH
  } else {
    stop('Problem with pH input (bee721)')
  }
  
  # Calculate rates etc. for output, from state variables
  # NTS: check units, use dens???
  dat$rNH3 <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))
  dat$NH3_emis_rate <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))
  dat$NH3_emis_rate_slurry <- dat$NH3_emis_rate / (dat$slurry_mass / 1000)
  dat$NH3_flux <- dat$NH3_emis_rate / pars$area
  dat$rCH4 <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  dat$rCH4_A <- c(0, diff(dat$CH4_A_emis_cum))/c(1, diff(dat$time))
  dat$CH4_emis_rate <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  dat$CH4_emis_rate_slurry <- dat$CH4_emis_rate / (dat$slurry_mass / 1000)
  dat$CH4_A_emis_rate <- c(0, diff(dat$CH4_A_emis_cum))/c(1, diff(dat$time))
  dat$CH4_flux <- dat$CH4_emis_rate / pars$area
  dat$CO2_emis_rate <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
  dat$CO2_emis_rate_slurry <- dat$CO2_emis_rate / (dat$slurry_mass / 1000)
  dat$CO2_flux <- dat$CO2_emis_rate / pars$area
  
  ## NTS: Add others, e.g., mu
  
  # Cut startup period before calculating cumulative flows
  dat <- dat[dat$time > startup, ]
  
  # Calculate COD/VS flows
  # First concentrations in g/kg
  dat$dCOD_conc_fresh <- dat$conc_fresh_VFA + dat$conc_fresh_xa_dead + dat$conc_fresh_VSd + rowSums(dat[, grepl("xa_fresh_", colnames(dat))])
  dat$COD_conc_fresh <- dat$conc_fresh_VFA + dat$conc_fresh_xa_dead + dat$conc_fresh_VSd + rowSums(dat[, grepl("xa_fresh_", colnames(dat))]) 
  dat$ndCOD_conc_fresh <- dat$COD_conc_fresh - dat$dCOD_conc_fresh
  dat$C_conc_fresh <- dat$conc_fresh_VFA * pars$COD_conv[['C_VFA']] + dat$conc_fresh_xa_dead * pars$COD_conv[['C_xa_dead']] +
                      dat$conc_fresh_VSd * pars$COD_conv[['C_VSd']] + rowSums(dat[, grepl("xa_fresh_", colnames(dat))]) * pars$COD_conv[['C_xa_dead']] +
                      dat$conc_fresh_urea * pars$COD_conv[['C_N_urea']]
  
  # ndCOD is almost conserved, same everywhere always
  dat$ndCOD_conc <- ndCOD_conc <- dat$COD_conc_fresh - dat$dCOD_conc_fresh 
  dat$dCOD_conc <- dCOD_conc <- dat$xa_dead_conc + dat$VFA_conc + dat$VSd_conc + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE])
  dat$COD_conc <- COD_conc <- ndCOD_conc + dCOD_conc
  dat$VS_conc <- pars$COD_conv[['VS']] * COD_conc # this needs to be updated. Note this is both deg and undeg VS, whereas if VSd is given it is only deg VS
  dat$C_conc <- dat$VFA_conc * pars$COD_conv[['C_VFA']] + dat$xa_dead_conc * pars$COD_conv[['C_xa_dead']] +
                dat$VSd_conc * pars$COD_conv[['C_VSd']]  + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE]) * pars$COD_conv[['C_xa_dead']] +
                dat$urea_conc * pars$COD_conv[['C_N_urea']]
  
  dat$C_eff_conc <- dat$VFA_eff_conc * pars$COD_conv[['C_VFA']] + dat$xa_dead_eff_conc * pars$COD_conv[['C_xa_dead']] +
                    dat$VSd_eff_conc * pars$COD_conv[['C_VSd']] + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE]) * pars$COD_conv[['C_xa_dead']] +
                    dat$urea_eff_conc * pars$COD_conv[['C_N_urea']]
  
  
  dat$Ninorg_conc_fresh <- dat$conc_fresh_urea + dat$conc_fresh_TAN
  dat$Norg_conc_fresh <- 0
  dat$N_conc_fresh <- dat$Ninorg_conc_fresh + dat$Norg_conc_fresh
  dat$Ninorg_conc <- Ninorg_conc <- dat$urea_conc + dat$TAN_conc
  dat$Norg_conc <- Norg_conc <- 0
  dat$N_conc <- N_conc <- dat$Ninorg_conc + dat$Norg_conc
  dat$Ninorg_eff_conc <- Ninorg_eff_conc <- dat$urea_eff_conc + dat$TAN_eff_conc
  dat$N_eff_conc <- Ninorg_eff_conc

  
  # And flows in g/d
  dat$COD_load_rate <- dat$COD_conc_fresh * dat$slurry_prod_rate
  dat$dCOD_load_rate <- dat$dCOD_conc_fresh * dat$slurry_prod_rate
  dat$ndCOD_load_rate <- dat$ndCOD_conc_fresh * dat$slurry_prod_rate
  dat$VS_load_rate <- pars$COD_conv[['VS']] * dat$COD_load_rate
  dat$C_load_rate <- dat$C_conc_fresh * dat$slurry_prod_rate
  
  dat$Ninorg_load_rate <- dat$Ninorg_conc_fresh * dat$slurry_prod_rate
  dat$Norg_load_rate <- dat$Norg_conc_fresh * dat$slurry_prod_rate
  dat$N_load_rate <- dat$N_conc_fresh * dat$slurry_prod_rate
  
  # Cumulative flow in g
  dat$COD_load_cum <- cumsum(dat$COD_load_rate * c(0, diff(dat$time))) + dat$COD_conc[1] * dat$slurry_mass[1]
  dat$dCOD_load_cum <- cumsum(dat$dCOD_load_rate * c(0, diff(dat$time))) + dat$dCOD_conc[1]* dat$slurry_mass[1]
  dat$ndCOD_load_cum <- cumsum(dat$ndCOD_load_rate * c(0, diff(dat$time))) + dat$ndCOD_conc[1]* dat$slurry_mass[1]
  dat$VS_load_cum <- cumsum(dat$VS_load_rate * c(0, diff(dat$time))) + dat$VS_conc[1]* dat$slurry_mass[1]
  dat$C_load_cum <- cumsum(dat$C_load_rate * c(0, diff(dat$time))) + dat$C_conc[1] * dat$slurry_mass[1]
  
  dat$Ninorg_load_cum <- cumsum(dat$Ninorg_load_rate * c(0, diff(dat$time))) + dat$Ninorg_conc[1] * dat$slurry_mass[1]
  dat$Norg_load_cum <- cumsum(dat$Norg_load_rate * c(0, diff(dat$time))) + dat$Norg_conc[1] * dat$slurry_mass[1]
  dat$N_load_cum <- cumsum(dat$N_load_rate * c(0, diff(dat$time))) + dat$N_conc[1] * dat$slurry_mass[1]

  # And relative emission
  # g CH4/g COD in
  
  dat$CH4_emis_rate_COD <- dat$CH4_emis_rate / dat$COD_load_rate
  dat$CH4_emis_rate_dCOD <- dat$CH4_emis_rate / dat$dCOD_load_rate
  dat$CH4_emis_rate_VS <- dat$CH4_emis_rate / dat$VS_load_rate
  dat$CH4_emis_rate_C <- dat$CH4_emis_rate / dat$C_load_rate
  dat$CH4_emis_cum_COD <- dat$CH4_emis_cum / dat$COD_load_cum
  dat$CH4_emis_cum_dCOD <- dat$CH4_emis_cum / dat$dCOD_load_cum
  dat$CH4_emis_cum_VS <- dat$CH4_emis_cum / dat$VS_load_cum
  dat$CH4_emis_cum_C <- dat$CH4_emis_cum / dat$C_load_cum
  
  # g NH3 N/g N
  dat$NH3_emis_rate_Ninorg <- dat$NH3_emis_rate / dat$Ninorg_load_rate
  dat$NH3_emis_rate_Norg <- dat$NH3_emis_rate / dat$Norg_load_rate
  dat$NH3_emis_rate_N <- dat$NH3_emis_rate / dat$N_load_rate
  dat$NH3_emis_cum_Ninorg <- dat$NH3_emis_cum / dat$Ninorg_load_cum
  dat$NH3_emis_cum_Norg <- dat$NH3_emis_cum / dat$Norg_load_cum
  dat$NH3_emis_cum_N <- dat$NH3_emis_cum / dat$N_load_cum
  
  # Same for CO2
  dat$CO2_emis_rate_COD <- dat$CO2_emis_rate / dat$COD_load_rate
  dat$CO2_emis_rate_dCOD <- dat$CO2_emis_rate / dat$dCOD_load_rate
  dat$CO2_emis_rate_VS <- dat$CO2_emis_rate / dat$VS_load_rate
  dat$CO2_emis_rate_C <- dat$CO2_emis_rate / dat$C_load_rate
  dat$CO2_emis_cum_COD <- dat$CO2_emis_cum / dat$COD_load_cum
  dat$CO2_emis_cum_dCOD <- dat$CO2_emis_cum / dat$dCOD_load_cum
  dat$CO2_emis_cum_VS <- dat$CO2_emis_cum / dat$VS_load_cum
  dat$CO2_emis_cum_C <- dat$CO2_emis_cum / dat$C_load_cum
  
  # Apparent COD conversion fraction
  dat$f_COD_CH4_rate <-  dat$CH4_emis_rate / pars$COD_conv[['CH4']] / dat$COD_load_rate
  dat$f_COD_CH4_cum <-  dat$COD_conv_cum_meth / dat$COD_load_cum
  dat$f_COD_respir_cum <-  dat$COD_conv_cum_respir / dat$COD_load_cum
  dat$f_COD_sr_cum <-  dat$COD_conv_cum_sr / dat$COD_load_cum
  
  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))
  
  Ninorg_load <- dat$Ninorg_load_cum[nrow(dat)]
  Norg_load <- dat$Norg_load_cum[nrow(dat)]
  N_load <- dat$N_load_cum[nrow(dat)]
  
  COD_load <- dat$COD_load_cum[nrow(dat)]
  dCOD_load <- dat$dCOD_load_cum[nrow(dat)]
  ndCOD_load <- dat$ndCOD_load_cum[nrow(dat)]
  VS_load <- dat$VS_load_cum[nrow(dat)]
  C_load <- dat$C_load_cum[nrow(dat)]
  
  NH3_emis_cum <- dat$NH3_emis_cum[nrow(dat)]
  NH3_emis_rate <- NH3_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  NH3_emis_Ninorg <- NH3_emis_cum / Ninorg_load
  NH3_emis_Norg <- NH3_emis_cum / Norg_load
  NH3_emis_N <- NH3_emis_cum / N_load
  
  CH4_emis_cum <- dat$CH4_emis_cum[nrow(dat)]
  CH4_A_emis_cum <- dat$CH4_A_emis_cum[nrow(dat)]
  CH4_emis_rate <- CH4_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CH4_A_emis_rate <- CH4_A_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CH4_emis_COD <- CH4_emis_cum / COD_load
  CH4_emis_dCOD <- CH4_emis_cum / dCOD_load
  CH4_emis_VS <- CH4_emis_cum / VS_load
  CH4_emis_C <- CH4_emis_cum / C_load
  
  CO2_emis_cum <- dat$CO2_emis_cum[nrow(dat)] - dat$CO2_emis_cum[1]
  CO2_emis_rate <- CO2_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CO2_emis_COD <- CO2_emis_cum / COD_load
  CO2_emis_dCOD <- CO2_emis_cum / dCOD_load
  CO2_emis_VS <- CO2_emis_cum / VS_load
  CO2_emis_C <- CO2_emis_cum / C_load
  
  COD_conv_meth <- dat$COD_conv_cum_meth[nrow(dat)] - dat$COD_conv_cum_meth[1]
  COD_conv_respir <- dat$COD_conv_cum_respir[nrow(dat)] - dat$COD_conv_cum_respir[1]
  COD_conv_sr <- dat$COD_conv_cum_sr[nrow(dat)] - dat$COD_conv_cum_sr[1]
  f_COD_CH4 <- COD_conv_meth / COD_load
  f_COD_respir <- COD_conv_respir / COD_load
  f_COD_sr <- COD_conv_sr / COD_load
  
  summ <- c(C_load = C_load, COD_load = COD_load, dCOD_load = dCOD_load, ndCOD_load = ndCOD_load, VS_load = VS_load, Ninorg_load = Ninorg_load, Norg_load = Norg_load, 
            N_load = N_load, NH3_emis_cum = NH3_emis_cum, NH3_emis_rate = NH3_emis_rate, NH3_emis_Ninorg = NH3_emis_Ninorg, NH3_emis_Norg = NH3_emis_Norg,
            NH3_emis_N = NH3_emis_N, CH4_emis_cum = CH4_emis_cum, CH4_emis_rate = CH4_emis_rate, CH4_A_emis_rate = CH4_A_emis_rate, CH4_A_emis_cum = CH4_A_emis_cum, 
            CH4_emis_COD = CH4_emis_COD, CH4_emis_dCOD = CH4_emis_dCOD, CH4_emis_VS = CH4_emis_VS, CH4_emis_C = CH4_emis_C,
            CO2_emis_cum = CO2_emis_cum, CO2_emis_rate = CO2_emis_rate, CO2_emis_COD = CO2_emis_COD, CO2_emis_dCOD = CO2_emis_dCOD, CO2_emis_VS = CO2_emis_VS, CO2_emis_C = CO2_emis_C,
            COD_conv_meth = COD_conv_meth, COD_conv_respir = COD_conv_respir, COD_conv_sr = COD_conv_sr,
            f_COD_CH4 = f_COD_CH4, f_COD_respir = f_COD_respir, f_COD_sr = f_COD_sr) 

  # Add slurry depth
  dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens

  # Check slurry mass, warn if needed
  if (any(dat$slurry_mass > pars$max_slurry_mass)) {
    warning('Maximum slurry mass exceeded.\nCheck output.')
  }
  
  # Return results
  # Average only
  if (substring(value, 1, 3) == 'sum') return(summ)
  # ts = time series
  if (value == 'ts') return(dat)
  # tsel = time series after startup period
  if (value == 'tsel') return(dat_sel)
  # Or everything
  return(list(pars = pars, ts = dat, tsel = dat_sel, summ = summ))
  
}

