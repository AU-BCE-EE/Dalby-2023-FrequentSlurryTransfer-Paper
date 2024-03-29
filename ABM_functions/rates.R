rates <- function(t, y, parms, temp_C_fun = temp_C_fun, temp_air_C_fun = temp_air_C_fun, pH_fun = pH_fun, H2SO4_inhibition_fun = H2SO4_inhibition_fun) {
  
     y[y < 1E-10] <- 1E-10
    
     days <- parms$days

     ks_SO4 <- parms$ks_SO4
     ki_H2S_meth <- parms$ki_H2S_meth
     ki_H2S_sr <- parms$ki_H2S_sr
     alpha_opt <- parms$alpha_opt
     alpha_T_opt <- parms$alpha_T_opt
     alpha_T_min <- parms$alpha_T_min
     alpha_T_max <- parms$alpha_T_max
     slurry_prod_rate <- parms$slurry_prod_rate
     floor_area <- parms$floor_area # for NH3 mass transfer equation
     area <- parms$area
     temp_C <- parms$temp_C
     temp_air_C <- parms$temp_air_C # for NH3 mass transfer equation
     wash_int <- parms$wash_int
     rest_d <- parms$rest_d
     rain <- parms$rain # kg/ m2 / day
     evap <- parms$evap # kg/ m2 / day
     EF_NH3 <- parms$EF_NH3
     pH <- parms$pH
    
     graze_int <- parms$graze_int
     graze_hours <- parms$graze_hours
    
     conc_fresh <- parms$conc_fresh
     slopes <- parms$slopes
     scale <- parms$scale # scale for optim qhat
     yield <- parms$yield
     xa_fresh <- parms$xa_fresh
     xa_init <- parms$xa_init
     decay_rate <- parms$decay_rate
     ks_coefficient <- parms$ks_coefficient
     resid_enrich <- parms$resid_enrich
     qhat_opt <- parms$qhat_opt
     T_opt <- parms$T_opt
     T_min <- parms$T_min
     T_max <- parms$T_max
     ki_NH3_min <- parms$ki_NH3_min
     ki_NH3_max <- parms$ki_NH3_max
     ki_NH4_min <- parms$ki_NH4_min
     ki_NH4_max <- parms$ki_NH4_max
     pH_upr <- parms$pH_upr
     pH_lwr <- parms$pH_lwr
     COD_conv <- parms$COD_conv
     kl<- parms$kl
    
     t_run <- parms$t_run
    
    # correct slurry production rate in periods with grazing
    
    if(!is.null(graze_int)) {
      slurry_prod_rate <- graze_fun(t,  t_run, days, slurry_prod_rate, graze_int, graze_hours)
    }
    
    
    # pH, numeric, variable, or from H2SO4
    if (is.numeric(pH) | is.data.frame(pH)) {
      pH <- pH_fun(t + t_run)
    } else if (pH == 'calc') {
      pH <- H2SO4_titrat(conc_SO4 = conc_fresh[['SO4']], class_anim = "pig")$pH
    } else {
      stop('pH problem (xi342)')
    }
    
    # calculate time of a batch
    if (!is.na(wash_int)){
      batches <- c(floor((t + t_run)/(wash_int + rest_d)))
      t_batch <- (t + t_run) - batches * (wash_int + rest_d)
      if (t_batch > wash_int) t_batch <- 0
    } else {
      t_batch <- 0
    }
    
    # if urea fresh increase during a batch 
    
    if (!is.na(slopes[['urea']]) & !is.data.frame(conc_fresh)) {
      start_urea <- conc_fresh[['urea']] - slopes[['urea']] * wash_int/2
      conc_fresh[['urea']] <- slopes[['urea']] * t_batch + start_urea
    }
    
    # if slurry production increase during a batch
    slurry_prod_rate_default <- slurry_prod_rate
    
    if (!is.na(slopes[['slurry_prod_rate']]) && slurry_prod_rate_default != 0) {
      start_slurry_prod_rate <- slurry_prod_rate_default - slopes[['slurry_prod_rate']] * wash_int/2
      slurry_prod_rate <- slopes[['slurry_prod_rate']] * t_batch + start_slurry_prod_rate
      if (t_batch > wash_int) slurry_prod_rate <- 0
    }

    # when variable fresh concentration is used
    if(is.data.frame(conc_fresh) | is.data.frame(xa_fresh)){
      conc_fresh <- variable_conc(conc_fresh, xa_fresh, t, t_run)$conc_fresh
      xa_fresh <- variable_conc(conc_fresh, xa_fresh, t, t_run)$xa_fresh
    }
    
    # Hard-wired temp settings settings
    temp_standard <- 298
    temp_zero <- 273
    
    #temp functions
    temp_air_C <- temp_air_C_fun(t + t_run)
    temp_C <- temp_C_fun(t + t_run)
    temp_K <- temp_C + 273.15
    temp_air_K <- temp_air_C + 273.15
    
    # Find methanogens and sulfate reducers
    i_meth <- which(grepl('^[mp]', names(qhat_opt)))
    i_sr <- which(grepl('^sr', names(qhat_opt)))
    n_mic <- length(qhat_opt)
    
    # Extract state variable values from y argument

    xa <- y[1:n_mic]
    y <- as.list(y[-c(1:n_mic)])
    list2env(y, envir = .GlobalEnv)
    
    # Hard-wired equilibrium constants
    log_ka <- c(NH3 = - 0.09046 - 2729.31/temp_K, 
                H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K),
                HAC = -4.8288 + 21.42/temp_K)
    
    kH_oxygen <- 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000
   
    # Hard-wire NH4+ activity coefficient
    g_NH4 <- 0.7
    
    # Hydrolysis rate with CTM
    alpha <- scale['alpha_opt'] * CTM(temp_K, alpha_T_opt, alpha_T_min, alpha_T_max, alpha_opt)

    # Microbial substrate utilization rate (vectorized calculation)
    qhat <- scale['qhat_opt'] * CTM(temp_K, T_opt, T_min, T_max, qhat_opt)
    
    # Ks temperature dependence
    ks <- ks_coefficient * (0.8157 * exp(-0.063 * temp_C)) 
    
    # NTS: Move this all out to a speciation function???
    # Chemical speciation (in rates() because is pH dependent)
    pH_surf <- pH + 1 # rough approximation from Rotz et al. IFSM 2012., "markfoged" thesis, "Petersen et al. 2014", "bilds?e et al. not published", "Elzing & Aarnik 1998", and own measurements..
    pH_surf_floor <- 8 # pH of manure on floor is not affected by acidification. Kept at 6.8
    
    HAC_frac <- 1-(1/(1 + 10^(-log_ka[['HAC']] - pH)))
    NH3_frac <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH)))) 
    NH3_frac_surf <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH_surf))))
    NH3_frac_floor <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH_surf_floor))))
    H2S_frac <- 1 - (1/(1 + 10^(- log_ka[['H2S']] - pH))) # H2S fraction of total sulfid
    # NTS: or just add H2CO3* here?
    # NTS: need TIC production too
    
    # HAC inhibition
    HAC_inhib <- ifelse(HAC_frac * VFA/slurry_mass >= 0.05, 1.16*0.31/(0.31 + (HAC_frac * VFA/slurry_mass)), 1) 
    
    # pH inhibition
    pH_inhib <- (1 + 2*10^(0.5*(pH_lwr - pH_upr)))/(1 + 10^(pH - pH_upr) + 10^(pH_lwr - pH))

    # NH3, NH4 inhibition
    NH3_inhib <- ifelse(NH3_frac * TAN/slurry_mass <= ki_NH3_min, 1, exp(-2.77259 * ((NH3_frac * (TAN/(slurry_mass)) - ki_NH3_min)/(ki_NH3_max - ki_NH3_min))^2))
    NH4_inhib <- ifelse((1 - NH3_frac) * (TAN/slurry_mass) <= ki_NH4_min, 1, exp(-2.77259*(((1 - NH3_frac) * (TAN/(slurry_mass)) - ki_NH4_min)/(ki_NH4_max - ki_NH4_min))^2))

    # H2S inhibition
    H2S_inhib_meth <-  1 - (H2S_frac*(sulfide/(slurry_mass))) / ki_H2S_meth # H2S _inhib factor for methanogens
    H2S_inhib_sr <-  1 - (H2S_frac*(sulfide/(slurry_mass))) / ki_H2S_sr     # H2S _inhib factor for sr
    if (H2S_frac*(sulfide/(slurry_mass)) > ki_H2S_meth) H2S_inhib_meth <-0
    if (H2S_frac*(sulfide/(slurry_mass))>ki_H2S_sr) H2S_inhib_sr <-0 
    
    cum_inhib_meth <- HAC_inhib * pH_inhib * NH3_inhib * NH4_inhib * H2S_inhib_meth
    cum_inhib_sr <- HAC_inhib * pH_inhib * NH3_inhib * NH4_inhib * H2S_inhib_sr
    
    # H2SO4 inhibition
    H2SO4_inhib <- H2SO4_inhibition_fun(sulfate/slurry_mass)
    
    if(any(H2SO4_inhib < cum_inhib_meth)) cum_inhib_meth[i_meth] <- H2SO4_inhib
    if(any(H2SO4_inhib < cum_inhib_sr)) cum_inhib_sr[i_sr] <- H2SO4_inhib
    
    # Henrys constant temp dependency
    H.NH3 <- 1431 * 1.053^(293 - temp_K)
    
    # Reduction from cover 
    kl[['NH3']] <- kl[['NH3']] * EF_NH3 
  
    # NH3 emission g(N) pr day
    NH3_emis_rate_floor <- kl[['NH3_floor']] * floor_area * ((NH3_frac_floor * TAN)/(slurry_mass)) * 1000 / H.NH3 # multiply by 1000 to get from g/kg to g/m3
    NH3_emis_rate_pit <- kl[['NH3']] * area * ((NH3_frac_surf * TAN)/(slurry_mass)) * 1000 / H.NH3 # multiply by 1000 to get from g/kg to g/m3
    
    H2S_emis_rate <- kl[['H2S']] * area * ((H2S_frac * sulfide)/(slurry_mass)) * 1000 # multiply by 1000 to get from g/kg to g/m3
    
    # Respiration gCOD/d, second-order reaction where kl applies for substrate concentration of 100 g COD / kg slurry
    
    kl.oxygen <- 0.0946 * temp_C ^ 1.2937 # from own lab experiments. 
    
    sub_respir <- VSd 
    if (sub_respir <= 0) sub_respir <- 1E-20
    respiration <- kl.oxygen * area * ((kH_oxygen * 0.208) - 0) * (sub_respir / slurry_mass) / 100 


    rut <- NA * qhat
    
    # VFA consumption rate by sulfate reducers (g/d). Slow speed
    rut[i_sr] <- ((qhat[i_sr] * VFA / (slurry_mass) * xa[i_sr] / (slurry_mass) / (scale['ks_coefficient'] * ks[i_sr] + VFA / (slurry_mass))) * (slurry_mass) *
                    (sulfate / (slurry_mass)) / (ks_SO4 + sulfate / (slurry_mass))) * cum_inhib_sr[i_sr]

    # Substrate utilization rate by methanogen groups in g/d affected by inhibition terms. Slow speed
    rut[i_meth] <- ((qhat[i_meth] * VFA / (slurry_mass) * xa[i_meth] / (slurry_mass)) / (scale['ks_coefficient'] * ks[i_meth] + VFA / (slurry_mass)) *
                      (slurry_mass)) * cum_inhib_meth[i_meth]
   
    # Some checks for safety
    if (any(rut < 0)) stop('In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)')

    # If there are no SRs...
    rutsr <- rut[i_sr]
    if (length(rutsr) == 0) rutsr <- 0

    # Derivatives, all in g/d except slurry_mass = kg/d
    # NTS: Some of these repeated calculations could be moved up
    derivatives <- c(
      xa = scale['yield'] * yield * rut + scale['xa_fresh'] * xa_fresh * slurry_prod_rate - decay_rate * xa, # expands to multiple elements with element for each mic group
      slurry_mass = slurry_prod_rate + (rain - evap) * area,
      xa_dead = slurry_prod_rate * conc_fresh[['xa_dead']] - alpha[['xa_dead']] * xa_dead + sum(decay_rate * xa),
      VSd = slurry_prod_rate * conc_fresh[['VSd']] - alpha[['VSd']] * VSd - respiration * VSd/sub_respir,
      VFA = alpha[['xa_dead']] * xa_dead + alpha[['VSd']] * VSd - sum(rut) + slurry_prod_rate * conc_fresh[['VFA']],
      urea = slurry_prod_rate * conc_fresh[['urea']] - alpha[['urea']] * urea,
      TAN = slurry_prod_rate * conc_fresh[['TAN']] + alpha[['urea']] * urea  - NH3_emis_rate_pit - NH3_emis_rate_floor,
      sulfate = slurry_prod_rate * conc_fresh[['SO4']] - sum(rutsr) * COD_conv[['S']],
      sulfide = slurry_prod_rate * conc_fresh[['S2']] + sum(rutsr) * COD_conv[['S']] - H2S_emis_rate,
      NH3_emis_cum = NH3_emis_rate_pit + NH3_emis_rate_floor,
      CH4_emis_cum = sum(rut[i_meth]) * COD_conv[['CH4']],
      CO2_emis_cum = sum(rut[i_meth]) * COD_conv[['CO2_anaer']] + sum(rutsr) * COD_conv[['CO2_sr']] + respiration * COD_conv[['CO2_aer']] + alpha[['urea']] * urea * COD_conv[['CO2_ureo']],
      COD_conv_cum = sum(rut[i_meth]) + respiration + sum(rutsr),
      COD_conv_cum_meth = sum(rut[i_meth]),
      COD_conv_cum_respir = respiration,
      COD_conv_cum_sr = rutsr
     )

    return(list(derivatives, c(H2S_emis_rate = H2S_emis_rate, NH3_emis_rate_pit = NH3_emis_rate_pit,
                               NH3_emis_rate_floor = NH3_emis_rate_floor, H2S_inhib_meth = H2S_inhib_meth, 
                               H2S_inhib_sr = H2S_inhib_sr, pH_inhib = pH_inhib, qhat = qhat, alpha = alpha, 
                               NH3_inhib = NH3_inhib, NH4_inhib = NH4_inhib, HAC_inhib = HAC_inhib, H2SO4_inhib = H2SO4_inhib, cum_inhib_meth = cum_inhib_meth, 
                               cum_inhib_sr = cum_inhib_sr, rut = rut, t_run = t_run, t_batch = t_batch, conc_fresh = conc_fresh, xa_init = xa_init, 
                               xa_fresh = xa_fresh * scale['xa_fresh'], area = area, t_batch = t_batch, slurry_prod_rate = slurry_prod_rate,
                               respiration = respiration, rain = rain, evap = evap)))
    
  }
