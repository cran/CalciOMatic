`fluo` <-
function(Ca=1,
                 R_min=0.136,
                 R_max=2.701,
                 K_eff=3.637,
                 K_d=0.58,
                 B_T=100,
                 phi=1.25,
                 S_B=10,
                 T_stim=0.015,
                 P=400,
                 P_B=400
                 ) {
  ## Function fluo
  ## Returns a vector of Fluorescence values according to Eq. 18, p8.
  ##
  ## Ca: a vector of [Ca] in muM
  ## R_min: R_min parameter from calibration experiments (dimensionless),
  ##        representing the minimum ratiometric measurement observable
  ##        [see Eq. 21, p10]
  ## R_max: R_max parameter from calibration experiments (dimensionless),
  ##        representing the maximum ratiometric measurement observable
  ##        [see Eq. 22, p10]
  ## K_eff: Effective fura dissociation constant from calibration experiments (in muM),
  ##      [see Eq. 25, p10]
  ## K_d: fura dissociation constant in muM (from calibration experiments),
  ##      [see Sec. 3.2]
  ## B_T: the total fura concentration (in muM)
  ## phi: dimensionless experiment specific parameter
  ## S_B: background (+ dark current) fluorescence intensity at a given wavelength,
  ##      in count/pixel/sec
  ## T_stim: exposure time at a given wavelength (usually 340 or 380 nm) (in s)
  ## P: number of pixels used for the fluorescence data binning
  ## P_B: number of pixels used for the background data binning
  
  pre_factor <- B_T * phi / (K_d+Ca)
  mid_factor <- R_min * K_eff + R_max*Ca
  post_factor <- P * T_stim

  fluo_bg <- S_B * T_stim * P_B
  
  ## Create the fluorescence transient
  result <- pre_factor * mid_factor * post_factor + fluo_bg
  
  return(result)
}

