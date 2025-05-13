#' Create data for stress-strain curve from tree winching data
#'
#' Function to calculate stress (xx) and strain (xx) information from tree winching
#' data as a pre-requisite for calculate modulus of elasticity.
#' @param moment_kNm numeric, vector of values representing turning moment in
#' kNm. Usually calcualted as applied force (kN) from tensiometer times height
#' of strap attachment (m) 
#' @param tilt_deg numeric, vector of values representing tilt measurements
#' from an angle gauge. Angles should be positive and  from vertical
#' @param ht_m numeric, height of tilt measurement
#' @param diam_cm numeric, diameter of tree at point of tilt measurement,
#' usually at the point of strap attachment
#' @param plot boolean, if `TRUE` plots a stress-strain diagram
#' @export
stress_strain = function(moment_kNm, tilt_deg, ht_m, diam_cm, plot=TRUE) {
  #ensure values are entered correctly
  stopifnot(is.numeric(moment_kNm))
  stopifnot(is.numeric(tilt_deg))
  stopifnot(is.numeric(ht_m))
  stopifnot(is.numeric(diam_cm))
  stopifnot(length(moment_kNm) == length(tilt_deg))
  stopifnot(length(ht_m)==1)
  stopifnot(length(diam_cm)==1)
  stopifnot(!is.null(plot))
  stopifnot(!is.na(plot))
  stopifnot(plot==TRUE|plot==FALSE)
  
  # calculate values
  c = (diam_cm/2) / 100 #distance from neutral axis (ie radius)
  I = pi/4 * (diam_cm/200)^2 # second moment of area for a circle
  stress = moment_kNm * c / I # in kN (Knm*m / m2)
  strain = ht_m * sin(pi*tilt_deg/180) #horizontal displacement of attachment point (strain)
  #organize and clean up data for s-s diagram
  df = data.frame(stress=stress, strain=strain)
  if(plot) {
    plot(stress~strain, df, type='l', lwd=1, col=grey(0.5),
         xlab = 'shear strain (m)',
         ylab = expression(shear~stress~(kN~m)))
    points(stress~strain, df, pch=16, cex=0.1)
  }
  return(df)
}

#' Calculate Modulus of elasticity from stress-strain curve
#'
#' Function to estimate modulus of elasticity from the linear portion of a 
#' stress-strain curve. Function looks to 75% of the values to the left of 
#' peak stress, and identifies the linear portion through iterative linear
#' modeling until the best fit is found
#' @param stress numeric, vector of values representing stress calculated from
#' function `stress_strain()`. 
#' @param strain numeric, vector of values representing strain calculated from
#' function `stress_strain()`. 
#' @param plot boolean, if `TRUE` plots a stress-strain diagram with fitted
#' line on linear portion representing modulus of elasticity
#' @param robustness, numeric between 0 and 1 representing the minimum proportion
#' of stress/strain readings to include in the linear fit. A value of 0.33 
#' represents that of values left of the peak, >=33% will be included in the model
#' fit.
#' @export
getMOE = function(stress, strain, plot=TRUE, robustness = 0.33, cleanup_width=5){
  stopifnot(robustness<1 & robustness > 0)
  stopifnot(is.numeric(cleanup_width))
  stopifnot(is.numeric(stress))
  stopifnot(is.numeric(strain))
  stopifnot(length(stress) == length(strain))
  stopifnot(!is.null(plot))
  stopifnot(!is.na(plot))
  stopifnot(plot==TRUE|plot==FALSE)
  
  # tidy input data
  df = data.frame(stress=stress, strain=strain)
  df = tidyr::drop_na(df)
  df = dplyr::filter(df, strain > 0 & stress > 0 )
  #fill in some missing values width = 5
  rap = c(zoo::rollapply(df$stress, width=cleanup_width, FUN = median), rep(NA,cleanup_width-1))
  head(rap)
  head(df$stress)
  
  # get MOE
  peak_stress = max(df$stress, na.rm=TRUE)
  peak_strain = df$strain[df$stress == peak_stress][1]
  
  # subset curve ahead of peak
  moe = dplyr::filter(df, strain < peak_strain)
  n = ceiling(nrow(df)*0.75)
  n_req = ceiling(n * robustness) # minimum number of points to include in curve
  out = list()
  
  # iteratively subset stress/strain data and find 
  for(i in n:n_req) {
    moe_tmp = dplyr::slice(moe, 1:i)
    mod = lm(stress ~ strain, data = moe_tmp)
    rmse = sqrt(mean(resid(mod)^2))
    out_tmp = data.frame(moe = coef(mod)[2],
                         intercept = coef(mod)[1],
                         rmse = rmse,
                         r2 = summary(mod)$r.squared)
    out[[length(out)+1]] = out_tmp
  }; out = do.call(rbind, out); 
  out = dplyr::filter(out, r2 == max(r2))
  rownames(out) = NULL
  
  # plot stress-strain + MOE
  if(plot) {
    plot(stress~strain, df, type='l', lwd=1, col=grey(0.5),
         xlab = 'shear strain (m)',
         ylab = expression(shear~stress~(kN~m)))
    points(stress~strain, df, pch=16, cex=0.1)
    abline(a = out$intercept, b = out$moe, col='red')
    r2 = round(out$r2,4)
    moe = round(out$moe,1)
    title(paste0('MOE: ', moe, '; R2: ', r2))
  }
  return(out)
}

