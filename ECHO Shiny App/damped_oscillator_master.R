# Extended Oscillations Function Source
# By Hannah De los Santos
# ECHO v 4.0
# Code description: Contains all the funcitons for extended harmonic oscillator work, in order to have less confusion between scripts.

# optimization support functions ----

# function to calculate optimization objective for weighted nls (with replicates)
# inputs:
#  p: paramters
#  t: time points
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
#  w: weights
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn <- function(p, t, s1, s2, w, y, fcall, jcall){
  sqrt(w)*(y - do.call("fcall", c(list(t = t, s1 =s1, s2 = s2), as.list(p))))
}

# function to calculate optimization objective for non-weighted nls (with no replicates)
# inputs:
#  p: paramters
#  t: time points
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn_one_rep <- function(p, t, s1, s2, y, fcall, jcall){
  (y - do.call("fcall", c(list(t = t, s1 =s1, s2 = s2), as.list(p))))
}

# function to calculate optimization objective for weighted nls (with replicates) for one omics type
# inputs:
#  p: paramters
#  t: time points
#  w: weights
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn_single <- function(p, t, w, y, fcall, jcall){
  sqrt(w)*(y - do.call("fcall", c(list(t = t), as.list(p))))
}

# function to calculate optimization objective for non-weighted nls (with no replicates) for one omics type
# inputs:
#  p: paramters
#  t: time points
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn_single_one_rep <- function(p, t, y, fcall, jcall){
  (y - do.call("fcall", c(list(t = t), as.list(p))))
}


# other ----

# function to represent damped oscillator with phase and equilibrium shift formula
# inputs:
#  a: Amplitude
#  gam: Amplitude.Change.Coefficient (amount of damping/driving)
#  omega: Radial frequency
#  phi: Phase Shift (radians)
#  y_shift: Equilibrium shift
#  t: time
# outputs:
#  result of inputs into formula
alt_form <- function(a,gam,omega,phi,y_shift,t){
  return (a*exp(-1*gam*(t)/2)*cos((omega*t)+phi)+y_shift)
}

# function to calculate the average of replicates
# inputs:
#  current_gene: row number of gene being examined
#  num_reps: number of replicates
# outputs:
#  a vector of means of gene expressions for replicates
avg_rep <- function(current_gene, num_reps){
  return(rbind(sapply(seq(2,ncol(genes), by = num_reps), function(x) mean(unlist(genes[current_gene,c(x:(num_reps-1+x))]), na.rm = TRUE))))
}

# function to calculate the average of all replicates. Used primarily for cases with multiple replicates.
# inputs:
#  num_reps: number of replicates
# outputs:
#  a matrix of means of gene expressions for replicates
avg_all_rep <- function(num_reps){
  # originally based on heat map code, but it will work fine here

  #get matrix of just the relative expression over time
  hm_mat <- as.matrix(genes[,2:ncol(genes)])

  #if there are replicates, average the relative expression for each replicate
  mtx_reps <- list() # to store actual replicate matrix
  mtx_count <- list() # to store how many are NA
  for (i in 1:num_reps){
    mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
    mtx_count[[i]] <- is.na(mtx_reps[[i]])
    mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
  }
  repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
  hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
  for (i in 1:num_reps){
    hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
    repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
  }
  repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
  hm_mat <- hm_mat/repmtx

  return(hm_mat)
}

# function to calculate the variance of replicates at a certain time point
# inputs:
#  x: time point
#  current_gene: row number of gene being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_var <- function(x, current_gene, num_reps){
  eps <- 1e-7 # slight adjustment for 0 variance case
  std2 <- var(unlist(genes[current_gene,c(x:(num_reps-1+x))]), na.rm = TRUE) # calc variance
  if (is.na(std2)){ # possible bug for NA
    std2 <- 1
  }
  return(1/(std2+eps))
}

# function to calculate the weights for replicate fitting (these weights are the inverse of the variance at each time point)
# inputs:
#  current_gene: row number of gene being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_weights <- function(current_gene, num_reps){
  return(sapply(seq(2,ncol(genes), by = num_reps), function(x) calc_var(x,current_gene,num_reps)))
}

# function to find confidence intervals by the standard jackknifing method (leave one out)
# inputs:
#  temp: data frame with time points, y values, and weights (if used)
#  parameters: parameters found by full nlsLM fit
#  num_reps: number of replicates
#  start_param: list of starting parameters for fits
# outputs:
#  confidence intervals for each parameter, low and then high
jackknife <- function(temp, parameters, num_reps, start_param){
  # preallocate spaces
  all_fits <- list()
  all_pars <- vector(mode = "list", 5)
  n <- nrow(temp)

  # go through and compute fit for each leave one out
  for (i in 1:nrow(temp)){
    temp_edit <- temp[-i,]

    if (num_reps==1){
      all_fits[[i]] <- suppressWarnings(
        nls.lm(par = start_param,
               fn = fcn_single_one_rep,
               fcall = alt_form,
               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
               upper=c(Inf, Inf, low, Inf, max(temp$y)),
               t = temp_edit$t,
               y = temp_edit$y,
               control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)
        ))
    } else {
      #fitting
      all_fits[[i]] <- suppressWarnings(
        nls.lm(par = start_param,
               fn = fcn_single,
               fcall = alt_form,
               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
               upper=c(Inf, Inf, low, Inf, max(temp$y)),
               t = temp_edit$t,
               y = temp_edit$y,
               w = temp_edit$w,
               control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)
        ))
    }

    # store all parameters
    for (p in 1:5){
      all_pars[[p]] <- c(all_pars[[p]],coef(all_fits[[i]])[p])
    }
  }

  # compute jackknifed confidence intervals
  par_names <- c("gam","a","omega","phi","y_shift")
  ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
  ci_int <- 1:(2*5)
  for (p in 1:5){
    ps <- (n*parameters[p])-((n-1)*all_pars[[p]])
    v <- 1/n/(n-1)*sum((ps-(sum(ps)/n))^2)
    ci_int[p] <- mean(ps) - qt(0.975, n-1)*sqrt(var(ps)/n)
    ci_int[p+5] <- mean(ps) + qt(0.975, n-1)*sqrt(var(ps)/n)
  }
  names(ci_int) <- ci_names

  return(ci_int)
}


# function to calculate the variance of replicates at a certain time point for bootstrapping
# inputs:
#  x: time point
#  boot_gene: gene time course being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_var_boot <- function(x, boot_gene, num_reps){
  eps <- 1e-7 # slight adjustment for 0 variance case
  std2 <- var(unlist(boot_gene[c(x:(num_reps+x))]), na.rm = TRUE) # calc variance
  if (is.na(std2)){ # possible bug for NA
    std2 <- 1
  }
  return(1/(std2+eps))
}

# function to calculate the weights for replicate fitting for bootstrapping (these weights are the inverse of the variance at each time point)
# inputs:
#  boot_gene: gene time course being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_weights_boot <- function(boot_gene, num_reps){
  return(sapply(seq(1,length(boot_gene), by = num_reps), function(x) calc_var_boot(x,boot_gene,num_reps)))
}

# function to find confidence intervals by the standard jackknifing method (leave one out)
# inputs:
#  temp: data frame with time points, y values, and weights (if used)
#  fit: original fit data
#  parameters: parameters found by full nlsLM fit
#  num_reps: number of replicates
#  start_param: list of starting parameters for fits
#  current_gene: row number of gene being examined
#  seed: random seed to set for reproducibility
# outputs:
#  confidence intervals for each parameter, low and then high
bootstrap <- function(temp, fit, start_param, num_reps, current_gene, seed){
  df <- temp
  df$fitted <- as.numeric(fitted(fit))
  df$resid <- as.numeric(temp$y) - as.numeric(df$fitted)
  fun <- function(df, inds) {
    # get the new bootstrapped values
    df$bootGPP <- df$fitted + df$resid[inds]
    # remove trycatch -- no errors
    if (num_reps > 1){ # with replicates
      # get new weights based on these
      # insert nas into bootcpp, to stagger weights appropriately
      boot_var <- as.numeric(genes[current_gene,-1])
      boot_var[!is.na(boot_var)] <- df$bootGPP
      w <- rep(calc_weights_boot(boot_var, num_reps), each = num_reps)
      df$w <- w[!is.na(genes[current_gene,-1])]
      
      suppressWarnings(coef(
        nls.lm(par = start_param,
               fn = fcn_single,
               fcall = alt_form,
               lower=c(-Inf, -Inf, high, -Inf, min(df$y)),
               upper=c(Inf, Inf, low, Inf, max(df$y)),
               t = df$t,
               y = df$bootGPP,
               w = df$w,
               control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)
        )))
    } else { # one replicate
      suppressWarnings(coef(nls.lm(par = start_param,
                                   fn = fcn_single_one_rep,
                                   fcall = alt_form,
                                   lower=c(-Inf, -Inf, high, -Inf, min(df$y)),
                                   upper=c(Inf, Inf, low, Inf, max(df$y)),
                                   t = df$t,
                                   y = df$bootGPP,
                                   control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                            ftol=1e-6, ptol=1e-6, gtol=1e-6)
      )))
    }
  }

  set.seed(seed)
  b <- boot(df, fun, R = 999)
  # get the 95% confidence intervals
  par_names <- c("gam","a","omega","phi","y_shift")
  ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
  ci_int <- 1:(2*5)
  for (p in 1:5){
    bci <- boot.ci(b, index=p, type = "perc")

    if (!is.null(bci)){
      ci_int[p] <- bci$percent[4]
      ci_int[p+5] <- bci$percent[5]
    } else { # if there's no variation in the parameter - fixed interval
      ci_int[p] <- ci_int[p+5] <- b$t0[3]
    }
  }
  names(ci_int) <- ci_names


  return(ci_int)
}

# Function to calculate the parameters for the extended harmonic oscillator equation for a specific gene.
#
#  current_gene: row number of current gene we want to calculate parameters for
#  timen: time points for dataset
#  resol: resolution of time points
#  num_reps: number of replicates
#  tied: if replicate data, whether the replicates are related (paired) or not (unpaired)
#  is_smooth: boolean that indicates whether data should be smoothed or not
#  is_weighted: if there is smoothing, is it weighted (1,2,1) smoothing, or unweighted smoothing (1,1,1)
#  low: the highest frequency we are looking for, in radians (lowest period)
#  high: the lowest frequency we are looking for, in radians (highest period)
#  rem_unexpr: logical indicating whether genes with less than rem_unexpr_amt % expression should not be considered
#  rem_unexpr_amt: percentage of expression for which genes should not be considered
#  run_conf: boolean of whether or not to run confidence intervals
#  which_conf: string of which type of confidence interval to compute
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#  seed: number for random seed to fix for bootstrapping for confidence intervals
#  results, a data frame which contains:
#   gene: gene name
#   conv: did the fit converge, or descriptor of type of data (constant, unexpressed, etc.)
#   iter: number of iterations
#   gamma: forcing coefficient value for fit
#   type_gam: Type of oscillation (damped, forced, etc.)
#   amplitude: Amplitude value for fit
#   omega: Radial frequency for fit
#   period: Period for fit (in time units)
#   phase.shift: Phase shift for fit (radians)
#   hours.shift: Phase shift for fit (hours)
#   tau: Kendall's tau between original and fitted values
#   y_shift: Equilibrium shift for fit
#   pval: P-value calculated based on Kendall's tau
#   (ci_int: confidence for each parameter)
#   original.values: original values for gene
#   fitted.values: fitted values for gene
# calculate_param <- function(current_gene,timen,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, which_conf = "Bootstrap", harm_cut = .03, over_cut = .15, seed = 30)

# function to smooth, with weighting, a given averaged expression
# inputs:
#  y_val: vector of averaged expressions
# outputs:
#  smoothed y_val
smooth_y_val <- function(y_val){
  # smoothing the starting averages - weighted or unweighted?
  #get matrix of just the relative expression over time
  all_reps <- matrix(c(y_val), nrow = 1, ncol = length(y_val))
  
  # weighted averaging
  
  #if there are replicates, average the relative expression for each replicate
  center_reps <- list() # to store actual replicate matrix
  mtx_count <- list() # to store how many are NA
  for (i in 1:1){
    center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=1)]*2
    mtx_count[[i]] <- is.na(center_reps[[i]])*2
    center_reps[[i]][is.na(center_reps[[i]])] <- 0
  }
  
  repmtx_l <- list() # store amount to divide by for each rep
  repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
  repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges
  
  # sum the replicates
  left <- c(rep(0,1),center_reps[[i]][-length(center_reps[[i]])]/2) # left shifted matrix
  right <- c(center_reps[[i]][-1]/2,rep(0,1)) # right shifted matrix
  center_reps[[i]] <- left + center_reps[[i]] + right
  
  # figure out how many replicates are actually available for each time point
  left_na <- c(rep(0,1),mtx_count[[i]][-length(mtx_count[[i]])]/2) # left shifted matrix
  right_na <- c(mtx_count[[i]][-1]/2,rep(0,1)) # right shifted matrix
  repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
  # to avoid division by 0 and induce NAs if there are no time points available
  repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  
  dat <- all_reps
  x <- 0
  # averaging to get smoothed replicates
  dat[,seq(1+x,ncol(all_reps),by=1)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  dat[is.na(all_reps)] <- NA # do not impute missing values
  
  return(dat[1,])
}

# function to calculate the "peaks" vector for starting points for echo-type models, "peaks" can be peaks or troughs
# inputs:
#  y_val: vector of averaged expressions
#  resol: resolution
#  y0: initial equilibrium shift
#  timen: numeric vector of time points
# outputs:
#  list with numeric vector of peaks, and the time points for those peaks
calc_peaks <- function(y_val, resol, y0, timen){
  # figure out the resolution modifier (must have at least 24 hour data)
  mod <- 0
  if (resol <= ((1/12)+10^-8)){ # 17 hour surround
    mod <- 102
  } else if (resol <= ((1/6)+10^-8)){ # 15 hour surround
    mod <- 45
  } else if (resol <= ((1/4)+10^-8)){ # 13 hour surround
    mod <- 26
  } else if (resol <= ((1/2)+10^-8)){ # 11 hour surround
    mod <- 11
  } else if (resol <= 1){ # 8 hour surround
    mod <- 4
  } else if (resol <=2){ # 6 hour surround
    mod <- 3
  } else if (resol <= 4){ # 4 hour surround
    mod <- 2
  } else{
    mod <- 1
  }
  
  # calculate the amount of peaks
  peaks <- c(); # vector of peak values
  peaks_time <- c(); # vector of peak timen
  counting <- 1; # counter
  # go through all time points and calculate peaks (max values)
  for(i in (mod+1):(length(y_val)-mod)){
    # deal with complete missingness
    if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
      next
    }
    # otherwise continue as normal
    if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
      peaks[counting] <- y_val[i]
      peaks_time[counting] <- timen[i]
      counting <- counting+1
    }
  }
  
  # calculate amount of troughs
  troughs <- c() # vector of trough values
  troughs_time <- c() # vector of trough timen
  counting <- 1 # counter
  # go through all time points and calculate peaks (max values)
  for(i in (mod+1):(length(y_val)-mod)){
    # deal with complete missingness
    if (suppressWarnings(min(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
      next
    }
    # otherwise continue as normal
    if (y_val[i] == min(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
      troughs[counting] <- y_val[i]
      troughs_time[counting] <- timen[i]
      counting <- counting+1
    }
  }
  
  peaks_per_diff <- 0
  troughs_per_diff <- 0
  if (length(peaks)>1 & length(troughs)>1){
    # new criterion: which peaks are more evenly distributed?
    # what should the period be, based on the amount of peaks/troughs seen in the time course
    peaks_per <- length(timen)*resol/length(peaks)
    troughs_per <- length(timen)*resol/length(troughs)
    
    # get avg period, calculated as the difference in time between successive peaks
    peaks_avg_per <- mean(sapply(seq(length(peaks_time),2,-1),
                                 function(x){peaks_time[x]-peaks_time[x-1]}))
    troughs_avg_per <- mean(sapply(seq(length(troughs_time),2,-1),
                                   function(x){troughs_time[x]-troughs_time[x-1]}))
    
    # calculate the difference between the theoretical and estimated period
    peaks_per_diff <- abs(peaks_per-peaks_avg_per)
    troughs_per_diff <- abs(troughs_per-troughs_avg_per)
  }
  
  # if there is only one peak, or if our troughs are more consistent, we choose troughs
  if ((length(peaks)<=1 & length(troughs)>length(peaks)) |
      (troughs_per_diff < peaks_per_diff)){
    # flip the troughs for accurate calculations later
    # absolute value the difference from the midline
    peaks <- abs(y0 - troughs)
    peaks_time <- troughs_time
  } else {# else we choose peaks
    # absolute value the difference from the midline
    peaks <- abs(peaks - y0)
  }
  # it's possible for this not to work if you don't have more than one
  # oscillation over your time course
  # }
  
  return(list("peaks"=peaks, "peaks_time"=peaks_time))
  
}

# function to calculate the starting points for echo models
# inputs:
#  y_val: vector of averaged expressions
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  list with all starting points for parameters (amplitude, ac coeff, frequency, phase shift, equilibrium shift)
calc_start_echo <- function(y_val, over_cut){
  timen <- timen
  # smooth the expression
  y_val <- smooth_y_val(y_val)
  
  # calculate starting y_shift
  y0 <- mean(y_val,na.rm = TRUE) # intial value for the equilibrium shift
  if (y0 < 10^-10 && y0 > -10^-10){
    y0 <- 10^-8 # avoiding 0 mean, which makes gradient singular
  }
  
  # calculate peaks
  peaks_list <- calc_peaks(y_val, resol, y0, timen)
  peaks <- peaks_list$peaks
  peaks_time <- peaks_list$peaks_time
  
  # calc starting amplitude
  if (length(peaks) >0){
    a0 <- abs(peaks[1]) # y0 already removed
  } else {
    a0 <- max(y_val,na.rm = TRUE) - y0
  }
  
  # intial value for gamma
  if (length(peaks)==0){ # if there are no peaks, we account for that
    gam0 <- 0
  } else if (which.max(peaks)==1){ # if the highest peak is first, then damping is likely
    if (length(peaks)>1){
      # trying to figure out gamma based on logarithmic decrement
      n <- 1
      # n <- peaks_time[2]-peaks_time[1]
      log_dec <- (1/n)*log(abs(peaks[1]/peaks[2]))
      gam0 <- 1/(sqrt(1+((2*pi/log_dec)^2)))/2
    } else{ # use a guess
      gam0 <- .01
    }
  } else{ # otherwise forcing is likely
    if (length(peaks)>1){
      # trying to figure out gamma based on logarithmic decrement
      n <- 1
      # n <- peaks_time[2]-peaks_time[1]
      log_dec <- (1/n)*log(abs(peaks[2]/peaks[1]))
      gam0 <- -1*1/(sqrt(1+((2*pi/log_dec)^2)))/2
    } else{ # use a guess
      gam0 <- -.01
    }
  }
  
  # restricting ac coeff to not be overexpressed/repressed
  if (gam0 < -over_cut){
    gam0 <- -over_cut
  } else if (gam0 > over_cut){
    gam0 <- over_cut
  }
  
  # let frequency depend on amount of peaks = (length(timen)*resol/(no of peaks+1 [accounts for phase shift])
  if (length(peaks) == 0){
    if (high == -Inf || low == Inf){
      w0 <- 2*pi/(length(timen)*resol/2)
    } else{
      # want to get their actual integer period values
      highfix <- (high/2/pi)^-1
      lowfix <- (low/2/pi)^-1
      w0 <- 2*pi/(length(timen)*resol/((highfix+lowfix)/2))
    }
  } else if (length(peaks) == 1){ # phase shift causes only one peak to appear
    w0 <- 2*pi/(length(timen)*resol/(length(peaks)+1))
  } else{
    w0 <- 2*pi/(length(timen)*resol/(length(peaks)))
  }
  
  # can't be outside the specified parameters
  if (w0 > low){
    w0 <- low
  } else if (w0 < high){
    w0 <- high
  }
  
  
  # we estimate our phase shift on the second and third nonmissing sets of points for accuracy
  # if you have less than 3 points nonmissing, I have no hope for you
  
  # higher fidelity guess for higher amounts of points
  # this still works for 3 points,
  beg1 <- which(!is.na(y_val))[2]
  beg2 <- which(!is.na(y_val))[3]
  # middle
  mid1 <- which(!is.na(y_val))[floor(length(timen)/2)]
  mid2 <- which(!is.na(y_val))[floor(length(timen)/2)+1]
  # end
  en1 <- which(!is.na(y_val))[length(timen)-2]
  en2 <- which(!is.na(y_val))[length(timen)-1]
  
  min_vect <- rep(Inf, length(0:11))
  for (i in 0:11){
    # if the phase shift at the beginning, middle, and end are the smallest value available for the fitted value
    min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg1])-y_val[beg1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid1])-y_val[mid1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en1])-y_val[en1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg2])-y_val[beg2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid2])-y_val[mid2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en2])-y_val[en2]))
  }
  phi0 <- (which.min(abs(min_vect))-1)*pi/6 # intial value for phase shift
  
  # form all starting parameters into a nice list
  start_param <- list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0)
  
  return(start_param)
}

# function to categorize ac coefficient into categories
# inputs:
#  gam: ac coefficient
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  string, ac coefficient categories
get_type_gam <- function(gam, harm_cut, over_cut){
  # categorize gamma in ac coeff categories
  if (gam < -over_cut){
    type_gam <- "Overexpressed"
  } else if (gam <= -harm_cut){
    type_gam <- "Forced"
  } else if (gam <= harm_cut){
    type_gam <- "Harmonic"
  } else if (gam <= over_cut){
    type_gam <- "Damped"
  } else{
    type_gam <- "Repressed"
  }
  return(type_gam)
}

# function to calculate hours shifted (phase shift, in hours) for ECHO-type rhythms
# inputs:
#  a: initial amplitude
#  phi: phase shift, radians
#  omega: radian frequency
# outputs:
#  numeric, phase shift, in hours, relative to the 0 specified in the time course
get_hs <- function(a,phi,omega){
  # calculating the phase shift in terms of period (omega inverse of period)
  if (!is.na(a)){ # all param will either be na or not
    if (a >= 0){ # positive amplitudes
      if (phi > 0){ # shift to the left
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        phase_hours <- (2*pi/omega)-dist_peak
      } else { # shift to the right
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        phase_hours <- abs(dist_peak)
      }
    } else { # negative ampltitudes
      if (phi > 0){ # shift to the left
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        if (abs(frac_part) < .5){
          phase_hours <- (2*pi/omega)-dist_peak - (2*pi/omega/2)
        } else {
          phase_hours <- (2*pi/omega)-dist_peak + (2*pi/omega/2)
        }
      } else { # shift to the right
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        if (abs(frac_part) < .5){
          phase_hours <- abs(dist_peak) + (2*pi/omega/2)
        } else {
          phase_hours <- abs(dist_peak) - (2*pi/omega/2)
        }
      }
    }
  } else {
    phase_hours <- NA
  }
  
  return(phase_hours)
}


# Function to calculate the parameters for the extended harmonic oscillator equation for a specific gene.
#
#  current_gene: row number of current gene we want to calculate parameters for
#  timen: time points for dataset
#  resol: resolution of time points
#  num_reps: number of replicates
#  tied: if replicate data, whether the replicates are related (paired) or not (unpaired)
#  is_smooth: boolean that indicates whether data should be smoothed or not
#  is_weighted: if there is smoothing, is it weighted (1,2,1) smoothing, or unweighted smoothing (1,1,1)
#  low: the highest frequency we are looking for, in radians (lowest period)
#  high: the lowest frequency we are looking for, in radians (highest period)
#  rem_unexpr: logical indicating whether genes with less than rem_unexpr_amt % expression should not be considered
#  rem_unexpr_amt: percentage of expression for which genes should not be considered
#  run_conf: boolean of whether or not to run confidence intervals
#  which_conf: string of which type of confidence interval to compute
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#  seed: number for random seed to fix for bootstrapping for confidence intervals
#  results, a data frame which contains:
#   gene: gene name
#   conv: did the fit converge, or descriptor of type of data (constant, unexpressed, etc.)
#   iter: number of iterations
#   gamma: forcing coefficient value for fit
#   type_gam: Type of oscillation (damped, forced, etc.)
#   amplitude: Amplitude value for fit
#   omega: Radial frequency for fit
#   period: Period for fit (in time units)
#   phase.shift: Phase shift for fit (radians)
#   hours.shift: Phase shift for fit (hours)
#   tau: Kendall's tau between original and fitted values
#   y_shift: Equilibrium shift for fit
#   pval: P-value calculated based on Kendall's tau
#   (ci_int: confidence for each parameter)
#   original.values: original values for gene
#   fitted.values: fitted values for gene
calculate_param <- function(current_gene,timen,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, which_conf = "Bootstrap", harm_cut = .03, over_cut = .15, seed = 30){

  if (run_conf){
    ci_int <- rep(NA, 10)
    par_names <- c("gam","a","omega","phi","y_shift")
    ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
    names(ci_int) <- ci_names
  }

  gene_n <- as.character(genes[current_gene,1]) # gene name
  # first we need to check whether or not the gene is just a straight line
  if (!is_deviating(current_gene)){ # one replicate
    results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)

    if (run_conf){
      results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    } else {
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    }
    return(results)
  }

  # then we need to check if < threshold % are expressed (if desired)
  if (rem_unexpr){
    if (rem_unexpr_vect[current_gene]){
      results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      if (run_conf){
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      }
      return(results)
    }
  }

  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- as.numeric(genes[current_gene,-1]) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  # REMOVING THE TRYCATCH -- THERE SHOULD NO LONGER BE ERRORS
  # tryCatch({ # throw exception upon error
    
  # calculate starting values
  start_param <- calc_start_echo(y_val, over_cut)
  
  temp <- data.frame()
  if (num_reps == 1){ # one replicate
    # put the timen and data points into a data frame
    temp <- data.frame(y=y_val,t=timen)
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    # fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single_one_rep,
                              fcall = alt_form,
                              lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                              upper=c(Inf, Inf, low, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    # get coefficients
    coeff <- coef(oscillator_init)
    
    # put coefficient estimates into an nls object
    oscillator.fit <- nls_edit(
      formula = y ~ alt_form(a, gam, omega, phi, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], y_shift=coeff[5]),
      control = nls.control(maxiter = 0)
    )
    # end
  } else{ # multiple replicates
    #put the timen and data points into a data frame, and weights
    weights <- calc_weights(current_gene,num_reps)
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    # fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single,
                              fcall = alt_form,
                              lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                              upper=c(Inf, Inf, low, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    # get coefficients
    coeff <- coef(oscillator_init)
    
    # put coefficient estimates into an nls object
    oscillator.fit <- nls_edit(
      formula = y ~ alt_form(a, gam, omega, phi, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], y_shift=coeff[5]),
      control = nls.control(maxiter = 0)
    )
  }

  did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
  num_iter <- oscillator.fit$convInfo$finIter # amount of iterations

  parameters <- coef(oscillator.fit) #extract parameter estimates

  # alt_form parameters:
  # the parameters go in the order of: gam,a,omega,phi,y_shift
  gam <- parameters[1]
  a <- parameters[2]
  omega <- parameters[3]
  phi <- parameters[4]
  y_shift <- parameters[5]

  if (run_conf){
    if (which_conf == "Bootstrap"){
      ci_int <- bootstrap(temp, oscillator.fit, start_param, num_reps, current_gene, seed)
    } else {
      ci_int <- jackknife(temp, parameters, num_reps, start_param)
    }
  }

  # calculating whether (over)damped, (over)forced, harmonic
  type_gam <- get_type_gam(gam, harm_cut, over_cut)
  
  # get hours shifted
  phase_hours <- get_hs(a,phi,omega)

  # calculate p-value
  ref_wave <- (alt_form(a,gam,omega,phi,y_shift,timen)) # fitted values
  all_pred <- rep(ref_wave, each=num_reps)[!is.na(unlist(genes[current_gene,-1]))]
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value
  tau <- testing$estimate

  # list of parameters and other resulting values
  results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
  if (!run_conf){
    if (num_reps == 1){
      results <- cbind(results, rbind(y_val), rbind(ref_wave))
    } else {
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
    }
  } else {
    if (num_reps == 1){
      results <- cbind(results, rbind(ci_int), rbind(y_val), rbind(ref_wave))
    } else {
      results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
    }
  }
  return (results)
  # }, error = function(e){ # if there's failure in convergence
  #   results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
  #   if (!run_conf){
  #       results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
  #   } else {
  #     results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
  #   }
  # 
  #   return (results)
  # })
}

#UNNEEDED
# function for determining whether gene values are deviating or constant for full matrix
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean vector if there is deviation within the gene values
is_deviating_all <- function(){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])

  # vector of standard deviations
  all_row_stdev <- sqrt(rowSums((all_reps - rowMeans(all_reps,na.rm = TRUE))^2, na.rm = TRUE)/(dim(all_reps)[2] - 1))

  # return vector if there is deviation within the gene values
  return (all_row_stdev != 0)
}

# function for determining whether gene values are unexpressed (less than 70% expressed (i.e., below specified cutoff)) for full matrix
# inputs:
#  rem_unexpr_amt: what percentage of the time course must be expressed
#  rem_unexpr_amt_below: cutoff for expression
#  boolean if there is 70% expression
genes_unexpressed_all <- function(rem_unexpr_amt, rem_unexpr_amt_below){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  # get how many genes are expressed for each gene
  tot_expressed <- rowSums(abs(all_reps) > abs(rem_unexpr_amt_below),na.rm = TRUE)

  # return true if amount is less than threshold
  return(tot_expressed <= (ncol(all_reps)*rem_unexpr_amt))
}


# function for determining whether gene values are deviating or constant (one replicate)
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean if there is deviation within the gene values
is_deviating <- function(current_gene){
  if (all(is.na(genes[current_gene,-1])) | sum(!is.na(genes[current_gene,-1]))==1){
    stdev <- 0
  } else {
    stdev <- sd(genes[current_gene,-1], na.rm = TRUE)
  }
  return (stdev != 0)
}

# DEPRECIATED
# function for determining whether gene values are deviating or constant (multiple replicates)
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
#  num_reps: number of replicates
# outputs:
#  boolean if there is deviation within the gene values
is_deviating_rep <- function(current_gene,num_reps){
  y_val <- avg_rep(current_gene,num_reps)
  stdev <- sd(y_val,na.rm = TRUE)
  return (stdev != 0)
}

# function for determining whether gene values are unexpressed (less than 70% expressed (i.e., not 0))
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean if there is 70% expression
genes_unexpressed <- function(current_gene,rem_unexpr_amt){
  y_val <- genes[current_gene,!as.logical(is.na(genes[current_gene,]))]
  y_val <- y_val[,-1]
  tot_expressed <- sum(y_val != 0,na.rm = TRUE)
  return(tot_expressed <= (length(y_val)*rem_unexpr_amt))
}

# function to adjust pvals according to the benjamini-hochberg criterion
# inputs:
#  pvals: vector of pvalues
# outputs:
#  BH adjusted pvalues
adjust_p_values <- function(pvals){
  return (p.adjust(unlist(pvals), method = "BH"))
}

# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for paired
# (tied) replicates
# inputs:
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_all_tied <- function(is_weighted, num_reps){
  # originally based on heat map code, but it will work fine here

  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])

  if (!is_weighted){

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are NA
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]
      mtx_count[[i]] <- is.na(center_reps[[i]])
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+3 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 2 at edges
  } else{ # weighted averaging

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are NA
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]*2
      mtx_count[[i]] <- is.na(center_reps[[i]])*2
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges

  }
  for (i in 1:num_reps){
    # sum the replicates
    left <- cbind(matrix(0,nrow(all_reps),1),center_reps[[i]][,-ncol(center_reps[[i]])]/2) # left shifted matrix
    right <- cbind(center_reps[[i]][,-1]/2,matrix(0,nrow(genes),1)) # right shifted matrix
    center_reps[[i]] <- left + center_reps[[i]] + right

    # figure out how many replicates are actually available for each time point
    left_na <- cbind(matrix(0,nrow(all_reps),1),mtx_count[[i]][,-ncol(mtx_count[[i]])]/2) # left shifted matrix
    right_na <- cbind(mtx_count[[i]][,-1]/2,matrix(0,nrow(genes),1)) # right shifted matrix
    repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
    # to avoid division by 0 and induce NAs if there are no time points available
    repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  }

  dat <- genes
  for (x in 0:(num_reps-1)){
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  }
  dat[is.na(genes)] <- NA # do not impute missing values
  return(dat)
}

# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for paired
# (tied) replicates
# inputs:
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_all_untied <- function(is_weighted, num_reps){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])

  if (!is_weighted){

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are

    side_reps <- avg_genes
    mtx_side_count <- is.na(side_reps)
    side_reps[is.na(side_reps)] <- 0
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]
      mtx_count[[i]] <- is.na(center_reps[[i]])
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+3 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges
  } else{ # weighted averaging

    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are

    side_reps <- avg_genes
    mtx_side_count <- is.na(side_reps)
    side_reps[is.na(side_reps)] <- 0

    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]*2
      mtx_count[[i]] <- is.na(center_reps[[i]])*2
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }

    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges

  }
  for (i in 1:num_reps){
    # sum the replicates
    left <- cbind(matrix(0,nrow(genes),1),side_reps[,-ncol(side_reps)]) # left shifted matrix
    right <- cbind(side_reps[,-1],matrix(0,nrow(genes),1)) # right shifted matrix
    center_reps[[i]] <- left + center_reps[[i]] + right

    # figure out how many replicates are actually available for each time point
    left_na <- cbind(matrix(0,nrow(all_reps),1),mtx_side_count[,-ncol(mtx_side_count)]) # left shifted matrix
    right_na <- cbind(mtx_side_count[,-1],matrix(0,nrow(genes),1)) # right shifted matrix
    repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
    # to avoid division by 0 and induce NAs if there are no time points available
    repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  }

  dat <- genes # assigning to dataframe to return
  for (x in 0:(num_reps-1)){
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  }
  dat[is.na(genes)] <- NA # do not impute missing values
  return(dat)
}

# function to normalize expressions in a matricized manner, by row
# thanks to: https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
# outputs:
#   normalized expression data
normalize_all <- function(){

  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])

  # vector of row means
  all_row_mean <- rowMeans(all_reps, na.rm = TRUE)
  # vector of standard deviations
  all_row_stdev <- sqrt(rowSums((all_reps - rowMeans(all_reps,na.rm = TRUE))^2, na.rm = TRUE)/(dim(all_reps)[2] - 1))

  # get the matrix of normalized expressions
  all_reps_normal <- (all_reps - all_row_mean)/all_row_stdev
  # if standard deviation is 0, imposes NA, so this expression shouldn't be considered anyway
  # and is now constant
  all_reps_normal[is.na(all_reps_normal)] <- 0

  # create dataframe with normalized expressions
  dat <- genes
  dat[,-1] <- all_reps_normal
  dat[is.na(genes)] <- NA # do not impute missing values

  return(list("dat"=dat, "means"=all_row_mean, "stdevs"=all_row_stdev))
}

# DEPRECIATED
# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for paired
# (tied) replicates
# smooths each replicate separately
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_tied <- function(current_gene, is_weighted, num_reps){
  ldat <- lapply(c(0:(num_reps - 1)), function (x) { # return each replicate in a list of smoothed data
    if (is_weighted){ # (1,2,1) weighted average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(genes[current_gene,seq(z-num_reps,z+num_reps,by = num_reps)],c(1,2,1),na.rm=TRUE))
      first <- weighted.mean(genes[current_gene,c(2+x,2+num_reps+x)],c(2,1),na.rm=TRUE)
      last <- weighted.mean(genes[current_gene,c(ncol(genes)-(2*num_reps)+1+x,ncol(genes)-num_reps+x+1)],c(1,2),na.rm=TRUE)
    }
    else{ # (1,1,1) average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(genes[current_gene,seq(z-num_reps,z+num_reps,by = num_reps)],c(1,1,1),na.rm=TRUE))
      first <- weighted.mean(genes[current_gene,c(2+x,2+num_reps+x)],c(1,1),na.rm=TRUE)
      last <- weighted.mean(genes[current_gene,c(ncol(genes)-(2*num_reps)+1+x,ncol(genes)-num_reps+x+1)],c(1,1),na.rm=TRUE)
    }
    return (c(first,center.dat,last))}
  )
  dat <- genes[current_gene,] # original data
  for (x in 0:(num_reps-1)){ # replace original data by smoothed data, keeping na values
    ldat[[x+1]][is.na(genes[current_gene,seq(2+x,ncol(genes),by=num_reps)])] <- NA
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- ldat[[x+1]]
  }

  return(dat)
}



# DEPRECIATED
# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for unpaired
# (untied) replicates
# smooths each replicate with the following scheme: left average, data of specific replicate, right average
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_untied <- function(current_gene, is_weighted, num_reps){
  y_val <- avg_rep(current_gene,num_reps) # calculating the average of the technical reps
  y_val <- c(genes[current_gene,1],y_val) # adding the gene name for similarity of index

  ldat <- lapply(c(0:(num_reps - 1)), function (x) { # return each replicate in a list of smoothed data
    if (is_weighted){ # (1,2,1) weighted average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(c(as.numeric(y_val[z-num_reps]),genes[current_gene,z],as.numeric(y_val[z+num_reps])),c(1,2,1),na.rm=TRUE))
      first <- weighted.mean(c(genes[current_gene,2+x],as.numeric(y_val[2+num_reps+x])),c(2,1),na.rm=TRUE)
      last <- weighted.mean(c(as.numeric(y_val[ncol(genes)-(2*num_reps)+1+x]),genes[current_gene,ncol(genes)-num_reps+x+1]),c(1,2),na.rm=TRUE)
    }
    else{ # (1,1,1) average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(c(as.numeric(y_val[z-num_reps]),genes[current_gene,z],as.numeric(y_val[z+num_reps])),c(1,1,1),na.rm=TRUE))
      first <- weighted.mean(c(genes[current_gene,2+x],as.numeric(y_val[2+num_reps+x])),c(1,1),na.rm=TRUE)
      last <- weighted.mean(c(as.numeric(y_val[ncol(genes)-(2*num_reps)+1+x]),genes[current_gene,ncol(genes)-num_reps+x+1]),c(1,1),na.rm=TRUE)
    }
    return (c(first,center.dat,last))}
  )
  dat <- genes[current_gene,] # original data
  for (x in 0:(num_reps-1)){ # replace original data by smoothed data, keeping na values
    ldat[[x+1]][is.na(genes[current_gene,seq(2+x,ncol(genes),by=num_reps)])] <- NA
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- ldat[[x+1]]
  }

  return(dat)
}

# DEPRECIATED
de_linear_trend <- function(current_gene, time_begin, time_end, resol,num_reps,timen,
                            rem_unexpr_vect){
  if (!rem_unexpr_vect[current_gene]){
    gene_n <- as.character(genes[current_gene,1]) # gene name
    rep_timen <- rep(timen,each=num_reps)
    y_val <- as.numeric(as.character(t(genes[current_gene,-1]))) # all the y values
    #do linear regression
    trend_test <- lm((y_val) ~ rep_timen)
    coeff <- trend_test$coefficients # resulting coefficients

    # detrend the data
    adjusted_y_val <- y_val - (coeff[1] + rep_timen*coeff[2])
    df2 <- cbind(data.frame("Gene.Name" = gene_n),rbind(adjusted_y_val))
    colnames(df2) <- colnames(genes)
    return (df2)
  } else {
    return (genes[current_gene,])
  }

}

# function to remove linear trend from all data
de_linear_trend_all <- function(timen,num_reps,tied){
  all_rep <- as.matrix(genes[,-1]) # y values for linear fit
  if (!tied){ # if they're not paired, we just fit an aggregate data model

    # x values for linear fit
    xrow <- rep(timen,each=num_reps)
    xmtx <- matrix(rep(xrow,each=nrow(all_rep)),nrow = nrow(all_rep))

    # covariance
    cov <- rowSums((all_rep-rowMeans(all_rep,na.rm = TRUE))*(xmtx-rowMeans(xmtx)),na.rm = TRUE)
    # variance
    var <- rowSums((xmtx - rowMeans(xmtx))^2,na.rm = TRUE)

    beta <- cov/var
    alph <- rowMeans(all_rep,na.rm = TRUE)-(beta*rowMeans(xmtx))

    df <- all_rep-(alph+(beta*xmtx)) # linear fit
  } else { # we have to do the models separately for each replicate
    # preallocate matrix where we put results
    df <- matrix(NA, nrow = dim(all_rep)[1], ncol = dim(all_rep)[2])

    # x values for linear fit
    xmtx <- matrix(rep(timen,each=nrow(all_rep)),nrow = nrow(all_rep))
    # preallocating to store slopes
    beta.df <- data.frame(matrix(NA, nrow(genes), num_reps))
    for (i in 1:num_reps){
      each_rep <- all_rep[,seq(i,ncol(all_rep),by=num_reps)]

      # covariance
      cov <- rowSums((each_rep-rowMeans(each_rep,na.rm = TRUE))*(xmtx-rowMeans(xmtx)),na.rm = TRUE)
      # variance
      var <- rowSums((xmtx - rowMeans(xmtx))^2,na.rm = TRUE)

      beta.df[,i] <- beta <- cov/var
      alph <- rowMeans(each_rep,na.rm = TRUE)-(beta*rowMeans(xmtx))

      df[,seq(i,ncol(all_rep),by=num_reps)] <- each_rep -(alph+(beta*xmtx)) # linear fit
    }
    # get average slope for each gene
    beta <- rowMeans(beta.df, na.rm = T)
  }
  # get the data frame correctly set up for returning
  res_df <- genes
  res_df[,-1] <- df

  # now we return the slope and the altered expressions
  res_list <- list("res_df" = res_df, "beta" = beta)

  return (res_list)
}

# function for getting an object of class nlsModel, EDITED VERY SPECIFICALLY FOR MOSAIC PROBLEM
# inputs:
#  form: formula specifying model
#  data: dataframe containing the data to be used
#  start: named initial values for parameters
#  wts: weights for data
# outputs:
#  object of class nlsModel
nlsModel_edit <- function (form, data, start, wts)
{
  thisEnv <- environment()
  env <- new.env(parent = environment(form))
  for (i in names(data)) {
    assign(i, data[[i]], envir = env)
  }
  ind <- as.list(start)
  parLength <- 0
  for (i in names(ind)) {
    temp <- start[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp, envir = env)
    ind[[i]] <- parLength + seq(along = start[[i]])
    parLength <- parLength + length(start[[i]])
  }
  useParams <- rep(TRUE, parLength)
  lhs <- eval(form[[2]], envir = env)
  rhs <- eval(form[[3]], envir = env)
  resid <- lhs - rhs
  dev <- sum(resid^2)
  if (is.null(attr(rhs, "gradient"))) {
    # the gradient doesn't matter, we're not evaluating anything
    # getRHS.noVarying <- function() numericDeriv(form[[3]],
    # names(ind), env)
    # getRHS <- getRHS.noVarying
    # rhs <- getRHS()
    attr(rhs, "gradient") <- matrix(1, nrow = length(data$y), ncol = length(names(ind)))
  }
  else {
    getRHS.noVarying <- function() eval(form[[3]], envir = env)
    getRHS <- getRHS.noVarying
  }
  dimGrad <- dim(attr(rhs, "gradient"))
  marg <- length(dimGrad)
  if (marg > 0) {
    gradSetArgs <- vector("list", marg + 1)
    for (i in 2:marg) gradSetArgs[[i]] <- rep(TRUE, dimGrad[i -
                                                              1])
    useParams <- rep(TRUE, dimGrad[marg])
  }
  else {
    gradSetArgs <- vector("list", 2)
    useParams <- rep(TRUE, length(attr(rhs, "gradient")))
  }
  npar <- length(useParams)
  gradSetArgs[[1]] <- (~attr(ans, "gradient"))[[2]]
  gradCall <- switch(length(gradSetArgs) - 1, call("[", gradSetArgs[[1]],
                                                   gradSetArgs[[2]]), call("[", gradSetArgs[[1]], gradSetArgs[[2]],
                                                                           gradSetArgs[[2]]), call("[", gradSetArgs[[1]], gradSetArgs[[2]],
                                                                                                   gradSetArgs[[2]], gradSetArgs[[3]]), call("[", gradSetArgs[[1]],
                                                                                                                                             gradSetArgs[[2]], gradSetArgs[[2]], gradSetArgs[[3]],
                                                                                                                                             gradSetArgs[[4]]))
  getRHS.varying <- function() {
    ans <- getRHS.noVarying()
    attr(ans, "gradient") <- eval(gradCall)
    ans
  }
  QR <- qr(attr(rhs, "gradient"))
  qrDim <- min(dim(QR$qr))
  # if (QR$rank < qrDim)
  #   stop("singular gradient matrix at initial parameter estimates")
  getPars.noVarying <- function() unlist(setNames(lapply(names(ind),
                                                         get, envir = env), names(ind)))
  getPars.varying <- function() unlist(setNames(lapply(names(ind),
                                                       get, envir = env), names(ind)))[useParams]
  getPars <- getPars.noVarying
  internalPars <- getPars()
  setPars.noVarying <- function(newPars) {
    assign("internalPars", newPars, envir = thisEnv)
    for (i in names(ind)) {
      assign(i, unname(newPars[ind[[i]]]), envir = env)
    }
  }
  setPars.varying <- function(newPars) {
    internalPars[useParams] <- newPars
    for (i in names(ind)) {
      assign(i, unname(internalPars[ind[[i]]]), envir = env)
    }
  }
  setPars <- setPars.noVarying
  on.exit(remove(i, data, parLength, start, temp, m))
  m <- list(resid = function() resid, fitted = function() rhs,
            formula = function() form, deviance = function() dev,
            gradient = function() attr(rhs, "gradient"), conv = function() {
              rr <- qr.qty(QR, resid)
              sqrt(sum(rr[1:npar]^2)/sum(rr[-(1:npar)]^2))
            }, incr = function() qr.coef(QR, resid), setVarying = function(vary = rep(TRUE,
                                                                                      length(useParams))) {
              assign("useParams", if (is.character(vary)) {
                temp <- logical(length(useParams))
                temp[unlist(ind[vary])] <- TRUE
                temp
              } else if (is.logical(vary) && length(vary) != length(useParams)) stop("setVarying : vary length must match length of parameters") else {
                vary
              }, envir = thisEnv)
              gradCall[[length(gradCall)]] <<- useParams
              if (all(useParams)) {
                assign("setPars", setPars.noVarying, envir = thisEnv)
                assign("getPars", getPars.noVarying, envir = thisEnv)
                assign("getRHS", getRHS.noVarying, envir = thisEnv)
                assign("npar", length(useParams), envir = thisEnv)
              } else {
                assign("setPars", setPars.varying, envir = thisEnv)
                assign("getPars", getPars.varying, envir = thisEnv)
                assign("getRHS", getRHS.varying, envir = thisEnv)
                assign("npar", length((1:length(useParams))[useParams]),
                       envir = thisEnv)
              }
            }, setPars = function(newPars) {
              setPars(newPars)
              assign("resid", lhs - assign("rhs", getRHS(), envir = thisEnv),
                     envir = thisEnv)
              assign("dev", sum(resid^2), envir = thisEnv)
              assign("QR", qr(attr(rhs, "gradient")), envir = thisEnv)
              return(QR$rank < min(dim(QR$qr)))
            }, getPars = function() getPars(), getAllPars = function() getPars(),
            getEnv = function() env, trace = function() cat(format(dev),
                                                            ": ", format(getPars()), "\n"), Rmat = function() qr.R(QR),
            predict = function(newdata = list(), qr = FALSE) {
              eval(form[[3]], as.list(newdata), env)
            })
  class(m) <- "nlsModel"
  m
}

# function for getting an object of class nls, EDITED VERY SPECIFICIALLY FOR MOSAIC PROBLEM
# inputs:
#  formula: a nonlinear model formula including variables and parameters. Will be coerced to a formula if necessary.
#  data: an optional data frame in which to evaluate the variables in formula and weights. Can also be a list or an environment, but not a matrix.
# start: a named list or named numeric vector of starting estimates. When start is missing (and formula is not a self-starting model, see selfStart), a very cheap guess for start is tried (if algorithm != "plinear").
# control: an optional list of control settings. See nls.control for the names of the settable control values and their effect.
# algorithm: character string specifying the algorithm to use. The default algorithm is a Gauss-Newton algorithm. Other possible values are "plinear" for the Golub-Pereyra algorithm for partially linear least-squares models and "port" for the 'nl2sol' algorithm from the Port library - see the references. Can be abbreviated.
# trace: logical value indicating if a trace of the iteration progress should be printed. Default is FALSE. If TRUE the residual (weighted) sum-of-squares and the parameter values are printed at the conclusion of each iteration. When the "plinear" algorithm is used, the conditional estimates of the linear parameters are printed after the nonlinear parameters. When the "port" algorithm is used the objective function value printed is half the residual (weighted) sum-of-squares.
# subset: an optional vector specifying a subset of observations to be used in the fitting process.
# weights: an optional numeric vector of (fixed) weights. When present, the objective function is weighted least squares.
# na.action: a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit. Value na.exclude can be useful.
# model: logical. If true, the model frame is returned as part of the object. Default is FALSE.
# lower, upper: vectors of lower and upper bounds, replicated to be as long as start. If unspecified, all parameters are assumed to be unconstrained. Bounds can only be used with the "port" algorithm. They are ignored, with a warning, if given for other algorithms.
# ...: Additional optional arguments. None are used at present.
# outputs:
#  a list containing:
# m: an nlsModel object incorporating the model.
# data: the expression that was passed to nls as the data argument. The actual data values are present in the environment of the m component.
# call: the matched call with several components, notably algorithm.
# na.action: the "na.action" attribute (if any) of the model frame.
# dataClasses: the "dataClasses" attribute (if any) of the "terms" attribute of the model frame.
# model: if model = TRUE, the model frame.
# weights: if weights is supplied, the weights.
# convInfo: a list with convergence information.
# control: the control list used, see the control argument.
# convergence, message: for an algorithm = "port" fit only, a convergence code (0 for convergence) and message. To use these is deprecated, as they are available from convInfo now.
nls_edit <- function (formula, data = parent.frame(), start, control = nls.control(),
                      algorithm = c("default", "plinear", "port"), trace = FALSE,
                      subset, weights, na.action, model = FALSE, lower = -Inf,
                      upper = Inf, ...)
{
  formula <- as.formula(formula)
  algorithm <- match.arg(algorithm)
  if (!is.list(data) && !is.environment(data))
    stop("'data' must be a list or an environment")
  mf <- cl <- match.call()
  varNames <- all.vars(formula)
  if (length(formula) == 2L) {
    formula[[3L]] <- formula[[2L]]
    formula[[2L]] <- 0
  }
  form2 <- formula
  form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)
  mWeights <- missing(weights)
  pnames <- if (missing(start)) {
    if (!is.null(attr(data, "parameters"))) {
      names(attr(data, "parameters"))
    }
    else {
      cll <- formula[[length(formula)]]
      if (is.symbol(cll)) {
        cll <- substitute(S + 0, list(S = cll))
      }
      fn <- as.character(cll[[1L]])
      if (is.null(func <- tryCatch(get(fn), error = function(e) NULL)))
        func <- get(fn, envir = parent.frame())
      if (!is.null(pn <- attr(func, "pnames")))
        as.character(as.list(match.call(func, call = cll))[-1L][pn])
    }
  }
  else names(start)
  env <- environment(formula)
  if (is.null(env))
    env <- parent.frame()
  if (length(pnames))
    varNames <- varNames[is.na(match(varNames, pnames))]
  lenVar <- function(var) tryCatch(length(eval(as.name(var),
                                               data, env)), error = function(e) -1L)
  if (length(varNames)) {
    n <- vapply(varNames, lenVar, 0)
    if (any(not.there <- n == -1L)) {
      nnn <- names(n[not.there])
      if (missing(start)) {
        if (algorithm == "plinear")
          stop("no starting values specified")
        warning("No starting values specified for some parameters.\n",
                "Initializing ", paste(sQuote(nnn), collapse = ", "),
                " to '1.'.\n", "Consider specifying 'start' or using a selfStart model",
                domain = NA)
        start <- setNames(as.list(rep_len(1, length(nnn))),
                          nnn)
        varNames <- varNames[i <- is.na(match(varNames,
                                              nnn))]
        n <- n[i]
      }
      else stop(gettextf("parameters without starting value in 'data': %s",
                         paste(nnn, collapse = ", ")), domain = NA)
    }
  }
  else {
    if (length(pnames) && any((np <- sapply(pnames, lenVar)) ==
                              -1)) {
      message(sprintf(ngettext(sum(np == -1), "fitting parameter %s without any variables",
                               "fitting parameters %s without any variables"),
                      paste(sQuote(pnames[np == -1]), collapse = ", ")),
              domain = NA)
      n <- integer()
    }
    else stop("no parameters to fit")
  }
  respLength <- length(eval(formula[[2L]], data, env))
  if (length(n) > 0L) {
    varIndex <- n%%respLength == 0
    if (is.list(data) && diff(range(n[names(n) %in% names(data)])) >
        0) {
      mf <- data
      if (!missing(subset))
        warning("argument 'subset' will be ignored")
      if (!missing(na.action))
        warning("argument 'na.action' will be ignored")
      if (missing(start))
        start <- getInitial(formula, mf)
      startEnv <- new.env(hash = FALSE, parent = environment(formula))
      for (i in names(start)) assign(i, start[[i]], envir = startEnv)
      rhs <- eval(formula[[3L]], data, startEnv)
      n <- NROW(rhs)
      wts <- if (mWeights)
        rep_len(1, n)
      else eval(substitute(weights), data, environment(formula))
    }
    else {
      vNms <- varNames[varIndex]
      if (any(nEQ <- vNms != make.names(vNms)))
        vNms[nEQ] <- paste0("`", vNms[nEQ], "`")
      mf$formula <- as.formula(paste("~", paste(vNms,
                                                collapse = "+")), env = environment(formula))
      mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
      mf$lower <- mf$upper <- NULL
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights)
        model.weights(mf)
      else rep_len(1, n)
    }
    if (any(wts < 0 | is.na(wts)))
      stop("missing or negative weights not allowed")
  }
  else {
    varIndex <- logical()
    mf <- list(0)
    wts <- numeric()
  }
  if (missing(start))
    start <- getInitial(formula, mf)
  for (var in varNames[!varIndex]) mf[[var]] <- eval(as.name(var),
                                                     data, env)
  varNamesRHS <- varNamesRHS[varNamesRHS %in% varNames[varIndex]]
  m <- switch(algorithm, plinear = nlsModel_edit.plinear(formula,
                                                         mf, start, wts), port = nlsModel_edit(formula, mf, start,
                                                                                               wts, upper), nlsModel_edit(formula, mf, start, wts))
  ctrl <- nls.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if (algorithm != "port") {
    if (!identical(lower, -Inf) || !identical(upper, +Inf)) {
      warning("upper and lower bounds ignored unless algorithm = \"port\"")
      cl$lower <- NULL
      cl$upper <- NULL
    }
    if (ctrl$maxiter > 0){
      convInfo <- .Call(C_nls_iter, m, ctrl, trace)
    } else {
      convInfo <- list("isConv" = F, "finIter" = 0, "finTol" = NA, "stopCode" = 9,
                       "stopMessage"="The number of iterations has reached maxiter.")
    }
    nls.out <- list(m = m, convInfo = convInfo, data = substitute(data),
                    call = cl)
  }
  else {
    pfit <- nls_port_fit(m, start, lower, upper, control,
                         trace, give.v = TRUE)
    iv <- pfit[["iv"]]
    msg.nls <- port_msg(iv[1L])
    conv <- (iv[1L] %in% 3:6)
    if (!conv) {
      msg <- paste("Convergence failure:", msg.nls)
      if (ctrl$warnOnly)
        warning(msg)
      else stop(msg)
    }
    v. <- port_get_named_v(pfit[["v"]])
    cInfo <- list(isConv = conv, finIter = iv[31L], finTol = v.[["NREDUC"]],
                  nEval = c(`function` = iv[6L], gradient = iv[30L]),
                  stopCode = iv[1L], stopMessage = msg.nls)
    cl$lower <- lower
    cl$upper <- upper
    nls.out <- list(m = m, data = substitute(data), call = cl,
                    convInfo = cInfo, convergence = as.integer(!conv),
                    message = msg.nls)
  }
  nls.out$call$algorithm <- algorithm
  nls.out$call$control <- ctrl
  nls.out$call$trace <- trace
  nls.out$na.action <- attr(mf, "na.action")
  nls.out$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  if (model)
    nls.out$model <- mf
  if (!mWeights)
    nls.out$weights <- wts
  nls.out$control <- control
  class(nls.out) <- "nls"
  nls.out
}

