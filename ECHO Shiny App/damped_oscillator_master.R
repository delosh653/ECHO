# Extended Oscillations Function Source
# By Hannah De los Santos
# ECHO v 3.0
# Code description: Contains all the funcitons for extended harmonic oscillator work, in order to have less confusion between scripts.

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
      all_fits[[i]] <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                              data=temp_edit,
                                              start=, start_param,
                                              lower=c(-Inf, -Inf, high, -Inf, min(temp_edit$y)),
                                              upper=c(Inf, Inf, low, Inf, max(temp_edit$y)),
                                              control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)))
    } else {
      #fitting
      all_fits[[i]] <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                              data=temp_edit,
                                              start=, start_param,
                                              lower=c(-Inf, -Inf, high, -Inf, min(temp_edit$y)),
                                              upper=c(Inf, Inf, low, Inf, max(temp_edit$y)),
                                              control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6),
                                              weights = w))
    }
    
    for (p in 1:5){
      all_pars[[p]] <- c(all_pars[[p]],all_fits[[i]]$m$getAllPars()[p])
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


# Function to calculate the parameters for the extended harmonic oscillator equation for a specific gene.
#
#  current_gene: row number of current gene we want to calculate parameters for
#  times: time points for dataset
#  resol: resolution of time points
#  num_reps: number of replicates
#  tied: if replicate data, whether the replicates are related (paired) or not (unpaired)
#  is_smooth: boolean that indicates whether data should be smoothed or not
#  is_weighted: if there is smoothing, is it weighted (1,2,1) smoothing, or unweighted smoothing (1,1,1)
#  low: the highest frequency we are looking for, in radians (lowest period)
#  high: the lowest frequency we are looking for, in radians (highest period)
#  rem_unexpr: logical indicating whether genes with less than rem_unexpr_amt % expression should not be considered
#  rem_unexpr_amt: percentage of expression for which genes should not be considered
#  jtklist: contains the exact p-value distribution for replicate data
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
#   original.values: original values for gene
#   fitted.values: fitted values for gene
calculate_param <- function(current_gene,times,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, harm_cut = .03, over_cut = .15){
  
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
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
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
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
      }
      return(results)
    }
  }
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- rbind(as.numeric(as.character(t(genes[current_gene,c(2:ncol(genes))])))) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  tryCatch({ # throw exception upon error
    # calculate the amount of peaks
    peaks <- c(); # vector of peak values
    peaks_time <- c(); # vector of peak times
    counting <- 1; # counter
    {
    if (resol <= ((1/12)+10^-8)){ # 17 hour surround
      mod <- 102
      for(i in (mod+1):(length(y_val)-mod)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
    } else if (resol <= ((1/6)+10^-8)){ # 15 hour surround
      mod <- 45
      for(i in (mod+1):(length(y_val)-mod)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
    } else if (resol <= ((1/4)+10^-8)){ # 13 hour surround
      mod <- 26
      for(i in (mod+1):(length(y_val)-mod)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
    } else if (resol <= ((1/2)+10^-8)){ # 11 hour surround
      mod <- 11
      for(i in (mod+1):(length(y_val)-mod)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
    } else if (resol <= 1){
      # go through gene values and find maximum as compared to 8 surrounding values
      # finding peaks for first 4 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:9], na.rm = TRUE)) != -Inf){
      #   for (i in 1:4){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:9], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
      
      for(i in 5:(length(y_val)-4)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[i-4],y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3],y_val[i+4], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[i-4],y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3],y_val[i+4], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
      
      # finding peaks for last 4 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-8):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in (length(y_val)-3):length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-8):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
    } else if (resol <=2){
      # go through gene values and find maximum as compared to six surrounding values
      # finding peaks for first 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:7], na.rm = TRUE)) != -Inf){
      #   for (i in 1:3){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:7], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
      for(i in 4:(length(y_val)-3)){
        # deal with complete missingness
        if (suppressWarnings(max(y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
          next
        }
        # otherwise continue as normal
        if (y_val[i] == max(y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
      # finding peaks for last 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-6):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in (length(y_val)-2):length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-6):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
    } else if (resol <= 4){
      # finding peaks for first 2 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:5], na.rm = TRUE)) != -Inf){
      #   for (i in 1:2){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:5], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
      # go through gene values and find maximum as compared to four surrounding values
      for(i in 3:(length(y_val)-2)){
        # to deal with complete missingness
        if(suppressWarnings(max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE))== -Inf  | is.na(y_val[i])){
          next
        }
        if (y_val[i] == max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
      
      # finding peaks for last 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-4):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in (length(y_val)-1):length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-4):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
      
    } else{
      # finding peaks for first point
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[1:3], na.rm = TRUE)) != -Inf){
      #   for (i in 1){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[1:3], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
      
      # go through gene values and find maximum as compared to two surrounding values
      for(i in 2:(length(y_val)-1)){
        # to deal with complete missingness
        if(suppressWarnings(max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE))== -Inf  | is.na(y_val[i])){
          next
        }
        if (y_val[i] == max(y_val[i-1],y_val[i],y_val[i+1], na.rm = TRUE)){
          peaks[counting] <- y_val[i]
          peaks_time[counting] <- times[i]
          counting <- counting+1
        }
      }
      
      # finding peaks for last 3 points
      # deal with complete missingness
      # if (suppressWarnings(max(y_val[(length(y_val)-2):length(y_val)], na.rm = TRUE)) != -Inf){
      #   for (i in length(y_val)){
      #     # otherwise continue as normal
      #     if (y_val[i] == max(y_val[(length(y_val)-2):length(y_val)], na.rm = TRUE)){
      #       peaks[counting] <- y_val[i]
      #       peaks_time[counting] <- times[i]
      #       counting <- counting+1
      #     }
      #   }
      # }
    }
    }
    # calculate starting amplitude, y_shift
    y0 <- mean(y_val,na.rm = TRUE) # intial value for the equilibrium shift
    if (y0 < 10^-10 && y0 > -10^-10){
      y0 <- 10^-8 # avoiding 0 mean, which makes gradient singular
    }
    x0 <- min(times) # the x start parameter
    a0 <- max(y_val,na.rm = TRUE) - y0 # mean(y_val) # initial guess for amplitude
    
    # intial value for gamma
    if (length(peaks)==0){ # if there are no peaks, we account for that
      gam0 <- 0
    } else if (which.max(peaks)==1){ # if the highest peak is first, then damping is likely
      if (length(peaks)>1){
        # trying to figure out gamma based on logarithmic decrement
        n <- peaks_time[2]-peaks_time[1]
        log_dec <- (1/n)*log(abs(peaks[1]/peaks[2]))
        gam0 <- 1/(sqrt(1+((2*pi/log_dec)^2)))
      } else{ # use a guess
        gam0 <- .01
      }
    } else{ # otherwise driving is likely
      if (length(peaks)>1){
        # trying to figure out gamma based on logarithmic decrement
        n <- peaks_time[2]-peaks_time[1]
        log_dec <- (1/n)*log(abs(peaks[2]/peaks[1]))
        gam0 <- -1*1/(sqrt(1+((2*pi/log_dec)^2)))
      } else{ # use a guess
        gam0 <- -.01
      }
    }
    
    # let frequency depend on amount of peaks = (length(times)*resol/(no of peaks+1 [accounts for phase shift])
    if (length(peaks) == 0){
      if (high == -Inf || low == Inf){
        w0 <- 2*pi/(length(times)*resol/2)
      } else{
        # want to get their actual integer period values
        highfix <- (high/2/pi)^-1
        lowfix <- (low/2/pi)^-1
        w0 <- 2*pi/(length(times)*resol/((highfix+lowfix)/2))
      }
    } else if (length(peaks) == 1){ # phase shift causes only one peak to appear
      w0 <- 2*pi/(length(times)*resol/(length(peaks)+1))
    } else{
      w0 <- 2*pi/(length(times)*resol/(length(peaks)))
    }
    
    # can't be outside the specified parameters
    if (w0 > low){
      w0 <- low
    } else if (w0 < high){
      w0 <- high
    }
    
    
    # we estimate our phase shift on the second and third nonmissing points for accuracy
    # if you have less than 3 points nonmissing, I have no hope for you
    second <- which(!is.na(y_val))[2]
    third <- which(!is.na(y_val))[3]
    min_i <- 0;
    min_vect <- rep(10000, length(0:11))
    for (i in 0:11){
      # if the phase shift at both the second and third gene values are the smallest value available for the fitted value
      min_vect[i+1] <- sum(alt_form(a0,gam0,w0,(i*pi/6),y0,times[second])-y_val[second],
                           alt_form(a0,gam0,w0,(i*pi/6),y0,times[third])-y_val[third])
    }
    # phi0 <- min_i*pi/6 # intial value for phase shift
    phi0 <- (which.min(abs(min_vect))-1)*pi/6 # intial value for phase shift
    
    start_param <- list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0)
    if (num_reps == 1){ # one replicate
      # put the times into a data frame
      temp <- data.frame(y=t(y_val),t=times)
      temp <- temp[!is.na(temp$y),] # remove any missing data points
      
      # fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=start_param,
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)))
    } else{ # multiple replicates
      #put the times and data point into a data frame
      weights <- calc_weights(current_gene,num_reps)
      temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(times,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
      temp <- temp[!is.na(temp$y),] # remove any missing data points
      
      #fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=start_param,
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6),
                                               weights = w))
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
      ci_int <- jackknife(temp, parameters, num_reps, start_param)
    }
    
    # calculating whether (over)damped, (over)forced, harmonic
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
    
    # calculate p-value
    ref_wave <- (alt_form(a,gam,omega,phi,y_shift,times)) # fitted values
    all_pred <- rep(ref_wave, each=num_reps)[!is.na(unlist(genes[current_gene,-1]))]
    pval <- cor.test(all_pred,temp$y, method = "kendall")$p.value
    tau <- cor.test(all_pred,temp$y, method = "kendall")$estimate
    
    # list of parameters and other resulting values
    results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
    if (!run_conf){
      if (num_reps == 1){
        results <- cbind(results, y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    } else {
      if (num_reps == 1){
        results <- cbind(results, rbind(ci_int), y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    }
    return (results)
  }, error = function(e){ # if there's failure in convergence
    results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
    if (!run_conf){
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
    } else {
      results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
    }
    
    return (results)
  })
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
