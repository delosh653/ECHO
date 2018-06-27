# Extended Oscillations Function Source
# By Hannah De los Santos
# ECHO v 1.83
# Code description: Contains all the funcitons for extended harmonic oscillator work, in order to have less confusion between scripts.

# function to represent damped oscillator with phase and equilibrium shift formula
# inputs:
#  a: Amplitude
#  gam: Forcing.Coefficient (amount of damping/driving)
#  omega: Radial frequency
#  phi: Phase Shift (radians)
#  y_shift: Equilibrium shift
#  t: time
# outputs:
#  result of inputs into formula
alt_form <- function(a,gam,omega,phi,y_shift,t){
  return (a*exp(-1*gam*(t)/2)*cos((omega*t)+phi)+y_shift)
}

# p value calculation functions

# function to calculate Kendall's tau value
# inputs:
#  ref_waveform: answer to damped oscillator
#  y_val: actual values
# outputs:
#  returns kendall's tau value
calc_kendalls_tau <- function(ref_waveform, y_val){
  # both the reference waveform and y_val are lists, so we need
  # to turn them into vectors first
  return (cor((unlist(y_val)), unlist(ref_waveform), method = "kendall"))
}

# function to calculate frequencies for CDF
# inputs:
#  n: number of time points
# outputs:
#  freq_list: frequencies for CDF
calc_freq <- function(n){
  # i'm not going to save the results at each step, but it's certainly a possibility
  # might be a misnomer -- think about changing this
  tot_inversions <- n*(n-1)/2 # amount of pairs for kendall's tau
  freq_n <- numeric(tot_inversions+1) # number of frequencies for CDF
  freq_n[1] <- 1.0

  last <- numeric(tot_inversions+1) # place to store intermediate steps of frequency vector

  # loops to go through the amount of inversions and collect frequencies
  # essentially we're creating a normal curve
  freq_list <- list()
  freq_list[[1]] <- 1.0
  l_index <- 2 # counter to keep track of index size
  for (i in 2:n-1){
    last <- freq_n;
    for (j in 1:(tot_inversions+1)){
      if ((j-1) >= max(1,j-i)){ # in order to make sure that it doesn't exceed index size
        freq_n[j] <- freq_n[j] + sum(last[max(1,j-i):(j-1)])
      }
    }
    freq_list[[l_index]] <- freq_n[1:(1+((i+1)*i/2))]
    l_index <- l_index + 1
  }

  return (freq_list)
}

# function to calculate CDF based on frequencies
# inputs:
#  n: number of time points
# outputs:
#  cdfs: CDF for number of time points (normal distribution)
calc_cdf <- function(n){
  freq_list <- calc_freq(n) # calculate frquencies

  cdfs <- list() # declaring the list of cdfs
  for (i in 1:n){ # calculating the cdfs
    cdfs[[i]] <- freq_list[[i]]
    if (length(freq_list[[i]]) >= 2){
      for (j in 2:length(freq_list[[i]])){
        cdfs[[i]][j] = cdfs[[i]][j] + cdfs[[i]][j-1]
      }
    }
    cdfs[[i]] = cdfs[[i]]/sum(freq_list[[i]])
  }

  return (cdfs)
}

# function to calculate p-value based on Kendall's tau (only used for 1 replicate)
# inputs:
#  n: number of time points
#  tau: Kendall's tau value
#  two_tail: boolean value indicating whether or not a two tailed p-value is wanted
# outputs:
#  pvalue based on kendall's tau
calc_p_value <- function(n, tau, two_tail){
  # error checking
  # epsilon <- .0000001 # machine accuracy
  # if (tau <= (-1-epsilon) + tau >= (1+epsilon)){
  #   print("invalid tau")
  #   return
  # }

  tot_inversions <- n*(n-1)/2 # total number of inversions
  #convert tau to value in terms of number of inversions
  inv_score <- round((tot_inversions - (tau*tot_inversions))/2.0) # automatic round to nearest int

  # for two tailed pvalue, we're getting tail from a symmetric distribution
  min_inv_score <- min(inv_score, tot_inversions - inv_score)

  # calcluate the cdfs
  cdfs <- calc_cdf(n)
  # calculate the pvalue from the CDF
  if (two_tail){
    pval = cdfs[[n]][min_inv_score+1]*2.0
  } else{ # one tailed: probability of getting that amount or fewer inversions
    pval = cdfs[[n]][min_inv_score+1]
  }

  # if inv_score is 0, might have larger than .5 prob
  return (min(pval,1.0))

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

# Shewchuk algorithms for adaptive precision summation used in jtkdist
# http://www.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
fast.two.sum <- function(a,b) {                   # known abs(a) >= abs(b)
  x <- a+b
  bv <- x-a
  y <- b-bv
  if(y==0) return(x)
  c(x,y)
}

two.sum <- function(a,b) {                        # unknown order
  x <- a+b
  bv <- x-a
  av <- x-bv
  br <- b-bv
  ar <- a-av
  y <- ar+br
  if(y==0) return(x)
  c(x,y)
}

expansion.sum <- function(g) {
  g <- g[order(abs(g))]
  z <- fast.two.sum(g[2],g[1])
  q <- z[1]
  h <- NULL
  if(length(z)!=1) h <- z[2]
  n <- length(g)
  if(n==2) return(c(h,q))
  for(i in 3:n) {
    z <- two.sum(q,g[i])
    q <- z[1]
    if(length(z)!=1) h <- c(h,z[2])
  }
  c(h,q)                                          # strongly non-overlapping values
}

# function to calculate exact distribution for jonckheere-terpsta-kendall's S using the Harding algorithm
#  details: http://www.jstor.org/pss/2347656
# inputs:
#  timepoints: number of total amount of timepoints
#  reps: number of replicates
#  normal: logical indicating whether a normal distribution is desired
#  alt: logical indicating whether there are missing time points
# outputs:
#  exact distribution of JTK, whether that be based on missing values or not
jtkdist <- function(timepoints,reps=1,normal=FALSE,alt=FALSE) {


  if(length(reps)==timepoints) {
    tim <- reps # support for unbalanced replication
  } else {
    tim <- rep(reps[1],timepoints) # balanced replication
  }

  maxnlp <- lfactorial(sum(tim))-sum(lfactorial(tim)) #maximum possible negative log p-value
  limit <- log(.Machine$double.xmax) #largest representable nlp
  normal <- normal | (maxnlp>limit-1) #switch to normal approximation if maxnlp is too large

  if(alt) { # if there are missing time points
    lab <- paste(sort(tim), collapse=",") # create label
    if(lab %in% names(jtk.alt)) { # if this distribution has already been created
      alt.id <- match(lab, names(jtk.alt))
      return(alt.id) # no need to create a new distribution, use the the old label
    }
    nn <- sum(tim) # number of time points
    M <- (nn^2-sum(tim^2))/2 # number of combinations
    jtk.alt[[lab]] <<- list() # save in the global space for efficiency
    jtk.alt[[lab]]$MAX <<- M
    if(normal) { # create normal distribution is desired or necessary
      var <- (nn^2*(2*nn+3) -
                sum(tim^2*(2*tim+3)))/72 # variance
      jtk.alt[[lab]]$SDV <<- sqrt(var) # standard deviation
      jtk.alt[[lab]]$EXV <<- M/2 # expected value2
      return(length(jtk.alt))
    }
  } else { # exact distribution
    jtk.grp.size <- tim # sizes of each replicate group
    jtk.num.grps <- length(tim) # timepoints = number of groups
    jtk.num.vals <- nn <- sum(tim) # number of data values (independent of period and lag)
    jtk.max <- M <- (nn^2-sum(tim^2))/2 # maximum possible jtk statistic
    jtk.grps <- rep(1:length(tim), ti=tim) # group labels
    jtk.dims <- c(nn*(nn-1)/2,1) # dimensions for distribution

    if(normal) { # unused, practically
      jtk.var <- (nn^2*(2*nn+3) -
                    sum(tim^2*(2*tim+3)))/72 # variance of jtk
      jtk.sdv <- sqrt(jtk.var) # standard deviation of jtk
      jtk.exv <- jtk.max/2 # expected value of jtk
      jtk.exact <- FALSE
      
      jtklist <- list("jtk.grp.size"=jtk.grp.size,"jtk.num.grps"=jtk.num.grps,"jtk.num.vals"=jtk.num.vals,"jtk.max"=jtk.max,"jtk.grps"=jtk.grps,"jtk.dims"=jtk.dims,"jtk.exact"=jtk.exact, "jtk.var"=jtk.var, "jtk.sdv" = jtk.sdv, "jtk.exv"=jtk.exv)
      
      return(jtklist) # omit calculation of exact distribution
    }
  }
  MM <- floor(M/2) # mode of this possibly alternative jtk distribution
  cf <- as.list(rep(1,MM+1)) # initial lower half cumulative frequency distribution

  size <- tim # sizes of each group of known replicate values
  size <- size[order(size)]  # ascending order for fastest calculation
  k <- length(tim) # number of groups of known replicate values

  # start building distribution
  N <- size[k]
  if(k>2) for(i in (k-1):2) {
    N <- c(size[i]+N[1],N)
  }
  for(i in 1:(k-1)) { # count permutations using the Harding algorithm
    m <- size[i]
    n <- N[i]

    if(n < MM) {
      P <- min(m+n,MM)
      for(t in (n+1):P) { # zero-based offset t
        for(u in 1+MM:t) { # one-based descending index u
          cf[[u]] <- expansion.sum( # Shewchuck algorithm
            c(cf[[u]],-cf[[u-t]]))
        }
      }
    }
    Q <- min(m,MM)
    for(s in 1:Q) { # zero-based offset s
      for(u in 1+s:MM) { # one-based ascending index u
        cf[[u]] <- expansion.sum( # Shewchuck algorithm
          c(cf[[u]],cf[[u-s]]))
      }
    }
  }
  cf <- sapply(cf,sum)

  # cf now contains the lower-half cumulative frequency distribution;
  # append the symmetric upper-half cumulative distribution to cf

  if(M %% 2) {
    cf <- c(cf,2*cf[MM+1]-c(cf[MM:1],0)) # if M is odd (mode is duplicated)
  } else {
    cf <- c(cf,cf[MM+1]+cf[MM]-c(cf[MM:2-1],0)) # if M is even (unique mode is in lower half)
  }
  jtkcf <- rev(cf) # upper-tail cumulative frequencies for all integer jtk
  ajtkcf <- (jtkcf[-length(cf)]+jtkcf[-1])/2 # interpolated cumulative frequency values for all half-integer jtk

  id <- 1+0:(2*M) # one-based indices for all jtk values
  cf <- id # container for the jtk frequency distribution
  cf[!!id%%2] <- jtkcf # odd indices for integer jtk
  cf[!id%%2] <- ajtkcf # even indices for half-integer jtk
  cp <- cf/jtkcf[1] # all upper-tail p-values

  if(alt) { # if missing values
    jtk.alt[[lab]]$CP <<- cp
    return(length(jtk.alt))
  }
  jtk.cp <- cp
  jtk.exact <- TRUE

  # list of relevant statistics
  jtklist <- list("jtk.grp.size"=jtk.grp.size,"jtk.num.grps"=jtk.num.grps,"jtk.num.vals"=jtk.num.vals,"jtk.max"=jtk.max,"jtk.grps"=jtk.grps,"jtk.dims"=jtk.dims,"jtk.cp"=jtk.cp,"jtk.exact"=jtk.exact)

  return(jtklist)
}

# function to calculate pvalues for replicates using the exact JTK distribution
# inputs:
#  ref_waveform: answer to damped oscillator
#  times: vector of time points
#  num_reps: number of replicates
#  current_gene: row number of gene being examined
#  jtklist: exact distribution for JTK
# outputs:
#  list containing pvalue, Kendall's tau, and Kendall's S statistic
calc_tau_and_p_jtk <- function(ref_waveform,times,num_reps,current_gene,jtklist){
  # forming the correlation matrix for the experimental values
  rep_numerator <- sign(outer(unlist(genes[current_gene,-1]),unlist(genes[current_gene,-1]),"-"))
  rep_numerator <- rep_numerator[lower.tri(rep_numerator)]

  # forming the correlation matrix for the fitted values
  ref_waveform_rep <- rep(ref_waveform,each=num_reps)
  ref_numerator <- sign(outer(ref_waveform_rep,ref_waveform_rep,"-"))
  ref_numerator <- ref_numerator[lower.tri(ref_numerator)]
  S <- sum(rep_numerator*ref_numerator,na.rm = TRUE) # Kendall's S statistic

  tim <- rep(times, each = num_reps)[!is.na(genes[current_gene,-1])] # time points examined (not including points with missing values)
  tim <- as.integer(table(tim)) # getting the number of occurences for each time point, in vector form
  nn <- sum(tim) # sum number of occurences

  alt <- any(is.na(unlist(genes[current_gene,-1]))) # are there any missing values?

  if(alt) { # if missing values
    tab <- table(jtklist$jtk.grps[is.finite(unlist(genes[current_gene,-1]))])
    alt.id <- jtkdist(length(tab),
                      as.integer(tab),
                      alt=alt)
  }
  M <- switch(1+alt, # maximum possible S score for this distribution
              jtklist$jtk.max, # no missing values
              jtk.alt[[alt.id]]$MAX) # if missing values

  if(!S){return (list("pval"=1,"S"=0,"tau"=0))} # if not calculated, must have a pvalue of 1
  jtk <- (abs(S)+M)/2 # two-tailed jtk statistic

  if(jtklist$jtk.exact) { # if we are using the exact distribution
    jtki <- 1+2*jtk # index into the exact upper-tail distribution
    p <- switch(1+alt, # deciding between missing values (2) or not (1)
                2*jtklist$jtk.cp[jtki],
                2*jtk.alt[[alt.id]]$CP[jtki])
  } else { # if we are using the normal distribution (unused, practically)
    p <- switch(1+alt,
                2*pnorm(-(jtk-1/2),
                        -jtklist$jtk.exv,jtklist$jtk.sdv),
                2*pnorm(-(jtk-1/2),
                        -jtk.alt[[alt.id]]$EXV,
                        jtk.alt[[alt.id]]$SDV))
  }
  # include tau = S/M for this distribution

  return(list("pval"=p,"S"=S,"tau"=S/M))
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
#   type_gam: Type of oscillation (damped, driven, etc.)
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
calculate_param <- function(current_gene,times,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,jtklist=list()){

  # if(is_smooth){ # smooth the data, if requested
  #   if (tied){ # if paired replicates
  #     genes[current_gene,] <- smoothing_tied(current_gene, is_weighted, num_reps)
  #   }
  #   else{ # if unpaired replicates
  #     genes[current_gene,] <- smoothing_untied(current_gene, is_weighted, num_reps)
  #   }
  # }

  gene_n <- as.character(genes[current_gene,1]) # gene name
  # first we need to check whether or not the gene is just a straight line
  if (!is_deviating(current_gene)){ # one replicate
    if (num_reps == 1){
      results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
      return (results)
    }
    else{ # multiple replicates
      results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
      return (results)
    }
  }

  # then we need to check if 70% are expressed (if desired)
  if (rem_unexpr){
    if (rem_unexpr_vect[current_gene]){
      if (num_reps == 1){ # one replicate

        results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
        return (results)

      } else{ # multiple replicates

        results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
        return (results)
      }
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
    for (i in 0:11){
      # if the phase shift at both the second and third gene values are the smallest value available for the fitted value
      if ((abs(alt_form(a0,gam0,w0,(i*pi/6),y0,times[second])-y_val[second]))<(abs(alt_form(a0,gam0,w0,(min_i*pi/6),y0,times[second])-y_val[second]))){
        if ((abs(alt_form(a0,gam0,w0,(i*pi/6),y0,times[third])-y_val[third]))<(abs(alt_form(a0,gam0,w0,(min_i*pi/6),y0,times[third])-y_val[third]))){
          min_i <- i
        }
      }
    }
    phi0 <- min_i*pi/6 # intial value for phase shift

    if (num_reps == 1){ # one replicate
      # put the times into a data frame
      temp <- data.frame(y=t(y_val),t=times)
      temp <- temp[!is.na(temp$y),] # remove any missing data points

      # fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0),
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
                                               start=list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0),
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

    # calculating whether (over)damped, (over)driven, harmonic
    if (gam < -.15){
      type_gam <- "Overexpressed"
    } else if (gam <= -.01){
      type_gam <- "Driven"
    } else if (gam <= .01){
      type_gam <- "Harmonic"
    } else if (gam <= .15){
      type_gam <- "Damped"
    } else{
      type_gam <- "Repressed"
    }

    # calculating the phase shift in terms of period (omega inverse of period)
    if (phi > 0){ # shift to the left
      frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
      dist_peak <- frac_part*(2*pi/omega) # distance from first peak
      phase_hours <- (2*pi/omega)-dist_peak
    } else { # shift to the right
      frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
      dist_peak <- frac_part*(2*pi/omega) # distance from first peak
      phase_hours <- abs(dist_peak)
    }
    # should be no negative shifts
    # if (phase_hours < 0){ # only output positive shifts
    #   phase_hours <- phase_hours + (2*pi/omega)
    # }

    # calculate p-value
    ref_wave <- (alt_form(a,gam,omega,phi,y_shift,times)) # fitted values

    # if (num_reps > 1 | (num_reps == 1 & use_exact)){ # multiple replicates
      # calculate pvalue by exact JTK distribution
      taulist <- calc_tau_and_p_jtk(ref_wave,times,num_reps,current_gene,jtklist)
      tau <- taulist$tau # Kendall's tau
      pval <- taulist$pval # pvalue
    # }
    # else { # one replicate
    #   tau <- calc_kendalls_tau(as.numeric(ref_wave),as.numeric(y_val)) # calculate Kendall's tau
    #   pval <- calc_p_value(24,tau, TRUE) # calculate pvalue
    # }

    # list of parameters and other resulting values
    if (num_reps == 1){
      results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
      results <- cbind(results, y_val, rbind(ref_wave))
    } else{
      results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
    }

    return (results)
  }, error = function(e){ # if there's failure in convergence

    if (num_reps == 1){
      results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, y_val, rbind(rep(NA,length(times))))
    } else{
      results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(times))))
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

# function for determining whether gene values are unexpressed (less than 70% expressed (i.e., not 0)) for full matrix
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean if there is 70% expression
genes_unexpressed_all <- function(rem_unexpr_amt){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  # get how many genes are expressed for each gene
  tot_expressed <- rowSums(all_reps != 0,na.rm = TRUE)

  # return false if amount is less than threshold
  return(tot_expressed <= (ncol(all_reps)*rem_unexpr_amt))
}


# function for determining whether gene values are deviating or constant (one replicate)
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean if there is deviation within the gene values
is_deviating <- function(current_gene){
  stdev <- sd(genes[current_gene,-1], na.rm = TRUE)
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
    left_na <- cbind(matrix(0,nrow(all_reps),1),mtx_count[[i]][,-ncol(mtx_count[[i]])]) # left shifted matrix
    right_na <- cbind(mtx_count[[i]][,-1],matrix(0,nrow(genes),1)) # right shifted matrix
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
    for (i in 1:num_reps){
      each_rep <- all_rep[,seq(i,ncol(all_rep),by=num_reps)]

      # covariance
      cov <- rowSums((each_rep-rowMeans(each_rep,na.rm = TRUE))*(xmtx-rowMeans(xmtx)),na.rm = TRUE)
      # variance
      var <- rowSums((xmtx - rowMeans(xmtx))^2,na.rm = TRUE)

      beta <- cov/var
      alph <- rowMeans(each_rep,na.rm = TRUE)-(beta*rowMeans(xmtx))

      df[,seq(i,ncol(all_rep),by=num_reps)] <- each_rep -(alph+(beta*xmtx)) # linear fit
    }
  }
  # get the data frame correctly set up for returning
  res_df <- genes
  res_df[,-1] <- df

  return (res_df)
}
