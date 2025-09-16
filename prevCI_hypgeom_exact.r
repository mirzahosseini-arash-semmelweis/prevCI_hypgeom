prevCI.hypgeom.exact <- function(N, n, n_pos, Se, Sp,
                                 conf.level = 0.95, alt = "two.sided") {
  # Description
  # This function calculates exact Clopper-Pearson and Blaker type confidence
  # intervals (CIs) for the true prevalence, assuming sampling without
  # replacement from a finite population, that is, using the hypergeometrical
  # model, and allowing for misclassification with known sensitivity and
  # specificity. Also provides point estimate of true prevalence.
  
  # Arguments
  # N: size of the population
  # n: sample size
  # n_pos: number of observed positives in the sample
  # Se: sensitivity of the diagnostic test (no uncertainty assumed)
  # Sp: specificity of the diagnostic test (no uncertainty assumed)
  # conf.level: prescribed confidence level, defaults to 0.95
  # alt: alternative - "two.sided" (default), "less", "greater"
  
  # Details
  # If a level (1 - alpha) one-sided CI is needed,
  # a two-sided CI with (1 - 2*alpha) is constructed
  # and (0, UCL) or (LCL, 1) given for "less" or "greater", respectively.
  # Binary search is used to narrow the relevant number of truly infected for efficiency.
  # Search for Clopper-Pearson CI starts from exact binomial CI estimates for speed.
  
  # Value
  # MLE:              The value of truly infected in the population
  #                   that maximizes the likelihood of observing the
  #                   actual n_pos test positives, accounting for
  #                   misclassification (sensitivity and specificity).
  # RoganGladen:      Rogan-Gladen point estimate.
  # ClopperPearson:   Lower and upper bounds of the exact
  #                   Clopper-Pearson confidence interval.
  # Blaker:           If alt is "two.sided" the lower and upper bounds
  #                   of Blaker's confidence interval are also given.
  
  alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
  
  # error handling
  if (!is.numeric(N) | !is.numeric(n) | !is.numeric(n_pos)) {stop("Population and sample data must be numeric.")}
  if (N <= 0) {stop("Population size must be positive.")}
  if (n <= 0 | n > N) {stop("Sample size (n) must be in (0, N].")}
  if (n_pos < 0 | n_pos > n) {stop("The number of positives (n_pos) must be in [0, n].")}
  if (!is.numeric(Se) | !is.numeric(Sp)) {stop("Se and Sp must be numeric.")}
  if (Se < 0.5 | Se > 1 | Sp < 0.5 | Sp > 1) {stop("Se and Sp must be in [0.5, 1].")}
  if (conf.level <= 0 | conf.level >= 1) {stop("Confidence level must be in (0, 1).")}
  
  # set alpha; for alt "two.sided" alpha will actually be alpha/2 for clearer coding
  if (alt == "two.sided") {
    alpha <- (1 - conf.level)/2
  } else {
    alpha <- (1 - conf.level)
  }
  
  # Rogan-Gladen point estimate
  RG <- (n_pos/n + Sp - 1)/(Se + Sp - 1)
  RG <- min(max(RG, 0), 1)
  
  # helper function to compute distribution of positives (pos_dist) in sample
  # conditional on true sick in population (Nsick)
  get_pos_dist <- function(Nsick, from, to) {
    # Support of X: true infected in sample
    x_vals <- max(0, n - (N - Nsick)):min(n, Nsick)
    # Precompute Hypergeometric pmf for X
    hx <- dhyper(x_vals, m = Nsick, n = N - Nsick, k = n)
    
    # Support of Y: observed positives
    y_vals <- 0:n
    pos_dist <- numeric(length(y_vals))
    
    # Convolution: sum_x P(X=x) * P(Y=y | X=x)
    for (i in seq_along(x_vals)) {
      x <- x_vals[i]
      # Y | X=x = Bin(x, Se) + Bin(n-x, 1-Sp)
      y1 <- 0:x
      y2 <- 0:(n - x)
      # Convolve the two binomials by summing matching pairs y1+y2 = y
      p1 <- dbinom(y1, size = x, prob = Se)
      p2 <- dbinom(y2, size = n - x, prob = 1 - Sp)
      conv <- convolve(p1, rev(p2), type = "open")  # length x+(n-x)+1 = n+1
      pos_dist <- pos_dist + hx[i] * conv
    }
    
    names(pos_dist) <- as.character(y_vals)
    return(pos_dist[(from + 1):(to + 1)])
  }
  
  # Likelihood: search only in vicinity of Rogan-Gladen estimate
  likelihoods <- numeric(N + 1)
  # look for MLE from (RG estimate - 2) to (RG estimate + 2)
  low <- max(floor(N*RG) - 2, 0)
  high <- min(ceiling(N*RG) + 2, N)
  for (Nsick in low:high) {
    likelihoods[Nsick + 1] <- get_pos_dist(Nsick, from = n_pos, to = n_pos)
  }
  current_max <- max(likelihoods)
  # if max not on the edge then accept as MLE
  if (max(likelihoods[low + 1], likelihoods[high + 1]) != current_max) {
    Nsick_MLE <- which.max(likelihoods) - 1
  } else {
    # if max on extremes of population then accept as MLE
    if (which.max(likelihoods) %in% c(1, N + 1)) {
      Nsick_MLE <- which.max(likelihoods) - 1
    } else {
      # if max on the edge but not on population extremes, keep looking for MLE
      low_index <- low + 1
      high_index <- high + 1
      low_value <- likelihoods[low_index]
      high_value <- likelihoods[high_index]
      direction <- sign(high_value - low_value) # direction to look further
      current_max_index <- which.max(likelihoods)
      next_inner_index <- current_max_index - direction
      next_inner_value <- likelihoods[next_inner_index]
      while (current_max > next_inner_value) {
        next_inner_index <- next_inner_index + direction
        next_inner_value <- current_max
        current_max_index <- current_max_index + direction
        current_max <- get_pos_dist(current_max_index - 1, from = n_pos, to = n_pos)
      }
      # current_max is no longer the max, accept next_inner_value as max
      Nsick_MLE <- next_inner_index - 1
    }
  }
  
  # MLE
  prev_MLE <- Nsick_MLE/N
  
  # Clopper-Pearson acceptability function
  CP_accept <- function(Nsick, from, to) {
    pos_dist <- get_pos_dist(Nsick, from, to)
    p_sum <- sum(pos_dist)
    
    # acceptability criterion
    return(p_sum > alpha)
  }
  
  # find exact binomial Clopper-Pearson CI estimates for faster binary search
  binomCI_lower <- 0
  binomCI_upper <- 1
  if (n_pos != 0) {
    binomCI_lower <- qbeta(alpha, n_pos, n - n_pos + 1)
  }
  if (n_pos != n) {
    binomCI_upper <- qbeta(1 - alpha, n_pos + 1, n - n_pos)
  }
  binomCI <- c(binomCI_lower, binomCI_upper)
  binomCI_SeSp <- (binomCI + Sp - 1)/(Se + Sp - 1)
  binomCI_SeSp[1] = min(1, max(0, binomCI_SeSp[1]))
  binomCI_SeSp[2] = min(1, max(0, binomCI_SeSp[2]))
  
  # binary search for LCL
  low <- max(floor(N*binomCI_SeSp[1]) - 1, 0)
  high <- min(Nsick_MLE + 1, N)
  if (CP_accept(low, from = n_pos, to = n)) {
    LCL <- low
  } else {
    while ((high - low) > 1) {
      mid <- floor((low + high)/2)
      if (CP_accept(mid, from = n_pos, to = n)) {
        high <- mid
      } else {
        low <- mid
      }
    }
    LCL <- high
  }
  
  # binary search for UCL
  low <- max(Nsick_MLE - 1, 0)
  high <- min(ceiling(N*binomCI_SeSp[2]) + 1, N)
  if (CP_accept(high, from = 0, to = n_pos)) {
    UCL <- high
  } else {
    while ((high - low) > 1) {
      mid <- ceiling((low + high)/2)
      if (CP_accept(mid, from = 0, to = n_pos)) {
        low <- mid
      } else {
        high <- mid
      }
    }
    UCL <- low
  }
  
  # Clopper-Pearson CIs
  if (alt == "two.sided") {
    CI_CP <- c(LCL, UCL)/N
  } else if (alt == "less") {
    CI_CP <- c(0, UCL)/N
  } else if (alt == "greater") {
    CI_CP <- c(LCL, N)/N
  }
  
  # Blaker acceptability function
  if (alt == "two.sided" & LCL == UCL) {CI_B = CI_CP}
  if (alt == "two.sided" & LCL != UCL) {
    blaker_accept <- function(Nsick) {
      pos_dist <- get_pos_dist(Nsick, from = 0, to = n)
      p_lower <- sum(pos_dist[as.integer(names(pos_dist)) <= n_pos])
      p_upper <- sum(pos_dist[as.integer(names(pos_dist)) >= n_pos])
      
      # p1: sum of probabilities from n_obs and closer tail
      p1 <- min(p_lower, p_upper)
      
      # p2: sum of all probabilities < p1 in the opposite tail
      cum <- cumsum(pos_dist)
      revcum <- rev(cumsum(rev(pos_dist)))
      if (p1 == p_lower) {
        p2_index <- which(revcum <= p1)
        if (length(p2_index) == 0) {p2 <- 0} else {p2 <- revcum[min(p2_index)]}
      } else {
        p2_index <- which(cum <= p1)
        if (length(p2_index) == 0) {p2 <- 0} else {p2 <- cum[max(p2_index)]}
      }
      
      # acceptability criterion
      return((p1 + p2) > 2*alpha)
    }
    
    # Blaker lower bound: start at LCL, move up
    for (Nsick in LCL:Nsick_MLE) {
      if (blaker_accept(Nsick)) {
        Blaker_LCL <- Nsick
        break
      }
    }
    if (!exists("Blaker_LCL")) {
      Blaker_LCL <- LCL # if none of values are accepted, use C-P as CI bound
    }
    
    # Blaker upper bound: start at UCL, move down
    for (Nsick in UCL:Nsick_MLE) {
      if (blaker_accept(Nsick)) {
        Blaker_UCL <- Nsick
        break
      }
    }
    if (!exists("Blaker_UCL")) {
      Blaker_UCL <- UCL # if none of values are accepted, use C-P as CI bound
    }
    
    # Blaker CI
    CI_B <- c(Blaker_LCL, Blaker_UCL)/N
  }
  
  # result
  result <- list(
    MLE = prev_MLE,
    RoganGladen = RG,
    ClopperPearson = CI_CP
  )
  if (alt == "two.sided") {
    result$Blaker <- CI_B
  }
  
  return(result)
}

