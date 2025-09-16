prevCI.hypgeom.simul <- function(N, n, n_pos, Se, Sp,
                                 conf.level = 0.95, alt = "two.sided",
                                 R = 1000){
  # Description
  # This function calculates approximate confidence interval (CI) from simulation
  # for the prevalence, assuming sampling without replacement from a
  # finite population, that is, using the hypergeometrical model, and
  # allowing for misclassification with known sensitivity and specificity.
  # Also provides point estimate of true prevalence.
  
  # Arguments
  # N: size of the population
  # n: sample size
  # n_pos: number of observed positives in the sample
  # Se: sensitivity of the diagnostic test (no uncertainty assumed)
  # Sp: specificity of the diagnostic test (no uncertainty assumed)
  # conf.level: prescribed confidence level, defaults to 0.95
  # alt: alternative - "two.sided" (default), "less", "greater"
  # R: number of simulations, defaults to 1000
  
  # Details
  # If a level (1 - alpha) one-sided CI is needed, 
  # a two-sided CI with (1 - 2*alpha) is constructed
  # and (0, UCL) or (LCL, 1) given for "less" or "greater", respectively.
  
  # Value
  # MLE:              The value of truly infected in the population
  #                   that maximizes the likelihood of observing the
  #                   actual n_pos test positives, accounting for
  #                   misclassification (sensitivity and specificity).
  # RoganGladen:      Rogan-Gladen point estimate.
  # ClopperPearson:   Lower and upper bounds of the exact
  #                   Clopper-Pearson confidence interval
  # Blaker:           If alt is "two.sided" the lower and upper
  #                   bounds of Blaker's confidence interval are also given.
  
  alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
  
  # error handling
  if (!is.numeric(N) | !is.numeric(n) | !is.numeric(n_pos)) {stop("Population and sample data must be numeric.")}
  if (N <= 0) {stop("Population size must be positive.")}
  if (n <= 0 | n > N) {stop("Sample size (n) must be in (0, N].")}
  if (n_pos < 0 | n_pos > n) {stop("The number of positives (n_pos) must be in [0, n].")}
  if (!is.numeric(Se) | !is.numeric(Sp)) {stop("Se and Sp must be numeric.")}
  if (Se < 0.5 | Se > 1 | Sp < 0.5 | Sp > 1) {stop("Se and Sp must be in [0.5, 1].")}
  if (conf.level <= 0 | conf.level >= 1) {stop("Confidence level must be in (0, 1).")}
  
  # set alpha, thresholds and safe initialization
  # for alt "two.sided" alpha will actually be alpha/2 for clearer coding
  if (alt == "two.sided") {
    alpha <- (1 - conf.level)/2
  } else {
    alpha <- (1 - conf.level)
  }
  threshold1 <- floor(R*alpha)
  threshold2 <- ceiling(R*(1 - alpha))
  LCL <- 0
  UCL <- N
  
  # for Blaker acceptance region
  accept <- logical(N + 1)
  
  # for MLE
  likelihoods <- numeric(N + 1)
  
  # CI construction for a misclassified hypergeometric
  for (Nsick in 0:N) {
    # simulate true positives in sample
    hyp <- rhyper(R, Nsick, N - Nsick, n)
    # vectorized misclassification
    false_negatives <- rbinom(R, hyp, 1 - Se)
    false_positives <- rbinom(R, n - hyp, 1 - Sp)
    observed_pos <- hyp + false_positives - false_negatives
    
    # Likelihood for MLE
    likelihoods[Nsick + 1] <- mean(observed_pos == n_pos)
    
    # Blaker
    if (alt == "two.sided") {
      probs <- prop.table(table(observed_pos))
      cum <- cumsum(probs)
      revcum <- rev(cumsum(rev(probs)))
      p_lower <- mean(observed_pos <= n_pos)
      p_upper <- mean(observed_pos >= n_pos)
      p1 <- min(p_lower, p_upper)
      if (p_lower < p_upper) {
        p2 <- max(revcum[revcum <= p1], 0)
      } else {
        p2 <- max(cum[cum <= p1], 0)
      }
      accept[Nsick + 1] <- ((p1 + p2) > 2*alpha)
    }
    
    # Clopper-Pearson and exit
    # tails of probability distribution
    lower <- sort(observed_pos)[threshold1]
    upper <- sort(observed_pos)[threshold2]
    # find Clopper-Pearson limits
    if (upper < n_pos) {LCL <- Nsick + 1}
    if (lower <= n_pos) {UCL <- Nsick}
    if (lower > n_pos) break
  }
  
  # maximum likelihood estimate
  Nsick_MLE <- which.max(likelihoods) - 1
  prev_MLE <- Nsick_MLE/N
  
  # Rogan-Gladen point estimate
  RG <- (n_pos/n + Sp - 1)/(Se + Sp - 1)
  RG <- min(max(RG, 0), 1)
  
  # determine CIs
  if (alt == "two.sided") {
    CI_CP <- c(LCL, UCL)/N
    if (any(accept)) {
      Nsick_accept <- which(accept) - 1
      CI_B <- c(min(Nsick_accept), max(Nsick_accept))/N
    } else {
      CI_B <- CI_CP
    }
  } else if (alt == "less") {
    CI_CP <- c(0, UCL)/N
  } else if (alt == "greater") {
    CI_CP <- c(LCL, N)/N
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

