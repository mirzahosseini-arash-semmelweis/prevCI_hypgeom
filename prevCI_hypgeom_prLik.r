prevCI.hypgeom.prLik <- function(N, n, n_pos,
                                 n_Se, k_Se, n_Sp, k_Sp,
                                 conf.level = 0.95, alt = "two.sided") {
  # Description
  # This function calculates profile likelihood confidence interval for prevalence 
  # when sensitivity and specificity are estimated from samples 
  # that are independent of the sample from which prevalence is estimated.
  # Boundary smoothing / continuity correction is applied.
  # For the prevalence estimation, sampling without replacement from a
  # finite population is assumed, that is, the hypergeometrical model is used.
  # The function also provides the Rogan-Gladen point estimate of true prevalence.
  
  # Arguments
  # N: size of the population
  # n: sample size
  # n_pos: number of observed positives in the sample
  # n_Se: sample size in the sensitivity validation study
  # k_Se: number of positives in the sensitivity validation study
  # n_Sp: sample size in the specificity validation study
  # k_Sp: number of positives in the specificity validation study
  # conf.level: prescribed confidence level, defaults to 0.95
  # alt: alternative - "two.sided" (default), "less", "greater"
  
  # Details
  # If a level (1 - alpha) one-sided CI is needed,
  # a two-sided CI with (1 - 2*alpha) is constructed
  # and (0, UCL) or (LCL, 1) given for "less" or "greater", respectively.
  # Binary search is used to narrow the relevant number of truly infected for efficiency.
  
  # Value
  # MLE:              The value of truly infected in the population
  #                   that maximizes the likelihood of observing the
  #                   actual n_pos test positives, accounting for
  #                   misclassification (sensitivity and specificity).
  # RoganGladen:      Rogan-Gladen point estimate.
  # ProfLikCI:        Lower and upper bounds of the profile likelihood
  #                   confidence interval.
  # ML.Se:            Joint maximum likelihood estimate of diagnostic sensitivity.
  # ML.Sp:            Joint maximum likelihood estimate of diagnostic specificity.
  #                   (jointly estimated with prevalence under the misclassification-adjusted hypergeometric model)
  
  alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
  
  # error handling
  if (!is.numeric(N) | !is.numeric(n) | !is.numeric(n_pos)) {stop("Population and sample data must be numeric.")}
  if (!is.numeric(n_Se) | !is.numeric(n_Sp) | !is.numeric(k_Se) | !is.numeric(k_Sp)) {stop("Diagnostic test data must be numeric.")}
  if (N <= 0) {stop("Population size must be positive.")}
  if (n <= 0 | n > N) {stop("Sample size (n) must be in (0, N].")}
  if (n_pos < 0 | n_pos > n) {stop("The number of positives (n_pos) must be in [0, n].")}
  if (k_Se < 0 | k_Se > n_Se) {stop("The number of positives in sensitivity study (k_Se) must be in [0, n_Se].")}
  if (k_Sp < 0 | k_Sp > n_Sp) {stop("The number of positives in specificity study (k_Sp) must be in [0, n_Sp].")}
  if (conf.level <= 0 | conf.level >= 1) {stop("Confidence level must be in (0, 1).")}
  stopifnot(N == as.integer(N),
            n == as.integer(n),
            n_pos == as.integer(n_pos),
            n_Se == as.integer(n_Se),
            k_Se == as.integer(k_Se),
            n_Sp == as.integer(n_Sp),
            k_Sp == as.integer(k_Sp))
  
  # set alpha; for alt "less" or "greater" alpha will actually be 2*signf.level for clearer coding
  if (alt == "two.sided") {
    alpha <- (1 - conf.level)
  } else {
    alpha <- 2*(1 - conf.level)
  }
  
  ##############################################################################
  # helper function for optim
  optim_with_fallback <- function(
    fn, par,
    gr = NULL,
    control = list(),
    max_retries_per_method = 2,
    jitter_scale = 0.1,
    verbose = FALSE,
    ...
  ) {
    stopifnot(is.numeric(par), length(par) >= 1)
    
    # small helper: run optim safely and check convergence
    run_optim <- function(method, par0) {
      if (verbose) message("Trying method = ", method, " ...")
      res <- tryCatch(
        optim(par = par0, fn = fn, gr = if (method == "BFGS") gr else NULL,
              method = method, control = control, ...),
        error = function(e) e
      )
      # on error, return list with error
      if (inherits(res, "error")) {
        if (verbose) message("  -> error: ", conditionMessage(res))
        return(list(ok = FALSE, res = NULL, err = conditionMessage(res)))
      }
      # check convergence: code 0 = success (per ?optim)
      if (!is.null(res$convergence) && res$convergence == 0 && is.finite(res$value)) {
        if (verbose) message("  -> success (convergence = 0)")
        return(list(ok = TRUE, res = res, err = NULL))
      } else {
        if (verbose) {
          message(sprintf("  -> non-convergence (code = %s, value = %s)",
                          as.character(res$convergence), as.character(res$value)))
        }
        return(list(ok = FALSE, res = res, err = sprintf("non-convergence: %s", res$convergence)))
      }
    }
    
    # jitter generator for retries
    jitter_par <- function(p, scale) {
      p + rnorm(length(p), sd = pmax(1e-8, abs(p)*scale + 1e-8))
    }
    
    methods <- c("Nelder-Mead", "BFGS", "SANN")
    last_err <- NULL
    tried <- list()
    
    for (meth in methods) {
      # try initial par first
      attempt <- run_optim(meth, par)
      tried[[meth]] <- attempt
      if (attempt$ok) {
        attempt$res$method_used <- meth
        return(attempt$res)
      }
      last_err <- attempt$err
      
      # retry with jittered starts (useful for NM/SANN, sometimes for BFGS)
      for (k in seq_len(max_retries_per_method)) {
        par_try <- jitter_par(par, jitter_scale)
        attempt <- run_optim(meth, par_try)
        tried[[paste0(meth, "_retry", k)]] <- attempt
        if (attempt$ok) {
          attempt$res$method_used <- meth
          attempt$res$par_start <- par_try
          return(attempt$res)
        }
        last_err <- attempt$err
      }
    }
    
    # if here, all methods failed: return the "best" among failed attempts if available
    if (verbose) message("All methods failed; returning best nonconverged result if any.")
    best <- NULL; best_val <- Inf; best_method <- NA
    for (nm in names(tried)) {
      att <- tried[[nm]]
      if (!is.null(att$res) && is.list(att$res) && is.finite(att$res$value)) {
        if (att$res$value < best_val) {
          best <- att$res; best_val <- att$res$value; best_method <- nm
        }
      }
    }
    if (!is.null(best)) {
      best$method_used <- paste0(best_method, " (nonconverged)")
      best$warning <- paste("All methods failed to converge. Last error:", last_err)
      return(best)
    }
    
    stop("optim_with_fallback: all methods failed; last error: ", last_err)
  }
  
  # helper function for extension of hypergeometric distribution to continuous values
  dhypgamma <- function(N, K, n, k, log = FALSE) {
    # basic checks
    if (N <= 0 || n < 0 || n > N) stop("Require 0 <= n <= N, N > 0.")
    if (any(K < 0 | K > N)) stop("Require 0 <= K <= N.")
    
    # support bounds for k
    lo <- pmax(0, n - (N - K))
    hi <- pmin(n, K)
    inside <- (k >= lo) & (k <= hi)
    
    # allocate output
    out <- rep(if (log) -Inf else 0, length(k))
    
    # log binomial coefficient via lgamma
    lC <- function(a, b) {
      lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)
    }
    
    # compute log pmf
    logv <- lC(K, k) + lC(N - K, n - k) - lC(N, n)
    
    if (log) {
      out[inside] <- logv[inside]
    } else {
      out[inside] <- exp(logv[inside])
    }
    out
  }
  
  # helper function to compute distribution of positives (pos_dist) in sample
  # conditional on true sick in population (Nsick)
  get_pos_dist <- function(Nsick, Se, Sp, from, to) {
    # Support of X: true infected in sample
    x_vals <- max(0, n - (N - round(Nsick))):min(n, round(Nsick))
    # Precompute Hypergeometric pmf for X
    hx <- dhypgamma(N, Nsick, n, x_vals)
    
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
    pos_dist[pos_dist <= 0] <- .Machine$double.eps
    return(pos_dist[(from + 1):(to + 1)])
  }
  ##############################################################################
  
  # main function
  compute_npos <- function(n_pos_local) {
    # observed relative frequencies
    prev_obs <- n_pos_local/n
    Se_obs <- k_Se/n_Se
    Sp_obs <- k_Sp/n_Sp
    
    # Rogan-Gladen point estimate of true prevalence
    RG <- (prev_obs + Sp_obs - 1)/(Se_obs + Sp_obs - 1)
    RG <- min(max(0, RG), 1)
    
    # Likelihood function for prev, Se, and Sp in logit scale
    logLik.s <- function(theta) {
      p_s <- plogis(theta[1])
      Se_s <- plogis(theta[2])
      Sp_s <- plogis(theta[3])
      Nsick_s <- p_s*N
      logLik_s <-
        log(get_pos_dist(Nsick_s, Se_s, Sp_s, from = n_pos_local, to = n_pos_local)) +
        dbinom(k_Se, n_Se, Se_s, log = TRUE) +
        dbinom(k_Sp, n_Sp, Sp_s, log = TRUE)
      return(-logLik_s)
    }
    
    # parameter initialization and logLikelihood maximization over all parameters
    init_prev <- qlogis(min(max(1e-6, RG), 1 - 1e-6))
    init_Se <- qlogis(min(max(1e-6, Se_obs), 1 - 1e-6))
    init_Sp <- qlogis(min(max(1e-6, Sp_obs), 1 - 1e-6))
    logL_s <- optim_with_fallback(par = c(init_prev, init_Se, init_Sp), fn = logLik.s)
    logL_allparams <- plogis(logL_s$par)
    logL_par <- logL_allparams[1]
    
    # Profile likelihood function
    profLik.p <- function(Nsick_i) {
      logLik.0 <- function(theta) {
        Se_0 <- plogis(theta[1])
        Sp_0 <- plogis(theta[2])
        logLik_0 <-
          log(get_pos_dist(Nsick_i, Se_0, Sp_0, from = n_pos_local, to = n_pos_local)) +
          dbinom(k_Se, n_Se, Se_0, log = TRUE) +
          dbinom(k_Sp, n_Sp, Sp_0, log = TRUE)
        return(-logLik_0)
      }
      logL_0 <- optim_with_fallback(par = c(init_Se, init_Sp), fn = logLik.0)
      return(logL_0$value)
    }
    
    # set critical value and threshold for LR test
    crit.val <- qchisq(1 - alpha, 1)
    thresh <- logL_s$value + crit.val/2
    
    # binary search for LCL
    low <- 0
    high <- ceiling(logL_par*N)
    if (profLik.p(low) < thresh) {
      LCL <- low
    } else {
      while ((high - low) > 1) {
        mid <- floor((low + high)/2)
        if (profLik.p(mid) < thresh) {
          high <- mid
        } else {
          low <- mid
        }
      }
      LCL <- high
    }
    
    # binary search for UCL
    low <- floor(logL_par*N)
    high <- N
    if (profLik.p(high) < thresh) {
      UCL <- high
    } else {
      while ((high - low) > 1) {
        mid <- ceiling((low + high)/2)
        if (profLik.p(mid) < thresh) {
          low <- mid
        } else {
          high <- mid
        }
      }
      UCL <- low
    }
    
    # set confidence interval
    if (alt == "two.sided") {
      CI_prLik <- c(LCL, UCL)/N
    } else if (alt == "less") {
      CI_prLik <- c(0, UCL)/N
    } else if (alt == "greater") {
      CI_prLik <- c(LCL, N)/N
    }
    
    # result
    list(
      MLE = logL_allparams[1],
      RoganGladen = RG,
      ProfLikCI = CI_prLik,
      ML.Se = logL_allparams[2],
      ML.Sp = logL_allparams[3],
      LCL_count = LCL,
      UCL_count = UCL
    )
  }
  
  # compute "base" result at the observed n_pos
  base <- compute_npos(n_pos)
  
  # boundary adjustment for CI
  if (n_pos == 0) {
    adj <- compute_npos(1)
    LCL_adj <- base$LCL_count
    UCL_adj <- as.integer(floor(mean(c(adj$UCL_count, base$UCL_count))))
    
    if (alt == "two.sided") {
      base$ProfLikCI <- c(LCL_adj, UCL_adj)/N
    } else if (alt == "less") {
      base$ProfLikCI <- c(0, UCL_adj)/N
    } else {
      base$ProfLikCI <- c(LCL_adj, N)/N
    }
  }
  if (n_pos == n) {
    adj <- compute_npos(n - 1)
    LCL_adj <- as.integer(ceiling(mean(c(adj$LCL_count, base$LCL_count))))
    UCL_adj <- base$UCL_count
    
    if (alt == "two.sided") {
      base$ProfLikCI <- c(LCL_adj, UCL_adj)/N
    } else if (alt == "less") {
      base$ProfLikCI <- c(0, UCL_adj)/N
    } else {
      base$ProfLikCI <- c(LCL_adj, N)/N
    }
  }
  
  # result
  result <- list(
    MLE = base$MLE,
    RoganGladen = base$RoganGladen,
    ProfLikCI = base$ProfLikCI,
    ML.Se = base$ML.Se,
    ML.Sp = base$ML.Sp
  )
  return(result)
}

# Example run
prevCI.hypgeom.prLik(
  N     = 100,
  n     = 50,
  n_pos = 25,
  n_Se  = 100,
  k_Se  = 90,
  n_Sp  = 100,
  k_Sp  = 99,
  conf.level = 0.95,
  alt = "two.sided"
)