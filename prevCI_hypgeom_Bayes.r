prevCI.hypgeom.Bayes <- function(N, n, n_pos, Se, Sp, cred.level = 0.95, alt = "two.sided",
                                 prior = rep(1, N + 1)) {
  # Description
  # Calculates approximate Bayesian credible interval (CI)
  # for the true prevalence, assuming sampling without replacement from a
  # finite population, that is, using the hypergeometrical model, and
  # allowing for misclassification with known sensitivity and specificity.
  # Also provides point estimate of true prevalence.
  
  # Arguments
  # N: size of the population
  # n: sample size
  # n_pos: number of observed positives in the sample;
  #        a noisy version of the true positives in the sample
  # Se: sensitivity of the diagnostic test (no uncertainty assumed)
  # Sp: specificity of the diagnostic test (no uncertainty assumed)
  # cred.level: credible level, defaults to 0.95
  # alt: alternative - "two.sided" (default), "less", "greater"
  # prior: prior over number of truly sick in population (default: flat)
  
  # Details
  # If a level (1 - alpha) one-sided CI is needed, 
  # a two-sided CI with (1 - 2*alpha) is constructed
  # and (0, UCL) or (LCL, 1) given for "less" or "greater", respectively.
  # Use non-uniform priors (e.g. Beta-binomial) if prior knowledge is available, e.g.
  # beta_binomial_prior <- function(N, alpha, beta) {
  #   k <- 0:N
  #   prior <- choose(N, k)*beta(k + alpha, N - k + beta)/beta(alpha, beta)
  #   return(prior/sum(prior))
  #   }
  
  # Value
  # MAP:              The Bayesian posterior mode estimate
  #                   (maximum a posteriori estimate - MAP) of prevalence.
  # RoganGladen:      Rogan-Gladen point estimate.
  # CredibleInterval: Lower (LCL) and upper (UCL) bounds of the
  #                   central credible interval.
  
  require(ggplot2)
  
  alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
  
  # error handling
  if (!is.numeric(N) | !is.numeric(n) | !is.numeric(n_pos)) {stop("Population and sample data must be numeric.")}
  if (N <= 0) {stop("Population size must be provided.")}
  if (n <= 0 | n > N) {stop("Sample size (n) must be between 0 and N.")}
  if (n_pos < 0 | n_pos > n) {stop("The number of positives (n_pos) must be between 0 and n.")}
  if (!is.numeric(Se) | !is.numeric(Sp)) {stop("Se and Sp must be numeric.")}
  if (Se < 0.5 | Se > 1 | Sp < 0.5 | Sp > 1) {stop("Se and Sep must be between 0.5 and 1.")}
  if (cred.level <= 0 | cred.level >= 1) {stop("Credible level must be between 0 and 1.")}
  if (length(prior) != N + 1) {stop("Length of prior must be N + 1.")}
  
  # Rogan-Gladen point estimate
  RG <- (n_pos/n + Sp - 1)/(Se + Sp - 1)
  
  # ensure prior sums to 1
  prior <- prior/sum(prior)
  
  # Likelihood function (vectorized)
  Nsick_vals <- 0:N
  likelihood <- sapply(Nsick_vals, function(Nsick) {
    p_truesick <- dhyper(0:n, Nsick, N - Nsick, n)
    obs_pos_prob <- sapply(0:n, function(ts) {
      sum(dbinom(0:ts, ts, 1 - Se)*
          dbinom(n_pos - (ts - 0:ts), n - ts, 1 - Sp)*
          (n_pos - (ts - 0:ts) >= 0)*
          (n_pos - (ts - 0:ts) <= (n - ts))
      )
    })
    sum(p_truesick*obs_pos_prob)
  })
  
  # posterior over Nsick
  unnorm_posterior <- prior*likelihood
  posterior <- unnorm_posterior/sum(unnorm_posterior)
  
  # map Nsick to prevalence
  prevalence_values <- 0:N/N
  
  # posterior mode (maximum a posteriori estimate)
  map_index <- which.max(posterior)
  prev_MAP <- prevalence_values[map_index]
  
  # compute central credible interval
  if (alt == "two.sided") {
    alpha <- (1 - cred.level)/2
  } else {
    alpha <- (1 - cred.level)
  }
  cum_p <- cumsum(posterior)
  lower_idx <- which(cum_p >= alpha)[1]
  upper_idx <- which(cum_p >= 1 - alpha)[1]
  if (alt == "two.sided") {
    CI <- c(prevalence_values[lower_idx], prevalence_values[upper_idx])
  } else if (alt == "less") {
    CI <- c(0, prevalence_values[upper_idx])
  } else if (alt == "greater") {
    CI <- c(prevalence_values[lower_idx], 1)
  }
  
  # plot posterior
  df <- data.frame(
    prevalence = prevalence_values,
    posterior = posterior
  )
  
  plot <- ggplot(df, aes(x = prevalence, y = posterior)) +
    geom_line(linewidth = 1.1) +
    geom_area(data = subset(df, prevalence >= CI[1] & prevalence <= CI[2]),
              aes(y = posterior), fill = "steelblue", alpha = 0.5) +
    geom_vline(xintercept = CI, linetype = "dashed", color = "red") +
    geom_vline(xintercept = prev_MAP, linetype = "dotted", color = "black") +
    annotate("text", x = CI, y = 0.9*max(posterior), label = round(CI, 2),
             col = "darkgrey", hjust = c(1, 0)) +
    annotate("text", x = prev_MAP, y = 0.95*max(posterior),
             label = paste("MAP =", round(prev_MAP, 3)),
             vjust = -1, col = "black", size = 3.5) +
    labs(title = "Posterior Distribution of Prevalence",
         subtitle = paste0(cred.level*100, "% Credible Interval \"", alt, "\""),
         x = "Prevalence",
         y = "Posterior Density") +
    theme_minimal()
  
  cat("MAP\n")
  print(prev_MAP)
  cat("\nRoganGladen\n")
  print(min(max(RG, 0), 1))
  cat("\nCredibleInterval\n")
  print(CI)
  print(plot)
  invisible(list(MAP = prev_MAP, RoganGladen = min(max(RG, 0), 1),
                 CredibleInterval = CI, plot = plot, posterior = df))
}

# example
prevCI.hypgeom.Bayes(50, 25, 10, 0.9, 0.9)
