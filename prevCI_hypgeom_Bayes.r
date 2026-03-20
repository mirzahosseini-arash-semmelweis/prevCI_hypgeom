prevCI.hypgeom.Bayes <- function(
    N, n, n_pos,
    Se = NULL, Sp = NULL,
    Se_study = NULL,
    Sp_study = NULL,
    Se_prior_ab = c(1, 1),
    Sp_prior_ab = c(1, 1),
    M = 100,
    cred.level = 0.95,
    prior = rep(1, N + 1),
    seed = NULL,
    show.plot = TRUE
) {
  # Description
  # Bayesian credible interval (HDI) for prevalence in a finite population
  # (hypergeometric sampling without replacement) with misclassification.
  #
  # Extension: If Se_study and Sp_study are provided, Se and Sp are treated as
  # uncertain and integrated out via Monte Carlo using Beta posteriors derived
  # from the independent studies.
  
  # Arguments
  # N: size of the population
  # n: sample size
  # n_pos: number of observed positives in the sample;
  #        a noisy version of the true positives in the sample
  # Either provide
  # Se: sensitivity of the diagnostic test (no uncertainty assumed)
  # Sp: specificity of the diagnostic test (no uncertainty assumed)
  # Or provide independent study results for Se/Sp uncertainty:
  # Se_study: c(k_tp, n_diseased)
  # Sp_study: c(k_tn, n_nondiseased)
  # Se_prior_ab: prior for Se ~ Beta(a,b)
  # Sp_prior_ab: prior for Sp ~ Beta(a,b)
  # M: Monte Carlo draws for marginalizing Se/Sp
  # cred.level: credible level, defaults to 0.95
  # prior: prior over number of truly sick in population (default: flat)
  # seed: for Monde Carlo simulation
  # show.plot: defaults to TRUE
  
  # Details
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
  #                   highest density credible interval.
  # Se posterior / Sp posterior: If Se_study / Sp_study were given.
  
  require(ggplot2)
  
  # error handling
  if (!is.numeric(N) | !is.numeric(n) | !is.numeric(n_pos)) stop("Population and sample data must be numeric.")
  if (N <= 0) stop("Population size must be positive.")
  if (n <= 0 | n > N) stop("Sample size (n) must be in (0, N].")
  if (n_pos < 0 | n_pos > n) stop("The number of positives (n_pos) must be in [0, n].")
  if (cred.level <= 0 | cred.level >= 1) stop("Credible level must be in (0, 1).")
  if (length(prior) != N + 1) stop("Length of prior must be N + 1.")
  if (!is.null(seed)) set.seed(seed)
  
  # determine mode: fixed Se/Sp vs uncertain from studies
  using_studies <- (!is.null(Se_study) || !is.null(Sp_study))
  
  if (using_studies) {
    if (is.null(Se_study) || is.null(Sp_study)) {
      stop("If using studies, provide BOTH Se_study and Sp_study.")
    }
    if (!is.numeric(Se_study) || length(Se_study) != 2) stop("Se_study must be numeric length-2: c(k_tp, n_diseased).")
    if (!is.numeric(Sp_study) || length(Sp_study) != 2) stop("Sp_study must be numeric length-2: c(k_tn, n_nondiseased).")
    kSe <- Se_study[1]; nSe <- Se_study[2]
    kSp <- Sp_study[1]; nSp <- Sp_study[2]
    if (nSe <= 0 || nSp <= 0) stop("Study totals must be positive.")
    if (kSe < 0 || kSe > nSe) stop("Se_study successes must be in [0, total].")
    if (kSp < 0 || kSp > nSp) stop("Sp_study successes must be in [0, total].")
    if (any(Se_prior_ab <= 0) || length(Se_prior_ab) != 2) stop("Se_prior_ab must be positive length-2.")
    if (any(Sp_prior_ab <= 0) || length(Sp_prior_ab) != 2) stop("Sp_prior_ab must be positive length-2.")
    
    aSe <- Se_prior_ab[1] + kSe
    bSe <- Se_prior_ab[2] + (nSe - kSe)
    aSp <- Sp_prior_ab[1] + kSp
    bSp <- Sp_prior_ab[2] + (nSp - kSp)
    
    # pre-draw Se/Sp once (re-used for every Nsick)
    Se_draws <- rbeta(M, aSe, bSe)
    Sp_draws <- rbeta(M, aSp, bSp)
    Se_draws <- pmin(pmax(Se_draws, 1e-6), 1 - 1e-6)
    Sp_draws <- pmin(pmax(Sp_draws, 1e-6), 1 - 1e-6)
    
    # for display / RG point estimate, posterior means
    Se_hat <- aSe / (aSe + bSe)
    Sp_hat <- aSp / (aSp + bSp)
  } else {
    # fixed Se/Sp path
    if (is.null(Se) || is.null(Sp)) stop("Provide either fixed Se and Sp, or Se_study and Sp_study.")
    if (!is.numeric(Se) | !is.numeric(Sp)) stop("Se and Sp must be numeric.")
    if (Se < 0.5 | Se > 1 | Sp < 0.5 | Sp > 1) stop("Se and Sp must be in [0.5, 1].")
    Se_hat <- Se
    Sp_hat <- Sp
  }
  
  # Rogan-Gladen point estimate (using Se_hat/Sp_hat)
  denom <- (Se_hat + Sp_hat - 1)
  if (abs(denom) < 1e-12) {
    RG <- NA_real_
    warning("Se_hat + Sp_hat is ~ 1; Rogan-Gladen is ill-defined.")
  } else {
    RG <- (n_pos / n + Sp_hat - 1) / denom
    RG <- min(max(RG, 0), 1)
  }
  
  # ensure prior sums to 1
  prior <- prior / sum(prior)
  
  # core likelihood helper:
  # p(n_pos | Nsick, Se, Sp) with hypergeom + misclass convolution
  lik_one <- function(Nsick, Se_val, Sp_val) {
    # distribution of true sick in sample: ts = 0..n
    p_truesick <- dhyper(0:n, Nsick, N - Nsick, n)
    
    # for each possible true sick count ts, compute P(obs positives = n_pos | ts, Se, Sp)
    # obs positives = TP + FP; TP ~ Binom(ts, Se); FP ~ Binom(n-ts, 1-Sp)
    # convolution: sum_{tp} P(TP=tp) P(FP = n_pos - tp)
    obs_pos_prob <- vapply(0:n, function(ts) {
      tp_vals <- 0:ts
      fp_needed <- n_pos - tp_vals
      ok <- (fp_needed >= 0) & (fp_needed <= (n - ts))
      if (!any(ok)) return(0)
      
      sum(
        dbinom(tp_vals[ok], size = ts, prob = Se_val) *
          dbinom(fp_needed[ok], size = (n - ts), prob = (1 - Sp_val))
      )
    }, numeric(1))
    
    sum(p_truesick * obs_pos_prob)
  }
  
  # likelihood over Nsick (0..N)
  Nsick_vals <- 0:N
  
  if (using_studies) {
    # Monte Carlo marginalization over Se,Sp
    likelihood <- vapply(Nsick_vals, function(Nsick) {
      # average over draws
      mean(vapply(seq_len(M), function(m) lik_one(Nsick, Se_draws[m], Sp_draws[m]), numeric(1)))
    }, numeric(1))
  } else {
    likelihood <- vapply(Nsick_vals, function(Nsick) lik_one(Nsick, Se_hat, Sp_hat), numeric(1))
  }
  
  # posterior over Nsick
  unnorm_posterior <- prior * likelihood
  if (!is.finite(sum(unnorm_posterior)) || sum(unnorm_posterior) <= 0) {
    stop("Posterior normalization failed (all zero/NA). Check inputs.")
  }
  posterior <- unnorm_posterior / sum(unnorm_posterior)
  
  # map to prevalence
  prevalence_values <- (0:N) / N
  
  # MAP
  map_index <- which.max(posterior)
  prev_MAP <- prevalence_values[map_index]
  
  # HDI (highest density interval) on discrete grid
  ord <- order(posterior, decreasing = TRUE)
  cs <- cumsum(posterior[ord])
  k <- which(cs >= cred.level)[1]
  thr <- posterior[ord][k]
  idx <- which(posterior >= thr - 1e-15)
  CI_hdi <- c(prevalence_values[min(idx)], prevalence_values[max(idx)])
  
  # plotting
  df <- data.frame(prevalence = prevalence_values, posterior = posterior)
  
  plot <- ggplot(df, aes(x = prevalence, y = posterior)) +
    geom_line(linewidth = 1.1) +
    geom_area(
      data = subset(df, prevalence >= CI_hdi[1] & prevalence <= CI_hdi[2]),
      aes(y = posterior), fill = "steelblue", alpha = 0.5
    ) +
    geom_vline(xintercept = CI_hdi, linetype = "dashed", color = "red") +
    geom_vline(xintercept = prev_MAP, linetype = "dotted", color = "black") +
    annotate("text", x = CI_hdi, y = 0.9*max(posterior), label = round(CI_hdi, 2),
             col = "darkgrey", hjust = c(1, 0)) +
    annotate("text", x = prev_MAP, y = 0.95*max(posterior),
             label = paste("MAP =", round(prev_MAP, 3)),
             vjust = -1, col = "black", size = 3.5) +
    labs(
      title = "Posterior Distribution of Prevalence",
      subtitle = paste0(round(cred.level * 100, 1), "% Highest Density Credible Interval"),
      x = "Prevalence",
      y = "Posterior Density"
    ) +
    theme_minimal()
  
  # output
  cat("MAP\n"); print(prev_MAP)
  cat("\nRoganGladen\n"); print(RG)
  cat("\nHDCredibleInterval\n"); print(CI_hdi)
  
  if (using_studies) {
    cat("\nSe posterior (Beta)\n")
    print(c(a = aSe, b = bSe, mean = Se_hat))
    cat("\nSp posterior (Beta)\n")
    print(c(a = aSp, b = bSp, mean = Sp_hat))
  }
  
  if (show.plot) print(plot)
  
  invisible(list(
    MAP = prev_MAP,
    RoganGladen = RG,
    HDCredibleInterval = CI_hdi,
    plot = plot,
    posterior = df,
    Se_hat = Se_hat,
    Sp_hat = Sp_hat,
    using_studies = using_studies
  ))
}

# Example run
prevCI.hypgeom.Bayes(
    N     = 100,
    n     = 50,
    n_pos = 25,
    Se = 0.9,
    Sp = 0.99,
    cred.level = 0.95,
    prior = rep(1, N + 1)
)
