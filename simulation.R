library(parallel)

################################################################################
# Function from Reiczigel et al. 2010, https://doi.org/10.1017/S0950268810000385

BinomialConfIntSeSp <- function(n, y, conf.level = 0.95, alt, met, Se, Sp){
  # Calculates a confidence interval for the true prevalence, assuming that
  # the sensitivity and specificity of the diagnostic test are known constants.
  # Based on Reiczigel et al. (2010), Epidemiol. Infect. 138, 1674–1678.
  
  # Args:
  #   n: The number of individuals that are tested.
  #   y: Number of individuals with positive test results.
  #   conf.level: Confidence level.
  #   alt: alternative - "two.sided", "less", "greater".
  #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"
  #   Se: sensitivity of the diagnostic procedure (a number between 0 and 1)
  #   Sp: specificity of the diagnostic procedure (a number between 0 and 1)
  #
  # Returns:
  #   A vector containing the LCL and UCL of the binomial CI.
  

	BinomialConfInt <- function(n, y, conf.level = 0.95, alt = "two.sided", met = "Clopper-Pearson") {
	  # Function to calculate a binomial confidence interval with several different methods
	  #  
	  # Args:
	  #   n: The number of individuals that are tested.
	  #   y: Number of individuals with positive test results.
	  #   conf.level: Confidence level.
	  #   alt: alternative - "two.sided", "less", "greater".
	  #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"
	  #
	  # Returns:
	  #   A vector containing the LCL and UCL of the binomial CI.
	  
	  p_obs <- y / n
	  alpha <- 1 - conf.level
	  z1s <- qnorm(1 - alpha)
	  z2s <- qnorm(1 - alpha / 2)
	  cilower <- 0
	  ciupper <- 1
	  
	  alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
	  met <- match.arg(met, choices = c("Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"))
	  
	  # Error message
	  if (p_obs < 0 | p_obs > 1) stop("The number of positives (y) must be between 0 and n.")
	  
	  # Methods
	  if (met == "Wald") {
		if (alt == "two.sided") {
		  CI <- c(max(p_obs - z2s * sqrt(p_obs * (1 - p_obs)/(n)), 0), min(p_obs + z2s * sqrt(p_obs * (1 - p_obs)/(n)), 1))
		} else {
		  if (alt == "less") {
			CI <- c(0, min(p_obs + z1s * sqrt(p_obs * (1 - p_obs)/(n)), 1))
		  } else {
			if (alt == "greater") {
			  CI = c(max(p_obs - z1s * sqrt(p_obs * (1 - p_obs)/(n)), 0), 1)
			} else {
			  stop("Argument alt must be specified!")
			}
		  }
		}
	  }
	  
	  
	  if (met == "Wilson") {
		if (alt == "two.sided") {
		  CI <- c((2 * n * p_obs + z2s^2 - z2s * sqrt(z2s^2 + 4 * n * p_obs * (1 - p_obs)))/(2 * (n + z2s^2)),
				  (2 * n * p_obs + z2s^2 + z2s * sqrt(z2s^2 + 4 * n * p_obs * (1-p_obs)))/(2 * (n + z2s^2)))
		} else {
		  if (alt == "less") {
			CI <- c(0, (2 * n * p_obs + z1s^2 + z1s * sqrt(z1s^2 + 4 * n * p_obs * (1 - p_obs)))/(2 * (n + z1s^2)))
		  } else {
			if (alt == "greater") {
			  CI <- c((2 * n * p_obs + z1s^2 - z1s * sqrt(z1s^2 + 4 * n * p_obs * (1-p_obs)))/(2 * (n + z1s^2)), 1)
			} else {
			  stop("Argument alt must be specified!")
			}
		  }
		}
	  }
	  
	  
	  if (met == "Agresti-Coull") {
		phat1s = (p_obs * n + (z1s^2) / 2) / (n + z1s^2)
		phat2s = (p_obs * n + (z2s^2) / 2) / (n + z2s^2)
		nhat1s = n + z1s^2
		nhat2s = n + z2s^2
		if (alt == "two.sided") {
		  CI <- c(max(phat2s - z2s * sqrt(phat2s * (1 - phat2s)/(nhat2s)), 0),
				  min(phat2s + z2s * sqrt(phat2s * (1 - phat2s)/(nhat2s)), 1))
		} else {
		  if (alt == "less") {
			CI <- c(0, min(phat1s + z1s * sqrt(phat1s * (1 - phat1s)/(nhat1s)), 1))
		  } else {
			if (alt == "greater") {
			  CI <- c(max(0, phat1s - z1s * sqrt(phat1s * (1 - phat1s)/(nhat1s))), 1)
			} else {
			  stop("Argument alt must be specified!")
			}
		  }
		}
	  }
	  
	  
	  if (met == "Clopper-Pearson") {
		if (alt == "two.sided") {
		  if (y != 0) {
			cilower <- qbeta((1 - conf.level)/2, y, n - y + 1)
		  }
		  if (y != n) {
			ciupper <- qbeta(1 - (1 - conf.level)/2, y + 1, n - y)
		  }
		}
		if (alt == "less") {
		  if (y != n) {
			ciupper <- qbeta(1 - (1 - conf.level), y + 1, n - y)
		  }
		}
		if (alt == "greater") {
		  if (y != 0) {
			cilower <- qbeta((1 - conf.level), y, n - y + 1)
		  }
		}
		CI <- c(cilower, ciupper)
	  }
	  
	  
	  if (met == "Blaker") {
		step.value = 1e-04
		blakeraccept <- function(y, n, p) {
		  p1 = 1 - pbinom(y - 1, n, p)
		  p2 = pbinom(y, n, p)
		  a1 = p1 + pbinom(qbinom(p1, n, p) - 1, n, p)
		  a2 = p2 + 1 - pbinom(qbinom(1 - p2, n, p), n, p)
		  return(min(a1, a2))
		}
		if (alt == "two.sided") {
		  if (y != 0) {
			cilower <- qbeta((1 - conf.level)/2, y, n - (y) + 1)
			{
			  while (blakeraccept(y, n, cilower + step.value) < (1 - 
					conf.level)) cilower = cilower + step.value
			}
		  }
		  if (y != n) {
			ciupper <- qbeta(1 - (1 - conf.level)/2, (y) + 1, n - (y))
			{
			  while (blakeraccept(y, n, ciupper - step.value) < (1 - 
					conf.level)) ciupper = ciupper - step.value
			}
		  }
		}
		if (alt == "less") {
		  if (y != n) {
			ciupper <- qbeta(1 - (1 - conf.level), (y) + 1, n - (y))
			{
			  while (blakeraccept(y, n, ciupper - step.value) < (1 - 
					conf.level)) ciupper = ciupper - step.value
			}
		  }
		}
		if (alt == "greater") {
		  if (y != 0) {
			cilower <- qbeta((1 - conf.level), y, n - (y) + 1)
			{
			  while (blakeraccept(y, n, cilower + step.value) < (1 - 
					conf.level)) cilower = cilower + step.value
			}
		  }
		}
		CI <- c(cilower, ciupper)
	  }
	  return(CI)
	}


	# Calculating a CI for the apparent prevalence
	CI.observed <- BinomialConfInt(n, y, conf.level, alt, met)
	
	# Tranforming its endpoints
	CI.true <- (CI.observed+Sp-1)/(Se+Sp-1)
	CI.true[1] = min(1,max(0,CI.true[1]))
	CI.true[2] = min(1,max(0,CI.true[2]))
	
	return(CI.true)
}

################################################################################
# Functions from Ge et al. 2024, https://doi.org/10.1080/00031305.2023.2250401
##############################################################################
#                                                                            #
#                                                                            #
#  Self-defined functions used for "Enhanced Inference for Finite Population #
#  Sampling-Based Prevalence Estimation with Misclassification Errors"       #
#                                                                            #
#                                                                            #
##############################################################################

## A function for the proposed Bayesian Credible interval approach

RS_BC = function(N, n.RS, n_pos.RS, Se, Sp){
  
  ## Parameters
  ## N: total population size (N=1, output the BC for prevalence estimator;
  ##    N=Ntot, output the BC for case count estimator)
  ## n_pos.RS: number of positive test results from the imperfect test
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  # set the number of posterior samples for BC interval
  m = 1000
  
  # calculate the estimate of test positive frequency ("pi" in paper) and the 
  # bias-corrected disease prevalence ("pi_c" in paper) with threshold [0, 1]
  p_star_RS = min(max(n_pos.RS/n.RS,0.0001), 0.9999)     # 1>= p_star_RS >= 0
  p_RS = min(max((p_star_RS+Sp-1)/(Se+Sp-1),0.0001), 0.9999)  # 1>= p_RS >= 0
  
  # V1(pi) in eqn.(1)
  V_p_star_RS = p_star_RS*(1-p_star_RS)/n.RS
  
  # V3(pi) in eqn.(12)
  fpc = n.RS*(N-n.RS)/(N*(n.RS-1)) 
  V_p_RS = fpc*V_p_star_RS+(p_RS*Se*(1-Se) + (1-p_RS)*Sp*(1-Sp))/N
  
  # calculate the scale and shift parameter a and b
  a=sqrt(V_p_RS/V_p_star_RS)
  b=p_star_RS*(1-a)
  
  # generate beta posterior samples and adjust by a and b
  p_star_post = rbeta(m,n_pos.RS+1/2,n.RS-n_pos.RS+1/2)
  p_true_post = a*p_star_post + b
  p_true_post = (p_true_post+Sp-1)/(Se+Sp-1)
  
  # calculate the BC interval
  # Ntot=1, we calculate the BC for prevalence estimator
  # Ntot=N, we calculate the BC for case count estimator
  Ntot = 1
  N_iter = Ntot*p_true_post
  
  N_iter_lower = max(quantile(N_iter,c(0.025)),0)
  N_iter_upper = max(quantile(N_iter,c(0.975)),0)
  
  N_iter_interval_width = N_iter_upper - N_iter_lower
  
  N_iter_median = max(quantile(N_iter,c(0.5)),0)
  
  return(list(BC_lower = N_iter_lower,BC_upper = N_iter_upper,BC_width = N_iter_interval_width,
              BC_median = N_iter_median))
}

RS_BC(N = 50, n_pos.RS = 10, n.RS = 25, Se = 0.9, Sp = 0.9)

## A function for the proposed FPC Wald interval approach

RS_Wald <- function(N, n.RS, n_pos.RS, Se, Sp) {
  pi_hat = n_pos.RS/n.RS #estimate of test positive frequency
  pi_hat_c = max((pi_hat + Sp - 1)/(Se + Sp - 1), 0) #bias-corrected disease prevalence
  
  # standard errors of pi_hat and pi_hat_c
  V1_pi = pi_hat*(1 - pi_hat)/n.RS
  fpc2 = n.RS*(N - n.RS)/(N*(n.RS - 1)) 
  V2_pi = V1_pi*fpc2
  V_pi_extra = (pi_hat_c*Se*(1 - Se) + (1 - pi_hat_c)*Sp*(1 - Sp))/N
  V3_pi = V2_pi + V_pi_extra
  se = 1/(Se + Sp - 1)*sqrt(V3_pi)
  
  # bias-corrected, FPC Wald-type CI with misclassification
  CI <- c(pi_hat_c - 1.96*se, pi_hat_c + 1.96*se)
  return(CI)
}

################################################################################
# Functions for simulation

evaluate_prevCI_hypgeom <- function(N, n, true_prev, Se, Sp, B = 1000) {
  true_sick <- round(N*true_prev)
  
  contains_true_exact_CP <- numeric(B)
  widths_exact_CP <- numeric(B)
  contains_true_exact_B <- numeric(B)
  widths_exact_B <- numeric(B)
  loc_exact <- numeric(B)
  contains_true_simul_CP <- numeric(B)
  widths_simul_CP <- numeric(B)
  contains_true_simul_B <- numeric(B)
  widths_simul_B <- numeric(B)
  loc_simul <- numeric(B)
  contains_true_Bayes <- numeric(B)
  widths_Bayes <- numeric(B)
  loc_Bayes <- numeric(B)
  RoganGladen <- numeric(B)
  contains_true_FPC <- numeric(B)
  widths_FPC <- numeric(B)
  contains_true_Ge <- numeric(B)
  widths_Ge <- numeric(B)
  loc_Ge <- numeric(B)
  contains_true_binom_CP <- numeric(B)
  widths_binom_CP <- numeric(B)
  contains_true_binom_B <- numeric(B)
  widths_binom_B <- numeric(B)
  
  for (b in 1:B) {
    # uniform random sampling of size n is taken from population of size N without replacement
    sample_ids <- sample(N, n, replace = FALSE)
    # assuming the first true_sick members of the population are sick,
    # the number of truly sick in sample is calculated
    sampled_sick <- sum(sample_ids <= true_sick)
    
    # calculating positives in sample considering misclassification
    TP <- rbinom(1, sampled_sick, Se)
    FP <- rbinom(1, n - sampled_sick, 1 - Sp)
    observed_pos <- TP + FP
    
    CI_exact <- prevCI.hypgeom.exact(N, n, observed_pos, Se, Sp)
    CI_simul <- prevCI.hypgeom.simul(N, n, observed_pos, Se, Sp)
    CI_Bayes <- prevCI.hypgeom.Bayes(N, n, observed_pos, Se, Sp)
    CI_FPC <- RS_Wald(N, n, observed_pos, Se, Sp)
    CI_Ge <- RS_BC(N, n, observed_pos, Se, Sp)
    CI_binom_CP <- BinomialConfIntSeSp(n = n, y = observed_pos, alt = "two.sided", met = "Clopper-Pearson", Se = Se, Sp = Sp)
    CI_binom_B <- BinomialConfIntSeSp(n = n, y = observed_pos, alt = "two.sided", met = "Blaker", Se = Se, Sp = Sp)
    
    contains_true_exact_CP[b] <- (true_prev >= CI_exact$ClopperPearson[1]) &
                                 (true_prev <= CI_exact$ClopperPearson[2])
    widths_exact_CP[b] <- CI_exact$ClopperPearson[2] - CI_exact$ClopperPearson[1]
    contains_true_exact_B[b] <- (true_prev >= CI_exact$Blaker[1]) &
                                (true_prev <= CI_exact$Blaker[2])
    widths_exact_B[b] <- CI_exact$Blaker[2] - CI_exact$Blaker[1]
    loc_exact[b] <- CI_exact$MLE
    contains_true_simul_CP[b] <- (true_prev >= CI_simul$ClopperPearson[1]) &
                                 (true_prev <= CI_simul$ClopperPearson[2])
    widths_simul_CP[b] <- CI_simul$ClopperPearson[2] - CI_simul$ClopperPearson[1]
    contains_true_simul_B[b] <- (true_prev >= CI_simul$Blaker[1]) &
                                (true_prev <= CI_simul$Blaker[2])
    widths_simul_B[b] <- CI_simul$Blaker[2] - CI_simul$Blaker[1]
    loc_simul[b] <- CI_simul$MLE
    contains_true_Bayes[b] <- (true_prev >= CI_Bayes$CredibleInterval[1]) &
                              (true_prev <= CI_Bayes$CredibleInterval[2])
    widths_Bayes[b] <- CI_Bayes$CredibleInterval[2] - CI_Bayes$CredibleInterval[1]
    loc_Bayes[b] <- CI_Bayes$MAP
    
    RoganGladen[b] <- min(max((observed_pos/n + Sp - 1)/(Se + Sp - 1), 0), 1)
    
    contains_true_FPC[b] <- (true_prev >= CI_FPC[1]) & (true_prev <= CI_FPC[2])
    widths_FPC[b] <- CI_FPC[2] - CI_FPC[1]
    contains_true_Ge[b] <- (true_prev >= CI_Ge$BC_lower) & (true_prev <= CI_Ge$BC_upper)
    widths_Ge[b] <- CI_Ge$BC_width
    loc_Ge[b] <- CI_Ge$BC_median
    
    contains_true_binom_CP[b] <- (true_prev >= CI_binom_CP[1]) & (true_prev <= CI_binom_CP[2])
    widths_binom_CP[b] <- CI_binom_CP[2] - CI_binom_CP[1]
    contains_true_binom_B[b] <- (true_prev >= CI_binom_B[1]) & (true_prev <= CI_binom_B[2])
    widths_binom_B[b] <- CI_binom_B[2] - CI_binom_B[1]
  }
  
  data.frame(
    N = N,
    true_prev = true_prev,
    n = n,
    Se = Se,
    Sp = Sp,
    repl = B,
    coverage_exact_CP_mean = mean(contains_true_exact_CP),
    width_exact_CP_mean = mean(widths_exact_CP),
    width_exact_CP_sd = sd(widths_exact_CP),
    coverage_exact_B_mean = mean(contains_true_exact_B),
    width_exact_B_mean = mean(widths_exact_B),
    width_exact_B_sd = sd(widths_exact_B),
    loc_exact_mean = mean(loc_exact),
    loc_exact_sd = sd(loc_exact),
    coverage_simul_CP_mean = mean(contains_true_simul_CP),
    width_simul_CP_mean = mean(widths_simul_CP),
    width_simul_CP_sd = sd(widths_simul_CP),
    coverage_simul_B_mean = mean(contains_true_simul_B),
    width_simul_B_mean = mean(widths_simul_B),
    width_simul_B_sd = sd(widths_simul_B),
    loc_simul_mean = mean(loc_simul),
    loc_simul_sd = sd(loc_simul),
    coverage_Bayes_mean = mean(contains_true_Bayes),
    width_Bayes_mean = mean(widths_Bayes),
    width_Bayes_sd = sd(widths_Bayes),
    loc_Bayes_mean = mean(loc_Bayes),
    loc_Bayes_sd = sd(loc_Bayes),
    loc_RG_mean = mean(RoganGladen),
    loc_RG_sd = sd(RoganGladen),
    coverage_FPC_mean = mean(contains_true_FPC),
    width_FPC_mean = mean(widths_FPC),
    width_FPC_sd = sd(widths_FPC),
    coverage_Ge_mean = mean(contains_true_Ge),
    width_Ge_mean = mean(widths_Ge),
    width_Ge_sd = sd(widths_Ge),
    loc_Ge_mean = mean(loc_Ge),
    loc_Ge_sd = sd(loc_Ge),
    coverage_binom_CP_mean = mean(contains_true_binom_CP),
    width_binom_CP_mean = mean(widths_binom_CP),
    width_binom_CP_sd = sd(widths_binom_CP),
    coverage_binom_B_mean = mean(contains_true_binom_B),
    width_binom_B_mean = mean(widths_binom_B),
    width_binom_B_sd = sd(widths_binom_B)
  )
}

################################################################################
# Parallelized data simulation and evaluation

run_sim_grid_parallel <- function(
    N = 300,
    prev_seq = seq(0.01, 0.5, length = 2),
    n_seq = c(10, 30),
    Se_seq = c(0.9, 0.99),
    Sp_seq = c(0.9, 0.99),
    B = 500,
    cores = detectCores() - 1
) {
  # Create all (n, prevalence) combinations
  grid <- expand.grid(
    true_prev = prev_seq,
    n = n_seq,
    Se = Se_seq,
    Sp = Sp_seq
  )
  
  # Function wrapper for parallel call
  wrapper <- function(params) {
    evaluate_prevCI_hypgeom(
      N = N,
      n = params["n"],
      true_prev = params["true_prev"],
      Se = params["Se"],
      Sp = params["Sp"],
      B = B
    )
  }
  
  cl <- makeCluster(min(detectCores() - 1, nrow(grid)))
  clusterExport(cl, varlist = c("evaluate_prevCI_hypgeom",
                                "prevCI.hypgeom.exact",
                                "prevCI.hypgeom.simul",
                                "prevCI.hypgeom.Bayes",
                                "RS_Wald",
                                "RS_BC",
                                "BinomialConfIntSeSp"), envir = environment())
  results <- parLapply(cl, split(grid, seq(nrow(grid))), function(row) wrapper(unlist(row)))
  stopCluster(cl)
  
  sim_data <- do.call(rbind, results)
  return(sim_data)

}
