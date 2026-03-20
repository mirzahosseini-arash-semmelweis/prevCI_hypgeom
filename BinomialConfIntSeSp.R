BinomialConfIntSeSp <- function(n, y, conf.level, alt, met, Se, Sp){
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

# Exmaple run
BinomialConfIntSeSp(
  n = 50,
  y = 25,
  conf.level = 0.95,
  alt = "two.sided",
  met = "Blaker",
  Se = 0.9,
  Sp = 0.99
)