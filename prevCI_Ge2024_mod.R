##############################################################################
#                                                                            #
#                                                                            #
#  Self-defined functions used for "Enhanced Inference for Finite Population #
#  Sampling-Based Prevalence Estimation with Misclassification Errors"       #
#                                                                            #
#                                                                            #
##############################################################################

## A function for the proposed Bayesian Credible interval approach

RS_BC_mod = function(N, n.RS, n_pos.RS, Se, Sp){
  
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
  
  N_iter_lower = min(max(quantile(N_iter,c(0.025)),0),1)
  N_iter_upper = min(max(quantile(N_iter,c(0.975)),0),1)
  
  N_iter_interval_width = N_iter_upper - N_iter_lower
  
  N_iter_median = min(max(quantile(N_iter,c(0.5)),0),1)
  
  return(list(BC_lower = N_iter_lower,BC_upper = N_iter_upper,BC_width = N_iter_interval_width,
              BC_median = N_iter_median))
}

# Example run
RS_BC_mod(N = 100, n.RS = 50, n_pos.RS = 25, Se = 0.9, Sp = 0.99)