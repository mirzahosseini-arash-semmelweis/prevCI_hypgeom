library(dplyr)
library(parallel)
library(data.table)
source("BinomialConfintSeSp.R")
source("prevCI_Ge2024.R")
source("prevCI_Ge2024_mod.R")
source("prevCI_hypgeom_exact.r")
source("prevCI_hypgeom_simul.r")
source("prevCI_hypgeom_Bayes.r")

evaluate_prevCI_hypgeom <- function(N, n, true_prev, Se, Sp, B = 5000) {
  true_sick <- round(N*true_prev)
  
  contains_true_exact_CP <- numeric(B)
  widths_exact_CP        <- numeric(B)
  contains_true_exact_B  <- numeric(B)
  widths_exact_B         <- numeric(B)
  loc_exact              <- numeric(B)
  contains_true_simul_CP <- numeric(B)
  widths_simul_CP        <- numeric(B)
  contains_true_simul_B  <- numeric(B)
  widths_simul_B         <- numeric(B)
  loc_simul              <- numeric(B)
  contains_true_Bayes    <- numeric(B)
  widths_Bayes           <- numeric(B)
  loc_Bayes              <- numeric(B)
  RoganGladen            <- numeric(B)
  contains_true_FPC      <- numeric(B)
  widths_FPC             <- numeric(B)
  contains_true_Ge       <- numeric(B)
  widths_Ge              <- numeric(B)
  loc_Ge                 <- numeric(B)
  contains_true_Ge_mod   <- numeric(B)
  widths_Ge_mod          <- numeric(B)
  loc_Ge_mod             <- numeric(B)
  contains_true_binom_CP <- numeric(B)
  widths_binom_CP        <- numeric(B)
  contains_true_binom_B  <- numeric(B)
  widths_binom_B         <- numeric(B)
  
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
    
    CI_exact    <- prevCI.hypgeom.exact(N, n, observed_pos, Se, Sp)
    CI_simul    <- prevCI.hypgeom.simul(N, n, observed_pos, Se, Sp)
    CI_Bayes    <- prevCI.hypgeom.Bayes(N, n, observed_pos, Se, Sp)
    CI_FPC      <- RS_Wald(N, n, observed_pos, Se, Sp)
    CI_Ge       <- RS_BC(N, n, observed_pos, Se, Sp)
    CI_Ge_mod   <- RS_BC_mod(N, n, observed_pos, Se, Sp)
    CI_binom_CP <- BinomialConfIntSeSp(
      n = n,
      y = observed_pos,
      Se = Se, Sp = Sp,
      alt = "two.sided", met = "Clopper-Pearson", conf.level = 0.95
    )
    CI_binom_B  <- BinomialConfIntSeSp(
      n = n,
      y = observed_pos,
      Se = Se, Sp = Sp,
      alt = "two.sided", met = "Blaker", conf.level = 0.95
    )
    
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
    contains_true_Bayes[b] <- (true_prev >= CI_Bayes$HDCredibleInterval[1]) &
                              (true_prev <= CI_Bayes$HDCredibleInterval[2])
    widths_Bayes[b] <- CI_Bayes$HDCredibleInterval[2] - CI_Bayes$HDCredibleInterval[1]
    loc_Bayes[b] <- CI_Bayes$MAP
    
    RoganGladen[b] <- min(max((observed_pos/n + Sp - 1)/(Se + Sp - 1), 0), 1)
    
    contains_true_FPC[b] <- (true_prev >= CI_FPC[1]) & (true_prev <= CI_FPC[2])
    widths_FPC[b] <- CI_FPC[2] - CI_FPC[1]
    contains_true_Ge[b] <- (true_prev >= CI_Ge$BC_lower) & (true_prev <= CI_Ge$BC_upper)
    widths_Ge[b] <- CI_Ge$BC_width
    loc_Ge[b] <- CI_Ge$BC_median
    contains_true_Ge_mod[b] <- (true_prev >= CI_Ge_mod$BC_lower) & (true_prev <= CI_Ge_mod$BC_upper)
    widths_Ge_mod[b] <- CI_Ge_mod$BC_width
    loc_Ge_mod[b] <- CI_Ge_mod$BC_median
    
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
    coverage_Ge_mod_mean = mean(contains_true_Ge_mod),
    width_Ge_mod_mean = mean(widths_Ge_mod),
    width_Ge_mod_sd = sd(widths_Ge_mod),
    loc_Ge_mod_mean = mean(loc_Ge_mod),
    loc_Ge_mod_sd = sd(loc_Ge_mod),
    coverage_binom_CP_mean = mean(contains_true_binom_CP),
    width_binom_CP_mean = mean(widths_binom_CP),
    width_binom_CP_sd = sd(widths_binom_CP),
    coverage_binom_B_mean = mean(contains_true_binom_B),
    width_binom_B_mean = mean(widths_binom_B),
    width_binom_B_sd = sd(widths_binom_B)
  )
}

run_sim_grid_parallel <- function(
    N = 100,
    prev_seq = seq(0.01, 0.5, length = 2),
    n_seq = c(10, 30),
    Se_seq = c(0.9, 0.99),
    Sp_seq = c(0.9, 0.99),
    B = 5000,
    cores = detectCores() - 1
) {
  library(parallel)
  
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
                                "RS_BC_mod",
                                "BinomialConfIntSeSp"), envir = environment())
  results <- parLapply(cl, split(grid, seq(nrow(grid))), function(row) wrapper(unlist(row)))
  stopCluster(cl)
  
  sim_data <- do.call(rbind, results)
  return(sim_data)
}

sim_data_N50 <- run_sim_grid_parallel(
  N = 50,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(5, 10, 25),
  Se_seq = c(0.8, 0.9, 0.95, 0.99),
  Sp_seq = c(0.8, 0.9, 0.95, 0.99),
  B = 5e3
)
sim_data_N100 <- run_sim_grid_parallel(
  N = 100,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(10, 20, 50),
  Se_seq = c(0.8, 0.9, 0.95, 0.99),
  Sp_seq = c(0.8, 0.9, 0.95, 0.99),
  B = 5e3
)
sim_data_N200 <- run_sim_grid_parallel(
  N = 200,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(20, 40, 100),
  Se_seq = c(0.8, 0.9, 0.95, 0.99),
  Sp_seq = c(0.8, 0.9, 0.95, 0.99),
  B = 5e3
)

sim_data_1_exact <- bind_rows(sim_data_N50, sim_data_N100, sim_data_N200)
fwrite(sim_data_1_exact, "sim_data_1_exact.csv")
