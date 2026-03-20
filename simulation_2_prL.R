library(dplyr)
library(parallel)
library(data.table)
source("prevCI_hypgeom_prLik.r")

evaluate_prevCI_hyp_prLik <- function(N, n, true_prev, Se, n_Se, Sp, n_Sp, B = 5000) {
  true_sick <- round(N*true_prev)
  
  contains_true_prLik <- numeric(B)
  widths_prLik        <- numeric(B)
  loc_prLik           <- numeric(B)
  RoganGladen         <- numeric(B)
  
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
    
    # generate Se and Sp evaluation samples
    k_Se <- rbinom(1, n_Se, Se)
    k_Sp <- rbinom(1, n_Sp, Sp)
    
    CI_prLik <- prevCI.hypgeom.prLik(N, n, observed_pos, n_Se, k_Se, n_Sp, k_Sp)
    
    contains_true_prLik[b] <- (true_prev >= CI_prLik$ProfLikCI[1]) &
                              (true_prev <= CI_prLik$ProfLikCI[2])
    widths_prLik[b] <- CI_prLik$ProfLikCI[2] - CI_prLik$ProfLikCI[1]
    loc_prLik[b]    <- CI_prLik$MLE
    RoganGladen[b]  <- CI_prLik$RoganGladen
  }
  
  data.frame(
    N = N,
    true_prev = true_prev,
    n = n,
    Se = Se,
    n_Se = n_Se,
    Sp = Sp,
    n_Sp = n_Sp,
    repl = B,
    coverage_prLik_mean = mean(contains_true_prLik),
    width_prLik_mean = mean(widths_prLik),
    width_prLik_sd = sd(widths_prLik),
    loc_prLik_mean = mean(loc_prLik),
    loc_prLik_sd = sd(loc_prLik),
    loc_RG_mean = mean(RoganGladen),
    loc_RG_sd = sd(RoganGladen)
  )
}

run_sim_grid_parallel <- function(
    N = 100,
    prev_seq = seq(0.01, 0.5, length = 2),
    n_seq = c(10, 30),
    Se_seq = c(0.9, 0.99),
    n_Se_seq = c(100),
    Sp_seq = c(0.9, 0.99),
    n_Sp_seq = c(100),
    B = 5000,
    cores = detectCores() - 1
) {
  # Create all (n, prevalence) combinations
  grid <- expand.grid(
    true_prev = prev_seq,
    n = n_seq,
    Se = Se_seq,
    n_Se = n_Se_seq,
    Sp = Sp_seq,
    n_Sp = n_Sp_seq
  )
  
  # Function wrapper for parallel call
  wrapper <- function(params) {
    evaluate_prevCI_hyp_prLik(
      N = N,
      n = params["n"],
      true_prev = params["true_prev"],
      Se = params["Se"],
      n_Se = params["n_Se"],
      Sp = params["Sp"],
      n_Sp = params["n_Sp"],
      B = B
    )
  }
  
  cl <- makeCluster(min(detectCores() - 2, nrow(grid)))
  clusterExport(cl, varlist = c("evaluate_prevCI_hyp_prLik",
                                "prevCI.hypgeom.prLik"), envir = environment())
  results <- parLapply(cl, split(grid, seq(nrow(grid))), function(row) wrapper(unlist(row)))
  stopCluster(cl)
  
  sim_data <- do.call(rbind, results)
  return(sim_data)
}

N50 <- run_sim_grid_parallel(
  N = 50,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(5, 10, 25),
  Se_seq = c(0.8, 0.9, 0.99),
  n_Se_seq = c(50, 100, 1000),
  Sp_seq = c(0.8, 0.9, 0.99),
  n_Sp_seq = c(50, 100, 1000),
  B = 5e3
)
N100 <- run_sim_grid_parallel(
  N = 100,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(10, 20, 50),
  Se_seq = c(0.8, 0.9, 0.99),
  n_Se_seq = c(50, 100, 1000),
  Sp_seq = c(0.8, 0.9, 0.99),
  n_Sp_seq = c(50, 100, 1000),
  B = 5e3
)
N200 <- run_sim_grid_parallel(
  N = 200,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(20, 40),
  Se_seq = c(0.8, 0.9, 0.99),
  n_Se_seq = c(50, 100, 1000),
  Sp_seq = c(0.8, 0.9, 0.99),
  n_Sp_seq = c(50, 100, 1000),
  B = 5e3
)
N200n100S99 <- run_sim_grid_parallel(
  N = 200,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(100),
  Se_seq = c(0.99),
  n_Se_seq = c(50, 100, 1000),
  Sp_seq = c(0.99),
  n_Sp_seq = c(50, 100, 1000),
  B = 5e3
)
N200n100S9 <- run_sim_grid_parallel(
  N = 200,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(100),
  Se_seq = c(0.9),
  n_Se_seq = c(50, 100, 1000),
  Sp_seq = c(0.9),
  n_Sp_seq = c(50, 100, 1000),
  B = 5e3
)
N200n100S8 <- run_sim_grid_parallel(
  N = 200,
  prev_seq = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  n_seq = c(100),
  Se_seq = c(0.8),
  n_Se_seq = c(50, 100, 1000),
  Sp_seq = c(0.8),
  n_Sp_seq = c(50, 100, 1000),
  B = 5e3
)

sim_data_2_prL <- bind_rows(
  N50, N100, N200,
  N200n100S99,
  N200n100S9,
  N200n100S8
)
fwrite(sim_data_2_prL, "sim_data_2_prL.csv")
