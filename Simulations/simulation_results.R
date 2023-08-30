setwd("C:/Users/spammer/Desktop/M.Sc. in Economics/So 2023/panel_data_hmethods/Simulations")

source("sim_dgp_fn.R")
source("sim_montecarlo_fn.R")
source("sim_final_fn.R")


set.seed(42)


################## Simulation Results ############################

n <- c(30, 100, 300)
T <- c(12, 30)
beta <- c(0.5, 0.5)
DGPs <- c("K1", "K2", "K3", "K4")



################# Baseline ########################

# 1. Exogeneous regressors, Gaussian errors

all_results_b <- list()


for (T_val in T) {
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    
    
    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = FALSE, error = "homo")
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_b[[label]] <- results
    }
  }
}



# 2. Endogenous regressors, Gaussian erros




all_results_bT <- list()


for (T_val in T) {
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = TRUE, error = "homo")
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_bT[[label]] <- results
    }
  }
}

######################### Heteroskedasticty ########################################

# 3. Exogenous regressors, time-level heteroskedasticity

all_results_hetero_time <- list()


for (T_val in T) {
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = FALSE, error = "hetero_time", error.type = 4)
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_hetero_time[[label]] <- results
    }
  }
}

# 4. Exogenous regressors, individual-level heteroskedasticity 


all_results_hetero_individual<- list()


for (T_val in T) {
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    

    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = FALSE, error = "hetero_individual", error.type = 2)
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_hetero_individual[[label]] <- results
    }
  }
}

# 5. Exogenous regressors, individual-level and time-level heteroskedastictiy


all_results_hetero_both <- list()


for (T_val in T) { 
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    

    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = FALSE, error = "hetero_both", error.type = 6)
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_hetero_both[[label]] <- results
    }
  }
}

###################################### Autocorrelation ################################################

# 6. Exogenous regressors, low autocorrelation in erros 

all_results_autocorr_low <- list()

for (T_val in T) {
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    
    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = FALSE, error = "autocorr_low", error.type = 5 )
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_autocorr_low[[label]] <- results
    }
  }
}


# 7. Exogenous regressors, high autocorrelation in erros


all_results_autocorr_high <- list()

for (T_val in T) {
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    
    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = FALSE, error = "autocorr_high", error.type = 5 )
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_autocorr_high[[label]] <- results
    }
  }
}

############################### Robustness #################################################

all_results_KSSdim <- list()

DGPs <- c("K3")

for (T_val in T) {
  for (n_val in n) {
    nsim_val <- ifelse(n_val == 300, 500, 1000)
    for (DGP_type in DGPs) {
      results <- doMonteCarlo(nsim_val, T_val, n_val, beta, DGP_type, endogenous = FALSE, error = "homo", true.dim = TRUE)
      label <- paste0("T", T_val, "_n", n_val, "_nsim", nsim_val, "_", DGP_type)
      all_results_KSSdim[[label]] <- results
    }
  }
}


############################# Tables ######################################################

DGPs <- c("K1", "K2", "K3", "K4")
scenarios <- list(
  b = all_results_b,
  bT = all_results_bT,
  hetero_time = all_results_hetero_time,
  hetero_individual = all_results_hetero_individual,
  hetero_both = all_results_hetero_both,
  autocorr_low = all_results_autocorr_low,
  autocorr_high = all_results_autocorr_high,
  KSSdim = all_results_KSSdim
)

for (scenario_name in names(scenarios)) {
  for (DGP in DGPs) {
    all_means <- calculate_means_for_paths_and_estimators(scenarios[[scenario_name]], DGP)
    generate_combined_latex_table(all_means, scenario_name, DGP)
  }
}





