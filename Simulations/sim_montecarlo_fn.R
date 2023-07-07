library(doParallel)
library(foreach)
source("sim_dgp_fn.R")


extract_simulation_results <- function(sim_results) {
  #' Extracts simulation results into a structured format for analysis.
  #'
  #' @param sim_results A list containing simulation results.
  #'
  #' @return A list with structured simulation results for different estimation methods.
  #' The list contains sublists for each estimation method (KSS, Eup, Within),
  #' with each sublist containing vectors for different metrics.
  
  # Initialize results
  results <- list(
    KSS = list(),
    Eup = list(),
    Within = list()
  )
  
  # Define metric names
  metrics_all <- c("b1", "b2", "used.dim", "sd1", "sd2", "mse", "b1.bias", "b2.bias", "mse.effect", "mse.coeff")
  metrics_within <- c("b1", "b2", "sd1", "sd2", "mse", "b1.bias", "b2.bias", "mse.coeff")
  
  # Number of simulations
  nsim <- length(sim_results)
  
  # Loop over each simulation result
  for (i in 1:nsim) {
    
    # Loop over each type (KSS, Eup, Within)
    for (type in c("KSS", "Eup", "Within")) {
      
      # Select metrics based on the type
      if (type == "Within") {
        metrics <- metrics_within
      } else {
        metrics <- metrics_all
      }
      
      # Loop over each metric
      for (metric in metrics) {
        
        # Append metric to the respective vector within results list
        results[[type]][[metric]] <- c(results[[type]][[metric]], sim_results[[i]][[type]][[metric]])
        
      }
    }
  }
  
  # Return the structured results
  return(results)
}


doMonteCarlo <- function(nsim, T, n, beta, DGP = "K3", endogenous = FALSE, error = "homo", rho = 0.5) {
  #' Run Monte Carlo simulation for panel data models
  #'
  #' This function performs a Monte Carlo simulation to estimate panel data models using various methods.
  #'
  #' @param nsim The number of simulations to run.
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param beta The true parameter values.
  #' @param DGP The data generating process. Default is "K3".
  #' @param endogenous Logical indicating whether the data is endogenous. Default is FALSE.
  #' @param error The type of panel errors to generate. Supported values are "homo", "hetero", and "autocorr". Default is "homo".
  #' @param rho Level of error-autocorrelation for error type "autocorr". 
  #'
  #' @return A list of estimation results for different methods and performance metrics
  
  # Register the parallel backend
  no_cores <- detectCores() # detect number of cores
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  # Define package dependencies and functions to be exported
  required_packages <- c("mvtnorm", "matlib", "phtt", "plm", "doParallel", "foreach")
  exported_functions <- c("generate_X", "shift_mean", 
                          "generate_V1", "generate_V2", 
                          "generate_V3", "generate_V4", 
                          "generate_V5", "generate_V6",
                          "generate_panel_errors", "endogenous_X", 
                          "create_plm_data", "create_Tn_mat",
                          "DataGeneratingFunction")
  
  # Run the Monte Carlo simulation using parallel processing
  sim_results <- foreach(i = 1:nsim,
                         .inorder = FALSE,
                         .packages = required_packages,
                         .export = exported_functions) %dopar% {
    
    # Generate data
    data <- DataGeneratingFunction(T, n, beta, DGP, endogenous, error, rho)
    data_plm <- data[["plm"]]
    data_Tn <-  data[["Tn"]]
    
    # Perform estimation using method 1 (KSS)
    KSS.fit <- KSS(data_Tn$Y ~ data_Tn$X_1 + data_Tn$X_2 - 1)
    
    # Store results
    KSS_results <- list(
      b1 = KSS.fit$slope.para[1],
      b2 = KSS.fit$slope.para[2],
      used.dim = KSS.fit$used.dim,
      sd1 = sqrt(KSS.fit$beta.V[1, 1]),
      sd2 = sqrt(KSS.fit$beta.V[2, 2]),
      b1.bias = KSS.fit$slope.para[1] - beta[1],
      b2.bias = KSS.fit$slope.para[2] - beta[2],
      mse.effect = (1/T) * (1/n) * sum((KSS.fit$unob.fact.stru - data_Tn$V)^2),
      mse = (1/T) * (1/n) * sum(KSS.fit$residuals^2)
    )
    # Now, add mse.coeff to the KSS_results list
    KSS_results$mse.coeff <- ((KSS_results$b1 - beta[1])^2 + (KSS_results$b2 - beta[2])^2)  / 2
    
    # Perform estimation using method 2 (Eup)
    Eup.fit <- Eup(data_Tn$Y ~ data_Tn$X_1 + data_Tn$X_2 - 1)
    Eup.fit.summary <- summary(Eup.fit)
    
    # Store results
    Eup_results <- list(
      b1 = Eup.fit$slope.para[1],
      b2 = Eup.fit$slope.para[2],
      used.dim = Eup.fit$used.dim,
      sd1 = Eup.fit.summary$coefficients[1, 2],
      sd2 = Eup.fit.summary$coefficients[2, 2],
      b1.bias = Eup.fit$slope.para[1] - beta[1],
      b2.bias = Eup.fit$slope.para[2] - beta[2],
      mse.effect = (1/T) * (1/n) * sum((Eup.fit$unob.fact.stru - data_Tn$V)^2),
      mse = (1/T) * (1/n) * sum(Eup.fit$residuals^2)
    )
    # Now, add mse.coeff to the Eup_results list
    Eup_results$mse.coeff <- ((Eup_results$b1 - beta[1])^2 + (Eup_results$b2 - beta[2])^2)  / 2
    
    # Perform estimation using method 3 (Within)
    Within.fit <- plm(data_plm$Y ~ data_plm$X_1 + data_plm$X_2 - 1, data = data_plm, model = "within", effects = "individual")
    
    # Store results
    Within_results <- list(
      b1 = Within.fit$coefficients[1],
      b2 = Within.fit$coefficients[2],
      sd1 = sqrt(Within.fit$vcov[1, 1]),
      sd2 = sqrt(Within.fit$vcov[2, 2]),
      b1.bias = Within.fit$coefficients[1] - beta[1],
      b2.bias = Within.fit$coefficients[2] - beta[2],
      mse = (1/T) * (1/n) * sum(Within.fit$residuals^2)
    )
    # Now, add mse.coeff to the Within_results list
    Within_results$mse.coeff <- ((Within_results$b1 - beta[1])^2 + (Within_results$b2 - beta[2])^2)  / 2
    
    # Return a list of results for this simulation
    list(KSS = KSS_results, Eup = Eup_results, Within = Within_results)
  }
  
  stopCluster(cl)
  
  # Extract and structure the results
  structured_results <- extract_simulation_results(sim_results)
  
  # Return the structured results
  return(structured_results)
}



