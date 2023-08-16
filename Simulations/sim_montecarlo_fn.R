library(doParallel)
library(foreach)
library(lmtest)
source("sim_dgp_fn.R")


extract_simulation_results <- function(sim_results) {
  # Initialize results
  results <- list(
    KSS = list(),
    Eup = list(),
    Within = list(),
    Bai = list()
  )
  
  # Define metric names
  metrics_all <- c("b1", "b2", "used.dim", "sd1", "sd2", "mse", "b1.bias", "b2.bias", "mse.effect", "mse.coeff", "power1", "power2")
  metrics_within <- c("b1", "b2", "sd1", "sd2", "mse", "b1.bias", "b2.bias", "mse.coeff", "power1", "power2")
  
  # Number of simulations
  nsim <- length(sim_results)
  
  # Loop over each simulation result
  for (i in 1:nsim) {
    
    # Loop over each type (KSS, Eup, Within, Bai)
    for (type in c("KSS", "Eup", "Within", "Bai")) {
      
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


doMonteCarlo <- function(nsim, T, n, beta, DGP = "K3", endogenous = FALSE, error = "homo", rho = 0.5, error.type = 1) {
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
  required_packages <- c("mvtnorm", "matlib", "phtt", "plm", "doParallel", "foreach", "lmtest")
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
                           KSS.fit.summary <- summary(KSS.fit)
                           
                           # Store results
                           KSS_results <- list(
                             b1 = KSS.fit$slope.para[1],
                             b2 = KSS.fit$slope.para[2],
                             used.dim = KSS.fit$used.dim,
                             sd1 = sqrt(KSS.fit$beta.V[1, 1]),
                             sd2 = sqrt(KSS.fit$beta.V[2, 2]),
                             b1.bias = KSS.fit$slope.para[1] - beta[1],
                             b2.bias = KSS.fit$slope.para[2] - beta[2],
                             mse.effect = (1/n) * (1/T) * sum((KSS.fit$unob.fact.stru - data_Tn$V)^2),
                             mse = (1/n) * (1/T) * sum(KSS.fit$residuals^2),
                             power1 = KSS.fit.summary$coefficients[1,4] < 0.05,
                             power2 = KSS.fit.summary$coefficients[2,4] < 0.05
                           )
                           KSS_results$mse.coeff <- mean((KSS_results$b1 - beta[1])^2, (KSS_results$b2 - beta[2])^2)

                           
                           # Perform estimation using method 2 (Eup)
                           Eup.fit <- Eup(data_Tn$Y ~ data_Tn$X_1 + data_Tn$X_2 - 1)
                           Eup.fit.summary <- summary(Eup.fit, error.type = error.type)
                           Eup_results <- list(
                             b1 = Eup.fit.summary$coefficients[1,1],
                             b2 = Eup.fit.summary$coefficients[2,1],
                             used.dim = Eup.fit$used.dim,
                             sd1 = Eup.fit.summary$coefficients[1, 2],
                             sd2 = Eup.fit.summary$coefficients[2, 2],
                             b1.bias = Eup.fit$slope.para[1] - beta[1],
                             b2.bias = Eup.fit$slope.para[2] - beta[2],
                             mse.effect = (1/n) * (1/T) * sum((Eup.fit$unob.fact.stru - data_Tn$V)^2),
                             mse = (1/n) * (1/T) * sum(Eup.fit$residuals^2),
                             power1 = Eup.fit.summary$coefficients[1,4] < 0.05,
                             power2 = Eup.fit.summary$coefficients[2,4] < 0.05
                           )
                           
                           Eup_results$mse.coeff <- mean((Eup_results$b1 - beta[1])^2, (Eup_results$b2 - beta[2])^2)

                           
                           # Perform estimation using Bai's method (similar to Eup)
                           Bai.fit <- Eup(data_Tn$Y ~ data_Tn$X_1 + data_Tn$X_2 - 1, factor.dim = 8)
                           Bai.fit.summary <- summary(Bai.fit, error.type = error.type)
                           Bai_results <- list(
                             b1 = Bai.fit.summary$coefficients[1,1],
                             b2 = Bai.fit.summary$coefficients[2,1],
                             used.dim = Bai.fit$used.dim,
                             sd1 = Bai.fit.summary$coefficients[1, 2],
                             sd2 = Bai.fit.summary$coefficients[2, 2],
                             b1.bias = Bai.fit$slope.para[1] - beta[1],
                             b2.bias = Bai.fit$slope.para[2] - beta[2],
                             mse.effect = (1/n) * (1/T) * sum((Bai.fit$unob.fact.stru - data_Tn$V)^2),
                             mse = (1/n) * (1/T) * sum(Bai.fit$residuals^2),
                             p1 = Bai.fit.summary$coefficients[1,4],
                             p2 = Bai.fit.summary$coefficients[2,4]
                           )
                           
                           Bai_results$mse.coeff <- mean((Bai_results$b1 - beta[1])^2, (Bai_results$b2 - beta[2])^2)
                           Bai_results$power1 <- Bai_results$p1 < 0.05
                           Bai_results$power2 <- Bai_results$p2 < 0.05
                        
                           
                           # Perform estimation using method 3 (Within)
                           Within.fit <- plm(data_plm$Y ~ data_plm$X_1 + data_plm$X_2 - 1, data = data_plm, model = "within", effects = "individual")
                           
                           
                           if (error %in% c("hetero_time", "hetero_individual", "hetero_both")) {
                             robust_se <- coeftest(Within.fit, vcov = vcovHC(Within.fit, type="HC0"))
                             
                             Within_results <- list(
                               b1 = Within.fit$coefficients[1],
                               b2 = Within.fit$coefficients[2],
                               sd1 = robust_se[1, 2], 
                               sd2 = robust_se[2, 2], 
                               b1.bias = Within.fit$coefficients[1] - beta[1],
                               b2.bias = Within.fit$coefficients[2] - beta[2],
                               mse =  (1/n) * (1/T) * sum(Within.fit$residuals^2),
                               power1 = robust_se[1, 4] < 0.05, 
                               power2 = robust_se[2, 4] < 0.05  
                             )
                             
                             
                           } else if (error == "autocorr") {
                             robust_se <- coeftest(Within.fit, vcov = vcovHC(Within.fit, type = "sss"))
                             
                             Within_results <- list(
                               b1 = Within.fit$coefficients[1],
                               b2 = Within.fit$coefficients[2],
                               sd1 = robust_se[1, 2], 
                               sd2 = robust_se[2, 2], 
                               b1.bias = Within.fit$coefficients[1] - beta[1],
                               b2.bias = Within.fit$coefficients[2] - beta[2],
                               mse =  (1/n) * (1/T) * sum(Within.fit$residuals^2),
                               power1 = robust_se[1, 4] < 0.05, 
                               power2 = robust_se[2, 4] < 0.05  
                             )
                             
                           } else {
                             robust_se <- summary(Within.fit)
                             
                             Within_results <- list(
                               b1 = Within.fit$coefficients[1],
                               b2 = Within.fit$coefficients[2],
                               sd1 = robust_se$coefficients[1, 2], 
                               sd2 = robust_se$coefficients[2, 2], 
                               b1.bias = Within.fit$coefficients[1] - beta[1],
                               b2.bias = Within.fit$coefficients[2] - beta[2],
                               mse =  (1/n) * (1/T) * sum(Within.fit$residuals^2),
                               power1 = robust_se$coefficients[1, 4] < 0.05, 
                               power2 = robust_se$coefficients[2, 4] < 0.05  
                             )
                           }
                           
                           
                           Within_results$mse.coeff <- mean((Within_results$b1 - beta[1])^2, (Within_results$b2 - beta[2])^2)
                           
                           # Return a list of results for this simulation
                           list(KSS = KSS_results, Eup = Eup_results, Within = Within_results, Bai = Bai_results)
                         }
  
  stopCluster(cl)
  
  # Extract and structure the results
  structured_results <- extract_simulation_results(sim_results)
  
  # Return the structured results
  return(structured_results)
}




