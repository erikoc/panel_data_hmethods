source("sim_dgp_fn.R")

calculate_MSE_effects <- function(T, n, estimate, true_value, type) {
  #' Calculate Mean Squared Error (MSE) or Normalized Mean Squared Error (NMSE) for effects
  #'
  #' This function calculates either Mean Squared Error (MSE) or Normalized Mean Squared Error (NMSE)
  #' for the estimated effects compared to the true values.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param estimate Vector of estimated effect values.
  #' @param true_value Vector of true effect values.
  #' @param type The type of error measure to calculate. Supported values are "MSE" and "NMSE".
  #'
  #' @return The calculated error measure (MSE or NMSE).
  
  if (type == "MSE") {
    mse <- (1/T) * (1/n) * sum((estimate - true_value)^2)
  } else if (type == "NMSE") {
    mse <- sum((estimate - true_value)^2) / sum(true_value^2)
  }
  return(mse)
}


doMonteCarlo <- function(nsim, T, n, beta, DGP = "K3", endogenous = FALSE, error = "homo", type = "MSE") {
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
  #' @param type The type of error measure to calculate. Supported values are "MSE" and "NMSE". Default is "MSE".
  #'
  #' @return A list of estimation results for different methods and performance metrics.
  
  # Create empty lists to store results
  results <- list(
    KSS = list(b1 = numeric(), b2 = numeric(), used.dim = numeric(),
               sd1 = numeric(), sd2 = numeric(),
               mse = numeric(), b1.bias = numeric(), b2.bias = numeric(),
               mse.effect = numeric(), mse.coeff= numeric()),
    Eup = list(b1 = numeric(), b2 = numeric(), used.dim = numeric(),
               sd1 = numeric(), sd2 = numeric(),
               mse = numeric(), b1.bias = numeric(), b2.bias = numeric(),
               mse.effect = numeric(), mse.coeff= numeric()),
    Within = list(b1 = numeric(), b2 = numeric(),
                  sd1 = numeric(), sd2 = numeric(),
                  mse = numeric(), b1.bias = numeric(), b2.bias = numeric(), mse.coeff= numeric())
  )
  
  # Run the Monte Carlo simulation
  for (i in 1:nsim) {
    # Generate data
    data <- DataGeneratingFunction(T, n, beta, DGP, endogenous, error)
    data_plm <- data[["plm"]]
    
    data_Tn <-  data[["Tn"]]
    
    # Perform estimation using method 1 (KSS)
    KSS.fit <- KSS(data_Tn$Y ~ data_Tn$X_1 + data_Tn$X_2 -1)
    
    # Save results
    results$KSS$b1[i] <- KSS.fit$slope.para[1]
    results$KSS$b2[i] <- KSS.fit$slope.para[2]
    results$KSS$used.dim[i] <- KSS.fit$used.dim
    
    # Preformance metrics 
    results$KSS$sd1[i] <- sqrt(KSS.fit$beta.V[1, 1])
    results$KSS$sd2[i] <- sqrt(KSS.fit$beta.V[2, 2])
    results$KSS$b1.bias[i] <- KSS.fit$slope.para[1] - beta[1]
    results$KSS$b2.bias[i] <- KSS.fit$slope.para[2] - beta[2]
    results$KSS$mse.effect[i] <- calculate_MSE_effects(T, n, KSS.fit$unob.fact.stru,
                                                       data$V, type)
    results$KSS$mse[i] <- (1/T) * (1/n) *sum(KSS.fit$residuals^2)
    results$KSS$mse.coeff[i] <- (results$KSS$b1.bias[i]^2 + results$KSS$sd1[i]^2 +
                                results$KSS$b2.bias[i]^2 + results$KSS$sd2[i]^2)/2 
    
    # Perform estimation using method 2 (Eup)
    Eup.fit <- Eup(data_Tn$Y ~ data_Tn$X_1 + data_Tn$X_2 -1)
    Eup.fit.summary <- summary(Eup.fit)
    
    # Save results
    results$Eup$b1[i] <- Eup.fit$slope.para[1]
    results$Eup$b2[i] <- Eup.fit$slope.para[2]
    results$Eup$used.dim[i] <- Eup.fit$used.dim
    
    # Performance metrics
    results$Eup$sd1[i] <- Eup.fit.summary$coefficients[1, 2]
    results$Eup$sd2[i] <- Eup.fit.summary$coefficients[2, 2]
    results$Eup$b1.bias[i] <- Eup.fit$slope.para[1] - beta[1]
    results$Eup$b2.bias[i] <- Eup.fit$slope.para[2] - beta[2]
    results$Eup$mse.effect[i] <- calculate_MSE_effects(T, n, Eup.fit$unob.fact.stru,
                                                       data$V, type)
    results$Eup$mse[i] <- (1/T) * (1/n) * sum(Eup.fit$residuals^2)
    results$Eup$mse.coeff[i] <- (results$Eup$b1.bias[i]^2 + results$Eup$sd1[i]^2 +
                                   results$Eup$b2.bias[i]^2 + results$Eup$sd2[i]^2)/2 
    
    # Perform estimation using method 3 (Within)
    Within.fit <- plm(data_plm$Y ~ data_plm$X_1 + data_plm$X_2 -1, data = data_plm,  model = "within", effects = "twoways")
    
    # Save results
    results$Within$b1[i] <- Within.fit$coefficients[1]
    results$Within$b2[i] <- Within.fit$coefficients[2]
    
    # Preformance metrics 
    results$Within$sd1[i] <- sqrt(Within.fit$vcov[1, 1])
    results$Within$sd2[i] <- sqrt(Within.fit$vcov[2, 2])
    results$Within$b1.bias[i] <- Within.fit$coefficients[1] - beta[1]
    results$Within$b2.bias[i] <- Within.fit$coefficients[2] - beta[2]
    results$Within$mse[i] <- (1/T) * (1/n) *sum(Within.fit$residuals^2)
    results$Within$mse.coeff[i] <- (results$Within$b1.bias[i]^2 + results$Within$sd1[i]^2 +
                                   results$Within$b2.bias[i]^2 + results$Within$sd2[i]^2)/2 
    
  }
  
  
  
  # Return the estimation results
  return(results)
}
