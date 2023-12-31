library(mvtnorm)
library(matlib)
library(phtt)
library(plm)

generate_X <- function(T, n) {
  #' Generates the exogenous variables X
  #'
  #' This function generates the exogenous variables X for a panel data model.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #'
  #' @return The matrix X representing the exogenous variables.
  
  R <- rbind(c(0.4, 0.05), c(0.05, 0.4))
  covariance_matrix <- inv(diag(2) - R^2)
  
  
  X_list <- vector("list", n)
  
  for (i in 1:n) {
    X_i <- matrix(NA, nrow = T, ncol = 2)
    X_i[1, ] <- rmvnorm(1, mean = rep(0, nrow(covariance_matrix)), sigma = covariance_matrix)
    
    for (t in 2:T) {
      nu_i <- rmvnorm(1, mean = rep(0, nrow(covariance_matrix)), sigma = diag(2))
      X_i[t, ] <- X_i[t - 1, ] %*% R + nu_i
    }
    X_list[[i]] <- X_i
  }
  
  X <- do.call(rbind, X_list)
  
  return(X)
}


  shift_mean <- function(X, T, n) {
    #' Shifts the mean of each group in the input matrix
    #'
    #' This function shifts the mean of each group in the input matrix to predefined target means.
    #'
    #' @param X The input matrix of dimension (n*T, 2), where n is the number of individuals and T is the number of time periods.
    #'
    #' @return The shifted matrix with the mean of each group adjusted to the target means.
    
    u_list <- list(
      c(5, 5),
      c(7.5, 7.5),
      c(10, 10)
    )
    
    group_size <- n / 3
    
    group_indices <- rep(1:3, each = group_size * T)
    
    X_grouped <- lapply(1:3, function(i) X[group_indices == i, , drop = FALSE])
    
    for (i in 1:3) {
      mean_diff <- u_list[[i]] - colMeans(X_grouped[[i]])
      if (any(mean_diff != 0)) {
        X_grouped[[i]] <- X_grouped[[i]] + matrix(rep(mean_diff, each = nrow(X_grouped[[i]])), ncol = 2, byrow = TRUE)
      }
    }
    
    X_shifted <- do.call(rbind, X_grouped)
    
    return(X_shifted)
  }
  
  
generate_V1 <- function(T, n) {
  #' Generates the individual effects V using DGP K1
  #'
  #' This function generates the individual effects V for a panel data model using DGP K1.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #'
  #' @return The matrix V representing the individual effects.
  

  theta_0 <- rnorm(n, mean = 0, sd = 1) * 5
  theta_1 <- rnorm(n, mean = 0, sd = 1) * 5
  theta_2 <- rnorm(n, mean = 0, sd = 1) * 5
  
  v_list <- vector("list", n)
  
  for (i in 1:n) {
    v_i <- matrix(NA, nrow = T, ncol = 1)
    for (t in 1:T) {
      v_i[t] <- theta_0[i] + theta_1[i] * (t / T) + theta_2[i] * (t / T) ^ 2
    }
    v_list[[i]] <- v_i
  }
  
  V <- do.call(rbind, v_list)
  
  return(V)
}


generate_V2 <- function(T, n) {
  #' Generates the individual effects V using DGP K2
  #'
  #' This function generates the individual effects V for a panel data model using DGP K2.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #'
  #' @return The matrix V representing the individual effects.
  
  # Creating common time factor
  r <- rep(NA, T)
  r[1] <- rnorm(1, mean = 0, sd = 1) * 3
  
  for (t in 2:T) {
    r[t] <- r[t - 1] + rnorm(1, mean = 0, sd = 1) * 3
  }
  
  
  # Creating individual loading parameters
  phi <- rnorm(n, mean = 0, sd = 1) * 3
  
  # Generating time-invariant error terms
  v_list <- vector("list", n)
  
  for (i in 1:n) {
    v_i <- matrix(NA, nrow = T, ncol = 1)
    for (t in 1:T) {
      v_i[t] <- phi[i] * r[t]
    }
    v_list[[i]] <- v_i
  }
  
  V <- do.call(rbind, v_list)
  
  return(V)
}

generate_V3 <- function(T, n) {
  #' Generates the individual effects V using DGP K3
  #'
  #' This function generates the individual effects V for a panel data model using DGP K3.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #'
  #' @return The matrix V representing the individual effects.

  v_1 <- rnorm(n, mean = 0, sd = 1) * 3
  v_2 <- rnorm(n, mean = 0, sd = 1) * 3
  
  v_list <- vector("list", n)
  
  for (i in 1:n) {
    v_i <- matrix(NA, nrow = T, ncol = 1)
    for (t in 1:T) {
      v_i[t] <- v_1[i] * sin(pi * t / 4) + v_2[i] * cos(pi * t / 4)
    }
    v_list[[i]] <- v_i
  }
  
  V <- do.call(rbind, v_list)
}



generate_V4 <- function(T, n) {
  #' Generates the individual effects V using DGP K4
  #'
  #' This function generates the individual effects V for a panel data model using DGP K4.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #'
  #' @return The matrix V representing the individual effects.
  
  xi <- rnorm(n, mean = 0, sd = 1) * 3
  V <- rep(xi, each = T)
  return(V)
}


endogenous_X <- function(T, n, X, V, rho = 0.5) {
  #' Generates the endogenous variables X with correlated shocks
  #'
  #' This function generates the endogenous variables X with correlated shocks for a panel data model.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param X The matrix representing the exogenous variables X.
  #' @param V The matrix representing the individual effects V.
  #' @param rho The correlation coefficient for the shocks.
  #'
  #' @return The matrix X representing the endogenous variables with correlated shocks.
  
  
  # Generate random noise
  eps <- rnorm(n * T, mean = 0, sd = 1)
  
  # Calculate W
  V <- V*10
  W <- rho * V+ sd(V) * sqrt(1 - rho^2) * eps
  
  # Update X$X_2 by adding W
  X[,2] <- X[,2] + W
  
  return(X)
}



generate_panel_errors <- function(T, n, X, error = "homo") {
  #' Generates the panel errors
  #'
  #' This function generates the panel errors for a panel data model.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param X The matrix representing the exogenous variables X.
  #' @param error The type of panel errors to generate. Supported values are "homo", "hetero_time", "hetero_individual", "hetero_both" and "autocorr".
  #' @param rho The correlation coefficient for the autocorrelated errors.
  #'
  #' @return The vector of panel errors.
  
  if (error == "homo") {
    e <- rnorm(n * T)
  } else if (error == "hetero_time") {
    e <-  rep(NA, times = n*T) 
    for (i in 1:n) {
      for (t in 1:T) {
        sigma <- sqrt(1 + t)  # variance increasing over time
        e[(i-1)*T + t] <- rnorm(1, mean = 0, sd = sigma)
      }
    }
  } else if (error == "hetero_individual") {
    e <-  rep(NA, times = n*T)
    for (i in 1:n) {
      sigma_i <- sqrt(1 + i/10)  # variance increasing with individual
      for (t in 1:T) {
        e[(i-1)*T + t] <- rnorm(1, mean = 0, sd = sigma_i)
      }
    }
  } else if (error == "hetero_both") {
    e <-  rep(NA, times = n*T)
    for (i in 1:n) {
      for (t in 1:T) {
        sigma <- sqrt(1 + t + i/10)  # variance increasing with time and individual
        e[(i-1)*T + t] <- rnorm(1, mean = 0, sd = sigma)
      }
    }
  } else if (error == "autocorr_low") {
    rho <- runif(1, min = -0.3, max = 0.3)
    e <- rnorm(n * T)
    for (i in 1:n) {
      for (t in 2:T) {
        e[(i-1)*T + t] <- rho * e[(i-1)*T + t-1] + sqrt(1 - rho^2) * e[(i-1)*T + t]
      }
    }
  } else if (error == "autocorr_high") {
    rho <- runif(1, min = 0.6, max = 0.8)
    e <- rnorm(n * T)
    for (i in 1:n) {
      for (t in 2:T) {
        e[(i-1)*T + t] <- rho * e[(i-1)*T + t-1] + sqrt(1 - rho^2) * e[(i-1)*T + t]
      }
    }
  } else {
    stop("Invalid error type. Supported values are 'homo', 'hetero_time', 'hetero_individual', 'hetero_both', and 'autocorr'.")
  }
  
  return(e)
}


create_plm_data <- function(T, n, Y, X, V) {
  #' Create a Panel Data Frame
  #' 
  #' This function constructs a panel data frame with given input vectors or matrices.
  #' It assembles the data in long format and defines the individual and time indices.
  #' 
  #' @param T An integer, the number of time periods.
  #' @param n An integer, the number of cross-sectional units (individuals).
  #' @param Y A vector or matrix of dependent variable observations.
  #' @param X A matrix of independent variable observations with two columns.
  #' @param V A vector or matrix of additional variable observations.
  #' 
  #' @return A pdata.frame object with the assembled data, including individual and time indices.
  
  individual_index <- rep(1:n, each = T)
  time_index <- rep(1:T, times = n)
  
  plm <- data.frame(X_1 = X[,1], X_2 = X[,2], V = V, Y = Y, individual = individual_index, time = time_index)
  plm <- pdata.frame(plm, index = c("individual","time"))
  
  return(plm)
}


create_Tn_mat <- function(T, n, Y, X, V) {
  #' Creates Tn matrices from the data
  #'
  #' This function creates Tn matrices from the panel data.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param Y The vector of dependent variables.
  #' @param X The matrix representing the exogenous variables X.
  #' @param V The matrix representing the individual effects V.
  #'
  #' @return A list containing Tn matrices for Y, X_1, X_2, and V.
  Y_model <- matrix(Y, T, n)
  X1_model <- matrix(X[,1], T, n)
  X2_model <- matrix(X[,2], T, n)
  V_model <- matrix(V, T, n)
  
  
  list <- list(Y = Y_model, X_1 = X1_model, X_2 = X2_model, V = V_model)
  
  return(list)
}



DataGeneratingFunction <- function(T, n, beta, DGP, endogenous = FALSE, error = "homo") {
  #' Generates panel data using the specified parameters
  #'
  #' This function generates panel data using the specified parameters.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param beta The vector of true coefficients.
  #' @param DGP The data generating process for individual effects (supported values are "K1", "K2", "K3", "K4").
  #' @param endogenous Whether to generate endogenous regressors.
  #' @param error The type of panel errors to generate. Supported values are "homo", "hetero", and "autocorr".
  #' @param rho Level of error-autocorrelation for error type "autocorr". 
  #'
  #' @return A list containing Tn matrices for Y, X_1, X_2, and V.


  X <- generate_X(T, n)
  X <- shift_mean(X, T, n)
  
  # Create individual effects according to specified DGP
  if (DGP == "K1") {
    V <- generate_V1(T, n)
  } else if (DGP == "K2") {
    V <- generate_V2(T, n)
  } else if (DGP == "K3") {
    V <- generate_V3(T, n)
  } else if (DGP == "K4") {
    V <- generate_V4(T, n)
  }
  
  # Create endogenous regressors, when specified
  if (endogenous) {
    X <- endogenous_X(T, n, X, V)
  }
  
  # Generate error terms, as specified
  e <- generate_panel_errors(T, n, X, error)
  
  # The Model

  Y <- X %*% beta + V + e
  
  plm <- create_plm_data(T, n, Y, X, V)
  
  # Convert to T,n matrixes
  Tn <- create_Tn_mat(T, n, Y, X, V)
  
  # Determine true dimension
  
  # Extract information about d
  
  if (DGP == "K1") {
    d <- 3
  } else if (DGP == "K2") {
    d <- 1
  } else if (DGP == "K3") {
    d <- 2
  } else if (DGP == "K4") {
    d <- 1
  } else {
    stop("Unknown DGP value") # This will stop execution if DGP does not match any known value
  }
  
  return(list(plm = plm, Tn = Tn, d = d))
}





