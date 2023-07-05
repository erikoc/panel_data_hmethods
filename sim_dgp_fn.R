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
      nu_i <- rmvnorm(T, mean = rep(0, nrow(covariance_matrix)), sigma = diag(2))
      X_i[t, ] <- X_i[t - 1, ] %*% R + nu_i[t, ]
    }
    X_list[[i]] <- X_i
  }
  
  X <- do.call(rbind, X_list)
  
  return(X)
}

shift_mean <- function(X) {
  #' Shifts the mean of each group in the input matrix
  #'
  #' This function shifts the mean of each group in the input matrix to predefined target means.
  #'
  #' @param X The input matrix of dimension (n, T), where n is the number of samples and T is the number of features.
  #'
  #' @return The shifted matrix with the mean of each group adjusted to the target means.

  u_list <- list(
    c(5, 5),
    c(7.5, 7.5),
    c(10, 10)
  )
  
  group_size <- n*T / 3
  
  group_indices <- rep(1:3, each = group_size)
  
  X_grouped <- lapply(1:3, function(i) X[group_indices == i, , drop = FALSE])
  
  for (i in 1:3) {
    mean_diff <- u_list[[i]] - colMeans(X_grouped[[i]])
    if (any(mean_diff != 0)) {
      X_grouped[[i]] <- X_grouped[[i]] + mean_diff
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
  

  theta_0 <- rnorm(n, mean = 0, sd = 1) * 3
  theta_1 <- rnorm(n, mean = 0, sd = 1) * 3
  theta_2 <- rnorm(n, mean = 0, sd = 1) * 3
  
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



generate_V5 <- function(T,n){
  #' Generates the individual effects V using DGP B1
  #'
  #' This function generates the individual effects V for a panel data model using DGP K4.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #'
  #' @return The matrix V representing the individual effects.
  
  f_1 <- rnorm(T, mean = 0, sd = 1)
  f_2 <- rnorm(T, mean = 0, sd = 1)
  
  lamb_1 <- rnorm(n, mean = 0, sd = 1)
  lamb_2 <- rnorm(n, mean = 0, sd = 1)
  
  # Generating time-invariant error terms
  v_list <- vector("list", n)
  
  for (i in 1:n) {
    v_i <- matrix(NA, nrow = T, ncol = 1)
    for (t in 1:T) {
      v_i[t] <- lamb_1[i]*f_1[t] + lamb_2[i] * f_2[t]
    }
    v_list[[i]] <- v_i
  }
  
  V <- do.call(rbind, v_list)
  
  return(V)
}


generate_V6 <- function(T, n, mu = 0.5, sigma = 3) {
  #' Generates the individual effects V using DGP E2 incorporating Brownian motion with drift
  #'
  #' This function generates the individual effects V for a panel data model using Brownian motion with drift.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param mu The drift term.
  #' @param sigma The volatility term.
  #'
  #' @return The matrix V representing the individual effects.
  
  # Creating common time factor (Brownian motion with drift)
  r <- rep(NA, T)
  r[1] <- rnorm(1, mean = 0, sd = 1) * sigma + mu
  
  for (t in 2:T) {
    r[t] <- r[t - 1] + rnorm(1, mean = 0, sd = 1) * sigma + mu
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
  
  # Combining the individual effects into a matrix
  V <- do.call(rbind, v_list)
  
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
  V <- V
  W <- rho * V + sd(V) * sqrt(1 - rho^2) * eps
  
  # Update X$X_2 by adding W
  X[,2] <- X[,2] + W
  
  return(X)
}



generate_panel_errors <- function(T, n, X, error = "homo", rho = 0.5) {
  #' Generates the panel errors
  #'
  #' This function generates the panel errors for a panel data model.
  #'
  #' @param T The number of time periods.
  #' @param n The number of individuals.
  #' @param X The matrix representing the exogenous variables X.
  #' @param error The type of panel errors to generate. Supported values are "homo", "hetero", and "autocorr".
  #' @param rho The correlation coefficient for the autocorrelated errors.
  #'
  #' @return The vector of panel errors.

  if (error == "homo") {
    e <- matrix(rnorm(n * T), nrow = n, ncol = T)
  } else if (error == "hetero") {
    sigma <- sqrt(1 + X[,1]^2 + X[,2]^2)
    e <- matrix(rnorm(n * T, mean = 0, sd = sigma), nrow = T, ncol = n)
  } else if (error == "autocorr") {
    e <- matrix(rnorm(n * T), nrow = T, ncol = n)
    for (i in 1:n) {
      for (t in 2:T) {
        e[t, i] <- rho * e[t-1, i] + sqrt(1 - rho^2) * e[t, i]
      }
    }
  } else {
    stop("Invalid error type. Supported values are 'homo', 'hetero', and 'autocorr'.")
  }
  return(as.vector(e))
}


create_plm_data <- function(T, n, Y, X, V){
  
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
  X_1 <- X[, 1]
  X_2 <- X[, 2]
  X_1_model <- matrix(X_1, T, n)
  X_2_model <- matrix(X_2, T, n)
  V_model <- matrix(V, T, n)
  
  
  list <- list(Y = Y_model, X_1 = X_1_model, X_2 = X_2_model, V = V_model)
  
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
  #'
  #' @return A list containing Tn matrices for Y, X_1, X_2, and V.


  X <- generate_X(T, n)
  X <- shift_mean(X)
  
  # Create individual effects according to specified DGP
  if (DGP == "K1") {
    V <- generate_V1(T, n)
  } else if (DGP == "K2") {
    V <- generate_V2(T, n)
  } else if (DGP == "K3") {
    V <- generate_V3(T, n)
  } else if (DGP == "K4") {
    V <- generate_V4(T, n)
  } else if (DGP == "B1") {
  V <- generate_V5(T,n)
  } else if (DGP == "E1") {
  V <- generate_V6(T,n)
  }
  
  # Create endogenous regressors, when specified
  if (endogenous) {
    X <- endogenous_X(T, n, X, V)
  }
  
  # Generate error terms, as specified
  e <- generate_panel_errors(T, n, error)
  
  # The Model

  Y <- X %*% beta + V + e
  
  plm <- create_plm_data(T, n, Y, X, V)
  
  # Convert to T,n matrixes
  Tn <- create_Tn_mat(T, n, Y, X, V)
  
  return(list(plm = plm, Tn = Tn))
}





