library(ggplot2)
library(gridExtra)
library(reshape2)

calculate_means <- function(results, path, estimator) {
  b1_mean <- mean(results[[path]][[estimator]]$b1.bias)
  b2_mean <- mean(results[[path]][[estimator]]$b2.bias)
  sd1_mean <- mean(results[[path]][[estimator]]$sd1)
  sd2_mean <- mean(results[[path]][[estimator]]$sd2)
  mse_mean <- mean(results[[path]][[estimator]]$mse)
  
  if (estimator %in% c("KSS", "Eup")) {
    mse_effects_mean <- mean(results[[path]][[estimator]]$mse.effect)
    used.dim_mean <- mean(results[[path]][[estimator]]$used.dim)
    return(list(b1_mean = b1_mean, b2_mean = b2_mean, sd1_mean = sd1_mean, sd2_mean = sd2_mean, mse_mean = mse_mean, mse_effects_mean = mse_effects_mean, used.dim_mean = used.dim_mean))
  } else {
    return(list(b1_mean = b1_mean, b2_mean = b2_mean, sd1_mean = sd1_mean, sd2_mean = sd2_mean, mse_mean = mse_mean))
  }
}

calculate_means_for_paths_and_estimators <- function(all_results, DGP, 
                                                     paths = c("T12_n30_nsim1000_", "T30_n30_nsim1000_", 
                                                               "T12_n100_nsim1000_", "T30_n100_nsim1000_", 
                                                               "T12_n300_nsim500_", "T30_n300_nsim500_"), 
                                                     estimators = c("Eup", "KSS", "Within")) {
  #' Calculate means for specified paths and estimatWors
  #'
  #' This function loops through each specified path and estimator,
  #' and calculates means using a custom function `calculate_means`.
  #'
  #' @param all_results_homo A data structure containing the results.
  #' @param paths A vector of strings representing the paths to be processed.
  #' @param DGP A string containing the DGP name
  #' @param estimators A vector of strings representing the estimators to be processed.
  #'
  #' @return A list containing the means for each combination of path and estimator.
  # Update the paths with DGP
  paths <- paste(paths, DGP, sep = "")
  
  # Rest of the function remains the same
  all_means_list <- list()
  for (path in paths) {
    means_list <- list()
    for (estimator in estimators) {
      means_list[[estimator]] <- calculate_means(all_results, path, estimator)
    }
    all_means_list[[path]] <- means_list
  }
  return(all_means_list)
}



generate_latex_table <- function(all_means_list, DGP){
  categories <- c('KSS', 'Eup', 'Within')
  T_values <- c('T12', 'T30')
  n_values <- c('n30', 'n100', 'n300')
  
  cat('\\begin{tabular}{ccccccc}', "\n")
  cat('\\hline', "\n")
  cat('\\multicolumn{7}{c}{MSE, Bias, Variance for Coefficients} \\\\ \\hline', "\n")
  cat('& \\multicolumn{3}{c}{$T=12$} & \\multicolumn{3}{c}{$T=30$} \\\\ \\cline{2-4} \\cline{5-7}', "\n")
  cat('& KSS & Eup & Within & KSS & Eup & Within \\\\')
  
  for (n in n_values){
    nsim_value <- ifelse(n == "n300", "nsim500", "nsim1000")
    n_display <- sub("n", "", n) # Remove 'n' prefix for display
    
    cat('\\multicolumn{7}{l}{$n =', n_display, '} \\\\')
    
    cat('MSE ')
    for (T in T_values){
      for (cat in categories){
        key <- paste(T, '_', n, '_', nsim_value, '_', DGP,  sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$mse_mean, 5), nsmall=5))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ BIAS1 ')
    for (T in T_values){
      for (cat in categories){
        key <-  paste(T, '_', n, '_', nsim_value, '_', DGP,  sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$b1_mean, 5), nsmall=4))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ BIAS2 ')
    for (T in T_values){
      for (cat in categories){
        key <-  paste(T, '_', n, '_', nsim_value, '_', DGP,  sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$b2_mean, 5), nsmall=5))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ VAR1 ')
    for (T in T_values){
      for (cat in categories){
        key <-  paste(T, '_', n, '_', nsim_value, '_', DGP,  sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$sd1_mean, 5), nsmall=5))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ VAR2 ')
    for (T in T_values){
      for (cat in categories){
        key <-  paste(T, '_', n, '_', nsim_value, '_', DGP,  sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$sd2_mean, 5), nsmall=5))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ \\hline', "\n")
  }
  
  cat('\\end{tabular}', "\n")
}

generate_latex_mse_effects_table <- function(results_type, DGP) {
  n_values <- c(30, 100, 300)
  T_values <- c(12, 30)
  
  cat('\\begin{tabular}{lcccccc}', "\n")
  cat('\\hline \\multicolumn{6}{c}{MSE of Effects} \\\\ \\hline', "\n")
  cat('$n$ & $T$ & KSS & Eup & $d_{KSS}$ & $d_{Eup}$ \\\\\n')
  cat('\\hline\n')
  
  for (n in n_values) {
    nsim <- if (n == 300) { "nsim500" } else { "nsim1000" }
    for (T in T_values) {
      key <- paste("T", T, "_n", n, "_", nsim, "_", DGP, sep = "")
      
      if (T == 12) {
        cat(n, "&", T)
      } else {
        cat("&", T)
      }
      
      cat(" & ", format(round(results_type[[key]]$KSS$mse_effects_mean, 4), nsmall=4),
          " & ", format(round(results_type[[key]]$Eup$mse_effects_mean, 4), nsmall=4),
          " & ", format(round(results_type[[key]]$KSS$used.dim_mean, 4), nsmall=4),
          " & ", format(round(results_type[[key]]$Eup$used.dim_mean, 4), nsmall=4),
          " \\\\\n")
    }
  }
  
  cat('\\hline\n')
  cat('\\end{tabular}', "\n")
}


create_plots <- function(results_type, DGP) {
  #' Create violin plots in a 2x3 grid for given results type and DGP.
  #' 
  #'
  #' @param results_type List with the Monte Carlo results.
  #' @param DGP Type of DGP.
  #'
  #' @return The matrix X representing the exogenous variables.

  # Construct the data_combine dataframe dynamically based on DGP
  data_combine <- data.frame(
    KSS_b_1 = results_type[[paste("T12_n30_nsim1000_", DGP, sep = "")]]$KSS$mse,
    Eup_b_1 = results_type[[paste("T12_n30_nsim1000_", DGP, sep = "")]]$Eup$mse,
    #Within_b_1 = results_type[[paste("T12_n30_nsim1000_", DGP, sep = "")]]$Within$mse.coeff,
    
    KSS_b_2 = results_type[[paste("T30_n30_nsim1000_", DGP, sep = "")]]$KSS$mse,
    Eup_b_2 = results_type[[paste("T30_n30_nsim1000_", DGP, sep = "")]]$Eup$mse,
    #Within_b_2 = results_type[[paste("T30_n30_nsim1000_", DGP, sep = "")]]$Within$mse.coeff,
    
    KSS_b_3 = results_type[[paste("T12_n100_nsim1000_", DGP, sep = "")]]$KSS$mse,
    Eup_b_3 = results_type[[paste("T12_n100_nsim1000_", DGP, sep = "")]]$Eup$mse,
    #Within_b_3 = results_type[[paste("T12_n100_nsim1000_", DGP, sep = "")]]$Within$mse.coeff,
    
    KSS_b_4 = results_type[[paste("T30_n100_nsim1000_", DGP, sep = "")]]$KSS$mse,
    Eup_b_4 = results_type[[paste("T30_n100_nsim1000_", DGP, sep = "")]]$Eup$mse,
   # Within_b_4 = results_type[[paste("T30_n100_nsim1000_", DGP, sep = "")]]$Within$mse.coeff,
    
    KSS_b_5 = results_type[[paste("T12_n300_nsim500_", DGP, sep = "")]]$KSS$mse,
    Eup_b_5 = results_type[[paste("T12_n300_nsim500_", DGP, sep = "")]]$Eup$mse,
    #Within_b_5 = results_type[[paste("T12_n300_nsim500_", DGP, sep = "")]]$Within$mse.coeff,
    
    KSS_b_6 = results_type[[paste("T30_n300_nsim500_", DGP, sep = "")]]$KSS$mse,
    Eup_b_6 = results_type[[paste("T30_n300_nsim500_", DGP, sep = "")]]$Eup$mse
    #Within_b_6 = results_type[[paste("T30_n300_nsim500_", DGP, sep = "")]]$Within$mse.coeff
  )
  
  # Extract data for each panel
  data_panel1 <- data.frame(
    KSS = data_combine$KSS_b_1,
    Eup = data_combine$Eup_b_1
    #Within = data_combine$Within_b_1
  )
  
  data_panel2 <- data.frame(
    KSS = data_combine$KSS_b_3,
    Eup = data_combine$Eup_b_3
    #Within = data_combine$Within_b_3
  )
  
  data_panel3 <- data.frame(
    KSS = data_combine$KSS_b_2,
    Eup = data_combine$Eup_b_2
    #Within = data_combine$Within_b_2
  )
  
  data_panel4 <- data.frame(
    KSS = data_combine$KSS_b_4,
    Eup = data_combine$Eup_b_4
    #Within = data_combine$Within_b_4
  )
  
  # New data for n = 300
  data_panel5 <- data.frame(
    KSS = data_combine$KSS_b_5,
    Eup = data_combine$Eup_b_5
    #Within = data_combine$Within_b_5
  )
  
  data_panel6 <- data.frame(
    KSS = data_combine$KSS_b_6,
    Eup = data_combine$Eup_b_6
    #Within = data_combine$Within_b_6
  )
  
  
  print(str(data_panel1))
  
  
  # Melt the data
  data_panel1 <- melt(data_panel1)
  data_panel2 <- melt(data_panel2)
  data_panel3 <- melt(data_panel3)
  data_panel4 <- melt(data_panel4)
  data_panel5 <- melt(data_panel5)
  data_panel6 <- melt(data_panel6)
  
  # Create plots
  plot1 <- ggplot(data_panel1, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    labs(title="n30 T12")
  
  plot2 <- ggplot(data_panel2, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    labs(title="n100 T12")
  
  plot3 <- ggplot(data_panel3, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    labs(title="n30 T30")
  
  plot4 <- ggplot(data_panel4, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    labs(title="n100 T30")
  
  # New plots for n = 300
  plot5 <- ggplot(data_panel5, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    labs(title="n300 T12")
  
  plot6 <- ggplot(data_panel6, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    labs(title="n300 T30")
  
  # Combine plots into 2x3 grid with a specified order
  grid.arrange(plot1, plot2, plot5, plot3, plot4, plot6, ncol=3)
}


