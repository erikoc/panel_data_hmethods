library(ggplot2)
library(gridExtra)
library(reshape2)

calculate_means <- function(results, path, estimator) {
  b1_mean <- mean(results[[path]][[estimator]]$b1.bias)
  b2_mean <- mean(results[[path]][[estimator]]$b2.bias)
  sd1_mean <- mean(results[[path]][[estimator]]$sd1)
  sd2_mean <- mean(results[[path]][[estimator]]$sd2)
  mse_mean <- mean(results[[path]][[estimator]]$mse.coeff)
  return(list(b1_mean = b1_mean, b2_mean = b2_mean, sd1_mean = sd1_mean, sd2_mean = sd2_mean, mse_mean = mse_mean))
}

calculate_means_for_paths_and_estimators <- function(all_results, paths = c("T12_n30_nsim1000_K1", "T30_n30_nsim1000_K1", "T12_n100_nsim1000_K1", "T30_n100_nsim1000_K1",
                                                                            "T12_n300_nsim500_K1","T30_n300_nsim500_K1" ), estimators = c("Eup", "Within", "KSS")) {
  #' Calculate means for specified paths and estimators
  #'
  #' This function loops through each specified path and estimator,
  #' and calculates means using a custom function `calculate_means`.
  #'
  #' @param all_results_homo A data structure containing the results.
  #' @param paths A vector of strings representing the paths to be processed.
  #' @param estimators A vector of strings representing the estimators to be processed.
  #'
  #' @return A list containing the means for each combination of path and estimator.
  
  # Create an empty list to store the results
  all_means_list <- list()
  
  # Loop through each path and estimator
  for (path in paths) {
    means_list <- list()
    for (estimator in estimators) {
      means_list[[estimator]] <- calculate_means(all_results, path, estimator)
    }
    all_means_list[[path]] <- means_list
  }
  
  return(all_means_list)
}


generate_latex_table <- function(all_means_list){
  categories <- c('KSS', 'Eup', 'Within')
  T_values <- c('T12', 'T30')
  n_values <- c('n30', 'n100', 'n300')
  
  cat('\\begin{tabular}{ccccccc}', "\n")
  cat('\\hline', "\n")
  cat('\\multicolumn{7}{c}{MSE, Bias, Variance for Coefficients} \\\\ \\hline', "\n")
  cat('& \\multicolumn{3}{c}{$T=12$} & \\multicolumn{3}{c}{$T=30$} \\\\ \\cline{2-4} \\cline{5-7}', "\n")
  cat('& KSS & Eup & Within & KSS & Eup & Within \\\\')
  
  for (n in n_values){
    nsim_value <- ifelse(n == "n300", "nsim50", "nsim100")
    n_display <- sub("n", "", n) # Remove 'n' prefix for display
    
    cat('\\multicolumn{7}{c}{$n =', n_display, '} \\\\')
    
    cat('MSE ')
    for (T in T_values){
      for (cat in categories){
        key <- paste(T, '_', n, '_', nsim_value, '_K1', sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$mse_mean, 4), nsmall=4))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ BIAS1 ')
    for (T in T_values){
      for (cat in categories){
        key <- paste(T, '_', n, '_', nsim_value, '_K1', sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$b1_mean, 4), nsmall=4))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ BIAS2 ')
    for (T in T_values){
      for (cat in categories){
        key <- paste(T, '_', n, '_', nsim_value, '_K1', sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$b2_mean, 4), nsmall=4))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ VAR1 ')
    for (T in T_values){
      for (cat in categories){
        key <- paste(T, '_', n, '_', nsim_value, '_K1', sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$sd1_mean, 4), nsmall=4))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ VAR2 ')
    for (T in T_values){
      for (cat in categories){
        key <- paste(T, '_', n, '_', nsim_value, '_K1', sep = '')
        if (cat %in% names(all_means_list[[key]])){
          cat(' &', format(round(all_means_list[[key]][[cat]]$sd2_mean, 4), nsmall=4))
        } else {
          cat(' & NULL')
        }
      }
    }
    
    cat('\\\\ \\hline', "\n")
  }
  
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
    KSS_b1_1 = results_type[[paste("T12_n30_nsim100_", DGP, sep = "")]]$KSS$b1,
    Eup_b1_1 = results_type[[paste("T12_n30_nsim100_", DGP, sep = "")]]$Eup$b1,
    Within_b1_1 = results_type[[paste("T12_n30_nsim100_", DGP, sep = "")]]$Within$b1,
    
    KSS_b1_2 = results_type[[paste("T30_n30_nsim100_", DGP, sep = "")]]$KSS$b1,
    Eup_b1_2 = results_type[[paste("T30_n30_nsim100_", DGP, sep = "")]]$Eup$b1,
    Within_b1_2 = results_type[[paste("T30_n30_nsim100_", DGP, sep = "")]]$Within$b1,
    
    KSS_b1_3 = results_type[[paste("T12_n100_nsim100_", DGP, sep = "")]]$KSS$b1,
    Eup_b1_3 = results_type[[paste("T12_n100_nsim100_", DGP, sep = "")]]$Eup$b1,
    Within_b1_3 = results_type[[paste("T12_n100_nsim100_", DGP, sep = "")]]$Within$b1,
    
    KSS_b1_4 = results_type[[paste("T30_n100_nsim100_", DGP, sep = "")]]$KSS$b1,
    Eup_b1_4 = results_type[[paste("T30_n100_nsim100_", DGP, sep = "")]]$Eup$b1,
    Within_b1_4 = results_type[[paste("T30_n100_nsim100_", DGP, sep = "")]]$Within$b1,
    
    KSS_b1_5 = results_type[[paste("T12_n300_nsim50_", DGP, sep = "")]]$KSS$b1,
    Eup_b1_5 = results_type[[paste("T12_n300_nsim50_", DGP, sep = "")]]$Eup$b1,
    Within_b1_5 = results_type[[paste("T12_n300_nsim50_", DGP, sep = "")]]$Within$b1,
    
    KSS_b1_6 = results_type[[paste("T30_n300_nsim50_", DGP, sep = "")]]$KSS$b1,
    Eup_b1_6 = results_type[[paste("T30_n300_nsim50_", DGP, sep = "")]]$Eup$b1,
    Within_b1_6 = results_type[[paste("T30_n300_nsim50_", DGP, sep = "")]]$Within$b1
  )
  
  # Extract data for each panel
  data_panel1 <- data.frame(
    KSS = data_combine$KSS_b1_1,
    Eup = data_combine$Eup_b1_1,
    Within = data_combine$Within_b1_1
  )
  
  data_panel2 <- data.frame(
    KSS = data_combine$KSS_b1_3,
    Eup = data_combine$Eup_b1_3,
    Within = data_combine$Within_b1_3
  )
  
  data_panel3 <- data.frame(
    KSS = data_combine$KSS_b1_2,
    Eup = data_combine$Eup_b1_2,
    Within = data_combine$Within_b1_2
  )
  
  data_panel4 <- data.frame(
    KSS = data_combine$KSS_b1_4,
    Eup = data_combine$Eup_b1_4,
    Within = data_combine$Within_b1_4
  )
  
  # New data for n = 300
  data_panel5 <- data.frame(
    KSS = data_combine$KSS_b1_5,
    Eup = data_combine$Eup_b1_5,
    Within = data_combine$Within_b1_5
  )
  
  data_panel6 <- data.frame(
    KSS = data_combine$KSS_b1_6,
    Eup = data_combine$Eup_b1_6,
    Within = data_combine$Within_b1_6
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
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
    labs(title="n30 T12")
  
  plot2 <- ggplot(data_panel2, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
    labs(title="n100 T12")
  
  plot3 <- ggplot(data_panel3, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
    labs(title="n30 T30")
  
  plot4 <- ggplot(data_panel4, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
    labs(title="n100 T30")
  
  # New plots for n = 300
  plot5 <- ggplot(data_panel5, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
    labs(title="n300 T12")
  
  plot6 <- ggplot(data_panel6, aes(x=variable, y=value, fill=variable)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, color="grey") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black") +
    labs(title="n300 T30")
  
  # Combine plots into 2x3 grid with a specified order
  grid.arrange(plot1, plot2, plot5, plot3, plot4, plot6, ncol=3)
}


library(xtable)

create_latex_table <- function(results_type, DGP) {
  # Calculate the mean values
  mean_values <- data.frame()
  
  combinations <- list(
    "T12_n30_nsim100_",
    "T30_n30_nsim100_",
    "T12_n100_nsim100_",
    "T30_n100_nsim100_",
    "T12_n300_nsim50_", 
    "T30_n300_nsim50_"
  )
  
  variables <- c("b1", "b2", "sd1", "sd2", "mse", "b1.bias", "b2.bias", "mse.effects")
  
  for (combination in combinations) {
    for (var in variables) {
      mean_values[paste(combination[1], var, sep="")] <- mean(results_type[[paste(combination[1], DGP, sep="")]]$KSS[[var]])
      mean_values[paste(combination[2], var, sep="")] <- mean(results_type[[paste(combination[2], DGP, sep="")]]$KSS[[var]])
    }
  }
  
  # Create LaTeX table
  latex_table <- xtable(mean_values)
  return(latex_table)
}

# Example usage
results_type <- all_results_homo  # Replace with your actual results_type
DGP <- "K3"  # Replace with your desired DGP

table <- create_latex_table(results_type, DGP)
print(table, type = "latex")

