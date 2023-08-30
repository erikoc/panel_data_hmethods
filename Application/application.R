setwd("C:/Users/spammer/Desktop/M.Sc. in Economics/So 2023/panel_data_hmethods/Application")
library(lmtest)
library(boot)
library(phtt)
library(ggplot2)
library(readstata13)
library(dplyr)
library(plm)
library(xtable)


set.seed(42)



####################### Data from AER ############################################

AER <- read.dta13("democracy-balanced-l4.dta")


######################## Dynamic Panel estimates #########################################################


# 0. Working the Data #
AER_panel <- pdata.frame(AER, index = c("id","year"))

AER_panel_Transitions <- AER_panel %>%
  arrange(country_name, year) %>%
  group_by(country_name) %>%
  mutate(lag_dem = dplyr::lag(dem, 1)) %>%
  ungroup()


AER_panel_Transitions$transition <- AER_panel_Transitions$dem - AER_panel_Transitions$lag_dem


transitions <- AER_panel_Transitions %>%
  filter(transition %in% c(-1, 1))



subset_df_1 <- AER_panel_Transitions %>% filter(transition == 1)


latex_table_1 <- xtable(subset_df_1)


print(latex_table_1, type = "latex")


subset_df_minus_1 <- AER_panel_Transitions %>% filter(transition == -1)


latex_table_minus_1 <- xtable(subset_df_minus_1)


print(latex_table_minus_1, type = "latex")




constant_values <- AER_panel %>%
  group_by(country_name) %>%
  summarize(overall_mean = mean(dem)) %>%
  ungroup() %>%
  filter(overall_mean %in% c(0,1))

always_1_countries <- constant_values %>%
  filter(overall_mean == 1) %>%
  pull(country_name)

always_0_countries <- constant_values %>%
  filter(overall_mean == 0) %>%
  pull(country_name)

######### 1. Fixed effects model##########
# 1.1 AER data - No differences



fe.fit1    <- plm(lgdp ~ dem + plm::lag(lgdp, 1) - 1, AER_panel, model = "within", effect = "twoways", index = c("id","year"))
fe.robust1 <- coeftest(fe.fit1 , vcov = vcovHC(fe.fit1, cluster = 'group'))

fe.fit2    <- plm(lgdp ~ dem + plm::lag(lgdp, 1:2) - 1, AER_panel, model = "within", effect = "twoways", index = c("id","year"))
fe.robust2 <- coeftest(fe.fit2 , vcov = vcovHC(fe.fit2, cluster = 'group'))

fe.fit3    <- plm(lgdp ~ dem + plm::lag(lgdp, 1:3) - 1, AER_panel, model = "within", effect = "twoways", index = c("id","year"))
fe.robust3 <- coeftest(fe.fit3 , vcov = vcovHC(fe.fit3, cluster = 'group'))

fe.fit4    <- plm(lgdp ~ dem + plm::lag(lgdp, 1:4) - 1, AER_panel, model = "within", effect = "twoways", index = c("id","year"))
fe.robust4 <- coeftest(fe.fit4 , vcov = vcovHC(fe.fit4, cluster = 'group'))

fe.robust1
fe.robust2
fe.robust3
fe.robust4


# 1.2 AER data - Differences

AER_diff <- AER %>%
  group_by(id) %>%
  mutate(dlgdp = c(NA, diff(lgdp))) %>%
  ungroup()

AER_diff <- AER_diff %>%
  filter(!is.na(dlgdp))


AER_panel <- pdata.frame(AER_diff, index = c("id","year"))


fe.fit1    <- plm(lgdp ~ dem + plm::lag(lgdp, 1) - 1, AER_diff, model = "within", effect = "twoways", index = c("id","year"))
fe.robust1 <- coeftest(fe.fit1 , vcov = vcovHC(fe.fit1, cluster = 'group'))

fe.fit2    <- plm(lgdp ~ dem + plm::lag(lgdp, 1:2) - 1, AER_diff, model = "within", effect = "twoways", index = c("id","year"))
fe.robust2 <- coeftest(fe.fit2 , vcov = vcovHC(fe.fit2, cluster = 'group'))

fe.fit3    <- plm(lgdp ~ dem + plm::lag(lgdp, 1:3) - 1, AER_diff, model = "within", effect = "twoways", index = c("id","year"))
fe.robust3 <- coeftest(fe.fit3 , vcov = vcovHC(fe.fit3, cluster = 'group'))

fe.fit4    <- plm(lgdp ~ dem + plm::lag(lgdp, 1:4) - 1, AER_diff, model = "within", effect = "twoways", index = c("id","year"))
fe.robust4 <- coeftest(fe.fit4 , vcov = vcovHC(fe.fit4, cluster = 'group'))

fe.robust1
fe.robust2
fe.robust3
fe.robust4

########## 2.Factor Models #################

#  2.1 AER data, Eup - PC1

n <- as.numeric(length(unique(AER$id)))
T <- as.numeric(length(unique(AER$year)))


lgdp.dy <- diff(matrix(AER$lgdp, T, n))
dem.dy <- matrix(AER$dem[-1], T-1, n)


Eup.fit.1 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) - 1,
                 additive.effects = "twoways", dim.criterion = "PC1")

Eup.fit.2 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) + dplyr::lag(lgdp.dy, 2)  - 1,
                 additive.effects = "twoways", dim.criterion = "PC1")

Eup.fit.3 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) + dplyr::lag(lgdp.dy, 2) + dplyr::lag(lgdp.dy, 3) - 1,
                 additive.effects = "twoways", dim.criterion = "PC1")

Eup.fit.4 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) + dplyr::lag(lgdp.dy, 2) + dplyr::lag(lgdp.dy, 3) + dplyr::lag(lgdp.dy, 4) - 1,
                 additive.effects = "twoways", dim.criterion = "PC1")


Eup.sum1 <-summary(Eup.fit.1, error.type = 5)
print(Eup.sum1)


Eup.sum2 <-summary(Eup.fit.2, error.type = 5)
print(Eup.sum2)


Eup.sum3 <-summary(Eup.fit.3, error.type = 5)
print(Eup.sum3)


Eup.sum4 <-summary(Eup.fit.4, error.type = 5)
print(Eup.sum4)
plot(Eup.sum4)

coef(Eup.fit.4)$Var.shares.of.loadings.param[1]
coef(Eup.fit.4)$Var.shares.of.loadings.param[2]
coef(Eup.fit.4)$Var.shares.of.loadings.param[3]
coef(Eup.fit.4)$Var.shares.of.loadings.param[4]



# 2.1.1 AER data, Eup - IPC1



Eup.fit.1 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) - 1,
                 additive.effects = "twoways", dim.criterion = "IPC1")

Eup.fit.2 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) + dplyr::lag(lgdp.dy, 2)  - 1,
                 additive.effects = "twoways", dim.criterion = "IPC1")

Eup.fit.3 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) + dplyr::lag(lgdp.dy, 2) + dplyr::lag(lgdp.dy, 3) - 1,
                 additive.effects = "twoways", dim.criterion = "IPC1")

Eup.fit.4 <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) + dplyr::lag(lgdp.dy, 2) + dplyr::lag(lgdp.dy, 3) + dplyr::lag(lgdp.dy, 4) - 1,
                 additive.effects = "twoways", dim.criterion = "IPC1")


Eup.sum1 <-summary(Eup.fit.1, error.type = 5)
print(Eup.sum1)


Eup.sum2 <-summary(Eup.fit.2, error.type = 5)
print(Eup.sum2)


Eup.sum3 <-summary(Eup.fit.3, error.type = 5)
print(Eup.sum3)


Eup.sum4 <-summary(Eup.fit.4, error.type = 5)
print(Eup.sum4)


# 2.1.2 Z-test

coef_estimate <- Eup.sum4$coefficients[1, 1]
standard_error <- Eup.sum4$coefficients[1, 2]
beta_upper <- fe.robust4[1,1]/4 # Define beta_upper here

z_value <- (coef_estimate - beta_upper) / standard_error
p_value <- pnorm(z_value)

alpha <- 0.05
z_alpha <- qnorm(alpha) # Get the z-value for the given alpha
lower_bound <- -Inf
upper_bound <- coef_estimate - z_alpha * standard_error

cat("z-value:", z_value, "\np-value:", p_value, "\n")
cat("95% One-sided Confidence Interval: [", lower_bound, ", ", upper_bound, "]", "\n")

# Decision based on p-value
if (p_value < alpha) {
  cat("Reject H0 in favor of H1: β is significantly less than", beta_upper, "\n")
} else {
  cat("Fail to reject H0: There's insufficient evidence to claim that β is less than", beta_upper, "\n")
}




# 2.1.3 Bootstrap regression

country_list <- split(AER, AER$id)


EupWrapper <- function(data_list, indices) {
  d <- do.call(rbind, data_list[indices])
  
  l.d.gdp.dy <- diff(matrix(d$lgdp, T, n))
  dem.dy <- matrix(d$dem[-1], T-1, n)
  
  Eup.fit <- Eup(lgdp.dy ~ dem.dy + dplyr::lag(lgdp.dy, 1) + dplyr::lag(lgdp.dy, 2) + dplyr::lag(lgdp.dy, 3) + dplyr::lag(lgdp.dy, 4) - 1,
                 additive.effects = "twoways", dim.criterion = "PC1")
  
  Eup.sum <-summary(Eup.fit, error.type = 5)
  
  return(Eup.sum$coefficients[,1]) 
}


boot_results_Eup <- boot(country_list, EupWrapper, R = 10000, stype = "i")

boot_cis <- lapply(1:5, function(i) boot.ci(boot.out = boot_results_Eup, index = i, type = c("norm", "perc", "bca")))

boot_results_Eup_summary <-summary(boot_results_Eup)
# Plots

bootstrap_results <- boot_results_Eup$t[,1]
median_value <- boot_results_Eup_summary[1,5]

ggplot(data.frame(BootstrapValues = bootstrap_results), aes(x = BootstrapValues)) +
  geom_histogram(aes(y = ..density..), fill = "blue", color = "black", binwidth = 0.0005, alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  geom_vline(aes(xintercept = median_value), color = "lightgray", linetype = "dashed", size = 1) + 
  annotate("text", x = median_value - 0.004, y = max(density(bootstrap_results)$y) + 0.1, label = paste("Median:", round(median_value, 4)), vjust = -0.5, color = "black") +
  labs(
    x = "Bootstrap Estimates",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )




## 2.2 AER - Data, KSS

l.gdp.dy <- matrix(AER$lgdp, T, n)
dem.dy.2 <- matrix(AER$dem, T, n)


KSS.fit4 <- KSS(l.gdp.dy ~ dem.dy.2 + dplyr::lag(l.gdp.dy, 1) + dplyr::lag(l.gdp.dy, 2) + dplyr::lag(l.gdp.dy, 3) + dplyr::lag(l.gdp.dy, 4) - 1,
               additive.effects = "twoways", consult.dim.crit = TRUE)

summary(KSS.fit4)
plot(summary(KSS.fit4))

coef(KSS.fit4)$Var.shares.of.loadings.param[1]
coef(KSS.fit4)$Var.shares.of.loadings.param[2]
coef(KSS.fit4)$Var.shares.of.loadings.param[3]
coef(KSS.fit4)$Var.shares.of.loadings.param[4]



country_list <- split(AER, AER$id)


KSSdyWrapper <- function(data_list, indices) {
  d <- do.call(rbind, data_list[indices])
  
  l.gdp.dy <- matrix(d$lgdp, T, n)
  dem.dy.2 <- matrix(d$dem, T, n)
  
  KSS.fit <- KSS(l.gdp.dy ~ dem.dy.2 + dplyr::lag(l.gdp.dy, 1) + dplyr::lag(l.gdp.dy, 2) + dplyr::lag(l.gdp.dy, 3) + dplyr::lag(l.gdp.dy, 4) - 1,
                 additive.effects = "twoways", factor.dim = 4)
  
  
  return(KSS.fit$slope.para) 
}

boot_results_KSSdy <- boot(country_list, KSSdyWrapper, R = 10000, stype = "i")

boot_cis <- lapply(1:5, function(i) boot.ci(boot.out = boot_results_KSSdy, index = i, type = c("norm", "perc", "bca")))

boot_results_KSSdy
boot_cis

# 2.2.1 Plot 
boot_results_KSS_summary <-summary(boot_results_KSSdy)

library(ggplot2)
library(scales)

bootstrap_results <- boot_results_KSSdy$t[,1]
median_value <- median(bootstrap_results)

ggplot(data.frame(BootstrapValues = bootstrap_results), aes(x = BootstrapValues)) +
  geom_histogram(aes(y = ..density..), fill = "blue", color = "black", binwidth = 0.0005, alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  geom_vline(aes(xintercept = median_value), color = "lightgray", linetype = "dashed", size = 1) + 
  annotate("text", x = median_value - 0.008, y = max(density(bootstrap_results)$y) + 0.1, label = paste("Median:", sprintf("%.4f", median_value)), vjust = -0.5, color = "black") +
  labs(
    x = "Bootstrap Estimates",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

