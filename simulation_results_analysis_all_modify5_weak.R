
rm(list = ls())
gc()

# load the simulated data
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_0weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_10weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_25weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_50weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_100weak.rda")


# load("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_0weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_10weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_25weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_50weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_100weak.rda")

# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")

##########################################

dim(sim_data_list[[1]]$u)
library(corrplot)

x <- sim_data_list[[1]]$x
# Correlation matrix
cor_x <- cor(x, use = "pairwise.complete.obs")
# Plot
corrplot(
  cor_x,
  method = "color",
  type = "lower",
  order = "hclust",
  tl.cex = 0.3,
  tl.col = "black",
  cl.cex = 0.7,
  diag = FALSE
)
plot(cor_x)




z <- sim_data_list[[1]]$z
# Correlation matrix
cor_z <- cor(z, use = "pairwise.complete.obs")
# Plot
corrplot(
  cor_z,
  method = "color",
  type = "lower",
  order = "hclust",
  tl.cex = 0.3,
  tl.col = "black",
  cl.cex = 0.7,
  diag = FALSE
)
plot(cor_z)


u <- sim_data_list[[1]]$u
# Correlation matrix
cor_u <- cor(u, use = "pairwise.complete.obs")
# Plot
corrplot(
  cor_u,
  method = "color",
  type = "lower",
  order = "hclust",
  tl.cex = 0.3,
  tl.col = "black",
  cl.cex = 0.7,
  diag = FALSE
)
plot(cor_u)



d <- sim_data_list[[1]]$d
# Correlation matrix
cor_d <- cor(d, use = "pairwise.complete.obs")
# Plot
corrplot(
  cor_d,
  method = "color",
  type = "lower",
  order = "hclust",
  tl.cex = 0.3,
  tl.col = "black",
  cl.cex = 0.7,
  diag = FALSE
)
plot(cor_d)



hist(x,
      breaks = 60,
      main = "Distribution of x",
      xlab = "x",
      ylab = "Frequency" ,
      col = "white")
abline(v = 0, lty = 1, col = "red", lwd = 3)


hist(z,
      breaks = 60,
      main = "Distribution of z",
      xlab = "z",
      ylab = "Frequency" ,
      col = "white")
abline(v = 0, lty = 1, col = "red", lwd = 3)


# Draw a separate histogram for each column of z
par(mfrow=c(2,3))
for (j in 1:ncol(z)) {
  hist(z[, j],
       breaks = 60,
       main = paste("Distribution of z[,", j, "]"),
       xlab = paste0("z[,", j, "]"),
       ylab = "Frequency",
       col = "white")

  abline(v = 0, lty = 1, col = "red", lwd = 3)
}




# Draw a separate histogram for each column of x
par(mfrow=c(2,3))
for (j in 1:ncol(x)) {
  hist(x[, j],
       breaks = 60,
       main = paste("Distribution of x[,", j, "]"),
       xlab = paste0("x[,", j, "]"),
       ylab = "Frequency",
       col = "white")

  abline(v = 0, lty = 1, col = "red", lwd = 3)
}


hist( u,
      breaks = 60,
      main = "Distribution of u",
      xlab = "u",
      ylab = "Frequency" ,
      col = "white")
abline(v = 0, lty = 1, col = "red", lwd = 3)


# Draw a separate histogram for each column of d
# m=13, j=2, implies jm=17+13=30
par(mfrow=c(2,3))
for (jm in 25:30) {
  hist(u[, jm],
       breaks = 60,
       main = paste("Distribution of u[,", jm, "]"),
       xlab = paste0("u[,", jm, "]"),
       ylab = "Frequency",
       col = "white")

  abline(v = 0, lty = 1, col = "red", lwd = 3)
}


hist( d,
      breaks = 60,
      main = "Distribution of d",
      xlab = "d",
      ylab = "Frequency" ,
      col = "white")
abline(v = 0, lty = 1, col = "red", lwd = 3)


#####################################################

#####################
# HHS model

# load("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_0weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_10weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_25weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_50weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_100weak.rda")


# load("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_0weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_10weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_25weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_50weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_100weak.rda")


# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")


nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]
J <- dim(sim_data_list[[1]]$x)[2]
M <- dim(sim_data_list[[1]]$z)[2]

alpha_cover_per_HHS <- c()
beta_cover_per_HHS <- c()
gamma_cover_per_HHS <- c()
varphi_cover_per_HHS <- c()
delta_cover_per_HHS <- c()

delta_RMSE_HHS <- c()
RMSE_all_main_effects_HHS <- c()

beta_cover_nonzero_per_HHS <- c()
gamma_cover_nonzero_per_HHS <- c()
varphi_cover_nonzero_per_HHS <- c()
delta_cover_nonzero_per_HHS <- c()

alpha_est_mat_HHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_HHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_HHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_HHS <- matrix(NA, nrow = nsample, ncol = K)

beta_est_list_HHS <- list()
beta_true_list_HHS <- list()
beta_lower_list_HHS <- list()
beta_upper_list_HHS <- list()

gamma_est_list_HHS <- list()
gamma_true_list_HHS <- list()
gamma_lower_list_HHS <- list()
gamma_upper_list_HHS <- list()

varphi_est_list_HHS <- list()
varphi_true_list_HHS <- list()
varphi_lower_list_HHS <- list()
varphi_upper_list_HHS <- list()

delta_est_list_HHS <- list()
delta_true_list_HHS <- list()
delta_lower_list_HHS <- list()
delta_upper_list_HHS <- list()

main_sensitivity_HHS <- c()
main_specificity_HHS <- c()
main_precision_HHS <- c()

interaction_sensitivity_HHS <- c()
interaction_specificity_HHS <- c()
interaction_precision_HHS <- c()

main_TP_HHS <- c()
main_TN_HHS <- c()
main_FP_HHS <- c()
main_FN_HHS <- c()

interaction_TP_HHS <- c()
interaction_TN_HHS <- c()
interaction_FP_HHS <- c()
interaction_FN_HHS <- c()


# ============================================================
# Additional diagnostics
# ============================================================

beta_bias_HHS <- c()
gamma_bias_HHS <- c()
varphi_bias_HHS <- c()
delta_bias_HHS <- c()

beta_bias_nonzero_HHS <- c()
gamma_bias_nonzero_HHS <- c()
varphi_bias_nonzero_HHS <- c()
delta_bias_nonzero_HHS <- c()

beta_bias_zero_HHS <- c()
gamma_bias_zero_HHS <- c()
varphi_bias_zero_HHS <- c()
delta_bias_zero_HHS <- c()

beta_RMSE_HHS <- c()
gamma_RMSE_HHS <- c()
varphi_RMSE_HHS <- c()
delta_RMSE_HHS <- c()

beta_RMSE_nonzero_HHS <- c()
gamma_RMSE_nonzero_HHS <- c()
varphi_RMSE_nonzero_HHS <- c()
delta_RMSE_nonzero_HHS <- c()

beta_RMSE_zero_HHS <- c()
gamma_RMSE_zero_HHS <- c()
varphi_RMSE_zero_HHS <- c()
delta_RMSE_zero_HHS <- c()

beta_TP_HHS <- c()
gamma_TP_HHS <- c()
varphi_TP_HHS <- c()

beta_TN_HHS <- c()
gamma_TN_HHS <- c()
varphi_TN_HHS <- c()

beta_FP_HHS <- c()
gamma_FP_HHS <- c()
varphi_FP_HHS <- c()

beta_FN_HHS <- c()
gamma_FN_HHS <- c()
varphi_FN_HHS <- c()

beta_sensitivity_HHS <- c()
gamma_sensitivity_HHS <- c()
varphi_sensitivity_HHS <- c()

beta_specificity_HHS <- c()
gamma_specificity_HHS <- c()
varphi_specificity_HHS <- c()

beta_precision_HHS <- c()
gamma_precision_HHS <- c()
varphi_precision_HHS <- c()

delta_RMSE_indices_mat_HHS <- matrix(NA, nrow = nsample, ncol = J*M)
delta_cover_indices_mat_HHS <- matrix(NA, nrow = nsample, ncol = J*M)
delta_var_indices_mat_HHS <- matrix(NA, nrow = nsample, ncol = J*M)

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  d <- sim_data_list[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_HHS[sample, ] <- alpha_true
  alpha_est_mat_HHS[sample, ] <- apply(theta_update_HHS_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_HHS_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_HHS[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_HHS_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_HHS[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_HHS[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_HHS[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_HHS_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_HHS_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_HHS_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_HHS[[sample]] <- beta_est
  beta_lower_list_HHS[[sample]] <- beta_lower
  beta_upper_list_HHS[[sample]] <- beta_upper
  beta_cover_per_HHS[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_HHS[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # Additional diagnostics
  nonzero_beta <- beta_true != 0
  zero_beta <- beta_true == 0

  beta_bias <- beta_est - beta_true

  beta_bias_HHS[[sample]] <- beta_bias
  beta_bias_nonzero_HHS[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_zero_HHS[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_HHS[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_nonzero_HHS[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_zero_HHS[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_HHS[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_HHS_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_HHS_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_HHS_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_HHS[[sample]] <- gamma_est
  gamma_lower_list_HHS[[sample]] <- gamma_lower
  gamma_upper_list_HHS[[sample]] <- gamma_upper
  gamma_cover_per_HHS[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_HHS[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  # Additional diagnostics
  nonzero_gamma <- gamma_true != 0
  zero_gamma <- gamma_true == 0

  gamma_bias <- gamma_est - gamma_true

  gamma_bias_HHS[[sample]] <- gamma_bias
  gamma_bias_nonzero_HHS[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_zero_HHS[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_HHS[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_nonzero_HHS[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_zero_HHS[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

  # varphi
  varphi_true <- sim_data_list[[sample]]$varphi_true
  varphi_true_list_HHS[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_HHS_list[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_HHS_list[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_HHS_list[[sample]][, o, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    varphi_cover[(o-1-J-M-J*M), ] <- (varphi_true[(o-1-J-M-J*M), ] >= varphi_lower[(o-1-J-M-J*M), ]) &
      (varphi_true[(o-1-J-M-J*M), ] <= varphi_upper[(o-1-J-M-J*M), ])
    varphi_MSE[o-1-J-M-J*M, ] <- (varphi_est[(o-1-J-M-J*M), ] - varphi_true[o-1-J-M-J*M, ])^2
  }
  varphi_est_list_HHS[[sample]] <- varphi_est
  varphi_lower_list_HHS[[sample]] <- varphi_lower
  varphi_upper_list_HHS[[sample]] <- varphi_upper
  varphi_cover_per_HHS[sample] <- (sum(varphi_cover)/(dim(varphi_true)[1]*dim(varphi_true)[2]))*100
  nonzero_idx_varphi <- varphi_true != 0
  varphi_cover_nonzero_per_HHS[sample] <- (sum(varphi_cover[nonzero_idx_varphi]) / sum(nonzero_idx_varphi)) * 100

  # Additional diagnostics
  nonzero_varphi <- varphi_true != 0
  zero_varphi <- varphi_true == 0

  varphi_bias <- varphi_est - varphi_true

  varphi_bias_HHS[[sample]] <- varphi_bias
  varphi_bias_nonzero_HHS[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_zero_HHS[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_HHS[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_nonzero_HHS[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_zero_HHS[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


  main_MSE <- c(beta_MSE, gamma_MSE, varphi_MSE)
  RMSE_all_main_effects_HHS[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_HHS[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HHS_list[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_HHS_list[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_HHS_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_HHS_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_HHS[[sample]] <- delta_est
  delta_lower_list_HHS[[sample]] <- delta_lower
  delta_upper_list_HHS[[sample]] <- delta_upper
  delta_cover_per_HHS[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_HHS[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_HHS[sample] <- sqrt(mean(delta_MSE))

  delta_cover_indices_mat_HHS[sample, ] <- rowMeans(delta_cover)*100
  delta_RMSE_indices_mat_HHS[sample, ] <- sqrt(rowMeans(delta_MSE))
  delta_var_indices_mat_HHS[sample, ] <- rowMeans(delta_var)

  # Additional diagnostics
  nonzero_delta <- delta_true != 0
  zero_delta <- delta_true == 0

  delta_bias <- delta_est - delta_true

  delta_bias_HHS[[sample]] <- delta_bias
  delta_bias_nonzero_HHS[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_zero_HHS[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_HHS[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_nonzero_HHS[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_zero_HHS[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_varphi <- as.vector(varphi_true != 0)
  selected_varphi <- as.vector((varphi_lower > 0) | (varphi_upper < 0))

  truth_main <- c(truth_beta, truth_gamma, truth_varphi)
  selected_main <- c(selected_beta, selected_gamma, selected_varphi)

  main_TP_HHS[sample] <- sum(truth_main == 1 & selected_main == 1)
  main_FN_HHS[sample] <- sum(truth_main == 1 & selected_main == 0)
  main_TN_HHS[sample] <- sum(truth_main == 0 & selected_main == 0)
  main_FP_HHS[sample] <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_HHS[sample] <- (main_TP_HHS[sample] /
                                        (main_TP_HHS[sample] + main_FN_HHS[sample]))
  main_specificity_HHS[sample] <- (main_TN_HHS[sample] /
                                        (main_TN_HHS[sample] + main_FP_HHS[sample]))
  main_precision_HHS[sample] <- (main_TP_HHS[sample] /
                                      (main_TP_HHS[sample] + main_FP_HHS[sample]))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP_HHS[sample] <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN_HHS[sample] <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN_HHS[sample] <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP_HHS[sample] <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_HHS[sample] <- (interaction_TP_HHS[sample] /
                                               (interaction_TP_HHS[sample] + interaction_FN_HHS[sample]))
  interaction_specificity_HHS[sample] <- (interaction_TN_HHS[sample] /
                                               (interaction_TN_HHS[sample] + interaction_FP_HHS[sample]))
  interaction_precision_HHS[sample] <- (interaction_TP_HHS[sample] /
                                             (interaction_TP_HHS[sample] + interaction_FP_HHS[sample]))

  # Additional diaagnostics
  beta_TP_HHS[sample] <- sum(truth_beta == 1 & selected_beta == 1)
  beta_FN_HHS[sample] <- sum(truth_beta == 1 & selected_beta == 0)
  beta_TN_HHS[sample] <- sum(truth_beta == 0 & selected_beta == 0)
  beta_FP_HHS[sample] <- sum(truth_beta == 0 & selected_beta == 1)

  gamma_TP_HHS[sample] <- sum(truth_gamma == 1 & selected_gamma == 1)
  gamma_FN_HHS[sample] <- sum(truth_gamma == 1 & selected_gamma == 0)
  gamma_TN_HHS[sample] <- sum(truth_gamma == 0 & selected_gamma == 0)
  gamma_FP_HHS[sample] <- sum(truth_gamma == 0 & selected_gamma == 1)

  varphi_TP_HHS[sample] <- sum(truth_varphi == 1 & selected_varphi == 1)
  varphi_FN_HHS[sample] <- sum(truth_varphi == 1 & selected_varphi == 0)
  varphi_TN_HHS[sample] <- sum(truth_varphi == 0 & selected_varphi == 0)
  varphi_FP_HHS[sample] <- sum(truth_varphi == 0 & selected_varphi == 1)

  beta_sensitivity_HHS[sample] <- (beta_TP_HHS[sample] /
                                        (beta_TP_HHS[sample] + beta_FN_HHS[sample]))
  beta_specificity_HHS[sample] <- (beta_TN_HHS[sample] /
                                        (beta_TN_HHS[sample] + beta_FP_HHS[sample]))
  beta_precision_HHS[sample] <- (beta_TP_HHS[sample] /
                                      (beta_TP_HHS[sample] + beta_FP_HHS[sample]))

  gamma_sensitivity_HHS[sample] <- (gamma_TP_HHS[sample] /
                                         (gamma_TP_HHS[sample] + gamma_FN_HHS[sample]))
  gamma_specificity_HHS[sample] <- (gamma_TN_HHS[sample] /
                                         (gamma_TN_HHS[sample] + gamma_FP_HHS[sample]))
  gamma_precision_HHS[sample] <- (gamma_TP_HHS[sample] /
                                       (gamma_TP_HHS[sample] + gamma_FP_HHS[sample]))

  varphi_sensitivity_HHS[sample] <- (varphi_TP_HHS[sample] /
                                          (varphi_TP_HHS[sample] + varphi_FN_HHS[sample]))
  varphi_specificity_HHS[sample] <- (varphi_TN_HHS[sample] /
                                          (varphi_TN_HHS[sample] + varphi_FP_HHS[sample]))
  varphi_precision_HHS[sample] <- (varphi_TP_HHS[sample] /
                                        (varphi_TP_HHS[sample] + varphi_FP_HHS[sample]))
}




round(mean(RMSE_all_main_effects_HHS)*100, 2)
round(mean(delta_RMSE_HHS)*100, 2)

round((mean(beta_cover_per_HHS) + mean(gamma_cover_per_HHS) +
         mean(varphi_cover_per_HHS))/3, 2)
round(mean(delta_cover_per_HHS), 2)

round((mean(beta_cover_nonzero_per_HHS) + mean(gamma_cover_nonzero_per_HHS) +
         mean(varphi_cover_nonzero_per_HHS))/3, 2)
round(mean(delta_cover_nonzero_per_HHS), 2)

round(mean(main_sensitivity_HHS) * 100, 2)
round(mean(interaction_sensitivity_HHS) * 100, 2)

round(mean(main_specificity_HHS) * 100, 2)
round(mean(interaction_specificity_HHS) * 100, 2)

round(mean(main_precision_HHS) * 100, 2)
round(mean(interaction_precision_HHS) * 100, 2)


###############################

round(mean(main_TP_HHS), 2)
round(mean(interaction_TP_HHS), 2)

round(mean(main_TN_HHS), 2)
round(mean(interaction_TN_HHS), 2)

round(mean(main_FP_HHS), 2)
round(mean(interaction_FP_HHS), 2)

round(mean(main_FN_HHS), 2)
round(mean(interaction_FN_HHS), 2)




round(mean(beta_TP_HHS), 2)
round(mean(gamma_TP_HHS), 2)
round(mean(varphi_TP_HHS), 2)

round(mean(beta_TN_HHS), 2)
round(mean(gamma_TN_HHS), 2)
round(mean(varphi_TN_HHS), 2)

round(mean(beta_FP_HHS), 2)
round(mean(gamma_FP_HHS), 2)
round(mean(varphi_FP_HHS), 2)


round(mean(beta_FN_HHS), 2)
round(mean(gamma_FN_HHS), 2)
round(mean(varphi_FN_HHS), 2)
######################

round(mean(alpha_est_mat_HHS), 2)
round(mean(alpha_cover_per_HHS), 2)

round(mean(beta_cover_per_HHS), 2)
round(mean(beta_cover_nonzero_per_HHS), 2)

round(mean(gamma_cover_per_HHS), 2)
round(mean(gamma_cover_nonzero_per_HHS), 2)

round(mean(varphi_cover_per_HHS), 2)
round(mean(varphi_cover_nonzero_per_HHS), 2)

round(mean(beta_RMSE_HHS)*100, 2)
round(mean(beta_RMSE_nonzero_HHS)*100, 2)
round(mean(beta_RMSE_zero_HHS)*100, 2)

round(mean(gamma_RMSE_HHS)*100, 2)
round(mean(gamma_RMSE_nonzero_HHS)*100, 2)
round(mean(gamma_RMSE_zero_HHS)*100, 2)

round(mean(varphi_RMSE_HHS)*100, 2)
round(mean(varphi_RMSE_nonzero_HHS)*100, 2)
round(mean(varphi_RMSE_zero_HHS)*100, 2)

round(mean(delta_RMSE_HHS)*100, 2)
round(mean(delta_RMSE_nonzero_HHS)*100, 2)
round(mean(delta_RMSE_zero_HHS)*100, 2)

round(mean(beta_sensitivity_HHS)*100, 2)
round(mean(gamma_sensitivity_HHS)*100, 2)
round(mean(varphi_sensitivity_HHS)*100, 2)

round(mean(beta_specificity_HHS)*100, 2)
round(mean(gamma_specificity_HHS)*100, 2)
round(mean(varphi_specificity_HHS)*100, 2)

round(mean(beta_precision_HHS)*100, 2)
round(mean(gamma_precision_HHS)*100, 2)
round(mean(varphi_precision_HHS)*100, 2)



beta_true_array <- simplify2array(beta_true_list_HHS)
beta_true_mean_HHS <- apply(beta_true_array, c(1, 2), mean)
beta_est_array <- simplify2array(beta_est_list_HHS)
beta_est_mean_HHS <- apply(beta_est_array, c(1, 2), mean)
dim(beta_est_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_true_mean_HHS[, k],
    beta_est_mean_HHS[, k],
    xlab = "True beta",
    ylab = "Mean estimated beta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



gamma_true_array <- simplify2array(gamma_true_list_HHS)
gamma_true_mean_HHS <- apply(gamma_true_array, c(1, 2), mean)
gamma_est_array <- simplify2array(gamma_est_list_HHS)
gamma_est_mean_HHS <- apply(gamma_est_array, c(1, 2), mean)
dim(gamma_est_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_true_mean_HHS[, k],
    gamma_est_mean_HHS[, k],
    xlab = "True gamma",
    ylab = "Mean estimated gamma",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



delta_true_array <- simplify2array(delta_true_list_HHS)
delta_true_mean_HHS <- apply(delta_true_array, c(1, 2), mean)
delta_est_array <- simplify2array(delta_est_list_HHS)
delta_est_mean_HHS <- apply(delta_est_array, c(1, 2), mean)
dim(delta_est_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_true_mean_HHS[, k],
    delta_est_mean_HHS[, k],
    xlab = "True delta",
    ylab = "Mean estimated delta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



beta_bias_array <- simplify2array(beta_bias_HHS)
beta_bias_mean_HHS <- apply(beta_bias_array, c(1, 2), mean)
mean(beta_bias_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_mean_HHS[, k],
    xlab = "Modifier",
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_nonzero_array <- simplify2array(beta_bias_nonzero_HHS)
beta_bias_nonzero_mean_HHS <- apply(beta_bias_nonzero_array, c(1, 2), mean)
mean(beta_bias_nonzero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_nonzero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_zero_array <- simplify2array(beta_bias_zero_HHS)
beta_bias_zero_mean_HHS <- apply(beta_bias_zero_array, c(1, 2), mean)
mean(beta_bias_zero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_zero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}





gamma_bias_array <- simplify2array(gamma_bias_HHS)
gamma_bias_mean_HHS <- apply(gamma_bias_array, c(1, 2), mean)
mean(gamma_bias_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_mean_HHS[, k],
    xlab = "Exposure",
    main = paste("Outcome", k),
    pch = 19
  )
}

gamma_bias_nonzero_array <- simplify2array(gamma_bias_nonzero_HHS)
gamma_bias_nonzero_mean_HHS <- apply(gamma_bias_nonzero_array, c(1, 2), mean)
mean(gamma_bias_nonzero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_nonzero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_zero_array <- simplify2array(gamma_bias_zero_HHS)
gamma_bias_zero_mean_HHS <- apply(gamma_bias_zero_array, c(1, 2), mean)
mean(gamma_bias_zero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_zero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_array <- simplify2array(varphi_bias_HHS)
varphi_bias_mean_HHS <- apply(varphi_bias_array, c(1, 2), mean)
mean(varphi_bias_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_mean_HHS[, k],
    xlab = "Covariate",
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_nonzero_array <- simplify2array(varphi_bias_nonzero_HHS)
varphi_bias_nonzero_mean_HHS <- apply(varphi_bias_nonzero_array, c(1, 2), mean)
mean(varphi_bias_nonzero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_nonzero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_zero_array <- simplify2array(varphi_bias_zero_HHS)
varphi_bias_zero_mean_HHS <- apply(varphi_bias_zero_array, c(1, 2), mean)
mean(varphi_bias_zero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_zero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_array <- simplify2array(delta_bias_HHS)
delta_bias_mean_HHS <- apply(delta_bias_array, c(1, 2), mean)
mean(delta_bias_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_mean_HHS[, k],
    xlab = "Interaction",
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_nonzero_array <- simplify2array(delta_bias_nonzero_HHS)
delta_bias_nonzero_mean_HHS <- apply(delta_bias_nonzero_array, c(1, 2), mean)
mean(delta_bias_nonzero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_nonzero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_zero_array <- simplify2array(delta_bias_zero_HHS)
delta_bias_zero_mean_HHS <- apply(delta_bias_zero_array, c(1, 2), mean)
mean(delta_bias_zero_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_zero_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}





#######################################################

# HHSSI model

# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_0weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_10weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_25weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_50weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_100weak.rda")


# load("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_0weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_10weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_25weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_50weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_100weak.rda")


# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")



nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]
J <- dim(sim_data_list[[1]]$x)[2]
M <- dim(sim_data_list[[1]]$z)[2]


alpha_cover_per_HHSSI <- c()
beta_cover_per_HHSSI <- c()
gamma_cover_per_HHSSI <- c()
varphi_cover_per_HHSSI <- c()
delta_cover_per_HHSSI <- c()

delta_RMSE_HHSSI <- c()
RMSE_all_main_effects_HHSSI <- c()

beta_cover_nonzero_per_HHSSI <- c()
gamma_cover_nonzero_per_HHSSI <- c()
varphi_cover_nonzero_per_HHSSI <- c()
delta_cover_nonzero_per_HHSSI <- c()

alpha_est_mat_HHSSI <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_HHSSI <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_HHSSI <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_HHSSI <- matrix(NA, nrow = nsample, ncol = K)

beta_est_list_HHSSI <- list()
beta_true_list_HHSSI <- list()
beta_lower_list_HHSSI <- list()
beta_upper_list_HHSSI <- list()

gamma_est_list_HHSSI <- list()
gamma_true_list_HHSSI <- list()
gamma_lower_list_HHSSI <- list()
gamma_upper_list_HHSSI <- list()

varphi_est_list_HHSSI <- list()
varphi_true_list_HHSSI <- list()
varphi_lower_list_HHSSI <- list()
varphi_upper_list_HHSSI <- list()

delta_est_list_HHSSI <- list()
delta_true_list_HHSSI <- list()
delta_lower_list_HHSSI <- list()
delta_upper_list_HHSSI <- list()

main_sensitivity_HHSSI <- c()
main_specificity_HHSSI <- c()
main_precision_HHSSI <- c()

interaction_sensitivity_HHSSI <- c()
interaction_specificity_HHSSI <- c()
interaction_precision_HHSSI <- c()

main_TP_HHSSI <- c()
main_TN_HHSSI <- c()
main_FP_HHSSI <- c()
main_FN_HHSSI <- c()

interaction_TP_HHSSI <- c()
interaction_TN_HHSSI <- c()
interaction_FP_HHSSI <- c()
interaction_FN_HHSSI <- c()


# ============================================================
# Additional diagnostics
# ============================================================

beta_bias_HHSSI <- c()
gamma_bias_HHSSI <- c()
varphi_bias_HHSSI <- c()
delta_bias_HHSSI <- c()

beta_bias_nonzero_HHSSI <- c()
gamma_bias_nonzero_HHSSI <- c()
varphi_bias_nonzero_HHSSI <- c()
delta_bias_nonzero_HHSSI <- c()

beta_bias_zero_HHSSI <- c()
gamma_bias_zero_HHSSI <- c()
varphi_bias_zero_HHSSI <- c()
delta_bias_zero_HHSSI <- c()

beta_RMSE_HHSSI <- c()
gamma_RMSE_HHSSI <- c()
varphi_RMSE_HHSSI <- c()
delta_RMSE_HHSSI <- c()

beta_RMSE_nonzero_HHSSI <- c()
gamma_RMSE_nonzero_HHSSI <- c()
varphi_RMSE_nonzero_HHSSI <- c()
delta_RMSE_nonzero_HHSSI <- c()

beta_RMSE_zero_HHSSI <- c()
gamma_RMSE_zero_HHSSI <- c()
varphi_RMSE_zero_HHSSI <- c()
delta_RMSE_zero_HHSSI <- c()

beta_TP_HHSSI <- c()
gamma_TP_HHSSI <- c()
varphi_TP_HHSSI <- c()

beta_TN_HHSSI <- c()
gamma_TN_HHSSI <- c()
varphi_TN_HHSSI <- c()

beta_FP_HHSSI <- c()
gamma_FP_HHSSI <- c()
varphi_FP_HHSSI <- c()

beta_FN_HHSSI <- c()
gamma_FN_HHSSI <- c()
varphi_FN_HHSSI <- c()

beta_sensitivity_HHSSI <- c()
gamma_sensitivity_HHSSI <- c()
varphi_sensitivity_HHSSI <- c()

beta_specificity_HHSSI <- c()
gamma_specificity_HHSSI <- c()
varphi_specificity_HHSSI <- c()

beta_precision_HHSSI <- c()
gamma_precision_HHSSI <- c()
varphi_precision_HHSSI <- c()

delta_RMSE_indices_mat_HHSSI <- matrix(NA, nrow = nsample, ncol = J*M)
delta_cover_indices_mat_HHSSI <- matrix(NA, nrow = nsample, ncol = J*M)
delta_var_indices_mat_HHSSI <- matrix(NA, nrow = nsample, ncol = J*M)

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  d <- sim_data_list[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_HHSSI[sample, ] <- alpha_true
  alpha_est_mat_HHSSI[sample, ] <- apply(theta_update_HHSSI_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_HHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_HHSSI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_HHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_HHSSI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_HHSSI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_HHSSI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_HHSSI_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_HHSSI_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_HHSSI_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_HHSSI[[sample]] <- beta_est
  beta_lower_list_HHSSI[[sample]] <- beta_lower
  beta_upper_list_HHSSI[[sample]] <- beta_upper
  beta_cover_per_HHSSI[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_HHSSI[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # Additional diagnostics
  nonzero_beta <- beta_true != 0
  zero_beta <- beta_true == 0

  beta_bias <- beta_est - beta_true

  beta_bias_HHSSI[[sample]] <- beta_bias
  beta_bias_nonzero_HHSSI[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_zero_HHSSI[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_HHSSI[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_nonzero_HHSSI[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_zero_HHSSI[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_HHSSI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_HHSSI_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_HHSSI_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_HHSSI_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_HHSSI[[sample]] <- gamma_est
  gamma_lower_list_HHSSI[[sample]] <- gamma_lower
  gamma_upper_list_HHSSI[[sample]] <- gamma_upper
  gamma_cover_per_HHSSI[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_HHSSI[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  # Additional diagnostics
  nonzero_gamma <- gamma_true != 0
  zero_gamma <- gamma_true == 0

  gamma_bias <- gamma_est - gamma_true

  gamma_bias_HHSSI[[sample]] <- gamma_bias
  gamma_bias_nonzero_HHSSI[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_zero_HHSSI[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_HHSSI[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_nonzero_HHSSI[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_zero_HHSSI[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

  # varphi
  varphi_true <- sim_data_list[[sample]]$varphi_true
  varphi_true_list_HHSSI[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list[[sample]][, o, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    varphi_cover[(o-1-J-M-J*M), ] <- (varphi_true[(o-1-J-M-J*M), ] >= varphi_lower[(o-1-J-M-J*M), ]) &
      (varphi_true[(o-1-J-M-J*M), ] <= varphi_upper[(o-1-J-M-J*M), ])
    varphi_MSE[o-1-J-M-J*M, ] <- (varphi_est[(o-1-J-M-J*M), ] - varphi_true[o-1-J-M-J*M, ])^2
  }
  varphi_est_list_HHSSI[[sample]] <- varphi_est
  varphi_lower_list_HHSSI[[sample]] <- varphi_lower
  varphi_upper_list_HHSSI[[sample]] <- varphi_upper
  varphi_cover_per_HHSSI[sample] <- (sum(varphi_cover)/(dim(varphi_true)[1]*dim(varphi_true)[2]))*100
  nonzero_idx_varphi <- varphi_true != 0
  varphi_cover_nonzero_per_HHSSI[sample] <- (sum(varphi_cover[nonzero_idx_varphi]) / sum(nonzero_idx_varphi)) * 100

  # Additional diagnostics
  nonzero_varphi <- varphi_true != 0
  zero_varphi <- varphi_true == 0

  varphi_bias <- varphi_est - varphi_true

  varphi_bias_HHSSI[[sample]] <- varphi_bias
  varphi_bias_nonzero_HHSSI[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_zero_HHSSI[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_HHSSI[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_nonzero_HHSSI[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_zero_HHSSI[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


  main_MSE <- c(beta_MSE, gamma_MSE, varphi_MSE)
  RMSE_all_main_effects_HHSSI[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_HHSSI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_HHSSI[[sample]] <- delta_est
  delta_lower_list_HHSSI[[sample]] <- delta_lower
  delta_upper_list_HHSSI[[sample]] <- delta_upper
  delta_cover_per_HHSSI[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_HHSSI[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_HHSSI[sample] <- sqrt(mean(delta_MSE))

  delta_cover_indices_mat_HHSSI[sample, ] <- rowMeans(delta_cover)*100
  delta_RMSE_indices_mat_HHSSI[sample, ] <- sqrt(rowMeans(delta_MSE))
  delta_var_indices_mat_HHSSI[sample, ] <- rowMeans(delta_var)

  # Additional diagnostics
  nonzero_delta <- delta_true != 0
  zero_delta <- delta_true == 0

  delta_bias <- delta_est - delta_true

  delta_bias_HHSSI[[sample]] <- delta_bias
  delta_bias_nonzero_HHSSI[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_zero_HHSSI[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_HHSSI[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_nonzero_HHSSI[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_zero_HHSSI[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_varphi <- as.vector(varphi_true != 0)
  selected_varphi <- as.vector((varphi_lower > 0) | (varphi_upper < 0))

  truth_main <- c(truth_beta, truth_gamma, truth_varphi)
  selected_main <- c(selected_beta, selected_gamma, selected_varphi)

  main_TP_HHSSI[sample] <- sum(truth_main == 1 & selected_main == 1)
  main_FN_HHSSI[sample] <- sum(truth_main == 1 & selected_main == 0)
  main_TN_HHSSI[sample] <- sum(truth_main == 0 & selected_main == 0)
  main_FP_HHSSI[sample] <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_HHSSI[sample] <- (main_TP_HHSSI[sample] /
                                        (main_TP_HHSSI[sample] + main_FN_HHSSI[sample]))
  main_specificity_HHSSI[sample] <- (main_TN_HHSSI[sample] /
                                        (main_TN_HHSSI[sample] + main_FP_HHSSI[sample]))
  main_precision_HHSSI[sample] <- (main_TP_HHSSI[sample] /
                                      (main_TP_HHSSI[sample] + main_FP_HHSSI[sample]))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP_HHSSI[sample] <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN_HHSSI[sample] <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN_HHSSI[sample] <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP_HHSSI[sample] <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_HHSSI[sample] <- (interaction_TP_HHSSI[sample] /
                                               (interaction_TP_HHSSI[sample] + interaction_FN_HHSSI[sample]))
  interaction_specificity_HHSSI[sample] <- (interaction_TN_HHSSI[sample] /
                                               (interaction_TN_HHSSI[sample] + interaction_FP_HHSSI[sample]))
  interaction_precision_HHSSI[sample] <- (interaction_TP_HHSSI[sample] /
                                             (interaction_TP_HHSSI[sample] + interaction_FP_HHSSI[sample]))

  # Additional diaagnostics
  beta_TP_HHSSI[sample] <- sum(truth_beta == 1 & selected_beta == 1)
  beta_FN_HHSSI[sample] <- sum(truth_beta == 1 & selected_beta == 0)
  beta_TN_HHSSI[sample] <- sum(truth_beta == 0 & selected_beta == 0)
  beta_FP_HHSSI[sample] <- sum(truth_beta == 0 & selected_beta == 1)

  gamma_TP_HHSSI[sample] <- sum(truth_gamma == 1 & selected_gamma == 1)
  gamma_FN_HHSSI[sample] <- sum(truth_gamma == 1 & selected_gamma == 0)
  gamma_TN_HHSSI[sample] <- sum(truth_gamma == 0 & selected_gamma == 0)
  gamma_FP_HHSSI[sample] <- sum(truth_gamma == 0 & selected_gamma == 1)

  varphi_TP_HHSSI[sample] <- sum(truth_varphi == 1 & selected_varphi == 1)
  varphi_FN_HHSSI[sample] <- sum(truth_varphi == 1 & selected_varphi == 0)
  varphi_TN_HHSSI[sample] <- sum(truth_varphi == 0 & selected_varphi == 0)
  varphi_FP_HHSSI[sample] <- sum(truth_varphi == 0 & selected_varphi == 1)

  beta_sensitivity_HHSSI[sample] <- (beta_TP_HHSSI[sample] /
                                        (beta_TP_HHSSI[sample] + beta_FN_HHSSI[sample]))
  beta_specificity_HHSSI[sample] <- (beta_TN_HHSSI[sample] /
                                        (beta_TN_HHSSI[sample] + beta_FP_HHSSI[sample]))
  beta_precision_HHSSI[sample] <- (beta_TP_HHSSI[sample] /
                                      (beta_TP_HHSSI[sample] + beta_FP_HHSSI[sample]))

  gamma_sensitivity_HHSSI[sample] <- (gamma_TP_HHSSI[sample] /
                                         (gamma_TP_HHSSI[sample] + gamma_FN_HHSSI[sample]))
  gamma_specificity_HHSSI[sample] <- (gamma_TN_HHSSI[sample] /
                                         (gamma_TN_HHSSI[sample] + gamma_FP_HHSSI[sample]))
  gamma_precision_HHSSI[sample] <- (gamma_TP_HHSSI[sample] /
                                       (gamma_TP_HHSSI[sample] + gamma_FP_HHSSI[sample]))

  varphi_sensitivity_HHSSI[sample] <- (varphi_TP_HHSSI[sample] /
                                          (varphi_TP_HHSSI[sample] + varphi_FN_HHSSI[sample]))
  varphi_specificity_HHSSI[sample] <- (varphi_TN_HHSSI[sample] /
                                          (varphi_TN_HHSSI[sample] + varphi_FP_HHSSI[sample]))
  varphi_precision_HHSSI[sample] <- (varphi_TP_HHSSI[sample] /
                                        (varphi_TP_HHSSI[sample] + varphi_FP_HHSSI[sample]))
}




round(mean(RMSE_all_main_effects_HHSSI)*100, 2)
round(mean(delta_RMSE_HHSSI)*100, 2)

round((mean(beta_cover_per_HHSSI) + mean(gamma_cover_per_HHSSI) +
         mean(varphi_cover_per_HHSSI))/3, 2)
round(mean(delta_cover_per_HHSSI), 2)

round((mean(beta_cover_nonzero_per_HHSSI) + mean(gamma_cover_nonzero_per_HHSSI) +
         mean(varphi_cover_nonzero_per_HHSSI))/3, 2)
round(mean(delta_cover_nonzero_per_HHSSI), 2)

round(mean(main_sensitivity_HHSSI) * 100, 2)
round(mean(interaction_sensitivity_HHSSI) * 100, 2)

round(mean(main_specificity_HHSSI) * 100, 2)
round(mean(interaction_specificity_HHSSI) * 100, 2)

round(mean(main_precision_HHSSI) * 100, 2)
round(mean(interaction_precision_HHSSI) * 100, 2)


dim(delta_RMSE_indices_mat_HHSSI)
###############################

round(mean(main_TP_HHSSI), 2)
round(mean(interaction_TP_HHSSI), 2)

round(mean(main_TN_HHSSI), 2)
round(mean(interaction_TN_HHSSI), 2)

round(mean(main_FP_HHSSI), 2)
round(mean(interaction_FP_HHSSI), 2)

round(mean(main_FN_HHSSI), 2)
round(mean(interaction_FN_HHSSI), 2)



round(mean(beta_TP_HHSSI), 2)
round(mean(gamma_TP_HHSSI), 2)
round(mean(varphi_TP_HHSSI), 2)

round(mean(beta_TN_HHSSI), 2)
round(mean(gamma_TN_HHSSI), 2)
round(mean(varphi_TN_HHSSI), 2)

round(mean(beta_FP_HHSSI), 2)
round(mean(gamma_FP_HHSSI), 2)
round(mean(varphi_FP_HHSSI), 2)


round(mean(beta_FN_HHSSI), 2)
round(mean(gamma_FN_HHSSI), 2)
round(mean(varphi_FN_HHSSI), 2)
######################

round(mean(alpha_est_mat_HHSSI), 2)
round(mean(alpha_cover_per_HHSSI), 2)

round(mean(beta_cover_per_HHSSI), 2)
round(mean(beta_cover_nonzero_per_HHSSI), 2)

round(mean(gamma_cover_per_HHSSI), 2)
round(mean(gamma_cover_nonzero_per_HHSSI), 2)

round(mean(varphi_cover_per_HHSSI), 2)
round(mean(varphi_cover_nonzero_per_HHSSI), 2)

round(mean(beta_RMSE_HHSSI)*100, 2)
round(mean(beta_RMSE_nonzero_HHSSI)*100, 2)
round(mean(beta_RMSE_zero_HHSSI)*100, 2)

round(mean(gamma_RMSE_HHSSI)*100, 2)
round(mean(gamma_RMSE_nonzero_HHSSI)*100, 2)
round(mean(gamma_RMSE_zero_HHSSI)*100, 2)

round(mean(varphi_RMSE_HHSSI)*100, 2)
round(mean(varphi_RMSE_nonzero_HHSSI)*100, 2)
round(mean(varphi_RMSE_zero_HHSSI)*100, 2)

round(mean(delta_RMSE_HHSSI)*100, 2)
round(mean(delta_RMSE_nonzero_HHSSI)*100, 2)
round(mean(delta_RMSE_zero_HHSSI)*100, 2)

round(mean(beta_sensitivity_HHSSI)*100, 2)
round(mean(gamma_sensitivity_HHSSI)*100, 2)
round(mean(varphi_sensitivity_HHSSI)*100, 2)

round(mean(beta_specificity_HHSSI)*100, 2)
round(mean(gamma_specificity_HHSSI)*100, 2)
round(mean(varphi_specificity_HHSSI)*100, 2)

round(mean(beta_precision_HHSSI)*100, 2)
round(mean(gamma_precision_HHSSI)*100, 2)
round(mean(varphi_precision_HHSSI)*100, 2)



beta_true_array <- simplify2array(beta_true_list_HHSSI)
beta_true_mean_HHSSI <- apply(beta_true_array, c(1, 2), mean)
beta_est_array <- simplify2array(beta_est_list_HHSSI)
beta_est_mean_HHSSI <- apply(beta_est_array, c(1, 2), mean)
dim(beta_est_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_true_mean_HHSSI[, k],
    beta_est_mean_HHSSI[, k],
    xlab = "True beta",
    ylab = "Mean estimated beta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



gamma_true_array <- simplify2array(gamma_true_list_HHSSI)
gamma_true_mean_HHSSI <- apply(gamma_true_array, c(1, 2), mean)
gamma_est_array <- simplify2array(gamma_est_list_HHSSI)
gamma_est_mean_HHSSI <- apply(gamma_est_array, c(1, 2), mean)
dim(gamma_est_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_true_mean_HHSSI[, k],
    gamma_est_mean_HHSSI[, k],
    xlab = "True gamma",
    ylab = "Mean estimated gamma",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



delta_true_array <- simplify2array(delta_true_list_HHSSI)
delta_true_mean_HHSSI <- apply(delta_true_array, c(1, 2), mean)
delta_est_array <- simplify2array(delta_est_list_HHSSI)
delta_est_mean_HHSSI <- apply(delta_est_array, c(1, 2), mean)
dim(delta_est_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_true_mean_HHSSI[, k],
    delta_est_mean_HHSSI[, k],
    xlab = "True delta",
    ylab = "Mean estimated delta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



beta_bias_array <- simplify2array(beta_bias_HHSSI)
beta_bias_mean_HHSSI <- apply(beta_bias_array, c(1, 2), mean)
mean(beta_bias_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_mean_HHSSI[, k],
    xlab = "Modifier",
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_nonzero_array <- simplify2array(beta_bias_nonzero_HHSSI)
beta_bias_nonzero_mean_HHSSI <- apply(beta_bias_nonzero_array, c(1, 2), mean)
mean(beta_bias_nonzero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_nonzero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_zero_array <- simplify2array(beta_bias_zero_HHSSI)
beta_bias_zero_mean_HHSSI <- apply(beta_bias_zero_array, c(1, 2), mean)
mean(beta_bias_zero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_zero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}





gamma_bias_array <- simplify2array(gamma_bias_HHSSI)
gamma_bias_mean_HHSSI <- apply(gamma_bias_array, c(1, 2), mean)
mean(gamma_bias_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_mean_HHSSI[, k],
    xlab = "Exposure",
    main = paste("Outcome", k),
    pch = 19
  )
}

gamma_bias_nonzero_array <- simplify2array(gamma_bias_nonzero_HHSSI)
gamma_bias_nonzero_mean_HHSSI <- apply(gamma_bias_nonzero_array, c(1, 2), mean)
mean(gamma_bias_nonzero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_nonzero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_zero_array <- simplify2array(gamma_bias_zero_HHSSI)
gamma_bias_zero_mean_HHSSI <- apply(gamma_bias_zero_array, c(1, 2), mean)
mean(gamma_bias_zero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_zero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_array <- simplify2array(varphi_bias_HHSSI)
varphi_bias_mean_HHSSI <- apply(varphi_bias_array, c(1, 2), mean)
mean(varphi_bias_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_mean_HHSSI[, k],
    xlab = "Covariate",
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_nonzero_array <- simplify2array(varphi_bias_nonzero_HHSSI)
varphi_bias_nonzero_mean_HHSSI <- apply(varphi_bias_nonzero_array, c(1, 2), mean)
mean(varphi_bias_nonzero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_nonzero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_zero_array <- simplify2array(varphi_bias_zero_HHSSI)
varphi_bias_zero_mean_HHSSI <- apply(varphi_bias_zero_array, c(1, 2), mean)
mean(varphi_bias_zero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_zero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_array <- simplify2array(delta_bias_HHSSI)
delta_bias_mean_HHSSI <- apply(delta_bias_array, c(1, 2), mean)
mean(delta_bias_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_mean_HHSSI[, k],
    xlab = "Interaction",
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_nonzero_array <- simplify2array(delta_bias_nonzero_HHSSI)
delta_bias_nonzero_mean_HHSSI <- apply(delta_bias_nonzero_array, c(1, 2), mean)
mean(delta_bias_nonzero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_nonzero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_zero_array <- simplify2array(delta_bias_zero_HHSSI)
delta_bias_zero_mean_HHSSI <- apply(delta_bias_zero_array, c(1, 2), mean)
mean(delta_bias_zero_mean_HHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_zero_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}





######################################
# NHHS model

# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_0weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_10weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_25weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_50weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_100weak.rda")


# load("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_0weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_10weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_25weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_50weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_100weak.rda")


# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")



nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]
J <- dim(sim_data_list[[1]]$x)[2]
M <- dim(sim_data_list[[1]]$z)[2]

alpha_cover_per_NHHS <- c()
beta_cover_per_NHHS <- c()
gamma_cover_per_NHHS <- c()
varphi_cover_per_NHHS <- c()
delta_cover_per_NHHS <- c()

delta_RMSE_NHHS <- c()
RMSE_all_main_effects_NHHS <- c()

beta_cover_nonzero_per_NHHS <- c()
gamma_cover_nonzero_per_NHHS <- c()
varphi_cover_nonzero_per_NHHS <- c()
delta_cover_nonzero_per_NHHS <- c()

alpha_est_mat_NHHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_NHHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_NHHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_NHHS <- matrix(NA, nrow = nsample, ncol = K)

beta_est_list_NHHS <- list()
beta_true_list_NHHS <- list()
beta_lower_list_NHHS <- list()
beta_upper_list_NHHS <- list()

gamma_est_list_NHHS <- list()
gamma_true_list_NHHS <- list()
gamma_lower_list_NHHS <- list()
gamma_upper_list_NHHS <- list()

varphi_est_list_NHHS <- list()
varphi_true_list_NHHS <- list()
varphi_lower_list_NHHS <- list()
varphi_upper_list_NHHS <- list()

delta_est_list_NHHS <- list()
delta_true_list_NHHS <- list()
delta_lower_list_NHHS <- list()
delta_upper_list_NHHS <- list()

main_sensitivity_NHHS <- c()
main_specificity_NHHS <- c()
main_precision_NHHS <- c()

interaction_sensitivity_NHHS <- c()
interaction_specificity_NHHS <- c()
interaction_precision_NHHS <- c()

main_TP_NHHS <- c()
main_TN_NHHS <- c()
main_FP_NHHS <- c()
main_FN_NHHS <- c()

interaction_TP_NHHS <- c()
interaction_TN_NHHS <- c()
interaction_FP_NHHS <- c()
interaction_FN_NHHS <- c()


# ============================================================
# Additional diagnostics
# ============================================================

beta_bias_NHHS <- c()
gamma_bias_NHHS <- c()
varphi_bias_NHHS <- c()
delta_bias_NHHS <- c()

beta_bias_nonzero_NHHS <- c()
gamma_bias_nonzero_NHHS <- c()
varphi_bias_nonzero_NHHS <- c()
delta_bias_nonzero_NHHS <- c()

beta_bias_zero_NHHS <- c()
gamma_bias_zero_NHHS <- c()
varphi_bias_zero_NHHS <- c()
delta_bias_zero_NHHS <- c()

beta_RMSE_NHHS <- c()
gamma_RMSE_NHHS <- c()
varphi_RMSE_NHHS <- c()
delta_RMSE_NHHS <- c()

beta_RMSE_nonzero_NHHS <- c()
gamma_RMSE_nonzero_NHHS <- c()
varphi_RMSE_nonzero_NHHS <- c()
delta_RMSE_nonzero_NHHS <- c()

beta_RMSE_zero_NHHS <- c()
gamma_RMSE_zero_NHHS <- c()
varphi_RMSE_zero_NHHS <- c()
delta_RMSE_zero_NHHS <- c()

beta_TP_NHHS <- c()
gamma_TP_NHHS <- c()
varphi_TP_NHHS <- c()

beta_TN_NHHS <- c()
gamma_TN_NHHS <- c()
varphi_TN_NHHS <- c()

beta_FP_NHHS <- c()
gamma_FP_NHHS <- c()
varphi_FP_NHHS <- c()

beta_FN_NHHS <- c()
gamma_FN_NHHS <- c()
varphi_FN_NHHS <- c()

beta_sensitivity_NHHS <- c()
gamma_sensitivity_NHHS <- c()
varphi_sensitivity_NHHS <- c()

beta_specificity_NHHS <- c()
gamma_specificity_NHHS <- c()
varphi_specificity_NHHS <- c()

beta_precision_NHHS <- c()
gamma_precision_NHHS <- c()
varphi_precision_NHHS <- c()

delta_RMSE_indices_mat_NHHS <- matrix(NA, nrow = nsample, ncol = J*M)
delta_cover_indices_mat_NHHS <- matrix(NA, nrow = nsample, ncol = J*M)
delta_var_indices_mat_NHHS <- matrix(NA, nrow = nsample, ncol = J*M)

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  d <- sim_data_list[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_NHHS[sample, ] <- alpha_true
  alpha_est_mat_NHHS[sample, ] <- apply(theta_update_NHHS_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NHHS_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NHHS[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NHHS_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NHHS[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NHHS[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_NHHS[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NHHS_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NHHS_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NHHS_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_NHHS[[sample]] <- beta_est
  beta_lower_list_NHHS[[sample]] <- beta_lower
  beta_upper_list_NHHS[[sample]] <- beta_upper
  beta_cover_per_NHHS[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_NHHS[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # Additional diagnostics
  nonzero_beta <- beta_true != 0
  zero_beta <- beta_true == 0

  beta_bias <- beta_est - beta_true

  beta_bias_NHHS[[sample]] <- beta_bias
  beta_bias_nonzero_NHHS[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_zero_NHHS[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_NHHS[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_nonzero_NHHS[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_zero_NHHS[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_NHHS[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_NHHS[[sample]] <- gamma_est
  gamma_lower_list_NHHS[[sample]] <- gamma_lower
  gamma_upper_list_NHHS[[sample]] <- gamma_upper
  gamma_cover_per_NHHS[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_NHHS[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  # Additional diagnostics
  nonzero_gamma <- gamma_true != 0
  zero_gamma <- gamma_true == 0

  gamma_bias <- gamma_est - gamma_true

  gamma_bias_NHHS[[sample]] <- gamma_bias
  gamma_bias_nonzero_NHHS[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_zero_NHHS[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_NHHS[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_nonzero_NHHS[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_zero_NHHS[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

  # varphi
  varphi_true <- sim_data_list[[sample]]$varphi_true
  varphi_true_list_NHHS[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_NHHS_list[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_NHHS_list[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_NHHS_list[[sample]][, o, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    varphi_cover[(o-1-J-M-J*M), ] <- (varphi_true[(o-1-J-M-J*M), ] >= varphi_lower[(o-1-J-M-J*M), ]) &
      (varphi_true[(o-1-J-M-J*M), ] <= varphi_upper[(o-1-J-M-J*M), ])
    varphi_MSE[o-1-J-M-J*M, ] <- (varphi_est[(o-1-J-M-J*M), ] - varphi_true[o-1-J-M-J*M, ])^2
  }
  varphi_est_list_NHHS[[sample]] <- varphi_est
  varphi_lower_list_NHHS[[sample]] <- varphi_lower
  varphi_upper_list_NHHS[[sample]] <- varphi_upper
  varphi_cover_per_NHHS[sample] <- (sum(varphi_cover)/(dim(varphi_true)[1]*dim(varphi_true)[2]))*100
  nonzero_idx_varphi <- varphi_true != 0
  varphi_cover_nonzero_per_NHHS[sample] <- (sum(varphi_cover[nonzero_idx_varphi]) / sum(nonzero_idx_varphi)) * 100

  # Additional diagnostics
  nonzero_varphi <- varphi_true != 0
  zero_varphi <- varphi_true == 0

  varphi_bias <- varphi_est - varphi_true

  varphi_bias_NHHS[[sample]] <- varphi_bias
  varphi_bias_nonzero_NHHS[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_zero_NHHS[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_NHHS[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_nonzero_NHHS[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_zero_NHHS[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


  main_MSE <- c(beta_MSE, gamma_MSE, varphi_MSE)
  RMSE_all_main_effects_NHHS[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NHHS[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_NHHS[[sample]] <- delta_est
  delta_lower_list_NHHS[[sample]] <- delta_lower
  delta_upper_list_NHHS[[sample]] <- delta_upper
  delta_cover_per_NHHS[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_NHHS[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_NHHS[sample] <- sqrt(mean(delta_MSE))

  delta_cover_indices_mat_NHHS[sample, ] <- rowMeans(delta_cover)*100
  delta_RMSE_indices_mat_NHHS[sample, ] <- sqrt(rowMeans(delta_MSE))
  delta_var_indices_mat_NHHS[sample, ] <- rowMeans(delta_var)

  # Additional diagnostics
  nonzero_delta <- delta_true != 0
  zero_delta <- delta_true == 0

  delta_bias <- delta_est - delta_true

  delta_bias_NHHS[[sample]] <- delta_bias
  delta_bias_nonzero_NHHS[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_zero_NHHS[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_NHHS[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_nonzero_NHHS[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_zero_NHHS[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_varphi <- as.vector(varphi_true != 0)
  selected_varphi <- as.vector((varphi_lower > 0) | (varphi_upper < 0))

  truth_main <- c(truth_beta, truth_gamma, truth_varphi)
  selected_main <- c(selected_beta, selected_gamma, selected_varphi)

  main_TP_NHHS[sample] <- sum(truth_main == 1 & selected_main == 1)
  main_FN_NHHS[sample] <- sum(truth_main == 1 & selected_main == 0)
  main_TN_NHHS[sample] <- sum(truth_main == 0 & selected_main == 0)
  main_FP_NHHS[sample] <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_NHHS[sample] <- (main_TP_NHHS[sample] /
                                        (main_TP_NHHS[sample] + main_FN_NHHS[sample]))
  main_specificity_NHHS[sample] <- (main_TN_NHHS[sample] /
                                        (main_TN_NHHS[sample] + main_FP_NHHS[sample]))
  main_precision_NHHS[sample] <- (main_TP_NHHS[sample] /
                                      (main_TP_NHHS[sample] + main_FP_NHHS[sample]))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP_NHHS[sample] <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN_NHHS[sample] <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN_NHHS[sample] <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP_NHHS[sample] <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NHHS[sample] <- (interaction_TP_NHHS[sample] /
                                               (interaction_TP_NHHS[sample] + interaction_FN_NHHS[sample]))
  interaction_specificity_NHHS[sample] <- (interaction_TN_NHHS[sample] /
                                               (interaction_TN_NHHS[sample] + interaction_FP_NHHS[sample]))
  interaction_precision_NHHS[sample] <- (interaction_TP_NHHS[sample] /
                                             (interaction_TP_NHHS[sample] + interaction_FP_NHHS[sample]))

  # Additional diaagnostics
  beta_TP_NHHS[sample] <- sum(truth_beta == 1 & selected_beta == 1)
  beta_FN_NHHS[sample] <- sum(truth_beta == 1 & selected_beta == 0)
  beta_TN_NHHS[sample] <- sum(truth_beta == 0 & selected_beta == 0)
  beta_FP_NHHS[sample] <- sum(truth_beta == 0 & selected_beta == 1)

  gamma_TP_NHHS[sample] <- sum(truth_gamma == 1 & selected_gamma == 1)
  gamma_FN_NHHS[sample] <- sum(truth_gamma == 1 & selected_gamma == 0)
  gamma_TN_NHHS[sample] <- sum(truth_gamma == 0 & selected_gamma == 0)
  gamma_FP_NHHS[sample] <- sum(truth_gamma == 0 & selected_gamma == 1)

  varphi_TP_NHHS[sample] <- sum(truth_varphi == 1 & selected_varphi == 1)
  varphi_FN_NHHS[sample] <- sum(truth_varphi == 1 & selected_varphi == 0)
  varphi_TN_NHHS[sample] <- sum(truth_varphi == 0 & selected_varphi == 0)
  varphi_FP_NHHS[sample] <- sum(truth_varphi == 0 & selected_varphi == 1)

  beta_sensitivity_NHHS[sample] <- (beta_TP_NHHS[sample] /
                                        (beta_TP_NHHS[sample] + beta_FN_NHHS[sample]))
  beta_specificity_NHHS[sample] <- (beta_TN_NHHS[sample] /
                                        (beta_TN_NHHS[sample] + beta_FP_NHHS[sample]))
  beta_precision_NHHS[sample] <- (beta_TP_NHHS[sample] /
                                      (beta_TP_NHHS[sample] + beta_FP_NHHS[sample]))

  gamma_sensitivity_NHHS[sample] <- (gamma_TP_NHHS[sample] /
                                         (gamma_TP_NHHS[sample] + gamma_FN_NHHS[sample]))
  gamma_specificity_NHHS[sample] <- (gamma_TN_NHHS[sample] /
                                         (gamma_TN_NHHS[sample] + gamma_FP_NHHS[sample]))
  gamma_precision_NHHS[sample] <- (gamma_TP_NHHS[sample] /
                                       (gamma_TP_NHHS[sample] + gamma_FP_NHHS[sample]))

  varphi_sensitivity_NHHS[sample] <- (varphi_TP_NHHS[sample] /
                                          (varphi_TP_NHHS[sample] + varphi_FN_NHHS[sample]))
  varphi_specificity_NHHS[sample] <- (varphi_TN_NHHS[sample] /
                                          (varphi_TN_NHHS[sample] + varphi_FP_NHHS[sample]))
  varphi_precision_NHHS[sample] <- (varphi_TP_NHHS[sample] /
                                        (varphi_TP_NHHS[sample] + varphi_FP_NHHS[sample]))
}




round(mean(RMSE_all_main_effects_NHHS)*100, 2)
round(mean(delta_RMSE_NHHS)*100, 2)

round((mean(beta_cover_per_NHHS) + mean(gamma_cover_per_NHHS) +
         mean(varphi_cover_per_NHHS))/3, 2)
round(mean(delta_cover_per_NHHS), 2)

round((mean(beta_cover_nonzero_per_NHHS) + mean(gamma_cover_nonzero_per_NHHS) +
         mean(varphi_cover_nonzero_per_NHHS))/3, 2)
round(mean(delta_cover_nonzero_per_NHHS), 2)

round(mean(main_sensitivity_NHHS) * 100, 2)
round(mean(interaction_sensitivity_NHHS) * 100, 2)

round(mean(main_specificity_NHHS) * 100, 2)
round(mean(interaction_specificity_NHHS) * 100, 2)

round(mean(main_precision_NHHS) * 100, 2)
round(mean(interaction_precision_NHHS) * 100, 2)


###############################

round(mean(main_TP_NHHS), 2)
round(mean(interaction_TP_NHHS), 2)

round(mean(main_TN_NHHS), 2)
round(mean(interaction_TN_NHHS), 2)

round(mean(main_FP_NHHS), 2)
round(mean(interaction_FP_NHHS), 2)

round(mean(main_FN_NHHS), 2)
round(mean(interaction_FN_NHHS), 2)



round(mean(beta_TP_NHHS), 2)
round(mean(gamma_TP_NHHS), 2)
round(mean(varphi_TP_NHHS), 2)

round(mean(beta_TN_NHHS), 2)
round(mean(gamma_TN_NHHS), 2)
round(mean(varphi_TN_NHHS), 2)

round(mean(beta_FP_NHHS), 2)
round(mean(gamma_FP_NHHS), 2)
round(mean(varphi_FP_NHHS), 2)


round(mean(beta_FN_NHHS), 2)
round(mean(gamma_FN_NHHS), 2)
round(mean(varphi_FN_NHHS), 2)
######################

round(mean(alpha_est_mat_NHHS), 2)
round(mean(alpha_cover_per_NHHS), 2)

round(mean(beta_cover_per_NHHS), 2)
round(mean(beta_cover_nonzero_per_NHHS), 2)

round(mean(gamma_cover_per_NHHS), 2)
round(mean(gamma_cover_nonzero_per_NHHS), 2)

round(mean(varphi_cover_per_NHHS), 2)
round(mean(varphi_cover_nonzero_per_NHHS), 2)

round(mean(beta_RMSE_NHHS)*100, 2)
round(mean(beta_RMSE_nonzero_NHHS)*100, 2)
round(mean(beta_RMSE_zero_NHHS)*100, 2)

round(mean(gamma_RMSE_NHHS)*100, 2)
round(mean(gamma_RMSE_nonzero_NHHS)*100, 2)
round(mean(gamma_RMSE_zero_NHHS)*100, 2)

round(mean(varphi_RMSE_NHHS)*100, 2)
round(mean(varphi_RMSE_nonzero_NHHS)*100, 2)
round(mean(varphi_RMSE_zero_NHHS)*100, 2)

round(mean(delta_RMSE_NHHS)*100, 2)
round(mean(delta_RMSE_nonzero_NHHS)*100, 2)
round(mean(delta_RMSE_zero_NHHS)*100, 2)

round(mean(beta_sensitivity_NHHS)*100, 2)
round(mean(gamma_sensitivity_NHHS)*100, 2)
round(mean(varphi_sensitivity_NHHS)*100, 2)

round(mean(beta_specificity_NHHS)*100, 2)
round(mean(gamma_specificity_NHHS)*100, 2)
round(mean(varphi_specificity_NHHS)*100, 2)

round(mean(beta_precision_NHHS)*100, 2)
round(mean(gamma_precision_NHHS)*100, 2)
round(mean(varphi_precision_NHHS)*100, 2)



beta_true_array <- simplify2array(beta_true_list_NHHS)
beta_true_mean_NHHS <- apply(beta_true_array, c(1, 2), mean)
beta_est_array <- simplify2array(beta_est_list_NHHS)
beta_est_mean_NHHS <- apply(beta_est_array, c(1, 2), mean)
dim(beta_est_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_true_mean_NHHS[, k],
    beta_est_mean_NHHS[, k],
    xlab = "True beta",
    ylab = "Mean estimated beta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



gamma_true_array <- simplify2array(gamma_true_list_NHHS)
gamma_true_mean_NHHS <- apply(gamma_true_array, c(1, 2), mean)
gamma_est_array <- simplify2array(gamma_est_list_NHHS)
gamma_est_mean_NHHS <- apply(gamma_est_array, c(1, 2), mean)
dim(gamma_est_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_true_mean_NHHS[, k],
    gamma_est_mean_NHHS[, k],
    xlab = "True gamma",
    ylab = "Mean estimated gamma",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



delta_true_array <- simplify2array(delta_true_list_NHHS)
delta_true_mean_NHHS <- apply(delta_true_array, c(1, 2), mean)
delta_est_array <- simplify2array(delta_est_list_NHHS)
delta_est_mean_NHHS <- apply(delta_est_array, c(1, 2), mean)
dim(delta_est_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_true_mean_NHHS[, k],
    delta_est_mean_NHHS[, k],
    xlab = "True delta",
    ylab = "Mean estimated delta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



beta_bias_array <- simplify2array(beta_bias_NHHS)
beta_bias_mean_NHHS <- apply(beta_bias_array, c(1, 2), mean)
mean(beta_bias_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_mean_NHHS[, k],
    xlab = "Modifier",
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_nonzero_array <- simplify2array(beta_bias_nonzero_NHHS)
beta_bias_nonzero_mean_NHHS <- apply(beta_bias_nonzero_array, c(1, 2), mean)
mean(beta_bias_nonzero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_nonzero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_zero_array <- simplify2array(beta_bias_zero_NHHS)
beta_bias_zero_mean_NHHS <- apply(beta_bias_zero_array, c(1, 2), mean)
mean(beta_bias_zero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_zero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}





gamma_bias_array <- simplify2array(gamma_bias_NHHS)
gamma_bias_mean_NHHS <- apply(gamma_bias_array, c(1, 2), mean)
mean(gamma_bias_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_mean_NHHS[, k],
    xlab = "Exposure",
    main = paste("Outcome", k),
    pch = 19
  )
}

gamma_bias_nonzero_array <- simplify2array(gamma_bias_nonzero_NHHS)
gamma_bias_nonzero_mean_NHHS <- apply(gamma_bias_nonzero_array, c(1, 2), mean)
mean(gamma_bias_nonzero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_nonzero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_zero_array <- simplify2array(gamma_bias_zero_NHHS)
gamma_bias_zero_mean_NHHS <- apply(gamma_bias_zero_array, c(1, 2), mean)
mean(gamma_bias_zero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_zero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_array <- simplify2array(varphi_bias_NHHS)
varphi_bias_mean_NHHS <- apply(varphi_bias_array, c(1, 2), mean)
mean(varphi_bias_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_mean_NHHS[, k],
    xlab = "Covariate",
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_nonzero_array <- simplify2array(varphi_bias_nonzero_NHHS)
varphi_bias_nonzero_mean_NHHS <- apply(varphi_bias_nonzero_array, c(1, 2), mean)
mean(varphi_bias_nonzero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_nonzero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_zero_array <- simplify2array(varphi_bias_zero_NHHS)
varphi_bias_zero_mean_NHHS <- apply(varphi_bias_zero_array, c(1, 2), mean)
mean(varphi_bias_zero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_zero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_array <- simplify2array(delta_bias_NHHS)
delta_bias_mean_NHHS <- apply(delta_bias_array, c(1, 2), mean)
mean(delta_bias_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_mean_NHHS[, k],
    xlab = "Interaction",
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_nonzero_array <- simplify2array(delta_bias_nonzero_NHHS)
delta_bias_nonzero_mean_NHHS <- apply(delta_bias_nonzero_array, c(1, 2), mean)
mean(delta_bias_nonzero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_nonzero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_zero_array <- simplify2array(delta_bias_zero_NHHS)
delta_bias_zero_mean_NHHS <- apply(delta_bias_zero_array, c(1, 2), mean)
mean(delta_bias_zero_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_zero_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}







###################
# NHHSSI model

# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_0weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_10weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_25weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_50weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_100weak.rda")


# load("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_0weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_10weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_25weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_50weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_100weak.rda")


# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")



nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]
J <- dim(sim_data_list[[1]]$x)[2]
M <- dim(sim_data_list[[1]]$z)[2]

alpha_cover_per_NHHSSI <- c()
beta_cover_per_NHHSSI <- c()
gamma_cover_per_NHHSSI <- c()
varphi_cover_per_NHHSSI <- c()
delta_cover_per_NHHSSI <- c()

delta_RMSE_NHHSSI <- c()
RMSE_all_main_effects_NHHSSI <- c()

beta_cover_nonzero_per_NHHSSI <- c()
gamma_cover_nonzero_per_NHHSSI <- c()
varphi_cover_nonzero_per_NHHSSI <- c()
delta_cover_nonzero_per_NHHSSI <- c()

alpha_est_mat_NHHSSI <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_NHHSSI <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_NHHSSI <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_NHHSSI <- matrix(NA, nrow = nsample, ncol = K)

beta_est_list_NHHSSI <- list()
beta_true_list_NHHSSI <- list()
beta_lower_list_NHHSSI <- list()
beta_upper_list_NHHSSI <- list()

gamma_est_list_NHHSSI <- list()
gamma_true_list_NHHSSI <- list()
gamma_lower_list_NHHSSI <- list()
gamma_upper_list_NHHSSI <- list()

varphi_est_list_NHHSSI <- list()
varphi_true_list_NHHSSI <- list()
varphi_lower_list_NHHSSI <- list()
varphi_upper_list_NHHSSI <- list()

delta_est_list_NHHSSI <- list()
delta_true_list_NHHSSI <- list()
delta_lower_list_NHHSSI <- list()
delta_upper_list_NHHSSI <- list()

main_sensitivity_NHHSSI <- c()
main_specificity_NHHSSI <- c()
main_precision_NHHSSI <- c()

interaction_sensitivity_NHHSSI <- c()
interaction_specificity_NHHSSI <- c()
interaction_precision_NHHSSI <- c()

main_TP_NHHSSI <- c()
main_TN_NHHSSI <- c()
main_FP_NHHSSI <- c()
main_FN_NHHSSI <- c()

interaction_TP_NHHSSI <- c()
interaction_TN_NHHSSI <- c()
interaction_FP_NHHSSI <- c()
interaction_FN_NHHSSI <- c()


# ============================================================
# Additional diagnostics
# ============================================================

beta_bias_NHHSSI <- c()
gamma_bias_NHHSSI <- c()
varphi_bias_NHHSSI <- c()
delta_bias_NHHSSI <- c()

beta_bias_nonzero_NHHSSI <- c()
gamma_bias_nonzero_NHHSSI <- c()
varphi_bias_nonzero_NHHSSI <- c()
delta_bias_nonzero_NHHSSI <- c()

beta_bias_zero_NHHSSI <- c()
gamma_bias_zero_NHHSSI <- c()
varphi_bias_zero_NHHSSI <- c()
delta_bias_zero_NHHSSI <- c()

beta_RMSE_NHHSSI <- c()
gamma_RMSE_NHHSSI <- c()
varphi_RMSE_NHHSSI <- c()
delta_RMSE_NHHSSI <- c()

beta_RMSE_nonzero_NHHSSI <- c()
gamma_RMSE_nonzero_NHHSSI <- c()
varphi_RMSE_nonzero_NHHSSI <- c()
delta_RMSE_nonzero_NHHSSI <- c()

beta_RMSE_zero_NHHSSI <- c()
gamma_RMSE_zero_NHHSSI <- c()
varphi_RMSE_zero_NHHSSI <- c()
delta_RMSE_zero_NHHSSI <- c()

beta_TP_NHHSSI <- c()
gamma_TP_NHHSSI <- c()
varphi_TP_NHHSSI <- c()

beta_TN_NHHSSI <- c()
gamma_TN_NHHSSI <- c()
varphi_TN_NHHSSI <- c()

beta_FP_NHHSSI <- c()
gamma_FP_NHHSSI <- c()
varphi_FP_NHHSSI <- c()

beta_FN_NHHSSI <- c()
gamma_FN_NHHSSI <- c()
varphi_FN_NHHSSI <- c()

beta_sensitivity_NHHSSI <- c()
gamma_sensitivity_NHHSSI <- c()
varphi_sensitivity_NHHSSI <- c()

beta_specificity_NHHSSI <- c()
gamma_specificity_NHHSSI <- c()
varphi_specificity_NHHSSI <- c()

beta_precision_NHHSSI <- c()
gamma_precision_NHHSSI <- c()
varphi_precision_NHHSSI <- c()

delta_RMSE_indices_mat_NHHSSI <- matrix(NA, nrow = nsample, ncol = J*M)
delta_cover_indices_mat_NHHSSI <- matrix(NA, nrow = nsample, ncol = J*M)
delta_var_indices_mat_NHHSSI <- matrix(NA, nrow = nsample, ncol = J*M)

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  d <- sim_data_list[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_NHHSSI[sample, ] <- alpha_true
  alpha_est_mat_NHHSSI[sample, ] <- apply(theta_update_NHHSSI_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NHHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NHHSSI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NHHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NHHSSI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NHHSSI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_NHHSSI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NHHSSI_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NHHSSI_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NHHSSI_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_NHHSSI[[sample]] <- beta_est
  beta_lower_list_NHHSSI[[sample]] <- beta_lower
  beta_upper_list_NHHSSI[[sample]] <- beta_upper
  beta_cover_per_NHHSSI[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_NHHSSI[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # Additional diagnostics
  nonzero_beta <- beta_true != 0
  zero_beta <- beta_true == 0

  beta_bias <- beta_est - beta_true

  beta_bias_NHHSSI[[sample]] <- beta_bias
  beta_bias_nonzero_NHHSSI[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_zero_NHHSSI[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_NHHSSI[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_nonzero_NHHSSI[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_zero_NHHSSI[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_NHHSSI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NHHSSI_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NHHSSI_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NHHSSI_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_NHHSSI[[sample]] <- gamma_est
  gamma_lower_list_NHHSSI[[sample]] <- gamma_lower
  gamma_upper_list_NHHSSI[[sample]] <- gamma_upper
  gamma_cover_per_NHHSSI[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_NHHSSI[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  # Additional diagnostics
  nonzero_gamma <- gamma_true != 0
  zero_gamma <- gamma_true == 0

  gamma_bias <- gamma_est - gamma_true

  gamma_bias_NHHSSI[[sample]] <- gamma_bias
  gamma_bias_nonzero_NHHSSI[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_zero_NHHSSI[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_NHHSSI[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_nonzero_NHHSSI[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_zero_NHHSSI[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

  # varphi
  varphi_true <- sim_data_list[[sample]]$varphi_true
  varphi_true_list_NHHSSI[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list[[sample]][, o, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    varphi_cover[(o-1-J-M-J*M), ] <- (varphi_true[(o-1-J-M-J*M), ] >= varphi_lower[(o-1-J-M-J*M), ]) &
      (varphi_true[(o-1-J-M-J*M), ] <= varphi_upper[(o-1-J-M-J*M), ])
    varphi_MSE[o-1-J-M-J*M, ] <- (varphi_est[(o-1-J-M-J*M), ] - varphi_true[o-1-J-M-J*M, ])^2
  }
  varphi_est_list_NHHSSI[[sample]] <- varphi_est
  varphi_lower_list_NHHSSI[[sample]] <- varphi_lower
  varphi_upper_list_NHHSSI[[sample]] <- varphi_upper
  varphi_cover_per_NHHSSI[sample] <- (sum(varphi_cover)/(dim(varphi_true)[1]*dim(varphi_true)[2]))*100
  nonzero_idx_varphi <- varphi_true != 0
  varphi_cover_nonzero_per_NHHSSI[sample] <- (sum(varphi_cover[nonzero_idx_varphi]) / sum(nonzero_idx_varphi)) * 100

  # Additional diagnostics
  nonzero_varphi <- varphi_true != 0
  zero_varphi <- varphi_true == 0

  varphi_bias <- varphi_est - varphi_true

  varphi_bias_NHHSSI[[sample]] <- varphi_bias
  varphi_bias_nonzero_NHHSSI[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_zero_NHHSSI[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_NHHSSI[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_nonzero_NHHSSI[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_zero_NHHSSI[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


  main_MSE <- c(beta_MSE, gamma_MSE, varphi_MSE)
  RMSE_all_main_effects_NHHSSI[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NHHSSI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_NHHSSI[[sample]] <- delta_est
  delta_lower_list_NHHSSI[[sample]] <- delta_lower
  delta_upper_list_NHHSSI[[sample]] <- delta_upper
  delta_cover_per_NHHSSI[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_NHHSSI[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_NHHSSI[sample] <- sqrt(mean(delta_MSE))

  delta_cover_indices_mat_NHHSSI[sample, ] <- rowMeans(delta_cover)*100
  delta_RMSE_indices_mat_NHHSSI[sample, ] <- sqrt(rowMeans(delta_MSE))
  delta_var_indices_mat_NHHSSI[sample, ] <- rowMeans(delta_var)

  # Additional diagnostics
  nonzero_delta <- delta_true != 0
  zero_delta <- delta_true == 0

  delta_bias <- delta_est - delta_true

  delta_bias_NHHSSI[[sample]] <- delta_bias
  delta_bias_nonzero_NHHSSI[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_zero_NHHSSI[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_NHHSSI[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_nonzero_NHHSSI[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_zero_NHHSSI[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_varphi <- as.vector(varphi_true != 0)
  selected_varphi <- as.vector((varphi_lower > 0) | (varphi_upper < 0))

  truth_main <- c(truth_beta, truth_gamma, truth_varphi)
  selected_main <- c(selected_beta, selected_gamma, selected_varphi)

  main_TP_NHHSSI[sample] <- sum(truth_main == 1 & selected_main == 1)
  main_FN_NHHSSI[sample] <- sum(truth_main == 1 & selected_main == 0)
  main_TN_NHHSSI[sample] <- sum(truth_main == 0 & selected_main == 0)
  main_FP_NHHSSI[sample] <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_NHHSSI[sample] <- (main_TP_NHHSSI[sample] /
                       (main_TP_NHHSSI[sample] + main_FN_NHHSSI[sample]))
  main_specificity_NHHSSI[sample] <- (main_TN_NHHSSI[sample] /
                                        (main_TN_NHHSSI[sample] + main_FP_NHHSSI[sample]))
  main_precision_NHHSSI[sample] <- (main_TP_NHHSSI[sample] /
                                      (main_TP_NHHSSI[sample] + main_FP_NHHSSI[sample]))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP_NHHSSI[sample] <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN_NHHSSI[sample] <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN_NHHSSI[sample] <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP_NHHSSI[sample] <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NHHSSI[sample] <- (interaction_TP_NHHSSI[sample] /
                          (interaction_TP_NHHSSI[sample] + interaction_FN_NHHSSI[sample]))
  interaction_specificity_NHHSSI[sample] <- (interaction_TN_NHHSSI[sample] /
                      (interaction_TN_NHHSSI[sample] + interaction_FP_NHHSSI[sample]))
  interaction_precision_NHHSSI[sample] <- (interaction_TP_NHHSSI[sample] /
                  (interaction_TP_NHHSSI[sample] + interaction_FP_NHHSSI[sample]))

  # Additional diaagnostics
  beta_TP_NHHSSI[sample] <- sum(truth_beta == 1 & selected_beta == 1)
  beta_FN_NHHSSI[sample] <- sum(truth_beta == 1 & selected_beta == 0)
  beta_TN_NHHSSI[sample] <- sum(truth_beta == 0 & selected_beta == 0)
  beta_FP_NHHSSI[sample] <- sum(truth_beta == 0 & selected_beta == 1)

  gamma_TP_NHHSSI[sample] <- sum(truth_gamma == 1 & selected_gamma == 1)
  gamma_FN_NHHSSI[sample] <- sum(truth_gamma == 1 & selected_gamma == 0)
  gamma_TN_NHHSSI[sample] <- sum(truth_gamma == 0 & selected_gamma == 0)
  gamma_FP_NHHSSI[sample] <- sum(truth_gamma == 0 & selected_gamma == 1)

  varphi_TP_NHHSSI[sample] <- sum(truth_varphi == 1 & selected_varphi == 1)
  varphi_FN_NHHSSI[sample] <- sum(truth_varphi == 1 & selected_varphi == 0)
  varphi_TN_NHHSSI[sample] <- sum(truth_varphi == 0 & selected_varphi == 0)
  varphi_FP_NHHSSI[sample] <- sum(truth_varphi == 0 & selected_varphi == 1)

  beta_sensitivity_NHHSSI[sample] <- (beta_TP_NHHSSI[sample] /
                              (beta_TP_NHHSSI[sample] + beta_FN_NHHSSI[sample]))
  beta_specificity_NHHSSI[sample] <- (beta_TN_NHHSSI[sample] /
                            (beta_TN_NHHSSI[sample] + beta_FP_NHHSSI[sample]))
  beta_precision_NHHSSI[sample] <- (beta_TP_NHHSSI[sample] /
                              (beta_TP_NHHSSI[sample] + beta_FP_NHHSSI[sample]))

  gamma_sensitivity_NHHSSI[sample] <- (gamma_TP_NHHSSI[sample] /
                          (gamma_TP_NHHSSI[sample] + gamma_FN_NHHSSI[sample]))
  gamma_specificity_NHHSSI[sample] <- (gamma_TN_NHHSSI[sample] /
                        (gamma_TN_NHHSSI[sample] + gamma_FP_NHHSSI[sample]))
  gamma_precision_NHHSSI[sample] <- (gamma_TP_NHHSSI[sample] /
                            (gamma_TP_NHHSSI[sample] + gamma_FP_NHHSSI[sample]))

  varphi_sensitivity_NHHSSI[sample] <- (varphi_TP_NHHSSI[sample] /
                              (varphi_TP_NHHSSI[sample] + varphi_FN_NHHSSI[sample]))
  varphi_specificity_NHHSSI[sample] <- (varphi_TN_NHHSSI[sample] /
                              (varphi_TN_NHHSSI[sample] + varphi_FP_NHHSSI[sample]))
  varphi_precision_NHHSSI[sample] <- (varphi_TP_NHHSSI[sample] /
                              (varphi_TP_NHHSSI[sample] + varphi_FP_NHHSSI[sample]))
}




round(mean(RMSE_all_main_effects_NHHSSI)*100, 2)
round(mean(delta_RMSE_NHHSSI)*100, 2)

round((mean(beta_cover_per_NHHSSI) + mean(gamma_cover_per_NHHSSI) +
         mean(varphi_cover_per_NHHSSI))/3, 2)
round(mean(delta_cover_per_NHHSSI), 2)

round((mean(beta_cover_nonzero_per_NHHSSI) + mean(gamma_cover_nonzero_per_NHHSSI) +
         mean(varphi_cover_nonzero_per_NHHSSI))/3, 2)
round(mean(delta_cover_nonzero_per_NHHSSI), 2)

round(mean(main_sensitivity_NHHSSI) * 100, 2)
round(mean(interaction_sensitivity_NHHSSI) * 100, 2)

round(mean(main_specificity_NHHSSI) * 100, 2)
round(mean(interaction_specificity_NHHSSI) * 100, 2)

round(mean(main_precision_NHHSSI) * 100, 2)
round(mean(interaction_precision_NHHSSI) * 100, 2)


###############################

round(mean(main_TP_NHHSSI), 2)
round(mean(interaction_TP_NHHSSI), 2)

round(mean(main_TN_NHHSSI), 2)
round(mean(interaction_TN_NHHSSI), 2)

round(mean(main_FP_NHHSSI), 2)
round(mean(interaction_FP_NHHSSI), 2)

round(mean(main_FN_NHHSSI), 2)
round(mean(interaction_FN_NHHSSI), 2)



round(mean(beta_TP_NHHSSI), 2)
round(mean(gamma_TP_NHHSSI), 2)
round(mean(varphi_TP_NHHSSI), 2)

round(mean(beta_TN_NHHSSI), 2)
round(mean(gamma_TN_NHHSSI), 2)
round(mean(varphi_TN_NHHSSI), 2)

round(mean(beta_FP_NHHSSI), 2)
round(mean(gamma_FP_NHHSSI), 2)
round(mean(varphi_FP_NHHSSI), 2)


round(mean(beta_FN_NHHSSI), 2)
round(mean(gamma_FN_NHHSSI), 2)
round(mean(varphi_FN_NHHSSI), 2)
######################

round(mean(alpha_est_mat_NHHSSI), 2)
round(mean(alpha_cover_per_NHHSSI), 2)

round(mean(beta_cover_per_NHHSSI), 2)
round(mean(beta_cover_nonzero_per_NHHSSI), 2)

round(mean(gamma_cover_per_NHHSSI), 2)
round(mean(gamma_cover_nonzero_per_NHHSSI), 2)

round(mean(varphi_cover_per_NHHSSI), 2)
round(mean(varphi_cover_nonzero_per_NHHSSI), 2)

round(mean(beta_RMSE_NHHSSI)*100, 2)
round(mean(beta_RMSE_nonzero_NHHSSI)*100, 2)
round(mean(beta_RMSE_zero_NHHSSI)*100, 2)

round(mean(gamma_RMSE_NHHSSI)*100, 2)
round(mean(gamma_RMSE_nonzero_NHHSSI)*100, 2)
round(mean(gamma_RMSE_zero_NHHSSI)*100, 2)

round(mean(varphi_RMSE_NHHSSI)*100, 2)
round(mean(varphi_RMSE_nonzero_NHHSSI)*100, 2)
round(mean(varphi_RMSE_zero_NHHSSI)*100, 2)

round(mean(delta_RMSE_NHHSSI)*100, 2)
round(mean(delta_RMSE_nonzero_NHHSSI)*100, 2)
round(mean(delta_RMSE_zero_NHHSSI)*100, 2)

round(mean(beta_sensitivity_NHHSSI)*100, 2)
round(mean(gamma_sensitivity_NHHSSI)*100, 2)
round(mean(varphi_sensitivity_NHHSSI)*100, 2)

round(mean(beta_specificity_NHHSSI)*100, 2)
round(mean(gamma_specificity_NHHSSI)*100, 2)
round(mean(varphi_specificity_NHHSSI)*100, 2)

round(mean(beta_precision_NHHSSI)*100, 2)
round(mean(gamma_precision_NHHSSI)*100, 2)
round(mean(varphi_precision_NHHSSI)*100, 2)



beta_true_array <- simplify2array(beta_true_list_NHHSSI)
beta_true_mean_NHHSSI <- apply(beta_true_array, c(1, 2), mean)
beta_est_array <- simplify2array(beta_est_list_NHHSSI)
beta_est_mean_NHHSSI <- apply(beta_est_array, c(1, 2), mean)
dim(beta_est_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_true_mean_NHHSSI[, k],
    beta_est_mean_NHHSSI[, k],
    xlab = "True beta",
    ylab = "Mean estimated beta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



gamma_true_array <- simplify2array(gamma_true_list_NHHSSI)
gamma_true_mean_NHHSSI <- apply(gamma_true_array, c(1, 2), mean)
gamma_est_array <- simplify2array(gamma_est_list_NHHSSI)
gamma_est_mean_NHHSSI <- apply(gamma_est_array, c(1, 2), mean)
dim(gamma_est_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_true_mean_NHHSSI[, k],
    gamma_est_mean_NHHSSI[, k],
    xlab = "True gamma",
    ylab = "Mean estimated gamma",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



delta_true_array <- simplify2array(delta_true_list_NHHSSI)
delta_true_mean_NHHSSI <- apply(delta_true_array, c(1, 2), mean)
delta_est_array <- simplify2array(delta_est_list_NHHSSI)
delta_est_mean_NHHSSI <- apply(delta_est_array, c(1, 2), mean)
dim(delta_est_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_true_mean_NHHSSI[, k],
    delta_est_mean_NHHSSI[, k],
    xlab = "True delta",
    ylab = "Mean estimated delta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



beta_bias_array <- simplify2array(beta_bias_NHHSSI)
beta_bias_mean_NHHSSI <- apply(beta_bias_array, c(1, 2), mean)
mean(beta_bias_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_mean_NHHSSI[, k],
    xlab = "Modifier",
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_nonzero_array <- simplify2array(beta_bias_nonzero_NHHSSI)
beta_bias_nonzero_mean_NHHSSI <- apply(beta_bias_nonzero_array, c(1, 2), mean)
mean(beta_bias_nonzero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_nonzero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_zero_array <- simplify2array(beta_bias_zero_NHHSSI)
beta_bias_zero_mean_NHHSSI <- apply(beta_bias_zero_array, c(1, 2), mean)
mean(beta_bias_zero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_zero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}





gamma_bias_array <- simplify2array(gamma_bias_NHHSSI)
gamma_bias_mean_NHHSSI <- apply(gamma_bias_array, c(1, 2), mean)
mean(gamma_bias_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_mean_NHHSSI[, k],
    xlab = "Exposure",
    main = paste("Outcome", k),
    pch = 19
  )
}

gamma_bias_nonzero_array <- simplify2array(gamma_bias_nonzero_NHHSSI)
gamma_bias_nonzero_mean_NHHSSI <- apply(gamma_bias_nonzero_array, c(1, 2), mean)
mean(gamma_bias_nonzero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_nonzero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_zero_array <- simplify2array(gamma_bias_zero_NHHSSI)
gamma_bias_zero_mean_NHHSSI <- apply(gamma_bias_zero_array, c(1, 2), mean)
mean(gamma_bias_zero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_zero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_array <- simplify2array(varphi_bias_NHHSSI)
varphi_bias_mean_NHHSSI <- apply(varphi_bias_array, c(1, 2), mean)
mean(varphi_bias_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_mean_NHHSSI[, k],
    xlab = "Covariate",
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_nonzero_array <- simplify2array(varphi_bias_nonzero_NHHSSI)
varphi_bias_nonzero_mean_NHHSSI <- apply(varphi_bias_nonzero_array, c(1, 2), mean)
mean(varphi_bias_nonzero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_nonzero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_zero_array <- simplify2array(varphi_bias_zero_NHHSSI)
varphi_bias_zero_mean_NHHSSI <- apply(varphi_bias_zero_array, c(1, 2), mean)
mean(varphi_bias_zero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_zero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_array <- simplify2array(delta_bias_NHHSSI)
delta_bias_mean_NHHSSI <- apply(delta_bias_array, c(1, 2), mean)
mean(delta_bias_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_mean_NHHSSI[, k],
    xlab = "Interaction",
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_nonzero_array <- simplify2array(delta_bias_nonzero_NHHSSI)
delta_bias_nonzero_mean_NHHSSI <- apply(delta_bias_nonzero_array, c(1, 2), mean)
mean(delta_bias_nonzero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_nonzero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_zero_array <- simplify2array(delta_bias_zero_NHHSSI)
delta_bias_zero_mean_NHHSSI <- apply(delta_bias_zero_array, c(1, 2), mean)
mean(delta_bias_zero_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_zero_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}







############################################################################
# NP model

# load("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_0weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_10weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_25weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_50weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_100weak.rda")

# load("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_0weak.rda")
# load("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_10weak.rda")
# load("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_25weak.rda")
# load("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_50weak.rda")
# load("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_100weak.rda")



# load("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")
# load("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")


nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]
J <- dim(sim_data_list[[1]]$x)[2]
M <- dim(sim_data_list[[1]]$z)[2]

alpha_cover_per_NP <- c()
beta_cover_per_NP <- c()
gamma_cover_per_NP <- c()
varphi_cover_per_NP <- c()
delta_cover_per_NP <- c()

delta_RMSE_NP <- c()
RMSE_all_main_effects_NP <- c()

beta_cover_nonzero_per_NP <- c()
gamma_cover_nonzero_per_NP <- c()
varphi_cover_nonzero_per_NP <- c()
delta_cover_nonzero_per_NP <- c()

alpha_est_mat_NP <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_NP <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_NP <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_NP <- matrix(NA, nrow = nsample, ncol = K)

beta_est_list_NP <- list()
beta_true_list_NP <- list()
beta_lower_list_NP <- list()
beta_upper_list_NP <- list()

gamma_est_list_NP <- list()
gamma_true_list_NP <- list()
gamma_lower_list_NP <- list()
gamma_upper_list_NP <- list()

varphi_est_list_NP <- list()
varphi_true_list_NP <- list()
varphi_lower_list_NP <- list()
varphi_upper_list_NP <- list()

delta_est_list_NP <- list()
delta_true_list_NP <- list()
delta_lower_list_NP <- list()
delta_upper_list_NP <- list()

main_sensitivity_NP <- c()
main_specificity_NP <- c()
main_precision_NP <- c()

interaction_sensitivity_NP <- c()
interaction_specificity_NP <- c()
interaction_precision_NP <- c()

main_TP_NP <- c()
main_TN_NP <- c()
main_FP_NP <- c()
main_FN_NP <- c()

interaction_TP_NP <- c()
interaction_TN_NP <- c()
interaction_FP_NP <- c()
interaction_FN_NP <- c()


# ============================================================
# Additional diagnostics
# ============================================================

beta_bias_NP <- c()
gamma_bias_NP <- c()
varphi_bias_NP <- c()
delta_bias_NP <- c()

beta_bias_nonzero_NP <- c()
gamma_bias_nonzero_NP <- c()
varphi_bias_nonzero_NP <- c()
delta_bias_nonzero_NP <- c()

beta_bias_zero_NP <- c()
gamma_bias_zero_NP <- c()
varphi_bias_zero_NP <- c()
delta_bias_zero_NP <- c()

beta_RMSE_NP <- c()
gamma_RMSE_NP <- c()
varphi_RMSE_NP <- c()
delta_RMSE_NP <- c()

beta_RMSE_nonzero_NP <- c()
gamma_RMSE_nonzero_NP <- c()
varphi_RMSE_nonzero_NP <- c()
delta_RMSE_nonzero_NP <- c()

beta_RMSE_zero_NP <- c()
gamma_RMSE_zero_NP <- c()
varphi_RMSE_zero_NP <- c()
delta_RMSE_zero_NP <- c()

beta_TP_NP <- c()
gamma_TP_NP <- c()
varphi_TP_NP <- c()

beta_TN_NP <- c()
gamma_TN_NP <- c()
varphi_TN_NP <- c()

beta_FP_NP <- c()
gamma_FP_NP <- c()
varphi_FP_NP <- c()

beta_FN_NP <- c()
gamma_FN_NP <- c()
varphi_FN_NP <- c()

beta_sensitivity_NP <- c()
gamma_sensitivity_NP <- c()
varphi_sensitivity_NP <- c()

beta_specificity_NP <- c()
gamma_specificity_NP <- c()
varphi_specificity_NP <- c()

beta_precision_NP <- c()
gamma_precision_NP <- c()
varphi_precision_NP <- c()

delta_RMSE_indices_mat_NP <- matrix(NA, nrow = nsample, ncol = J*M)
delta_cover_indices_mat_NP <- matrix(NA, nrow = nsample, ncol = J*M)
delta_var_indices_mat_NP <- matrix(NA, nrow = nsample, ncol = J*M)

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  d <- sim_data_list[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_NP[sample, ] <- alpha_true
  alpha_est_mat_NP[sample, ] <- apply(theta_update_NP_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NP_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NP[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NP_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NP[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NP[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_NP[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NP_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NP_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NP_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_NP[[sample]] <- beta_est
  beta_lower_list_NP[[sample]] <- beta_lower
  beta_upper_list_NP[[sample]] <- beta_upper
  beta_cover_per_NP[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_NP[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # Additional diagnostics
  nonzero_beta <- beta_true != 0
  zero_beta <- beta_true == 0

  beta_bias <- beta_est - beta_true

  beta_bias_NP[[sample]] <- beta_bias
  beta_bias_nonzero_NP[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_zero_NP[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_NP[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_nonzero_NP[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_zero_NP[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_NP[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NP_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NP_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NP_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_NP[[sample]] <- gamma_est
  gamma_lower_list_NP[[sample]] <- gamma_lower
  gamma_upper_list_NP[[sample]] <- gamma_upper
  gamma_cover_per_NP[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_NP[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  # Additional diagnostics
  nonzero_gamma <- gamma_true != 0
  zero_gamma <- gamma_true == 0

  gamma_bias <- gamma_est - gamma_true

  gamma_bias_NP[[sample]] <- gamma_bias
  gamma_bias_nonzero_NP[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_zero_NP[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_NP[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_nonzero_NP[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_zero_NP[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

  # varphi
  varphi_true <- sim_data_list[[sample]]$varphi_true
  varphi_true_list_NP[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_NP_list[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_NP_list[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_NP_list[[sample]][, o, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    varphi_cover[(o-1-J-M-J*M), ] <- (varphi_true[(o-1-J-M-J*M), ] >= varphi_lower[(o-1-J-M-J*M), ]) &
      (varphi_true[(o-1-J-M-J*M), ] <= varphi_upper[(o-1-J-M-J*M), ])
    varphi_MSE[o-1-J-M-J*M, ] <- (varphi_est[(o-1-J-M-J*M), ] - varphi_true[o-1-J-M-J*M, ])^2
  }
  varphi_est_list_NP[[sample]] <- varphi_est
  varphi_lower_list_NP[[sample]] <- varphi_lower
  varphi_upper_list_NP[[sample]] <- varphi_upper
  varphi_cover_per_NP[sample] <- (sum(varphi_cover)/(dim(varphi_true)[1]*dim(varphi_true)[2]))*100
  nonzero_idx_varphi <- varphi_true != 0
  varphi_cover_nonzero_per_NP[sample] <- (sum(varphi_cover[nonzero_idx_varphi]) / sum(nonzero_idx_varphi)) * 100

  # Additional diagnostics
  nonzero_varphi <- varphi_true != 0
  zero_varphi <- varphi_true == 0

  varphi_bias <- varphi_est - varphi_true

  varphi_bias_NP[[sample]] <- varphi_bias
  varphi_bias_nonzero_NP[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_zero_NP[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_NP[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_nonzero_NP[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_zero_NP[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


  main_MSE <- c(beta_MSE, gamma_MSE, varphi_MSE)
  RMSE_all_main_effects_NP[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NP[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NP_list[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_NP_list[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NP_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NP_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_NP[[sample]] <- delta_est
  delta_lower_list_NP[[sample]] <- delta_lower
  delta_upper_list_NP[[sample]] <- delta_upper
  delta_cover_per_NP[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_NP[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_NP[sample] <- sqrt(mean(delta_MSE))

  delta_cover_indices_mat_NP[sample, ] <- rowMeans(delta_cover)*100
  delta_RMSE_indices_mat_NP[sample, ] <- sqrt(rowMeans(delta_MSE))
  delta_var_indices_mat_NP[sample, ] <- rowMeans(delta_var)

  # Additional diagnostics
  nonzero_delta <- delta_true != 0
  zero_delta <- delta_true == 0

  delta_bias <- delta_est - delta_true

  delta_bias_NP[[sample]] <- delta_bias
  delta_bias_nonzero_NP[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_zero_NP[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_NP[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_nonzero_NP[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_zero_NP[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_varphi <- as.vector(varphi_true != 0)
  selected_varphi <- as.vector((varphi_lower > 0) | (varphi_upper < 0))

  truth_main <- c(truth_beta, truth_gamma, truth_varphi)
  selected_main <- c(selected_beta, selected_gamma, selected_varphi)

  main_TP_NP[sample] <- sum(truth_main == 1 & selected_main == 1)
  main_FN_NP[sample] <- sum(truth_main == 1 & selected_main == 0)
  main_TN_NP[sample] <- sum(truth_main == 0 & selected_main == 0)
  main_FP_NP[sample] <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_NP[sample] <- (main_TP_NP[sample] /
                                        (main_TP_NP[sample] + main_FN_NP[sample]))
  main_specificity_NP[sample] <- (main_TN_NP[sample] /
                                        (main_TN_NP[sample] + main_FP_NP[sample]))
  main_precision_NP[sample] <- (main_TP_NP[sample] /
                                      (main_TP_NP[sample] + main_FP_NP[sample]))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP_NP[sample] <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN_NP[sample] <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN_NP[sample] <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP_NP[sample] <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NP[sample] <- (interaction_TP_NP[sample] /
                                               (interaction_TP_NP[sample] + interaction_FN_NP[sample]))
  interaction_specificity_NP[sample] <- (interaction_TN_NP[sample] /
                                               (interaction_TN_NP[sample] + interaction_FP_NP[sample]))
  interaction_precision_NP[sample] <- (interaction_TP_NP[sample] /
                                             (interaction_TP_NP[sample] + interaction_FP_NP[sample]))

  # Additional diaagnostics
  beta_TP_NP[sample] <- sum(truth_beta == 1 & selected_beta == 1)
  beta_FN_NP[sample] <- sum(truth_beta == 1 & selected_beta == 0)
  beta_TN_NP[sample] <- sum(truth_beta == 0 & selected_beta == 0)
  beta_FP_NP[sample] <- sum(truth_beta == 0 & selected_beta == 1)

  gamma_TP_NP[sample] <- sum(truth_gamma == 1 & selected_gamma == 1)
  gamma_FN_NP[sample] <- sum(truth_gamma == 1 & selected_gamma == 0)
  gamma_TN_NP[sample] <- sum(truth_gamma == 0 & selected_gamma == 0)
  gamma_FP_NP[sample] <- sum(truth_gamma == 0 & selected_gamma == 1)

  varphi_TP_NP[sample] <- sum(truth_varphi == 1 & selected_varphi == 1)
  varphi_FN_NP[sample] <- sum(truth_varphi == 1 & selected_varphi == 0)
  varphi_TN_NP[sample] <- sum(truth_varphi == 0 & selected_varphi == 0)
  varphi_FP_NP[sample] <- sum(truth_varphi == 0 & selected_varphi == 1)

  beta_sensitivity_NP[sample] <- (beta_TP_NP[sample] /
                                        (beta_TP_NP[sample] + beta_FN_NP[sample]))
  beta_specificity_NP[sample] <- (beta_TN_NP[sample] /
                                        (beta_TN_NP[sample] + beta_FP_NP[sample]))
  beta_precision_NP[sample] <- (beta_TP_NP[sample] /
                                      (beta_TP_NP[sample] + beta_FP_NP[sample]))

  gamma_sensitivity_NP[sample] <- (gamma_TP_NP[sample] /
                                         (gamma_TP_NP[sample] + gamma_FN_NP[sample]))
  gamma_specificity_NP[sample] <- (gamma_TN_NP[sample] /
                                         (gamma_TN_NP[sample] + gamma_FP_NP[sample]))
  gamma_precision_NP[sample] <- (gamma_TP_NP[sample] /
                                       (gamma_TP_NP[sample] + gamma_FP_NP[sample]))

  varphi_sensitivity_NP[sample] <- (varphi_TP_NP[sample] /
                                          (varphi_TP_NP[sample] + varphi_FN_NP[sample]))
  varphi_specificity_NP[sample] <- (varphi_TN_NP[sample] /
                                          (varphi_TN_NP[sample] + varphi_FP_NP[sample]))
  varphi_precision_NP[sample] <- (varphi_TP_NP[sample] /
                                        (varphi_TP_NP[sample] + varphi_FP_NP[sample]))
}




round(mean(RMSE_all_main_effects_NP)*100, 2)
round(mean(delta_RMSE_NP)*100, 2)

round((mean(beta_cover_per_NP) + mean(gamma_cover_per_NP) +
         mean(varphi_cover_per_NP))/3, 2)
round(mean(delta_cover_per_NP), 2)

round((mean(beta_cover_nonzero_per_NP) + mean(gamma_cover_nonzero_per_NP) +
         mean(varphi_cover_nonzero_per_NP))/3, 2)
round(mean(delta_cover_nonzero_per_NP), 2)

round(mean(main_sensitivity_NP) * 100, 2)
round(mean(interaction_sensitivity_NP) * 100, 2)

round(mean(main_specificity_NP) * 100, 2)
round(mean(interaction_specificity_NP) * 100, 2)

round(mean(main_precision_NP) * 100, 2)
round(mean(interaction_precision_NP) * 100, 2)


###############################

round(mean(main_TP_NP), 2)
round(mean(interaction_TP_NP), 2)

round(mean(main_TN_NP), 2)
round(mean(interaction_TN_NP), 2)

round(mean(main_FP_NP), 2)
round(mean(interaction_FP_NP), 2)

round(mean(main_FN_NP), 2)
round(mean(interaction_FN_NP), 2)



round(mean(beta_TP_NP), 2)
round(mean(gamma_TP_NP), 2)
round(mean(varphi_TP_NP), 2)

round(mean(beta_TN_NP), 2)
round(mean(gamma_TN_NP), 2)
round(mean(varphi_TN_NP), 2)

round(mean(beta_FP_NP), 2)
round(mean(gamma_FP_NP), 2)
round(mean(varphi_FP_NP), 2)


round(mean(beta_FN_NP), 2)
round(mean(gamma_FN_NP), 2)
round(mean(varphi_FN_NP), 2)
######################

round(mean(alpha_est_mat_NP), 2)
round(mean(alpha_cover_per_NP), 2)

round(mean(beta_cover_per_NP), 2)
round(mean(beta_cover_nonzero_per_NP), 2)

round(mean(gamma_cover_per_NP), 2)
round(mean(gamma_cover_nonzero_per_NP), 2)

round(mean(varphi_cover_per_NP), 2)
round(mean(varphi_cover_nonzero_per_NP), 2)

round(mean(beta_RMSE_NP)*100, 2)
round(mean(beta_RMSE_nonzero_NP)*100, 2)
round(mean(beta_RMSE_zero_NP)*100, 2)

round(mean(gamma_RMSE_NP)*100, 2)
round(mean(gamma_RMSE_nonzero_NP)*100, 2)
round(mean(gamma_RMSE_zero_NP)*100, 2)

round(mean(varphi_RMSE_NP)*100, 2)
round(mean(varphi_RMSE_nonzero_NP)*100, 2)
round(mean(varphi_RMSE_zero_NP)*100, 2)

round(mean(delta_RMSE_NP)*100, 2)
round(mean(delta_RMSE_nonzero_NP)*100, 2)
round(mean(delta_RMSE_zero_NP)*100, 2)

round(mean(beta_sensitivity_NP)*100, 2)
round(mean(gamma_sensitivity_NP)*100, 2)
round(mean(varphi_sensitivity_NP)*100, 2)

round(mean(beta_specificity_NP)*100, 2)
round(mean(gamma_specificity_NP)*100, 2)
round(mean(varphi_specificity_NP)*100, 2)

round(mean(beta_precision_NP)*100, 2)
round(mean(gamma_precision_NP)*100, 2)
round(mean(varphi_precision_NP)*100, 2)



beta_true_array <- simplify2array(beta_true_list_NP)
beta_true_mean_NP <- apply(beta_true_array, c(1, 2), mean)
beta_est_array <- simplify2array(beta_est_list_NP)
beta_est_mean_NP <- apply(beta_est_array, c(1, 2), mean)
dim(beta_est_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_true_mean_NP[, k],
    beta_est_mean_NP[, k],
    xlab = "True beta",
    ylab = "Mean estimated beta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



gamma_true_array <- simplify2array(gamma_true_list_NP)
gamma_true_mean_NP <- apply(gamma_true_array, c(1, 2), mean)
gamma_est_array <- simplify2array(gamma_est_list_NP)
gamma_est_mean_NP <- apply(gamma_est_array, c(1, 2), mean)
dim(gamma_est_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_true_mean_NP[, k],
    gamma_est_mean_NP[, k],
    xlab = "True gamma",
    ylab = "Mean estimated gamma",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



delta_true_array <- simplify2array(delta_true_list_NP)
delta_true_mean_NP <- apply(delta_true_array, c(1, 2), mean)
delta_est_array <- simplify2array(delta_est_list_NP)
delta_est_mean_NP <- apply(delta_est_array, c(1, 2), mean)
dim(delta_est_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_true_mean_NP[, k],
    delta_est_mean_NP[, k],
    xlab = "True delta",
    ylab = "Mean estimated delta",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}



beta_bias_array <- simplify2array(beta_bias_NP)
beta_bias_mean_NP <- apply(beta_bias_array, c(1, 2), mean)
mean(beta_bias_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_mean_NP[, k],
    xlab = "Modifier",
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_nonzero_array <- simplify2array(beta_bias_nonzero_NP)
beta_bias_nonzero_mean_NP <- apply(beta_bias_nonzero_array, c(1, 2), mean)
mean(beta_bias_nonzero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_nonzero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_zero_array <- simplify2array(beta_bias_zero_NP)
beta_bias_zero_mean_NP <- apply(beta_bias_zero_array, c(1, 2), mean)
mean(beta_bias_zero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_zero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}





gamma_bias_array <- simplify2array(gamma_bias_NP)
gamma_bias_mean_NP <- apply(gamma_bias_array, c(1, 2), mean)
mean(gamma_bias_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_mean_NP[, k],
    xlab = "Exposure",
    main = paste("Outcome", k),
    pch = 19
  )
}

gamma_bias_nonzero_array <- simplify2array(gamma_bias_nonzero_NP)
gamma_bias_nonzero_mean_NP <- apply(gamma_bias_nonzero_array, c(1, 2), mean)
mean(gamma_bias_nonzero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_nonzero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_zero_array <- simplify2array(gamma_bias_zero_NP)
gamma_bias_zero_mean_NP <- apply(gamma_bias_zero_array, c(1, 2), mean)
mean(gamma_bias_zero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_zero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_array <- simplify2array(varphi_bias_NP)
varphi_bias_mean_NP <- apply(varphi_bias_array, c(1, 2), mean)
mean(varphi_bias_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_mean_NP[, k],
    xlab = "Covariate",
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_nonzero_array <- simplify2array(varphi_bias_nonzero_NP)
varphi_bias_nonzero_mean_NP <- apply(varphi_bias_nonzero_array, c(1, 2), mean)
mean(varphi_bias_nonzero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_nonzero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_zero_array <- simplify2array(varphi_bias_zero_NP)
varphi_bias_zero_mean_NP <- apply(varphi_bias_zero_array, c(1, 2), mean)
mean(varphi_bias_zero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_zero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_array <- simplify2array(delta_bias_NP)
delta_bias_mean_NP <- apply(delta_bias_array, c(1, 2), mean)
mean(delta_bias_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_mean_NP[, k],
    xlab = "Interaction",
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_nonzero_array <- simplify2array(delta_bias_nonzero_NP)
delta_bias_nonzero_mean_NP <- apply(delta_bias_nonzero_array, c(1, 2), mean)
mean(delta_bias_nonzero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_nonzero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_zero_array <- simplify2array(delta_bias_zero_NP)
delta_bias_zero_mean_NP <- apply(delta_bias_zero_array, c(1, 2), mean)
mean(delta_bias_zero_mean_NP)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_zero_mean_NP[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}




###################################################################
library(ggplot2)


# rmse heatmap with simulated predictors
rmse_hhssi <- colMeans(delta_RMSE_indices_mat_HHSSI)
rmse_nhhssi <- colMeans(delta_RMSE_indices_mat_NHHSSI)
rmse_hhs <- colMeans(delta_RMSE_indices_mat_HHS)
rmse_nhhs <- colMeans(delta_RMSE_indices_mat_NHHS)
rmse_np <- colMeans(delta_RMSE_indices_mat_NP)

df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$RMSE_HHS <- rmse_hhs
df$RMSE_HHSSI <- rmse_hhssi
df$RMSE_NHHS <- rmse_nhhs
df$RMSE_NHHSSI <- rmse_nhhssi
df$RMSE_NP <- rmse_np
rmse_range <- range(c(rmse_hhs, rmse_hhssi,
                      rmse_nhhs, rmse_nhhssi,
                      rmse_np))
mid_value <- mean(rmse_range)
df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_HHS,
                            Method = "HHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_NHHS,
                            Method = "NHHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_NP,
                            Method = "NP"))

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.2) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = rmse_range,
    name = "RMSE"
  ) +

  labs(
    title = "RMSE Heatmap of Modifier–Exposure Interaction Estimates",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))











rmse_range <- range(c(rmse_hhs, rmse_hhssi,
                      rmse_nhhs, rmse_nhhssi))
mid_value <- mean(rmse_range)
df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_HHS,
                            Method = "HHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_NHHS,
                            Method = "NHHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"))

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.2) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = rmse_range,
    name = "RMSE"
  ) +

  labs(
    title = "RMSE Heatmap of Modifier–Exposure Interaction Estimates",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9)
  )


###################################################################

# coverage heatmap with simulated predictors
cover_hhs <- colMeans(delta_cover_indices_mat_HHS)
cover_hhssi <- colMeans(delta_cover_indices_mat_HHSSI)
cover_nhhs <- colMeans(delta_cover_indices_mat_NHHS)
cover_nhhssi <- colMeans(delta_cover_indices_mat_NHHSSI)
cover_np <- colMeans(delta_cover_indices_mat_NP)

df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$cover_HHS <- cover_hhs
df$cover_HHSSI <- cover_hhssi
df$cover_NHHS <- cover_nhhs
df$cover_NHHSSI <- cover_nhhssi
df$cover_NP <- cover_np
cover_range <- range(c(cover_hhs,
                       cover_hhssi,
                       cover_nhhs,
                       cover_nhhssi,
                       cover_np))
mid_value <- mean(cover_range)
df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_HHS,
                            Method = "HHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_NHHS,
                            Method = "NHHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_NHHSSI,
                            Method = "NHHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_NP,
                            Method = "NP"))

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

ggplot(df_long, aes(x = exposure, y = modifier, fill = cover)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.2) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = cover_range,
    name = "cover"
  ) +

  labs(
    title = "Coverage Heatmap of Modifier–Exposure Interaction Estimates",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))





cover_range <- range(c(cover_hhs,
                       cover_hhssi,
                       cover_nhhs,
                       cover_nhhssi))
mid_value <- mean(cover_range)
df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_HHS,
                            Method = "HHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_NHHS,
                            Method = "NHHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_NHHSSI,
                            Method = "NHHSSI"))

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

ggplot(df_long, aes(x = exposure, y = modifier, fill = cover)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.2) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = cover_range,
    name = "cover"
  ) +

  labs(
    title = "Coverage Heatmap of Modifier–Exposure Interaction Estimates",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9))

########################################################################

# variance heatmap with simulated predictors
var_hhs <- colMeans(delta_var_indices_mat_HHS)
var_hhssi <- colMeans(delta_var_indices_mat_HHSSI)
var_nhhs <- colMeans(delta_var_indices_mat_NHHS)
var_nhhssi <- colMeans(delta_var_indices_mat_NHHSSI)
var_np <- colMeans(delta_var_indices_mat_NP)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$var_HHS <- var_hhs
df$var_HHSSI <- var_hhssi
df$var_NHHS <- var_nhhs
df$var_NHHSSI <- var_nhhssi
df$var_NP <- var_np
var_range <- range(c(var_hhs,
                     var_hhssi,
                     var_nhhs,
                     var_nhhssi,
                     var_np))
mid_value <- mean(var_range)
df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_HHS,
                            Method = "HHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_NHHS,
                            Method = "NHHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_NHHSSI,
                            Method = "NHHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_NP,
                            Method = "NP"))

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

ggplot(df_long, aes(x = exposure, y = modifier, fill = var)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.2) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = var_range,
    name = "var"
  ) +

  labs(
    title = "Variance Heatmap of Modifier–Exposure Interaction Estimates",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))









var_range <- range(c(var_hhs,
                     var_hhssi,
                     var_nhhs,
                     var_nhhssi))
mid_value <- mean(var_range)
df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_HHS,
                            Method = "HHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_NHHS,
                            Method = "NHHS"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_NHHSSI,
                            Method = "NHHSSI"))

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

ggplot(df_long, aes(x = exposure, y = modifier, fill = var)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.2) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = var_range,
    name = "var"
  ) +

  labs(
    title = "Variance Heatmap of Modifier–Exposure Interaction Estimates",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

####################################################################################

# Convergence of NHHSSI for weak heredity
sample <- 1

alpha_post_mean <- apply(theta_update_NHHSSI_list[[sample]][, 1, ], 2, mean)  # 1 × K
alpha_true <- sim_data_list[[sample]]$alpha_true
alpha_lower <- apply(theta_update_NHHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
alpha_upper <- apply(theta_update_NHHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)

beta_post_mean <- apply(theta_update_NHHSSI_list[[sample]][, 2:(1 + J), ], c(2, 3), mean)  # J × K
beta_lower <- apply(theta_update_NHHSSI_list[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.025)  # J × K
beta_upper <- apply(theta_update_NHHSSI_list[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.975)  # J × K
beta_true <- sim_data_list[[sample]]$beta_true

gamma_post_mean <- apply(theta_update_NHHSSI_list[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), mean)  # M × K
gamma_lower <- apply(theta_update_NHHSSI_list[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.025)  # M × K
gamma_upper <- apply(theta_update_NHHSSI_list[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.975)  # M × K
gamma_true <- sim_data_list[[sample]]$gamma_true

delta_post_mean <- apply(theta_update_NHHSSI_list[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), mean)  # JM × K
delta_lower <- apply(theta_update_NHHSSI_list[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.025)  # JM × K
delta_upper <- apply(theta_update_NHHSSI_list[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.975)  # JM × K
delta_true <- sim_data_list[[sample]]$delta_true


alpha_df <- data.frame(
  coef = paste0("A", "_k", rep(1:K)),
  mean = as.vector(alpha_post_mean),
  lower = as.vector(alpha_lower),
  upper = as.vector(alpha_upper),
  true  = as.vector(alpha_true))
ggplot(alpha_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(y = "Alpha coefficient", x = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of alpha")


beta_df <- data.frame(
  coef = paste0("B", rep(1:J, times = K), "_k", rep(1:K, each = J)),
  mean = as.vector(beta_post_mean),
  lower = as.vector(beta_lower),
  upper = as.vector(beta_upper),
  true  = as.vector(beta_true))
ggplot(beta_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(y = "Beta coefficient", x = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of beta")


gamma_df <- data.frame(
  coef = paste0("g", rep(1:M, times = K), "_k", rep(1:K, each = M)),
  mean = as.vector(gamma_post_mean),
  lower = as.vector(gamma_lower),
  upper = as.vector(gamma_upper),
  true  = as.vector(gamma_true))
ggplot(gamma_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(x = "gamma coefficient", y = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of gamma")


delta_df <- data.frame(
  coef = paste0("d", rep(1:(J*M), times = K), "_k", rep(1:K, each = (J*M))),
  mean = as.vector(delta_post_mean),
  lower = as.vector(delta_lower),
  upper = as.vector(delta_upper),
  true  = as.vector(delta_true))
ggplot(delta_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(x = "delta coefficient", y = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of delta")
#####################################################################################


####################################################################################

# Convergence of HHSSI for weak heredity

sample <- 1

alpha_post_mean <- apply(theta_update_HHSSI_list[[sample]][, 1, ], 2, mean)  # 1 × K
alpha_true <- sim_data_list[[sample]]$alpha_true
alpha_lower <- apply(theta_update_HHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
alpha_upper <- apply(theta_update_HHSSI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)

beta_post_mean <- apply(theta_update_HHSSI_list[[sample]][, 2:(1 + J), ], c(2, 3), mean)  # J × K
beta_lower <- apply(theta_update_HHSSI_list[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.025)  # J × K
beta_upper <- apply(theta_update_HHSSI_list[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.975)  # J × K
beta_true <- sim_data_list[[sample]]$beta_true

gamma_post_mean <- apply(theta_update_HHSSI_list[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), mean)  # M × K
gamma_lower <- apply(theta_update_HHSSI_list[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.025)  # M × K
gamma_upper <- apply(theta_update_HHSSI_list[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.975)  # M × K
gamma_true <- sim_data_list[[sample]]$gamma_true

delta_post_mean <- apply(theta_update_HHSSI_list[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), mean)  # JM × K
delta_lower <- apply(theta_update_HHSSI_list[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.025)  # JM × K
delta_upper <- apply(theta_update_HHSSI_list[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.975)  # JM × K
delta_true <- sim_data_list[[sample]]$delta_true


alpha_df <- data.frame(
  coef = paste0("A", "_k", rep(1:K)),
  mean = as.vector(alpha_post_mean),
  lower = as.vector(alpha_lower),
  upper = as.vector(alpha_upper),
  true  = as.vector(alpha_true))
ggplot(alpha_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(y = "Alpha coefficient", x = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of alpha")


beta_df <- data.frame(
  coef = paste0("B", rep(1:J, times = K), "_k", rep(1:K, each = J)),
  mean = as.vector(beta_post_mean),
  lower = as.vector(beta_lower),
  upper = as.vector(beta_upper),
  true  = as.vector(beta_true))
ggplot(beta_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(y = "Beta coefficient", x = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of beta")


gamma_df <- data.frame(
  coef = paste0("g", rep(1:M, times = K), "_k", rep(1:K, each = M)),
  mean = as.vector(gamma_post_mean),
  lower = as.vector(gamma_lower),
  upper = as.vector(gamma_upper),
  true  = as.vector(gamma_true))
ggplot(gamma_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(x = "gamma coefficient", y = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of gamma")


delta_df <- data.frame(
  coef = paste0("d", rep(1:(J*M), times = K), "_k", rep(1:K, each = (J*M))),
  mean = as.vector(delta_post_mean),
  lower = as.vector(delta_lower),
  upper = as.vector(delta_upper),
  true  = as.vector(delta_true))
ggplot(delta_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(x = "delta coefficient", y = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of delta")

#####################################################################################





# delta true vs estimated heatmap
library(ggplot2)

# Select one dataset and one outcome, if your matrices contain K outcomes
k <- 1

delta_est_array <- simplify2array(delta_est_list_HHSSI)
delta_est_HHSSI <- apply(delta_est_array, c(1, 2), mean)



delta_est_array <- simplify2array(delta_est_list_NHHSSI)
delta_est_NHHSSI <- apply(delta_est_array, c(1, 2), mean)




# If these are J*M x K matrices:
true_delta <- delta_true[, k]
est_hhssi <- delta_est_HHSSI[, k]
est_nhhssi <- delta_est_NHHSSI[, k]

# Calculate estimation errors
error_hhssi <- est_hhssi - true_delta
error_nhhssi <- est_nhhssi - true_delta

# Create heatmap data
df <- expand.grid(
  exposure = 1:M,
  modifier = 1:J
)

df$error_HHSSI <- error_hhssi
df$error_NHHSSI <- error_nhhssi


df_long <- rbind(
  data.frame(
    exposure = df$exposure,
    modifier = df$modifier,
    error = df$error_HHSSI,
    Method = "HHSSI"
  ),
  data.frame(
    exposure = df$exposure,
    modifier = df$modifier,
    error = df$error_NHHSSI,
    Method = "NHHSSI"
  )
)

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

# Use a common color scale centered at zero
error_range <- max(abs(df_long$error), na.rm = TRUE)

ggplot(df_long, aes(x = exposure, y = modifier, fill = error)) +
  geom_tile() +
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  facet_wrap(~Method) +
  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = 1:J) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-error_range, error_range),
    name = "Estimate - True"
  ) +
  labs(
    title = "Heatmap of Interaction Estimation Errors",
    subtitle = "Black borders indicate active interactions",
    x = "Exposure index",
    y = "Modifier index"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )




par(mfrow = c(1, 2))
matplot(delta_var_indices_mat_HHSSI, type = "l")
matplot(delta_var_indices_mat_NHHSSI, type = "l")

par(mfrow = c(1, 1))
apply(delta_var_indices_mat_HHSSI, 2, mean)
apply(delta_var_indices_mat_NHHSSI, 2, mean)
plot(sqrt(apply(delta_var_indices_mat_HHSSI, 2, mean)),
     sqrt(apply(delta_var_indices_mat_NHHSSI, 2, mean)))


plot(error_hhssi, error_nhhssi)
plot(apply(delta_est_HHSSI, 1, mean), apply(delta_est_NHHSSI, 1, mean))

plot(apply(delta_true, 1, mean), apply(delta_est_NHHSSI, 1, mean))
plot(apply(delta_true, 1, mean), apply(delta_est_HHSSI, 1, mean))

#########################

library(ggplot2)

# Select one outcome
for(k in 1:4)
{
# Extract true and estimated interaction coefficients
true_delta <- delta_true[, k]
est_hhssi <- delta_est_HHSSI[, k]
est_nhhssi <- delta_est_NHHSSI[, k]

# Create heatmap data
df <- expand.grid(
  exposure = 1:M,
  modifier = 1:J
)

df$True <- true_delta
df$HHSSI <- est_hhssi
df$NHHSSI <- est_nhhssi

# Convert to long format
df_long <- rbind(
  data.frame(
    exposure = df$exposure,
    modifier = df$modifier,
    delta = df$True,
    Method = "True"
  ),
  data.frame(
    exposure = df$exposure,
    modifier = df$modifier,
    delta = df$HHSSI,
    Method = "HHSSI"
  ),
  data.frame(
    exposure = df$exposure,
    modifier = df$modifier,
    delta = df$NHHSSI,
    Method = "NHHSSI"
  )
)

# Identify active interactions: strong + weak heredity
df_long$Interaction_idx <- (df_long$modifier - 1) * M + df_long$exposure
strong_idx <- (1:2 - 1) * M + rep(1:3, each = 2)
fixed_weak_rows <- sim_data_list[[1]]$active_weak_idx
active_delta_idx <- c(strong_idx, fixed_weak_rows)
df_long$Active <- df_long$Interaction_idx %in% active_delta_idx

# Use the same color scale for all panels
delta_range <- range(df_long$delta, na.rm = TRUE)
delta_limit <- max(abs(delta_range))

g1 <- ggplot(df_long, aes(x = exposure, y = modifier, fill = delta)) +
  geom_tile() +
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  facet_wrap(~Method, nrow = 1) +
  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = 1:J) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-delta_limit, delta_limit),
    name = expression(delta)
  ) +
  labs(
    title = "True and Estimated Modifier–Exposure Interactions",
    subtitle = "Black borders indicate active interactions",
    x = "Exposure index",
    y = "Modifier index"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(g1)
}



################################################

colnames(x) # "h_mbmi_None", "e3_gac_None" are first two modifiers
colnames(z) # "hs_cu_c_Log2" is 13th exposure

u13_1 <- z[, 13] * x[, 1]
u13_2 <- z[, 13] * x[, 2]
par(mfrow = c(1, 2))
hist(u13_1,
     main = "Exposure 13 × Modifier 1",
     xlab = "Interaction value",
     col = "lightblue",
     breaks = 30)
hist(u13_2,
     main = "Exposure 13 × Modifier 2",
     xlab = "Interaction value",
     col = "lightpink",
     breaks = 30)
par(mfrow = c(1, 1))




interaction_summary <- data.frame(
  Interaction = c("Exposure 13 × Modifier 1",
                  "Exposure 13 × Modifier 2"),
  Mean = c(mean(u13_1), mean(u13_2)),
  SD = c(sd(u13_1), sd(u13_2)),
  Min = c(min(u13_1), min(u13_2)),
  Max = c(max(u13_1), max(u13_2)),
  Near_zero = c(mean(abs(u13_1) < 0.01),
                mean(abs(u13_2) < 0.01))
)
interaction_summary





active_rows_1 <- which(rowSums(sim_data_list[[1]]$delta_true != 0) > 0)
active_rows_2 <- which(rowSums(sim_data_list[[2]]$delta_true != 0) > 0)

active_rows_1==active_rows_2


sim_data_list[[1]]$delta_true[active_rows_1,]


> load("~/Desktop/Packages/BVSMEInteractionMVN/practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_100weak.rda")
>

  sim_data_list[[1]]$delta_true[active_rows_1,]
[,1]        [,2]         [,3]         [,4]
[1,] -2.252166920 -0.20222937  0.160336953 -0.336736500
[2,] -1.221610967  0.07745647  0.757535286  0.501317751
[3,] -1.523022677 -0.47460808  0.399381202 -0.816431361
[4,] -0.199301069  0.02112817 -1.095849791  0.976950979
[5,]  0.680208029 -1.44302315  0.018852536 -0.025120198
[6,]  0.922643544  0.84898056 -0.764372218  0.822848954
[7,]  1.010059141 -0.46032589 -0.790950120 -1.760285239
[8,] -0.702918669 -0.27313316  0.858575621  1.606054493
[9,] -0.528797448  0.14868192 -0.486028608 -0.258187102
[10,] -0.117604028 -0.23374565  0.467275258  0.833895880
[11,]  0.295110299  0.43023417 -1.131580128  0.032067906
[12,]  2.501703451 -1.76230689  0.569926204  0.643144320
[13,] -1.669046890 -2.06379967  0.020707514  1.038281965
[14,] -1.039413023  0.59786117 -1.782537164 -0.472364943
[15,] -0.824409608 -0.70407318  1.137913898  2.310675514
[16,] -0.958409297 -0.80523339  0.055705848  0.423924556
[17,] -0.126475143 -0.77531433 -1.512409318 -2.238780972
[18,]  1.269887915 -1.00522531 -0.185366036 -0.505531074
[19,]  0.445386316  0.26003309 -0.005102893 -0.460773400
[20,] -0.175340009  0.98780410 -1.612186691 -0.675167127
[21,]  1.111680527 -1.09612789  1.490501033 -1.786839032
[22,]  0.175211753 -0.76252598  0.962129583  0.674059403
[23,]  1.800413267 -0.30009211  0.710042617  0.873797956
[24,] -1.210836140  0.59966597 -1.254792771 -1.446873446
[25,]  2.091496146  1.49945031  0.112740117 -2.305004970
[26,]  0.561966146  1.08386078  0.432585690 -0.499982120
[27,]  1.436922850  0.20814754  1.438529423  1.435991895
[28,] -0.442853976  0.59569307  1.148663817 -1.822068011
[29,]  0.757077854  1.63076447  1.239806524 -0.059770301
[30,] -0.984055349 -1.65551175 -0.740810868  0.309353365
[31,]  0.423393947 -1.67121060 -0.683758735 -0.354303298
[32,]  0.365576375 -1.58016771  1.117359066  0.167805700
[33,]  0.502816090  0.22141529 -1.367068160 -1.879045943
[34,] -0.281287105 -0.95332345 -0.590760733 -0.869379708
[35,]  0.416047834  0.73955005  0.639801812 -0.380807048
[36,]  0.626422178  0.76705397  0.870979603 -0.176324211
[37,]  1.493352599  0.72288010  0.131485590  1.958905313
[38,]  0.122343348 -1.60733338  0.843373445  0.418221015
[39,] -1.515491717  0.08211289 -0.781604179  1.739228962
[40,] -0.409628593 -0.03282725  0.496145645  1.078099072
[41,]  1.470116393  0.65753928  0.594379410 -0.791862310
[42,] -0.780406855  0.03417712  1.403910882  0.513951726
[43,]  2.642018174 -1.05831999  1.162548743 -0.390920799
[44,]  1.363655579  0.49198953 -0.106733374  0.546712005
[45,]  0.748717201  0.73454793  0.910247651 -0.613827450
[46,]  0.796223152 -0.42196123  0.035971604  0.231480731
[47,]  0.719866623 -0.36665232  0.722324417 -1.573184023
[48,]  1.534439659  1.07406608  1.509609990 -0.731265404
[49,] -0.002296398 -1.16076327 -0.921947553  0.418028924
[50,] -1.034194753 -1.57398600  0.270412104  0.692521111
[51,] -1.700385914 -1.53430899  1.668648337 -0.194317869
[52,]  0.061841061  0.02633444  1.255174517  0.241327127
[53,] -0.690356744  0.31483623  1.280163323 -0.677253317
[54,] -1.525897344  1.86873044  0.466958527  0.243660072
[55,]  0.434434599  1.36473323  0.723247652 -0.387526177
[56,]  0.574678067 -0.65390323 -1.271858248 -0.343877484
[57,] -0.513234688 -0.58964855  1.395451854  1.592051516
[58,] -1.576183135 -0.69118475 -1.286563034 -1.366094827
[59,]  0.814238758 -0.45455818  0.985527046 -0.551306921
[60,] -1.239379414 -0.10289446 -0.006642433  0.678580687
[61,]  0.969998133 -0.38479788 -0.718293633 -0.237769168
[62,]  0.035269383  1.05080551  0.248541508 -0.840521564
[63,] -0.610548776  0.01623400 -0.478778112  1.568112369
[64,]  1.900211027  0.04907783  1.229110312 -0.869086683
[65,] -0.508408595  1.07225884 -0.840967957  0.493298616
[66,] -0.714406733 -0.14669563  0.593968143  2.529347500
[67,] -0.388983167 -1.23872145 -0.346405434 -0.506915363
[68,] -1.567809123 -0.45510203  0.956482787  0.068595149
[69,] -0.141686116 -0.12384074 -1.627885653 -0.542823713
[70,]  1.143857240  1.17452491 -0.735074835  0.952547682
[71,] -1.047465799 -1.76302587  1.386285827  1.335563858
[72,]  1.336870970 -0.76467398  0.071092183  0.699623370
[73,] -0.392116583  0.57781711 -0.366590320  0.521054522
[74,]  1.937157584 -0.71142084  0.396970165 -0.648352199
[75,] -0.823345394  0.26883815 -1.230530778 -0.931249653
[76,] -1.080664124  0.99066082  0.279554738  0.735120526
[77,] -0.648895561  0.38797422 -1.581850000  0.440465335
[78,]  0.241565519 -0.22470259  1.507137413 -0.009413817
[79,]  0.106935157 -0.04594134  0.509591826  0.125642717
[80,] -0.546465144  0.73661823 -0.084261586 -0.387327494
[81,]  0.465389541  1.51813029 -1.327976910 -1.325377942
[82,] -1.188306405  0.66548368  0.607438438 -1.399337244
[83,]  0.561697818  0.38660740  0.302437489  0.885130915
[84,] -0.281474703  1.92726767 -0.475504068 -0.047648598
[85,]  0.962156642  0.30222523  0.815188549 -1.934732345
[86,]  1.630763322  0.38732655  0.626712381 -0.106839785
[87,] -3.026819465  1.05853975 -0.863710260 -2.172060792
[88,]  0.157215460  0.90501340 -0.604345672 -0.956440724
[89,] -0.773697963 -1.03674921  0.717205506 -0.087568068
[90,]  0.491168333  0.01073137 -0.262233885  0.110812791
[91,]  1.113372856 -0.84267590  0.627781357 -0.169461758
[92,]  0.377278137 -1.17381053 -0.838872413 -0.954337473
[93,]  0.641981577 -0.75416318 -1.554703955 -0.715168631
[94,] -0.995570016 -3.55946128 -0.989958768 -0.586893105
[95,]  0.018251568 -1.22624037 -0.182828253  1.119331614
[96,] -1.486050869 -0.99974896  1.928209557 -0.717549626
[97,]  0.828759799  0.88974147  0.366522404 -1.255865303
[98,]  2.515457280 -0.55171780  0.910028734 -0.377813637
[99,] -1.956208950  0.04030571  1.114538300  0.100029264
[100,] -0.340540402  1.47072734  0.286311414  0.599341891
[101,]  0.349898547  0.23711243 -0.783784089  1.598149755
[102,]  0.689846454  0.43042219  0.372553060  1.291913389
[103,] -0.683941688 -0.22937695 -0.781629022  0.546980236
[104,]  0.598281860 -0.51627179 -1.404806738 -0.560852555
[105,]  0.796726807  0.27369852 -2.222426320  0.646081630
[106,] -0.696190981  0.73382868 -0.087028227 -1.762156554
[107,]  1.476636232 -0.60445304  0.183933408  0.565352400
[108,] -0.903177268 -1.03797686 -0.755931817  1.915188831
[109,]  0.331542742  0.43695162 -0.834528530  1.559285898
[110,]  0.075530537 -0.21202756 -0.274559124 -0.866631919
[111,] -0.011144129 -0.45327174  0.879230223  1.663311569
[112,]  0.906328582  1.38105255 -0.316966848  0.408397765
[113,] -1.060990743 -0.21543762 -0.358512129  1.484962929
[114,] -0.877137466 -1.86029536  1.408139915 -1.000016601
[115,]  0.132569877  0.05593941  0.565550170 -1.045799809
[116,]  0.207043300 -1.08751710  0.439087197 -0.110629848
[117,]  1.415784132 -0.75207447  0.138764334 -1.704350027
[118,] -1.259188094 -0.08190110 -0.890491599  0.046688554
[119,] -0.455298039 -0.83789083 -0.705883963  0.171900583
[120,]  3.046300789  1.03123264 -0.670190930 -1.529976821
[121,] -0.053587695 -0.29083165 -1.815763587 -1.321555841
[122,] -0.731411818  1.72223157  0.315233317 -1.332891855
[123,]  0.945679138  0.16422592 -0.313386585 -1.815565243
[124,] -0.065308703  0.78492201 -1.373259028  0.873851421
[125,] -0.304040198  0.12571232 -1.305712810  0.064117331
[126,]  0.850391422  1.35713290  0.193680018  0.198351416
[127,] -0.573243388 -0.88639220  0.915111486  0.308976764
[128,] -1.215993581  0.35053104 -1.610724830  0.231735697
[129,]  0.655348477  0.19202183 -0.160343066 -0.445569733
[130,]  1.444910754 -1.93538729  1.130204242  1.375977624
[131,] -1.764280558  2.54676937 -0.041256987  0.551648985
[132,]  0.627553595 -1.32768238  0.737025563  0.129585427
[133,]  1.061272330  0.37919512  1.602924836  0.053555798
[134,]  0.540950235 -3.25636754  0.258322081 -1.767789775
[135,] -1.977507104  1.20842132 -0.658097047  0.369922887
[136,]  0.876762597  0.65228005  0.136884749  0.927505618
[137,] -1.740029197 -0.90913289  0.602793506  0.431304703
[138,]  1.310890216  0.22324192 -0.084753774  1.757455384
[139,]  0.418252297  0.96440028  0.520997065 -0.370822297
[140,] -0.322444143  0.51438645 -0.160957770  0.574035676
[141,] -0.330924792  0.32280154  0.671650261 -0.066625513
[142,]  0.119714653  1.99166204 -0.339151913  2.915143567
[143,] -0.580424994  0.21485983  1.021537002 -0.728263246
[144,]  0.034911001  0.83042249  1.287330509 -1.008030364
[145,] -0.459372903 -0.60603872 -0.416769682  0.762973417
[146,] -0.478908072  0.76483575 -0.046922524 -1.407947981
[147,]  0.079498021  1.11989328 -0.890193272  1.302286886
[148,] -0.830029434 -0.72902389 -0.844373423  0.602902830
[149,] -1.222182206 -0.36218347 -1.145549202 -1.541704901
[150,]  1.922200059  0.38857583 -1.084884769 -1.731259010
[151,]  0.094805220  0.09970004  1.007438418  0.474187276
[152,] -1.284990455  0.20420664  1.602995348  0.696468546
[153,]  1.212936374 -1.34748728  0.434210581  0.503858495
[154,] -1.646053905  1.13783393  0.953630280  1.972679627
[155,]  0.938129278 -1.34603272  0.280387298 -0.437322080
[156,]  1.253338669  0.28395954 -0.427129347 -0.136282012
[157,] -0.894409872  1.48303598 -0.557350367  0.825354871
[158,] -1.160775967 -0.97413227  0.157607998  1.789634990
[159,]  1.460739507 -2.21180602  0.452721508 -0.951894385
[160,] -0.020203773 -0.70041595  0.553959731 -0.203659506
>
