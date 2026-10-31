
rm(list = ls())
gc()

# load the simulated data
# load("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")

# load("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")

# load("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")

# load("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")

# load("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")


# load("practice_data/results_datasimulate_simvar_regcoef_normal0_1_modify5.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5.rda")


# load("practice_data/results_datasimulate_simvar_regcoef_unif1_2_modify5.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5.rda")



# load("practice_data/results_datasimulate_simvar_regcoef_unif2_3_modify5.rda")
# load("practice_data/results_datasimulate_realvar_regcoef_unif2_3_modify5.rda")
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

# load("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")

# load("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")

# load("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")

# load("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")


# load("practice_data/results_HHS_simvar_regcoef_normal0_1_modify5.rda")
# load("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5.rda")


# load("practice_data/results_HHS_simvar_regcoef_unif1_2_modify5.rda")
# load("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5.rda")


# load("practice_data/results_HHS_simvar_regcoef_unif2_3_modify5.rda")
# load("practice_data/results_HHS_realvar_regcoef_unif2_3_modify5.rda")




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

beta_bias_active_HHS <- c()
gamma_bias_active_HHS <- c()
varphi_bias_active_HHS <- c()
delta_bias_active_HHS <- c()

beta_bias_inactive_HHS <- c()
gamma_bias_inactive_HHS <- c()
varphi_bias_inactive_HHS <- c()
delta_bias_inactive_HHS <- c()

beta_RMSE_HHS <- c()
gamma_RMSE_HHS <- c()
varphi_RMSE_HHS <- c()
delta_RMSE_HHS <- c()

beta_RMSE_active_HHS <- c()
gamma_RMSE_active_HHS <- c()
varphi_RMSE_active_HHS <- c()
delta_RMSE_active_HHS <- c()

beta_RMSE_inactive_HHS <- c()
gamma_RMSE_inactive_HHS <- c()
varphi_RMSE_inactive_HHS <- c()
delta_RMSE_inactive_HHS <- c()

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
  beta_bias_active_HHS[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_inactive_HHS[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_HHS[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_active_HHS[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_inactive_HHS[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



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
  gamma_bias_active_HHS[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_inactive_HHS[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_HHS[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_active_HHS[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_inactive_HHS[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

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
  varphi_bias_active_HHS[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_inactive_HHS[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_HHS[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_active_HHS[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_inactive_HHS[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


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
  delta_bias_active_HHS[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_inactive_HHS[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_HHS[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_active_HHS[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_inactive_HHS[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

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

round(mean(alpha_cover_per_HHS), 2)

round(mean(varphi_cover_per_HHS), 2)
round(mean(varphi_cover_nonzero_per_HHS), 2)
round(mean(delta_cover_nonzero_per_HHS), 2)

round(mean(alpha_est_mat_HHS), 2)

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

cor(beta_est_mean_HHS, beta_true_mean_HHS)

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

cor(gamma_est_mean_HHS, gamma_true_mean_HHS)

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

cor(delta_est_mean_HHS, delta_true_mean_HHS)


varphi_true_array <- simplify2array(varphi_true_list_HHS)
varphi_true_mean_HHS <- apply(varphi_true_array, c(1, 2), mean)
varphi_est_array <- simplify2array(varphi_est_list_HHS)
varphi_est_mean_HHS <- apply(varphi_est_array, c(1, 2), mean)
dim(varphi_est_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_true_mean_HHS[, k],
    varphi_est_mean_HHS[, k],
    xlab = "True varphi",
    ylab = "Mean estimated varphi",
    main = paste("Outcome", k),
    pch = 19
  )
  abline(0, 1, lty = 2)
}

cor(varphi_est_mean_HHS, varphi_true_mean_HHS)


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

beta_bias_active_array <- simplify2array(beta_bias_active_HHS)
beta_bias_active_mean_HHS <- apply(beta_bias_active_array, c(1, 2), mean)
mean(beta_bias_active_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_active_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_inactive_array <- simplify2array(beta_bias_inactive_HHS)
beta_bias_inactive_mean_HHS <- apply(beta_bias_inactive_array, c(1, 2), mean)
mean(beta_bias_inactive_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_inactive_mean_HHS[, k],
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

gamma_bias_active_array <- simplify2array(gamma_bias_active_HHS)
gamma_bias_active_mean_HHS <- apply(gamma_bias_active_array, c(1, 2), mean)
mean(gamma_bias_active_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_active_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_inactive_array <- simplify2array(gamma_bias_inactive_HHS)
gamma_bias_inactive_mean_HHS <- apply(gamma_bias_inactive_array, c(1, 2), mean)
mean(gamma_bias_inactive_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_inactive_mean_HHS[, k],
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


varphi_bias_active_array <- simplify2array(varphi_bias_active_HHS)
varphi_bias_active_mean_HHS <- apply(varphi_bias_active_array, c(1, 2), mean)
mean(varphi_bias_active_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_active_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_inactive_array <- simplify2array(varphi_bias_inactive_HHS)
varphi_bias_inactive_mean_HHS <- apply(varphi_bias_inactive_array, c(1, 2), mean)
mean(varphi_bias_inactive_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_inactive_mean_HHS[, k],
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



delta_bias_active_array <- simplify2array(delta_bias_active_HHS)
delta_bias_active_mean_HHS <- apply(delta_bias_active_array, c(1, 2), mean)
mean(delta_bias_active_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_active_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_inactive_array <- simplify2array(delta_bias_inactive_HHS)
delta_bias_inactive_mean_HHS <- apply(delta_bias_inactive_array, c(1, 2), mean)
mean(delta_bias_inactive_mean_HHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_inactive_mean_HHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


round(mean(beta_RMSE_HHS)*100, 2)
round(mean(beta_RMSE_active_HHS)*100, 2)
round(mean(beta_RMSE_inactive_HHS)*100, 2)

round(mean(gamma_RMSE_HHS)*100, 2)
round(mean(gamma_RMSE_active_HHS)*100, 2)
round(mean(gamma_RMSE_inactive_HHS)*100, 2)

round(mean(varphi_RMSE_HHS)*100, 2)
round(mean(varphi_RMSE_active_HHS)*100, 2)
round(mean(varphi_RMSE_inactive_HHS)*100, 2)

round(mean(delta_RMSE_HHS)*100, 2)
round(mean(delta_RMSE_active_HHS)*100, 2)
round(mean(delta_RMSE_inactive_HHS)*100, 2)

round(mean(beta_sensitivity_HHS)*100, 2)
round(mean(gamma_sensitivity_HHS)*100, 2)
round(mean(varphi_sensitivity_HHS)*100, 2)

round(mean(beta_specificity_HHS)*100, 2)
round(mean(gamma_specificity_HHS)*100, 2)
round(mean(varphi_specificity_HHS)*100, 2)

round(mean(beta_precision_HHS)*100, 2)
round(mean(gamma_precision_HHS)*100, 2)
round(mean(varphi_precision_HHS)*100, 2)





#######################################################

# HHSSI model

# load("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")

# load("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")

# load("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")

# load("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")

# load("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")


# load("practice_data/results_HHSSI_simvar_regcoef_normal0_1_modify5.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5.rda")


# load("practice_data/results_HHSSI_simvar_regcoef_unif1_2_modify5.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5.rda")


# load("practice_data/results_HHSSI_simvar_regcoef_unif2_3_modify5.rda")
# load("practice_data/results_HHSSI_realvar_regcoef_unif2_3_modify5.rda")




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

beta_bias_active_HHSSI <- c()
gamma_bias_active_HHSSI <- c()
varphi_bias_active_HHSSI <- c()
delta_bias_active_HHSSI <- c()

beta_bias_inactive_HHSSI <- c()
gamma_bias_inactive_HHSSI <- c()
varphi_bias_inactive_HHSSI <- c()
delta_bias_inactive_HHSSI <- c()

beta_RMSE_HHSSI <- c()
gamma_RMSE_HHSSI <- c()
varphi_RMSE_HHSSI <- c()
delta_RMSE_HHSSI <- c()

beta_RMSE_active_HHSSI <- c()
gamma_RMSE_active_HHSSI <- c()
varphi_RMSE_active_HHSSI <- c()
delta_RMSE_active_HHSSI <- c()

beta_RMSE_inactive_HHSSI <- c()
gamma_RMSE_inactive_HHSSI <- c()
varphi_RMSE_inactive_HHSSI <- c()
delta_RMSE_inactive_HHSSI <- c()

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
  beta_bias_active_HHSSI[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_inactive_HHSSI[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_HHSSI[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_active_HHSSI[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_inactive_HHSSI[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



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
  gamma_bias_active_HHSSI[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_inactive_HHSSI[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_HHSSI[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_active_HHSSI[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_inactive_HHSSI[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

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
  varphi_bias_active_HHSSI[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_inactive_HHSSI[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_HHSSI[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_active_HHSSI[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_inactive_HHSSI[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


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
  delta_bias_active_HHSSI[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_inactive_HHSSI[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_HHSSI[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_active_HHSSI[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_inactive_HHSSI[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

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

round(mean(alpha_cover_per_HHSSI), 2)

round(mean(varphi_cover_per_HHSSI), 2)
round(mean(varphi_cover_nonzero_per_HHSSI), 2)
round(mean(delta_cover_nonzero_per_HHSSI), 2)

round(mean(alpha_est_mat_HHSSI), 2)

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


beta_bias_active_array <- simplify2array(beta_bias_active_HHSSI)
beta_bias_active_mean_HHSSI <- apply(beta_bias_active_array, c(1, 2), mean)
mean(beta_bias_active_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_active_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_inactive_array <- simplify2array(beta_bias_inactive_HHSSI)
beta_bias_inactive_mean_HHSSI <- apply(beta_bias_inactive_array, c(1, 2), mean)
mean(beta_bias_inactive_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_inactive_mean_HHSSI[, k],
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

gamma_bias_active_array <- simplify2array(gamma_bias_active_HHSSI)
gamma_bias_active_mean_HHSSI <- apply(gamma_bias_active_array, c(1, 2), mean)
mean(gamma_bias_active_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_active_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_inactive_array <- simplify2array(gamma_bias_inactive_HHSSI)
gamma_bias_inactive_mean_HHSSI <- apply(gamma_bias_inactive_array, c(1, 2), mean)
mean(gamma_bias_inactive_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_inactive_mean_HHSSI[, k],
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


varphi_bias_active_array <- simplify2array(varphi_bias_active_HHSSI)
varphi_bias_active_mean_HHSSI <- apply(varphi_bias_active_array, c(1, 2), mean)
mean(varphi_bias_active_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_active_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_inactive_array <- simplify2array(varphi_bias_inactive_HHSSI)
varphi_bias_inactive_mean_HHSSI <- apply(varphi_bias_inactive_array, c(1, 2), mean)
mean(varphi_bias_inactive_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_inactive_mean_HHSSI[, k],
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



delta_bias_active_array <- simplify2array(delta_bias_active_HHSSI)
delta_bias_active_mean_HHSSI <- apply(delta_bias_active_array, c(1, 2), mean)
mean(delta_bias_active_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_active_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_inactive_array <- simplify2array(delta_bias_inactive_HHSSI)
delta_bias_inactive_mean_HHSSI <- apply(delta_bias_inactive_array, c(1, 2), mean)
mean(delta_bias_inactive_mean_HHSSI)

par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_inactive_mean_HHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


round(mean(beta_RMSE_HHSSI)*100, 2)
round(mean(beta_RMSE_active_HHSSI)*100, 2)
round(mean(beta_RMSE_inactive_HHSSI)*100, 2)

round(mean(gamma_RMSE_HHSSI)*100, 2)
round(mean(gamma_RMSE_active_HHSSI)*100, 2)
round(mean(gamma_RMSE_inactive_HHSSI)*100, 2)

round(mean(varphi_RMSE_HHSSI)*100, 2)
round(mean(varphi_RMSE_active_HHSSI)*100, 2)
round(mean(varphi_RMSE_inactive_HHSSI)*100, 2)

round(mean(delta_RMSE_HHSSI)*100, 2)
round(mean(delta_RMSE_active_HHSSI)*100, 2)
round(mean(delta_RMSE_inactive_HHSSI)*100, 2)

round(mean(beta_sensitivity_HHSSI)*100, 2)
round(mean(gamma_sensitivity_HHSSI)*100, 2)
round(mean(varphi_sensitivity_HHSSI)*100, 2)

round(mean(beta_specificity_HHSSI)*100, 2)
round(mean(gamma_specificity_HHSSI)*100, 2)
round(mean(varphi_specificity_HHSSI)*100, 2)

round(mean(beta_precision_HHSSI)*100, 2)
round(mean(gamma_precision_HHSSI)*100, 2)
round(mean(varphi_precision_HHSSI)*100, 2)




######################################
# NHHS model

# load("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")

# load("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")

# load("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")

# load("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")


# load("practice_data/results_NHHS_simvar_regcoef_normal0_1_modify5.rda")
# load("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5.rda")


# load("practice_data/results_NHHS_simvar_regcoef_unif1_2_modify5.rda")
# load("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5.rda")



# load("practice_data/results_NHHS_simvar_regcoef_unif2_3_modify5.rda")
# load("practice_data/results_NHHS_realvar_regcoef_unif2_3_modify5.rda")





nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

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

beta_bias_active_NHHS <- c()
gamma_bias_active_NHHS <- c()
varphi_bias_active_NHHS <- c()
delta_bias_active_NHHS <- c()

beta_bias_inactive_NHHS <- c()
gamma_bias_inactive_NHHS <- c()
varphi_bias_inactive_NHHS <- c()
delta_bias_inactive_NHHS <- c()

beta_RMSE_NHHS <- c()
gamma_RMSE_NHHS <- c()
varphi_RMSE_NHHS <- c()
delta_RMSE_NHHS <- c()

beta_RMSE_active_NHHS <- c()
gamma_RMSE_active_NHHS <- c()
varphi_RMSE_active_NHHS <- c()
delta_RMSE_active_NHHS <- c()

beta_RMSE_inactive_NHHS <- c()
gamma_RMSE_inactive_NHHS <- c()
varphi_RMSE_inactive_NHHS <- c()
delta_RMSE_inactive_NHHS <- c()

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
  beta_bias_active_NHHS[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_inactive_NHHS[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_NHHS[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_active_NHHS[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_inactive_NHHS[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



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
  gamma_bias_active_NHHS[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_inactive_NHHS[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_NHHS[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_active_NHHS[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_inactive_NHHS[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

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
  varphi_bias_active_NHHS[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_inactive_NHHS[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_NHHS[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_active_NHHS[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_inactive_NHHS[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


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
  delta_bias_active_NHHS[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_inactive_NHHS[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_NHHS[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_active_NHHS[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_inactive_NHHS[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

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

round(mean(alpha_cover_per_NHHS), 2)

round(mean(beta_cover_per_NHHS), 2)
round(mean(beta_cover_nonzero_per_NHHS), 2)

round(mean(gamma_cover_per_NHHS), 2)
round(mean(gamma_cover_nonzero_per_NHHS), 2)

round(mean(varphi_cover_per_NHHS), 2)
round(mean(varphi_cover_nonzero_per_NHHS), 2)


round(mean(alpha_est_mat_NHHS), 2)

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



varphi_true_array <- simplify2array(varphi_true_list_NHHS)
varphi_true_mean_NHHS <- apply(varphi_true_array, c(1, 2), mean)
varphi_est_array <- simplify2array(varphi_est_list_NHHS)
varphi_est_mean_NHHS <- apply(varphi_est_array, c(1, 2), mean)
dim(varphi_est_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_true_mean_NHHS[, k],
    varphi_est_mean_NHHS[, k],
    xlab = "True varphi",
    ylab = "Mean estimated varphi",
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


beta_bias_active_array <- simplify2array(beta_bias_active_NHHS)
beta_bias_active_mean_NHHS <- apply(beta_bias_active_array, c(1, 2), mean)
mean(beta_bias_active_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_active_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_inactive_array <- simplify2array(beta_bias_inactive_NHHS)
beta_bias_inactive_mean_NHHS <- apply(beta_bias_inactive_array, c(1, 2), mean)
mean(beta_bias_inactive_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_inactive_mean_NHHS[, k],
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

gamma_bias_active_array <- simplify2array(gamma_bias_active_NHHS)
gamma_bias_active_mean_NHHS <- apply(gamma_bias_active_array, c(1, 2), mean)
mean(gamma_bias_active_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_active_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_inactive_array <- simplify2array(gamma_bias_inactive_NHHS)
gamma_bias_inactive_mean_NHHS <- apply(gamma_bias_inactive_array, c(1, 2), mean)
mean(gamma_bias_inactive_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_inactive_mean_NHHS[, k],
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


varphi_bias_active_array <- simplify2array(varphi_bias_active_NHHS)
varphi_bias_active_mean_NHHS <- apply(varphi_bias_active_array, c(1, 2), mean)
mean(varphi_bias_active_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_active_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_inactive_array <- simplify2array(varphi_bias_inactive_NHHS)
varphi_bias_inactive_mean_NHHS <- apply(varphi_bias_inactive_array, c(1, 2), mean)
mean(varphi_bias_inactive_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_inactive_mean_NHHS[, k],
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



delta_bias_active_array <- simplify2array(delta_bias_active_NHHS)
delta_bias_active_mean_NHHS <- apply(delta_bias_active_array, c(1, 2), mean)
mean(delta_bias_active_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_active_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_inactive_array <- simplify2array(delta_bias_inactive_NHHS)
delta_bias_inactive_mean_NHHS <- apply(delta_bias_inactive_array, c(1, 2), mean)
mean(delta_bias_inactive_mean_NHHS)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_inactive_mean_NHHS[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


round(mean(beta_RMSE_NHHS)*100, 2)
round(mean(beta_RMSE_active_NHHS)*100, 2)
round(mean(beta_RMSE_inactive_NHHS)*100, 2)

round(mean(gamma_RMSE_NHHS)*100, 2)
round(mean(gamma_RMSE_active_NHHS)*100, 2)
round(mean(gamma_RMSE_inactive_NHHS)*100, 2)

round(mean(varphi_RMSE_NHHS)*100, 2)
round(mean(varphi_RMSE_active_NHHS)*100, 2)
round(mean(varphi_RMSE_inactive_NHHS)*100, 2)

round(mean(delta_RMSE_NHHS)*100, 2)
round(mean(delta_RMSE_active_NHHS)*100, 2)
round(mean(delta_RMSE_inactive_NHHS)*100, 2)

round(mean(beta_sensitivity_NHHS)*100, 2)
round(mean(gamma_sensitivity_NHHS)*100, 2)
round(mean(varphi_sensitivity_NHHS)*100, 2)

round(mean(beta_specificity_NHHS)*100, 2)
round(mean(gamma_specificity_NHHS)*100, 2)
round(mean(varphi_specificity_NHHS)*100, 2)

round(mean(beta_precision_NHHS)*100, 2)
round(mean(gamma_precision_NHHS)*100, 2)
round(mean(varphi_precision_NHHS)*100, 2)







###################
# NHHSSI model

# load("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")

# load("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")

# load("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")

# load("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")

# load("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")


# load("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_modify5.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5.rda")


# load("practice_data/results_NHHSSI_simvar_regcoef_unif1_2_modify5.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5.rda")


# load("practice_data/results_NHHSSI_simvar_regcoef_unif2_3_modify5.rda")
# load("practice_data/results_NHHSSI_realvar_regcoef_unif2_3_modify5.rda")


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

beta_bias_active_NHHSSI <- c()
gamma_bias_active_NHHSSI <- c()
varphi_bias_active_NHHSSI <- c()
delta_bias_active_NHHSSI <- c()

beta_bias_inactive_NHHSSI <- c()
gamma_bias_inactive_NHHSSI <- c()
varphi_bias_inactive_NHHSSI <- c()
delta_bias_inactive_NHHSSI <- c()

beta_RMSE_NHHSSI <- c()
gamma_RMSE_NHHSSI <- c()
varphi_RMSE_NHHSSI <- c()
delta_RMSE_NHHSSI <- c()

beta_RMSE_active_NHHSSI <- c()
gamma_RMSE_active_NHHSSI <- c()
varphi_RMSE_active_NHHSSI <- c()
delta_RMSE_active_NHHSSI <- c()

beta_RMSE_inactive_NHHSSI <- c()
gamma_RMSE_inactive_NHHSSI <- c()
varphi_RMSE_inactive_NHHSSI <- c()
delta_RMSE_inactive_NHHSSI <- c()

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
  beta_bias_active_NHHSSI[[sample]] <- beta_bias[nonzero_beta]
  beta_bias_inactive_NHHSSI[[sample]] <- beta_bias[zero_beta]

  beta_RMSE_NHHSSI[sample] <- sqrt(mean(beta_bias^2))
  beta_RMSE_active_NHHSSI[sample] <- sqrt(mean((beta_bias[nonzero_beta])^2))
  beta_RMSE_inactive_NHHSSI[sample] <- sqrt(mean((beta_bias[zero_beta])^2))



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
  gamma_bias_active_NHHSSI[[sample]] <- gamma_bias[nonzero_gamma]
  gamma_bias_inactive_NHHSSI[[sample]] <- gamma_bias[zero_gamma]

  gamma_RMSE_NHHSSI[sample] <- sqrt(mean(gamma_bias^2))
  gamma_RMSE_active_NHHSSI[sample] <- sqrt(mean((gamma_bias[nonzero_gamma])^2))
  gamma_RMSE_inactive_NHHSSI[sample] <- sqrt(mean((gamma_bias[zero_gamma])^2))

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
  varphi_bias_active_NHHSSI[[sample]] <- varphi_bias[nonzero_varphi]
  varphi_bias_inactive_NHHSSI[[sample]] <- varphi_bias[zero_varphi]

  varphi_RMSE_NHHSSI[sample] <- sqrt(mean(varphi_bias^2))
  varphi_RMSE_active_NHHSSI[sample] <- sqrt(mean((varphi_bias[nonzero_varphi])^2))
  varphi_RMSE_inactive_NHHSSI[sample] <- sqrt(mean((varphi_bias[zero_varphi])^2))


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
  delta_bias_active_NHHSSI[[sample]] <- delta_bias[nonzero_delta]
  delta_bias_inactive_NHHSSI[[sample]] <- delta_bias[zero_delta]

  delta_RMSE_NHHSSI[sample] <- sqrt(mean(delta_bias^2))
  delta_RMSE_active_NHHSSI[sample] <- sqrt(mean((delta_bias[nonzero_delta])^2))
  delta_RMSE_inactive_NHHSSI[sample] <- sqrt(mean((delta_bias[zero_delta])^2))

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

round(mean(alpha_cover_per_NHHSSI), 2)

round(mean(varphi_cover_per_NHHSSI), 2)
round(mean(varphi_cover_nonzero_per_NHHSSI), 2)
round(mean(delta_cover_nonzero_per_NHHSSI), 2)

round(mean(alpha_est_mat_NHHSSI), 2)

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


beta_bias_active_array <- simplify2array(beta_bias_active_NHHSSI)
beta_bias_active_mean_NHHSSI <- apply(beta_bias_active_array, c(1, 2), mean)
mean(beta_bias_active_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_active_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


beta_bias_inactive_array <- simplify2array(beta_bias_inactive_NHHSSI)
beta_bias_inactive_mean_NHHSSI <- apply(beta_bias_inactive_array, c(1, 2), mean)
mean(beta_bias_inactive_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    beta_bias_inactive_mean_NHHSSI[, k],
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

gamma_bias_active_array <- simplify2array(gamma_bias_active_NHHSSI)
gamma_bias_active_mean_NHHSSI <- apply(gamma_bias_active_array, c(1, 2), mean)
mean(gamma_bias_active_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_active_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


gamma_bias_inactive_array <- simplify2array(gamma_bias_inactive_NHHSSI)
gamma_bias_inactive_mean_NHHSSI <- apply(gamma_bias_inactive_array, c(1, 2), mean)
mean(gamma_bias_inactive_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    gamma_bias_inactive_mean_NHHSSI[, k],
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


varphi_bias_active_array <- simplify2array(varphi_bias_active_NHHSSI)
varphi_bias_active_mean_NHHSSI <- apply(varphi_bias_active_array, c(1, 2), mean)
mean(varphi_bias_active_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_active_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


varphi_bias_inactive_array <- simplify2array(varphi_bias_inactive_NHHSSI)
varphi_bias_inactive_mean_NHHSSI <- apply(varphi_bias_inactive_array, c(1, 2), mean)
mean(varphi_bias_inactive_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    varphi_bias_inactive_mean_NHHSSI[, k],
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



delta_bias_active_array <- simplify2array(delta_bias_active_NHHSSI)
delta_bias_active_mean_NHHSSI <- apply(delta_bias_active_array, c(1, 2), mean)
mean(delta_bias_active_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_active_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}



delta_bias_inactive_array <- simplify2array(delta_bias_inactive_NHHSSI)
delta_bias_inactive_mean_NHHSSI <- apply(delta_bias_inactive_array, c(1, 2), mean)
mean(delta_bias_inactive_mean_NHHSSI)
par(mfrow = c(2, 2))
for (k in 1:K) {
  plot(
    delta_bias_inactive_mean_NHHSSI[, k],
    main = paste("Outcome", k),
    pch = 19
  )
}


round(mean(beta_RMSE_NHHSSI)*100, 2)
round(mean(beta_RMSE_active_NHHSSI)*100, 2)
round(mean(beta_RMSE_inactive_NHHSSI)*100, 2)

round(mean(gamma_RMSE_NHHSSI)*100, 2)
round(mean(gamma_RMSE_active_NHHSSI)*100, 2)
round(mean(gamma_RMSE_inactive_NHHSSI)*100, 2)

round(mean(varphi_RMSE_NHHSSI)*100, 2)
round(mean(varphi_RMSE_active_NHHSSI)*100, 2)
round(mean(varphi_RMSE_inactive_NHHSSI)*100, 2)

round(mean(delta_RMSE_NHHSSI)*100, 2)
round(mean(delta_RMSE_active_NHHSSI)*100, 2)
round(mean(delta_RMSE_inactive_NHHSSI)*100, 2)

round(mean(beta_sensitivity_NHHSSI)*100, 2)
round(mean(gamma_sensitivity_NHHSSI)*100, 2)
round(mean(varphi_sensitivity_NHHSSI)*100, 2)

round(mean(beta_specificity_NHHSSI)*100, 2)
round(mean(gamma_specificity_NHHSSI)*100, 2)
round(mean(varphi_specificity_NHHSSI)*100, 2)

round(mean(beta_precision_NHHSSI)*100, 2)
round(mean(gamma_precision_NHHSSI)*100, 2)
round(mean(varphi_precision_NHHSSI)*100, 2)




###################################################

# NHRHS_SI model

# load("practice_data/results_NHHSSISI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")



nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_NHRHS_SI <- c()
beta_cover_per_NHRHS_SI <- c()
gamma_cover_per_NHRHS_SI <- c()
delta_cover_per_NHRHS_SI <- c()
delta_RMSE_NHRHS_SI <- c()
RMSE_all_main_effects_NHRHS_SI <- c()
beta_cover_nonzero_per_NHRHS_SI <- c()
gamma_cover_nonzero_per_NHRHS_SI <- c()
delta_cover_nonzero_per_NHRHS_SI <- c()
alpha_est_mat_NHRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_NHRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_NHRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_NHRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
beta_est_list_NHRHS_SI <- list()
beta_true_list_NHRHS_SI <- list()
beta_lower_list_NHRHS_SI <- list()
beta_upper_list_NHRHS_SI <- list()
gamma_est_list_NHRHS_SI <- list()
gamma_true_list_NHRHS_SI <- list()
gamma_lower_list_NHRHS_SI <- list()
gamma_upper_list_NHRHS_SI <- list()
delta_est_list_NHRHS_SI <- list()
delta_true_list_NHRHS_SI <- list()
delta_lower_list_NHRHS_SI <- list()
delta_upper_list_NHRHS_SI <- list()

main_sensitivity_NHRHS_SI <- c()
main_specificity_NHRHS_SI <- c()
main_precision_NHRHS_SI <- c()

interaction_sensitivity_NHRHS_SI <- c()
interaction_specificity_NHRHS_SI <- c()
interaction_precision_NHRHS_SI <- c()

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_NHRHS_SI[sample, ] <- alpha_true
  alpha_est_mat_NHRHS_SI[sample, ] <- apply(theta_update_NHRHS_SI_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NHRHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NHRHS_SI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NHRHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NHRHS_SI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NHRHS_SI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_NHRHS_SI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NHRHS_SI_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NHRHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NHRHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_NHRHS_SI[[sample]] <- beta_est
  beta_lower_list_NHRHS_SI[[sample]] <- beta_lower
  beta_upper_list_NHRHS_SI[[sample]] <- beta_upper
  beta_cover_per_NHRHS_SI[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_NHRHS_SI[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_NHRHS_SI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NHRHS_SI_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NHRHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NHRHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_NHRHS_SI[[sample]] <- gamma_est
  gamma_lower_list_NHRHS_SI[[sample]] <- gamma_lower
  gamma_upper_list_NHRHS_SI[[sample]] <- gamma_upper
  gamma_cover_per_NHRHS_SI[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_NHRHS_SI[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_NHRHS_SI[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NHRHS_SI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHRHS_SI_list[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NHRHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NHRHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_NHRHS_SI[[sample]] <- delta_est
  delta_lower_list_NHRHS_SI[[sample]] <- delta_lower
  delta_upper_list_NHRHS_SI[[sample]] <- delta_upper
  delta_cover_per_NHRHS_SI[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_NHRHS_SI[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_NHRHS_SI[sample] <- sqrt(mean(delta_MSE))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_main <- c(truth_beta, truth_gamma)
  selected_main <- c(selected_beta, selected_gamma)

  main_TP <- sum(truth_main == 1 & selected_main == 1)
  main_FN <- sum(truth_main == 1 & selected_main == 0)
  main_TN <- sum(truth_main == 0 & selected_main == 0)
  main_FP <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_NHRHS_SI[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_NHRHS_SI[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_NHRHS_SI[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NHRHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_NHRHS_SI[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_NHRHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_NHRHS_SI)*100, 2)
round(mean(delta_RMSE_NHRHS_SI)*100, 2)
round((mean(beta_cover_per_NHRHS_SI) + mean(gamma_cover_per_NHRHS_SI))/2, 2)
round((mean(beta_cover_nonzero_per_NHRHS_SI) + mean(gamma_cover_nonzero_per_NHRHS_SI))/2, 2)
round(mean(delta_cover_per_NHRHS_SI), 2)
round(mean(delta_cover_nonzero_per_NHRHS_SI), 2)




round(mean(main_sensitivity_NHRHS_SI) * 100, 2)
round(mean(interaction_sensitivity_NHRHS_SI) * 100, 2)

round(mean(main_specificity_NHRHS_SI) * 100, 2)
round(mean(interaction_specificity_NHRHS_SI) * 100, 2)

round(mean(main_precision_NHRHS_SI) * 100, 2)
round(mean(interaction_precision_NHRHS_SI) * 100, 2)








###################
# HRHS_SI model

# load("practice_data/theta_update_HRHS_SI_list.rda")

# load("practice_data/theta_update_HRHS_SI_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_HRHS_SI <- c()
beta_cover_per_HRHS_SI <- c()
gamma_cover_per_HRHS_SI <- c()
delta_cover_per_HRHS_SI <- c()
delta_RMSE_HRHS_SI <- c()
RMSE_all_main_effects_HRHS_SI <- c()
beta_cover_nonzero_per_HRHS_SI <- c()
gamma_cover_nonzero_per_HRHS_SI <- c()
delta_cover_nonzero_per_HRHS_SI <- c()
alpha_est_mat_HRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_HRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_HRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_HRHS_SI <- matrix(NA, nrow = nsample, ncol = K)
beta_est_list_HRHS_SI <- list()
beta_true_list_HRHS_SI <- list()
beta_lower_list_HRHS_SI <- list()
beta_upper_list_HRHS_SI <- list()
gamma_est_list_HRHS_SI <- list()
gamma_true_list_HRHS_SI <- list()
gamma_lower_list_HRHS_SI <- list()
gamma_upper_list_HRHS_SI <- list()
delta_est_list_HRHS_SI <- list()
delta_true_list_HRHS_SI <- list()
delta_lower_list_HRHS_SI <- list()
delta_upper_list_HRHS_SI <- list()

main_sensitivity_HRHS_SI <- c()
main_specificity_HRHS_SI <- c()
main_precision_HRHS_SI <- c()

interaction_sensitivity_HRHS_SI <- c()
interaction_specificity_HRHS_SI <- c()
interaction_precision_HRHS_SI <- c()

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_HRHS_SI[sample, ] <- alpha_true
  alpha_est_mat_HRHS_SI[sample, ] <- apply(theta_update_HRHS_SI_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_HRHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_HRHS_SI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_HRHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_HRHS_SI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_HRHS_SI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_HRHS_SI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_HRHS_SI_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_HRHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_HRHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_HRHS_SI[[sample]] <- beta_est
  beta_lower_list_HRHS_SI[[sample]] <- beta_lower
  beta_upper_list_HRHS_SI[[sample]] <- beta_upper
  beta_cover_per_HRHS_SI[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_HRHS_SI[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_HRHS_SI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_HRHS_SI_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_HRHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_HRHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_HRHS_SI[[sample]] <- gamma_est
  gamma_lower_list_HRHS_SI[[sample]] <- gamma_lower
  gamma_upper_list_HRHS_SI[[sample]] <- gamma_upper
  gamma_cover_per_HRHS_SI[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_HRHS_SI[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_HRHS_SI[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_HRHS_SI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HRHS_SI_list[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_HRHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_HRHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_HRHS_SI[[sample]] <- delta_est
  delta_lower_list_HRHS_SI[[sample]] <- delta_lower
  delta_upper_list_HRHS_SI[[sample]] <- delta_upper
  delta_cover_per_HRHS_SI[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_HRHS_SI[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_HRHS_SI[sample] <- sqrt(mean(delta_MSE))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_main <- c(truth_beta, truth_gamma)
  selected_main <- c(selected_beta, selected_gamma)

  main_TP <- sum(truth_main == 1 & selected_main == 1)
  main_FN <- sum(truth_main == 1 & selected_main == 0)
  main_TN <- sum(truth_main == 0 & selected_main == 0)
  main_FP <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_HRHS_SI[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_HRHS_SI[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_HRHS_SI[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_HRHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_HRHS_SI[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_HRHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_HRHS_SI)*100, 2)
round(mean(delta_RMSE_HRHS_SI)*100, 2)
round((mean(beta_cover_per_HRHS_SI) + mean(gamma_cover_per_HRHS_SI))/2, 2)
round((mean(beta_cover_nonzero_per_HRHS_SI) + mean(gamma_cover_nonzero_per_HRHS_SI))/2, 2)
round(mean(delta_cover_per_HRHS_SI), 2)
round(mean(delta_cover_nonzero_per_HRHS_SI), 2)




round(mean(main_sensitivity_HRHS_SI) * 100, 2)
round(mean(interaction_sensitivity_HRHS_SI) * 100, 2)

round(mean(main_specificity_HRHS_SI) * 100, 2)
round(mean(interaction_specificity_HRHS_SI) * 100, 2)

round(mean(main_precision_HRHS_SI) * 100, 2)
round(mean(interaction_precision_HRHS_SI) * 100, 2)








###################################
# HRHS model

# load("practice_data/theta_update_HRHS_list.rda")

# load("practice_data/theta_update_HRHS_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_HRHS <- c()
beta_cover_per_HRHS <- c()
gamma_cover_per_HRHS <- c()
delta_cover_per_HRHS <- c()
delta_RMSE_HRHS <- c()
RMSE_all_main_effects_HRHS <- c()
beta_cover_nonzero_per_HRHS <- c()
gamma_cover_nonzero_per_HRHS <- c()
delta_cover_nonzero_per_HRHS <- c()
alpha_est_mat_HRHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_HRHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_HRHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_HRHS <- matrix(NA, nrow = nsample, ncol = K)
beta_est_list_HRHS <- list()
beta_true_list_HRHS <- list()
beta_lower_list_HRHS <- list()
beta_upper_list_HRHS <- list()
gamma_est_list_HRHS <- list()
gamma_true_list_HRHS <- list()
gamma_lower_list_HRHS <- list()
gamma_upper_list_HRHS <- list()
delta_est_list_HRHS <- list()
delta_true_list_HRHS <- list()
delta_lower_list_HRHS <- list()
delta_upper_list_HRHS <- list()

main_sensitivity_HRHS <- c()
main_specificity_HRHS <- c()
main_precision_HRHS <- c()

interaction_sensitivity_HRHS <- c()
interaction_specificity_HRHS <- c()
interaction_precision_HRHS <- c()

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_HRHS[sample, ] <- alpha_true
  alpha_est_mat_HRHS[sample, ] <- apply(theta_update_HRHS_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_HRHS_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_HRHS[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_HRHS_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_HRHS[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_HRHS[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_HRHS[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_HRHS_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_HRHS_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_HRHS_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_HRHS[[sample]] <- beta_est
  beta_lower_list_HRHS[[sample]] <- beta_lower
  beta_upper_list_HRHS[[sample]] <- beta_upper
  beta_cover_per_HRHS[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_HRHS[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_HRHS[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_HRHS_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_HRHS_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_HRHS_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_HRHS[[sample]] <- gamma_est
  gamma_lower_list_HRHS[[sample]] <- gamma_lower
  gamma_upper_list_HRHS[[sample]] <- gamma_upper
  gamma_cover_per_HRHS[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_HRHS[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_HRHS[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_HRHS[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HRHS_list[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_HRHS_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_HRHS_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_HRHS[[sample]] <- delta_est
  delta_lower_list_HRHS[[sample]] <- delta_lower
  delta_upper_list_HRHS[[sample]] <- delta_upper
  delta_cover_per_HRHS[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_HRHS[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_HRHS[sample] <- sqrt(mean(delta_MSE))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_main <- c(truth_beta, truth_gamma)
  selected_main <- c(selected_beta, selected_gamma)

  main_TP <- sum(truth_main == 1 & selected_main == 1)
  main_FN <- sum(truth_main == 1 & selected_main == 0)
  main_TN <- sum(truth_main == 0 & selected_main == 0)
  main_FP <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_HRHS[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_HRHS[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_HRHS[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_HRHS[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_HRHS[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_HRHS[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_HRHS)*100, 2)
round(mean(delta_RMSE_HRHS)*100, 2)
round((mean(beta_cover_per_HRHS) + mean(gamma_cover_per_HRHS))/2, 2)
round((mean(beta_cover_nonzero_per_HRHS) + mean(gamma_cover_nonzero_per_HRHS))/2, 2)
round(mean(delta_cover_per_HRHS), 2)
round(mean(delta_cover_nonzero_per_HRHS), 2)




round(mean(main_sensitivity_HRHS) * 100, 2)
round(mean(interaction_sensitivity_HRHS) * 100, 2)

round(mean(main_specificity_HRHS) * 100, 2)
round(mean(interaction_specificity_HRHS) * 100, 2)

round(mean(main_precision_HRHS) * 100, 2)
round(mean(interaction_precision_HRHS) * 100, 2)








#############################
# NHRHS model

# load("practice_data/theta_update_NHRHS_list.rda")

# load("practice_data/theta_update_NHRHS_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_NHRHS <- c()
beta_cover_per_NHRHS <- c()
gamma_cover_per_NHRHS <- c()
delta_cover_per_NHRHS <- c()
delta_RMSE_NHRHS <- c()
RMSE_all_main_effects_NHRHS <- c()
beta_cover_nonzero_per_NHRHS <- c()
gamma_cover_nonzero_per_NHRHS <- c()
delta_cover_nonzero_per_NHRHS <- c()
alpha_est_mat_NHRHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_NHRHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_NHRHS <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_NHRHS <- matrix(NA, nrow = nsample, ncol = K)
beta_est_list_NHRHS <- list()
beta_true_list_NHRHS <- list()
beta_lower_list_NHRHS <- list()
beta_upper_list_NHRHS <- list()
gamma_est_list_NHRHS <- list()
gamma_true_list_NHRHS <- list()
gamma_lower_list_NHRHS <- list()
gamma_upper_list_NHRHS <- list()
delta_est_list_NHRHS <- list()
delta_true_list_NHRHS <- list()
delta_lower_list_NHRHS <- list()
delta_upper_list_NHRHS <- list()

main_sensitivity_NHRHS <- c()
main_specificity_NHRHS <- c()
main_precision_NHRHS <- c()

interaction_sensitivity_NHRHS <- c()
interaction_specificity_NHRHS <- c()
interaction_precision_NHRHS <- c()

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_NHRHS[sample, ] <- alpha_true
  alpha_est_mat_NHRHS[sample, ] <- apply(theta_update_NHRHS_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NHRHS_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NHRHS[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NHRHS_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NHRHS[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NHRHS[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_NHRHS[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NHRHS_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NHRHS_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NHRHS_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_NHRHS[[sample]] <- beta_est
  beta_lower_list_NHRHS[[sample]] <- beta_lower
  beta_upper_list_NHRHS[[sample]] <- beta_upper
  beta_cover_per_NHRHS[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_NHRHS[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_NHRHS[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NHRHS_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NHRHS_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NHRHS_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_NHRHS[[sample]] <- gamma_est
  gamma_lower_list_NHRHS[[sample]] <- gamma_lower
  gamma_upper_list_NHRHS[[sample]] <- gamma_upper
  gamma_cover_per_NHRHS[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_NHRHS[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_NHRHS[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NHRHS[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHRHS_list[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NHRHS_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NHRHS_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_NHRHS[[sample]] <- delta_est
  delta_lower_list_NHRHS[[sample]] <- delta_lower
  delta_upper_list_NHRHS[[sample]] <- delta_upper
  delta_cover_per_NHRHS[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_NHRHS[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_NHRHS[sample] <- sqrt(mean(delta_MSE))

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_main <- c(truth_beta, truth_gamma)
  selected_main <- c(selected_beta, selected_gamma)

  main_TP <- sum(truth_main == 1 & selected_main == 1)
  main_FN <- sum(truth_main == 1 & selected_main == 0)
  main_TN <- sum(truth_main == 0 & selected_main == 0)
  main_FP <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_NHRHS[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_NHRHS[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_NHRHS[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NHRHS[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_NHRHS[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_NHRHS[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_NHRHS)*100, 2)
round(mean(delta_RMSE_NHRHS)*100, 2)
round((mean(beta_cover_per_NHRHS) + mean(gamma_cover_per_NHRHS))/2, 2)
round((mean(beta_cover_nonzero_per_NHRHS) + mean(gamma_cover_nonzero_per_NHRHS))/2, 2)
round(mean(delta_cover_per_NHRHS), 2)
round(mean(delta_cover_nonzero_per_NHRHS), 2)




round(mean(main_sensitivity_NHRHS) * 100, 2)
round(mean(interaction_sensitivity_NHRHS) * 100, 2)

round(mean(main_specificity_NHRHS) * 100, 2)
round(mean(interaction_specificity_NHRHS) * 100, 2)

round(mean(main_precision_NHRHS) * 100, 2)
round(mean(interaction_precision_NHRHS) * 100, 2)







####################################
# NP model

# load("practice_data/theta_update_NP_list.rda")

# load("practice_data/theta_update_NP_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_NP <- c()
beta_cover_per_NP <- c()
gamma_cover_per_NP <- c()
delta_cover_per_NP <- c()
delta_RMSE_NP <- c()
RMSE_all_main_effects_NP <- c()
beta_cover_nonzero_per_NP <- c()
gamma_cover_nonzero_per_NP <- c()
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

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


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

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_NP[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NP[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NP_list[[sample]][, jm, ], 2, mean)
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

  # For sensitivity, specificity, precision
  truth_beta <- as.vector(beta_true != 0)
  selected_beta <- as.vector((beta_lower > 0) | (beta_upper < 0))
  truth_gamma <- as.vector(gamma_true != 0)
  selected_gamma <- as.vector((gamma_lower > 0) | (gamma_upper < 0))
  truth_main <- c(truth_beta, truth_gamma)
  selected_main <- c(selected_beta, selected_gamma)

  main_TP <- sum(truth_main == 1 & selected_main == 1)
  main_FN <- sum(truth_main == 1 & selected_main == 0)
  main_TN <- sum(truth_main == 0 & selected_main == 0)
  main_FP <- sum(truth_main == 0 & selected_main == 1)

  main_sensitivity_NP[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_NP[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_NP[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NP[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_NP[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_NP[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_NP)*100, 2)
round(mean(delta_RMSE_NP)*100, 2)
round((mean(beta_cover_per_NP) + mean(gamma_cover_per_NP))/2, 2)
round((mean(beta_cover_nonzero_per_NP) + mean(gamma_cover_nonzero_per_NP))/2, 2)
round(mean(delta_cover_per_NP), 2)
round(mean(delta_cover_nonzero_per_NP), 2)




round(mean(main_sensitivity_NP) * 100, 2)
round(mean(interaction_sensitivity_NP) * 100, 2)

round(mean(main_specificity_NP) * 100, 2)
round(mean(interaction_specificity_NP) * 100, 2)

round(mean(main_precision_NP) * 100, 2)
round(mean(interaction_precision_NP) * 100, 2)




######################################################################################

# Credible interval plot
beta_lower_list_HHSSI
beta_upper_list_HHSSI

beta_lower_list_NHHSSI
beta_upper_list_NHHSSI

dim(beta_lower_list_HHSSI[[1]])

beta_est_list_HHSSI
beta_est_list_NHHSSI

dim(beta_est_list_HHSSI[[1]])



# Initialize matrices: rows = beta coefficients, columns = outcomes
beta_est_avg_HHSSI   <- matrix(NA, J, K)
beta_lower_avg_HHSSI <- matrix(NA, J, K)
beta_upper_avg_HHSSI <- matrix(NA, J, K)

beta_est_avg_NHHSSI   <- matrix(NA, J, K)
beta_lower_avg_NHHSSI <- matrix(NA, J, K)
beta_upper_avg_NHHSSI <- matrix(NA, J, K)

beta_true_avg  <- matrix(NA, J, K)

for (j in 1:J) {
  for (k in 1:K) {

    beta_est_avg_HHSSI[j, k] <- mean(
      sapply(1:nsample, function(i)
        beta_est_list_HHSSI[[i]][j, k]),
      na.rm = TRUE
    )

    beta_lower_avg_HHSSI[j, k] <- mean(
      sapply(1:nsample, function(i)
        beta_lower_list_HHSSI[[i]][j, k]),
      na.rm = TRUE
    )

    beta_upper_avg_HHSSI[j, k] <- mean(
      sapply(1:nsample, function(i)
        beta_upper_list_HHSSI[[i]][j, k]),
      na.rm = TRUE
    )

    beta_est_avg_NHHSSI[j, k] <- mean(
      sapply(1:nsample, function(i)
        beta_est_list_NHHSSI[[i]][j, k]),
      na.rm = TRUE
    )

    beta_lower_avg_NHHSSI[j, k] <- mean(
      sapply(1:nsample, function(i)
        beta_lower_list_NHHSSI[[i]][j, k]),
      na.rm = TRUE
    )

    beta_upper_avg_NHHSSI[j, k] <- mean(
      sapply(1:nsample, function(i)
        beta_upper_list_NHHSSI[[i]][j, k]),
      na.rm = TRUE
    )

    beta_true_avg[j, k] <- mean(
      sapply(1:nsample, function(i)
        sim_data_list[[i]]$beta_true[j, k]),
      na.rm = TRUE
    )
  }
}







par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))

for (k in 1:K) {

  # Average HHSSI estimates and bounds
  est_h <- beta_est_avg_HHSSI[, k]
  low_h <- beta_lower_avg_HHSSI[, k]
  up_h  <- beta_upper_avg_HHSSI[, k]

  # Average NHHSSI estimates and bounds
  est_n <- beta_est_avg_NHHSSI[, k]
  low_n <- beta_lower_avg_NHHSSI[, k]
  up_n  <- beta_upper_avg_NHHSSI[, k]

  ylim <- range(low_h, up_h, low_n, up_n, 0,
                na.rm = TRUE)

  # HHSSI estimates
  plot(
    1:J, est_h,
    type = "p",
    pch = 16,
    col = "red",
    ylim = ylim,
    xlim = c(0.5, J + 0.5),
    xaxt = "n",
    xlab = expression(beta ~ "coefficient index"),
    ylab = "Average estimate and interval",
    main = paste("Outcome", k)
  )

  axis(1, at = seq(0, J, by = 5))

  # HHSSI intervals
  segments(
    x0 = (1:J) - 0.15, y0 = low_h,
    x1 = (1:J) - 0.15, y1 = up_h,
    col = "red", lwd = 1.5
  )

  # NHHSSI estimates
  points(1:J, est_n, pch = 17, col = "blue")

  # NHHSSI intervals
  segments(
    x0 = (1:J) + 0.15, y0 = low_n,
    x1 = (1:J) + 0.15, y1 = up_n,
    col = "blue", lwd = 1.5
  )

  # Zero reference line
  abline(h = 0, lty = 2, col = "gray40")

  legend(
    "topright",
    legend = c("HHSSI", "NHHSSI"),
    col = c("red", "blue"),
    pch = c(16, 17),
    bty = "n",
    cex = 0.8
  )
}








####################################

# Initialize matrices: rows = gamma coefficients, columns = outcomes
gamma_est_avg_HHSSI   <- matrix(NA, M, K)
gamma_lower_avg_HHSSI <- matrix(NA, M, K)
gamma_upper_avg_HHSSI <- matrix(NA, M, K)

gamma_est_avg_NHHSSI   <- matrix(NA, M, K)
gamma_lower_avg_NHHSSI <- matrix(NA, M, K)
gamma_upper_avg_NHHSSI <- matrix(NA, M, K)

gamma_true_avg  <- matrix(NA, M, K)

for (m in 1:M) {
  for (k in 1:K) {

    gamma_est_avg_HHSSI[m, k] <- mean(
      sapply(1:nsample, function(i)
        gamma_est_list_HHSSI[[i]][m, k]),
      na.rm = TRUE
    )

    gamma_lower_avg_HHSSI[m, k] <- mean(
      sapply(1:nsample, function(i)
        gamma_lower_list_HHSSI[[i]][m, k]),
      na.rm = TRUE
    )

    gamma_upper_avg_HHSSI[m, k] <- mean(
      sapply(1:nsample, function(i)
        gamma_upper_list_HHSSI[[i]][m, k]),
      na.rm = TRUE
    )

    gamma_est_avg_NHHSSI[m, k] <- mean(
      sapply(1:nsample, function(i)
        gamma_est_list_NHHSSI[[i]][m, k]),
      na.rm = TRUE
    )

    gamma_lower_avg_NHHSSI[m, k] <- mean(
      sapply(1:nsample, function(i)
        gamma_lower_list_NHHSSI[[i]][m, k]),
      na.rm = TRUE
    )

    gamma_upper_avg_NHHSSI[m, k] <- mean(
      sapply(1:nsample, function(i)
        gamma_upper_list_NHHSSI[[i]][m, k]),
      na.rm = TRUE
    )


    gamma_true_avg[m, k] <- mean(
      sapply(1:nsample, function(i)
        sim_data_list[[i]]$gamma_true[m, k]),
      na.rm = TRUE
    )
  }
}







par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))

for (k in 1:K) {

  # Average HHSSI estimates and bounds
  est_h <- gamma_est_avg_HHSSI[, k]
  low_h <- gamma_lower_avg_HHSSI[, k]
  up_h  <- gamma_upper_avg_HHSSI[, k]

  # Average NHHSSI estimates and bounds
  est_n <- gamma_est_avg_NHHSSI[, k]
  low_n <- gamma_lower_avg_NHHSSI[, k]
  up_n  <- gamma_upper_avg_NHHSSI[, k]

  ylim <- range(low_h, up_h, low_n, up_n, 0,
                na.rm = TRUE)

  # HHSSI estimates
  plot(
    1:M, est_h,
    type = "p",
    pch = 16,
    col = "red",
    ylim = ylim,
    xlim = c(0.5, M + 0.5),
    xaxt = "n",
    xlab = expression(gamma ~ "coefficient index"),
    ylab = "Average estimate and interval",
    main = paste("Outcome", k)
  )
  axis(1, at = seq(0, M, by = 5))

  # HHSSI intervals
  segments(x0 = (1:M) - 0.15, y0 = low_h, x1 = (1:M) - 0.15, y1 = up_h,
    col = "red", lwd = 1.5)

  # NHHSSI estimates
  points(1:M, est_n, pch = 17, col = "blue")
  # NHHSSI intervals
  segments(x0 = (1:M) + 0.15, y0 = low_n, x1 = (1:M) + 0.15, y1 = up_n,
    col = "blue", lwd = 1.5)



  # Zero reference line
  abline(h = 0, lty = 2, col = "gray40")

  legend(
    "topright",
    legend = c("HHSSI", "NHHSSI"),
    col = c("red", "blue"),
    pch = c(16, 17),
    bty = "n",
    cex = 0.8
  )
}



####################################

# Initialize matrices: rows = delta coefficients, columns = outcomes
delta_est_avg_HHSSI   <- matrix(NA, J*M, K)
delta_lower_avg_HHSSI <- matrix(NA, J*M, K)
delta_upper_avg_HHSSI <- matrix(NA, J*M, K)

delta_est_avg_NHHSSI   <- matrix(NA, J*M, K)
delta_lower_avg_NHHSSI <- matrix(NA, J*M, K)
delta_upper_avg_NHHSSI <- matrix(NA, J*M, K)

delta_true_avg  <- matrix(NA, J*M, K)

for (j in 1:J) {
for (m in 1:M) {
  for (k in 1:K) {

    jm <- (j - 1) * M + m

    delta_est_avg_HHSSI[jm, k] <- mean(
      sapply(1:nsample, function(i)
        delta_est_list_HHSSI[[i]][jm, k]),
      na.rm = TRUE
    )

    delta_lower_avg_HHSSI[jm, k] <- mean(
      sapply(1:nsample, function(i)
        delta_lower_list_HHSSI[[i]][jm, k]),
      na.rm = TRUE
    )

    delta_upper_avg_HHSSI[jm, k] <- mean(
      sapply(1:nsample, function(i)
        delta_upper_list_HHSSI[[i]][jm, k]),
      na.rm = TRUE
    )

    delta_est_avg_NHHSSI[jm, k] <- mean(
      sapply(1:nsample, function(i)
        delta_est_list_NHHSSI[[i]][jm, k]),
      na.rm = TRUE
    )

    delta_lower_avg_NHHSSI[jm, k] <- mean(
      sapply(1:nsample, function(i)
        delta_lower_list_NHHSSI[[i]][jm, k]),
      na.rm = TRUE
    )

    delta_upper_avg_NHHSSI[jm, k] <- mean(
      sapply(1:nsample, function(i)
        delta_upper_list_NHHSSI[[i]][jm, k]),
      na.rm = TRUE
    )


    delta_true_avg[jm, k] <- mean(
      sapply(1:nsample, function(i)
        sim_data_list[[i]]$delta_true[jm, k]),
      na.rm = TRUE
    )
  }
 }
}







par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))

for (k in 1:K) {

  est_h <- delta_est_avg_HHSSI[, k]
  low_h <- delta_lower_avg_HHSSI[, k]
  up_h  <- delta_upper_avg_HHSSI[, k]

  est_n <- delta_est_avg_NHHSSI[, k]
  low_n <- delta_lower_avg_NHHSSI[, k]
  up_n  <- delta_upper_avg_NHHSSI[, k]

  idx <- seq_len(J * M)

  ylim <- range(low_h, up_h, low_n, up_n, 0, na.rm = TRUE)

  # Plot HHSSI estimates
  plot(
    idx, est_h,
    type = "p",
    pch = 16,
    cex = 0.5,
    col = "red",
    ylim = ylim,
    xlim = c(0.5, J * M + 0.5),
    xaxt = "n",
    xlab = expression(delta ~ "coefficient index"),
    ylab = "Average estimate and 95% CrI",
    main = paste("Outcome", k)
  )

  # Display fewer x-axis labels for readability
  axis(1, at = seq(0, J * M, by = 50))

  # HHSSI credible intervals
  segments(
    x0 = idx - 0.15, y0 = low_h,
    x1 = idx - 0.15, y1 = up_h,
    col = "red", lwd = 0.5
  )

  # NHHSSI estimates
  points(idx, est_n, pch = 17, cex = 0.5, col = "blue")

  # NHHSSI credible intervals
  segments(
    x0 = idx + 0.15, y0 = low_n,
    x1 = idx + 0.15, y1 = up_n,
    col = "blue", lwd = 0.5
  )

  # Zero reference line
  abline(h = 0, lty = 2, col = "gray40")

  legend(
    "topright",
    legend = c("HHSSI", "NHHSSI"),
    col = c("red", "blue"),
    pch = c(16, 17),
    bty = "n",
    cex = 0.8
  )
}





# coverage
# Coverage for each outcome
coverage_delta_avg <- data.frame(
  Outcome = 1:K,
  HHSSI_Covered = NA_integer_,
  HHSSI_Total = NA_integer_,
  HHSSI_Percentage = NA_real_,
  NHHSSI_Covered = NA_integer_,
  NHHSSI_Total = NA_integer_,
  NHHSSI_Percentage = NA_real_
)

for (k in 1:K) {

  # HHSSI
  valid_h <- complete.cases(
    delta_lower_avg_HHSSI[, k],
    delta_upper_avg_HHSSI[, k],
    delta_true_avg[, k]
  )

  covered_h <- (
    delta_lower_avg_HHSSI[valid_h, k] <= delta_true_avg[valid_h, k] &
      delta_upper_avg_HHSSI[valid_h, k] >= delta_true_avg[valid_h, k]
  )

  coverage_delta_avg$HHSSI_Covered[k] <- sum(covered_h)
  coverage_delta_avg$HHSSI_Total[k] <- sum(valid_h)
  coverage_delta_avg$HHSSI_Percentage[k] <-
    mean(covered_h) * 100

  # NHHSSI
  valid_n <- complete.cases(
    delta_lower_avg_NHHSSI[, k],
    delta_upper_avg_NHHSSI[, k],
    delta_true_avg[, k]
  )

  covered_n <- (
    delta_lower_avg_NHHSSI[valid_n, k] <= delta_true_avg[valid_n, k] &
      delta_upper_avg_NHHSSI[valid_n, k] >= delta_true_avg[valid_n, k]
  )

  coverage_delta_avg$NHHSSI_Covered[k] <- sum(covered_n)
  coverage_delta_avg$NHHSSI_Total[k] <- sum(valid_n)
  coverage_delta_avg$NHHSSI_Percentage[k] <-
    mean(covered_n) * 100
}
print(coverage_delta_avg)




###################################################################

# rmse heatmap with simulated predictors
rmse_hhssi <- colMeans(delta_RMSE_indices_mat_HHSSI)
rmse_nhhssi <- colMeans(delta_RMSE_indices_mat_NHHSSI)

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$RMSE_HHSSI <- rmse_hhssi
df$RMSE_NHHSSI <- rmse_nhhssi

rmse_range <- range(c(rmse_hhssi, rmse_nhhssi))
mid_value <- mean(rmse_range)

df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "RMSE") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = rmse_range,
    name = "RMSE"
  ) +

  labs(
    title = "RMSE Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of simulated predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))




# log transformed
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$RMSE_HHSSI <- log(rmse_hhssi)
df$RMSE_NHHSSI <- log(rmse_nhhssi)

rmse_range <- range(log(c(rmse_hhssi, rmse_nhhssi)))
mid_value <- mean(rmse_range)

df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "log RMSE") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = rmse_range,
    name = "log RMSE"
  ) +

  labs(
    title = "Log of RMSE Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of simulated predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))



###################################################################

# coverage heatmap with simulated predictors
cover_hhssi <- colMeans(delta_cover_indices_mat_HHSSI)
cover_nhhssi <- colMeans(delta_cover_indices_mat_NHHSSI)

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$cover_HHSSI <- cover_hhssi
df$cover_NHHSSI <- cover_nhhssi

cover_range <- range(c(cover_hhssi, cover_nhhssi))
mid_value <- mean(cover_range)

df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = cover)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "cover") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = cover_range,
    name = "cover"
  ) +

  labs(
    title = "Coverage Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of simulated predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))




########################################################################

# variance heatmap with simulated predictors
var_hhssi <- colMeans(delta_var_indices_mat_HHSSI)
var_nhhssi <- colMeans(delta_var_indices_mat_NHHSSI)

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$var_HHSSI <- var_hhssi
df$var_NHHSSI <- var_nhhssi

var_range <- range(c(var_hhssi, var_nhhssi))
mid_value <- mean(var_range)

df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            var = df$var_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, var = df$var_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = var)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "var") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = var_range,
    name = "var"
  ) +

  labs(
    title = "Variance Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of simulated predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))




# log transformed
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$var_HHSSI <- log(var_hhssi)
df$var_NHHSSI <- log(var_nhhssi)

var_range <- range(log(c(var_hhssi, var_nhhssi)))
mid_value <- mean(var_range)

df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = var)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "log var") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = var_range,
    name = "log var"
  ) +

  labs(
    title = "Log of Variance Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of simulated predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))



###################################################################

# rmse heatmap with real predictors
rmse_hhssi <- colMeans(delta_RMSE_indices_mat_HHSSI)
rmse_nhhssi <- colMeans(delta_RMSE_indices_mat_NHHSSI)

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$RMSE_HHSSI <- rmse_hhssi
df$RMSE_NHHSSI <- rmse_nhhssi

rmse_range <- range(c(rmse_hhssi, rmse_nhhssi))
mid_value <- mean(rmse_range)

df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "RMSE") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = rmse_range,
    name = "RMSE"
  ) +

  labs(
    title = "RMSE Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of real predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))




# log transformed
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$RMSE_HHSSI <- log(rmse_hhssi)
df$RMSE_NHHSSI <- log(rmse_nhhssi)

rmse_range <- range(log(c(rmse_hhssi, rmse_nhhssi)))
mid_value <- mean(rmse_range)

df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "log RMSE") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = rmse_range,
    name = "log RMSE"
  ) +

  labs(
    title = "Log of RMSE Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of real predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))



###################################################################

# coverage heatmap with real predictors
cover_hhssi <- colMeans(delta_cover_indices_mat_HHSSI)
cover_nhhssi <- colMeans(delta_cover_indices_mat_NHHSSI)

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$cover_HHSSI <- cover_hhssi
df$cover_NHHSSI <- cover_nhhssi

cover_range <- range(c(cover_hhssi, cover_nhhssi))
mid_value <- mean(cover_range)

df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            cover = df$cover_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = cover)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "cover") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = cover_range,
    name = "cover"
  ) +

  labs(
    title = "Coverage Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of real predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))



########################################################################

# variance heatmap with real predictors
var_hhssi <- colMeans(delta_var_indices_mat_HHSSI)
var_nhhssi <- colMeans(delta_var_indices_mat_NHHSSI)

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$var_HHSSI <- var_hhssi
df$var_NHHSSI <- var_nhhssi

var_range <- range(c(var_hhssi, var_nhhssi))
mid_value <- mean(var_range)

df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            var = df$var_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, var = df$var_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = var)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "var") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = var_range,
    name = "var"
  ) +

  labs(
    title = "Variance Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of real predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))




# log transformed
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$var_HHSSI <- log(var_hhssi)
df$var_NHHSSI <- log(var_nhhssi)

var_range <- range(log(c(var_hhssi, var_nhhssi)))
mid_value <- mean(var_range)

df_long <- rbind(data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure,
                            modifier = df$modifier,
                            var = df$var_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = var)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  # scale_fill_gradient(
  #   low = "white",
  #   high = "red",
  #   name = "log var") +

  scale_fill_gradient2(
    low = "white",
    mid = "green",
    high = "red",
    midpoint = mid_value,
    limits = var_range,
    name = "log var"
  ) +

  labs(
    title = "Log of Variance Heatmap of Modifier–Exposure Interaction Estimates (with 100 datasets of real predictors)",
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
sim_data_list[[1]]$delta_true


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

# Convergence of NHHSSI for weak heredity
dim(theta_update_HHSSI_list[[1]][, 63:810, ])

sim_data_list[[1]]$delta_true





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

# Identify active interactions
df_long$Active <- with(
  df_long,
  modifier %in% 1:2 & exposure %in% 1:3
)

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

# Highlight active interactions
df_long$Active <- with(
  df_long,
  modifier %in% 1:2 & exposure %in% 1:3
)

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







