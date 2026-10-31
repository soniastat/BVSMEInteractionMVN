
rm(list = ls())
gc()

# load the simulated data
# load("results/sample_data_list.rda")

# load("results/sample_data_list_Sigma.rda")

#####################################################

# HHS_SI model

# load("results/theta_update_HHS_SI_list.rda")

# load("results/theta_update_HHS_SI_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_HHS_SI <- c()
beta_cover_per_HHS_SI <- c()
gamma_cover_per_HHS_SI <- c()
delta_cover_per_HHS_SI <- c()
delta_RMSE_HHS_SI <- c()
RMSE_all_main_effects_HHS_SI <- c()
beta_cover_nonzero_per_HHS_SI <- c()
gamma_cover_nonzero_per_HHS_SI <- c()
delta_cover_nonzero_per_HHS_SI <- c()
alpha_est_mat_HHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_HHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_HHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_HHS_SI <- matrix(NA, nrow = nsample, ncol = K)
beta_est_list_HHS_SI <- list()
beta_true_list_HHS_SI <- list()
beta_lower_list_HHS_SI <- list()
beta_upper_list_HHS_SI <- list()
gamma_est_list_HHS_SI <- list()
gamma_true_list_HHS_SI <- list()
gamma_lower_list_HHS_SI <- list()
gamma_upper_list_HHS_SI <- list()
delta_est_list_HHS_SI <- list()
delta_true_list_HHS_SI <- list()
delta_lower_list_HHS_SI <- list()
delta_upper_list_HHS_SI <- list()

main_sensitivity_HHS_SI <- c()
main_specificity_HHS_SI <- c()
main_precision_HHS_SI <- c()

interaction_sensitivity_HHS_SI <- c()
interaction_specificity_HHS_SI <- c()
interaction_precision_HHS_SI <- c()

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$E
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_HHS_SI[sample, ] <- alpha_true
  alpha_est_mat_HHS_SI[sample, ] <- apply(theta_update_HHS_SI_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_HHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_HHS_SI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_HHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_HHS_SI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_HHS_SI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_HHS_SI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_HHS_SI_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_HHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_HHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_HHS_SI[[sample]] <- beta_est
  beta_lower_list_HHS_SI[[sample]] <- beta_lower
  beta_upper_list_HHS_SI[[sample]] <- beta_upper
  beta_cover_per_HHS_SI[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_HHS_SI[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_HHS_SI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_HHS_SI_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_HHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_HHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_HHS_SI[[sample]] <- gamma_est
  gamma_lower_list_HHS_SI[[sample]] <- gamma_lower
  gamma_upper_list_HHS_SI[[sample]] <- gamma_upper
  gamma_cover_per_HHS_SI[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_HHS_SI[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_HHS_SI[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_HHS_SI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HHS_SI_list[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_HHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_HHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_HHS_SI[[sample]] <- delta_est
  delta_lower_list_HHS_SI[[sample]] <- delta_lower
  delta_upper_list_HHS_SI[[sample]] <- delta_upper
  delta_cover_per_HHS_SI[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_HHS_SI[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_HHS_SI[sample] <- sqrt(mean(delta_MSE))

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

  main_sensitivity_HHS_SI[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_HHS_SI[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_HHS_SI[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_HHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_HHS_SI[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_HHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}



round(mean(RMSE_all_main_effects_HHS_SI)*100, 2)
round(mean(delta_RMSE_HHS_SI)*100, 2)
round((mean(beta_cover_per_HHS_SI) + mean(gamma_cover_per_HHS_SI ))/2, 2)
round((mean(beta_cover_nonzero_per_HHS_SI) + mean(gamma_cover_nonzero_per_HHS_SI ))/2, 2)
round(mean(delta_cover_per_HHS_SI), 2)
round(mean(delta_cover_nonzero_per_HHS_SI), 2)



round(mean(main_sensitivity_HHS_SI) * 100, 2)
round(mean(interaction_sensitivity_HHS_SI) * 100, 2)

round(mean(main_specificity_HHS_SI) * 100, 2)
round(mean(interaction_specificity_HHS_SI) * 100, 2)

round(mean(main_precision_HHS_SI) * 100, 2)
round(mean(interaction_precision_HHS_SI) * 100, 2)







###################
# NHHS_SI model

# load("results/theta_update_NHHS_SI_list.rda")

# load("results/theta_update_NHHS_SI_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_NHHS_SI <- c()
beta_cover_per_NHHS_SI <- c()
gamma_cover_per_NHHS_SI <- c()
delta_cover_per_NHHS_SI <- c()
delta_RMSE_NHHS_SI <- c()
RMSE_all_main_effects_NHHS_SI <- c()
beta_cover_nonzero_per_NHHS_SI <- c()
gamma_cover_nonzero_per_NHHS_SI <- c()
delta_cover_nonzero_per_NHHS_SI <- c()
alpha_est_mat_NHHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_NHHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_NHHS_SI <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_NHHS_SI <- matrix(NA, nrow = nsample, ncol = K)
beta_est_list_NHHS_SI <- list()
beta_true_list_NHHS_SI <- list()
beta_lower_list_NHHS_SI <- list()
beta_upper_list_NHHS_SI <- list()
gamma_est_list_NHHS_SI <- list()
gamma_true_list_NHHS_SI <- list()
gamma_lower_list_NHHS_SI <- list()
gamma_upper_list_NHHS_SI <- list()
delta_est_list_NHHS_SI <- list()
delta_true_list_NHHS_SI <- list()
delta_lower_list_NHHS_SI <- list()
delta_upper_list_NHHS_SI <- list()

main_sensitivity_NHHS_SI <- c()
main_specificity_NHHS_SI <- c()
main_precision_NHHS_SI <- c()

interaction_sensitivity_NHHS_SI <- c()
interaction_specificity_NHHS_SI <- c()
interaction_precision_NHHS_SI <- c()

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$E
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list[[sample]]$alpha_true
  alpha_true_mat_NHHS_SI[sample, ] <- alpha_true
  alpha_est_mat_NHHS_SI[sample, ] <- apply(theta_update_NHHS_SI_list[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NHHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NHHS_SI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NHHS_SI_list[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NHHS_SI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NHHS_SI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list[[sample]]$beta_true
  beta_true_list_NHHS_SI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NHHS_SI_list[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NHHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NHHS_SI_list[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_NHHS_SI[[sample]] <- beta_est
  beta_lower_list_NHHS_SI[[sample]] <- beta_lower
  beta_upper_list_NHHS_SI[[sample]] <- beta_upper
  beta_cover_per_NHHS_SI[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_NHHS_SI[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list[[sample]]$gamma_true
  gamma_true_list_NHHS_SI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NHHS_SI_list[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NHHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NHHS_SI_list[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_NHHS_SI[[sample]] <- gamma_est
  gamma_lower_list_NHHS_SI[[sample]] <- gamma_lower
  gamma_upper_list_NHHS_SI[[sample]] <- gamma_upper
  gamma_cover_per_NHHS_SI[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_NHHS_SI[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_NHHS_SI[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NHHS_SI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHHS_SI_list[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NHHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NHHS_SI_list[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_NHHS_SI[[sample]] <- delta_est
  delta_lower_list_NHHS_SI[[sample]] <- delta_lower
  delta_upper_list_NHHS_SI[[sample]] <- delta_upper
  delta_cover_per_NHHS_SI[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_NHHS_SI[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_NHHS_SI[sample] <- sqrt(mean(delta_MSE))

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

  main_sensitivity_NHHS_SI[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_NHHS_SI[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_NHHS_SI[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NHHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_NHHS_SI[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_NHHS_SI[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_NHHS_SI)*100, 2)
round(mean(delta_RMSE_NHHS_SI)*100, 2)
round((mean(beta_cover_per_NHHS_SI) + mean(gamma_cover_per_NHHS_SI ))/2, 2)
round((mean(beta_cover_nonzero_per_NHHS_SI) + mean(gamma_cover_nonzero_per_NHHS_SI ))/2, 2)
round(mean(delta_cover_per_NHHS_SI), 2)
round(mean(delta_cover_nonzero_per_NHHS_SI), 2)




round(mean(main_sensitivity_NHHS_SI) * 100, 2)
round(mean(interaction_sensitivity_NHHS_SI) * 100, 2)

round(mean(main_specificity_NHHS_SI) * 100, 2)
round(mean(interaction_specificity_NHHS_SI) * 100, 2)

round(mean(main_precision_NHHS_SI) * 100, 2)
round(mean(interaction_precision_NHHS_SI) * 100, 2)







####################
# NHRHS_SI model

# load("results/theta_update_NHRHS_SI_list.rda")

# load("results/theta_update_NHRHS_SI_list_Sigma.rda")

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
  u <- sim_data_list[[sample]]$E
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

# load("results/theta_update_HRHS_SI_list.rda")

# load("results/theta_update_HRHS_SI_list_Sigma.rda")

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
  u <- sim_data_list[[sample]]$E
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






#####################
# HHS model

# load("results/theta_update_HHS_list.rda")

# load("results/theta_update_HHS_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_HHS <- c()
beta_cover_per_HHS <- c()
gamma_cover_per_HHS <- c()
delta_cover_per_HHS <- c()
delta_RMSE_HHS <- c()
RMSE_all_main_effects_HHS <- c()
beta_cover_nonzero_per_HHS <- c()
gamma_cover_nonzero_per_HHS <- c()
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

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$E
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


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

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_HHS[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_HHS[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HHS_list[[sample]][, jm, ], 2, mean)
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

  main_sensitivity_HHS[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_HHS[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_HHS[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_HHS[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_HHS[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_HHS[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_HHS)*100, 2)
round(mean(delta_RMSE_HHS)*100, 2)
round((mean(beta_cover_per_HHS) + mean(gamma_cover_per_HHS))/2, 2)
round((mean(beta_cover_nonzero_per_HHS) + mean(gamma_cover_nonzero_per_HHS))/2, 2)
round(mean(delta_cover_per_HHS), 2)
round(mean(delta_cover_nonzero_per_HHS), 2)




round(mean(main_sensitivity_HHS) * 100, 2)
round(mean(interaction_sensitivity_HHS) * 100, 2)

round(mean(main_specificity_HHS) * 100, 2)
round(mean(interaction_specificity_HHS) * 100, 2)

round(mean(main_precision_HHS) * 100, 2)
round(mean(interaction_precision_HHS) * 100, 2)







###################################
# HRHS model

# load("results/theta_update_HRHS_list.rda")

# load("results/theta_update_HRHS_list_Sigma.rda")

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
  u <- sim_data_list[[sample]]$E
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






######################################
# NHHS model

# load("results/theta_update_NHHS_list.rda")

# load("results/theta_update_NHHS_list_Sigma.rda")

nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

alpha_cover_per_NHHS <- c()
beta_cover_per_NHHS <- c()
gamma_cover_per_NHHS <- c()
delta_cover_per_NHHS <- c()
delta_RMSE_NHHS <- c()
RMSE_all_main_effects_NHHS <- c()
beta_cover_nonzero_per_NHHS <- c()
gamma_cover_nonzero_per_NHHS <- c()
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

for(sample in 1:nsample)
{
  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$E
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


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

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_NHHS[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list[[sample]]$delta_true
  delta_true_list_NHHS[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHHS_list[[sample]][, jm, ], 2, mean)
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

  main_sensitivity_NHHS[sample] <- (main_TP / (main_TP + main_FN))
  main_specificity_NHHS[sample] <- (main_TN / (main_TN + main_FP))
  main_precision_NHHS[sample] <- (main_TP / (main_TP + main_FP))

  truth_delta <- as.vector(delta_true != 0)
  selected_delta <- as.vector((delta_lower > 0) | (delta_upper < 0))

  interaction_TP <- sum(truth_delta == 1 & selected_delta == 1)
  interaction_FN <- sum(truth_delta == 1 & selected_delta == 0)
  interaction_TN <- sum(truth_delta == 0 & selected_delta == 0)
  interaction_FP <- sum(truth_delta == 0 & selected_delta == 1)

  interaction_sensitivity_NHHS[sample] <- (interaction_TP / (interaction_TP + interaction_FN))
  interaction_specificity_NHHS[sample] <- (interaction_TN / (interaction_TN + interaction_FP))
  interaction_precision_NHHS[sample] <- (interaction_TP / (interaction_TP + interaction_FP))

}




round(mean(RMSE_all_main_effects_NHHS)*100, 2)
round(mean(delta_RMSE_NHHS)*100, 2)
round((mean(beta_cover_per_NHHS) + mean(gamma_cover_per_NHHS))/2, 2)
round((mean(beta_cover_nonzero_per_NHHS) + mean(gamma_cover_nonzero_per_NHHS))/2, 2)
round(mean(delta_cover_per_NHHS), 2)
round(mean(delta_cover_nonzero_per_NHHS), 2)




round(mean(main_sensitivity_NHHS) * 100, 2)
round(mean(interaction_sensitivity_NHHS) * 100, 2)

round(mean(main_specificity_NHHS) * 100, 2)
round(mean(interaction_specificity_NHHS) * 100, 2)

round(mean(main_precision_NHHS) * 100, 2)
round(mean(interaction_precision_NHHS) * 100, 2)







#############################
# NHRHS model

# load("results/theta_update_NHRHS_list.rda")

# load("results/theta_update_NHRHS_list_Sigma.rda")

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
  u <- sim_data_list[[sample]]$E
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

# load("results/theta_update_NP_list.rda")

# load("results/theta_update_NP_list_Sigma.rda")

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
  u <- sim_data_list[[sample]]$E
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








