
rm(list = ls())
gc()


########################
Sim_data_MVN <- function(K = 4, n = 1301, J = 44, M = 17, O = 6)
{
  # modifiers variables
  x <- matrix(rnorm(n * J), nrow = n)
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- rnorm(K)
  beta_true[2, ] <- rnorm(K)

  # environmental variables
  z <- matrix(rnorm(n * M), nrow = n)
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different beta for different k
  gamma_true[1, ] <- rnorm(K)
  gamma_true[2, ] <- rnorm(K)
  gamma_true[3, ] <- rnorm(K)


  # gene-environment interaction terms
  u <- matrix(0, n, J * M)
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      u[, col_idx] <- x[, j] * z[, m]
      col_idx <- col_idx + 1
    }
  }
  delta_true <-  matrix(0, nrow = J*M, ncol = K) # Different delta for different k
  active_beta  <- which(rowSums(abs(beta_true)) > 0)
  active_gamma <- which(rowSums(abs(gamma_true)) > 0)
  for (j in active_beta) {
    for (m in active_gamma) {
      idx <- (j - 1) * M + m
      delta_true[idx, ] <- rnorm(K)
    }
  }

  # covariates
  d <- matrix(rnorm(n * O), nrow = n)
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- rnorm(K)
  varphi_true[2, ] <- rnorm(K)

  alpha_true <- rnorm(K)

  theta_true <- rbind(matrix(alpha_true, nrow = 1),
                      beta_true,
                      gamma_true,
                      delta_true,
                      varphi_true)
  W <- cbind(1, x, z, u, d)
  W_theta_true <- W %*% theta_true # n*K dimension

  # Wishart-based PD matrix
  A <- matrix(rnorm(K*K), K, K)
  Sigma <- crossprod(A)   # AᵀA is automatically positive definite

  # Cholesky of Sigma
  R_Sigma <- chol(Sigma)

  # E ~ iid N(0,1)
  E_mat <- matrix(rnorm(n * K), nrow = n, ncol = K)

  # Y = M + E * chol(Sigma)
  Y <- W_theta_true + E_mat %*% R_Sigma

  return(list(
    x = x, z = z, u = u, d = d, W = W, Y = Y,
    alpha_true = alpha_true,
    beta_true = beta_true,
    gamma_true = gamma_true,
    delta_true = delta_true,
    varphi_true = varphi_true,
    W_theta_true = W_theta_true,
    Sigma = Sigma,
    theta_true = theta_true
  ))
}



#########################################################################

Sim_data_MVN_weak <- function(K = 4, n = 1301, J = 44, M = 17, O = 6)
{
  # modifiers variables
  x <- matrix(rnorm(n * J), nrow = n)
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- rnorm(K)
  beta_true[2, ] <- rnorm(K)

  # environmental variables
  z <- matrix(rnorm(n * M), nrow = n)
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different beta for different k
  gamma_true[1, ] <- rnorm(K)
  gamma_true[2, ] <- rnorm(K)
  gamma_true[3, ] <- rnorm(K)


  # gene-environment interaction terms
  u <- matrix(0, n, J * M)
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      u[, col_idx] <- x[, j] * z[, m]
      col_idx <- col_idx + 1
    }
  }

  delta_true <- matrix(0, nrow = J * M, ncol = K)
  active_beta  <- which(rowSums(abs(beta_true)) > 0)
  active_gamma <- which(rowSums(abs(gamma_true)) > 0)
  # Weak heredity: at least one corresponding main effect is active
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      if (j %in% active_beta || m %in% active_gamma)
      {
        idx <- (j - 1) * M + m
        delta_true[idx, ] <- rnorm(K)
      }
    }
  }


  # covariates
  d <- matrix(rnorm(n * O), nrow = n)
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- rnorm(K)
  varphi_true[2, ] <- rnorm(K)

  alpha_true <- rnorm(K)

  theta_true <- rbind(matrix(alpha_true, nrow = 1),
                      beta_true,
                      gamma_true,
                      delta_true,
                      varphi_true)
  W <- cbind(1, x, z, u, d)
  W_theta_true <- W %*% theta_true # n*K dimension

  # Wishart-based PD matrix
  A <- matrix(rnorm(K*K), K, K)
  Sigma <- crossprod(A)   # AᵀA is automatically positive definite

  # Cholesky of Sigma
  R_Sigma <- chol(Sigma)

  # E ~ iid N(0,1)
  E_mat <- matrix(rnorm(n * K), nrow = n, ncol = K)

  # Y = M + E * chol(Sigma)
  Y <- W_theta_true + E_mat %*% R_Sigma

  return(list(
    x = x, z = z, u = u, d = d, W = W, Y = Y,
    alpha_true = alpha_true,
    beta_true = beta_true,
    gamma_true = gamma_true,
    delta_true = delta_true,
    varphi_true = varphi_true,
    W_theta_true = W_theta_true,
    Sigma = Sigma,
    theta_true = theta_true
  ))
}



#############################################

# J = 10
sim_data_list <- list()
for(sample in 1:10)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_MVN(K = 4, n = 200, J = 10, M = 5, O = 3)
  sim_data_list[[sample]] <- sim_data
}
# Took 5 minutes to complete


# J = 20
sim_data_list <- list()
for(sample in 1:10)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_MVN(K = 4, n = 200, J = 20, M = 5, O = 3)
  sim_data_list[[sample]] <- sim_data
}




# J = 20
sim_data_list <- list()
for(sample in 1:10)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_MVN_weak(K = 4, n = 200, J = 20, M = 5, O = 3)
  sim_data_list[[sample]] <- sim_data
}




# J = 44
sim_data_list <- list()
for(sample in 1:10)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_MVN(K = 4, n = 200, J = 44, M = 5, O = 3)
  sim_data_list[[sample]] <- sim_data
}



# J = 44, O = 6
sim_data_list <- list()
for(sample in 1:10)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_MVN(K = 4, n = 200, J = 44, M = 5, O = 6)
  sim_data_list[[sample]] <- sim_data
}


# J = 44, M = 17, O = 6
sim_data_list <- list()
for(sample in 1:10)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_MVN(K = 4, n = 200, J = 44, M = 17, O = 6)
  sim_data_list[[sample]] <- sim_data
}

############################
# HHSSI model
############################

# Run all the functions from the package
devtools::load_all()

library(parallel)
library(foreach)
library(doParallel)
library(doRNG)

n_cores <- 6

cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export objects/functions needed by workers
clusterExport(cl,
  varlist = c(
    "fit_HHSSI_MVN_cov_modify3",
    "sim_data_list"
  ),
  envir = environment()
)

start_time_HHSSI <- Sys.time()

theta_update_HHSSI_list <- foreach(
  sample = 1:10,
  .packages = c("LaplacesDemon", "Matrix", "actuar"),
  .combine = "list",
  .multicombine = TRUE
) %dorng% {

  cat("Running sample", sample, "\n")

  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  d <- sim_data_list[[sample]]$d

  W <- cbind(1, x, z, u, d)
  WTW <- t(W) %*% W

  K <- ncol(Y)
  n <- nrow(Y)
  J <- ncol(x)
  M <- ncol(z)
  O <- ncol(d)
  n_all_par <- ncol(W)

  nu_0 <- K + 2
  Psi_0 <- diag(K)

  Sigma_init <- rinvwishart(nu_0, Psi_0)

  tryCatch({

    res <- fit_HHSSI_MVN_cov_modify3(
      niter = 2000,
      burn_in = 500,
      thin = 5,
      n = n,
      K = K,
      Y = Y,
      W = W,
      WTW = WTW,
      n_all_par = n_all_par,
      J = J,
      M = M,
      O = O,
      sigmasq_alpha = 100,
      theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
      lambdasq_beta_init = rep(0.5, J),
      tausq_beta_init = 1,
      lambdasq_gamma_init = rep(0.5, M),
      tausq_gamma_init = 1,
      lambdasq_delta_init = rep(0.5, J * M),
      tausq_delta_init = 1,
      psi_beta_init = rep(0.5, J),
      psi_gamma_init = rep(0.5, M),
      psi_delta_init = rep(0.5, J * M),
      xi_beta_init = 1,
      xi_gamma_init = 1,
      xi_delta_init = 1,
      Sigma_init = Sigma_init,
      nu_0 = nu_0,
      Psi_0 = Psi_0,
      sigmasq_varphi = 10
    )
    res$theta_update

  }, error = function(e) {

    list(error = TRUE, message = e$message, sample_id = sample)
  })
}

end_time_HHSSI <- Sys.time()
stopCluster(cl)

cat("Completed in", round(difftime(end_time_HHSSI,
    start_time_HHSSI, units = "mins"), 2), "minutes\n")


###############################################################


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

  delta_RMSE_indices_mat_HHSSI[sample, ] <- sqrt(rowMeans(delta_MSE))
  delta_cover_indices_mat_HHSSI[sample, ] <- rowMeans(delta_cover)*100
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




#######################################################################

# NHHSSI model

########################################################################

# Run all the functions from the package
devtools::load_all()

library(parallel)
library(foreach)
library(doParallel)
library(doRNG)

n_cores <- 6

cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export objects/functions needed by workers
clusterExport(cl,
              varlist = c(
                "fit_NHHSSI_MVN_cov_modify3",
                "sim_data_list"
              ),
              envir = environment()
)

start_time_NHHSSI <- Sys.time()

theta_update_NHHSSI_list <- foreach(
  sample = 1:10,
  .packages = c("LaplacesDemon", "Matrix", "actuar"),
  .combine = "list",
  .multicombine = TRUE
) %dorng% {

  cat("Running sample", sample, "\n")

  Y <- sim_data_list[[sample]]$Y
  x <- sim_data_list[[sample]]$x
  z <- sim_data_list[[sample]]$z
  u <- sim_data_list[[sample]]$u
  d <- sim_data_list[[sample]]$d

  W <- cbind(1, x, z, u, d)
  WTW <- t(W) %*% W

  K <- ncol(Y)
  n <- nrow(Y)
  J <- ncol(x)
  M <- ncol(z)
  O <- ncol(d)
  n_all_par <- ncol(W)

  nu_0 <- K + 2
  Psi_0 <- diag(K)

  Sigma_init <- rinvwishart(nu_0, Psi_0)

  tryCatch({

    res <- fit_NHHSSI_MVN_cov_modify3(
      niter = 2000, burn_in = 500, thin = 5,
      n=n, K=K, Y=Y, W=W,  WTW = WTW,
      n_all_par = n_all_par, J = J, M = M, O = O,
      sigmasq_alpha = 100,
      theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
      lambdasq_beta_init = rep(0.5, J),
      tausq_beta_init = 1,
      lambdasq_gamma_init = rep(0.5, M),
      tausq_gamma_init = 1,
      lambdasq_delta_init = rep(0.5, J*M),
      tausq_delta_init = 1,
      psi_beta_init = rep(0.5, J),
      psi_gamma_init = rep(0.5, M),
      psi_delta_init = rep(0.5, J*M),
      xi_beta_init = 1, xi_gamma_init = 1,
      xi_delta_init = 1,
      Sigma_init = Sigma_init,
      nu_0 = nu_0, Psi_0 = Psi_0,
      sigmasq_varphi = 10
    )
    res$theta_update

  }, error = function(e) {

    list(error = TRUE, message = e$message, sample_id = sample)
  })
}

end_time_NHHSSI <- Sys.time()
stopCluster(cl)

cat("Completed in", round(difftime(end_time_NHHSSI,
                                   start_time_NHHSSI, units = "mins"), 2), "minutes\n")


###############################################################


nsample <- length(sim_data_list)
K <- dim(sim_data_list[[1]]$Y)[2]

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

  delta_RMSE_indices_mat_NHHSSI[sample, ] <- sqrt(rowMeans(delta_MSE))
  delta_cover_indices_mat_NHHSSI[sample, ] <- rowMeans(delta_cover)*100
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





##########################################################
# separate data analysis to understand why does nhhssi perform better than hhssi
set.seed(23444)
sim_data_normal01_simpred1 <- Sim_data_MVN(K = 4, n = 200, J = 20, M = 5, O = 3)




############################
# HHSSI model
############################

# Run all the functions from the package
devtools::load_all()

  Y <- sim_data_normal01_simpred1$Y
  x <- sim_data_normal01_simpred1$x
  z <- sim_data_normal01_simpred1$z
  u <- sim_data_normal01_simpred1$u
  d <- sim_data_normal01_simpred1$d

  W <- cbind(1, x, z, u, d)
  WTW <- t(W) %*% W

  K <- ncol(Y)
  n <- nrow(Y)
  J <- ncol(x)
  M <- ncol(z)
  O <- ncol(d)
  n_all_par <- ncol(W)

  nu_0 <- K + 2
  Psi_0 <- diag(K)

  Sigma_init <- rinvwishart(nu_0, Psi_0)

  start_time <- Sys.time()
  set.seed(456557)
  res_normal01_simpred_hhssi1 <- fit_HHSSI_MVN_cov_modify3(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
    sigmasq_alpha = 100,
    theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
    lambdasq_beta_init = rep(0.5, J),
    tausq_beta_init = 1,
    lambdasq_gamma_init = rep(0.5, M),
    tausq_gamma_init = 1,
    lambdasq_delta_init = rep(0.5, J*M),
    tausq_delta_init = 1,
    psi_beta_init = rep(0.5, J),
    psi_gamma_init = rep(0.5, M),
    psi_delta_init = rep(0.5, J*M),
    xi_beta_init = 1, xi_gamma_init = 1,
    xi_delta_init = 1,
    Sigma_init = Sigma_init,
    nu_0 = nu_0, Psi_0 = Psi_0,
    sigmasq_varphi = 10
  )
  end_time <- Sys.time()
  end_time - start_time




  start_time <- Sys.time()
  set.seed(456557)
  res_normal01_simpred_nhhssi1 <- fit_NHHSSI_MVN_cov_modify3(
    niter = 6000, burn_in = 1000, thin = 5,
    n=n, K=K, Y=Y, W=W,  WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
    sigmasq_alpha = 100,
    theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
    lambdasq_beta_init = rep(0.5, J),
    tausq_beta_init = 1,
    lambdasq_gamma_init = rep(0.5, M),
    tausq_gamma_init = 1,
    lambdasq_delta_init = rep(0.5, J*M),
    tausq_delta_init = 1,
    psi_beta_init = rep(0.5, J),
    psi_gamma_init = rep(0.5, M),
    psi_delta_init = rep(0.5, J*M),
    xi_beta_init = 1, xi_gamma_init = 1,
    xi_delta_init = 1,
    Sigma_init = Sigma_init,
    nu_0 = nu_0, Psi_0 = Psi_0,
    sigmasq_varphi = 10
  )
  end_time <- Sys.time()
  end_time - start_time




###############################################################################################
# graphs
###############################################################################################

  #fs sim
  beta_lam_simpred_hhssi1  <- res_normal01_simpred_hhssi1$lambdasq_beta_update
  gamma_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_gamma_update
  delta_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_delta_update
  tausq_delta_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update

  delta_lam_simpred_nhhssi1  <- res_normal01_simpred_nhhssi1$lambdasq_delta_update
  tausq_delta_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update




  n_iter <- dim(beta_lam_simpred_hhssi1)[1]



####################################################################
# for normal01_simpred_hhssi1, normal01_simpred_nhhssi1

  # Construct beta x gamma component
  beta_gamma_lam_simpred_hhssi1 <- matrix(NA, nrow = n_iter, ncol = J * M)
  for (t in 1:n_iter)
  {
    col_idx <- 1
    for (j in 1:J) {
      for (m in 1:M) {
        beta_gamma_lam_simpred_hhssi1[t, col_idx] <- (beta_lam_simpred_hhssi1[t, j] *
                                                        gamma_lam_simpred_hhssi1[t, m])
        col_idx <- col_idx + 1
      }
    }
  }
  beta_gamma_lam_mean_simpred_hhssi1 <- apply(beta_gamma_lam_simpred_hhssi1, 2, mean)

  # prior variance from hhssi
  delta_var_simpred_hhssi1 <- (delta_lam_simpred_hhssi1 *
                                 beta_gamma_lam_simpred_hhssi1 *
                                 tausq_delta_simpred_hhssi1)
  delta_var_mean_simpred_hhssi1 <- apply(delta_var_simpred_hhssi1, 2, mean)
  dim(delta_var_simpred_hhssi1)
  # prior variance from nhhssi
  delta_var_simpred_nhhssi1 <- (delta_lam_simpred_nhhssi1 *
                                  tausq_delta_simpred_nhhssi1)
  delta_var_mean_simpred_nhhssi1 <- apply(delta_var_simpred_nhhssi1, 2, mean)


  # prior var estimates from hhssi vs nhhssi
  plot(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1),
       pch = 19,
       xlab = "log HHS-SI delta var mean",
       ylab = "log NHHS-SI delta var mean",
       main ="log transformed",
       xlim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))),
       ylim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)),
       ylim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
  plot(delta_var_mean_simpred_hhssi1[idx],
       delta_var_mean_simpred_nhhssi1[idx],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(delta_var_mean_simpred_hhssi1[idx],
                      delta_var_mean_simpred_nhhssi1[idx])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[idx],
                      delta_var_mean_simpred_nhhssi1[idx])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")

  plot(delta_var_mean_simpred_hhssi1[!idx],
       delta_var_mean_simpred_nhhssi1[!idx],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                      delta_var_mean_simpred_nhhssi1[!idx])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                      delta_var_mean_simpred_nhhssi1[!idx])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  # for active delta
  plot(delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       main="for corresponding active deltas",
       xlim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(0, 1, lty = 2, lwd = 2, col = "red")



  # for inactive delta
  plot(delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       main="for corresponding inactive deltas",
       xlim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)])))
  abline(0, 1, lty = 2, lwd = 2, col = "red")



  inactive_simpred_hhssi1 <- delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)]
  inactive_simpred_nhhssi1 <- delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)]
  inxinc <- inactive_simpred_nhhssi1 > inactive_simpred_hhssi1
  plot(inactive_simpred_hhssi1[inxinc],
       inactive_simpred_nhhssi1[inxinc],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(inactive_simpred_hhssi1[inxinc],
                      inactive_simpred_hhssi1[inxinc])),
       ylim = range(c(inactive_simpred_hhssi1[inxinc],
                      inactive_simpred_hhssi1[inxinc])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")

  plot(inactive_simpred_hhssi1[!inxinc],
       inactive_simpred_nhhssi1[!inxinc],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(inactive_simpred_hhssi1[!inxinc],
                      inactive_simpred_hhssi1[!inxinc])),
       ylim = range(c(inactive_simpred_hhssi1[!inxinc],
                      inactive_simpred_hhssi1[!inxinc])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



  # so the prior variance of lambdasq which are bigger from nhhssi model than hhssi model,
  # comes from the corresponding inactive delta values


  # other variance terms
  # lambdasq_beta
  lambdsq_beta_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update, 2, mean)
  lambdsq_beta_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update, 2, mean)
  plot(lambdsq_beta_mean_simpred_hhssi1,
       lambdsq_beta_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI lambdasq beta mean",
       ylab = "NHHS-SI lambdasq beta mean",
       main="overall",
       xlim = range(c(lambdsq_beta_mean_simpred_hhssi1,
                      lambdsq_beta_mean_simpred_nhhssi1)),
       ylim = range(c(lambdsq_beta_mean_simpred_hhssi1,
                      lambdsq_beta_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_beta_mean_simpred_hhssi1[1:2],
       lambdsq_beta_mean_simpred_nhhssi1[1:2],
       pch = 19,
       xlab = "HHS-SI lambdasq beta mean",
       ylab = "NHHS-SI lambdasq beta mean",
       main="for corresponding active parameters",
       xlim = range(c(lambdsq_beta_mean_simpred_hhssi1[1:2],
                      lambdsq_beta_mean_simpred_nhhssi1[1:2])),
       ylim = range(c(lambdsq_beta_mean_simpred_hhssi1[1:2],
                      lambdsq_beta_mean_simpred_nhhssi1[1:2])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_beta_mean_simpred_hhssi1[-(1:2)],
       lambdsq_beta_mean_simpred_nhhssi1[-(1:2)],
       pch = 19,
       xlab = "HHS-SI lambdasq beta mean",
       ylab = "NHHS-SI lambdasq beta mean",
       main="for corresponding inactive parameters",
       xlim = range(c(lambdsq_beta_mean_simpred_hhssi1[-(1:2)],
                      lambdsq_beta_mean_simpred_nhhssi1[-(1:2)])),
       ylim = range(c(lambdsq_beta_mean_simpred_hhssi1[-(1:2)],
                      lambdsq_beta_mean_simpred_nhhssi1[-(1:2)])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  # lambdasq_gamma
  lambdsq_gamma_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
  lambdsq_gamma_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update, 2, mean)
  plot(lambdsq_gamma_mean_simpred_hhssi1,
       lambdsq_gamma_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI lambdasq gamma mean",
       ylab = "NHHS-SI lambdasq gamma mean",
       main="overall",
       xlim = range(c(lambdsq_gamma_mean_simpred_hhssi1,
                      lambdsq_gamma_mean_simpred_nhhssi1)),
       ylim = range(c(lambdsq_gamma_mean_simpred_hhssi1,
                      lambdsq_gamma_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_gamma_mean_simpred_hhssi1[1:3],
       lambdsq_gamma_mean_simpred_nhhssi1[1:3],
       pch = 19,
       xlab = "HHS-SI lambdasq gamma mean",
       ylab = "NHHS-SI lambdasq gamma mean",
       main="for corresponding active parameters",
       xlim = range(c(lambdsq_gamma_mean_simpred_hhssi1[1:3],
                      lambdsq_gamma_mean_simpred_nhhssi1[1:3])),
       ylim = range(c(lambdsq_gamma_mean_simpred_hhssi1[1:3],
                      lambdsq_gamma_mean_simpred_nhhssi1[1:3])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_gamma_mean_simpred_hhssi1[-(1:3)],
       lambdsq_gamma_mean_simpred_nhhssi1[-(1:3)],
       pch = 19,
       xlab = "HHS-SI lambdasq gamma mean",
       ylab = "NHHS-SI lambdasq gamma mean",
       main="for corresponding inactive parameters",
       xlim = range(c(lambdsq_gamma_mean_simpred_hhssi1[-(1:3)],
                      lambdsq_gamma_mean_simpred_nhhssi1[-(1:3)])),
       ylim = range(c(lambdsq_gamma_mean_simpred_hhssi1[-(1:3)],
                      lambdsq_gamma_mean_simpred_nhhssi1[-(1:3)])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")





  # lambdasq_delta
  lambdsq_delta_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_delta_update, 2, mean)
  lambdsq_delta_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_delta_update, 2, mean)
  plot(lambdsq_delta_mean_simpred_hhssi1,
       lambdsq_delta_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="overall",
       xlim = range(c(lambdsq_delta_mean_simpred_hhssi1,
                      lambdsq_delta_mean_simpred_nhhssi1)),
       ylim = range(c(lambdsq_delta_mean_simpred_hhssi1,
                      lambdsq_delta_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
       lambdsq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)],
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="for corresponding active parameters",
       xlim = range(c(lambdsq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
                      lambdsq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(lambdsq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
                      lambdsq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
       lambdsq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))],
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="for corresponding inactive parameters",
       xlim = range(c(lambdsq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
                      lambdsq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))])),
       ylim = range(c(lambdsq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
                      lambdsq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  # tausq_beta
  tausq_beta_mean_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_beta_update
  tausq_beta_mean_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_beta_update
  plot(tausq_beta_mean_simpred_hhssi1,
       tausq_beta_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI tausq beta",
       ylab = "NHHS-SI tausq beta",
       main="overall",
       xlim = range(c(tausq_beta_mean_simpred_hhssi1,
                      tausq_beta_mean_simpred_nhhssi1)),
       ylim = range(c(tausq_beta_mean_simpred_hhssi1,
                      tausq_beta_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(tausq_beta_mean_simpred_hhssi1[1:2],
       tausq_beta_mean_simpred_nhhssi1[1:2],
       pch = 19,
       xlab = "HHS-SI tausq beta",
       ylab = "NHHS-SI tausq beta",
       main="for corresponding active parameters",
       xlim = range(c(tausq_beta_mean_simpred_hhssi1[1:2],
                      tausq_beta_mean_simpred_nhhssi1[1:2])),
       ylim = range(c(tausq_beta_mean_simpred_hhssi1[1:2],
                      tausq_beta_mean_simpred_nhhssi1[1:2])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(tausq_beta_mean_simpred_hhssi1[-(1:2)],
       tausq_beta_mean_simpred_nhhssi1[-(1:2)],
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="for corresponding inactive parameters",
       xlim = range(c(tausq_beta_mean_simpred_hhssi1[-(1:2)],
                      tausq_beta_mean_simpred_nhhssi1[-(1:2)])),
       ylim = range(c(tausq_beta_mean_simpred_hhssi1[-(1:2)],
                      tausq_beta_mean_simpred_nhhssi1[-(1:2)])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



  # tausq_gamma
  tausq_gamma_mean_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_gamma_update
  tausq_gamma_mean_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_gamma_update
  plot(tausq_gamma_mean_simpred_hhssi1,
       tausq_gamma_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI tausq gamma",
       ylab = "NHHS-SI tausq gamma",
       main="overall",
       xlim = range(c(tausq_gamma_mean_simpred_hhssi1,
                      tausq_gamma_mean_simpred_nhhssi1)),
       ylim = range(c(tausq_gamma_mean_simpred_hhssi1,
                      tausq_gamma_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(tausq_gamma_mean_simpred_hhssi1[1:3],
       tausq_gamma_mean_simpred_nhhssi1[1:3],
       pch = 19,
       xlab = "HHS-SI tausq gamma",
       ylab = "NHHS-SI tausq gamma",
       main="for corresponding active parameters",
       xlim = range(c(tausq_gamma_mean_simpred_hhssi1[1:3],
                      tausq_gamma_mean_simpred_nhhssi1[1:3])),
       ylim = range(c(tausq_gamma_mean_simpred_hhssi1[1:3],
                      tausq_gamma_mean_simpred_nhhssi1[1:3])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(tausq_gamma_mean_simpred_hhssi1[-(1:3)],
       tausq_gamma_mean_simpred_nhhssi1[-(1:3)],
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="for corresponding inactive parameters",
       xlim = range(c(tausq_gamma_mean_simpred_hhssi1[-(1:3)],
                      tausq_gamma_mean_simpred_nhhssi1[-(1:3)])),
       ylim = range(c(tausq_gamma_mean_simpred_hhssi1[-(1:3)],
                      tausq_gamma_mean_simpred_nhhssi1[-(1:3)])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")







  # tausq_delta
  tausq_delta_mean_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update
  tausq_delta_mean_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update
  plot(tausq_delta_mean_simpred_hhssi1,
       tausq_delta_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI tausq delta",
       ylab = "NHHS-SI tausq delta",
       main="overall",
       xlim = range(c(tausq_delta_mean_simpred_hhssi1,
                      tausq_delta_mean_simpred_nhhssi1)),
       ylim = range(c(tausq_delta_mean_simpred_hhssi1,
                      tausq_delta_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(tausq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
       tausq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)],
       pch = 19,
       xlab = "HHS-SI tausq delta",
       ylab = "NHHS-SI tausq delta",
       main="for corresponding active parameters",
       xlim = range(c(tausq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
                      tausq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(tausq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
                      tausq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(tausq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
       tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))],
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="for corresponding inactive parameters",
       xlim = range(c(tausq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
                      tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))])),
       ylim = range(c(tausq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
                      tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")






  # lambdasq_delta*tausq_delta
  lambdsq_tausq_delta_mean_simpred_hhssi1 <- apply((res_normal01_simpred_hhssi1$lambdasq_delta_update*
                                                      res_normal01_simpred_hhssi1$tausq_delta_update),
                                                   2, mean)
  lambdsq_tausq_delta_mean_simpred_nhhssi1 <- apply((res_normal01_simpred_nhhssi1$lambdasq_delta_update*
                                                       res_normal01_simpred_nhhssi1$tausq_delta_update),
                                                    2, mean)
  plot(lambdsq_tausq_delta_mean_simpred_hhssi1,
       lambdsq_tausq_delta_mean_simpred_nhhssi1,
       pch = 19,
       xlab = "HHS-SI lambdasq*tausq delta mean",
       ylab = "NHHS-SI lambdasq*tausq delta mean",
       main="overall",
       xlim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1,
                      lambdsq_tausq_delta_mean_simpred_nhhssi1)),
       ylim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1,
                      lambdsq_tausq_delta_mean_simpred_nhhssi1)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_tausq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
       lambdsq_tausq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)],
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="for corresponding active parameters",
       xlim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
                      lambdsq_tausq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[c(1:3, 6:8)],
                      lambdsq_tausq_delta_mean_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


  plot(lambdsq_tausq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
       lambdsq_tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))],
       pch = 19,
       xlab = "HHS-SI lambdasq delta mean",
       ylab = "NHHS-SI lambdasq delta mean",
       main="for corresponding inactive parameters",
       xlim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
                      lambdsq_tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))])),
       ylim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[-(c(1:3, 6:8))],
                      lambdsq_tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 6:8))])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


####################################################################################
  # delta prior variance estimates
####################################################################
  # for normal01_simpred_hhssi1, normal01_simpred_nhhssi1

  #fs sim
  beta_lam_simpred_hhssi1  <- res_normal01_simpred_hhssi1$lambdasq_beta_update
  gamma_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_gamma_update
  delta_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_delta_update
  tausq_delta_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update

  delta_lam_simpred_nhhssi1  <- res_normal01_simpred_nhhssi1$lambdasq_delta_update
  tausq_delta_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update


  # Construct beta x gamma component
  beta_gamma_lam_simpred_hhssi1 <- matrix(NA, nrow = n_iter, ncol = J * M)
  for (t in 1:n_iter)
  {
    col_idx <- 1
    for (j in 1:J) {
      for (m in 1:M) {
        beta_gamma_lam_simpred_hhssi1[t, col_idx] <- (beta_lam_simpred_hhssi1[t, j] *
                                                        gamma_lam_simpred_hhssi1[t, m])
        col_idx <- col_idx + 1
      }
    }
  }
  beta_gamma_lam_mean_simpred_hhssi1 <- apply(beta_gamma_lam_simpred_hhssi1, 2, mean)

  # prior variance from hhssi
  delta_var_simpred_hhssi1 <- (delta_lam_simpred_hhssi1 *
                                 beta_gamma_lam_simpred_hhssi1 *
                                 tausq_delta_simpred_hhssi1)
  delta_var_mean_simpred_hhssi1 <- apply(delta_var_simpred_hhssi1, 2, mean)
  dim(delta_var_simpred_hhssi1)
  # prior variance from nhhssi
  delta_var_simpred_nhhssi1 <- (delta_lam_simpred_nhhssi1 *
                                  tausq_delta_simpred_nhhssi1)
  delta_var_mean_simpred_nhhssi1 <- apply(delta_var_simpred_nhhssi1, 2, mean)


  # prior var estimates from hhssi vs nhhssi
  idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
  point_col <- ifelse(idx, "blue", "red")

  par(mfrow=c(1,1))
  plot(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1),
       pch = 19,
       col = point_col,
       xlab = "log HHS-SI delta var mean",
       ylab = "log NHHS-SI delta var mean",
       main ="log transformed",
       xlim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))),
       ylim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))))
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  plot(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1,
       pch = 19,
       col = point_col,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)),
       ylim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)))
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
  plot(delta_var_mean_simpred_hhssi1[idx],
       delta_var_mean_simpred_nhhssi1[idx],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(delta_var_mean_simpred_hhssi1[idx],
                      delta_var_mean_simpred_nhhssi1[idx])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[idx],
                      delta_var_mean_simpred_nhhssi1[idx])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



  plot(delta_var_mean_simpred_hhssi1[!idx],
       delta_var_mean_simpred_nhhssi1[!idx],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       xlim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                      delta_var_mean_simpred_nhhssi1[!idx])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                      delta_var_mean_simpred_nhhssi1[!idx])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  # for active delta
  plot(delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)],
       pch = 19,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       main="for corresponding active deltas",
       xlim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(0, 1, lty = 2, lwd = 2, col = "red")



  # for inactive delta
  idx_inct <- delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)] > delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)]
  point_col_inct <- ifelse(idx_inct, "blue", "red")
  plot(delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)],
       pch = 19,
       col= point_col_inct,
       xlab = "HHS-SI delta var mean",
       ylab = "NHHS-SI delta var mean",
       main="for corresponding inactive deltas",
       xlim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)])),
       ylim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)])))
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0, 1, lty = 2, lwd = 2, col = "red")


  #######################################################################################
  # sd of prior variance estimate
  sqrt_delta_var_mean_simpred_hhssi1 <- sqrt(delta_var_mean_simpred_hhssi1)
  sqrt_delta_var_mean_simpred_nhhssi1 <- sqrt(delta_var_mean_simpred_nhhssi1)

  plot(sqrt_delta_var_mean_simpred_hhssi1, sqrt_delta_var_mean_simpred_nhhssi1,
       pch = 19,
       col = point_col,
       xlab = "HHS-SI sqrt of delta var",
       ylab = "NHHS-SI sqrt of delta var",
       xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1, sqrt_delta_var_mean_simpred_nhhssi1)),
       ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1, sqrt_delta_var_mean_simpred_nhhssi1)))
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  idx <- sqrt_delta_var_mean_simpred_nhhssi1 > sqrt_delta_var_mean_simpred_hhssi1
  plot(sqrt_delta_var_mean_simpred_hhssi1[idx],
       sqrt_delta_var_mean_simpred_nhhssi1[idx],
       pch = 19,
       xlab = "HHS-SI sqrt of delta var",
       ylab = "NHHS-SI sqrt of delta var",
       xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[idx],
                      sqrt_delta_var_mean_simpred_nhhssi1[idx])),
       ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[idx],
                      sqrt_delta_var_mean_simpred_nhhssi1[idx])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



  plot(sqrt_delta_var_mean_simpred_hhssi1[!idx],
       sqrt_delta_var_mean_simpred_nhhssi1[!idx],
       pch = 19,
       xlab = "HHS-SI sqrt of delta var",
       ylab = "NHHS-SI sqrt of delta var",
       xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[!idx],
                      sqrt_delta_var_mean_simpred_nhhssi1[!idx])),
       ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[!idx],
                      sqrt_delta_var_mean_simpred_nhhssi1[!idx])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




  # for active delta
  plot(sqrt_delta_var_mean_simpred_hhssi1[c(1:3, 6:8)],
       sqrt_delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)],
       pch = 19,
       xlab = "HHS-SI sqrt of delta var",
       ylab = "NHHS-SI sqrt of delta var",
       main="for corresponding active deltas",
       xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], sqrt_delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[c(1:3, 6:8)], sqrt_delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(0, 1, lty = 2, lwd = 2, col = "red")



  # for inactive delta
  plot(sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)],
       sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)],
       pch = 19,
       col = point_col_inct,
       xlab = "HHS-SI sqrt of delta var",
       ylab = "NHHS-SI sqrt of delta var",
       main="for corresponding inactive deltas",
       xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)])),
       ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)], sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)])))
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0, 1, lty = 2, lwd = 2, col = "red")



  sqrt_inactive_simpred_hhssi1 <- sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)]
  sqrt_inactive_simpred_nhhssi1 <- sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)]
  inxinc <- sqrt_inactive_simpred_nhhssi1 > sqrt_inactive_simpred_hhssi1
  plot(sqrt_inactive_simpred_hhssi1[inxinc],
       sqrt_inactive_simpred_nhhssi1[inxinc],
       pch = 19,
       xlab = "HHS-SI sqrt of delta var",
       ylab = "NHHS-SI sqrt of delta var",
       xlim = range(c(sqrt_inactive_simpred_hhssi1[inxinc],
                      sqrt_inactive_simpred_hhssi1[inxinc])),
       ylim = range(c(sqrt_inactive_simpred_hhssi1[inxinc],
                      sqrt_inactive_simpred_hhssi1[inxinc])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")

  plot(sqrt_inactive_simpred_hhssi1[!inxinc],
       sqrt_inactive_simpred_nhhssi1[!inxinc],
       pch = 19,
       xlab = "HHS-SI sqrt of delta var",
       ylab = "NHHS-SI sqrt of delta var",
       xlim = range(c(sqrt_inactive_simpred_hhssi1[!inxinc],
                      sqrt_inactive_simpred_hhssi1[!inxinc])),
       ylim = range(c(sqrt_inactive_simpred_hhssi1[!inxinc],
                      sqrt_inactive_simpred_hhssi1[!inxinc])))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")





  #############################################################################################
  # Posterior sd estimate of delta
  delta_sd_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 27:126, ], c(2,3), sd)
  delta_sd_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 27:126, ], c(2,3), sd)

  delta_est_var_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 27:126, ], c(2,3), var)
  delta_est_var_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 27:126, ], c(2,3), var)



  # Posterior mean estimate of delta
  delta_est_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 27:126, ], c(2,3), mean)
  delta_est_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 27:126, ], c(2,3), mean)

  # Bias for each outcome
  delta_bias_normal01_simpred_hhssi1 <- (delta_est_normal01_simpred_hhssi1 - sim_data_normal01_simpred1$delta_true)
  delta_bias_normal01_simpred_nhhssi1 <- (delta_est_normal01_simpred_nhhssi1 - sim_data_normal01_simpred1$delta_true)



  mean(delta_est_var_normal01_simpred_hhssi1)
  mean(delta_est_var_normal01_simpred_nhhssi1)

  mean(delta_bias_normal01_simpred_hhssi1)
  mean(delta_bias_normal01_simpred_nhhssi1)

  ##########################################################################################


  # sd of Delta Estimates
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_sd_normal01_simpred_hhssi1[, k],
      delta_sd_normal01_simpred_nhhssi1[, k],
      col =point_col,
      xlab = "HHSSI delta sd",
      ylab = "NHHSSI delta sd",
      main = paste("Overall sd: Outcome", k),
      xlim = range(c(delta_sd_normal01_simpred_hhssi1[,k],
                     delta_sd_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(delta_sd_normal01_simpred_hhssi1[,k],
                     delta_sd_normal01_simpred_nhhssi1[,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")



  # active delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_sd_normal01_simpred_hhssi1[c(1:3, 6:8), k],
      delta_sd_normal01_simpred_nhhssi1[c(1:3, 6:8), k],
      xlab = "HHSSI delta sd",
      ylab = "NHHSSI delta sd",
      main = paste("active: Outcome", k),
      xlim = range(c(delta_sd_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_sd_normal01_simpred_nhhssi1[c(1:3, 6:8),k])),
      ylim = range(c(delta_sd_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_sd_normal01_simpred_nhhssi1[c(1:3, 6:8),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8), k],
      delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8), k],
      xlab = "HHSSI delta sd",
      ylab = "NHHSSI delta sd",
      col = point_col_inct,
      main = paste("inactive: Outcome", k),
      xlim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])),
      ylim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
    legend("topleft",
           legend = c(
             "NHHSSI variance > HHSSI variance",
             "NHHSSI variance <= HHSSI variance"),
           col = c("blue", "red"),
           pch = 19,
           bty = "n")
  }






  idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
  # idx
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_sd_normal01_simpred_hhssi1[idx, k],
      delta_sd_normal01_simpred_nhhssi1[idx, k],
      xlab = "HHSSI delta sd",
      ylab = "NHHSSI delta sd",
      main = paste("Overall: Outcome", k),
      xlim = range(c(delta_sd_normal01_simpred_hhssi1[idx,k],
                     delta_sd_normal01_simpred_nhhssi1[idx,k])),
      ylim = range(c(delta_sd_normal01_simpred_hhssi1[idx,k],
                     delta_sd_normal01_simpred_nhhssi1[idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }


  # !idx
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_sd_normal01_simpred_hhssi1[!idx, k],
      delta_sd_normal01_simpred_nhhssi1[!idx, k],
      xlab = "HHSSI delta sd",
      ylab = "NHHSSI delta sd",
      main = paste("Overall: Outcome", k),
      xlim = range(c(delta_sd_normal01_simpred_hhssi1[!idx,k],
                     delta_sd_normal01_simpred_nhhssi1[!idx,k])),
      ylim = range(c(delta_sd_normal01_simpred_hhssi1[!idx,k],
                     delta_sd_normal01_simpred_nhhssi1[!idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }





  idx_act <- delta_var_mean_simpred_hhssi1[c(1:3, 6:8)] < delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)]
  # no active delta with nhhssi var bigger than hhssi

  # inactive delta
  idx_inact <- delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)]
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8), ][idx_inact, k],
      delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8), ][idx_inact, k],
      xlab = "HHSSI delta sd",
      ylab = "NHHSSI delta sd",
      main = paste("Outcome", k),
      xlim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8),][idx_inact, k],
                     delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8),][idx_inact, k])),
      ylim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8),][idx_inact, k],
                     delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8),][idx_inact, k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }


  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8), ][!idx_inact, k],
      delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8), ][!idx_inact, k],
      xlab = "HHSSI delta sd",
      ylab = "NHHSSI delta sd",
      main = paste("Outcome", k),
      xlim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8),][!idx_inact, k],
                     delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8),][!idx_inact, k])),
      ylim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 6:8),][!idx_inact, k],
                     delta_sd_normal01_simpred_nhhssi1[-c(1:3, 6:8),][!idx_inact, k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }





  ##########################################################################################

  # variance of Delta Estimates
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_var_normal01_simpred_hhssi1[, k],
      delta_est_var_normal01_simpred_nhhssi1[, k],
      col =point_col,
      xlab = "HHSSI delta var",
      ylab = "NHHSSI delta var",
      main = paste("Overall sd: Outcome", k),
      xlim = range(c(delta_est_var_normal01_simpred_hhssi1[,k],
                     delta_est_var_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(delta_est_var_normal01_simpred_hhssi1[,k],
                     delta_est_var_normal01_simpred_nhhssi1[,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")



  # active delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_var_normal01_simpred_hhssi1[c(1:3, 6:8), k],
      delta_est_var_normal01_simpred_nhhssi1[c(1:3, 6:8), k],
      xlab = "HHSSI delta var",
      ylab = "NHHSSI delta var",
      main = paste("active: Outcome", k),
      xlim = range(c(delta_est_var_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_est_var_normal01_simpred_nhhssi1[c(1:3, 6:8),k])),
      ylim = range(c(delta_est_var_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_est_var_normal01_simpred_nhhssi1[c(1:3, 6:8),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_var_normal01_simpred_hhssi1[-c(1:3, 6:8), k],
      delta_est_var_normal01_simpred_nhhssi1[-c(1:3, 6:8), k],
      xlab = "HHSSI delta var",
      ylab = "NHHSSI delta var",
      col = point_col_inct,
      main = paste("inactive: Outcome", k),
      xlim = range(c(delta_est_var_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_est_var_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])),
      ylim = range(c(delta_est_var_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_est_var_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
    legend("topleft",
           legend = c(
             "NHHSSI variance > HHSSI variance",
             "NHHSSI variance <= HHSSI variance"),
           col = c("blue", "red"),
           pch = 19,
           bty = "n")
  }





  ##########################################################################################


  # Delta Estimates
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_normal01_simpred_hhssi1[, k],
      delta_est_normal01_simpred_nhhssi1[, k],
      xlab = "HHSSI delta estimates",
      ylab = "NHHSSI delta estimates",
      col = point_col,
      main = paste("Overall estimates: Outcome", k),
      xlim = range(c(delta_est_normal01_simpred_hhssi1[,k],
                     delta_est_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(delta_est_normal01_simpred_hhssi1[,k],
                     delta_est_normal01_simpred_nhhssi1[,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
    legend("topleft",
           legend = c(
             "NHHSSI variance > HHSSI variance",
             "NHHSSI variance <= HHSSI variance"),
           col = c("blue", "red"),
           pch = 19,
           bty = "n")
  }

  # active delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_normal01_simpred_hhssi1[c(1:3, 6:8), k],
      delta_est_normal01_simpred_nhhssi1[c(1:3, 6:8), k],
      xlab = "HHSSI delta estimates",
      ylab = "NHHSSI delta estimates",
      main = paste("active: Outcome", k),
      xlim = range(c(delta_est_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_est_normal01_simpred_nhhssi1[c(1:3, 6:8),k])),
      ylim = range(c(delta_est_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_est_normal01_simpred_nhhssi1[c(1:3, 6:8),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8), k],
      delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8), k],
      xlab = "HHSSI delta estimates",
      ylab = "NHHSSI delta estimates",
      col = point_col_inct,
      main = paste("inactive: Outcome", k),
      xlim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])),
      ylim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])))
    legend("topleft",
           legend = c(
             "NHHSSI variance > HHSSI variance",
             "NHHSSI variance <= HHSSI variance"),
           col = c("blue", "red"),
           pch = 19,
           bty = "n")
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }






  idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_normal01_simpred_hhssi1[idx, k],
      delta_est_normal01_simpred_nhhssi1[idx, k],
      xlab = "HHSSI delta estimates",
      ylab = "NHHSSI delta estimates",
      main = paste("Overall: Outcome", k),
      xlim = range(c(delta_est_normal01_simpred_hhssi1[idx,k],
                     delta_est_normal01_simpred_nhhssi1[idx,k])),
      ylim = range(c(delta_est_normal01_simpred_hhssi1[idx,k],
                     delta_est_normal01_simpred_nhhssi1[idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }


  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_normal01_simpred_hhssi1[!idx, k],
      delta_est_normal01_simpred_nhhssi1[!idx, k],
      xlab = "HHSSI delta estimates",
      ylab = "NHHSSI delta estimates",
      main = paste("Overall: Outcome", k),
      xlim = range(c(delta_est_normal01_simpred_hhssi1[!idx,k],
                     delta_est_normal01_simpred_nhhssi1[!idx,k])),
      ylim = range(c(delta_est_normal01_simpred_hhssi1[!idx,k],
                     delta_est_normal01_simpred_nhhssi1[!idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }





  idx <- delta_var_mean_simpred_hhssi1[c(1:3, 6:8)] < delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)]
  # no active delta with nhhssi var bigger than hhssi

  # inactive delta
  idx <- delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)]
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8), ][idx, k],
      delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8), ][idx, k],
      xlab = "HHSSI delta estimates",
      ylab = "NHHSSI delta estimates",
      main = paste("inactive: Outcome", k),
      xlim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8),][idx, k],
                     delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8),][idx, k])),
      ylim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8),][idx, k],
                     delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8),][idx, k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }


  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8), ][!idx, k],
      delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8), ][!idx, k],
      xlab = "HHSSI delta estimates",
      ylab = "NHHSSI delta estimates",
      main = paste("inactive: Outcome", k),
      xlim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8),][!idx, k],
                     delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8),][!idx, k])),
      ylim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 6:8),][!idx, k],
                     delta_est_normal01_simpred_nhhssi1[-c(1:3, 6:8),][!idx, k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }





  ##########################################################################################
  # Delta Bias
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_bias_normal01_simpred_hhssi1[, k],
      delta_bias_normal01_simpred_nhhssi1[, k],
      col = point_col,
      xlab = "HHSSI delta bias",
      ylab = "NHHSSI delta bias",
      main = paste("Overall bias: Outcome", k),
      xlim = range(c(delta_bias_normal01_simpred_hhssi1[,k],
                     delta_bias_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(delta_bias_normal01_simpred_hhssi1[,k],
                     delta_bias_normal01_simpred_nhhssi1[,k])))
    legend("topleft",
           legend = c(
             "NHHSSI variance > HHSSI variance",
             "NHHSSI variance <= HHSSI variance"),
           col = c("blue", "red"),
           pch = 19,
           bty = "n")
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }




  # active delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_bias_normal01_simpred_hhssi1[c(1:3, 6:8), k],
      delta_bias_normal01_simpred_nhhssi1[c(1:3, 6:8), k],
      xlab = "HHSSI delta bias",
      ylab = "NHHSSI delta bias",
      main = paste("active: Outcome", k),
      xlim = range(c(delta_bias_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_bias_normal01_simpred_nhhssi1[c(1:3, 6:8),k])),
      ylim = range(c(delta_bias_normal01_simpred_hhssi1[c(1:3, 6:8),k],
                     delta_bias_normal01_simpred_nhhssi1[c(1:3, 6:8),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive delta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8), k],
      delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8), k],
      xlab = "HHSSI delta bias",
      ylab = "NHHSSI delta bias",
      col = point_col_inct,
      main = paste("inactive: Outcome", k),
      xlim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])),
      ylim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8),k],
                     delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
    legend("topleft",
           legend = c(
             "NHHSSI variance > HHSSI variance",
             "NHHSSI variance <= HHSSI variance"),
           col = c("blue", "red"),
           pch = 19,
           bty = "n")
  }






  idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_bias_normal01_simpred_hhssi1[idx, k],
      delta_bias_normal01_simpred_nhhssi1[idx, k],
      xlab = "HHSSI delta bias",
      ylab = "NHHSSI delta bias",
      main = paste("Overall: Outcome", k),
      xlim = range(c(delta_bias_normal01_simpred_hhssi1[idx,k],
                     delta_bias_normal01_simpred_nhhssi1[idx,k])),
      ylim = range(c(delta_bias_normal01_simpred_hhssi1[idx,k],
                     delta_bias_normal01_simpred_nhhssi1[idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }


  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_bias_normal01_simpred_hhssi1[!idx, k],
      delta_bias_normal01_simpred_nhhssi1[!idx, k],
      xlab = "HHSSI delta bias",
      ylab = "NHHSSI delta bias",
      main = paste("Overall: Outcome", k),
      xlim = range(c(delta_bias_normal01_simpred_hhssi1[!idx,k],
                     delta_bias_normal01_simpred_nhhssi1[!idx,k])),
      ylim = range(c(delta_bias_normal01_simpred_hhssi1[!idx,k],
                     delta_bias_normal01_simpred_nhhssi1[!idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }






  idx <- delta_var_mean_simpred_hhssi1[c(1:3, 6:8)] < delta_var_mean_simpred_nhhssi1[c(1:3, 6:8)]
  # no active delta with nhhssi var bigger than hhssi

  # inactive delta
  idx <- delta_var_mean_simpred_hhssi1[-c(1:3, 6:8)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 6:8)]
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8), ][idx, k],
      delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8), ][idx, k],
      xlab = "HHSSI delta bias",
      ylab = "NHHSSI delta bias",
      main = paste("Outcome", k),
      xlim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8),][idx, k],
                     delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8),][idx, k])),
      ylim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8),][idx, k],
                     delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8),][idx, k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }


  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8), ][!idx, k],
      delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8), ][!idx, k],
      xlab = "HHSSI delta bias",
      ylab = "NHHSSI delta bias",
      main = paste("Outcome", k),
      xlim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8),][!idx, k],
                     delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8),][!idx, k])),
      ylim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 6:8),][!idx, k],
                     delta_bias_normal01_simpred_nhhssi1[-c(1:3, 6:8),][!idx, k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }



  ######################################################################################

  idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
  # delta_idx <- which(!idx)
  delta_idx <- which(idx)
  # delta_idx <- 1:748
  beta_idx <- ceiling(delta_idx / M)
  gamma_idx <- ((delta_idx - 1) %% M) + 1
  data.frame(
    delta_index = delta_idx,
    beta_index = beta_idx,
    gamma_index = gamma_idx
  )



  length(beta_idx)
  length(gamma_idx)

  unique_beta_idx <- unique(beta_idx) # they are all inactive beta for NHHSSI>HHSSI var
  unique_gamma_idx <- unique(gamma_idx) # they are all inactive gamma for NHHSSI>HHSSI var


  length(unique_beta_idx)
  length(unique_gamma_idx)
  ###########################

  # # this diag_dat for inactive with nhhssi smaller hhssi only
  # # Keep only terms with inactive beta and gamma indices
  # keep <- !(beta_idx %in% c(1, 2) |
  #             gamma_idx %in% c(1, 2, 3))
  # delta_idx <- delta_idx[keep]
  # beta_idx  <- beta_idx[keep]
  # gamma_idx <- gamma_idx[keep]
  # diag_dat <- data.frame(
  #   delta = delta_idx,
  #   beta = beta_idx,
  #   gamma = gamma_idx
  # )


  diag_dat <- data.frame(
    delta = delta_idx,
    beta  = beta_idx,
    gamma = gamma_idx
  )
  beta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$theta_update[, 2:21, ]
  gamma_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$theta_update[, 22:26, ]
  delta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$theta_update[, 27:126, ]


  beta_mean_hhssi <- apply(beta_normal01_simpred_hhssi1, c(2, 3), mean)
  gamma_mean_hhssi <- apply(gamma_normal01_simpred_hhssi1, c(2, 3), mean)
  delta_mean_hhssi <- apply(delta_normal01_simpred_hhssi1, c(2, 3), mean)

  diag_dat$delta_est1 <- delta_mean_hhssi[delta_idx, 1]
  diag_dat$delta_est2 <- delta_mean_hhssi[delta_idx, 2]
  diag_dat$delta_est3 <- delta_mean_hhssi[delta_idx, 3]
  diag_dat$delta_est4 <- delta_mean_hhssi[delta_idx, 4]

  diag_dat$beta_est1 <- beta_mean_hhssi[beta_idx, 1]
  diag_dat$beta_est2 <- beta_mean_hhssi[beta_idx, 2]
  diag_dat$beta_est3 <- beta_mean_hhssi[beta_idx, 3]
  diag_dat$beta_est4 <- beta_mean_hhssi[beta_idx, 4]


  diag_dat$gamma_est1 <- gamma_mean_hhssi[gamma_idx, 1]
  diag_dat$gamma_est2 <- gamma_mean_hhssi[gamma_idx, 2]
  diag_dat$gamma_est3 <- gamma_mean_hhssi[gamma_idx, 3]
  diag_dat$gamma_est4 <- gamma_mean_hhssi[gamma_idx, 4]

  dim(diag_dat)
  head(diag_dat)

  library(dplyr)
  library(tidyr)

  diag_long <- diag_dat %>%
    pivot_longer(
      cols = matches("^(delta|beta|gamma)_est[1-4]$"),
      names_to = c("parameter", "outcome"),
      names_pattern = "(delta|beta|gamma)_est([1-4])",
      values_to = "estimate"
    )

  diag_long <- data.frame(
    delta = rep(diag_dat$delta, 4),
    beta  = rep(diag_dat$beta, 4),
    gamma = rep(diag_dat$gamma, 4),
    outcome = rep(1:4, each = nrow(diag_dat)),
    delta_est = c(
      diag_dat$delta_est1,
      diag_dat$delta_est2,
      diag_dat$delta_est3,
      diag_dat$delta_est4
    ),
    beta_est = c(
      diag_dat$beta_est1,
      diag_dat$beta_est2,
      diag_dat$beta_est3,
      diag_dat$beta_est4
    ),
    gamma_est = c(
      diag_dat$gamma_est1,
      diag_dat$gamma_est2,
      diag_dat$gamma_est3,
      diag_dat$gamma_est4
    )
  )


  diag_long$beta_gamma_product <-
    diag_long$beta_est * diag_long$gamma_est
  par(mfrow=c(1,1))
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    tmp <- subset(diag_long, outcome == k)
    plot(
      tmp$beta_est * tmp$gamma_est,
      tmp$delta_est,
      xlab = expression(hat(beta) %*% hat(gamma)),
      ylab = expression(hat(delta)),
      main = paste("Outcome", k),
      pch = 16,
      cex = 0.6
    )
    abline(h = 0, v = 0, lty = 2)
  }


  ##########################################################################################

  # Posterior sd estimate of beta
  beta_sd_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 2:21, ], c(2,3), sd)
  beta_sd_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 2:21, ], c(2,3), sd)

  # Posterior mean estimate of beta
  beta_est_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 2:21, ], c(2,3), mean)
  beta_est_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 2:21, ], c(2,3), mean)

  # Bias for each outcome
  beta_bias_normal01_simpred_hhssi1 <- (beta_est_normal01_simpred_hhssi1 - sim_data_normal01_simpred1$beta_true)
  beta_bias_normal01_simpred_nhhssi1 <- (beta_est_normal01_simpred_nhhssi1 - sim_data_normal01_simpred1$beta_true)

  # prior var estimate of beta
  beta_var_normal01_simpred_hhssi1 <- apply((res_normal01_simpred_hhssi1$lambdasq_beta_update*
                                               res_normal01_simpred_hhssi1$tausq_beta_update), 2, mean)
  beta_var_normal01_simpred_nhhssi1 <- apply((res_normal01_simpred_nhhssi1$lambdasq_beta_update*
                                                res_normal01_simpred_nhhssi1$tausq_beta_update), 2, mean)


  # beta Estimates
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_est_normal01_simpred_hhssi1[, k],
      beta_est_normal01_simpred_nhhssi1[, k],
      xlab = "HHSSI beta estimates",
      ylab = "NHHSSI beta estimates",
      main = paste("Overall estimates: Outcome", k),
      xlim = range(c(beta_est_normal01_simpred_hhssi1[,k],
                     beta_est_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(beta_est_normal01_simpred_hhssi1[,k],
                     beta_est_normal01_simpred_nhhssi1[,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # active beta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_est_normal01_simpred_hhssi1[1:2, k],
      beta_est_normal01_simpred_nhhssi1[1:2, k],
      xlab = "HHSSI beta estimates",
      ylab = "NHHSSI beta estimates",
      main = paste("active: Outcome", k),
      xlim = range(c(beta_est_normal01_simpred_hhssi1[1:2,k],
                     beta_est_normal01_simpred_nhhssi1[1:2,k])),
      ylim = range(c(beta_est_normal01_simpred_hhssi1[1:2,k],
                     beta_est_normal01_simpred_nhhssi1[1:2,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive beta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_est_normal01_simpred_hhssi1[-c(1:2), k],
      beta_est_normal01_simpred_nhhssi1[-c(1:2), k],
      xlab = "HHSSI beta estimates",
      ylab = "NHHSSI beta estimates",
      main = paste("inactive: Outcome", k),
      xlim = range(c(beta_est_normal01_simpred_hhssi1[-c(1:2),k],
                     beta_est_normal01_simpred_nhhssi1[-c(1:2),k])),
      ylim = range(c(beta_est_normal01_simpred_hhssi1[-c(1:2),k],
                     beta_est_normal01_simpred_nhhssi1[-c(1:2),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }







  # all the indices are all inactive beta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_est_normal01_simpred_hhssi1[unique_beta_idx, k],
      beta_est_normal01_simpred_nhhssi1[unique_beta_idx, k],
      xlab = "HHSSI beta estimates",
      ylab = "NHHSSI beta estimates",
      main = paste("Outcome", k),
      xlim = range(c(beta_est_normal01_simpred_hhssi1[unique_beta_idx,k],
                     beta_est_normal01_simpred_nhhssi1[unique_beta_idx,k])),
      ylim = range(c(beta_est_normal01_simpred_hhssi1[unique_beta_idx,k],
                     beta_est_normal01_simpred_nhhssi1[unique_beta_idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }




  ##########################################################################################
  # beta Bias
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_bias_normal01_simpred_hhssi1[, k],
      beta_bias_normal01_simpred_nhhssi1[, k],
      xlab = "HHSSI beta bias",
      ylab = "NHHSSI beta bias",
      main = paste("Overall bias: Outcome", k),
      xlim = range(c(beta_bias_normal01_simpred_hhssi1[,k],
                     beta_bias_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(beta_bias_normal01_simpred_hhssi1[,k],
                     beta_bias_normal01_simpred_nhhssi1[,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # active beta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_bias_normal01_simpred_hhssi1[c(1:2), k],
      beta_bias_normal01_simpred_nhhssi1[c(1:2), k],
      xlab = "HHSSI beta bias",
      ylab = "NHHSSI beta bias",
      main = paste("active: Outcome", k),
      xlim = range(c(beta_bias_normal01_simpred_hhssi1[c(1:2),k],
                     beta_bias_normal01_simpred_nhhssi1[c(1:2),k])),
      ylim = range(c(beta_bias_normal01_simpred_hhssi1[c(1:2),k],
                     beta_bias_normal01_simpred_nhhssi1[c(1:2),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive beta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_bias_normal01_simpred_hhssi1[-c(1:2), k],
      beta_bias_normal01_simpred_nhhssi1[-c(1:2), k],
      xlab = "HHSSI beta bias",
      ylab = "NHHSSI beta bias",
      main = paste("inactive: Outcome", k),
      xlim = range(c(beta_bias_normal01_simpred_hhssi1[-c(1:2),k],
                     beta_bias_normal01_simpred_nhhssi1[-c(1:2),k])),
      ylim = range(c(beta_bias_normal01_simpred_hhssi1[-c(1:2),k],
                     beta_bias_normal01_simpred_nhhssi1[-c(1:2),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }



  # all the indices are all inactive beta
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      beta_bias_normal01_simpred_hhssi1[unique_beta_idx, k],
      beta_bias_normal01_simpred_nhhssi1[unique_beta_idx, k],
      xlab = "HHSSI beta bias",
      ylab = "NHHSSI beta bias",
      main = paste("Overall: Outcome", k),
      xlim = range(c(beta_bias_normal01_simpred_hhssi1[unique_beta_idx,k],
                     beta_bias_normal01_simpred_nhhssi1[unique_beta_idx,k])),
      ylim = range(c(beta_bias_normal01_simpred_hhssi1[unique_beta_idx,k],
                     beta_bias_normal01_simpred_nhhssi1[unique_beta_idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }




  ######################################################################################

  # beta prior variance estimate
  # overall
  par(mfrow = c(1,1))
  plot(beta_var_normal01_simpred_hhssi1,
       beta_var_normal01_simpred_nhhssi1,
       xlab = "HHSSI beta var",
       ylab = "NHHSSI beta var",
       xlim = range(c(beta_var_normal01_simpred_hhssi1,
                      beta_var_normal01_simpred_nhhssi1)),
       ylim = range(c(beta_var_normal01_simpred_hhssi1,
                      beta_var_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active beta
  plot(beta_var_normal01_simpred_hhssi1[c(1:2)],
       beta_var_normal01_simpred_nhhssi1[c(1:2)],
       xlab = "HHSSI beta var",
       ylab = "NHHSSI beta var",
       xlim = range(c(beta_var_normal01_simpred_hhssi1[c(1:2)],
                      beta_var_normal01_simpred_nhhssi1[c(1:2)])),
       ylim = range(c(beta_var_normal01_simpred_hhssi1[c(1:2)],
                      beta_var_normal01_simpred_nhhssi1[c(1:2)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive beta
  plot(beta_var_normal01_simpred_hhssi1[-c(1:2)],
       beta_var_normal01_simpred_nhhssi1[-c(1:2)],
       xlab = "HHSSI beta var",
       ylab = "NHHSSI beta var",
       xlim = range(c(beta_var_normal01_simpred_hhssi1[-c(1:2)],
                      beta_var_normal01_simpred_nhhssi1[-c(1:2)])),
       ylim = range(c(beta_var_normal01_simpred_hhssi1[-c(1:2)],
                      beta_var_normal01_simpred_nhhssi1[-c(1:2)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # all are infact inactive beta
  plot(beta_var_normal01_simpred_hhssi1[unique_beta_idx],
       beta_var_normal01_simpred_nhhssi1[unique_beta_idx],
       xlab = "HHSSI beta var",
       ylab = "NHHSSI beta var",
       xlim = range(c(beta_var_normal01_simpred_hhssi1[unique_beta_idx],
                      beta_var_normal01_simpred_nhhssi1[unique_beta_idx])),
       ylim = range(c(beta_var_normal01_simpred_hhssi1[unique_beta_idx],
                      beta_var_normal01_simpred_nhhssi1[unique_beta_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")





  ######################################################################################
  # Posterior sd estimate of gamma
  gamma_sd_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 22:26, ], c(2,3), sd)
  gamma_sd_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 22:26, ], c(2,3), sd)

  # Posterior mean estimate of gamma
  gamma_est_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 22:26, ], c(2,3), mean)
  gamma_est_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 22:26, ], c(2,3), mean)

  # Bias for each outcome
  gamma_bias_normal01_simpred_hhssi1 <- (gamma_est_normal01_simpred_hhssi1 - sim_data_normal01_simpred1$gamma_true)
  gamma_bias_normal01_simpred_nhhssi1 <- (gamma_est_normal01_simpred_nhhssi1 - sim_data_normal01_simpred1$gamma_true)

  # prior var estimate of gamma
  gamma_var_normal01_simpred_hhssi1 <- apply((res_normal01_simpred_hhssi1$lambdasq_gamma_update*
                                                res_normal01_simpred_hhssi1$tausq_gamma_update), 2, mean)
  gamma_var_normal01_simpred_nhhssi1 <- apply((res_normal01_simpred_nhhssi1$lambdasq_gamma_update*
                                                 res_normal01_simpred_nhhssi1$tausq_gamma_update), 2, mean)


  # gamma Estimates
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_est_normal01_simpred_hhssi1[, k],
      gamma_est_normal01_simpred_nhhssi1[, k],
      xlab = "HHSSI gamma estimates",
      ylab = "NHHSSI gamma estimates",
      main = paste("Overall estimates: Outcome", k),
      xlim = range(c(gamma_est_normal01_simpred_hhssi1[,k],
                     gamma_est_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(gamma_est_normal01_simpred_hhssi1[,k],
                     gamma_est_normal01_simpred_nhhssi1[,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # active gamma
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_est_normal01_simpred_hhssi1[1:3, k],
      gamma_est_normal01_simpred_nhhssi1[1:3, k],
      xlab = "HHSSI gamma estimates",
      ylab = "NHHSSI gamma estimates",
      main = paste("active: Outcome", k),
      xlim = range(c(gamma_est_normal01_simpred_hhssi1[1:3,k],
                     gamma_est_normal01_simpred_nhhssi1[1:3,k])),
      ylim = range(c(gamma_est_normal01_simpred_hhssi1[1:3,k],
                     gamma_est_normal01_simpred_nhhssi1[1:3,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive gamma
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_est_normal01_simpred_hhssi1[-c(1:3), k],
      gamma_est_normal01_simpred_nhhssi1[-c(1:3), k],
      xlab = "HHSSI gamma estimates",
      ylab = "NHHSSI gamma estimates",
      main = paste("inactive: Outcome", k),
      xlim = range(c(gamma_est_normal01_simpred_hhssi1[-c(1:3),k],
                     gamma_est_normal01_simpred_nhhssi1[-c(1:3),k])),
      ylim = range(c(gamma_est_normal01_simpred_hhssi1[-c(1:3),k],
                     gamma_est_normal01_simpred_nhhssi1[-c(1:3),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }







  # overall (infact all are inactive gamma)
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_est_normal01_simpred_hhssi1[unique_gamma_idx, k],
      gamma_est_normal01_simpred_nhhssi1[unique_gamma_idx, k],
      xlab = "HHSSI gamma estimates",
      ylab = "NHHSSI gamma estimates",
      main = paste("Outcome", k),
      xlim = range(c(gamma_est_normal01_simpred_hhssi1[unique_gamma_idx,k],
                     gamma_est_normal01_simpred_nhhssi1[unique_gamma_idx,k])),
      ylim = range(c(gamma_est_normal01_simpred_hhssi1[unique_gamma_idx,k],
                     gamma_est_normal01_simpred_nhhssi1[unique_gamma_idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }




  ##########################################################################################
  # gamma Bias
  # overall
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_bias_normal01_simpred_hhssi1[, k],
      gamma_bias_normal01_simpred_nhhssi1[, k],
      xlab = "HHSSI gamma bias",
      ylab = "NHHSSI gamma bias",
      main = paste("Overall bias: Outcome", k),
      xlim = range(c(gamma_bias_normal01_simpred_hhssi1[,k],
                     gamma_bias_normal01_simpred_nhhssi1[,k])),
      ylim = range(c(gamma_bias_normal01_simpred_hhssi1[,k],
                     gamma_bias_normal01_simpred_nhhssi1[,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # active gamma
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_bias_normal01_simpred_hhssi1[c(1:3), k],
      gamma_bias_normal01_simpred_nhhssi1[c(1:3), k],
      xlab = "HHSSI gamma bias",
      ylab = "NHHSSI gamma bias",
      main = paste("active: Outcome", k),
      xlim = range(c(gamma_bias_normal01_simpred_hhssi1[c(1:3),k],
                     gamma_bias_normal01_simpred_nhhssi1[c(1:3),k])),
      ylim = range(c(gamma_bias_normal01_simpred_hhssi1[c(1:3),k],
                     gamma_bias_normal01_simpred_nhhssi1[c(1:3),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }

  # inactive gamma
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_bias_normal01_simpred_hhssi1[-c(1:3), k],
      gamma_bias_normal01_simpred_nhhssi1[-c(1:3), k],
      xlab = "HHSSI gamma bias",
      ylab = "NHHSSI gamma bias",
      main = paste("inactive: Outcome", k),
      xlim = range(c(gamma_bias_normal01_simpred_hhssi1[-c(1:3),k],
                     gamma_bias_normal01_simpred_nhhssi1[-c(1:3),k])),
      ylim = range(c(gamma_bias_normal01_simpred_hhssi1[-c(1:3),k],
                     gamma_bias_normal01_simpred_nhhssi1[-c(1:3),k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }



  # overall (infact all are inactive gamma)
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    plot(
      gamma_bias_normal01_simpred_hhssi1[unique_gamma_idx, k],
      gamma_bias_normal01_simpred_nhhssi1[unique_gamma_idx, k],
      xlab = "HHSSI gamma bias",
      ylab = "NHHSSI gamma bias",
      main = paste("Overall: Outcome", k),
      xlim = range(c(gamma_bias_normal01_simpred_hhssi1[unique_gamma_idx,k],
                     gamma_bias_normal01_simpred_nhhssi1[unique_gamma_idx,k])),
      ylim = range(c(gamma_bias_normal01_simpred_hhssi1[unique_gamma_idx,k],
                     gamma_bias_normal01_simpred_nhhssi1[unique_gamma_idx,k])))
    abline(0,1, lty = 2, lwd = 1, col = "red")
  }




  ######################################################################################
  # gamma prior variance estimate
  # overall
  par(mfrow = c(1,1))
  plot(gamma_var_normal01_simpred_hhssi1,
       gamma_var_normal01_simpred_nhhssi1,
       xlab = "HHSSI gamma var",
       ylab = "NHHSSI gamma var",
       xlim = range(c(gamma_var_normal01_simpred_hhssi1,
                      gamma_var_normal01_simpred_nhhssi1)),
       ylim = range(c(gamma_var_normal01_simpred_hhssi1,
                      gamma_var_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active gamma
  plot(gamma_var_normal01_simpred_hhssi1[c(1:3)],
       gamma_var_normal01_simpred_nhhssi1[c(1:3)],
       xlab = "HHSSI gamma var",
       ylab = "NHHSSI gamma var",
       xlim = range(c(gamma_var_normal01_simpred_hhssi1[c(1:3)],
                      gamma_var_normal01_simpred_nhhssi1[c(1:3)])),
       ylim = range(c(gamma_var_normal01_simpred_hhssi1[c(1:3)],
                      gamma_var_normal01_simpred_nhhssi1[c(1:3)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive gamma
  plot(gamma_var_normal01_simpred_hhssi1[-c(1:3)],
       gamma_var_normal01_simpred_nhhssi1[-c(1:3)],
       xlab = "HHSSI gamma var",
       ylab = "NHHSSI gamma var",
       xlim = range(c(gamma_var_normal01_simpred_hhssi1[-c(1:3)],
                      gamma_var_normal01_simpred_nhhssi1[-c(1:3)])),
       ylim = range(c(gamma_var_normal01_simpred_hhssi1[-c(1:3)],
                      gamma_var_normal01_simpred_nhhssi1[-c(1:3)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # overall (infact all are inactive gamma)
  plot(gamma_var_normal01_simpred_hhssi1[unique_gamma_idx],
       gamma_var_normal01_simpred_nhhssi1[unique_gamma_idx],
       xlab = "HHSSI gamma var",
       ylab = "NHHSSI gamma var",
       xlim = range(c(gamma_var_normal01_simpred_hhssi1[unique_gamma_idx],
                      gamma_var_normal01_simpred_nhhssi1[unique_gamma_idx])),
       ylim = range(c(gamma_var_normal01_simpred_hhssi1[unique_gamma_idx],
                      gamma_var_normal01_simpred_nhhssi1[unique_gamma_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  ######################################################################################

  # lambdasq terms
  lambdasq_beta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update, 2, mean)
  lambdasq_beta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update, 2, mean)

  lambdasq_gamma_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
  lambdasq_gamma_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update, 2, mean)

  lambdasq_delta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_delta_update, 2, mean)
  lambdasq_delta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_delta_update, 2, mean)

  ######################################################################################

  # lambdasq_beta estimate
  # overall
  par(mfrow = c(1,1))
  plot(lambdasq_beta_normal01_simpred_hhssi1,
       lambdasq_beta_normal01_simpred_nhhssi1,
       xlab = "HHSSI lambdasq beta",
       ylab = "NHHSSI lambdasq beta",
       xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1,
                      lambdasq_beta_normal01_simpred_nhhssi1)),
       ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1,
                      lambdasq_beta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active beta
  plot(lambdasq_beta_normal01_simpred_hhssi1[c(1:2)],
       lambdasq_beta_normal01_simpred_nhhssi1[c(1:2)],
       xlab = "HHSSI lambdasq beta",
       ylab = "NHHSSI lambdasq beta",
       xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1[c(1:2)],
                      lambdasq_beta_normal01_simpred_nhhssi1[c(1:2)])),
       ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1[c(1:2)],
                      lambdasq_beta_normal01_simpred_nhhssi1[c(1:2)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive beta
  plot(lambdasq_beta_normal01_simpred_hhssi1[-c(1:2)],
       lambdasq_beta_normal01_simpred_nhhssi1[-c(1:2)],
       xlab = "HHSSI lambdasq beta",
       ylab = "NHHSSI lambdasq beta",
       xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1[-c(1:2)],
                      lambdasq_beta_normal01_simpred_nhhssi1[-c(1:2)])),
       ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1[-c(1:2)],
                      lambdasq_beta_normal01_simpred_nhhssi1[-c(1:2)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # all are infact inactive beta
  plot(lambdasq_beta_normal01_simpred_hhssi1[unique_beta_idx],
       lambdasq_beta_normal01_simpred_nhhssi1[unique_beta_idx],
       xlab = "HHSSI lambdasq beta",
       ylab = "NHHSSI lambdasq beta",
       xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1[unique_beta_idx],
                      lambdasq_beta_normal01_simpred_nhhssi1[unique_beta_idx])),
       ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1[unique_beta_idx],
                      lambdasq_beta_normal01_simpred_nhhssi1[unique_beta_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")




  ######################################################################################

  # lambdasq_gamma estimate
  # overall
  par(mfrow = c(1,1))
  plot(lambdasq_gamma_normal01_simpred_hhssi1,
       lambdasq_gamma_normal01_simpred_nhhssi1,
       xlab = "HHSSI lambdasq gamma",
       ylab = "NHHSSI lambdasq gamma",
       xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1,
                      lambdasq_gamma_normal01_simpred_nhhssi1)),
       ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1,
                      lambdasq_gamma_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active gamma
  plot(lambdasq_gamma_normal01_simpred_hhssi1[c(1:3)],
       lambdasq_gamma_normal01_simpred_nhhssi1[c(1:3)],
       xlab = "HHSSI lambdasq gamma",
       ylab = "NHHSSI lambdasq gamma",
       xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[c(1:3)],
                      lambdasq_gamma_normal01_simpred_nhhssi1[c(1:3)])),
       ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[c(1:3)],
                      lambdasq_gamma_normal01_simpred_nhhssi1[c(1:3)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive gamma
  plot(lambdasq_gamma_normal01_simpred_hhssi1[-c(1:3)],
       lambdasq_gamma_normal01_simpred_nhhssi1[-c(1:3)],
       xlab = "HHSSI lambdasq gamma",
       ylab = "NHHSSI lambdasq gamma",
       xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[-c(1:3)],
                      lambdasq_gamma_normal01_simpred_nhhssi1[-c(1:3)])),
       ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[-c(1:3)],
                      lambdasq_gamma_normal01_simpred_nhhssi1[-c(1:3)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # all are infact inactive gamma
  plot(lambdasq_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
       lambdasq_gamma_normal01_simpred_nhhssi1[unique_gamma_idx],
       xlab = "HHSSI lambdasq gamma",
       ylab = "NHHSSI lambdasq gamma",
       xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                      lambdasq_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])),
       ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                      lambdasq_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")




  ######################################################################################

  # lambdasq_delta estimate
  # overall
  par(mfrow = c(1,1))
  plot(lambdasq_delta_normal01_simpred_hhssi1,
       lambdasq_delta_normal01_simpred_nhhssi1,
       xlab = "HHSSI lambdasq delta",
       ylab = "NHHSSI lambdasq delta",
       xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1,
                      lambdasq_delta_normal01_simpred_nhhssi1)),
       ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1,
                      lambdasq_delta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active delta
  plot(lambdasq_delta_normal01_simpred_hhssi1[c(1:3, 6:8)],
       lambdasq_delta_normal01_simpred_nhhssi1[c(1:3, 6:8)],
       xlab = "HHSSI lambdasq delta",
       ylab = "NHHSSI lambdasq delta",
       xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1[c(1:3, 6:8)],
                      lambdasq_delta_normal01_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1[c(1:3, 6:8)],
                      lambdasq_delta_normal01_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive delta
  plot(lambdasq_delta_normal01_simpred_hhssi1[-c(1:3, 6:8)],
       lambdasq_delta_normal01_simpred_nhhssi1[-c(1:3, 6:8)],
       xlab = "HHSSI lambdasq delta",
       ylab = "NHHSSI lambdasq delta",
       xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1[-c(1:3, 6:8)],
                      lambdasq_delta_normal01_simpred_nhhssi1[-c(1:3, 6:8)])),
       ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1[-c(1:3, 6:8)],
                      lambdasq_delta_normal01_simpred_nhhssi1[-c(1:3, 6:8)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # all are infact inactive delta
  plot(lambdasq_delta_normal01_simpred_hhssi1[delta_idx],
       lambdasq_delta_normal01_simpred_nhhssi1[delta_idx],
       xlab = "HHSSI lambdasq delta",
       ylab = "NHHSSI lambdasq delta",
       xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1[delta_idx],
                      lambdasq_delta_normal01_simpred_nhhssi1[delta_idx])),
       ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1[delta_idx],
                      lambdasq_delta_normal01_simpred_nhhssi1[delta_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")









  ######################################################################################

  # psi terms
  psi_beta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$psi_beta_update, 2, mean)
  psi_beta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$psi_beta_update, 2, mean)

  psi_gamma_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$psi_gamma_update, 2, mean)
  psi_gamma_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$psi_gamma_update, 2, mean)

  psi_delta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$psi_delta_update, 2, mean)
  psi_delta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$psi_delta_update, 2, mean)

  ######################################################################################

  # psi_beta estimate
  # overall
  par(mfrow = c(1,1))
  plot(psi_beta_normal01_simpred_hhssi1,
       psi_beta_normal01_simpred_nhhssi1,
       xlab = "HHSSI psi beta",
       ylab = "NHHSSI psi beta",
       xlim = range(c(psi_beta_normal01_simpred_hhssi1,
                      psi_beta_normal01_simpred_nhhssi1)),
       ylim = range(c(psi_beta_normal01_simpred_hhssi1,
                      psi_beta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active beta
  plot(psi_beta_normal01_simpred_hhssi1[c(1:2)],
       psi_beta_normal01_simpred_nhhssi1[c(1:2)],
       xlab = "HHSSI psi beta",
       ylab = "NHHSSI psi beta",
       xlim = range(c(psi_beta_normal01_simpred_hhssi1[c(1:2)],
                      psi_beta_normal01_simpred_nhhssi1[c(1:2)])),
       ylim = range(c(psi_beta_normal01_simpred_hhssi1[c(1:2)],
                      psi_beta_normal01_simpred_nhhssi1[c(1:2)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive beta
  plot(psi_beta_normal01_simpred_hhssi1[-c(1:2)],
       psi_beta_normal01_simpred_nhhssi1[-c(1:2)],
       xlab = "HHSSI psi beta",
       ylab = "NHHSSI psi beta",
       xlim = range(c(psi_beta_normal01_simpred_hhssi1[-c(1:2)],
                      psi_beta_normal01_simpred_nhhssi1[-c(1:2)])),
       ylim = range(c(psi_beta_normal01_simpred_hhssi1[-c(1:2)],
                      psi_beta_normal01_simpred_nhhssi1[-c(1:2)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # all are infact inactive beta
  plot(psi_beta_normal01_simpred_hhssi1[unique_beta_idx],
       psi_beta_normal01_simpred_nhhssi1[unique_beta_idx],
       xlab = "HHSSI psi beta",
       ylab = "NHHSSI psi beta",
       xlim = range(c(psi_beta_normal01_simpred_hhssi1[unique_beta_idx],
                      psi_beta_normal01_simpred_nhhssi1[unique_beta_idx])),
       ylim = range(c(psi_beta_normal01_simpred_hhssi1[unique_beta_idx],
                      psi_beta_normal01_simpred_nhhssi1[unique_beta_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")




  ######################################################################################

  # psi_gamma estimate
  # overall
  par(mfrow = c(1,1))
  plot(psi_gamma_normal01_simpred_hhssi1,
       psi_gamma_normal01_simpred_nhhssi1,
       xlab = "HHSSI psi gamma",
       ylab = "NHHSSI psi gamma",
       xlim = range(c(psi_gamma_normal01_simpred_hhssi1,
                      psi_gamma_normal01_simpred_nhhssi1)),
       ylim = range(c(psi_gamma_normal01_simpred_hhssi1,
                      psi_gamma_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active gamma
  plot(psi_gamma_normal01_simpred_hhssi1[c(1:3)],
       psi_gamma_normal01_simpred_nhhssi1[c(1:3)],
       xlab = "HHSSI psi gamma",
       ylab = "NHHSSI psi gamma",
       xlim = range(c(psi_gamma_normal01_simpred_hhssi1[c(1:3)],
                      psi_gamma_normal01_simpred_nhhssi1[c(1:3)])),
       ylim = range(c(psi_gamma_normal01_simpred_hhssi1[c(1:3)],
                      psi_gamma_normal01_simpred_nhhssi1[c(1:3)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive gamma
  plot(psi_gamma_normal01_simpred_hhssi1[-c(1:3)],
       psi_gamma_normal01_simpred_nhhssi1[-c(1:3)],
       xlab = "HHSSI psi gamma",
       ylab = "NHHSSI psi gamma",
       xlim = range(c(psi_gamma_normal01_simpred_hhssi1[-c(1:3)],
                      psi_gamma_normal01_simpred_nhhssi1[-c(1:3)])),
       ylim = range(c(psi_gamma_normal01_simpred_hhssi1[-c(1:3)],
                      psi_gamma_normal01_simpred_nhhssi1[-c(1:3)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # all are infact inactive gamma
  plot(psi_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
       psi_gamma_normal01_simpred_nhhssi1[unique_gamma_idx],
       xlab = "HHSSI psi gamma",
       ylab = "NHHSSI psi gamma",
       xlim = range(c(psi_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                      psi_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])),
       ylim = range(c(psi_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                      psi_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")




  ######################################################################################

  # psi_delta estimate
  # overall
  par(mfrow = c(1,1))
  plot(psi_delta_normal01_simpred_hhssi1,
       psi_delta_normal01_simpred_nhhssi1,
       xlab = "HHSSI psi delta",
       ylab = "NHHSSI psi delta",
       xlim = range(c(psi_delta_normal01_simpred_hhssi1,
                      psi_delta_normal01_simpred_nhhssi1)),
       ylim = range(c(psi_delta_normal01_simpred_hhssi1,
                      psi_delta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # active delta
  plot(psi_delta_normal01_simpred_hhssi1[c(1:3, 6:8)],
       psi_delta_normal01_simpred_nhhssi1[c(1:3, 6:8)],
       xlab = "HHSSI psi delta",
       ylab = "NHHSSI psi delta",
       xlim = range(c(psi_delta_normal01_simpred_hhssi1[c(1:3, 6:8)],
                      psi_delta_normal01_simpred_nhhssi1[c(1:3, 6:8)])),
       ylim = range(c(psi_delta_normal01_simpred_hhssi1[c(1:3, 6:8)],
                      psi_delta_normal01_simpred_nhhssi1[c(1:3, 6:8)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")

  # inactive delta
  plot(psi_delta_normal01_simpred_hhssi1[-c(1:3, 6:8)],
       psi_delta_normal01_simpred_nhhssi1[-c(1:3, 6:8)],
       xlab = "HHSSI psi delta",
       ylab = "NHHSSI psi delta",
       xlim = range(c(psi_delta_normal01_simpred_hhssi1[-c(1:3, 6:8)],
                      psi_delta_normal01_simpred_nhhssi1[-c(1:3, 6:8)])),
       ylim = range(c(psi_delta_normal01_simpred_hhssi1[-c(1:3, 6:8)],
                      psi_delta_normal01_simpred_nhhssi1[-c(1:3, 6:8)])))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # all are infact inactive delta
  plot(psi_delta_normal01_simpred_hhssi1[delta_idx],
       psi_delta_normal01_simpred_nhhssi1[delta_idx],
       xlab = "HHSSI psi delta",
       ylab = "NHHSSI psi delta",
       xlim = range(c(psi_delta_normal01_simpred_hhssi1[delta_idx],
                      psi_delta_normal01_simpred_nhhssi1[delta_idx])),
       ylim = range(c(psi_delta_normal01_simpred_hhssi1[delta_idx],
                      psi_delta_normal01_simpred_nhhssi1[delta_idx])))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  ######################################################################################

  # tausq terms
  tausq_beta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_beta_update
  tausq_beta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_beta_update

  tausq_gamma_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_gamma_update
  tausq_gamma_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_gamma_update

  tausq_delta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update
  tausq_delta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update

  ######################################################################################

  # tausq_beta estimate
  plot(tausq_beta_normal01_simpred_hhssi1,
       tausq_beta_normal01_simpred_nhhssi1,
       xlab = "HHSSI tausq beta",
       ylab = "NHHSSI tausq beta",
       xlim = range(c(tausq_beta_normal01_simpred_hhssi1,
                      tausq_beta_normal01_simpred_nhhssi1)),
       ylim = range(c(tausq_beta_normal01_simpred_hhssi1,
                      tausq_beta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # tausq_gamma estimate
  plot(tausq_gamma_normal01_simpred_hhssi1,
       tausq_gamma_normal01_simpred_nhhssi1,
       xlab = "HHSSI tausq gamma",
       ylab = "NHHSSI tausq gamma",
       xlim = range(c(tausq_gamma_normal01_simpred_hhssi1,
                      tausq_gamma_normal01_simpred_nhhssi1)),
       ylim = range(c(tausq_gamma_normal01_simpred_hhssi1,
                      tausq_gamma_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # tausq_delta estimate
  plot(tausq_delta_normal01_simpred_hhssi1,
       tausq_delta_normal01_simpred_nhhssi1,
       xlab = "HHSSI tausq delta",
       ylab = "NHHSSI tausq delta",
       xlim = range(c(tausq_delta_normal01_simpred_hhssi1,
                      tausq_delta_normal01_simpred_nhhssi1)),
       ylim = range(c(tausq_delta_normal01_simpred_hhssi1,
                      tausq_delta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")





  ######################################################################################

  # xi terms
  xi_beta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$xi_beta_update
  xi_beta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$xi_beta_update

  xi_gamma_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$xi_gamma_update
  xi_gamma_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$xi_gamma_update

  xi_delta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$xi_delta_update
  xi_delta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$xi_delta_update

  ######################################################################################

  # xi_beta estimate
  plot(xi_beta_normal01_simpred_hhssi1,
       xi_beta_normal01_simpred_nhhssi1,
       xlab = "HHSSI xi beta",
       ylab = "NHHSSI xi beta",
       xlim = range(c(xi_beta_normal01_simpred_hhssi1,
                      xi_beta_normal01_simpred_nhhssi1)),
       ylim = range(c(xi_beta_normal01_simpred_hhssi1,
                      xi_beta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  # xi_gamma estimate
  plot(xi_gamma_normal01_simpred_hhssi1,
       xi_gamma_normal01_simpred_nhhssi1,
       xlab = "HHSSI xi gamma",
       ylab = "NHHSSI xi gamma",
       xlim = range(c(xi_gamma_normal01_simpred_hhssi1,
                      xi_gamma_normal01_simpred_nhhssi1)),
       ylim = range(c(xi_gamma_normal01_simpred_hhssi1,
                      xi_gamma_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")



  # xi_delta estimate
  plot(xi_delta_normal01_simpred_hhssi1,
       xi_delta_normal01_simpred_nhhssi1,
       xlab = "HHSSI xi delta",
       ylab = "NHHSSI xi delta",
       xlim = range(c(xi_delta_normal01_simpred_hhssi1,
                      xi_delta_normal01_simpred_nhhssi1)),
       ylim = range(c(xi_delta_normal01_simpred_hhssi1,
                      xi_delta_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


  #############################################################################################
  # prior var og=f beta, gamma
  beta_var_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update *
                                          res_normal01_simpred_hhssi1$tausq_beta_update, 2, mean)
  beta_var_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update *
                                           res_normal01_simpred_nhhssi1$tausq_beta_update, 2, mean)


  gamma_var_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update *
                                           res_normal01_simpred_hhssi1$tausq_gamma_update, 2, mean)
  gamma_var_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update *
                                            res_normal01_simpred_nhhssi1$tausq_gamma_update, 2, mean)


  point_col_beta <- ifelse(seq_along(beta_var_mean_simpred_hhssi1) %in% 1:2, "blue", "red")
  plot(
    log(beta_var_mean_simpred_hhssi1),
    log(beta_var_mean_simpred_nhhssi1),
    pch = 19,
    col = point_col_beta,
    main= "log transformed",
    xlab = "HHSSI beta prior var",
    ylab = "NHHSSI beta prior var",
    xlim=range(c(log(beta_var_mean_simpred_hhssi1), log(beta_var_mean_simpred_nhhssi1))),
    ylim=range(c(log(beta_var_mean_simpred_hhssi1), log(beta_var_mean_simpred_nhhssi1)))
  )
  # Legend
  legend("topleft",
         legend = c(
           "active beta",
           "inactive"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0, 1, lty = 2, col="red")


  plot(
    beta_var_mean_simpred_hhssi1[1:2],
    beta_var_mean_simpred_nhhssi1[1:2],
    pch = 19,
    col = "blue",
    xlab = "HHSSI beta prior var",
    ylab = "NHHSSI beta prior var",
    xlim=range(c(beta_var_mean_simpred_hhssi1[1:2], beta_var_mean_simpred_nhhssi1[1:2])),
    ylim=range(c(beta_var_mean_simpred_hhssi1[1:2], beta_var_mean_simpred_nhhssi1[1:2]))
  )
  abline(0, 1, lty = 2, col="red")


  plot(
    beta_var_mean_simpred_hhssi1[-(1:2)],
    beta_var_mean_simpred_nhhssi1[-(1:2)],
    pch = 19,
    col = "red",
    xlab = "HHSSI beta prior var",
    ylab = "NHHSSI beta prior var",
    xlim=range(c(beta_var_mean_simpred_hhssi1[-(1:2)], beta_var_mean_simpred_nhhssi1[-(1:2)])),
    ylim=range(c(beta_var_mean_simpred_hhssi1[-(1:2)], beta_var_mean_simpred_nhhssi1[-(1:2)]))
  )
  abline(0, 1, lty = 2, col="red")






  point_col_gamma <- ifelse(seq_along(gamma_var_mean_simpred_hhssi1) %in% 1:3, "blue", "red")
  plot(
    log(gamma_var_mean_simpred_hhssi1),
    log(gamma_var_mean_simpred_nhhssi1),
    pch = 19,
    col = point_col_gamma,
    main= "log transformed",
    xlab = "HHSSI gamma prior var",
    ylab = "NHHSSI gamma prior var",
    xlim=range(c(log(gamma_var_mean_simpred_hhssi1), log(gamma_var_mean_simpred_nhhssi1))),
    ylim=range(c(log(gamma_var_mean_simpred_hhssi1), log(gamma_var_mean_simpred_nhhssi1)))
  )
  # Legend
  legend("topleft",
         legend = c(
           "active gamma",
           "inactive"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0, 1, lty = 2, col="red")


  plot(
    gamma_var_mean_simpred_hhssi1[1:3],
    gamma_var_mean_simpred_nhhssi1[1:3],
    pch = 19,
    col = "blue",
    xlab = "HHSSI gamma prior var",
    ylab = "NHHSSI gamma prior var",
    xlim=range(c(gamma_var_mean_simpred_hhssi1[1:3], gamma_var_mean_simpred_nhhssi1[1:3])),
    ylim=range(c(gamma_var_mean_simpred_hhssi1[1:3], gamma_var_mean_simpred_nhhssi1[1:3]))
  )
  abline(0, 1, lty = 2, col="red")


  plot(
    gamma_var_mean_simpred_hhssi1[-(1:3)],
    gamma_var_mean_simpred_nhhssi1[-(1:3)],
    pch = 19,
    col = "red",
    xlab = "HHSSI gamma prior var",
    ylab = "NHHSSI gamma prior var",
    xlim=range(c(gamma_var_mean_simpred_hhssi1[-(1:3)], gamma_var_mean_simpred_nhhssi1[-(1:3)])),
    ylim=range(c(gamma_var_mean_simpred_hhssi1[-(1:3)], gamma_var_mean_simpred_nhhssi1[-(1:3)]))
  )
  abline(0, 1, lty = 2, col="red")

  ##############################################################################################

  # ratio of prior variance estimate from NHHSSI and HHSSI vs lambdasq_beta*lambdasq_gamma

  beta_gamma_lam_mean_simpred_hhssi1
  delta_tau_simpred_hhssi1 <- (delta_lam_simpred_hhssi1 *
                                 tausq_delta_simpred_hhssi1)
  delta_tau_mean_simpred_hhssi1 <- apply(delta_tau_simpred_hhssi1, 2, mean)


  plot(
    log(beta_gamma_lam_mean_simpred_hhssi1),
    log(delta_tau_mean_simpred_hhssi1),
    pch = 19,
    main= "log transformed",
    ylab = expression(hat(lambda)[delta]^2 * hat(tau)[delta]^2),
    xlab = expression(hat(lambda)[beta]^2 * hat(lambda)[gamma]^2),
    xlim=range(c(log(beta_gamma_lam_mean_simpred_hhssi1), log(delta_tau_mean_simpred_hhssi1))),
    ylim=range(c(log(beta_gamma_lam_mean_simpred_hhssi1), log(delta_tau_mean_simpred_hhssi1)))
  )
  abline(0, 1, lty = 2, col="red")


  # Identify whether NHHSSI variance is greater than HHSSI variance
  idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
  # Colors
  point_col <- ifelse(idx, "blue", "red")
  # Common limits
  x_vals <- log(beta_gamma_lam_mean_simpred_hhssi1)
  y_vals <- log(delta_tau_mean_simpred_hhssi1)
  plot(x_vals, y_vals,
       pch = 19,
       col = point_col,
       main = "HHSSI prior variance structure",
       ylab = expression(log(hat(lambda)[delta]^2 * hat(tau)[delta]^2)),
       xlab = expression(log(hat(lambda)[beta]^2 * hat(lambda)[gamma]^2)),
       xlim = range(c(x_vals, y_vals)),
       ylim = range(c(x_vals, y_vals))
  )
  # 45-degree reference line
  abline(0, 1, lty = 2, col = "black")
  # Legend
  legend("topright",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n"
  )



  #################################################

  # Beta-gamma shrinkage term
  bg <- beta_gamma_lam_mean_simpred_hhssi1
  # Log transform
  log_bg <- log(bg)
  # Rescale to a reasonable point-size range
  cex_vals <- 0.5 +
    2.5 * (log_bg - min(log_bg)) /
    (max(log_bg) - min(log_bg))
  # Prior variances
  x_vals <- log(delta_var_mean_simpred_hhssi1)
  y_vals <- log(delta_var_mean_simpred_nhhssi1)
  idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
  # Colors
  point_col <- ifelse(idx, "blue", "red")
  # Plot
  plot(x_vals, y_vals,
       pch = 21,
       cex = cex_vals,
       bg = point_col,
       col = "black",
       main = "Log-transformed prior variance",
       ylab = "log(NHHSSI prior variance)",
       xlab = "log(HHSSI prior variance)",
       xlim = range(c(x_vals, y_vals)),
       ylim = range(c(x_vals, y_vals))
  )
  # 45-degree reference line
  abline(
    0, 1,
    lty = 2,
    col = "black"
  )
  # Legend for variance comparison
  legend("topleft",
         legend = c("NHHSSI variance > HHSSI variance",
                    "NHHSSI variance <= HHSSI variance"
         ),
         pt.bg = c("blue", "red"),
         col = "black",
         pch = 21,
         bty = "n"
  )





  # Colors
  point_col <- ifelse(idx, "blue", "red")
  # Common limits
  x_vals <- log(delta_tau_mean_simpred_hhssi1)
  y_vals <- log(delta_var_mean_simpred_nhhssi1)
  plot(x_vals, y_vals,
       pch = 19,
       col = point_col,
       main = "HHSSI prior variance structure",
       ylab = expression(NHHSSI_log(hat(lambda)[delta]^2 * hat(tau)[delta]^2)),
       xlab = expression(HHSSI_log(hat(lambda)[delta]^2 * hat(tau)[delta]^2)),
       xlim = range(c(x_vals, y_vals)),
       ylim = range(c(x_vals, y_vals))
  )
  # 45-degree reference line
  abline(0, 1, lty = 2, col = "black")
  # Legend
  legend("topright",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n"
  )


  dat <- round(cbind(log(delta_tau_mean_simpred_hhssi1),
                     log(delta_var_mean_simpred_hhssi1),
                     log(delta_var_mean_simpred_nhhssi1),
                     log(beta_gamma_lam_mean_simpred_hhssi1)), 5)
  dat_sorted <- dat[order(dat[, 3], decreasing = TRUE), ]

  dat_sorted




  ##########################################################

  active_idx <- c(1:3, 6:8)

  J <- 20
  M <- 5

  active_info <- data.frame(
    index = active_idx,
    j = ((active_idx - 1) %/% M) + 1,
    m = ((active_idx - 1) %% M) + 1,
    hhssi_var = delta_var_mean_simpred_hhssi1[active_idx],
    nhhssi_var = delta_var_mean_simpred_nhhssi1[active_idx],
    bg = beta_gamma_lam_mean_simpred_hhssi1[active_idx]
  )
  active_info



  dim(delta_est_normal01_simpred_hhssi1)
  delta_est_normal01_simpred_nhhssi1
  delta_true <- sim_data_normal01_simpred1$delta_true

  rmse_hhssi <- sqrt(rowMeans((delta_est_normal01_simpred_hhssi1 - delta_true)^2))
  rmse_nhhssi <- sqrt(rowMeans((delta_est_normal01_simpred_nhhssi1 - delta_true)^2 ))
  active_rmse <- data.frame(
    index = active_idx,
    j = ((active_idx - 1) %/% M) + 1,
    m = ((active_idx - 1) %% M) + 1,
    HHSSI = rmse_hhssi[active_idx],
    NHHSSI = rmse_nhhssi[active_idx]
  )
  active_rmse




  bias_hhssi <- rowMeans(delta_est_normal01_simpred_hhssi1 - delta_true)
  bias_nhhssi <- rowMeans(delta_est_normal01_simpred_nhhssi1 - delta_true)
  active_bias <- data.frame(
    index = active_idx,
    j = ((active_idx - 1) %/% M) + 1,
    m = ((active_idx - 1) %% M) + 1,
    HHSSI = bias_hhssi[active_idx],
    NHHSSI = bias_nhhssi[active_idx]
  )
  active_bias




  active_absbias <- data.frame(
    index = active_idx,
    j = ((active_idx - 1) %/% M) + 1,
    m = ((active_idx - 1) %% M) + 1,
    HHSSI = abs(bias_hhssi[active_idx]),
    NHHSSI = abs(bias_nhhssi[active_idx])
  )
  active_absbias


  active_results <- data.frame(
    index = active_idx,
    j = ((active_idx - 1) %/% M) + 1,
    m = ((active_idx - 1) %% M) + 1,

    log_bg = log(beta_gamma_lam_mean_simpred_hhssi1[active_idx]),

    log_var_hhssi =
      log(delta_var_mean_simpred_hhssi1[active_idx]),

    log_var_nhhssi =
      log(delta_var_mean_simpred_nhhssi1[active_idx]),

    RMSE_HHSSI = rmse_hhssi[active_idx],
    RMSE_NHHSSI = rmse_nhhssi[active_idx],

    Bias_HHSSI = bias_hhssi[active_idx],
    Bias_NHHSSI = bias_nhhssi[active_idx],

    AbsBias_HHSSI = abs(bias_hhssi[active_idx]),
    AbsBias_NHHSSI = abs(bias_nhhssi[active_idx])
  )

  round(active_results, 4)




  plot(active_rmse$HHSSI, active_rmse$NHHSSI,
       pch = 19,
       xlab = "HHSSI RMSE",
       ylab = "NHHSSI RMSE",
       main = "RMSE for active interaction coefficients",
       xlim = range(c(active_rmse$HHSSI, active_rmse$NHHSSI)),
       ylim = range(c(active_rmse$HHSSI, active_rmse$NHHSSI))
  )
  abline(0, 1, lty = 2)
  text(active_rmse$HHSSI, active_rmse$NHHSSI,
       labels = paste0("(", active_rmse$j, ",", active_rmse$m, ")"),
       pos = 4,
       cex = 0.8)


  plot(active_bias$HHSSI, active_bias$NHHSSI,
       pch = 19,
       xlab = "HHSSI bias",
       ylab = "NHHSSI bias",
       main = "Bias for active interaction coefficients")
  abline(h = 0, v = 0, lty = 2)
  abline(0, 1, lty = 2)
  text(active_bias$HHSSI, active_bias$NHHSSI,
       labels = paste0("(", active_bias$j, ",", active_bias$m, ")"),
       pos = 4,
       cex = 0.8)




  point_col <- ifelse(idx, "blue", "red")
  plot(log(delta_var_mean_simpred_hhssi1), log(rmse_hhssi),
       pch = 19,
       cex = 0.7,
       col = point_col,
       xlab = "log(HHSSI prior variance)",
       ylab = "log(HHSSI RMSE)",
       main = "HHSSI prior variance vs RMSE")
  points(log(delta_var_mean_simpred_hhssi1[active_idx]),
         log(rmse_hhssi[active_idx]),
         pch = 19,
         cex = 1.5,
         col="red")
  text(log(delta_var_mean_simpred_hhssi1[active_idx]),
       log(rmse_hhssi[active_idx]),
       labels = active_idx,
       pos = 4,
       cex = 0.8)
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0,1)




  plot(log(delta_var_mean_simpred_nhhssi1),
       log(rmse_nhhssi),
       pch = 19,
       cex = 0.7,
       col = point_col,
       xlab = "log(NHHSSI prior variance)",
       ylab = "log(NHHSSI RMSE)",
       main = "NHHSSI prior variance vs RMSE")
  points(log(delta_var_mean_simpred_nhhssi1[active_idx]),
         log(rmse_nhhssi[active_idx]),
         pch = 19,
         cex = 1.5,
         col="red")
  text(log(delta_var_mean_simpred_nhhssi1[active_idx]),
       log(rmse_nhhssi[active_idx]),
       labels = active_idx,
       pos = 4,
       cex = 0.8)
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0,1)




  plot(log(delta_var_mean_simpred_hhssi1),
       bias_hhssi,
       pch = 19,
       cex = 0.7,
       col = point_col,
       xlab = "log(HHSSI prior variance)",
       ylab = "HHSSI bias",
       main = "HHSSI prior variance vs bias")
  points(log(delta_var_mean_simpred_hhssi1[active_idx]),
         bias_hhssi[active_idx],
         pch = 19,
         cex = 1.5,
         col="red")
  text(log(delta_var_mean_simpred_hhssi1[active_idx]),
       bias_hhssi[active_idx],
       labels = active_idx,
       pos = 4,
       cex = 0.8)
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(h=0)





  plot(log(delta_var_mean_simpred_nhhssi1),
       bias_nhhssi,
       pch = 19,
       cex = 0.7,
       col = point_col,
       xlab = "log(NHHSSI prior variance)",
       ylab = "NHHSSI bias",
       main = "NHHSSI prior variance vs bias")
  points(log(delta_var_mean_simpred_nhhssi1[active_idx]),
         bias_nhhssi[active_idx],
         pch = 19,
         cex = 1.5,
         col="red")
  text(log(delta_var_mean_simpred_nhhssi1[active_idx]),
       bias_nhhssi[active_idx],
       labels = active_idx,
       pos = 4,
       cex = 0.8)
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(h=0)




  point_col <- ifelse(idx, "blue", "red")
  plot(log(rmse_hhssi), log(rmse_nhhssi),
       pch = 19,
       cex = 0.7,
       col = point_col,
       xlab = "log(HHSSI RMSE)",
       ylab = "log(NHHSSI RMSE)",
       main = "RMSE of HHSSI vs NHHSSI")
  points(log(rmse_hhssi)[active_idx], log(rmse_nhhssi)[active_idx],
         pch = 19,
         cex = 1.5,
         col="red")
  text(log(rmse_hhssi[active_idx]),
       log(rmse_nhhssi[active_idx]),
       labels = active_idx,
       pos = 4,
       cex = 0.8)
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n"
  )
  abline(0,1)





  point_col <- ifelse(idx, "blue", "red")
  plot(bias_hhssi, bias_nhhssi,
       pch = 19,
       cex = 0.7,
       col = point_col,
       xlab = "HHSSI bias",
       ylab = "NHHSSI bias",
       main = "Bias of HHSSI vs NHHSSI")
  points(bias_hhssi[active_idx], bias_nhhssi[active_idx],
         pch = 19,
         cex = 1.5,
         col="red")
  text(bias_hhssi[active_idx],
       bias_nhhssi[active_idx],
       labels = active_idx,
       pos = 4,
       cex = 0.8)
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n"
  )
  abline(h=0)



  ##########################################################################

  delta_est_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 27:126, ], c(2,3), mean)
  delta_est_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 27:126, ], c(2,3), mean)

  delta_true <- sim_data_normal01_simpred1$delta_true

  rmse_hhssi1 <- sqrt(rowMeans((delta_est_normal01_simpred_hhssi1 - delta_true)^2))
  rmse_nhhssi1 <- sqrt(rowMeans((delta_est_normal01_simpred_nhhssi1 - delta_true)^2 ))

  library(ggplot2)
  df <- expand.grid(exposure = 1:M, modifier = 1:J)
  df$RMSE_HHSSI <- rmse_hhssi1
  df$RMSE_NHHSSI <- rmse_nhhssi1
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

    scale_fill_gradient(
      low = "white",
      high = "red",
      name = "RMSE") +

    labs(
      title = "RMSE of Modifier–Exposure Interaction Estimates (with smaller number of predictors)",
      subtitle = "(Black borders indicate active interactions)",
      x = "Exposure index",
      y = "Modifier index") +

    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5))





#####################################################################
# Heatmap of 10 datasets
#####################################################################
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
      title = "RMSE Heatmap of Modifier–Exposure Interaction Estimates (with 10 datasets of simulated predictors)",
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
      title = "Log of RMSE Heatmap of Modifier–Exposure Interaction Estimates (with 10 datasets of simulated predictors)",
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
      title = "Coverage Heatmap of Modifier–Exposure Interaction Estimates (with 10 datasets of simulated predictors)",
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
      title = "Variance Heatmap of Modifier–Exposure Interaction Estimates (with 10 datasets of simulated predictors)",
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
      title = "Log of Variance Heatmap of Modifier–Exposure Interaction Estimates (with 10 datasets of simulated predictors)",
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




