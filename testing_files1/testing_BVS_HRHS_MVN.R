

# Testing codes for Bayesian variable selection with hierarchical regularized
# Horseshoe priors for multivariate normal response (BVS-HRHS-MVN)

rm(list = ls())
gc()

###############################

# Run all the functions from the package
devtools::load_all()

library(doParallel) # for parallelization
library(foreach) # for parallelization
library(LaplacesDemon) # IW

library(doRNG) # used for instead of set.seed

library(ggplot2) # plot

#######################################################

# Simulate data
# J = 10, M = 5 (for example)
nsample <- 100
sim_data_list1 <- list()
for(sample in 1:nsample)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 10, M = 5)
  sim_data_list1[[sample]] <- sim_data
}

#######################################################
# Fit the BVS-HRHS-MVN model

# Detect and register cluster
n_cores <- 4
cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# For reproducibility
registerDoRNG(3567)

# Export needed data and functions
clusterExport(cl, varlist = "fit_BVS_HRHS_MVN", envir = environment())


# Load required packages on each worker
clusterEvalQ(cl, {
  devtools::load_all()

  library(MASS)
  library(mvtnorm)
  library(BayesLogit)
  library(actuar)
  library(Matrix)
  library(greybox)
  library(LaplacesDemon)
  NULL
})

# Run in parallel
start_time_BVS_HRHS_MVN <- Sys.time()

theta_update_list_BVS_HRHS_MVN <- foreach(
  sample = 1:100,
  .combine = "list",
  .multicombine = TRUE,
  .errorhandling = "pass",
  .packages = c("MASS",
                "mvtnorm",
                "BayesLogit",
                "actuar",
                "Matrix",
                "greybox",
                "LaplacesDemon")) %dopar%
  {
    cat("Running sample", sample, "...\n")

    Y <- sim_data_list1[[sample]]$Y
    x <- sim_data_list1[[sample]]$x
    z <- sim_data_list1[[sample]]$z
    u <- sim_data_list1[[sample]]$u
    W <- cbind(1, x, z, u)
    K <- ncol(Y)
    n <- dim(Y)[1]
    J <- dim(x)[2]
    M <- dim(z)[2]
    n_all_par <- dim(W)[2]

    Psi_0 <- diag(K)          # K x K
    nu_0 <- K + 2 # nu_0 should be greater that K+1
    Sigma_init <- rinvwishart(nu_0, Psi_0)

    res_all_par_BVS_HRHS_MVN <- tryCatch({
      fit_BVS_HRHS_MVN(
        niter = 6000, burn_in = 1000, thin = 5,
        n, K, Y, W, n_all_par, J, M,
        c = 2.5,
        theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
        lambdasq_beta_init = matrix(0.5, nrow = J, ncol = K),
        tausq_beta_init = rep(1, J),
        lambdasq_gamma_init = matrix(0.5, nrow = M, ncol = K),
        tausq_gamma_init = rep(1, M),
        lambdasq_delta_init = matrix(0.5, nrow = J*M, ncol = K),
        tausq_delta_init = rep(1, J*M),
        psi_beta_init = matrix(0.5, nrow = J, ncol = K),
        psi_gamma_init = matrix(0.5, nrow = M, ncol = K),
        psi_delta_init = matrix(0.5, nrow = J*M, ncol = K),
        xi_beta_init = 1, xi_gamma_init = 1,
        xi_delta_init = 1,
        Sigma_init = Sigma_init,
        nu_0 = nu_0, Psi_0 = Psi_0,
        omegasq_beta_lambdasq = 0.25, omegasq_beta_tausq = 0.25,
        omegasq_gamma_lambdasq = 0.25, omegasq_gamma_tausq = 0.25,
        omegasq_delta_lambdasq = 0.25, omegasq_delta_tausq = 0.25,
        accept_lambdasq_beta_init = matrix(1, nrow = J, ncol = K),
        accept_tausq_beta_init = rep(1, times = J),
        accept_lambdasq_gamma_init = matrix(1, nrow = M, ncol = K),
        accept_tausq_gamma_init = rep(1, times = M),
        accept_lambdasq_delta_init = matrix(1, nrow = J*M, ncol = K),
        accept_tausq_delta_init = rep(1, times = M)
      )

    }, error = function(e) {
      cat("⚠️ Error in sample", sample, ":", e$message, "\n")
      return(NULL)
    })
    if (!is.null(res_all_par_BVS_HRHS_MVN)) {
      return(res_all_par_BVS_HRHS_MVN$theta_update)
    } else {
      return(NA)
    }
  }

end_time_BVS_HRHS_MVN <- Sys.time()
stopCluster(cl)

cat("✅ All samples completed in", round(difftime(end_time_BVS_HRHS_MVN, start_time_BVS_HRHS_MVN, units = "mins"), 2), "minutes.\n")
cat("Length of theta_update_list_BVS_HRHS_MVN:", length(theta_update_list_BVS_HRHS_MVN), "\n")

################################################################################

K <- dim(sim_data_list1[[1]]$Y)[2]

alpha_cover_per_BVS_HRHS_MVN <- c()
beta_cover_per_BVS_HRHS_MVN <- c()
gamma_cover_per_BVS_HRHS_MVN <- c()
delta_cover_per_BVS_HRHS_MVN <- c()
delta_RMSE_BVS_HRHS_MVN <- c()
RMSE_all_main_effects_BVS_HRHS_MVN <- c()
beta_cover_nonzero_per_BVS_HRHS_MVN <- c()
gamma_cover_nonzero_per_BVS_HRHS_MVN <- c()
delta_cover_nonzero_per_BVS_HRHS_MVN <- c()
alpha_est_mat_BVS_HRHS_MVN <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_mat_BVS_HRHS_MVN <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_mat_BVS_HRHS_MVN <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_mat_BVS_HRHS_MVN <- matrix(NA, nrow = nsample, ncol = K)
beta_est_list_BVS_HRHS_MVN <- list()
beta_true_list_BVS_HRHS_MVN <- list()
beta_lower_list_BVS_HRHS_MVN <- list()
beta_upper_list_BVS_HRHS_MVN <- list()
gamma_est_list_BVS_HRHS_MVN <- list()
gamma_true_list_BVS_HRHS_MVN <- list()
gamma_lower_list_BVS_HRHS_MVN <- list()
gamma_upper_list_BVS_HRHS_MVN <- list()
delta_est_list_BVS_HRHS_MVN <- list()
delta_true_list_BVS_HRHS_MVN <- list()
delta_lower_list_BVS_HRHS_MVN <- list()
delta_upper_list_BVS_HRHS_MVN <- list()

for(sample in 1:nsample)
{
  Y <- sim_data_list1[[sample]]$Y
  x <- sim_data_list1[[sample]]$x
  z <- sim_data_list1[[sample]]$z
  u <- sim_data_list1[[sample]]$E
  W <- cbind(1, x, z, u)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list1[[sample]]$alpha_true
  alpha_true_mat_BVS_HRHS_MVN[sample, ] <- alpha_true
  alpha_est_mat_BVS_HRHS_MVN[sample, ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_BVS_HRHS_MVN[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_BVS_HRHS_MVN[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_BVS_HRHS_MVN[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list1[[sample]]$beta_true
  beta_true_list_BVS_HRHS_MVN[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_list_BVS_HRHS_MVN[[sample]] <- beta_est
  beta_lower_list_BVS_HRHS_MVN[[sample]] <- beta_lower
  beta_upper_list_BVS_HRHS_MVN[[sample]] <- beta_upper
  beta_cover_per_BVS_HRHS_MVN[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_per_BVS_HRHS_MVN[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list1[[sample]]$gamma_true
  gamma_true_list_BVS_HRHS_MVN[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_list_BVS_HRHS_MVN[[sample]] <- gamma_est
  gamma_lower_list_BVS_HRHS_MVN[[sample]] <- gamma_lower
  gamma_upper_list_BVS_HRHS_MVN[[sample]] <- gamma_upper
  gamma_cover_per_BVS_HRHS_MVN[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_per_BVS_HRHS_MVN[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(beta_MSE, gamma_MSE)
  RMSE_all_main_effects_BVS_HRHS_MVN[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list1[[sample]]$delta_true
  delta_true_list_BVS_HRHS_MVN[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_list_BVS_HRHS_MVN[[sample]] <- delta_est
  delta_lower_list_BVS_HRHS_MVN[[sample]] <- delta_lower
  delta_upper_list_BVS_HRHS_MVN[[sample]] <- delta_upper
  delta_cover_per_BVS_HRHS_MVN[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_per_BVS_HRHS_MVN[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RMSE_BVS_HRHS_MVN[sample] <- sqrt(mean(delta_MSE))


}




(mean(beta_cover_per_BVS_HRHS_MVN ) + mean(gamma_cover_per_BVS_HRHS_MVN ))/2
(mean(beta_cover_nonzero_per_BVS_HRHS_MVN ) + mean(gamma_cover_nonzero_per_BVS_HRHS_MVN ))/2
mean(delta_cover_per_BVS_HRHS_MVN)
mean(delta_cover_nonzero_per_BVS_HRHS_MVN )
mean(RMSE_all_main_effects_BVS_HRHS_MVN)
mean(delta_RMSE_BVS_HRHS_MVN)

###########################################################################################

sample <- 1

alpha_post_mean <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 1, ], 2, mean)  # 1 × K
alpha_true <- sim_data_list1[[sample]]$alpha_true
alpha_lower <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 1, ], 2, quantile, probs = 0.025)
alpha_upper <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 1, ], 2, quantile, probs = 0.975)

beta_post_mean <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 2:(1 + J), ], c(2, 3), mean)  # J × K
beta_lower <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.025)  # J × K
beta_upper <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.975)  # J × K
beta_true <- sim_data_list1[[sample]]$beta_true

gamma_post_mean <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), mean)  # M × K
gamma_lower <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.025)  # M × K
gamma_upper <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.975)  # M × K
gamma_true <- sim_data_list1[[sample]]$gamma_true

delta_post_mean <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), mean)  # JM × K
delta_lower <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.025)  # JM × K
delta_upper <- apply(theta_update_list_BVS_HRHS_MVN[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.975)  # JM × K
delta_true <- sim_data_list1[[sample]]$delta_true

###############

alpha_df <- data.frame(
  coef = paste0("A", "_k", rep(1:K)),
  mean = as.vector(alpha_post_mean),
  lower = as.vector(alpha_lower),
  upper = as.vector(alpha_upper),
  true  = as.vector(alpha_true)
)

ggplot(alpha_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(y = "Alpha coefficient", x = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of alpha (BVS-HRHS-MVN model)")

################

# put into long form
beta_df <- data.frame(
  coef = paste0("B", rep(1:J, times = K), "_k", rep(1:K, each = J)),
  mean = as.vector(beta_post_mean),
  lower = as.vector(beta_lower),
  upper = as.vector(beta_upper),
  true  = as.vector(beta_true)
)

ggplot(beta_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(y = "Beta coefficient", x = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of beta (BVS-HRHS-MVN model)")


###############
# put into long form
gamma_df <- data.frame(
  coef = paste0("g", rep(1:M, times = K), "_k", rep(1:K, each = M)),
  mean = as.vector(gamma_post_mean),
  lower = as.vector(gamma_lower),
  upper = as.vector(gamma_upper),
  true  = as.vector(gamma_true)
)

ggplot(gamma_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(x = "gamma coefficient", y = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of gamma (BVS-HRHS-MVN model)")


###############
# put into long form
delta_df <- data.frame(
  coef = paste0("d", rep(1:(J*M), times = K), "_k", rep(1:K, each = (J*M))),
  mean = as.vector(delta_post_mean),
  lower = as.vector(delta_lower),
  upper = as.vector(delta_upper),
  true  = as.vector(delta_true)
)

ggplot(delta_df, aes(x = mean, y = coef)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), shape = 16, size = 0.3) +
  geom_point(aes(x = true), shape = 15, color = "red", size = 1.25) +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = rel(0.9))) +
  labs(x = "delta coefficient", y = "Posterior mean with 95% CrI",
       title = "Posterior estimates vs True values of delta (BVS-HRHS-MVN model)")



