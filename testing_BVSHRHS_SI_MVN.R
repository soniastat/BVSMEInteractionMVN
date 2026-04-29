
devtools::load_all()

library(doParallel) # for parallelization
library(foreach) # for parallelization
library(LaplacesDemon) # IW


#######################################################
# Detect and register cluster
n_cores <- 4
cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export needed data and functions
clusterExport(cl, varlist = "fit_BVSHRHS_SI_MVN", envir = environment())


# Load required packages on each worker
clusterEvalQ(cl, {
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
start_time2 <- Sys.time()

theta_update_RHS_list2 <- foreach(
  sample = 1:100,
  .combine = "list",
  .multicombine = TRUE,
  .errorhandling = "stop",
  .packages = c("MASS",
                "mvtnorm",
                "BayesLogit",
                "actuar",
                "Matrix",
                "greybox",
                "LaplacesDemon")) %dopar%
  {
    cat("Running sample", sample, "...\n")
    set.seed(3567 + sample)

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
    nu_0 <- K+2
    Psi_0 <- diag(K)
    Sigma_init <- rinvwishart(nu_0, Psi_0)


    set.seed(3567 + sample)
    res_all_par2 <- tryCatch({
      fit_BVSHRHS_SI_MVN(
        niter = 6000, burn_in = 1000, thin = 5,
        n, K, Y, W, n_all_par, J, M,
        c = 2.5,
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
        omegasq_beta_lambdasq = 0.25, omegasq_beta_tausq = 0.25,
        omegasq_gamma_lambdasq = 0.25, omegasq_gamma_tausq = 0.25,
        omegasq_delta_lambdasq = 0.25, omegasq_delta_tausq = 0.25,
        accept_lambdasq_beta_init = rep(1, times = J),
        accept_tausq_beta_init = 1,
        accept_lambdasq_gamma_init = rep(1, times = M),
        accept_tausq_gamma_init = 1,
        accept_lambdasq_delta_init = rep(1, times = J*M),
        accept_tausq_delta_init = 1
      )

    }, error = function(e) {
      cat("⚠️ Error in sample", sample, ":", e$message, "\n")
      return(NULL)
    })
    if (!is.null(res_all_par2)) {
      return(res_all_par2$theta_update)
    } else {
      return(NA)
    }
  }

end_time2 <- Sys.time()
stopCluster(cl)

cat("✅ All samples completed in", round(difftime(end_time2, start_time2, units = "mins"), 2), "minutes.\n")
cat("Length of theta_update_RHS_list2:", length(theta_update_RHS_list2), "\n")

################################################################################

nsample <- length(sim_data_list1)
K <- dim(sim_data_list1[[1]]$Y)[2]

alpha_cover_RHS_per2 <- c()
beta_cover_RHS_per2 <- c()
gamma_cover_RHS_per2 <- c()
delta_cover_RHS_per2 <- c()
delta_RHS_RMSE2 <- c()
RMSE_all_main_RHS_effects2 <- c()
beta_cover_nonzero_RHS_per2 <- c()
gamma_cover_nonzero_RHS_per2 <- c()
delta_cover_nonzero_RHS_per2 <- c()
alpha_est_RHS_mat2 <- matrix(NA, nrow = nsample, ncol = K)
alpha_true_RHS_mat2 <- matrix(NA, nrow = nsample, ncol = K)
alpha_lower_RHS_mat2 <- matrix(NA, nrow = nsample, ncol = K)
alpha_upper_RHS_mat2 <- matrix(NA, nrow = nsample, ncol = K)
beta_est_RHS_list2 <- list()
beta_true_RHS_list2 <- list()
beta_lower_RHS_list2 <- list()
beta_upper_RHS_list2 <- list()
gamma_est_RHS_list2<- list()
gamma_true_RHS_list2 <- list()
gamma_lower_RHS_list2 <- list()
gamma_upper_RHS_list2 <- list()
delta_est_RHS_list2 <- list()
delta_true_RHS_list2 <- list()
delta_lower_RHS_list2 <- list()
delta_upper_RHS_list2 <- list()

for(sample in 1:nsample)
{
  Y <- sim_data_list1[[sample]]$Y
  x <- sim_data_list1[[sample]]$x
  z <- sim_data_list1[[sample]]$z
  E <- sim_data_list1[[sample]]$E
  W <- cbind(1, x, z, E)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]


  # alpha
  alpha_true <- sim_data_list1[[sample]]$alpha_true
  alpha_true_RHS_mat2[sample, ] <- alpha_true
  alpha_est_RHS_mat2[sample, ] <- apply(theta_update_RHS_list2[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_RHS_list2[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_RHS_mat2[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_RHS_list2[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_RHS_mat2[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_RHS_per2[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list1[[sample]]$beta_true
  beta_true_RHS_list2[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_RHS_list2[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_RHS_list2[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_RHS_list2[[sample]][, j, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    beta_cover[j-1, ] <- (beta_true[j-1, ] >= beta_lower[j-1, ]) &
      (beta_true[j-1, ] <= beta_upper[j-1, ])
    beta_MSE[j-1, ] <- (beta_est[j-1, ] - beta_true[j-1, ])^2
  }
  beta_est_RHS_list2[[sample]] <- beta_est
  beta_lower_RHS_list2[[sample]] <- beta_lower
  beta_upper_RHS_list2[[sample]] <- beta_upper
  beta_cover_RHS_per2[sample] <- (sum(beta_cover)/(dim(beta_true)[1]*dim(beta_true)[2]))*100
  nonzero_idx_beta <- beta_true != 0
  beta_cover_nonzero_RHS_per2[sample] <- (sum(beta_cover[nonzero_idx_beta]) / sum(nonzero_idx_beta)) * 100


  # gamma
  gamma_true <- sim_data_list1[[sample]]$gamma_true
  gamma_true_RHS_list2[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_RHS_list2[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_RHS_list2[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_RHS_list2[[sample]][, m, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    gamma_cover[(m-J-1), ] <- (gamma_true[(m-J-1), ] >= gamma_lower[(m-J-1), ]) &
      (gamma_true[(m-J-1), ] <= gamma_upper[(m-J-1), ])
    gamma_MSE[m-J-1, ] <- (gamma_est[(m-J-1), ] - gamma_true[m-J-1, ])^2
  }
  gamma_est_RHS_list2[[sample]] <- gamma_est
  gamma_lower_RHS_list2[[sample]] <- gamma_lower
  gamma_upper_RHS_list2[[sample]] <- gamma_upper
  gamma_cover_RHS_per2[sample] <- (sum(gamma_cover)/(dim(gamma_true)[1]*dim(gamma_true)[2]))*100
  nonzero_idx_gamma <- gamma_true != 0
  gamma_cover_nonzero_RHS_per2[sample] <- (sum(gamma_cover[nonzero_idx_gamma]) / sum(nonzero_idx_gamma)) * 100

  main_MSE <- c(as.vector(beta_MSE), as.vector(gamma_MSE))
  RMSE_all_main_RHS_effects2[sample] <- sqrt(mean(main_MSE))

  # delta
  delta_true <- sim_data_list1[[sample]]$delta_true
  delta_true_RHS_list2[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_RHS_list2[[sample]][, jm, ], 2, mean)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_RHS_list2[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_RHS_list2[[sample]][, jm, ], 2, quantile, probs = 0.975)
    # coverage: 1 if true value inside [lower, upper], else 0
    delta_cover[(jm-M-J-1), ] <- (delta_true[(jm-M-J-1), ] >= delta_lower[(jm-M-J-1), ]) &
      (delta_true[(jm-M-J-1), ] <= delta_upper[(jm-M-J-1), ])
    delta_MSE[(jm-M-J-1), ] <- (delta_est[(jm-M-J-1), ] - delta_true[jm-M-J-1, ])^2
  }
  delta_est_RHS_list2[[sample]] <- delta_est
  delta_lower_RHS_list2[[sample]] <- delta_lower
  delta_upper_RHS_list2[[sample]] <- delta_upper
  delta_cover_RHS_per2[sample] <- (sum(delta_cover)/(dim(delta_true)[1]*dim(delta_true)[2]))*100
  nonzero_idx_delta <- delta_true != 0
  delta_cover_nonzero_RHS_per2[sample] <- (sum(delta_cover[nonzero_idx_delta]) / sum(nonzero_idx_delta)) * 100
  delta_RHS_RMSE2[sample] <- sqrt(mean(delta_MSE))


}



(mean(beta_cover_RHS_per2) + mean(gamma_cover_RHS_per2))/2
(mean(beta_cover_nonzero_RHS_per2) + mean(gamma_cover_nonzero_RHS_per2))/2
mean(delta_cover_RHS_per2)
mean(delta_cover_nonzero_RHS_per2)
mean(RMSE_all_main_RHS_effects2)
mean(delta_RHS_RMSE2)
###########################################################################################

sample <- 1

alpha_post_mean <- apply(theta_update_RHS_list2[[sample]][, 1, ], 2, mean)  # 1 × K
alpha_true <- sim_data_list1[[sample]]$alpha_true
alpha_lower <- apply(theta_update_RHS_list2[[sample]][, 1, ], 2, quantile, probs = 0.025)
alpha_upper <- apply(theta_update_RHS_list2[[sample]][, 1, ], 2, quantile, probs = 0.975)

beta_post_mean <- apply(theta_update_RHS_list2[[sample]][, 2:(1 + J), ], c(2, 3), mean)  # J × K
beta_lower <- apply(theta_update_RHS_list2[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.025)  # J × K
beta_upper <- apply(theta_update_RHS_list2[[sample]][, 2:(1 + J), ], c(2, 3), quantile, probs = 0.975)  # J × K
beta_true <- sim_data_list1[[sample]]$beta_true

gamma_post_mean <- apply(theta_update_RHS_list2[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), mean)  # M × K
gamma_lower <- apply(theta_update_RHS_list2[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.025)  # M × K
gamma_upper <- apply(theta_update_RHS_list2[[sample]][, (1 + J + 1):(1 + J + M), ], c(2, 3), quantile, probs = 0.975)  # M × K
gamma_true <- sim_data_list1[[sample]]$gamma_true

delta_post_mean <- apply(theta_update_RHS_list2[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), mean)  # JM × K
delta_lower <- apply(theta_update_RHS_list2[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.025)  # JM × K
delta_upper <- apply(theta_update_RHS_list2[[sample]][, (1 + J + M + 1):(1 + J + M + J*M), ], c(2, 3), quantile, probs = 0.975)  # JM × K
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
       title = "Posterior estimates vs True values of alpha (Model 2 RHS)")

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
       title = "Posterior estimates vs True values of beta (Model 2 RHS)")

#################

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
       title = "Posterior estimates vs True values of gamma (Model 2 RHS)")

#################


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
       title = "Posterior estimates vs True values of delta (Model 2 RHS)")




