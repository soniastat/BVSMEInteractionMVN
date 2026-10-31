rm(list = ls())
gc()


########################


# Run all the functions from the package
devtools::load_all()


load("exposome.RData")
load("modifiers_covariates.rda")
load("exposures.rda")
load("cov_base_post.rda")



raw_Y <- cbind(
  phenotype[,5:6],
  cov_base_post
)

Y0 <- scale(raw_Y)

Z <- cbind(
  Organochlorines_Postnatal,
  Metals_Postnatal
)

Sigma <- cov(Y0)
K <- dim(Y0)[2]
n <- nrow(Y0)
J <- ncol(X)
M <- ncol(Z)
O <- ncol(D)

############################
# Simulate one dataset (first approach to simulate data from weak heredity (100%))
############################

# Simulate data (with weak heredity)
Sim_data_BVS_weak_real1 <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                                   x, z, d, Sigma)
{
  # modifiers
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- rnorm(K)
  beta_true[2, ] <- rnorm(K)

  # environmental exposures
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different gamma for different k
  gamma_true[1, ] <- rnorm(K)
  gamma_true[2, ] <- rnorm(K)
  gamma_true[3, ] <- rnorm(K)

  # modifier-environment interaction terms
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
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- rnorm(K)
  varphi_true[2, ] <- rnorm(K)

  alpha_true <- runif(K, -3.5, -2.5)

  theta_true <- rbind(matrix(alpha_true, nrow = 1),
                      beta_true, gamma_true, delta_true,
                      varphi_true)
  W <- cbind(1, x, z, u, d)
  W_theta_true <- W %*% theta_true # n*K dimension

  # Cholesky of Sigma
  R_Sigma <- chol(Sigma)

  # E ~ iid N(0,1) error
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

set.seed(23444)
sim_data1 <- Sim_data_BVS_weak_real1(K=K, n=n, J=J, M=M,
                                   O=O, x=X, z=Z, d=D,
                                   Sigma=Sigma)
Y1 <- sim_data1$Y
x1 <- sim_data1$x
z1 <- sim_data1$z
u1 <- sim_data1$u
d1 <- sim_data1$d

W1 <- cbind(1, x1, z1, u1, d1)
WTW <- t(W1) %*% W1

delta_true1 <- sim_data1$delta_true
beta_true1 <- sim_data1$beta_true
gamma_true1 <- sim_data1$gamma_true

########################################################################

############################
# Simulate one dataset (second approach to simulate data from weak heredity (100%))
############################

# Simulate data (with weak heredity)
Sim_data_BVS_weak_real2 <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                              x, z, d, Sigma, weak_per = 0.1)
{
  # modifiers
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- rnorm(K)
  beta_true[2, ] <- rnorm(K)
  active_beta  <- which(rowSums(abs(beta_true)) > 0)

  # environmental exposures
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different gamma for different k
  gamma_true[1, ] <- rnorm(K)
  gamma_true[2, ] <- rnorm(K)
  gamma_true[3, ] <- rnorm(K)
  active_gamma <- which(rowSums(abs(gamma_true)) > 0)

  # modifier-environment interaction terms
  u <- matrix(0, n, J * M)
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      u[, col_idx] <- x[, j] * z[, m]
      col_idx <- col_idx + 1
    }
  }

  delta_true <- matrix(0, nrow = J * M, ncol = K)
  # Identify strong-heredity and weak-heredity-only interactions
  strong_idx <- integer(0)
  weak_only_idx <- integer(0)
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      idx <- (j - 1) * M + m
      beta_active <- j %in% active_beta
      gamma_active <- m %in% active_gamma
      if (beta_active && gamma_active)
      {
        # Both corresponding main effects are active
        strong_idx <- c(strong_idx, idx)
      } else if (beta_active || gamma_active) {
        # Exactly one corresponding main effect is active
        weak_only_idx <- c(weak_only_idx, idx)
      }
    }
  }
  # Activate all strong-heredity interactions
  active_strong <- strong_idx
  # Randomly activate some% of weak-heredity-only interactions
  n_active_weak <- round(weak_per * length(weak_only_idx))
  active_weak <- sample(weak_only_idx, n_active_weak)

  # Generate coefficients for active interactions
  active_delta <- c(active_strong, active_weak)
  delta_true[active_delta, ] <- matrix(
    rnorm(length(active_delta) * K),
    nrow = length(active_delta),
    ncol = K
  )

  # covariates
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- rnorm(K)
  varphi_true[2, ] <- rnorm(K)

  alpha_true <- runif(K, -3.5, -2.5)

  theta_true <- rbind(matrix(alpha_true, nrow = 1),
                      beta_true, gamma_true, delta_true,
                      varphi_true)
  W <- cbind(1, x, z, u, d)
  W_theta_true <- W %*% theta_true # n*K dimension

  # Cholesky of Sigma
  R_Sigma <- chol(Sigma)

  # E ~ iid N(0,1) error
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




set.seed(23444)
sim_data2 <- Sim_data_BVS_weak_real2(K=K, n=n, J=J, M=M,
                                    O=O, x=X, z=Z, d=D,
                                    Sigma=Sigma,
                                    weak_per=1)
Y2 <- sim_data2$Y
x2 <- sim_data2$x
z2 <- sim_data2$z
u2 <- sim_data2$u
d2 <- sim_data2$d

W2 <- cbind(1, x2, z2, u2, d2)
WTW <- t(W2) %*% W2

delta_true2 <- sim_data2$delta_true
beta_true2 <- sim_data2$beta_true
gamma_true2 <- sim_data2$gamma_true

########################################################################



############################
# Simulate one dataset (third approach to simulate data from weak heredity (100%))
############################

# Simulate data (with weak heredity)
Sim_data_BVS_weak_real3 <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                               x, z, d, Sigma, weak_per = 0.1,
                               active_weak_idx = NULL)
{
  # modifiers
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- rnorm(K)
  beta_true[2, ] <- rnorm(K)
  active_beta  <- which(rowSums(abs(beta_true)) > 0)

  # environmental exposures
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different gamma for different k
  gamma_true[1, ] <- rnorm(K)
  gamma_true[2, ] <- rnorm(K)
  gamma_true[3, ] <- rnorm(K)
  active_gamma <- which(rowSums(abs(gamma_true)) > 0)

  # modifier-environment interaction terms
  u <- matrix(0, n, J * M)
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      u[, col_idx] <- x[, j] * z[, m]
      col_idx <- col_idx + 1
    }
  }

  # Identify strong and weak-heredity interactions
  strong_idx <- integer(0)
  weak_only_idx <- integer(0)
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      idx <- (j - 1) * M + m
      beta_active <- j %in% active_beta
      gamma_active <- m %in% active_gamma
      if (beta_active && gamma_active)
      {
        # Both main effects are active
        strong_idx <- c(strong_idx, idx)
      } else if (beta_active || gamma_active) {
        # Exactly one main effect is active
        weak_only_idx <- c(weak_only_idx, idx)
      }
    }
  }
  # Select fixed weak-heredity rows
  if (is.null(active_weak_idx))
  {
    n_active_weak <- round(weak_per * length(weak_only_idx))
    active_weak_idx <- sample(weak_only_idx, n_active_weak)
  }
  # Ensure the supplied rows are valid
  stopifnot(all(active_weak_idx %in% weak_only_idx))
  # All strong-heredity rows plus fixed weak-only rows
  active_delta <- c(strong_idx, active_weak_idx)

  # Generate delta coefficients
  delta_true <- matrix(0, nrow = J * M, ncol = K)
  delta_true[active_delta, ] <- matrix(rnorm(length(active_delta) * K),
                                       nrow = length(active_delta),
                                       ncol = K)

  # covariates
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- rnorm(K)
  varphi_true[2, ] <- rnorm(K)

  alpha_true <- runif(K, -3.5, -2.5)

  theta_true <- rbind(matrix(alpha_true, nrow = 1),
                      beta_true, gamma_true, delta_true,
                      varphi_true)
  W <- cbind(1, x, z, u, d)
  W_theta_true <- W %*% theta_true # n*K dimension

  # Cholesky of Sigma
  R_Sigma <- chol(Sigma)

  # E ~ iid N(0,1) error
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
    theta_true = theta_true,
    active_weak_idx = active_weak_idx,
    active_delta = active_delta
  ))
}




grid <- expand.grid(m = 1:M, j = 1:J)
weak_only_idx <- with(grid, (j - 1) * M + m)[xor(grid$m %in% 1:3, grid$j %in% 1:2)]
weak_per <- 1
set.seed(87665467)
fixed_weak_rows <- sample(weak_only_idx, round(weak_per * length(weak_only_idx)))


set.seed(23444)
sim_data3 <- Sim_data_BVS_weak_real3(K=K, n=n, J=J, M=M,
                                     O=O, x=X, z=Z, d=D,
                                     Sigma=Sigma,
                                     weak_per=1,
                                     active_weak_idx=fixed_weak_rows)
Y3 <- sim_data3$Y
x3 <- sim_data3$x
z3 <- sim_data3$z
u3 <- sim_data3$u
d3 <- sim_data3$d

W3 <- cbind(1, x3, z3, u3, d3)
WTW <- t(W3) %*% W3

delta_true3 <- sim_data3$delta_true
beta_true3 <- sim_data3$beta_true
gamma_true3 <- sim_data3$gamma_true

########################################################################



beta_true1==beta_true2
beta_true2==beta_true3


gamma_true1==gamma_true2
gamma_true2==gamma_true3


delta_true1==delta_true2
delta_true1==delta_true3
delta_true2==delta_true3



active1 <- which(rowSums(abs(delta_true1)) > 0)
active2 <- which(rowSums(abs(delta_true2)) > 0)
active3 <- which(rowSums(abs(delta_true3)) > 0)

identical(active1, active2)
identical(active2, active3)







# Extract nonzero interaction coefficients
d1 <- delta_true1[delta_true1 != 0]
d2 <- delta_true2[delta_true2 != 0]
d3 <- delta_true3[delta_true3 != 0]
# Summary statistics
summarize_delta <- function(x) {
  c(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    min = min(x),
    max = max(x)
  )
}
rbind(
  Method1 = summarize_delta(d1),
  Method2 = summarize_delta(d2),
  Method3 = summarize_delta(d3)
)




par(mfrow = c(1, 3))
hist(d1, breaks = 30, probability = TRUE,
     main = "Method 1", xlab = "Nonzero delta")
lines(density(d1), lwd = 2)
hist(d2, breaks = 30, probability = TRUE,
     main = "Method 2", xlab = "Nonzero delta")
lines(density(d2), lwd = 2)
hist(d3, breaks = 30, probability = TRUE,
     main = "Method 3", xlab = "Nonzero delta")
lines(density(d3), lwd = 2)
par(mfrow = c(1, 1))



# Compare coefficients only on common active rows
common_active <- Reduce(intersect, list(active1, active2, active3))
v1 <- as.vector(delta_true1[common_active, ])
v2 <- as.vector(delta_true2[common_active, ])
v3 <- as.vector(delta_true3[common_active, ])
rbind(
  Method1 = summarize_delta(v1),
  Method2 = summarize_delta(v2),
  Method3 = summarize_delta(v3)
)
ks.test(v1, v2)
ks.test(v1, v3)
ks.test(v2, v3)




ks.test(
  as.vector(delta_true1),
  as.vector(delta_true2)
)
ks.test(
  as.vector(delta_true1),
  as.vector(delta_true3)
)
ks.test(
  as.vector(delta_true2),
  as.vector(delta_true3)
)





# Are the active interaction rows identical?
identical(active1, active2)
identical(active1, active3)
identical(active2, active3)
# Are the true coefficient matrices identical?
identical(delta_true1, delta_true2)
identical(delta_true1, delta_true3)
identical(delta_true2, delta_true3)
# Check the expected marginal distribution
# of the nonzero coefficients
par(mfrow = c(1, 1))
qqnorm(d1, main = "Method 1: Normal Q-Q plot")
qqline(d1)
qqnorm(d2, main = "Method 2: Normal Q-Q plot")
qqline(d2)
qqnorm(d3, main = "Method 3: Normal Q-Q plot")
qqline(d3)






# Number of entries that differ
sum(delta_true1 != delta_true2)
sum(delta_true1 != delta_true3)
sum(delta_true2 != delta_true3)

# Compare only the active coefficients
sum(delta_true1[active1, ] != delta_true2[active2, ])
sum(delta_true1[active1, ] != delta_true3[active3, ])
sum(delta_true2[active2, ] != delta_true3[active3, ])

# Maximum absolute difference
max(abs(delta_true1 - delta_true2))
max(abs(delta_true1 - delta_true3))
max(abs(delta_true2 - delta_true3))

# Are the active coefficients equal after sorting?
sort(as.vector(delta_true1[active1, ])) |> head()
sort(as.vector(delta_true2[active2, ])) |> head()
sort(as.vector(delta_true3[active3, ])) |> head()





#######################################################################################
# run for first data
#######################################################################################

#############################################


sim_data_list1 <- list()
for(sample in 1:5)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS_weak_real1(K=K, n=n, J=J, M=M,
                                      O=O, x=X, z=Z, d=D,
                                      Sigma=Sigma)
  sim_data_list1[[sample]] <- sim_data
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

n_cores <- 3

cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export objects/functions needed by workers
clusterExport(cl,
              varlist = c(
                "fit_HHSSI_MVN_cov_modify3",
                "sim_data_list1"
              ),
              envir = environment()
)

start_time_HHSSI <- Sys.time()

theta_update_HHSSI_list1 <- foreach(
  sample = 1:5,
  .packages = c("LaplacesDemon", "Matrix", "actuar"),
  .combine = "list",
  .multicombine = TRUE
) %dorng% {

  cat("Running sample", sample, "\n")

  Y <- sim_data_list1[[sample]]$Y
  x <- sim_data_list1[[sample]]$x
  z <- sim_data_list1[[sample]]$z
  u <- sim_data_list1[[sample]]$u
  d <- sim_data_list1[[sample]]$d

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
      niter = 100,
      burn_in = 10,
      thin = 2,
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


nsample <- length(sim_data_list1)
K <- dim(sim_data_list1[[1]]$Y)[2]
J <- dim(sim_data_list1[[1]]$x)[2]
M <- dim(sim_data_list1[[1]]$z)[2]


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
  Y <- sim_data_list1[[sample]]$Y
  x <- sim_data_list1[[sample]]$x
  z <- sim_data_list1[[sample]]$z
  u <- sim_data_list1[[sample]]$u
  d <- sim_data_list1[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list1[[sample]]$alpha_true
  alpha_true_mat_HHSSI[sample, ] <- alpha_true
  alpha_est_mat_HHSSI[sample, ] <- apply(theta_update_HHSSI_list1[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_HHSSI_list1[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_HHSSI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_HHSSI_list1[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_HHSSI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_HHSSI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list1[[sample]]$beta_true
  beta_true_list_HHSSI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_HHSSI_list1[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_HHSSI_list1[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_HHSSI_list1[[sample]][, j, ], 2, quantile, probs = 0.975)
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
  gamma_true <- sim_data_list1[[sample]]$gamma_true
  gamma_true_list_HHSSI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_HHSSI_list1[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_HHSSI_list1[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_HHSSI_list1[[sample]][, m, ], 2, quantile, probs = 0.975)
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
  varphi_true <- sim_data_list1[[sample]]$varphi_true
  varphi_true_list_HHSSI[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list1[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list1[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list1[[sample]][, o, ], 2, quantile, probs = 0.975)
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
  delta_true <- sim_data_list1[[sample]]$delta_true
  delta_true_list_HHSSI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list1[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list1[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list1[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list1[[sample]][, jm, ], 2, quantile, probs = 0.975)
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

n_cores <- 3

cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export objects/functions needed by workers
clusterExport(cl,
              varlist = c(
                "fit_NHHSSI_MVN_cov_modify3",
                "sim_data_list1"
              ),
              envir = environment()
)

start_time_NHHSSI <- Sys.time()

theta_update_NHHSSI_list1 <- foreach(
  sample = 1:5,
  .packages = c("LaplacesDemon", "Matrix", "actuar"),
  .combine = "list",
  .multicombine = TRUE
) %dorng% {

  cat("Running sample", sample, "\n")

  Y <- sim_data_list1[[sample]]$Y
  x <- sim_data_list1[[sample]]$x
  z <- sim_data_list1[[sample]]$z
  u <- sim_data_list1[[sample]]$u
  d <- sim_data_list1[[sample]]$d

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
      niter = 100, burn_in = 10, thin = 2,
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


nsample <- length(sim_data_list1)
K <- dim(sim_data_list1[[1]]$Y)[2]

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
  Y <- sim_data_list1[[sample]]$Y
  x <- sim_data_list1[[sample]]$x
  z <- sim_data_list1[[sample]]$z
  u <- sim_data_list1[[sample]]$u
  d <- sim_data_list1[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list1[[sample]]$alpha_true
  alpha_true_mat_NHHSSI[sample, ] <- alpha_true
  alpha_est_mat_NHHSSI[sample, ] <- apply(theta_update_NHHSSI_list1[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NHHSSI_list1[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NHHSSI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NHHSSI_list1[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NHHSSI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NHHSSI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list1[[sample]]$beta_true
  beta_true_list_NHHSSI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NHHSSI_list1[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NHHSSI_list1[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NHHSSI_list1[[sample]][, j, ], 2, quantile, probs = 0.975)
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
  gamma_true <- sim_data_list1[[sample]]$gamma_true
  gamma_true_list_NHHSSI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NHHSSI_list1[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NHHSSI_list1[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NHHSSI_list1[[sample]][, m, ], 2, quantile, probs = 0.975)
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
  varphi_true <- sim_data_list1[[sample]]$varphi_true
  varphi_true_list_NHHSSI[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list1[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list1[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list1[[sample]][, o, ], 2, quantile, probs = 0.975)
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
  delta_true <- sim_data_list1[[sample]]$delta_true
  delta_true_list_NHHSSI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list1[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list1[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list1[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list1[[sample]][, jm, ], 2, quantile, probs = 0.975)
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




#####################################################################
# Heatmap of 5 datasets
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




#######################################################################################
# run for third data
#######################################################################################

grid <- expand.grid(m = 1:M, j = 1:J)
weak_only_idx <- with(grid, (j - 1) * M + m)[xor(grid$m %in% 1:3, grid$j %in% 1:2)]
weak_per <- 1
set.seed(87665467)
fixed_weak_rows <- sample(weak_only_idx, round(weak_per * length(weak_only_idx)))


sim_data_list3 <- list()
for(sample in 1:5)
{
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS_weak_real3(K=K, n=n, J=J, M=M,
                                      O=O, x=X, z=Z, d=D,
                                      Sigma=Sigma,
                                      weak_per=1,
                                      active_weak_idx=fixed_weak_rows)
  sim_data_list3[[sample]] <- sim_data
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

n_cores <- 3

cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export objects/functions needed by workers
clusterExport(cl,
              varlist = c(
                "fit_HHSSI_MVN_cov_modify3",
                "sim_data_list3"
              ),
              envir = environment()
)

start_time_HHSSI <- Sys.time()

theta_update_HHSSI_list3 <- foreach(
  sample = 1:5,
  .packages = c("LaplacesDemon", "Matrix", "actuar"),
  .combine = "list",
  .multicombine = TRUE
) %dorng% {

  cat("Running sample", sample, "\n")

  Y <- sim_data_list3[[sample]]$Y
  x <- sim_data_list3[[sample]]$x
  z <- sim_data_list3[[sample]]$z
  u <- sim_data_list3[[sample]]$u
  d <- sim_data_list3[[sample]]$d

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
      niter = 100,
      burn_in = 10,
      thin = 2,
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


nsample <- length(sim_data_list3)
K <- dim(sim_data_list3[[1]]$Y)[2]
J <- dim(sim_data_list3[[1]]$x)[2]
M <- dim(sim_data_list3[[1]]$z)[2]


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
  Y <- sim_data_list3[[sample]]$Y
  x <- sim_data_list3[[sample]]$x
  z <- sim_data_list3[[sample]]$z
  u <- sim_data_list3[[sample]]$u
  d <- sim_data_list3[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list3[[sample]]$alpha_true
  alpha_true_mat_HHSSI[sample, ] <- alpha_true
  alpha_est_mat_HHSSI[sample, ] <- apply(theta_update_HHSSI_list3[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_HHSSI_list3[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_HHSSI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_HHSSI_list3[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_HHSSI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_HHSSI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list3[[sample]]$beta_true
  beta_true_list_HHSSI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_HHSSI_list3[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_HHSSI_list3[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_HHSSI_list3[[sample]][, j, ], 2, quantile, probs = 0.975)
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
  gamma_true <- sim_data_list3[[sample]]$gamma_true
  gamma_true_list_HHSSI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_HHSSI_list3[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_HHSSI_list3[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_HHSSI_list3[[sample]][, m, ], 2, quantile, probs = 0.975)
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
  varphi_true <- sim_data_list3[[sample]]$varphi_true
  varphi_true_list_HHSSI[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list3[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list3[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_HHSSI_list3[[sample]][, o, ], 2, quantile, probs = 0.975)
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
  delta_true <- sim_data_list3[[sample]]$delta_true
  delta_true_list_HHSSI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list3[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list3[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list3[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_HHSSI_list3[[sample]][, jm, ], 2, quantile, probs = 0.975)
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

n_cores <- 3

cat("Using", n_cores, "cores for parallel computation.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export objects/functions needed by workers
clusterExport(cl,
              varlist = c(
                "fit_NHHSSI_MVN_cov_modify3",
                "sim_data_list3"
              ),
              envir = environment()
)

start_time_NHHSSI <- Sys.time()

theta_update_NHHSSI_list3 <- foreach(
  sample = 1:5,
  .packages = c("LaplacesDemon", "Matrix", "actuar"),
  .combine = "list",
  .multicombine = TRUE
) %dorng% {

  cat("Running sample", sample, "\n")

  Y <- sim_data_list3[[sample]]$Y
  x <- sim_data_list3[[sample]]$x
  z <- sim_data_list3[[sample]]$z
  u <- sim_data_list3[[sample]]$u
  d <- sim_data_list3[[sample]]$d

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
      niter = 100, burn_in = 10, thin = 2,
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


nsample <- length(sim_data_list3)
K <- dim(sim_data_list3[[1]]$Y)[2]

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
  Y <- sim_data_list3[[sample]]$Y
  x <- sim_data_list3[[sample]]$x
  z <- sim_data_list3[[sample]]$z
  u <- sim_data_list3[[sample]]$u
  d <- sim_data_list3[[sample]]$d
  W <- cbind(1, x, z, u, d)
  K <- dim(Y)[2] # Number of groups
  n <- dim(Y)[1]
  J <- dim(x)[2]
  M <- dim(z)[2]
  O <- dim(d)[2]


  # alpha
  alpha_true <- sim_data_list3[[sample]]$alpha_true
  alpha_true_mat_NHHSSI[sample, ] <- alpha_true
  alpha_est_mat_NHHSSI[sample, ] <- apply(theta_update_NHHSSI_list3[[sample]][, 1, ], 2, mean)
  alpha_lower <- apply(theta_update_NHHSSI_list3[[sample]][, 1, ], 2, quantile, probs = 0.025)
  alpha_lower_mat_NHHSSI[sample, ] <- alpha_lower
  alpha_upper <- apply(theta_update_NHHSSI_list3[[sample]][, 1, ], 2, quantile, probs = 0.975)
  alpha_upper_mat_NHHSSI[sample, ] <- alpha_upper
  # coverage: 1 if true value inside [lower, upper], else 0
  alpha_cover <- (alpha_true >= alpha_lower) &
    (alpha_true <= alpha_upper)
  alpha_cover_per_NHHSSI[sample] <- (sum(alpha_cover)/length(alpha_true))*100

  # beta
  beta_true <- sim_data_list3[[sample]]$beta_true
  beta_true_list_NHHSSI[[sample]] <- beta_true
  beta_est <- matrix(NA, nrow = J, ncol = K)
  beta_lower <- matrix(NA, nrow = J, ncol = K)
  beta_upper <- matrix(NA, nrow = J, ncol = K)
  beta_cover <- matrix(NA, nrow = J, ncol = K)  # coverage indicator
  beta_MSE <- matrix(NA, nrow = J, ncol = K)
  for (j in 2:(1 + J))
  {
    beta_est[j-1, ] <- apply(theta_update_NHHSSI_list3[[sample]][, j, ], 2, mean)
    beta_lower[j-1, ] <- apply(theta_update_NHHSSI_list3[[sample]][, j, ], 2, quantile, probs = 0.025)
    beta_upper[j-1, ] <- apply(theta_update_NHHSSI_list3[[sample]][, j, ], 2, quantile, probs = 0.975)
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
  gamma_true <- sim_data_list3[[sample]]$gamma_true
  gamma_true_list_NHHSSI[[sample]] <- gamma_true
  gamma_est <- matrix(NA, nrow = M, ncol = K)
  gamma_lower <- matrix(NA, nrow = M, ncol = K)
  gamma_upper <- matrix(NA, nrow = M, ncol = K)
  gamma_cover <- matrix(NA, nrow = M, ncol = K)  # coverage indicator
  gamma_MSE <- matrix(NA, nrow = M, ncol = K)
  for (m in (1 + J + 1):(1 + J + M))
  {
    gamma_est[(m-J-1), ] <- apply(theta_update_NHHSSI_list3[[sample]][, m, ], 2, mean)
    gamma_lower[(m-J-1), ] <- apply(theta_update_NHHSSI_list3[[sample]][, m, ], 2, quantile, probs = 0.025)
    gamma_upper[(m-J-1), ] <- apply(theta_update_NHHSSI_list3[[sample]][, m, ], 2, quantile, probs = 0.975)
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
  varphi_true <- sim_data_list3[[sample]]$varphi_true
  varphi_true_list_NHHSSI[[sample]] <- varphi_true
  varphi_est <- matrix(NA, nrow = O, ncol = K)
  varphi_lower <- matrix(NA, nrow = O, ncol = K)
  varphi_upper <- matrix(NA, nrow = O, ncol = K)
  varphi_cover <- matrix(NA, nrow = O, ncol = K)  # coverage indicator
  varphi_MSE <- matrix(NA, nrow = O, ncol = K)
  for (o in (1 + J + M + J*M + 1):(1 + J + M + J*M + O))
  {
    varphi_est[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list3[[sample]][, o, ], 2, mean)
    varphi_lower[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list3[[sample]][, o, ], 2, quantile, probs = 0.025)
    varphi_upper[(o-1-J-M-J*M), ] <- apply(theta_update_NHHSSI_list3[[sample]][, o, ], 2, quantile, probs = 0.975)
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
  delta_true <- sim_data_list3[[sample]]$delta_true
  delta_true_list_NHHSSI[[sample]] <- delta_true
  delta_est <- matrix(NA, nrow = J*M, ncol = K)
  delta_var <- matrix(NA, nrow = J*M, ncol = K)
  delta_lower <- matrix(NA, nrow = J*M, ncol = K)
  delta_upper <- matrix(NA, nrow = J*M, ncol = K)
  delta_cover <- matrix(NA, nrow = J*M, ncol = K)  # coverage indicator
  delta_MSE <- matrix(NA, nrow = J*M, ncol = K)
  for (jm in (1 + J + M + 1):(1 + J + M + J*M))
  {
    delta_est[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list3[[sample]][, jm, ], 2, mean)
    delta_var[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list3[[sample]][, jm, ], 2, var)
    delta_lower[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list3[[sample]][, jm, ], 2, quantile, probs = 0.025)
    delta_upper[(jm-M-J-1), ] <- apply(theta_update_NHHSSI_list3[[sample]][, jm, ], 2, quantile, probs = 0.975)
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




#####################################################################
# Heatmap of 5 datasets
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














