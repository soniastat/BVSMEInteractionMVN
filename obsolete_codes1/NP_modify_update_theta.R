
# Functions for simulating data and fitting BVS-NP-MVN-cov model


# Simulate data
set.seed(23444)
sim_data <- Sim_data_BVS(K = 4, n = 100, J = 44, M = 17)


Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
n <- dim(Y)[1]
O <- 6
d <- matrix(runif(n * O), nrow = n)
W <- cbind(1, x, z, u, d)
K <- ncol(Y)
J <- dim(x)[2]
M <- dim(z)[2]
n_all_par <- dim(W)[2]
WTW <- t(W) %*% W


library(Matrix)
library(MASS)
library(LaplacesDemon)


nu_0 <- K+2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


#########################################################
# Fit NP model
########################################################

# Update theta
update_theta_NP_MVN_cov_block_rnorm <- function(Y, K, W, WTW, n_all_par,
                       J, M, O,
                       sigmasq_varphi, Sigma_update)
{
  # Build V_theta (same for all k)
  var_alpha_0 <- 100 # dim 1*1
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  # vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
  vec_sigmasq_varphi <- rep(10, times = O)
  mat_var_varphi <- diag(vec_sigmasq_varphi)

  # Inverse using blocks separately
  V_theta_inv_block <- as.matrix(
    bdiag(
      solve(var_alpha_0),
      solve(var_beta),
      solve(var_gamma),
      solve(var_delta),
      solve(mat_var_varphi)
    )
  )

  # prior precision
  prior_prec <- kronecker(Diagonal(K), V_theta_inv_block)  # I_K ⊗ Vθ⁻¹

  # Kronecker term from likelihood
  Sigma_inv <- as.matrix(chol2inv(chol(Sigma_update)))
  # Sigma_inv <- as.matrix(chol2inv(chol(Sigma_init)))
  kron_term <- kronecker(Sigma_inv, WTW)  # Σ⁻¹ ⊗ WᵀW

  kr_pr <- kron_term + prior_prec
  # Full conditional covariance
  Sigma_theta <- chol2inv(chol(kr_pr))
  R <- chol(Sigma_theta) # right triangular R^TR = Sigma_theta

  # Full conditional mean
  vecY <- as.vector(Y)
  mean_theta_vec <- Sigma_theta %*% (kronecker(Sigma_inv, t(W)) %*% vecY)


  z_norm <- rnorm(length(mean_theta_vec))
  theta_vec <- as.vector(mean_theta_vec + t(R) %*% z_norm)
  # Reshape back to n_all_par x K
  theta_update_s <- matrix(theta_vec, nrow = n_all_par, ncol = K)



  return(theta_update_s)
}


#######################################################################################

# Update Sigma
update_Sigma_NP_MVN_cov <- function(Y, n, W, theta_update,
                                    Psi_0, nu_0)
{
  # Residuals
  Resid <- Y - W %*% theta_update

  # Posterior scale and df
  Psi_n <- Psi_0 + t(Resid) %*% Resid
  nu_n <- nu_0 + n

  # Draw from IW
  Sigma_update_s <- rinvwishart(nu_n, Psi_n)

  return(Sigma_update_s)
}


#######################################################################################

# Update parameters for NP-MVN-cov model
fit_NP_MVN_cov_block_rnorm <- function(niter = 20, burn_in = 2, thin = 1,
                           n, K, Y, W,
                           WTW, n_all_par,
                           J, M, O,
                           theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                           Sigma_init = diag(K),
                           nu_0, Psi_0,
                           sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))

  # Initialize
  theta_update[1, , ] <- theta_init
  Sigma_update[1, , ] <- Sigma_init

  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_NP_MVN_cov_block_rnorm(Y, K, W, WTW, n_all_par,
                                              J, M, O,
                                              sigmasq_varphi = sigmasq_varphi,
                                              Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s


    # Sigma_update
    Sigma_update_s <- update_Sigma_NP_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ]
  )
  return(results)
}


#######################################################

start_time <- Sys.time()
res_block_rnorm <- fit_NP_MVN_cov_block_rnorm(niter = 10, burn_in = 2, thin = 1,
                            n = n, K = K, Y = Y, W = W,
                            WTW = WTW, n_all_par = n_all_par,
                            J = J, M = M, O = O,
                            theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                            Sigma_init = Sigma_init,
                            nu_0 = nu_0, Psi_0 = Psi_0,
                            sigmasq_varphi = 10)
end_time <- Sys.time()
time_block_rnorm <- end_time - start_time #6.842803 secs
#################################################################

# Update theta
update_theta_NP_MVN_cov_rnorm <- function(Y, K, W, WTW, n_all_par,
                                          J, M, O,
                                          sigmasq_varphi, Sigma_update)
{
  # Build V_theta (same for all k)
  var_alpha_0 <- 100 # dim 1*1
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  # vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
  vec_sigmasq_varphi <- rep(10, times = O)
  mat_var_varphi <- diag(vec_sigmasq_varphi)

  V_theta <- as.matrix(bdiag(var_alpha_0, var_beta, var_gamma,
                             var_delta, mat_var_varphi))

  # prior precision
  V_theta_inv <- chol2inv(chol(V_theta))
  prior_prec <- kronecker(Diagonal(K), V_theta_inv)  # I_K ⊗ Vθ⁻¹

  # Kronecker term from likelihood
  Sigma_inv <- as.matrix(chol2inv(chol(Sigma_update)))
  # Sigma_inv <- as.matrix(chol2inv(chol(Sigma_init)))
  kron_term <- kronecker(Sigma_inv, WTW)  # Σ⁻¹ ⊗ WᵀW

  kr_pr <- kron_term + prior_prec
  # Full conditional covariance
  Sigma_theta <- chol2inv(chol(kr_pr))
  R <- chol(Sigma_theta) # right triangular R^TR = Sigma_theta

  # Full conditional mean
  vecY <- as.vector(Y)
  mean_theta_vec <- Sigma_theta %*% (kronecker(Sigma_inv, t(W)) %*% vecY)


  z_norm <- rnorm(length(mean_theta_vec))
  theta_vec <- as.vector(mean_theta_vec + t(R) %*% z_norm)
  # Reshape back to n_all_par x K
  theta_update_s <- matrix(theta_vec, nrow = n_all_par, ncol = K)



  return(theta_update_s)
}


#######################################################################################

# Update Sigma
update_Sigma_NP_MVN_cov <- function(Y, n, W, theta_update,
                                    Psi_0, nu_0)
{
  # Residuals
  Resid <- Y - W %*% theta_update

  # Posterior scale and df
  Psi_n <- Psi_0 + t(Resid) %*% Resid
  nu_n <- nu_0 + n

  # Draw from IW
  Sigma_update_s <- rinvwishart(nu_n, Psi_n)

  return(Sigma_update_s)
}


#######################################################################################

# Update parameters for NP-MVN-cov model
fit_NP_MVN_cov_rnorm <- function(niter = 20, burn_in = 2, thin = 1,
                                 n, K, Y, W,
                                 WTW, n_all_par,
                                 J, M, O,
                                 theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                                 Sigma_init = diag(K),
                                 nu_0, Psi_0,
                                 sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))

  # Initialize
  theta_update[1, , ] <- theta_init
  Sigma_update[1, , ] <- Sigma_init

  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_NP_MVN_cov_rnorm(Y, K, W, WTW, n_all_par,
                                                    J, M, O,
                                                    sigmasq_varphi = sigmasq_varphi,
                                                    Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s


    # Sigma_update
    Sigma_update_s <- update_Sigma_NP_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ]
  )
  return(results)
}


#######################################################

start_time <- Sys.time()
res_rnorm <- fit_NP_MVN_cov_rnorm(niter = 10, burn_in = 2, thin = 1,
                                  n = n, K = K, Y = Y, W = W,
                                  WTW = WTW, n_all_par = n_all_par,
                                  J = J, M = M, O = O,
                                  theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                                  Sigma_init = Sigma_init,
                                  nu_0 = nu_0, Psi_0 = Psi_0,
                                  sigmasq_varphi = 10)
end_time <- Sys.time()
time_rnorm <- end_time - start_time #6.781877 secs


#################################################################

# Update theta
update_theta_NP_MVN_cov <- function(Y, K, W, WTW, n_all_par,
                                          J, M, O,
                                          sigmasq_varphi, Sigma_update)
{
  # Build V_theta (same for all k)
  var_alpha_0 <- 100 # dim 1*1
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  # vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
  vec_sigmasq_varphi <- rep(10, times = O)
  mat_var_varphi <- diag(vec_sigmasq_varphi)

  V_theta <- as.matrix(bdiag(var_alpha_0, var_beta, var_gamma,
                             var_delta, mat_var_varphi))

  # prior precision
  V_theta_inv <- chol2inv(chol(V_theta))
  prior_prec <- kronecker(Diagonal(K), V_theta_inv)  # I_K ⊗ Vθ⁻¹

  # Kronecker term from likelihood
  Sigma_inv <- as.matrix(chol2inv(chol(Sigma_update)))
  # Sigma_inv <- as.matrix(chol2inv(chol(Sigma_init)))
  kron_term <- kronecker(Sigma_inv, WTW)  # Σ⁻¹ ⊗ WᵀW

  kr_pr <- kron_term + prior_prec
  # Full conditional covariance
  Sigma_theta <- chol2inv(chol(kr_pr))

  # Full conditional mean
  vecY <- as.vector(Y)
  mean_theta_vec <- Sigma_theta %*% (kronecker(Sigma_inv, t(W)) %*% vecY)


  # Sample from MVN
  theta_vec <- mvrnorm(1, mu = mean_theta_vec, Sigma = Sigma_theta)
  # Reshape back to n_all_par x K
  theta_update_s <- matrix(theta_vec, nrow = n_all_par, ncol = K)



  return(theta_update_s)
}


#######################################################################################

# Update Sigma
update_Sigma_NP_MVN_cov <- function(Y, n, W, theta_update,
                                    Psi_0, nu_0)
{
  # Residuals
  Resid <- Y - W %*% theta_update

  # Posterior scale and df
  Psi_n <- Psi_0 + t(Resid) %*% Resid
  nu_n <- nu_0 + n

  # Draw from IW
  Sigma_update_s <- rinvwishart(nu_n, Psi_n)

  return(Sigma_update_s)
}


#######################################################################################

# Update parameters for NP-MVN-cov model
fit_NP_MVN_cov <- function(niter = 20, burn_in = 2, thin = 1,
                                 n, K, Y, W,
                                 WTW, n_all_par,
                                 J, M, O,
                                 theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                                 Sigma_init = diag(K),
                                 nu_0, Psi_0,
                                 sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))

  # Initialize
  theta_update[1, , ] <- theta_init
  Sigma_update[1, , ] <- Sigma_init

  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_NP_MVN_cov(Y, K, W, WTW, n_all_par,
                                                    J, M, O,
                                                    sigmasq_varphi = sigmasq_varphi,
                                                    Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s


    # Sigma_update
    Sigma_update_s <- update_Sigma_NP_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ]
  )
  return(results)
}


#######################################################

start_time <- Sys.time()
res <- fit_NP_MVN_cov(niter = 10, burn_in = 2, thin = 1,
                                  n = n, K = K, Y = Y, W = W,
                                  WTW = WTW, n_all_par = n_all_par,
                                  J = J, M = M, O = O,
                                  theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                                  Sigma_init = Sigma_init,
                                  nu_0 = nu_0, Psi_0 = Psi_0,
                                  sigmasq_varphi = 10)
end_time <- Sys.time()
time_cov <- end_time - start_time #6.781877 secs



########################################################

# # Update theta
# # Build V_theta (same for all k)
# var_alpha_0 <- 100 # dim 1*1
# var_beta <- diag(rep(1, J)) # dim J*J
# var_gamma <- diag(rep(1, M)) # dim M*M
# var_delta <- diag(rep(1, J*M)) # dim JM*JM
#
# # vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
# vec_sigmasq_varphi <- rep(10, times = O)
# mat_var_varphi <- diag(vec_sigmasq_varphi)
#
# # Inverse using blocks separately
# V_theta_inv_block <- as.matrix(
#   bdiag(
#     solve(var_alpha_0),
#     solve(var_beta),
#     solve(var_gamma),
#     solve(var_delta),
#     solve(mat_var_varphi)
#   )
# )
#
#
# prior_prec <- kronecker(Diagonal(K), V_theta_inv_block)  # I_K ⊗ Vθ⁻¹
#
# # Sigma_inv <- as.matrix(chol2inv(chol(Sigma_update)))
# Sigma_inv <- as.matrix(chol2inv(chol(Sigma_init)))
# kron_term <- kronecker(Sigma_inv, WTW)  # Σ⁻¹ ⊗ WᵀW
# kr_pr <- kron_term + prior_prec
#
# vecY <- as.vector(Y)
#
# Q <- kr_pr
# Sigma_theta <- solve(Q)
#
# b <- kronecker(Sigma_inv, t(W)) %*% vecY
# mu1 <- Sigma_theta %*% b
#
# # Method 1
# theta1 <- MASS::mvrnorm(1, mu = mu1, Sigma = Sigma_theta)
#
# # Method 2
# R <- chol(Q)
# mu2 <- backsolve(R, forwardsolve(t(R), b))
# theta2 <- mu2 + backsolve(R, rnorm(length(b)))
#
# max(abs(mu1 - mu2))
# abs(mean(theta1) - mean(theta2))
# max(abs(Q %*% mu2 - b))
#
#
#
#
# nsim <- 1000
#
# draws1 <- replicate(
#   nsim,
#   MASS::mvrnorm(
#     1,
#     mu = as.vector(mu1),
#     Sigma = Sigma_theta
#   ),
#   simplify = "matrix"
# )
#
# dim(draws1)
#
# draws2 <- replicate(
#   nsim,
#   as.vector(mu2 + backsolve(R, rnorm(length(mu2))))
# )
#
# dim(draws2)
#
# max(abs(rowMeans(draws1) - rowMeans(draws2)))
# max(abs(cov(t(draws1)) - cov(t(draws2))))
#
#
#
#
# A <- backsolve(R, diag(nrow(R)))
# Sigma_from_Q <- chol2inv(chol(kr_pr))
# max(abs(A %*% t(A) - Sigma_from_Q))
#
#
#
# max(abs(kr_pr %*% mu2 - b))

##############################################################

# Update theta
update_theta_NP_MVN_cov_Q <- function(Y, K, W, WTW, n_all_par,
                                                J, M, O,
                                                sigmasq_varphi, Sigma_update)
{
  # Build V_theta (same for all k)
  var_alpha_0 <- 100 # dim 1*1
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  # vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
  vec_sigmasq_varphi <- rep(10, times = O)
  mat_var_varphi <- diag(vec_sigmasq_varphi)

  # Inverse using blocks separately
  V_theta_inv_block <- as.matrix(
    bdiag(
      solve(var_alpha_0),
      solve(var_beta),
      solve(var_gamma),
      solve(var_delta),
      solve(mat_var_varphi)
    )
  )

prior_prec <- kronecker(Diagonal(K), V_theta_inv_block)  # I_K ⊗ Vθ⁻¹

Sigma_inv <- as.matrix(chol2inv(chol(Sigma_update)))
kron_term <- kronecker(Sigma_inv, WTW)  # Σ⁻¹ ⊗ WᵀW

vecY <- as.vector(Y)
Q <- kron_term + prior_prec
b <- kronecker(Sigma_inv, t(W)) %*% vecY

R <- chol(Q)
mu_theta <- backsolve(R, forwardsolve(t(R), b))
theta_vec <- mu_theta + backsolve(R, rnorm(length(b)))
# Reshape back to n_all_par x K
theta_update_s <- matrix(theta_vec, nrow = n_all_par, ncol = K)

return(theta_update_s)
}

#######################################################################################

# Update Sigma
update_Sigma_NP_MVN_cov <- function(Y, n, W, theta_update,
                                    Psi_0, nu_0)
{
  # Residuals
  Resid <- Y - W %*% theta_update

  # Posterior scale and df
  Psi_n <- Psi_0 + t(Resid) %*% Resid
  nu_n <- nu_0 + n

  # Draw from IW
  Sigma_update_s <- rinvwishart(nu_n, Psi_n)

  return(Sigma_update_s)
}


#######################################################################################

# Update parameters for NP-MVN-cov model
fit_NP_MVN_cov_Q <- function(niter = 20, burn_in = 2, thin = 1,
                                       n, K, Y, W,
                                       WTW, n_all_par,
                                       J, M, O,
                                       theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                                       Sigma_init = diag(K),
                                       nu_0, Psi_0,
                                       sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))

  # Initialize
  theta_update[1, , ] <- theta_init
  Sigma_update[1, , ] <- Sigma_init

  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_NP_MVN_cov_Q(Y, K, W, WTW, n_all_par,
                                                          J, M, O,
                                                          sigmasq_varphi = sigmasq_varphi,
                                                          Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s


    # Sigma_update
    Sigma_update_s <- update_Sigma_NP_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ]
  )
  return(results)
}


#######################################################

start_time <- Sys.time()
res_Q <- fit_NP_MVN_cov_Q(niter = 10, burn_in = 2, thin = 1,
                                              n = n, K = K, Y = Y, W = W,
                                              WTW = WTW, n_all_par = n_all_par,
                                              J = J, M = M, O = O,
                                              theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                                              Sigma_init = Sigma_init,
                                              nu_0 = nu_0, Psi_0 = Psi_0,
                                              sigmasq_varphi = 10)
end_time <- Sys.time()
time_Q <- end_time - start_time #6.842803 secs

#################################################################




time_cov #23.22538 mins
time_rnorm #4.568626 mins
time_block_rnorm # 4.580743 mins
time_Q #1.855866 mins


#################################################################

# Posterior mean
mean_theta_cov  <- apply(res$theta_update, c(2,3), mean)
mean_theta_rnorm <- apply(res_rnorm$theta_update, c(2,3), mean)
mean_theta_block <- apply(res_block_rnorm$theta_update, c(2,3), mean)
mean_theta_Q     <- apply(res_Q$theta_update, c(2,3), mean)

max(abs(mean_theta_cov - mean_theta_rnorm))
max(abs(mean_theta_cov - mean_theta_block))
max(abs(mean_theta_cov - mean_theta_Q))





set.seed(1)
s1 <- update_theta_NP_MVN_cov(Y, K, W, WTW, n_all_par, J, M, O,
                              sigmasq_varphi = 10,
                              Sigma_update = Sigma_init)

set.seed(1)
s2 <- update_theta_NP_MVN_cov_rnorm(Y, K, W, WTW, n_all_par, J, M, O,
                                    sigmasq_varphi = 10,
                                    Sigma_update = Sigma_init)

set.seed(1)
s3 <- update_theta_NP_MVN_cov_Q(Y, K, W, WTW, n_all_par, J, M, O,
                                sigmasq_varphi = 10,
                                Sigma_update = Sigma_init)

max(abs(s1 - s2))
max(abs(s1 - s3))



acf(res$theta_update[,1,1])
acf(res_rnorm$theta_update[,1,1])





