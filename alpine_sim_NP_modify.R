
# Functions for simulating data and fitting BVS-NP-MVN-cov model


# Simulate data
Sim_data_BVS_real <- function(K = 4, n = 1000, J = 44, M = 17, O = 6,
                              x, z, d)
{
  # modifiers
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- seq(0.25, 1.25, length.out = K)
  beta_true[2, ] <- seq(-1.25, -.25, length.out = K)

  # environmental exposures
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different gamma for different k
  gamma_true[1, ] <- seq(0.5, 1.5, length.out = K)
  gamma_true[2, ] <- seq(-1.5, -0.5, length.out = K)
  gamma_true[3, ] <- seq(-1.5, 2, length.out = K)

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
  delta_true[1:3, ] <- rnorm(3*K)
  delta_true[(M+1):(M+3), ] <- rnorm(3*K)

  # covariates
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- seq(0.25, 1.25, length.out = K)
  varphi_true[2, ] <- seq(-1.25, -.25, length.out = K)

  alpha_true <- runif(K)

  theta_true <- rbind(matrix(alpha_true, nrow = 1),
                      beta_true, gamma_true, delta_true,
                      varphi_true)
  W <- cbind(1, x, z, u, d)
  W_theta_true <- W %*% theta_true # n*K dimension

  # Wishart-based PD matrix
  A <- matrix(rnorm(K*K), K, K)
  Sigma <- crossprod(A)   # AᵀA is automatically positive definite


  # Cholesky of Sigma
  R_Sigma <- chol(Sigma)

  # E ~ iid N(0,1) error
  E_mat <- matrix(rnorm(n * K), nrow = n, ncol = K)

  # Y = M + E * chol(Sigma)
  Y <- W_theta_true + E_mat %*% R_Sigma

  return(list(
    x = x, z = z, u = u, d = d, W = W, Y = Y,
    alpha_true = alpha_true, beta_true = beta_true,
    gamma_true = gamma_true, delta_true = delta_true,
    varphi_true, W_theta_true = W_theta_true,
    Sigma = Sigma, theta_true = theta_true
  ))
}







#########################################################
# Fit NP model
########################################################

# Update theta
update_theta_NP_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par, J, M, O,
                                    sigmasq_varphi, Sigma_update)
{

  # Build V_theta (same for all k)
  var_alpha_0 <- 100 # dim 1*1
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
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
update_Sigma_NP_MVN_cov <- function(Y, n, W, theta_update, Psi_0, nu_0)
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
fit_NP_MVN_cov_modify <- function(niter = 20, burn_in = 2, thin = 1,
                           n, K, Y, W, WTW, n_all_par, J, M, O,
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
    theta_update_s <- update_theta_NP_MVN_cov_modify(Y, K, W, WTW, n_all_par, J, M, O,
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







