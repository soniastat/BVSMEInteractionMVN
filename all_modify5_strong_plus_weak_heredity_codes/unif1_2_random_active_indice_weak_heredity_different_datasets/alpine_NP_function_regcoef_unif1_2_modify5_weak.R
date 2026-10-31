
# Functions for simulating data and fitting HHS-MVN-cov model

###################################################################

# Simulate data (with weak heredity)
Sim_data_BVS_real <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                              x, z, d, Sigma, weak_per = 0.1)
{
  # modifiers
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- rnorm(K)
  beta_true[2, ] <- rnorm(K)
  active_beta  <- which(rowSums(abs(beta_true)) > 0)

  # environmental exposures
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different gamma for different k
  gamma_true[1, ] <- runif(K, 1, 2)
  gamma_true[2, ] <- runif(K, 1, 2)
  gamma_true[3, ] <- runif(K, 1, 2)
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
    runif(length(active_delta) * K, 1, 2),
    nrow = length(active_delta),
    ncol = K
  )

  # covariates
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- runif(K, 1, 2)
  varphi_true[2, ] <- runif(K, 1, 2)

  alpha_true <- runif(K, 1, 2)

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





#########################################################
# Fit NP model
########################################################

# Update theta
update_theta_NP_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par, J, M, O,
                                           sigmasq_alpha = 100,
                                           sigmasq_varphi, Sigma_update)
{

  # Build V_theta (same for all k)
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)

  # Inverse using blocks separately
  V_theta_inv <- as.matrix(
    bdiag(
      1/sigmasq_alpha,
      diag(1/diag(var_beta)),
      diag(1/diag(var_gamma)),
      diag(1/diag(var_delta)),
      diag(1/vec_sigmasq_varphi)
    )
  )

  prior_prec <- kronecker(Diagonal(K), V_theta_inv)  # I_K ⊗ Vθ⁻¹
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

#' Update parameters for NP-MVN-cov model
#' @export
fit_NP_MVN_cov_modify3 <- function(niter = 6000, burn_in = 1000, thin = 5,
                                   n, K, Y, W, WTW, n_all_par, J, M, O,
                                   sigmasq_alpha = 100,
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
                                                     sigmasq_alpha = sigmasq_alpha,
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


