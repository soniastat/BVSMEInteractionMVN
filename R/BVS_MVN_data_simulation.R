
# Data simulation for the Bayesian variable selection model for
# modifier-interaction model for multivariate normal (MVN) response


################################################################################
# Data simulation
################################################################################
#' Smulate data
#' @export
Sim_data_BVS <- function(K = 5, n = 1000, J = 17, M = 50)
{
  # modifiers
  x <- matrix(runif(n * J), nrow = n)
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- seq(0.25, 1.25, length.out = K)
  beta_true[2, ] <- seq(-1.25, -.25, length.out = K)

  # environmental variables
  z <- matrix(rnorm(n * M), nrow = n)
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different beta for different k
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

  alpha_true <- runif(K)

  theta_true <- rbind(matrix(alpha_true, nrow = 1), beta_true, gamma_true, delta_true)
  W <- cbind(1, x, z, u)
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
    x = x, z = z, u = u, W = W, Y = Y,
    alpha_true = alpha_true, beta_true = beta_true,
    gamma_true = gamma_true, delta_true = delta_true,
    W_theta_true = W_theta_true, Sigma = Sigma, theta_true = theta_true
  ))
}




####################################################################

#' Smulate data
#' @export
Sim_data_BVS_real <- function(K = 5, n = 1000, J = 44, M = 17, O = 6,
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

