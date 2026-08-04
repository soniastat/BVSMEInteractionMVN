

# Simulate data
Sim_data_BVS_real <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                              x, z, d, Sigma)
{
  # Hyperpriors
  tau_beta  <- 0.7
  tau_gamma <- 0.8
  tau_delta <- 0.5

  # Local shrinkage parameters
  lambda_beta  <- rep(0, J)
  lambda_gamma <- rep(0, M)

  # Active modifier main effects
  lambda_beta[1:2] <- c(1.0, 1.2)

  # Active exposure main effects
  lambda_gamma[1:3] <- c(1.0, 1.1, 1.3)

  # Modifier effects
  beta_true <- matrix(0, J, K)
  for(j in 1:J)
  {
   beta_true[j, ] <- rnorm(K, mean = 0, sd = tau_beta * lambda_beta[j])
  }

  # Exposure effects
  gamma_true <- matrix(0, M, K)
  for(m in 1:M)
  {
   gamma_true[m, ] <- rnorm(K, mean = 0, sd = tau_gamma * lambda_gamma[m])
  }

  # Interaction design matrix
  u <- matrix(0, n, J * M)
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      u[, col_idx] <- x[, j] * z[, m]
      col_idx <- col_idx + 1
    }
  }

 # Interaction coefficients
 delta_true <- matrix(0, J*M, K)
 # Interaction-specific local scales
 lambda_delta <- rep(0, J*M)
 for(j in 1:J)
 {
  beta_active <- lambda_beta[j] > 0
  for(m in 1:M)
  {
   gamma_active <- lambda_gamma[m] > 0
   idx <- (j-1)*M + m
   # Exact strong heredity
   if(beta_active & gamma_active)
   {
    lambda_delta[idx] <- rgamma(1, 2, 3)

    sd_delta <- (tau_delta * lambda_beta[j] *
                   lambda_gamma[m] *
                   lambda_delta[idx])
    delta_true[idx, ] <- rnorm(K, mean = 0, sd = sd_delta)
    }
   }
  }

  # covariates
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- seq(0.25, 1.25, length.out = K)
  varphi_true[2, ] <- seq(-1.25, -0.25, length.out = K)

  alpha_true <- runif(K)

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
  raw_Y <- W_theta_true + E_mat %*% R_Sigma

  # scaled the simulated Y
  Y <- scale(raw_Y)

  scale_sd <- apply(raw_Y, 2, sd)
  alpha_scaled <- alpha_true/scale_sd
  beta_scaled <- beta_true / scale_sd
  gamma_scaled <- gamma_true / scale_sd
  delta_scaled <- delta_true / scale_sd
  varphi_scaled <- varphi_true / scale_sd
  theta_scaled <- theta_true / scale_sd

  return(list(
    x = x, z = z, u = u, d = d, W = W, Y = Y,
    alpha_true = alpha_scaled, beta_true = beta_scaled,
    gamma_true = delta_scaled, delta_true = delta_scaled,
    varphi_true = varphi_scaled, W_theta_true = W_theta_true,
    Sigma = Sigma, theta_true = theta_scaled
  ))
}


