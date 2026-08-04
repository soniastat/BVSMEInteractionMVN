
# Simulate data
Sim_data_BVS_real <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                              x, z, d, Sigma)
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

  # modifier-environment interaction terms
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
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- seq(0.25, 1.25, length.out = K)
  varphi_true[2, ] <- seq(-1.25, -.25, length.out = K)

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
  Y <- W_theta_true + E_mat %*% R_Sigma

  return(list(
    x = x, z = z, u = u, d = d, W = W, Y = Y,
    alpha_true = alpha_true, beta_true = beta_true,
    gamma_true = gamma_true, delta_true = delta_true,
    varphi_true, W_theta_true = W_theta_true,
    Sigma = Sigma, theta_true = theta_true
  ))
}


####################################################


# Simulate data
Sim_data_BVS_real <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                              x, z, d, Sigma)
{
  # Hyperpriors
  tau_beta  <- 0.5
  tau_gamma <- 0.6
  tau_delta <- 0.4

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
  beta_active <- max(abs(beta_true[j,])) > 0
  for(m in 1:M)
  {
   gamma_active <- max(abs(gamma_true[m,])) > 0
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
  varphi_true[2, ] <- seq(-1.25, -.25, length.out = K)

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
  Y <- W_theta_true + E_mat %*% R_Sigma

  return(list(
    x = x, z = z, u = u, d = d, W = W, Y = Y,
    alpha_true = alpha_true, beta_true = beta_true,
    gamma_true = gamma_true, delta_true = delta_true,
    varphi_true, W_theta_true = W_theta_true,
    Sigma = Sigma, theta_true = theta_true
  ))
}



############################
# Load data
############################

load("exposome.RData")
load("modifiers_covariates.rda")
load("exposures.rda")
load("cov_base_post.rda")



raw_Y <- cbind(
  phenotype[,5:6],
  cov_base_post
)

Y0 <- scale(raw_Y)

colMeans(Y0)

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
# Simulate one dataset
############################

set.seed(23444)
sim_data <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D,
                              Sigma=Sigma)

Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
d <- sim_data$d

W <- cbind(1, x, z, u, d)
WTW <- t(W) %*% W

n_all_par <- ncol(W)

matplot(Y, type = "l")
colMeans(Y) # 1.008969  4.796786 11.129089 14.016009
colMeans(Y) # -1.801501  2.569311 -3.499702 -1.004731





