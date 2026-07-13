




library(Matrix)
library(LaplacesDemon)



####
#data
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

K <- dim(Y0)[2]
n <- nrow(Y0)
J <- ncol(X)
M <- ncol(Z)
O <- ncol(D)


lambdasq_beta_update <- rinvgamma(J, 2, 3)
lambdasq_gamma_update <- rinvgamma(M, 2, 3)
lambdasq_delta_update <- rinvgamma(J*M, 2, 3)
tausq_beta_update <- rgamma(1, 2, 3)
tausq_gamma_update <- rgamma(1, 2, 3)
tausq_delta_update <- rgamma(1, 2, 3)
sigmasq_varphi <- 10


sigmasq_beta_update <- lambdasq_beta_update * tausq_beta_update
sigmasq_gamma_update <- lambdasq_gamma_update * tausq_gamma_update
sigmasq_delta_update <- rep(NA, J * M)
for (j in 1:J)
{
  for (m in 1:M)
  {
    jm <- (j - 1) * M + m   # row index for (j,m)
    sigmasq_delta_update[jm] <- (lambdasq_delta_update[jm] *
                                   lambdasq_beta_update[j] *
                                   lambdasq_gamma_update[m] *
                                   tausq_delta_update)
  }
}

# Build V_theta (same for all k)
var_alpha_0 <- 100
var_beta <- diag(sigmasq_beta_update)
var_gamma <- diag(sigmasq_gamma_update)
var_delta <- diag(sigmasq_delta_update)

vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
mat_var_varphi <- diag(vec_sigmasq_varphi)



# compare different functions
V_theta_inv1 <- as.matrix(bdiag(
  chol2inv(chol(var_alpha_0)),
  chol2inv(chol(var_beta)),
  chol2inv(chol(var_gamma)),
  chol2inv(chol(var_delta)),
  chol2inv(chol(mat_var_varphi)) ))


V_theta_inv2 <- as.matrix(bdiag(
  1/var_alpha_0,
  chol2inv(chol(var_beta)),
  chol2inv(chol(var_gamma)),
  chol2inv(chol(var_delta)),
  chol2inv(chol(mat_var_varphi)) ))


# this may provide infinite number
V_theta_inv3 <- as.matrix(bdiag(
  1/var_alpha_0,
  diag(1/diag(var_beta)),
  diag(1/diag(var_gamma)),
  diag(1/diag(var_delta)),
  diag(1/diag(mat_var_varphi))))



V_theta_inv4 <- as.matrix(bdiag(
  solve(var_alpha_0),
  solve(var_beta),
  solve(var_gamma),
  solve(var_delta),
  solve(mat_var_varphi) ))




V_theta_inv5 <- chol2inv(chol(as.matrix(
  bdiag(var_alpha_0, var_beta,
        var_gamma, var_delta,
        mat_var_varphi))))






all.equal(V_theta_inv1, V_theta_inv2)
all.equal(V_theta_inv1, V_theta_inv3)
all.equal(V_theta_inv1, V_theta_inv4)
all.equal(V_theta_inv1, V_theta_inv5)

all.equal(V_theta_inv2, V_theta_inv3)
all.equal(V_theta_inv2, V_theta_inv4)
all.equal(V_theta_inv2, V_theta_inv5)

all.equal(V_theta_inv3, V_theta_inv4)
all.equal(V_theta_inv3, V_theta_inv5)

all.equal(V_theta_inv4, V_theta_inv5)


max(abs(V_theta_inv1 - V_theta_inv2))
max(abs(V_theta_inv1 - V_theta_inv3))
max(abs(V_theta_inv1 - V_theta_inv4))
max(abs(V_theta_inv1 - V_theta_inv5)) #0

max(abs(V_theta_inv2 - V_theta_inv3))
max(abs(V_theta_inv2 - V_theta_inv4))
max(abs(V_theta_inv2 - V_theta_inv5))

max(abs(V_theta_inv3 - V_theta_inv4))
max(abs(V_theta_inv3 - V_theta_inv5))

max(abs(V_theta_inv4 - V_theta_inv5))



range(eigen(V_theta_inv1)$values -
        eigen(V_theta_inv2)$values)
range(eigen(V_theta_inv1)$values -
        eigen(V_theta_inv3)$values)
range(eigen(V_theta_inv1)$values -
        eigen(V_theta_inv4)$values)
range(eigen(V_theta_inv1)$values -
        eigen(V_theta_inv5)$values) #0

range(eigen(V_theta_inv2)$values -
        eigen(V_theta_inv3)$values)
range(eigen(V_theta_inv2)$values -
        eigen(V_theta_inv4)$values)
range(eigen(V_theta_inv2)$values -
        eigen(V_theta_inv5)$values)

range(eigen(V_theta_inv3)$values -
        eigen(V_theta_inv4)$values)
range(eigen(V_theta_inv3)$values -
        eigen(V_theta_inv5)$values)

range(eigen(V_theta_inv4)$values -
        eigen(V_theta_inv5)$values)



# timing
install.packages("microbenchmark")
library(microbenchmark)



mb <- microbenchmark(

  # chol = {
  #   V_theta_inv1 <- as.matrix(bdiag(
  #     chol2inv(chol(var_alpha_0)),
  #     chol2inv(chol(var_beta)),
  #     chol2inv(chol(var_gamma)),
  #     chol2inv(chol(var_delta)),
  #     chol2inv(chol(mat_var_varphi))
  #   ))
  # },

  reciprocal_chol <- {
    V_theta_inv2 <- as.matrix(bdiag(
      1/var_alpha_0,
      chol2inv(chol(var_beta)),
      chol2inv(chol(var_gamma)),
      chol2inv(chol(var_delta)),
      chol2inv(chol(mat_var_varphi))
    ))
  },

  reciprocal = {
    V_theta_inv3 <- as.matrix(bdiag(
      1/var_alpha_0,
      diag(1/diag(var_beta)),
      diag(1/diag(var_gamma)),
      diag(1/diag(var_delta)),
      diag(1/diag(mat_var_varphi))
    ))
  },


  # solve = {
  #   V_theta_inv4 <- as.matrix(bdiag(
  #     solve(var_alpha_0),
  #     solve(var_beta),
  #     solve(var_gamma),
  #     solve(var_delta),
  #     solve(mat_var_varphi)
  #   ))
  # },
  #
  # combined_chol = {
  #   V_theta_inv5 <- chol2inv(chol(as.matrix(
  #     bdiag(var_alpha_0, var_beta,
  #           var_gamma, var_delta,
  #           mat_var_varphi))))
  # },

  times = 6000
)

print(mb)



##############################################



library(tidyverse)

J_vals <- c(10, 15, 20, 25)

#----------------------------
# Main effects
#----------------------------

df_main <- tribble(
  ~Method, ~J, ~Sensitivity, ~Specificity, ~Precision,

  "HHS",10,48.00,99.38,97.74,
  "HHS",15,43.68,99.73,97.99,
  "HHS",20,40.08,99.77,97.84,
  "HHS",25,35.48,99.90,98.16,

  "HHS-SI",10,60.20,99.94,99.81,
  "HHS-SI",15,57.16,99.91,99.53,
  "HHS-SI",20,53.76,99.95,99.67,
  "HHS-SI",25,52.40,99.94,99.32,

  "HRHS",10,50.72,99.34,97.65,
  "HRHS",15,45.80,99.73,98.32,
  "HRHS",20,41.96,99.77,97.85,
  "HRHS",25,37.40,99.87,98.33,

  "HRHS-SI",10,59.20,99.86,99.51,
  "HRHS-SI",15,55.72,99.91,99.52,
  "HRHS-SI",20,52.32,99.93,99.45,
  "HRHS-SI",25,49.36,99.90,98.98,

  "NHHS",10,52.60,99.34,97.79,
  "NHHS",15,50.60,99.69,98.25,
  "NHHS",20,49.16,99.83,98.66,
  "NHHS",25,47.48,99.88,98.62,

  "NHHS-SI",10,59.12,99.90,99.68,
  "NHHS-SI",15,58.36,99.92,99.63,
  "NHHS-SI",20,55.92,99.97,99.81,
  "NHHS-SI",25,56.20,99.93,99.32,

  "NHRHS",10,54.24,99.34,97.80,
  "NHRHS",15,51.76,99.72,98.26,
  "NHRHS",20,50.28,99.75,98.14,
  "NHRHS",25,47.60,99.82,98.09,

  "NHRHS-SI",10,59.48,99.84,99.45,
  "NHRHS-SI",15,57.96,99.88,99.39,
  "NHRHS-SI",20,56.08,99.94,99.63,
  "NHRHS-SI",25,54.12,99.92,99.28,

  "NP",10,40.84,97.82,91.25,
  "NP",15,31.84,97.71,82.08,
  "NP",20,26.36,97.75,75.64,
  "NP",25,NA,NA,NA
)

df_main_long <-
  pivot_longer(
    df_main,
    cols = c(Sensitivity, Specificity, Precision),
    names_to = "Metric",
    values_to = "Value"
  )




library(ggplot2)

p_main <-
  ggplot(
    df_main_long,
    aes(
      x = J,
      y = Value,
      color = Method,
      group = Method
    )
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_x_continuous(breaks = c(10,15,20,25)) +
  labs(
    x = "Number of Modifiers (J)",
    y = "Percentage",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

p_main



#########

df_interaction <- tribble(
  ~Method, ~J, ~Sensitivity, ~Specificity, ~Precision,

  "HHS",10,30.13,99.98,99.47,
  "HHS",15,30.53,99.98,99.30,
  "HHS",20,28.17,99.99,99.41,
  "HHS",25,27.27,99.99,99.28,

  "HHS-SI",10,42.17,100.00,99.93,
  "HHS-SI",15,44.50,99.99,99.80,
  "HHS-SI",20,41.10,100.00,99.85,
  "HHS-SI",25,41.93,99.92,99.92,

  "HRHS",10,21.23,100.00,100.00,
  "HRHS",15,22.80,99.99,99.58,
  "HRHS",20,21.20,100.00,99.92,
  "HRHS",25,19.63,99.99,99.59,

  "HRHS-SI",10,40.83,100.00,100.00,
  "HRHS-SI",15,42.93,100.00,100.00,
  "HRHS-SI",20,40.00,99.99,99.79,
  "HRHS-SI",25,40.73,100.00,99.92,

  "NHHS",10,30.00,99.99,99.69,
  "NHHS",15,30.23,100.00,99.89,
  "NHHS",20,28.83,99.99,99.54,
  "NHHS",25,27.13,99.99,99.56,

  "NHHS-SI",10,41.57,99.98,99.69,
  "NHHS-SI",15,43.90,99.99,99.72,
  "NHHS-SI",20,40.57,99.99,99.46,
  "NHHS-SI",25,41.00,100.00,99.92,

  "NHRHS",10,33.13,99.95,98.86,
  "NHRHS",15,34.73,99.94,98.15,
  "NHRHS",20,33.00,99.94,97.62,
  "NHRHS",25,32.07,99.96,97.64,

  "NHRHS-SI",10,42.27,99.98,99.64,
  "NHRHS-SI",15,43.87,99.99,99.67,
  "NHRHS-SI",20,41.30,99.97,98.78,
  "NHRHS-SI",25,40.80,99.98,99.31,

  "NP",10,41.33,97.48,69.48,
  "NP",15,42.40,97.50,60.62,
  "NP",20,36.83,97.79,52.72,
  "NP",25,NA,NA,NA
)

df_interaction_long <-
  pivot_longer(
    df_interaction,
    cols = c(Sensitivity, Specificity, Precision),
    names_to = "Metric",
    values_to = "Value"
  )



p_interaction <-
  ggplot(
    df_interaction_long,
    aes(
      x = J,
      y = Value,
      color = Method,
      group = Method
    )
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_x_continuous(breaks = c(10,15,20,25)) +
  labs(
    x = "Number of Modifiers (J)",
    y = "Percentage",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

p_interaction



#################################################



# Functions for simulating data and fitting BVS-HHS-SI-MVN-cov model

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
# Fit HHS-SI model
########################################################

############################################################


# Update theta
update_theta_BVS_HHS_SI_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par,
                                                   J, M, O,
                                                   lambdasq_beta_update, tausq_beta_update,
                                                   lambdasq_gamma_update, tausq_gamma_update,
                                                   lambdasq_delta_update, tausq_delta_update,
                                                   sigmasq_varphi, Sigma_update
)
{
  sigmasq_beta_update <- lambdasq_beta_update * tausq_beta_update
  sigmasq_gamma_update <- lambdasq_gamma_update * tausq_gamma_update

  sigmasq_delta_update <- rep(NA, J * M)
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m   # row index for (j,m)
      sigmasq_delta_update[jm] <- (lambdasq_delta_update[jm] *
                                     lambdasq_beta_update[j] *
                                     lambdasq_gamma_update[m] *
                                     tausq_delta_update)
    }
  }


  # Build V_theta (same for all k)
  var_alpha_0 <- 100
  var_beta <- diag(sigmasq_beta_update)
  var_gamma <- diag(sigmasq_gamma_update)
  var_delta <- diag(sigmasq_delta_update)

  vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
  mat_var_varphi <- diag(vec_sigmasq_varphi)

  # Inverse using blocks separately
  V_theta_inv <- as.matrix(
    bdiag(
      1/var_alpha_0,
      diag(1/diag(var_beta)),
      diag(1/diag(var_gamma)),
      diag(1/diag(var_delta)),
      diag(1/diag(mat_var_varphi))
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
update_Sigma_BVS_HHS_SI_MVN_cov <- function(Y, n, W, theta_update, Psi_0, nu_0)
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

# Update lambdasq_beta (for all j = 1, 2, ..., J)
update_lambdasq_beta_BVS_HHS_SI_MVN_cov1 <- function(J, K, beta_update, psi_beta_update, tausq_beta_update)
{
  lambdasq_beta_update_s <- rep(NA, J)
  shape_lambdasq_beta <- (K+1)/2

  for(j in 1:J)
  {
    scale_lambdasq_beta <- ((1/psi_beta_update[j]) +
                              (sum((beta_update[j, ])^2) / (2 * tausq_beta_update)))
    lambdasq_beta_update_s[j] <- rinvgamma(1, shape = shape_lambdasq_beta, scale = scale_lambdasq_beta) # IG distribution
  }
  return(lambdasq_beta_update_s)
}





update_lambdasq_beta_BVS_HHS_SI_MVN_cov2 <- function(J, K, beta_update,
                                                     psi_beta_update,
                                                     tausq_beta_update)
{
  shape_lambdasq_beta <- (K + 1) / 2

  scale_lambdasq_beta <- (1 / psi_beta_update) +
    rowSums(beta_update^2) / (2 * tausq_beta_update)

  rinvgamma(J,
            shape = shape_lambdasq_beta,
            scale = scale_lambdasq_beta)
}
#######################################################################################

# Update lambdasq_gamma (for all m = 1, 2, ..., M)
update_lambdasq_gamma_BVS_HHS_SI_MVN_cov1 <- function(M, K, gamma_update, psi_gamma_update, tausq_gamma_update)
{
  lambdasq_gamma_update_s <- rep(NA, M)
  shape_lambdasq_gamma <- (K+1)/2

  for(m in 1:M)
  {
    scale_lambdasq_gamma <- ((1/psi_gamma_update[m]) +
                               (sum((gamma_update[m, ])^2) / (2 * tausq_gamma_update)))
    lambdasq_gamma_update_s[m] <- rinvgamma(1, shape = shape_lambdasq_gamma, scale = scale_lambdasq_gamma) # IG distribution
  }
  return(lambdasq_gamma_update_s)
}




update_lambdasq_gamma_BVS_HHS_SI_MVN_cov2 <- function(M, K, gamma_update,
                                                      psi_gamma_update,
                                                      tausq_gamma_update)
{
  shape_lambdasq_gamma <- (K + 1) / 2

  scale_lambdasq_gamma <- (1 / psi_gamma_update) +
    rowSums(gamma_update^2) / (2 * tausq_gamma_update)

  rinvgamma(M,
            shape = shape_lambdasq_gamma,
            scale = scale_lambdasq_gamma)
}
#####################################################################################

# Update lambdasq_delta (for all j and m)
update_lambdasq_delta_BVS_HHS_SI_MVN_cov <- function(J, M, K, delta_update, lambdasq_beta_update,
                                                     lambdasq_gamma_update, psi_delta_update, tausq_delta_update)
{
  lambdasq_delta_update_s <- rep(NA, J*M)
  shape_lambdasq_delta <- (K+1)/2

  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction
      scale_lambdasq_delta <- ((1/psi_delta_update[jm]) + (sum((delta_update[jm, ])^2) /
                                                             (2 * lambdasq_beta_update[j] * lambdasq_gamma_update[m] * tausq_delta_update)))

      lambdasq_delta_update_s[jm] <- rinvgamma(1, shape = shape_lambdasq_delta, scale = scale_lambdasq_delta) # IG distribution
    }
  }
  return(lambdasq_delta_update_s)
}

#######################################################################################

# Update tausq_beta
update_tausq_beta_BVS_HHS_SI_MVN_cov1 <- function(J, K, beta_update, lambdasq_beta_update, xi_beta_update)
{

  sum_term_tausq_beta <- 0
  for (j in 1:J)
  {
    for (k in 1:K)
    {
      sum_term_tausq_beta <- sum_term_tausq_beta + (beta_update[j, k]^2 / (2 * lambdasq_beta_update[j]))
    }
  }
  shape_tausq_beta <- (J*K + 1) / 2
  scale_tausq_beta <- ((1/xi_beta_update) + sum_term_tausq_beta)

  tausq_beta_update_s <- rinvgamma(1, shape = shape_tausq_beta, scale = scale_tausq_beta) # IG distribution

  return(tausq_beta_update_s)
}




update_tausq_beta_BVS_HHS_SI_MVN_cov2 <- function(
    J, K, beta_update,
    lambdasq_beta_update,
    xi_beta_update)
{
  rinvgamma(
    1,
    shape = (J * K + 1) / 2,
    scale = 1 / xi_beta_update +
      sum(rowSums(beta_update^2) / (2 * lambdasq_beta_update))
  )
}



#######################################################################################

# Update tausq_gamma
update_tausq_gamma_BVS_HHS_SI_MVN_cov1 <- function(M, K, gamma_update, lambdasq_gamma_update, xi_gamma_update)
{
  sum_term_tausq_gamma <- 0
  for(m in 1:M)
  {
    for (k in 1:K)
    {
      sum_term_tausq_gamma <- sum_term_tausq_gamma + (gamma_update[m, k]^2 /
                                                        (2 * lambdasq_gamma_update[m]))
    }
  }
  shape_tausq_gamma <- (M*K + 1) / 2
  scale_tausq_gamma <- ((1/xi_gamma_update) + sum_term_tausq_gamma)
  tausq_gamma_update_s <- rinvgamma(1, shape = shape_tausq_gamma, scale = scale_tausq_gamma) # IG distribution

  return(tausq_gamma_update_s)
}




update_tausq_gamma_BVS_HHS_SI_MVN_cov2 <- function(
    M, K, gamma_update,
    lambdasq_gamma_update,
    xi_gamma_update)
{
  rinvgamma(
    1,
    shape = (M * K + 1) / 2,
    scale = 1 / xi_gamma_update +
      sum(rowSums(gamma_update^2) / (2 * lambdasq_gamma_update))
  )
}


#######################################################################################

# Update tausq_delta
update_tausq_delta_BVS_HHS_SI_MVN_cov <- function(J, M, K, delta_update, lambdasq_delta_update,
                                                  lambdasq_beta_update, lambdasq_gamma_update, xi_delta_update)
{
  sum_term_tausq_delta <- 0
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction

      for (k in 1:K)
      {
        sum_term_tausq_delta <- sum_term_tausq_delta + ((delta_update[jm, k])^2 / (2 *
                                                                                     lambdasq_delta_update[jm] * lambdasq_beta_update[j] * lambdasq_gamma_update[m]))
      }
    }
  }
  shape_tausq_delta <- (J*M*K + 1) / 2
  scale_tausq_delta <- ((1/xi_delta_update) + sum_term_tausq_delta)
  tausq_delta_update_s <- rinvgamma(1, shape = shape_tausq_delta, scale = scale_tausq_delta) # IG distribution

  return(tausq_delta_update_s)
}

#######################################################################################

# Update psi_beta (for all j = 1, 2, ..., J)
update_psi_beta_BVS_HHS_SI_MVN_cov1 <- function(J, lambdasq_beta_update)
{
  psi_beta_update_s <- rep(NA, J)

  for(j in 1:J)
  {
    scale_psi_beta <- (1 + (1/lambdasq_beta_update[j]))
    psi_beta_update_s[j] <- rinvgamma(1, shape = 1, scale = scale_psi_beta) # IG distribution
  }
  return(psi_beta_update_s)
}




update_psi_beta_BVS_HHS_SI_MVN_cov2 <- function(J, lambdasq_beta_update)
{
  rinvgamma(
    n = J,
    shape = 1,
    scale = 1 + 1 / lambdasq_beta_update
  )
}


#######################################################################################

# Update psi_gamma (for all m = 1, 2, ..., M)
update_psi_gamma_BVS_HHS_SI_MVN_cov1 <- function(M, lambdasq_gamma_update)
{
  psi_gamma_update_s <- rep(NA, M)

  for(m in 1:M)
  {
    scale_psi_gamma <- (1 + (1/lambdasq_gamma_update[m]))
    psi_gamma_update_s[m] <- rinvgamma(1, shape = 1, scale = scale_psi_gamma) # IG distribution
  }
  return(psi_gamma_update_s)
}



update_psi_gamma_BVS_HHS_SI_MVN_cov2 <- function(M, lambdasq_gamma_update)
{
  rinvgamma(
    n = M,
    shape = 1,
    scale = 1 + 1 / lambdasq_gamma_update
  )
}


#######################################################################################

# Update psi_delta (for all j , m)
update_psi_delta_BVS_HHS_SI_MVN_cov1 <- function(J, M, K, lambdasq_delta_update)
{
  psi_delta_update_s <- rep(NA, J*M)

  for(jm in 1:(J*M))
  {
    scale_psi_delta <- (1 + (1/lambdasq_delta_update[jm]))
    psi_delta_update_s[jm] <- rinvgamma(1, shape = 1, scale = scale_psi_delta) # IG distribution
  }
  return(psi_delta_update_s)
}




update_psi_delta_BVS_HHS_SI_MVN_cov2 <- function(
    J, M,
    lambdasq_delta_update)
{
  rinvgamma(
    n = J * M,
    shape = 1,
    scale = 1 + 1 / lambdasq_delta_update
  )
}
#######################################################################################

# Update xi_beta
update_xi_beta_BVS_HHS_SI_MVN_cov <- function(J, tausq_beta_update)
{
  # Use more informative prior (xi_beta ~ IG(15, 3))
  shape_xi_beta <- ((1/2) + 15)
  scale_xi_beta <- (3 + (1 / tausq_beta_update))

  xi_beta_update_s <- rinvgamma(1, shape = shape_xi_beta, scale = scale_xi_beta) # IG distribution
  return(xi_beta_update_s)
}



#######################################################################################

# Update xi_gamma
update_xi_gamma_BVS_HHS_SI_MVN_cov <- function(M, tausq_gamma_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_gamma <- ((1/2) + 15)
  scale_xi_gamma <- (3 + (1 / tausq_gamma_update))

  xi_gamma_update_s <- rinvgamma(1, shape = shape_xi_gamma, scale = scale_xi_gamma) # IG distribution
  return(xi_gamma_update_s)
}



#######################################################################################

# Update xi_delta
update_xi_delta_BVS_HHS_SI_MVN_cov <- function(J, M, tausq_delta_update)
{
  # Use more informative prior (xi_delta ~ IG(15, 3))
  shape_xi_delta <- ((1/2) + 15)
  scale_xi_delta <- (3 + (1 / tausq_delta_update))

  xi_delta_update_s <- rinvgamma(1, shape = shape_xi_delta, scale = scale_xi_delta) # IG distribution
  return(xi_delta_update_s)
}


##################################################################################

# Update parameters for BVS-HHS-SI-MVN-modify model
fit_BVS_HHS_SI_MVN_cov_modify1 <- function(niter = 6000, burn_in = 1000, thin = 5,
                                           n, K, Y, W, WTW, n_all_par, J, M, O,
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
                                           Sigma_init = diag(K),
                                           nu_0, Psi_0,
                                           sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))
  lambdasq_beta_update <- matrix(NA, nrow = niter, ncol = J)
  tausq_beta_update <- rep(NA, niter)
  lambdasq_gamma_update <- matrix(NA, nrow = niter, ncol = M)
  tausq_gamma_update <- rep(NA, niter)
  lambdasq_delta_update <- matrix(NA, nrow = niter, ncol = J*M)
  tausq_delta_update <- rep(NA, niter)
  psi_beta_update <- matrix(NA, nrow = niter, ncol = J)
  psi_gamma_update <- matrix(NA, nrow = niter, ncol = M)
  psi_delta_update <- matrix(NA, nrow = niter, ncol = J*M)
  xi_beta_update <- rep(NA, niter)
  xi_gamma_update <- rep(NA, niter)
  xi_delta_update <- rep(NA, niter)

  alpha0_update <- matrix(NA, nrow = niter, ncol = K)
  beta_update <- array(NA, dim = c(niter, J, K))
  gamma_update <- array(NA, dim = c(niter, M, K))
  delta_update <- array(NA, dim = c(niter, J*M, K))
  varphi_update <- array(NA, dim = c(niter, O, K))

  # Initialize
  theta_update[1, , ] <- theta_init
  Sigma_update[1, , ] <- Sigma_init
  lambdasq_beta_update[1, ] <- lambdasq_beta_init
  tausq_beta_update[1] <- tausq_beta_init
  lambdasq_gamma_update[1, ] <- lambdasq_gamma_init
  tausq_gamma_update[1] <- tausq_gamma_init
  lambdasq_delta_update[1, ] <- lambdasq_delta_init
  tausq_delta_update[1] <- tausq_delta_init
  psi_beta_update[1, ] <- psi_beta_init
  psi_gamma_update[1, ] <- psi_gamma_init
  psi_delta_update[1, ] <- psi_delta_init
  xi_beta_update[1] <- xi_beta_init
  xi_gamma_update[1] <- xi_gamma_init
  xi_delta_update[1] <- xi_delta_init

  alpha0_update[1, ] <- theta_update[1, 1, ]
  beta_update[1, , ] <- theta_update[1, (1 + 1):(1 + J), ]
  gamma_update[1, , ] <- theta_update[1, (1 + J + 1):(1 + J + M), ]
  delta_update[1, , ] <- theta_update[1, (1 + J + M + 1):(1 + J + M + J*M), ]
  varphi_update[1, , ] <- theta_update[1, (1 + J + M + J*M + 1):(1 + J + M + J*M + O), ]



  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_BVS_HHS_SI_MVN_cov_modify(Y, K, W, WTW, n_all_par,
                                                             J, M, O,
                                                             lambdasq_beta_update = lambdasq_beta_update[(s-1), ],
                                                             tausq_beta_update = tausq_beta_update[(s-1)],
                                                             lambdasq_gamma_update = lambdasq_gamma_update[(s-1), ],
                                                             tausq_gamma_update = tausq_gamma_update[(s-1)],
                                                             lambdasq_delta_update = lambdasq_delta_update[(s-1), ],
                                                             tausq_delta_update = tausq_delta_update[(s-1)],
                                                             sigmasq_varphi = sigmasq_varphi,
                                                             Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s

    alpha0_update_s <- theta_update[s, 1, ]
    beta_update_s <- theta_update[s, (1 + 1):(1 + J), ]
    gamma_update_s <- theta_update[s, (1 + J + 1):(1 + J + M), ]
    delta_update_s <- theta_update[s, (1 + J + M + 1):(1 + J + M + J*M), ]
    varphi_update_s <- theta_update[s, (1 + J + M + J*M + 1):(1 + J + M + J*M + O), ]

    alpha0_update[s, ] <- alpha0_update_s
    beta_update[s, , ] <- beta_update_s
    gamma_update[s, , ] <- gamma_update_s
    delta_update[s, , ] <- delta_update_s
    varphi_update[s, , ] <- varphi_update_s

    # Sigma_update
    Sigma_update_s <- update_Sigma_BVS_HHS_SI_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s


    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_BVS_HHS_SI_MVN_cov1(J, K,
                                                                       beta_update = beta_update[s, , ],
                                                                       psi_beta_update = psi_beta_update[(s-1), ],
                                                                       tausq_beta_update = tausq_beta_update[(s-1)])
    lambdasq_beta_update[s, ] <- lambdasq_beta_update_s


    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_BVS_HHS_SI_MVN_cov1(M, K,
                                                                         gamma_update = gamma_update[s, , ],
                                                                         psi_gamma_update = psi_gamma_update[(s-1), ],
                                                                         tausq_gamma_update = tausq_gamma_update[(s-1)])
    lambdasq_gamma_update[s, ] <- lambdasq_gamma_update_s


    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_BVS_HHS_SI_MVN_cov(J, M, K,
                                                                        delta_update = delta_update[s, , ],
                                                                        lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                                        lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                                        psi_delta_update = psi_delta_update[(s-1), ],
                                                                        tausq_delta_update = tausq_delta_update[(s-1)])
    lambdasq_delta_update[s, ] <- lambdasq_delta_update_s



    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_BVS_HHS_SI_MVN_cov1(J, K,
                                                                 beta_update = beta_update[s, , ],
                                                                 lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                                 xi_beta_update = xi_beta_update[(s-1)])
    tausq_beta_update[s] <- tausq_beta_update_s


    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_BVS_HHS_SI_MVN_cov1(M, K,
                                                                   gamma_update = gamma_update[s, , ],
                                                                   lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                                   xi_gamma_update = xi_gamma_update[(s-1)])
    tausq_gamma_update[s] <- tausq_gamma_update_s


    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_BVS_HHS_SI_MVN_cov(J, M, K,
                                                                  delta_update = delta_update[s, , ],
                                                                  lambdasq_delta_update = lambdasq_delta_update[s, ],
                                                                  lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                                  lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                                  xi_delta_update = xi_delta_update[(s-1)])
    tausq_delta_update[s] <- tausq_delta_update_s



    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_BVS_HHS_SI_MVN_cov1(J,
                                                             lambdasq_beta_update = lambdasq_beta_update[s, ])
    psi_beta_update[s, ] <- psi_beta_update_s


    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_BVS_HHS_SI_MVN_cov1(M,
                                                               lambdasq_gamma_update = lambdasq_gamma_update[s, ])
    psi_gamma_update[s, ] <- psi_gamma_update_s


    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_BVS_HHS_SI_MVN_cov1(J, M,
                                                               lambdasq_delta_update = lambdasq_delta_update[s, ])
    psi_delta_update[s, ] <- psi_delta_update_s



    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_BVS_HHS_SI_MVN_cov(J, tausq_beta_update = tausq_beta_update[s])
    xi_beta_update[s] <- xi_beta_update_s


    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_BVS_HHS_SI_MVN_cov(M, tausq_gamma_update = tausq_gamma_update[s])
    xi_gamma_update[s] <- xi_gamma_update_s


    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_BVS_HHS_SI_MVN_cov(J, M, tausq_delta_update = tausq_delta_update[s])
    xi_delta_update[s] <- xi_delta_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ],
    lambdasq_beta_update = lambdasq_beta_update[indices_to_save, ],
    tausq_beta_update = tausq_beta_update[indices_to_save],
    lambdasq_gamma_update = lambdasq_gamma_update[indices_to_save, ],
    tausq_gamma_update = tausq_gamma_update[indices_to_save],
    lambdasq_delta_update = lambdasq_delta_update[indices_to_save, ],
    tausq_delta_update = tausq_delta_update[indices_to_save],
    psi_beta_update = psi_beta_update[indices_to_save, ],
    psi_gamma_update = psi_gamma_update[indices_to_save, ],
    psi_delta_update = psi_delta_update[indices_to_save, ],
    xi_beta_update = xi_beta_update[indices_to_save],
    xi_gamma_update = xi_gamma_update[indices_to_save],
    xi_delta_update = xi_delta_update[indices_to_save]
  )
  return(results)
}
################################################################


# Update parameters for BVS-HHS-SI-MVN-modify model
fit_BVS_HHS_SI_MVN_cov_modify2 <- function(niter = 6000, burn_in = 1000, thin = 5,
                                           n, K, Y, W, WTW, n_all_par, J, M, O,
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
                                           Sigma_init = diag(K),
                                           nu_0, Psi_0,
                                           sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))
  lambdasq_beta_update <- matrix(NA, nrow = niter, ncol = J)
  tausq_beta_update <- rep(NA, niter)
  lambdasq_gamma_update <- matrix(NA, nrow = niter, ncol = M)
  tausq_gamma_update <- rep(NA, niter)
  lambdasq_delta_update <- matrix(NA, nrow = niter, ncol = J*M)
  tausq_delta_update <- rep(NA, niter)
  psi_beta_update <- matrix(NA, nrow = niter, ncol = J)
  psi_gamma_update <- matrix(NA, nrow = niter, ncol = M)
  psi_delta_update <- matrix(NA, nrow = niter, ncol = J*M)
  xi_beta_update <- rep(NA, niter)
  xi_gamma_update <- rep(NA, niter)
  xi_delta_update <- rep(NA, niter)

  alpha0_update <- matrix(NA, nrow = niter, ncol = K)
  beta_update <- array(NA, dim = c(niter, J, K))
  gamma_update <- array(NA, dim = c(niter, M, K))
  delta_update <- array(NA, dim = c(niter, J*M, K))
  varphi_update <- array(NA, dim = c(niter, O, K))

  # Initialize
  theta_update[1, , ] <- theta_init
  Sigma_update[1, , ] <- Sigma_init
  lambdasq_beta_update[1, ] <- lambdasq_beta_init
  tausq_beta_update[1] <- tausq_beta_init
  lambdasq_gamma_update[1, ] <- lambdasq_gamma_init
  tausq_gamma_update[1] <- tausq_gamma_init
  lambdasq_delta_update[1, ] <- lambdasq_delta_init
  tausq_delta_update[1] <- tausq_delta_init
  psi_beta_update[1, ] <- psi_beta_init
  psi_gamma_update[1, ] <- psi_gamma_init
  psi_delta_update[1, ] <- psi_delta_init
  xi_beta_update[1] <- xi_beta_init
  xi_gamma_update[1] <- xi_gamma_init
  xi_delta_update[1] <- xi_delta_init

  alpha0_update[1, ] <- theta_update[1, 1, ]
  beta_update[1, , ] <- theta_update[1, (1 + 1):(1 + J), ]
  gamma_update[1, , ] <- theta_update[1, (1 + J + 1):(1 + J + M), ]
  delta_update[1, , ] <- theta_update[1, (1 + J + M + 1):(1 + J + M + J*M), ]
  varphi_update[1, , ] <- theta_update[1, (1 + J + M + J*M + 1):(1 + J + M + J*M + O), ]



  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_BVS_HHS_SI_MVN_cov_modify(Y, K, W, WTW, n_all_par,
                                                             J, M, O,
                                                             lambdasq_beta_update = lambdasq_beta_update[(s-1), ],
                                                             tausq_beta_update = tausq_beta_update[(s-1)],
                                                             lambdasq_gamma_update = lambdasq_gamma_update[(s-1), ],
                                                             tausq_gamma_update = tausq_gamma_update[(s-1)],
                                                             lambdasq_delta_update = lambdasq_delta_update[(s-1), ],
                                                             tausq_delta_update = tausq_delta_update[(s-1)],
                                                             sigmasq_varphi = sigmasq_varphi,
                                                             Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s

    alpha0_update_s <- theta_update[s, 1, ]
    beta_update_s <- theta_update[s, (1 + 1):(1 + J), ]
    gamma_update_s <- theta_update[s, (1 + J + 1):(1 + J + M), ]
    delta_update_s <- theta_update[s, (1 + J + M + 1):(1 + J + M + J*M), ]
    varphi_update_s <- theta_update[s, (1 + J + M + J*M + 1):(1 + J + M + J*M + O), ]

    alpha0_update[s, ] <- alpha0_update_s
    beta_update[s, , ] <- beta_update_s
    gamma_update[s, , ] <- gamma_update_s
    delta_update[s, , ] <- delta_update_s
    varphi_update[s, , ] <- varphi_update_s

    # Sigma_update
    Sigma_update_s <- update_Sigma_BVS_HHS_SI_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s


    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_BVS_HHS_SI_MVN_cov2(J, K,
                                                                       beta_update = beta_update[s, , ],
                                                                       psi_beta_update = psi_beta_update[(s-1), ],
                                                                       tausq_beta_update = tausq_beta_update[(s-1)])
    lambdasq_beta_update[s, ] <- lambdasq_beta_update_s


    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_BVS_HHS_SI_MVN_cov2(M, K,
                                                                         gamma_update = gamma_update[s, , ],
                                                                         psi_gamma_update = psi_gamma_update[(s-1), ],
                                                                         tausq_gamma_update = tausq_gamma_update[(s-1)])
    lambdasq_gamma_update[s, ] <- lambdasq_gamma_update_s


    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_BVS_HHS_SI_MVN_cov(J, M, K,
                                                                        delta_update = delta_update[s, , ],
                                                                        lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                                        lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                                        psi_delta_update = psi_delta_update[(s-1), ],
                                                                        tausq_delta_update = tausq_delta_update[(s-1)])
    lambdasq_delta_update[s, ] <- lambdasq_delta_update_s



    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_BVS_HHS_SI_MVN_cov2(J, K,
                                                                 beta_update = beta_update[s, , ],
                                                                 lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                                 xi_beta_update = xi_beta_update[(s-1)])
    tausq_beta_update[s] <- tausq_beta_update_s


    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_BVS_HHS_SI_MVN_cov2(M, K,
                                                                   gamma_update = gamma_update[s, , ],
                                                                   lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                                   xi_gamma_update = xi_gamma_update[(s-1)])
    tausq_gamma_update[s] <- tausq_gamma_update_s


    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_BVS_HHS_SI_MVN_cov(J, M, K,
                                                                  delta_update = delta_update[s, , ],
                                                                  lambdasq_delta_update = lambdasq_delta_update[s, ],
                                                                  lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                                  lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                                  xi_delta_update = xi_delta_update[(s-1)])
    tausq_delta_update[s] <- tausq_delta_update_s



    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_BVS_HHS_SI_MVN_cov2(J,
                                                             lambdasq_beta_update = lambdasq_beta_update[s, ])
    psi_beta_update[s, ] <- psi_beta_update_s


    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_BVS_HHS_SI_MVN_cov2(M,
                                                               lambdasq_gamma_update = lambdasq_gamma_update[s, ])
    psi_gamma_update[s, ] <- psi_gamma_update_s


    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_BVS_HHS_SI_MVN_cov2(J, M,
                                                               lambdasq_delta_update = lambdasq_delta_update[s, ])
    psi_delta_update[s, ] <- psi_delta_update_s



    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_BVS_HHS_SI_MVN_cov(J, tausq_beta_update = tausq_beta_update[s])
    xi_beta_update[s] <- xi_beta_update_s


    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_BVS_HHS_SI_MVN_cov(M, tausq_gamma_update = tausq_gamma_update[s])
    xi_gamma_update[s] <- xi_gamma_update_s


    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_BVS_HHS_SI_MVN_cov(J, M, tausq_delta_update = tausq_delta_update[s])
    xi_delta_update[s] <- xi_delta_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ],
    lambdasq_beta_update = lambdasq_beta_update[indices_to_save, ],
    tausq_beta_update = tausq_beta_update[indices_to_save],
    lambdasq_gamma_update = lambdasq_gamma_update[indices_to_save, ],
    tausq_gamma_update = tausq_gamma_update[indices_to_save],
    lambdasq_delta_update = lambdasq_delta_update[indices_to_save, ],
    tausq_delta_update = tausq_delta_update[indices_to_save],
    psi_beta_update = psi_beta_update[indices_to_save, ],
    psi_gamma_update = psi_gamma_update[indices_to_save, ],
    psi_delta_update = psi_delta_update[indices_to_save, ],
    xi_beta_update = xi_beta_update[indices_to_save],
    xi_gamma_update = xi_gamma_update[indices_to_save],
    xi_delta_update = xi_delta_update[indices_to_save]
  )
  return(results)
}


###########################################################





library(LaplacesDemon)
library(Matrix)



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

Z <- cbind(
  Organochlorines_Postnatal,
  Metals_Postnatal
)

K <- dim(Y0)[2]
n <- nrow(Y0)
J <- ncol(X)
M <- ncol(Z)
O <- ncol(D)


# devtools::install_github("soniastat/BVSMEInteractionMVN")
# library(BVSMEInteractionMVN)

set.seed(23444)
sim_data <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D)

Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
d <- sim_data$d

W <- cbind(1, x, z, u, d)
WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)

library(LaplacesDemon)
Sigma_init <- rinvwishart(nu_0, Psi_0)



############################
# Fit model
############################

start_time1 <- Sys.time()
set.seed(765878)
res1 <- tryCatch({
  fit_BVS_HHS_SI_MVN_cov_modify1(
    niter = 1000, burn_in = 500, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
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
}, error=function(e){
  list(error=TRUE, message=e$message)
})
end_time1 <- Sys.time()
end_time1 - start_time1







start_time2 <- Sys.time()
set.seed(765878)
res2 <- tryCatch({
  fit_BVS_HHS_SI_MVN_cov_modify2(
    niter = 2, burn_in = 1, thin = 1,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
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
}, error=function(e){
  list(error=TRUE, message=e$message)
})
end_time2 <- Sys.time()
end_time2 - start_time2





all.equal(res1$lambdasq_beta_update, res2$lambdasq_beta_update)
all.equal(res1$lambdasq_gamma_update, res2$lambdasq_gamma_update)
all.equal(res1$lambdasq_delta_update, res2$lambdasq_delta_update)
all.equal(res1$tausq_beta_update, res2$tausq_beta_update)
all.equal(res1$tausq_gamma_update, res2$tausq_gamma_update)
all.equal(res1$tausq_delta_update, res2$tausq_delta_update)
all.equal(res1$psi_beta_update, res2$psi_beta_update)
all.equal(res1$psi_gamma_update, res2$psi_gamma_update)
all.equal(res1$psi_delta_update, res2$psi_delta_update)
all.equal(res1$theta_update[,1,], res2$theta_update[,1,])
all.equal(res1$theta_update[,5,], res2$theta_update[,5,])

max(abs(res1$lambdasq_beta_update - res2$lambdasq_beta_update))
max(abs(res1$lambdasq_gamma_update - res2$lambdasq_gamma_update))
max(abs(res1$lambdasq_delta_update - res2$lambdasq_delta_update))
max(abs(res1$tausq_beta_update - res2$tausq_beta_update))
max(abs(res1$tausq_gamma_update - res2$tausq_gamma_update))
max(abs(res1$tausq_delta_update - res2$tausq_delta_update))
max(abs(res1$psi_beta_update - res2$psi_beta_update))
max(abs(res1$psi_gamma_update - res2$psi_gamma_update))
max(abs(res1$psi_delta_update - res2$psi_delta_update))
max(abs(res1$xi_beta_update - res2$xi_beta_update))
max(abs(res1$xi_gamma_update - res2$xi_gamma_update))
max(abs(res1$xi_delta_update - res2$xi_delta_update))
max(abs(res1$theta_update[,1,] - res2$theta_update[,1,]))
max(abs(res1$theta_update[,5,] - res2$theta_update[,5,]))
max(abs(res1$theta_update[,1,] - res2$theta_update[,1,]))
max(abs(res1$theta_update[,5,] - res2$theta_update[,5,]))






lambdasq_beta_update <- c(0.5, 1, 2, 4, 10)
scale_vec <- 1 + 1 / lambdasq_beta_update
scale_vec
set.seed(123)
x1 <- sapply(scale_vec, function(sc)
  rinvgamma(1, shape = 1, scale = sc))

set.seed(123)
x2 <- rinvgamma(length(scale_vec), shape = 1, scale = scale_vec)

max(abs(x1 - x2))







which.max(abs(res1$lambdasq_delta_update -
                res2$lambdasq_delta_update))

idx <- which.max(abs(res1$lambdasq_delta_update -
                       res2$lambdasq_delta_update))

res1$lambdasq_delta_update[idx]
res2$lambdasq_delta_update[idx]





rel_diff <- abs(res1$lambdasq_delta_update -
                  res2$lambdasq_delta_update) /
  pmax(abs(res1$lambdasq_delta_update),
       abs(res2$lambdasq_delta_update))

max(rel_diff)



summary(res1$lambdasq_delta_update)
quantile(res1$lambdasq_delta_update,
         probs = c(.5,.9,.95,.99,.999,.9999))
max(res1$lambdasq_delta_update)



sum(res1$lambdasq_delta_update > 1e6)
sum(res1$lambdasq_delta_update > 1e8)
sum(res1$lambdasq_delta_update > 1e10)



###########################################

sigmasq_delta_update <- matrix(NA, J * M, K)
for (j in 1:J) {
  for (m in 1:M) {
    jm <- (j - 1) * M + m   # row index for (j,m)
    for (k in 1:K)
    {
      sigmasq_delta_update[jm, k] <- (lambdasq_delta_update[jm, k] *
                                        lambdasq_beta_update[j, k] *
                                        lambdasq_gamma_update[m, k] *
                                        tausq_delta_update[jm])
    }
  }
}


beta_rep  <- lambdasq_beta_update[rep(1:J, each = M), ]
gamma_rep <- lambdasq_gamma_update[rep(1:M, times = J), ]
sigmasq_delta_update2 <- (lambdasq_delta_update * beta_rep *
                            gamma_rep *
                            matrix(
                              tausq_delta_update[1:(J*M)],
                              nrow = J*M,
                              ncol = K
                            ))


J <- 2
M <- 3
K <- 2
lambdasq_beta_update <- matrix(c(1, 2, 2, 3), ncol = 2)
lambdasq_gamma_update <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
lambdasq_delta_update <- matrix(c(1, 2, 3, 4, 5, 6,
                                  3, 4, 4, 1, 1, 3), ncol = 2)
tausq_delta_update <- c(1, 2, 3, 4, 5, 6)



#################################################


# Update theta
update_theta_BVS_HHS_MVN_cov <- function(Y, K, W, WTW, n_all_par,
                                         J, M, O,
                                         lambdasq_beta_update, tausq_beta_update,
                                         lambdasq_gamma_update, tausq_gamma_update,
                                         lambdasq_delta_update, tausq_delta_update,
                                         sigmasq_varphi, Sigma_update
)
{
  sigmasq_beta_update <- lambdasq_beta_update * tausq_beta_update
  sigmasq_gamma_update <- lambdasq_gamma_update * tausq_gamma_update

  sigmasq_delta_update <- matrix(NA, J * M, K)
  for (j in 1:J) {
    for (m in 1:M) {
      jm <- (j - 1) * M + m   # row index for (j,m)
      for (k in 1:K)
      {
        sigmasq_delta_update[jm, k] <- (lambdasq_delta_update[jm, k] *
                                          lambdasq_beta_update[j, k] *
                                          lambdasq_gamma_update[m, k] *
                                          tausq_delta_update[jm])
      }
    }
  }

  # Build list of prior covariance matrices V_theta_k for each k
  V_theta_list <- vector("list", K)
  for (k in 1:K)
  {
    # Gibbs update of theta
    var_alpha_0k <- 100
    mat_var_alpha_0k <- matrix(var_alpha_0k, nrow = 1, ncol = 1) #make it matrix
    var_beta_k <- diag(sigmasq_beta_update[, k]) # dim J*J
    var_gamma_k <- diag(sigmasq_gamma_update[, k]) # dim M*M
    var_delta_k <- diag(sigmasq_delta_update[, k]) # dim JM*JM

    vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
    mat_var_varphi <- diag(vec_sigmasq_varphi)

    V_theta_list[[k]] <- as.matrix(bdiag(mat_var_alpha_0k, var_beta_k, var_gamma_k,
                                         var_delta_k, mat_var_varphi))
  }
  # Build block-diagonal prior precision: blockdiag(V_theta1^{-1},...,V_thetaK^{-1})
  prior_prec_list <- lapply(V_theta_list, function(Vk) {
    chol2inv(chol(Vk))
  })

  prior_prec <- as.matrix(bdiag(prior_prec_list))   # pK x pK

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


###########

# Update theta
update_theta_BVS_HHS_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par,
                                                J, M, O,
                                                lambdasq_beta_update, tausq_beta_update,
                                                lambdasq_gamma_update, tausq_gamma_update,
                                                lambdasq_delta_update, tausq_delta_update,
                                                sigmasq_varphi, Sigma_update
)
{
  sigmasq_beta_update <- lambdasq_beta_update * tausq_beta_update
  sigmasq_gamma_update <- lambdasq_gamma_update * tausq_gamma_update

  beta_rep  <- lambdasq_beta_update[rep(1:J, each = M), ]
  gamma_rep <- lambdasq_gamma_update[rep(1:M, times = J), ]
  sigmasq_delta_update <- (lambdasq_delta_update * beta_rep *
                             gamma_rep *
                             matrix(tausq_delta_update[1:(J*M)],
                                    nrow = J*M,
                                    ncol = K))


  # Build list of prior covariance matrices V_theta_k for each k
  V_theta_k_inv_list <- vector("list", K)
  for (k in 1:K)
  {
    # Gibbs update of theta
    var_alpha_0k <- 100

    vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
    mat_var_varphi <- diag(vec_sigmasq_varphi)

    V_theta_k_inv_list[[k]] <- as.matrix(
      bdiag(
        1/var_alpha_0k,
        diag(1/sigmasq_beta_update[, k]),
        diag(1/sigmasq_gamma_update[, k]),
        diag(1/sigmasq_delta_update[, k]),
        diag(1/diag(mat_var_varphi))
      )
    )
  }
  prior_prec <- as.matrix(bdiag(V_theta_k_inv_list))   # pK x pK

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

####################

lambdasq_beta_update <- matrix(rinvgamma(J*K, 2, 3), ncol = K)
lambdasq_gamma_update <- matrix(rinvgamma(M*K, 2, 3), ncol = K)
lambdasq_delta_update <- matrix(rinvgamma(J*M*K, 2, 3), ncol = K)
tausq_beta_update <- rgamma(J, 2, 3)
tausq_gamma_update <- rgamma(M, 2, 3)
tausq_delta_update <- rgamma(J*M, 2, 3)
sigmasq_varphi <- 10

nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_update <- rinvwishart(nu_0, Psi_0)


set.seed(76789)
update_theta_modify <- update_theta_BVS_HHS_MVN_cov_modify(Y, K, W, WTW, n_all_par,
                                    J, M, O,
                                    lambdasq_beta_update, tausq_beta_update,
                                    lambdasq_gamma_update, tausq_gamma_update,
                                    lambdasq_delta_update, tausq_delta_update,
                                    sigmasq_varphi, Sigma_update)


set.seed(76789)
update_theta <- update_theta_BVS_HHS_MVN_cov(Y, K, W, WTW, n_all_par,
                                             J, M, O,
                                             lambdasq_beta_update, tausq_beta_update,
                                             lambdasq_gamma_update, tausq_gamma_update,
                                             lambdasq_delta_update, tausq_delta_update,
                                             sigmasq_varphi, Sigma_update)
all.equal(update_theta_modify, update_theta)
max(abs(update_theta_modify - update_theta))

all.equal(Q, Q_modify)
max(abs(Q-Q_modify))


############
# Update tausq_beta (for all j = 1, 2, ..., J)
update_tausq_beta_BVS_HHS_MVN_cov2 <- function(
    J, K,
    beta_update,
    lambdasq_beta_update,
    xi_beta_update)
{
  shape_tausq_beta <- (K + 1) / 2

  scale_tausq_beta <- (1 / xi_beta_update) +
    rowSums(beta_update^2 / (2 * lambdasq_beta_update))

  rinvgamma(J,
            shape = shape_tausq_beta,
            scale = scale_tausq_beta)
}





update_tausq_beta_BVS_HHS_MVN_cov <- function(J, K, beta_update, lambdasq_beta_update, xi_beta_update)
{
  tausq_beta_update_s <- rep(NA, J)

  shape_tausq_beta <- (K + 1) / 2
  for(j in 1:J)
  {
    sum_term_tausq_beta <- 0
    for (k in 1:K)
    {
      sum_term_tausq_beta <- sum_term_tausq_beta + (beta_update[j, k]^2 / (2 *
                                                                             lambdasq_beta_update[j, k]))
    }

    scale_tausq_beta <- ((1/xi_beta_update) + sum_term_tausq_beta)
    tausq_beta_update_s[j] <- rinvgamma(1, shape = shape_tausq_beta, scale = scale_tausq_beta) # IG distribution
  }
  return(tausq_beta_update_s)
}


xi_beta_update <- 1
beta_update <- update_theta_modify[1:J, ]

set.seed(456)
t1 <- update_tausq_beta_BVS_HHS_MVN_cov(J, K, beta_update,
             lambdasq_beta_update, xi_beta_update)
set.seed(456)
t2 <- update_tausq_beta_BVS_HHS_MVN_cov2(J, K, beta_update,
                      lambdasq_beta_update, xi_beta_update)


all.equal(t1, t2)
max(abs(t1-t2))


##########################
lambdasq_beta_update <- matrix(rinvgamma(J*K, 2, 3), ncol = K)
lambdasq_gamma_update <- matrix(rinvgamma(M*K, 2, 3), ncol = K)
lambdasq_delta_update <- matrix(rinvgamma(J*M*K, 2, 3), ncol = K)
tausq_beta_update <- rgamma(J, 2, 3)
tausq_gamma_update <- rgamma(M, 2, 3)
tausq_delta_update <- rgamma(J*M, 2, 3)
c <- 2.5

lambdasq_beta_update_tilde <- ((c^2 * lambdasq_beta_update) /
                                 (c^2 + matrix(tausq_beta_update, J, K) * lambdasq_beta_update))
sigmasq_beta_update <- lambdasq_beta_update_tilde * matrix(tausq_beta_update, J, K)

lambdasq_gamma_update_tilde <- ((c^2 * lambdasq_gamma_update) /
                                  (c^2 + matrix(tausq_gamma_update, M, K) * lambdasq_gamma_update))
sigmasq_gamma_update <- lambdasq_gamma_update_tilde * matrix(tausq_gamma_update, M, K)

sigmasq_delta_update <- matrix(NA, J * M, K)
for (j in 1:J) {
  for (m in 1:M) {
    jm <- (j - 1) * M + m   # row index for (j,m)
    for (k in 1:K) {
      lambdasq_delta_update_tilde <- ((c^2 * lambdasq_delta_update[jm, k]) /
                                        (c^2 + tausq_delta_update[jm] * lambdasq_delta_update[jm, k]))
      sigmasq_delta_update[jm, k] <-
        (lambdasq_delta_update_tilde *
           lambdasq_beta_update_tilde[j, k] *
           lambdasq_gamma_update_tilde[m, k] *
           tausq_delta_update[jm])
    }
  }
}



beta_rep  <- lambdasq_beta_update_tilde[rep(1:J, each = M), ]
gamma_rep <- lambdasq_gamma_update_tilde[rep(1:M, times = J), ]
tausq_rep <- matrix(tausq_delta_update, nrow = J * M, ncol = K)
lambdasq_delta_update_tilde <- ((c^2 * lambdasq_delta_update) /
                                  (c^2 + tausq_rep * lambdasq_delta_update))

sigmasq_delta_update2 <-
  lambdasq_delta_update_tilde *
  beta_rep *
  gamma_rep *
  tausq_rep

all.equal(sigmasq_delta_update, sigmasq_delta_update2)
max(abs(sigmasq_delta_update - sigmasq_delta_update2))

################################################

# Update theta
update_theta_BVS_HRHS_MVN_cov <- function(Y, K, W, WTW, n_all_par,
                                          J, M, O,
                                          c,
                                          lambdasq_beta_update, tausq_beta_update,
                                          lambdasq_gamma_update, tausq_gamma_update,
                                          lambdasq_delta_update, tausq_delta_update,
                                          sigmasq_varphi, Sigma_update
)
{
  lambdasq_beta_update_tilde <- ((c^2 * lambdasq_beta_update) /
                                   (c^2 + matrix(tausq_beta_update, J, K) * lambdasq_beta_update))
  sigmasq_beta_update <- lambdasq_beta_update_tilde * matrix(tausq_beta_update, J, K)

  lambdasq_gamma_update_tilde <- ((c^2 * lambdasq_gamma_update) /
                                    (c^2 + matrix(tausq_gamma_update, M, K) * lambdasq_gamma_update))
  sigmasq_gamma_update <- lambdasq_gamma_update_tilde * matrix(tausq_gamma_update, M, K)


  sigmasq_delta_update <- matrix(NA, J * M, K)
  for (j in 1:J) {
    for (m in 1:M) {
      jm <- (j - 1) * M + m   # row index for (j,m)
      for (k in 1:K) {
        lambdasq_delta_update_tilde <- ((c^2 * lambdasq_delta_update[jm, k]) /
                                          (c^2 + tausq_delta_update[jm] * lambdasq_delta_update[jm, k]))
        sigmasq_delta_update[jm, k] <-
          (lambdasq_delta_update_tilde *
             lambdasq_beta_update_tilde[j, k] *
             lambdasq_gamma_update_tilde[m, k] *
             tausq_delta_update[jm])
      }
    }
  }


  # Build list of prior covariance matrices V_theta_k for each k
  V_theta_list <- vector("list", K)
  for (k in 1:K)
  {
    # Gibbs update of theta
    var_alpha_0k <- 100
    mat_var_alpha_0k <- matrix(var_alpha_0k, nrow = 1, ncol = 1) #make it matrix
    var_beta_k <- diag(sigmasq_beta_update[, k]) # dim J*J
    var_gamma_k <- diag(sigmasq_gamma_update[, k]) # dim M*M
    var_delta_k <- diag(sigmasq_delta_update[, k]) # dim JM*JM

    vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
    mat_var_varphi <- diag(vec_sigmasq_varphi)

    V_theta_list[[k]] <- as.matrix(bdiag(mat_var_alpha_0k, var_beta_k, var_gamma_k,
                                         var_delta_k, mat_var_varphi))
  }
  # Build block-diagonal prior precision: blockdiag(V_theta1^{-1},...,V_thetaK^{-1})
  prior_prec_list <- lapply(V_theta_list, function(Vk) {
    chol2inv(chol(Vk))
  })
  prior_prec <- as.matrix(bdiag(prior_prec_list))   # pK x pK

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

#############
# Update theta
update_theta_BVS_HRHS_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par,
                                                 J, M, O,
                                                 c,
                                                 lambdasq_beta_update, tausq_beta_update,
                                                 lambdasq_gamma_update, tausq_gamma_update,
                                                 lambdasq_delta_update, tausq_delta_update,
                                                 sigmasq_varphi, Sigma_update
)
{
  lambdasq_beta_update_tilde <- ((c^2 * lambdasq_beta_update) /
                                   (c^2 + matrix(tausq_beta_update, J, K) * lambdasq_beta_update))
  sigmasq_beta_update <- lambdasq_beta_update_tilde * matrix(tausq_beta_update, J, K)

  lambdasq_gamma_update_tilde <- ((c^2 * lambdasq_gamma_update) /
                                    (c^2 + matrix(tausq_gamma_update, M, K) * lambdasq_gamma_update))
  sigmasq_gamma_update <- lambdasq_gamma_update_tilde * matrix(tausq_gamma_update, M, K)

  beta_rep  <- lambdasq_beta_update_tilde[rep(1:J, each = M), ]
  gamma_rep <- lambdasq_gamma_update_tilde[rep(1:M, times = J), ]
  tausq_rep <- matrix(tausq_delta_update, nrow = J * M, ncol = K)
  lambdasq_delta_update_tilde <- ((c^2 * lambdasq_delta_update) /
                                    (c^2 + tausq_rep * lambdasq_delta_update))
  sigmasq_delta_update <- (lambdasq_delta_update_tilde * beta_rep *
                             gamma_rep * tausq_rep)


  # Build list of prior covariance matrices V_theta_k for each k
  V_theta_list <- vector("list", K)
  for (k in 1:K)
  {
    # Gibbs update of theta
    var_alpha_0k <- 100
    var_delta_k <- diag(sigmasq_delta_update[, k]) # dim JM*JM

    vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
    mat_var_varphi <- diag(vec_sigmasq_varphi)

    V_theta_k_inv_list[[k]] <- as.matrix(
      bdiag(
        1/var_alpha_0k,
        diag(1/sigmasq_beta_update[, k]),
        diag(1/sigmasq_gamma_update[, k]),
        diag(1/sigmasq_delta_update[, k]),
        diag(1/diag(mat_var_varphi))
      )
    )
  }

  prior_prec <- as.matrix(bdiag(V_theta_k_inv_list))   # pK x pK

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

##############

lambdasq_beta_update <- matrix(rinvgamma(J*K, 2, 3), ncol = K)
lambdasq_gamma_update <- matrix(rinvgamma(M*K, 2, 3), ncol = K)
lambdasq_delta_update <- matrix(rinvgamma(J*M*K, 2, 3), ncol = K)
tausq_beta_update <- rgamma(J, 2, 3)
tausq_gamma_update <- rgamma(M, 2, 3)
tausq_delta_update <- rgamma(J*M, 2, 3)
sigmasq_varphi <- 10

nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_update <- rinvwishart(nu_0, Psi_0)

set.seed(232)
theta_update <- update_theta_BVS_HRHS_MVN_cov(Y, K, W, WTW, n_all_par,
                                                          J, M, O,
                                                          c,
                                                          lambdasq_beta_update, tausq_beta_update,
                                                          lambdasq_gamma_update, tausq_gamma_update,
                                                          lambdasq_delta_update, tausq_delta_update,
                                                          sigmasq_varphi, Sigma_update)


set.seed(232)
theta_update_modify <- update_theta_BVS_HRHS_MVN_cov_modify(Y, K, W, WTW, n_all_par,
                                                                        J, M, O,
                                                                        c,
                                                                        lambdasq_beta_update, tausq_beta_update,
                                                                        lambdasq_gamma_update, tausq_gamma_update,
                                                                        lambdasq_delta_update, tausq_delta_update,
                                                                        sigmasq_varphi, Sigma_update)
all.equal(theta_update, theta_update_modify)
max(abs(theta_update-theta_update_modify))

#################









