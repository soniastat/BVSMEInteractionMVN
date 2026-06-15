
# Functions for simulating data and fitting BVS-NHHS-MVN-cov model

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
# Fit NHHS model
########################################################

# Update theta
update_theta_BVS_NHHS_MVN_cov <- function(Y, K, W, n_all_par,
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
      for (k in 1:K) {
        sigmasq_delta_update[jm, k] <- (lambdasq_delta_update[jm, k] *
                                          tausq_delta_update[jm])
      }
    }
  }


  # Build list of prior covariance matrices V_theta_k for each k
  V_theta_list <- vector("list", K)
  for (k in 1:K)
  {
    # Gibbs update of theta
    var_alpha_0k <- 100 # dim 1*1
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
    chol2inv(chol(Vk))  # ensure Vk is PD before inversion
  })
  prior_prec <- as.matrix(bdiag(prior_prec_list))   # pK x pK

  # Kronecker product term: Σ^{-1} ⊗ (W^T W)
  Sigma_inv <- as.matrix(chol2inv(chol(Sigma_update)))  # K x K
  kron_term <- kronecker(Sigma_inv, t(W) %*% W) # pK x pK

  # Full conditional covariance
  Sigma_theta <- chol2inv(chol(kron_term + prior_prec)) # pK x pK

  # Full conditional mean
  vecY <- as.vector(Y)   # vec(Y), column-stacked
  mean_theta_vec <- Sigma_theta %*% (kronecker(Sigma_inv, t(W)) %*% vecY)

  # Draw from MVN
  theta_vec <- mvrnorm(1, mu = mean_theta_vec, Sigma = Sigma_theta)

  # Reshape back to n_all_par x K
  theta_update_s <- matrix(theta_vec, nrow = n_all_par, ncol = K)

  return(theta_update_s)
}

#######################################################################################
# Update Sigma
update_Sigma_BVS_NHHS_MVN_cov <- function(Y, n, W, theta_update, Psi_0, nu_0)
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

# Update lambdasq_beta (for all j = 1, 2, ..., J, and k = 1, 2, ..., K)
update_lambdasq_beta_BVS_NHHS_MVN_cov <- function(J, K, beta_update, psi_beta_update, tausq_beta_update)
{
  lambdasq_beta_update_s <- matrix(NA, nrow = J, ncol = K)

  for(j in 1:J)
  {
    for(k in 1:K)
    {
      scale_lambdasq_beta <- ((1/psi_beta_update[j, k]) +
                                ((beta_update[j, k])^2 / (2 * tausq_beta_update[j])))
      lambdasq_beta_update_s[j, k] <- rinvgamma(1, shape = 1, scale = scale_lambdasq_beta) # IG distribution
    }
  }
  return(lambdasq_beta_update_s)
}



#######################################################################################

# Update lambdasq_gamma (for all m = 1, 2, ..., M, and k = 1, 2, ..., K)
update_lambdasq_gamma_BVS_NHHS_MVN_cov <- function(M, K, gamma_update, psi_gamma_update, tausq_gamma_update)
{
  lambdasq_gamma_update_s <- matrix(NA, nrow = M, ncol = K)

  for(m in 1:M)
  {
    for(k in 1:K)
    {
      scale_lambdasq_gamma <- ((1/psi_gamma_update[m, k]) +
                                 ((gamma_update[m, k])^2 / (2 * tausq_gamma_update[m])))
      lambdasq_gamma_update_s[m, k] <- rinvgamma(1, shape = 1, scale = scale_lambdasq_gamma) # IG distribution
    }
  }
  return(lambdasq_gamma_update_s)
}



#####################################################################################

# Update lambdasq_delta (for all j and m, and k = 1, 2, ..., K)
update_lambdasq_delta_BVS_NHHS_MVN_cov <- function(J, M, K, delta_update,
                                                   psi_delta_update, tausq_delta_update)
{
  lambdasq_delta_update_s <- matrix(NA, nrow = J*M, ncol = K)

  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction
      for(k in 1:K)
      {
        scale_lambdasq_delta <- ((1/psi_delta_update[jm, k]) + ((delta_update[jm, k])^2 / (2 *
                                                                                             tausq_delta_update[jm])))

        lambdasq_delta_update_s[jm, k] <- rinvgamma(1, shape = 1, scale = scale_lambdasq_delta) # IG distribution
      }
    }
  }
  return(lambdasq_delta_update_s)
}


#######################################################################################

# Update tausq_beta (for all j = 1, 2, ..., J)
update_tausq_beta_BVS_NHHS_MVN_cov <- function(J, K, beta_update, lambdasq_beta_update, xi_beta_update)
{
  tausq_beta_update_s <- rep(NA, J)

  shape_tausq_beta <- (K+1) / 2
  for(j in 1:J)
  {
    sum_term_tausq_beta <- 0
    for (k in 1:K) {
      sum_term_tausq_beta <- sum_term_tausq_beta + (beta_update[j, k]^2 /
                                                      (2 * lambdasq_beta_update[j, k]))
    }

    scale_tausq_beta <- ((1/xi_beta_update) + sum_term_tausq_beta)
    tausq_beta_update_s[j] <- rinvgamma(1, shape = shape_tausq_beta, scale = scale_tausq_beta) # IG distribution
  }
  return(tausq_beta_update_s)
}


#######################################################################################

# Update tausq_gamma (for all m = 1, 2, ..., M)
update_tausq_gamma_BVS_NHHS_MVN_cov <- function(M, K, gamma_update, lambdasq_gamma_update, xi_gamma_update)
{
  tausq_gamma_update_s <- rep(NA, M)

  shape_tausq_gamma <- (K+1) / 2
  for(m in 1:M)
  {
    sum_term_tausq_gamma <- 0
    for (k in 1:K) {
      sum_term_tausq_gamma <- sum_term_tausq_gamma + (gamma_update[m, k]^2 / (2 * lambdasq_gamma_update[m, k]))
    }

    scale_tausq_gamma <- ((1/xi_gamma_update) + sum_term_tausq_gamma)
    tausq_gamma_update_s[m] <- rinvgamma(1, shape = shape_tausq_gamma, scale = scale_tausq_gamma) # IG distribution
  }
  return(tausq_gamma_update_s)
}


#######################################################################################

# Update tausq_delta (for all j, and m )
update_tausq_delta_BVS_NHHS_MVN_cov <- function(J, M, K, delta_update, lambdasq_delta_update,
                                                xi_delta_update)
{
  tausq_delta_update_s <- rep(NA, J * M)

  shape_tausq_delta <- (K+1) / 2
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction

      sum_term_tausq_delta <- 0
      for (k in 1:K)
      {
        sum_term_tausq_delta <- sum_term_tausq_delta + ((delta_update[jm, k])^2 / (2 *
                                                                                     lambdasq_delta_update[jm, k]))
      }

      scale_tausq_delta <- ((1/xi_delta_update) + sum_term_tausq_delta)
      tausq_delta_update_s[jm] <- rinvgamma(1, shape = shape_tausq_delta, scale = scale_tausq_delta) # IG distribution
    }
  }
  return(tausq_delta_update_s)
}


#######################################################################################

# Update psi_beta (for all j = 1, 2, ..., J, and k = 1, 2, ..., K)
update_psi_beta_BVS_NHHS_MVN_cov <- function(J, K, lambdasq_beta_update)
{
  psi_beta_update_s <- matrix(NA, nrow = J, ncol = K)

  for(j in 1:J)
  {
    for(k in 1:K)
    {
      scale_psi_beta <- (1 + (1/lambdasq_beta_update[j, k]))
      psi_beta_update_s[j, k] <- rinvgamma(1, shape = 1, scale = scale_psi_beta) # IG distribution
    }
  }
  return(psi_beta_update_s)
}



#######################################################################################

# Update psi_gamma (for all m = 1, 2, ..., M, and k = 1, 2, ..., K)
update_psi_gamma_BVS_NHHS_MVN_cov <- function(M, K, lambdasq_gamma_update)
{
  psi_gamma_update_s <- matrix(NA, nrow = M, ncol = K)

  for(m in 1:M)
  {
    for(k in 1:K)
    {
      scale_psi_gamma <- (1 + (1/lambdasq_gamma_update[m, k]))
      psi_gamma_update_s[m, k] <- rinvgamma(1, shape = 1, scale = scale_psi_gamma) # IG distribution
    }
  }
  return(psi_gamma_update_s)
}



#######################################################################################

# Update psi_delta (for all j , m, and k = 1, 2, ..., K)
update_psi_delta_BVS_NHHS_MVN_cov <- function(J, M, K, lambdasq_delta_update)
{
  psi_delta_update_s <- matrix(NA, nrow = J*M, ncol = K)

  for(jm in 1:(J*M))
  {
    for(k in 1:K)
    {
      scale_psi_delta <- (1 + (1/lambdasq_delta_update[jm, k]))
      psi_delta_update_s[jm, k] <- rinvgamma(1, shape = 1, scale = scale_psi_delta) # IG distribution
    }
  }
  return(psi_delta_update_s)
}




#######################################################################################

# Update xi_beta
update_xi_beta_BVS_NHHS_MVN_cov <- function(J, tausq_beta_update)
{
  # Use more informative prior (xi_beta ~ IG(15, 3))
  shape_xi_beta <- ((J/2) + 15)
  scale_xi_beta <- (3 + sum(1 / tausq_beta_update))

  xi_beta_update_s <- rinvgamma(1, shape = shape_xi_beta, scale = scale_xi_beta) # IG distribution
  return(xi_beta_update_s)
}



#######################################################################################

# Update xi_gamma
update_xi_gamma_BVS_NHHS_MVN_cov <- function(M, tausq_gamma_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_gamma <- ((M/2) + 15)
  scale_xi_gamma <- (3 + sum(1 / tausq_gamma_update))

  xi_gamma_update_s <- rinvgamma(1, shape = shape_xi_gamma, scale = scale_xi_gamma) # IG distribution
  return(xi_gamma_update_s)
}





#######################################################################################

# Update xi_delta (for j and m)
update_xi_delta_BVS_NHHS_MVN_cov <- function(J, M, tausq_delta_update)
{
  # Use more informative prior (xi_delta ~ IG(15, 3))
  shape_xi_delta <- ((J*M/2) + 15)
  scale_xi_delta <- (3 + sum(1 / tausq_delta_update))

  xi_delta_update_s <- rinvgamma(1, shape = shape_xi_delta, scale = scale_xi_delta) # IG distribution
  return(xi_delta_update_s)
}



##################################################################################

# Update parameters for BVS-NHHS-MVN-cov model
fit_BVS_NHHS_MVN_cov <- function(niter = 20, burn_in = 2, thin = 1,
                                 n, K, Y, W, n_all_par, J, M, O,
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
                                 Sigma_init,
                                 nu_0, Psi_0,
                                 sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))
  lambdasq_beta_update <- array(NA, dim = c(niter, J, K))
  tausq_beta_update <- matrix(NA, nrow = niter, ncol = J)
  lambdasq_gamma_update <- array(NA, dim = c(niter, M, K))
  tausq_gamma_update <- matrix(NA, nrow = niter, ncol = M)
  lambdasq_delta_update <- array(NA, dim = c(niter, J*M, K))
  tausq_delta_update <- matrix(NA, nrow = niter, ncol = J*M)
  psi_beta_update <- array(NA, dim = c(niter, J, K))
  psi_gamma_update <- array(NA, dim = c(niter, M, K))
  psi_delta_update <- array(NA, dim = c(niter, J*M, K))
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
  lambdasq_beta_update[1, , ] <- lambdasq_beta_init
  tausq_beta_update[1, ] <- tausq_beta_init
  lambdasq_gamma_update[1, , ] <- lambdasq_gamma_init
  tausq_gamma_update[1, ] <- tausq_gamma_init
  lambdasq_delta_update[1, , ] <- lambdasq_delta_init
  tausq_delta_update[1, ] <- tausq_delta_init
  psi_beta_update[1, , ] <- psi_beta_init
  psi_gamma_update[1, , ] <- psi_gamma_init
  psi_delta_update[1, , ] <- psi_delta_init
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
    theta_update_s <- update_theta_BVS_NHHS_MVN_cov(Y, K, W, n_all_par,
                                                    J, M, O,
                                                    lambdasq_beta_update = lambdasq_beta_update[(s-1), , ],
                                                    tausq_beta_update = tausq_beta_update[(s-1), ],
                                                    lambdasq_gamma_update = lambdasq_gamma_update[(s-1), , ],
                                                    tausq_gamma_update = tausq_gamma_update[(s-1), ],
                                                    lambdasq_delta_update = lambdasq_delta_update[(s-1), , ],
                                                    tausq_delta_update = tausq_delta_update[(s-1), ],
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
    Sigma_update_s <- update_Sigma_BVS_NHHS_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s



    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_BVS_NHHS_MVN_cov(J, K,
                                                                    beta_update = beta_update[s, , ],
                                                                    psi_beta_update = psi_beta_update[(s-1), , ],
                                                                    tausq_beta_update = tausq_beta_update[(s-1), ])
    lambdasq_beta_update[s, , ] <- lambdasq_beta_update_s



    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_BVS_NHHS_MVN_cov(M, K,
                                                                      gamma_update = gamma_update[s, , ],
                                                                      psi_gamma_update = psi_gamma_update[(s-1), , ],
                                                                      tausq_gamma_update = tausq_gamma_update[(s-1), ])
    lambdasq_gamma_update[s, , ] <- lambdasq_gamma_update_s



    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_BVS_NHHS_MVN_cov(J, M, K,
                                                                      delta_update = delta_update[s, , ],
                                                                      psi_delta_update = psi_delta_update[(s-1), , ],
                                                                      tausq_delta_update = tausq_delta_update[(s-1), ])
    lambdasq_delta_update[s, , ] <- lambdasq_delta_update_s



    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_BVS_NHHS_MVN_cov(J, K,
                                                              beta_update = beta_update[s, , ],
                                                              lambdasq_beta_update = lambdasq_beta_update[s, , ],
                                                              xi_beta_update = xi_beta_update[(s-1)])
    tausq_beta_update[s, ] <- tausq_beta_update_s



    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_BVS_NHHS_MVN_cov(M, K,
                                                                gamma_update = gamma_update[s, , ],
                                                                lambdasq_gamma_update = lambdasq_gamma_update[s, , ],
                                                                xi_gamma_update = xi_gamma_update[(s-1)])
    tausq_gamma_update[s, ] <- tausq_gamma_update_s



    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_BVS_NHHS_MVN_cov(J, M, K,
                                                                delta_update = delta_update[s, , ],
                                                                lambdasq_delta_update = lambdasq_delta_update[s, , ],
                                                                xi_delta_update = xi_delta_update[(s-1)])
    tausq_delta_update[s, ] <- tausq_delta_update_s



    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_BVS_NHHS_MVN_cov(J, K,
                                                          lambdasq_beta_update = lambdasq_beta_update[s, , ])
    psi_beta_update[s, , ] <- psi_beta_update_s



    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_BVS_NHHS_MVN_cov(M, K,
                                                            lambdasq_gamma_update = lambdasq_gamma_update[s, , ])
    psi_gamma_update[s, , ] <- psi_gamma_update_s



    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_BVS_NHHS_MVN_cov(J, M, K,
                                                            lambdasq_delta_update = lambdasq_delta_update[s, , ])
    psi_delta_update[s, , ] <- psi_delta_update_s



    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_BVS_NHHS_MVN_cov(J, tausq_beta_update = tausq_beta_update[s, ])
    xi_beta_update[s] <- xi_beta_update_s



    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_BVS_NHHS_MVN_cov(M, tausq_gamma_update = tausq_gamma_update[s, ])
    xi_gamma_update[s] <- xi_gamma_update_s



    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_BVS_NHHS_MVN_cov(J, M, tausq_delta_update = tausq_delta_update[s, ])
    xi_delta_update[s] <- xi_delta_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ],
    lambdasq_beta_update = lambdasq_beta_update[indices_to_save, , ],
    tausq_beta_update = tausq_beta_update[indices_to_save, ],
    lambdasq_gamma_update = lambdasq_gamma_update[indices_to_save, , ],
    tausq_gamma_update = tausq_gamma_update[indices_to_save, ],
    lambdasq_delta_update = lambdasq_delta_update[indices_to_save, , ],
    tausq_delta_update = tausq_delta_update[indices_to_save, ],
    psi_beta_update = psi_beta_update[indices_to_save, , ],
    psi_gamma_update = psi_gamma_update[indices_to_save, , ],
    psi_delta_update = psi_delta_update[indices_to_save, , ],
    xi_beta_update = xi_beta_update[indices_to_save],
    xi_gamma_update = xi_gamma_update[indices_to_save],
    xi_delta_update = xi_delta_update[indices_to_save]
  )
  return(results)
}

################################################################

