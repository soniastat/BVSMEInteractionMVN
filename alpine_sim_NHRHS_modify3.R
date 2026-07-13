
# Functions for simulating data and fitting NHRHS-MVN-cov model
# considered more variables to be active compared to the modify2 files
###################################################################


# Simulate data
Sim_data_BVS_real <- function(K = 4, n = 1301, J = 44, M = 17, O = 6,
                              x, z, d, Sigma)
{
  # modifiers
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- seq(0.25, 1.25, length.out = K)
  beta_true[2, ] <- seq(-1.25, -.25, length.out = K)
  beta_true[3, ] <- seq(0.5, 1.5, length.out = K)
  beta_true[4, ] <- seq(-1.25, 1.25, length.out = K)
  beta_true[5, ] <- seq(0.25, 2.25, length.out = K)
  beta_true[6, ] <- seq(1.5, 3.5, length.out = K)

  # environmental exposures
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different gamma for different k
  gamma_true[1, ] <- seq(0.5, 1.5, length.out = K)
  gamma_true[2, ] <- seq(-1.5, -0.5, length.out = K)
  gamma_true[3, ] <- seq(-1.5, 2, length.out = K)
  gamma_true[4, ] <- seq(1.25, 0.5, length.out = K)
  gamma_true[5, ] <- seq(1.5, 2, length.out = K)

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
  varphi_true[3, ] <- seq(-1, 2.25, length.out = K)

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
# all variables are simulated, mimicking real data by number of variable
# Sigma=error variance is same as the real data

Sim_data_BVS_real_sim_var <- function(K = 4, n = 1301, J = 44, M = 17, O = 6, Sigma)
{
  # genetic variables
  x <- matrix(rnorm(n * J), nrow = n)
  beta_true <- matrix(0, nrow = J, ncol = K) # Different beta for different k
  beta_true[1, ] <- seq(0.25, 1.25, length.out = K)
  beta_true[2, ] <- seq(-1.25, -.25, length.out = K)
  beta_true[3, ] <- seq(0.5, 1.5, length.out = K)
  beta_true[4, ] <- seq(-1.25, 1.25, length.out = K)
  beta_true[5, ] <- seq(0.25, 2.25, length.out = K)
  beta_true[6, ] <- seq(1.5, 3.5, length.out = K)

  # environmental variables
  z <- matrix(rnorm(n * M), nrow = n)
  gamma_true <- matrix(0, nrow = M, ncol = K) # Different beta for different k
  gamma_true[1, ] <- seq(0.5, 1.5, length.out = K)
  gamma_true[2, ] <- seq(-1.5, -0.5, length.out = K)
  gamma_true[3, ] <- seq(-1.5, 2, length.out = K)
  gamma_true[4, ] <- seq(1.25, 0.5, length.out = K)
  gamma_true[5, ] <- seq(1.5, 2, length.out = K)


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
  active_beta  <- which(rowSums(abs(beta_true)) > 0)
  active_gamma <- which(rowSums(abs(gamma_true)) > 0)
  for (j in active_beta) {
    for (m in active_gamma) {
      idx <- (j - 1) * M + m
      delta_true[idx, ] <- rnorm(K)
    }
  }

  # covariates
  d <- matrix(rnorm(n * O), nrow = n)
  varphi_true <- matrix(0, nrow = O, ncol = K) # Different beta for different k
  varphi_true[1, ] <- seq(0.25, 1.25, length.out = K)
  varphi_true[2, ] <- seq(-1.25, -.25, length.out = K)
  varphi_true[3, ] <- seq(-1, 2.25, length.out = K)

  alpha_true <- runif(K)

  theta_true <- rbind(matrix(alpha_true, nrow = 1),
                      beta_true,
                      gamma_true,
                      delta_true,
                      varphi_true)
  W <- cbind(1, x, z, u, d)
  W_theta_true <- W %*% theta_true # n*K dimension

  # Cholesky of Sigma
  R_Sigma <- chol(Sigma)

  # E ~ iid N(0,1)
  E_mat <- matrix(rnorm(n * K), nrow = n, ncol = K)

  # Y = M + E * chol(Sigma)
  Y <- W_theta_true + E_mat %*% R_Sigma

  return(list(
    x = x, z = z, u = u, d = d, W = W,
    alpha_true = alpha_true,
    beta_true = beta_true,
    gamma_true = gamma_true,
    delta_true = delta_true,
    varphi_true = varphi_true,
    W_theta_true = W_theta_true,
    Y = Y, Sigma = Sigma,
    theta_true = theta_true
  ))
}



#########################################################
# Fit NHRHS model
########################################################

# Update theta
update_theta_NHRHS_MVN_cov_modify <- function(Y, K, W, WTW,
                                           n_all_par,
                                           J, M, O,
                                           c,
                                           sigmasq_alpha = 100,
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


  tausq_rep <- matrix(tausq_delta_update, nrow = J * M, ncol = K)
  lambdasq_delta_update_tilde <- ((c^2 * lambdasq_delta_update) /
                                    (c^2 + tausq_rep * lambdasq_delta_update))
  sigmasq_delta_update <- (lambdasq_delta_update_tilde * tausq_rep)


  # Build list of prior covariance matrices V_theta_k for each k
  V_theta_k_inv_list <- vector("list", K)
  for (k in 1:K)
  {
    # Gibbs update of theta
    vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)

    V_theta_k_inv_list[[k]] <- as.matrix(
      bdiag(
        1/sigmasq_alpha,
        diag(1/sigmasq_beta_update[, k]),
        diag(1/sigmasq_gamma_update[, k]),
        diag(1/sigmasq_delta_update[, k]),
        diag(1/vec_sigmasq_varphi)
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

#######################################################################################
# Update Sigma
update_Sigma_NHRHS_MVN_cov <- function(Y, n, W, theta_update, Psi_0, nu_0)
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
update_lambdasq_beta_NHRHS_MVN_cov <- function(c, J, K, beta_update, psi_beta_update, tausq_beta_update,
                                                   lambdasq_beta_update_s_1, omegasq_beta_lambdasq)
{
  # MH update for lambdasq_beta
  lambdasq_beta_update_s <- lambdasq_beta_update_s_1
  accept_lambdasq_beta_update_s <- matrix(0, nrow = J, ncol = K)

  for(j in 1:J)
  {
    for(k in 1:K)
    {
      lambdasq_beta_prop <- rfnorm(1, lambdasq_beta_update_s[j, k], sqrt(omegasq_beta_lambdasq)) # Folded normal

      # value for j,k (single not a matrix)
      lambdasq_beta_prop_tilde <- ((c^2 * lambdasq_beta_prop) /
                                     (c^2 + tausq_beta_update[j] * lambdasq_beta_prop))
      lambdasq_beta_update_s_1_tilde <- ((c^2 * lambdasq_beta_update_s[j, k]) /
                                           (c^2 + tausq_beta_update[j] * lambdasq_beta_update_s[j, k]))

      log_accept_ratio_lambdasq_beta <- (
        -0.5 * (log(lambdasq_beta_prop_tilde) - log(lambdasq_beta_update_s_1_tilde))
        - (3/2) * (log(lambdasq_beta_prop) - log(lambdasq_beta_update_s[j, k]))
        - ((beta_update[j, k])^2 / (2 * tausq_beta_update[j])) *
          ((1/lambdasq_beta_prop_tilde) - (1/lambdasq_beta_update_s_1_tilde))
        - (1/psi_beta_update[j, k]) * ((1/lambdasq_beta_prop) - (1/lambdasq_beta_update_s[j, k]))
      )

      # Accept/reject the proposed lambdasq_beta
      if (log(runif(1)) < log_accept_ratio_lambdasq_beta)
      {
        lambdasq_beta_update_s[j, k] <- lambdasq_beta_prop
        accept_lambdasq_beta_update_s[j, k] <- 1
      }
    }
  }
  return(list(lambdasq_beta_update = lambdasq_beta_update_s,
              accept_lambdasq_beta_update = accept_lambdasq_beta_update_s))
}



#######################################################################################

# Update lambdasq_gamma (for all m = 1, 2, ..., M, and k = 1, 2, ..., K)
update_lambdasq_gamma_NHRHS_MVN_cov <- function(c, M, K, gamma_update, psi_gamma_update,
                                                    tausq_gamma_update,
                                                    lambdasq_gamma_update_s_1, omegasq_gamma_lambdasq)
{
  # MH update for lambdasq_gamma
  lambdasq_gamma_update_s <- lambdasq_gamma_update_s_1
  accept_lambdasq_gamma_update_s <- matrix(0, nrow = M, ncol = K)

  for(m in 1:M)
  {
    for(k in 1:K)
    {
      lambdasq_gamma_prop <- rfnorm(1, lambdasq_gamma_update_s[m, k], sqrt(omegasq_gamma_lambdasq)) # Folded normal

      # value for m,k (single not a matrix)
      lambdasq_gamma_prop_tilde <- ((c^2 * lambdasq_gamma_prop) /
                                      (c^2 + tausq_gamma_update[m] * lambdasq_gamma_prop))
      lambdasq_gamma_update_s_1_tilde <- ((c^2 * lambdasq_gamma_update_s[m, k]) /
                                            (c^2 + tausq_gamma_update[m] * lambdasq_gamma_update_s[m, k]))

      log_accept_ratio_lambdasq_gamma <- (
        -0.5 * (log(lambdasq_gamma_prop_tilde) - log(lambdasq_gamma_update_s_1_tilde))
        - (3/2) * (log(lambdasq_gamma_prop) - log(lambdasq_gamma_update_s[m, k]))
        - ((gamma_update[m, k])^2 / (2 * tausq_gamma_update[m])) *
          ((1/lambdasq_gamma_prop_tilde) - (1/lambdasq_gamma_update_s_1_tilde))
        - (1/psi_gamma_update[m, k]) * ((1/lambdasq_gamma_prop) - (1/lambdasq_gamma_update_s[m, k]))
      )

      # Accept/reject the proposed lambdasq_gamma
      if (log(runif(1)) < log_accept_ratio_lambdasq_gamma)
      {
        lambdasq_gamma_update_s[m, k] <- lambdasq_gamma_prop
        accept_lambdasq_gamma_update_s[m, k] <- 1
      }
    }
  }
  return(list(lambdasq_gamma_update = lambdasq_gamma_update_s,
              accept_lambdasq_gamma_update = accept_lambdasq_gamma_update_s))
}



#####################################################################################

# Update lambdasq_delta (for all j and m, and k = 1, 2, ..., K)
update_lambdasq_delta_NHRHS_MVN_cov <- function(c, J, M, K, delta_update, psi_delta_update,
                                                    tausq_delta_update,
                                                    lambdasq_delta_update_s_1, omegasq_delta_lambdasq)
{
  # MH update for lambdasq_delta
  lambdasq_delta_update_s <- lambdasq_delta_update_s_1
  accept_lambdasq_delta_update_s <- matrix(0, nrow = J*M, ncol = K)

  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction
      for(k in 1:K)
      {
        lambdasq_delta_prop <- rfnorm(1, lambdasq_delta_update_s[jm, k], sqrt(omegasq_delta_lambdasq)) # Folded normal

        # value for j,k (single not a matrix)
        lambdasq_delta_prop_tilde <- ((c^2 * lambdasq_delta_prop) /
                                        (c^2 + tausq_delta_update[jm] * lambdasq_delta_prop))
        lambdasq_delta_update_s_1_tilde <- ((c^2 * lambdasq_delta_update_s[jm, k]) /
                                              (c^2 + tausq_delta_update[jm] * lambdasq_delta_update_s[jm, k]))

        log_accept_ratio_lambdasq_delta <- (
          -0.5 * (log(lambdasq_delta_prop_tilde) - log(lambdasq_delta_update_s_1_tilde))
          - (3/2) * (log(lambdasq_delta_prop) - log(lambdasq_delta_update_s[jm, k]))
          - ((delta_update[jm, k])^2 / (2 * tausq_delta_update[jm])) *
            ((1/lambdasq_delta_prop_tilde) - (1/lambdasq_delta_update_s_1_tilde))
          - (1/psi_delta_update[jm, k]) * ((1/lambdasq_delta_prop) - (1/lambdasq_delta_update_s[jm, k]))
        )

        # Accept/reject the proposed lambdasq_delta
        if (log(runif(1)) < log_accept_ratio_lambdasq_delta)
        {
          lambdasq_delta_update_s[jm, k] <- lambdasq_delta_prop
          accept_lambdasq_delta_update_s[jm, k] <- 1
        }
      }
    }
  }
  return(list(lambdasq_delta_update = lambdasq_delta_update_s,
              accept_lambdasq_delta_update = accept_lambdasq_delta_update_s))
}



#######################################################################################

# Update tausq_beta (for all j = 1, 2, ..., J)
update_tausq_beta_NHRHS_MVN_cov <- function(c, J, K, beta_update,
                                                lambdasq_beta_update, xi_beta_update,
                                                tausq_beta_update_s_1, omegasq_beta_tausq)
{
  # MH update for tausq_beta
  tausq_beta_update_s <- tausq_beta_update_s_1
  accept_tausq_beta_update_s <- rep(0, times = J)

  for (j in 1:J)
  {
    tausq_beta_prop <- rfnorm(1, tausq_beta_update_s[j], sqrt(omegasq_beta_tausq)) # Folded normal

    lambdasq_beta_prop_tilde <- ((c^2 * lambdasq_beta_update[j, ]) /
                                   (c^2 + tausq_beta_prop * lambdasq_beta_update[j, ]))
    lambdasq_beta_update_s_1_tilde <- ((c^2 * lambdasq_beta_update[j, ]) /
                                         (c^2 + tausq_beta_update_s[j] * lambdasq_beta_update[j, ]))

    log_accept_ratio_tausq_beta <- (
      - ((K+3)/2) * (log(tausq_beta_prop) - log(tausq_beta_update_s[j]))
      - (1/tausq_beta_prop) * ((1/xi_beta_update) +
                                 sum(beta_update[j, ]^2 / (2*lambdasq_beta_prop_tilde)))
      + (1/tausq_beta_update_s[j]) * ((1/xi_beta_update) +
                                        sum(beta_update[j, ]^2 / (2*lambdasq_beta_update_s_1_tilde)))
    )

    # Accept/reject the proposed tausq_beta
    if (log(runif(1)) < log_accept_ratio_tausq_beta)
    {
      tausq_beta_update_s[j] <- tausq_beta_prop
      accept_tausq_beta_update_s[j] <- 1
    }
  }
  return(list(tausq_beta_update = tausq_beta_update_s,
              accept_tausq_beta_update = accept_tausq_beta_update_s))
}




#######################################################################################

# Update tausq_gamma (for all m = 1, 2, ..., M)
update_tausq_gamma_NHRHS_MVN_cov <- function(c, M, K, gamma_update,
                                                 lambdasq_gamma_update, xi_gamma_update,
                                                 tausq_gamma_update_s_1, omegasq_gamma_tausq)
{
  # MH update for tausq_gamma
  tausq_gamma_update_s <- tausq_gamma_update_s_1
  accept_tausq_gamma_update_s <- rep(0, times = M)

  for (m in 1:M)
  {
    tausq_gamma_prop <- rfnorm(1, tausq_gamma_update_s[m], sqrt(omegasq_gamma_tausq)) # Folded normal

    # a vector
    lambdasq_gamma_prop_tilde <- ((c^2 * lambdasq_gamma_update[m, ]) /
                                    (c^2 + tausq_gamma_prop * lambdasq_gamma_update[m, ]))
    lambdasq_gamma_update_s_1_tilde <- ((c^2 * lambdasq_gamma_update[m, ]) /
                                          (c^2 + tausq_gamma_update_s[m] * lambdasq_gamma_update[m, ]))

    log_accept_ratio_tausq_gamma <- (
      - ((K+3)/2) * (log(tausq_gamma_prop) - log(tausq_gamma_update_s[m]))
      - (1/tausq_gamma_prop) * ((1/xi_gamma_update) +
                                  sum(gamma_update[m, ]^2 / (2*lambdasq_gamma_prop_tilde)))
      + (1/tausq_gamma_update_s[m]) * ((1/xi_gamma_update) +
                                         sum(gamma_update[m, ]^2 / (2*lambdasq_gamma_update_s_1_tilde)))
    )

    # Accept/reject the proposed tausq_gamma
    if (log(runif(1)) < log_accept_ratio_tausq_gamma)
    {
      tausq_gamma_update_s[m] <- tausq_gamma_prop
      accept_tausq_gamma_update_s[m] <- 1
    }
  }
  return(list(tausq_gamma_update = tausq_gamma_update_s,
              accept_tausq_gamma_update = accept_tausq_gamma_update_s))
}


#######################################################################################

# Update tausq_delta (for all j, and m )
update_tausq_delta_NHRHS_MVN_cov <- function(c, J, M, K, delta_update, xi_delta_update,
                                                 lambdasq_delta_update,
                                                 tausq_delta_update_s_1,
                                                 omegasq_delta_tausq)
{
  # MH update for tausq_delta
  tausq_delta_update_s <- tausq_delta_update_s_1
  accept_tausq_delta_update_s <- rep(0, times = J*M)

  for (j in 1:J)
  {
    for(m in 1:M)
    {
      jm <- (j - 1) * M + m

      tausq_delta_prop <- rfnorm(1, tausq_delta_update_s[jm], sqrt(omegasq_delta_tausq)) # Folded normal

      # a vector
      lambdasq_delta_prop_tilde <- ((c^2 * lambdasq_delta_update[jm, ]) /
                                      (c^2 + tausq_delta_prop * lambdasq_delta_update[jm, ]))
      lambdasq_delta_update_s_1_tilde <- ((c^2 * lambdasq_delta_update[jm, ]) /
                                            (c^2 + tausq_delta_update_s[jm] * lambdasq_delta_update[jm, ]))

      log_accept_ratio_tausq_delta <- (
        - ((K+3)/2) * (log(tausq_delta_prop) - log(tausq_delta_update_s[jm]))
        - (1/tausq_delta_prop) * ((1/xi_delta_update) +
                                    sum(delta_update[jm, ]^2 / (2 * lambdasq_delta_prop_tilde)))
        + (1/tausq_delta_update_s[jm]) * ((1/xi_delta_update) +
                                            sum(delta_update[jm, ]^2 / (2 * lambdasq_delta_update_s_1_tilde)))
      )

      # Accept/reject the proposed tausq_delta
      if (log(runif(1)) < log_accept_ratio_tausq_delta)
      {
        tausq_delta_update_s[jm] <- tausq_delta_prop
        accept_tausq_delta_update_s[jm] <- 1
      }
    }
  }
  return(list(tausq_delta_update = tausq_delta_update_s,
              accept_tausq_delta_update = accept_tausq_delta_update_s))
}


#######################################################################################

# Update psi_beta (for all j = 1, 2, ..., J, and k = 1, 2, ..., K)
update_psi_beta_NHRHS_MVN_cov <- function(J, K, lambdasq_beta_update)
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
update_psi_gamma_NHRHS_MVN_cov <- function(M, K, lambdasq_gamma_update)
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
update_psi_delta_NHRHS_MVN_cov <- function(J, M, K, lambdasq_delta_update)
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
update_xi_beta_NHRHS_MVN_cov <- function(J, tausq_beta_update)
{
  # Use more informative prior (xi_beta ~ IG(15, 3))
  shape_xi_beta <- ((J/2) + 15)
  scale_xi_beta <- (3 + sum(1 / tausq_beta_update))

  xi_beta_update_s <- rinvgamma(1, shape = shape_xi_beta, scale = scale_xi_beta) # IG distribution
  return(xi_beta_update_s)
}



#######################################################################################

# Update xi_gamma
update_xi_gamma_NHRHS_MVN_cov <- function(M, tausq_gamma_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_gamma <- ((M/2) + 15)
  scale_xi_gamma <- (3 + sum(1 / tausq_gamma_update))

  xi_gamma_update_s <- rinvgamma(1, shape = shape_xi_gamma, scale = scale_xi_gamma) # IG distribution
  return(xi_gamma_update_s)
}





#######################################################################################

# Update xi_delta
update_xi_delta_NHRHS_MVN_cov <- function(J, M, tausq_delta_update)
{
  # Use more informative prior (xi_delta ~ IG(15, 3))
  shape_xi_delta <- ((J*M/2) + 15)
  scale_xi_delta <- (3 + sum(1 / tausq_delta_update))

  xi_delta_update_s <- rinvgamma(1, shape = shape_xi_delta, scale = scale_xi_delta) # IG distribution
  return(xi_delta_update_s)
}



##################################################################################


# Update parameters for NHRHS-MVN model
fit_NHRHS_MVN_cov_modify3 <- function(niter = 6000, burn_in = 1000, thin = 5,
                                  n, K, Y, W, WTW, n_all_par, J, M, O,
                                  c = 2.5,
                                  sigmasq_alpha = 100,
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
                                  omegasq_beta_lambdasq = 0.25, omegasq_beta_tausq = 0.25,
                                  omegasq_gamma_lambdasq = 0.25, omegasq_gamma_tausq = 0.25,
                                  omegasq_delta_lambdasq = 0.25, omegasq_delta_tausq = 0.25,
                                  accept_lambdasq_beta_init = matrix(1, nrow = J, ncol = K),
                                  accept_tausq_beta_init = rep(1, times = J),
                                  accept_lambdasq_gamma_init = matrix(1, nrow = M, ncol = K),
                                  accept_tausq_gamma_init = rep(1, times = M),
                                  accept_lambdasq_delta_init = matrix(1, nrow = J*M, ncol = K),
                                  accept_tausq_delta_init = rep(1, times = J*M),
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

  accept_lambdasq_beta_update <- array(NA, dim = c(niter, J, K))
  accept_tausq_beta_update <- matrix(NA, nrow = niter, ncol = J)
  accept_lambdasq_gamma_update <- array(NA, dim = c(niter, M, K))
  accept_tausq_gamma_update <- matrix(NA, nrow = niter, ncol = M)
  accept_lambdasq_delta_update <- array(NA, dim = c(niter, J*M, K))
  accept_tausq_delta_update <- matrix(NA, nrow = niter, ncol = J*M)

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


  accept_lambdasq_beta_update[1,,] <- accept_lambdasq_beta_init
  accept_tausq_beta_update[1,] <- accept_tausq_beta_init
  accept_lambdasq_gamma_update[1,,] <- accept_lambdasq_gamma_init
  accept_tausq_gamma_update[1,] <- accept_tausq_gamma_init
  accept_lambdasq_delta_update[1,,] <- accept_lambdasq_delta_init
  accept_tausq_delta_update[1,] <- accept_tausq_delta_init

  alpha0_update[1, ] <- theta_update[1, 1, ]
  beta_update[1, , ] <- theta_update[1, (1 + 1):(1 + J), ]
  gamma_update[1, , ] <- theta_update[1, (1 + J + 1):(1 + J + M), ]
  delta_update[1, , ] <- theta_update[1, (1 + J + M + 1):(1 + J + M + J*M), ]
  varphi_update[1, , ] <- theta_update[1, (1 + J + M + J*M + 1):(1 + J + M + J*M + O), ]



  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_NHRHS_MVN_cov_modify(Y, K, W, WTW, n_all_par,
                                                     J, M, O,
                                                     c,
                                                     sigmasq_alpha = 100,
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
    Sigma_update_s <- update_Sigma_NHRHS_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s


    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_NHRHS_MVN_cov(c, J, K,
                                                                     beta_update = beta_update[s, , ],
                                                                     psi_beta_update = psi_beta_update[(s-1), , ],
                                                                     tausq_beta_update = tausq_beta_update[(s-1), ],
                                                                     lambdasq_beta_update_s_1 = lambdasq_beta_update[(s-1),,],
                                                                     omegasq_beta_lambdasq)
    lambdasq_beta_update[s, , ] <- lambdasq_beta_update_s$lambdasq_beta_update
    accept_lambdasq_beta_update[s, ,] <- lambdasq_beta_update_s$accept_lambdasq_beta_update


    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_NHRHS_MVN_cov(c, M, K,
                                                                       gamma_update = gamma_update[s, , ],
                                                                       psi_gamma_update = psi_gamma_update[(s-1), , ],
                                                                       tausq_gamma_update = tausq_gamma_update[(s-1), ],
                                                                       lambdasq_gamma_update_s_1 = lambdasq_gamma_update[(s-1),,],
                                                                       omegasq_gamma_lambdasq)
    lambdasq_gamma_update[s, , ] <- lambdasq_gamma_update_s$lambdasq_gamma_update
    accept_lambdasq_gamma_update[s, ,] <- lambdasq_gamma_update_s$accept_lambdasq_gamma_update

    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_NHRHS_MVN_cov(c, J, M, K,
                                                                       delta_update = delta_update[s, , ],
                                                                       psi_delta_update = psi_delta_update[(s-1), , ],
                                                                       tausq_delta_update = tausq_delta_update[(s-1), ],
                                                                       lambdasq_delta_update_s_1 = lambdasq_delta_update[(s-1),,],
                                                                       omegasq_delta_lambdasq)
    lambdasq_delta_update[s, , ] <- lambdasq_delta_update_s$lambdasq_delta_update
    accept_lambdasq_delta_update[s, ,] <- lambdasq_delta_update_s$accept_lambdasq_delta_update



    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_NHRHS_MVN_cov(c, J, K,
                                                               beta_update = beta_update[s, , ],
                                                               lambdasq_beta_update = lambdasq_beta_update[s, , ],
                                                               xi_beta_update = xi_beta_update[(s-1)],
                                                               tausq_beta_update_s_1 = tausq_beta_update[(s-1), ],
                                                               omegasq_beta_tausq)
    tausq_beta_update[s, ] <- tausq_beta_update_s$tausq_beta_update
    accept_tausq_beta_update[s,] <- tausq_beta_update_s$accept_tausq_beta_update


    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_NHRHS_MVN_cov(c, M, K,
                                                                 gamma_update = gamma_update[s, , ],
                                                                 lambdasq_gamma_update = lambdasq_gamma_update[s, , ],
                                                                 xi_gamma_update = xi_gamma_update[(s-1)],
                                                                 tausq_gamma_update_s_1 = tausq_gamma_update[(s-1), ],
                                                                 omegasq_gamma_tausq)
    tausq_gamma_update[s, ] <- tausq_gamma_update_s$tausq_gamma_update
    accept_tausq_gamma_update[s,] <- tausq_gamma_update_s$accept_tausq_gamma_update


    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_NHRHS_MVN_cov(c, J, M, K,
                                                                 delta_update = delta_update[s, , ],
                                                                 xi_delta_update = xi_delta_update[(s-1)],
                                                                 lambdasq_delta_update = lambdasq_delta_update[s, , ],
                                                                 tausq_delta_update_s_1 = tausq_delta_update[(s-1),],
                                                                 omegasq_delta_tausq)
    tausq_delta_update[s, ] <- tausq_delta_update_s$tausq_delta_update
    accept_tausq_delta_update[s,] <- tausq_delta_update_s$accept_tausq_delta_update


    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_NHRHS_MVN_cov(J, K,
                                                           lambdasq_beta_update = lambdasq_beta_update[s, , ])
    psi_beta_update[s, , ] <- psi_beta_update_s


    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_NHRHS_MVN_cov(M, K,
                                                             lambdasq_gamma_update = lambdasq_gamma_update[s, , ])
    psi_gamma_update[s, , ] <- psi_gamma_update_s


    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_NHRHS_MVN_cov(J, M, K,
                                                             lambdasq_delta_update = lambdasq_delta_update[s, , ])
    psi_delta_update[s, , ] <- psi_delta_update_s



    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_NHRHS_MVN_cov(J, tausq_beta_update = tausq_beta_update[s, ])
    xi_beta_update[s] <- xi_beta_update_s


    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_NHRHS_MVN_cov(M, tausq_gamma_update = tausq_gamma_update[s, ])
    xi_gamma_update[s] <- xi_gamma_update_s


    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_NHRHS_MVN_cov(J, M, tausq_delta_update = tausq_delta_update[s, ])
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
    xi_delta_update = xi_delta_update[indices_to_save],

    accept_lambdasq_beta_update = accept_lambdasq_beta_update[indices_to_save,,],
    accept_tausq_beta_update = accept_tausq_beta_update[indices_to_save,],
    accept_lambdasq_gamma_update = accept_lambdasq_gamma_update[indices_to_save,,],
    accept_tausq_gamma_update = accept_tausq_gamma_update[indices_to_save,],
    accept_lambdasq_delta_update = accept_lambdasq_delta_update[indices_to_save,,],
    accept_tausq_delta_update = accept_tausq_delta_update[indices_to_save,]
  )
  return(results)
}

##################################################################



