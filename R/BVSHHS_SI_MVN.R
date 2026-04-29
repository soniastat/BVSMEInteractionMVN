
# Bayesian variable selection model for MVN response
# Parallelization

# Bayesian variable selection with hierarchical Horseshoe priors and
# Shared information across outcomes for multivariate normal response
# (BVSHHS-SI-MVN)


##############################################################
# Matrix
#' @importFrom Matrix bdiag Diagonal chol2inv

# MASS
#' @importFrom MASS mvrnorm

# mvtnorm
#' @importFrom mvtnorm dmvnorm

# actuar
#' @importFrom actuar rinvgamma

NULL
####################################################################################

# Update theta
update_theta_BVSHHS_SI_MVN <- function(Y, K, W, n_all_par,
                         J, M,
                         lambdasq_beta_update, tausq_beta_update,
                         lambdasq_gamma_update, tausq_gamma_update,
                         lambdasq_delta_update, tausq_delta_update,
                         Sigma_update
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
  mat_var_alpha <- matrix(100, nrow = 1, ncol = 1)
  var_beta <- diag(sigmasq_beta_update)
  var_gamma <- diag(sigmasq_gamma_update)
  var_delta <- diag(sigmasq_delta_update)

  V_theta <- as.matrix(bdiag(mat_var_alpha, var_beta, var_gamma, var_delta))

  # Prior precision
  V_theta_inv <- chol2inv(chol(V_theta))
  prior_prec <- kronecker(Diagonal(K), V_theta_inv)  # I_K ⊗ Vθ⁻¹

  # Kronecker term from likelihood
  Sigma_inv <- as.matrix(chol2inv(chol(Sigma_update)))
  kron_term <- kronecker(Sigma_inv, t(W) %*% W)  # Σ⁻¹ ⊗ WᵀW

  # Full conditional covariance
  Sigma_theta <- chol2inv(chol(kron_term + prior_prec))

  # Full conditional mean
  vecY <- as.vector(Y)
  mean_theta_vec <- Sigma_theta %*% (kronecker(Sigma_inv, t(W)) %*% vecY)

  # Sample from MVN
  theta_vec <- mvrnorm(1, mu = mean_theta_vec, Sigma = Sigma_theta)
  # Reshape back to n_all_par x K
  theta_update_s <- matrix(theta_vec, nrow = n_all_par, ncol = K)


  return(theta_update = theta_update_s)
}

#######################################################################################
# Update Sigma
update_Sigma_BVSHHS_SI_MVN <- function(Y, n, W, theta_update, Psi_0, nu_0)
{
  # Residuals
  Resid <- Y - W %*% theta_update

  # Posterior scale and df
  Psi_n <- Psi_0 + t(Resid) %*% Resid
  nu_n <- nu_0 + n

  # Draw from IW
  Sigma_update_s <- rinvwishart(nu_n, Psi_n)

  return(Sigma_update = Sigma_update_s)
}

#######################################################################################

# Update lambdasq_beta (for all j = 1, 2, ..., J)
update_lambdasq_beta_BVSHHS_SI_MVN <- function(J, K, beta_update, psi_beta_update, tausq_beta_update)
{
  lambdasq_beta_update_s <- rep(NA, J)
  shape_lambdasq_beta <- (K+1)/2

  for(j in 1:J)
  {
    scale_lambdasq_beta <- ((1/psi_beta_update[j]) +
                              (sum((beta_update[j, ])^2) / (2 * tausq_beta_update)))
    lambdasq_beta_update_s[j] <- rinvgamma(1, shape = shape_lambdasq_beta, scale = scale_lambdasq_beta) # IG distribution
  }
  return(lambdasq_beta_update = lambdasq_beta_update_s)
}


#######################################################################################

# Update lambdasq_gamma (for all m = 1, 2, ..., M)
update_lambdasq_gamma_BVSHHS_SI_MVN <- function(M, K, gamma_update, psi_gamma_update, tausq_gamma_update)
{
  lambdasq_gamma_update_s <- rep(NA, M)
  shape_lambdasq_gamma <- (K+1)/2

  for(m in 1:M)
  {
    scale_lambdasq_gamma <- ((1/psi_gamma_update[m]) +
                               (sum((gamma_update[m, ])^2) / (2 * tausq_gamma_update)))
    lambdasq_gamma_update_s[m] <- rinvgamma(1, shape = shape_lambdasq_gamma, scale = scale_lambdasq_gamma) # IG distribution
  }
  return(lambdasq_gamma_update = lambdasq_gamma_update_s)
}

#####################################################################################

# Update lambdasq_delta (for all j and m)
update_lambdasq_delta_BVSHHS_SI_MVN <- function(J, M, K, delta_update, lambdasq_beta_update,
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
  return(lambdasq_delta_update = lambdasq_delta_update_s)
}

#######################################################################################

# Update tausq_beta
update_tausq_beta_BVSHHS_SI_MVN <- function(J, K, beta_update, lambdasq_beta_update, xi_beta_update)
{

  sum_term_tausq_beta <- 0
  for(j in 1:J)
  {
    for (k in 1:K)
    {
  sum_term_tausq_beta <- sum_term_tausq_beta + (beta_update[j, k]^2 / (2 * lambdasq_beta_update[j]))
    }
  }
  shape_tausq_beta <- (J*K + 1) / 2
  scale_tausq_beta <- ((1/xi_beta_update) + sum_term_tausq_beta)

  tausq_beta_update_s <- rinvgamma(1, shape = shape_tausq_beta, scale = scale_tausq_beta) # IG distribution

  return(tausq_beta_update = tausq_beta_update_s)
}

#######################################################################################

# Update tausq_gamma
update_tausq_gamma_BVSHHS_SI_MVN <- function(M, K, gamma_update, lambdasq_gamma_update, xi_gamma_update)
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

  return(tausq_gamma_update = tausq_gamma_update_s)
}


#######################################################################################

# Update tausq_delta
update_tausq_delta_BVSHHS_SI_MVN <- function(J, M, K, delta_update, lambdasq_delta_update,
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

  return(tausq_delta_update = tausq_delta_update_s)
}

#######################################################################################

# Update psi_beta (for all j = 1, 2, ..., J)
update_psi_beta_BVSHHS_SI_MVN <- function(J, K, lambdasq_beta_update)
{
  psi_beta_update_s <- rep(NA, J)

  for(j in 1:J)
  {
    scale_psi_beta <- (1 + (1/lambdasq_beta_update[j]))
    psi_beta_update_s[j] <- rinvgamma(1, shape = 1, scale = scale_psi_beta) # IG distribution
  }
  return(psi_beta_update = psi_beta_update_s)
}


#######################################################################################

# Update psi_gamma (for all m = 1, 2, ..., M)
update_psi_gamma_BVSHHS_SI_MVN <- function(M, K, lambdasq_gamma_update)
{
  psi_gamma_update_s <- rep(NA, M)

  for(m in 1:M)
  {
    scale_psi_gamma <- (1 + (1/lambdasq_gamma_update[m]))
    psi_gamma_update_s[m] <- rinvgamma(1, shape = 1, scale = scale_psi_gamma) # IG distribution
  }
  return(psi_gamma_update = psi_gamma_update_s)
}



#######################################################################################

# Update psi_delta (for all j , m)
update_psi_delta_BVSHHS_SI_MVN <- function(J, M, K, lambdasq_delta_update)
{
  psi_delta_update_s <- rep(NA, J*M)

  for(jm in 1:(J*M))
  {
    scale_psi_delta <- (1 + (1/lambdasq_delta_update[jm]))
    psi_delta_update_s[jm] <- rinvgamma(1, shape = 1, scale = scale_psi_delta) # IG distribution
  }
  return(psi_delta_update = psi_delta_update_s)
}

#######################################################################################

# Update xi_beta
update_xi_beta_BVSHHS_SI_MVN <- function(J, tausq_beta_update)
{
  # Use more informative prior (xi_beta ~ IG(15, 3))
  shape_xi_beta <- ((1/2) + 15)
  scale_xi_beta <- (3 + (1 / tausq_beta_update))

  xi_beta_update_s <- rinvgamma(1, shape = shape_xi_beta, scale = scale_xi_beta) # IG distribution
  return(xi_beta_update = xi_beta_update_s)
}



#######################################################################################

# Update xi_gamma
update_xi_gamma_BVSHHS_SI_MVN <- function(M, tausq_gamma_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_gamma <- ((1/2) + 15)
  scale_xi_gamma <- (3 + (1 / tausq_gamma_update))

  xi_gamma_update_s <- rinvgamma(1, shape = shape_xi_gamma, scale = scale_xi_gamma) # IG distribution
  return(xi_gamma_update = xi_gamma_update_s)
}



#######################################################################################

# Update xi_delta
update_xi_delta_BVSHHS_SI_MVN <- function(J, M, tausq_delta_update)
{
  # Use more informative prior (xi_delta ~ IG(15, 3))
  shape_xi_delta <- ((1/2) + 15)
  scale_xi_delta <- (3 + (1 / tausq_delta_update))

  xi_delta_update_s <- rinvgamma(1, shape = shape_xi_delta, scale = scale_xi_delta) # IG distribution
  return(xi_delta_update = xi_delta_update_s)
}

##################################################################################

#' Update parameters for BVSHHS-SI-MVN model
#' @export
fit_BVSHHS_SI_MVN <- function(niter = 6000, burn_in = 1000, thin = 5,
                            n, K, Y, W, n_all_par, J, M,
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
                            Sigma_init,
                            nu_0, Psi_0)
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



  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_BVSHHS_SI_MVN(Y, K, W, n_all_par,
                                   J, M,
                                   lambdasq_beta_update = lambdasq_beta_update[(s-1), ],
                                   tausq_beta_update = tausq_beta_update[(s-1)],
                                   lambdasq_gamma_update = lambdasq_gamma_update[(s-1), ],
                                   tausq_gamma_update = tausq_gamma_update[(s-1)],
                                   lambdasq_delta_update = lambdasq_delta_update[(s-1), ],
                                   tausq_delta_update = tausq_delta_update[(s-1)]
                                   , Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s

    alpha0_update_s <- theta_update[s, 1, ]
    beta_update_s <- theta_update[s, (1 + 1):(1 + J), ]
    gamma_update_s <- theta_update[s, (1 + J + 1):(1 + J + M), ]
    delta_update_s <- theta_update[s, (1 + J + M + 1):(1 + J + M + J*M), ]

    alpha0_update[s, ] <- alpha0_update_s
    beta_update[s, , ] <- beta_update_s
    gamma_update[s, , ] <- gamma_update_s
    delta_update[s, , ] <- delta_update_s


    # Sigma_update
    Sigma_update_s <- update_Sigma_BVSHHS_SI_MVN(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s


    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_BVSHHS_SI_MVN(J, K,
                                                   beta_update = beta_update[s, , ],
                                                   psi_beta_update = psi_beta_update[(s-1), ],
                                                   tausq_beta_update = tausq_beta_update[(s-1)])
    lambdasq_beta_update[s, ] <- lambdasq_beta_update_s

    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_BVSHHS_SI_MVN(M, K,
                                                     gamma_update = gamma_update[s, , ],
                                                     psi_gamma_update = psi_gamma_update[(s-1), ],
                                                     tausq_gamma_update = tausq_gamma_update[(s-1)])
    lambdasq_gamma_update[s, ] <- lambdasq_gamma_update_s

    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_BVSHHS_SI_MVN(J, M, K,
                                                     delta_update = delta_update[s, , ],
                                                     lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                     lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                     psi_delta_update = psi_delta_update[(s-1), ],
                                                     tausq_delta_update = tausq_delta_update[(s-1)])
    lambdasq_delta_update[s, ] <- lambdasq_delta_update_s



    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_BVSHHS_SI_MVN(J, K,
                                             beta_update = beta_update[s, , ],
                                             lambdasq_beta_update = lambdasq_beta_update[s, ],
                                             xi_beta_update = xi_beta_update[(s-1)])
    tausq_beta_update[s] <- tausq_beta_update_s


    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_BVSHHS_SI_MVN(M, K,
                                               gamma_update = gamma_update[s, , ],
                                               lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                               xi_gamma_update = xi_gamma_update[(s-1)])
    tausq_gamma_update[s] <- tausq_gamma_update_s


    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_BVSHHS_SI_MVN(J, M, K,
                                               delta_update = delta_update[s, , ],
                                               lambdasq_delta_update = lambdasq_delta_update[s, ],
                                               lambdasq_beta_update = lambdasq_beta_update[s, ],
                                               lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                               xi_delta_update = xi_delta_update[(s-1)])
    tausq_delta_update[s] <- tausq_delta_update_s



    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_BVSHHS_SI_MVN(J, K,
                                         lambdasq_beta_update = lambdasq_beta_update[s, ])
    psi_beta_update[s, ] <- psi_beta_update_s


    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_BVSHHS_SI_MVN(M, K,
                                           lambdasq_gamma_update = lambdasq_gamma_update[s, ])
    psi_gamma_update[s, ] <- psi_gamma_update_s


    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_BVSHHS_SI_MVN(J, M, K,
                                           lambdasq_delta_update = lambdasq_delta_update[s, ])
    psi_delta_update[s, ] <- psi_delta_update_s



    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_BVSHHS_SI_MVN(J, tausq_beta_update = tausq_beta_update[s])
    xi_beta_update[s] <- xi_beta_update_s


    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_BVSHHS_SI_MVN(M, tausq_gamma_update = tausq_gamma_update[s])
    xi_gamma_update[s] <- xi_gamma_update_s


    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_BVSHHS_SI_MVN(J, M, tausq_delta_update = tausq_delta_update[s])
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


