
# Bayesian variable selection with hierarchical Horseshoe priors and
# Shared information across outcomes for multivariate normal response
# (HHSSI-MVN-cov)
# seperate covariates


##############################################################
# Matrix
#' @importFrom Matrix bdiag Diagonal chol2inv

# invwishart
#' @importFrom LaplacesDemon rinvwishart


NULL

####################################################################################

# Update theta
update_theta_HHSSI_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par,
                                              J, M, O,
                                              sigmasq_alpha = 100,
                                              lambdasq_beta_update, tausq_beta_update,
                                              lambdasq_gamma_update, tausq_gamma_update,
                                              lambdasq_delta_update, tausq_delta_update,
                                              sigmasq_varphi, Sigma_update
)
{
  sigmasq_beta_update <- lambdasq_beta_update * tausq_beta_update
  sigmasq_gamma_update <- lambdasq_gamma_update * tausq_gamma_update

  beta_rep  <- rep(lambdasq_beta_update, each = M)
  gamma_rep <- rep(lambdasq_gamma_update, times = J)
  sigmasq_delta_update <- (lambdasq_delta_update *
                             beta_rep *
                             gamma_rep *
                             tausq_delta_update)

  vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)

  # Inverse using blocks separately
  V_theta_inv <- as.matrix(
    bdiag(
      1/sigmasq_alpha,
      diag(1/sigmasq_beta_update),
      diag(1/sigmasq_gamma_update),
      diag(1/sigmasq_delta_update),
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
update_Sigma_HHSSI_MVN_cov <- function(Y, n, W, theta_update, Psi_0, nu_0)
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
update_lambdasq_beta_HHSSI_MVN_cov2 <- function(J, K, beta_update,
                                                psi_beta_update,
                                                tausq_beta_update)
{
  shape_lambdasq_beta <- (K + 1) / 2

  scale_lambdasq_beta <- (1 / psi_beta_update) +
    rowSums(beta_update^2) / (2 * tausq_beta_update)

  lambdasq_beta_update_s <- rinvgamma(J,
                                      shape = shape_lambdasq_beta,
                                      scale = scale_lambdasq_beta)
  return(lambdasq_beta_update_s)
}

#######################################################################################

# Update lambdasq_gamma (for all m = 1, 2, ..., M)
update_lambdasq_gamma_HHSSI_MVN_cov2 <- function(M, K, gamma_update,
                                                 psi_gamma_update,
                                                 tausq_gamma_update)
{
  shape_lambdasq_gamma <- (K + 1) / 2

  scale_lambdasq_gamma <- (1 / psi_gamma_update) +
    rowSums(gamma_update^2) / (2 * tausq_gamma_update)

  lambdasq_gamma_update_s <- rinvgamma(M,
                                       shape = shape_lambdasq_gamma,
                                       scale = scale_lambdasq_gamma)
  return(lambdasq_gamma_update_s)
}

#####################################################################################

# Update lambdasq_delta (for all j and m)
update_lambdasq_delta_HHSSI_MVN_cov2 <- function(
    J, M, K, delta_update,
    lambdasq_beta_update,
    lambdasq_gamma_update,
    psi_delta_update,
    tausq_delta_update)
{
  shape_lambdasq_delta <- (K + 1) / 2

  # Sum of squares for each (j,m)
  delta_ss <- rowSums(delta_update^2)

  # Repeat beta and gamma to match the ordering of delta_update
  beta_rep  <- rep(lambdasq_beta_update, each = M)
  gamma_rep <- rep(lambdasq_gamma_update, times = J)

  scale_lambdasq_delta <-
    1 / psi_delta_update +
    delta_ss / (2 * beta_rep * gamma_rep * tausq_delta_update)

  lambdasq_delta_update_s <- rinvgamma(
    J * M,
    shape = shape_lambdasq_delta,
    scale = scale_lambdasq_delta
  )
  return(lambdasq_delta_update_s)
}

#######################################################################################

# Update tausq_beta
update_tausq_beta_HHSSI_MVN_cov2 <- function(
    J, K, beta_update,
    lambdasq_beta_update,
    xi_beta_update)
{
  tausq_beta_update_s <- rinvgamma(
    1,
    shape = (J * K + 1) / 2,
    scale = 1 / xi_beta_update +
      sum(rowSums(beta_update^2) / (2 * lambdasq_beta_update))
  )
  return(tausq_beta_update_s)
}



#######################################################################################

# Update tausq_gamma
update_tausq_gamma_HHSSI_MVN_cov2 <- function(
    M, K, gamma_update,
    lambdasq_gamma_update,
    xi_gamma_update)
{
  tausq_gamma_update_s <- rinvgamma(
    1,
    shape = (M * K + 1) / 2,
    scale = 1 / xi_gamma_update +
      sum(rowSums(gamma_update^2) / (2 * lambdasq_gamma_update))
  )
  return(tausq_gamma_update_s)
}


#######################################################################################

# Update tausq_delta
update_tausq_delta_HHSSI_MVN_cov2 <- function(
    J, M, K,
    delta_update,
    lambdasq_delta_update,
    lambdasq_beta_update,
    lambdasq_gamma_update,
    xi_delta_update)
{
  # Repeat beta and gamma to match rows of delta_update
  beta_rep  <- rep(lambdasq_beta_update, each = M)
  gamma_rep <- rep(lambdasq_gamma_update, times = J)

  # Sum of squares for each (j,m)
  delta_ss <- rowSums(delta_update * delta_update)

  sum_term_tausq_delta <- sum(
    delta_ss /
      (2 * lambdasq_delta_update * beta_rep * gamma_rep)
  )

  shape_tausq_delta <- (J * M * K + 1) / 2
  scale_tausq_delta <- 1 / xi_delta_update + sum_term_tausq_delta

  tausq_delta_update_s <- rinvgamma(
    1,
    shape = shape_tausq_delta,
    scale = scale_tausq_delta
  )
  return(tausq_delta_update_s)
}

#######################################################################################

# Update psi_beta (for all j = 1, 2, ..., J)
update_psi_beta_HHSSI_MVN_cov2 <- function(J, lambdasq_beta_update)
{
  psi_beta_update_s <- rinvgamma(
    n = J,
    shape = 1,
    scale = 1 + 1 / lambdasq_beta_update
  )
  return(psi_beta_update_s)
}


#######################################################################################

# Update psi_gamma (for all m = 1, 2, ..., M)
update_psi_gamma_HHSSI_MVN_cov2 <- function(M, lambdasq_gamma_update)
{
  psi_gamma_update_s <- rinvgamma(
    n = M,
    shape = 1,
    scale = 1 + 1 / lambdasq_gamma_update
  )
  return(psi_gamma_update_s)
}


#######################################################################################

# Update psi_delta (for all j , m)
update_psi_delta_HHSSI_MVN_cov2 <- function(
    J, M,
    lambdasq_delta_update)
{
  psi_delta_update_s <- rinvgamma(
    n = J * M,
    shape = 1,
    scale = 1 + 1 / lambdasq_delta_update
  )
  return(psi_delta_update_s)
}

#######################################################################################

# Update xi_beta
update_xi_beta_HHSSI_MVN_cov <- function(J, tausq_beta_update)
{
  # Use more informative prior (xi_beta ~ IG(15, 3))
  shape_xi_beta <- ((1/2) + 15)
  scale_xi_beta <- (3 + (1 / tausq_beta_update))

  xi_beta_update_s <- rinvgamma(1, shape = shape_xi_beta, scale = scale_xi_beta) # IG distribution
  return(xi_beta_update_s)
}



#######################################################################################

# Update xi_gamma
update_xi_gamma_HHSSI_MVN_cov <- function(M, tausq_gamma_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_gamma <- ((1/2) + 15)
  scale_xi_gamma <- (3 + (1 / tausq_gamma_update))

  xi_gamma_update_s <- rinvgamma(1, shape = shape_xi_gamma, scale = scale_xi_gamma) # IG distribution
  return(xi_gamma_update_s)
}



#######################################################################################

# Update xi_delta
update_xi_delta_HHSSI_MVN_cov <- function(J, M, tausq_delta_update)
{
  # Use more informative prior (xi_delta ~ IG(15, 3))
  shape_xi_delta <- ((1/2) + 15)
  scale_xi_delta <- (3 + (1 / tausq_delta_update))

  xi_delta_update_s <- rinvgamma(1, shape = shape_xi_delta, scale = scale_xi_delta) # IG distribution
  return(xi_delta_update_s)
}


##################################################################################

#' Update parameters for HHSSI-MVN model
#' @export
fit_HHSSI_MVN_cov_modify3 <- function(niter = 6000, burn_in = 1000, thin = 5,
                                      n, K, Y, W, WTW, n_all_par, J, M, O,
                                      sigmasq_alpha = 100,
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
    theta_update_s <- update_theta_HHSSI_MVN_cov_modify(Y, K, W, WTW, n_all_par,
                                                        J, M, O,
                                                        sigmasq_alpha = 100,
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
    Sigma_update_s <- update_Sigma_HHSSI_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s


    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_HHSSI_MVN_cov2(J, K,
                                                                  beta_update = beta_update[s, , ],
                                                                  psi_beta_update = psi_beta_update[(s-1), ],
                                                                  tausq_beta_update = tausq_beta_update[(s-1)])
    lambdasq_beta_update[s, ] <- lambdasq_beta_update_s


    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_HHSSI_MVN_cov2(M, K,
                                                                    gamma_update = gamma_update[s, , ],
                                                                    psi_gamma_update = psi_gamma_update[(s-1), ],
                                                                    tausq_gamma_update = tausq_gamma_update[(s-1)])
    lambdasq_gamma_update[s, ] <- lambdasq_gamma_update_s


    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_HHSSI_MVN_cov2(J, M, K,
                                                                    delta_update = delta_update[s, , ],
                                                                    lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                                    lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                                    psi_delta_update = psi_delta_update[(s-1), ],
                                                                    tausq_delta_update = tausq_delta_update[(s-1)])
    lambdasq_delta_update[s, ] <- lambdasq_delta_update_s



    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_HHSSI_MVN_cov2(J, K,
                                                            beta_update = beta_update[s, , ],
                                                            lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                            xi_beta_update = xi_beta_update[(s-1)])
    tausq_beta_update[s] <- tausq_beta_update_s


    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_HHSSI_MVN_cov2(M, K,
                                                              gamma_update = gamma_update[s, , ],
                                                              lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                              xi_gamma_update = xi_gamma_update[(s-1)])
    tausq_gamma_update[s] <- tausq_gamma_update_s


    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_HHSSI_MVN_cov2(J, M, K,
                                                              delta_update = delta_update[s, , ],
                                                              lambdasq_delta_update = lambdasq_delta_update[s, ],
                                                              lambdasq_beta_update = lambdasq_beta_update[s, ],
                                                              lambdasq_gamma_update = lambdasq_gamma_update[s, ],
                                                              xi_delta_update = xi_delta_update[(s-1)])
    tausq_delta_update[s] <- tausq_delta_update_s



    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_HHSSI_MVN_cov2(J,
                                                        lambdasq_beta_update = lambdasq_beta_update[s, ])
    psi_beta_update[s, ] <- psi_beta_update_s


    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_HHSSI_MVN_cov2(M,
                                                          lambdasq_gamma_update = lambdasq_gamma_update[s, ])
    psi_gamma_update[s, ] <- psi_gamma_update_s


    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_HHSSI_MVN_cov2(J, M,
                                                          lambdasq_delta_update = lambdasq_delta_update[s, ])
    psi_delta_update[s, ] <- psi_delta_update_s



    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_HHSSI_MVN_cov(J, tausq_beta_update = tausq_beta_update[s])
    xi_beta_update[s] <- xi_beta_update_s


    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_HHSSI_MVN_cov(M, tausq_gamma_update = tausq_gamma_update[s])
    xi_gamma_update[s] <- xi_gamma_update_s


    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_HHSSI_MVN_cov(J, M, tausq_delta_update = tausq_delta_update[s])
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

