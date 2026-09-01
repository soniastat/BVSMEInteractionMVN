
# Bayesian variable selection with hierarchical Horseshoe priors
# for multivariate normal response (HHS-MVN-cov)
# seperate covariates


##############################################################
# Matrix
#' @importFrom Matrix bdiag Diagonal chol2inv

# MASS
#' @importFrom MASS mvrnorm

# mvtnorm
#' @importFrom mvtnorm dmvnorm

# invwishart
#' @importFrom LaplacesDemon rinvwishart

# folded normal
#' @importFrom greybox rfnorm

NULL
########################################################

# Update theta
update_theta_HHS_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par,
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
update_Sigma_HHS_MVN_cov <- function(Y, n, W, theta_update, Psi_0, nu_0)
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
update_lambdasq_beta_HHS_MVN_cov2 <- function(beta_update,
                                              psi_beta_update,
                                              tausq_beta_update)
{
  scale_lambdasq_beta <- (1 / psi_beta_update) +
    beta_update^2 / (2 * tausq_beta_update)

  lambdasq_beta_update_s <- matrix(rinvgamma(length(scale_lambdasq_beta),
                                             shape = 1,
                                             scale = c(t(scale_lambdasq_beta))),
                                   nrow = nrow(beta_update),
                                   byrow = TRUE
  )

  return(lambdasq_beta_update_s)
}



#######################################################################################

# Update lambdasq_gamma (for all m = 1, 2, ..., M, and k = 1, 2, ..., K)
update_lambdasq_gamma_HHS_MVN_cov2 <- function(gamma_update,
                                               psi_gamma_update,
                                               tausq_gamma_update)
{
  scale_lambdasq_gamma <-
    1 / psi_gamma_update +
    gamma_update^2 / (2 * tausq_gamma_update)

  lambdasq_gamma_update_s <- matrix(
    rinvgamma(length(scale_lambdasq_gamma),
              shape = 1,
              scale = c(t(scale_lambdasq_gamma))),
    nrow = nrow(gamma_update),
    byrow = TRUE
  )

   return(lambdasq_gamma_update_s)
}


#####################################################################################

# Update lambdasq_delta (for all j and m, and k = 1, 2, ..., K)
update_lambdasq_delta_HHS_MVN_cov2 <- function(J, M, K,
                                               delta_update,
                                               lambdasq_beta_update,
                                               lambdasq_gamma_update,
                                               psi_delta_update,
                                               tausq_delta_update)
{
  beta_part <- lambdasq_beta_update[rep(1:J, each = M), , drop = FALSE]

  gamma_part <- lambdasq_gamma_update[rep(1:M, times = J), , drop = FALSE]

  scale_lambdasq_delta <-
    1 / psi_delta_update +
    delta_update^2 /
    (2 * beta_part * gamma_part *
       tausq_delta_update)

  lambdasq_delta_update_s <- matrix(
    rinvgamma(length(scale_lambdasq_delta),
              shape = 1,
              scale = as.vector(t(scale_lambdasq_delta))),
    nrow = J*M,
    byrow = TRUE
  )

  return(lambdasq_delta_update_s)
}


#######################################################################################

# Update tausq_beta (for all j = 1, 2, ..., J)
update_tausq_beta_HHS_MVN_cov2 <- function(
    J, K,
    beta_update,
    lambdasq_beta_update,
    xi_beta_update)
{
  shape_tausq_beta <- (K + 1) / 2

  scale_tausq_beta <- (1 / xi_beta_update) +
    rowSums(beta_update^2 / (2 * lambdasq_beta_update))

  tausq_beta_update_s <- rinvgamma(J,
                                   shape = shape_tausq_beta,
                                   scale = scale_tausq_beta)

  return(tausq_beta_update_s)
}


#######################################################################################

# Update tausq_gamma (for all m = 1, 2, ..., M)
update_tausq_gamma_HHS_MVN_cov2 <- function(M, K, gamma_update,
                                            lambdasq_gamma_update,
                                            xi_gamma_update)
{
  shape_tausq_gamma <- (K + 1) / 2

  scale_tausq_gamma <- (1 / xi_gamma_update) +
    rowSums(gamma_update^2 / (2 * lambdasq_gamma_update))

  tausq_gamma_update_s <- rinvgamma(M,
                                    shape = shape_tausq_gamma,
                                    scale = scale_tausq_gamma)

   return(tausq_gamma_update_s)
}


#######################################################################################

# Update tausq_delta (for all j, and m )
update_tausq_delta_HHS_MVN_cov2 <- function(J, M, K,
                                            delta_update,
                                            lambdasq_delta_update,
                                            lambdasq_beta_update,
                                            lambdasq_gamma_update,
                                            xi_delta_update)
{
  shape_tausq_delta <- (K + 1) / 2

  ## Expand beta and gamma according to jm = (j-1)M + m
  beta_part  <- lambdasq_beta_update[rep(seq_len(J), each = M),
                                     , drop = FALSE]

  gamma_part <- lambdasq_gamma_update[rep(seq_len(M), times = J),
                                      , drop = FALSE]

  ## Compute all K terms simultaneously
  sum_term_tausq_delta <- rowSums(
    delta_update^2 /
      (2 * lambdasq_delta_update *
         beta_part *
         gamma_part)
  )

  scale_tausq_delta <- 1 / xi_delta_update + sum_term_tausq_delta

  tausq_delta_update_s <- rinvgamma(J * M,
                                    shape = shape_tausq_delta,
                                    scale = scale_tausq_delta)

  return(tausq_delta_update_s)
}



#######################################################################################

# Update psi_beta (for all j = 1, 2, ..., J, and k = 1, 2, ..., K)
update_psi_beta_HHS_MVN_cov2 <- function(lambdasq_beta_update)
{
  scale_psi_beta <- 1 + 1 / lambdasq_beta_update

  psi_beta_update_s <- matrix(
    rinvgamma(length(scale_psi_beta),
              shape = 1,
              scale = as.vector(t(scale_psi_beta))),
    nrow = nrow(lambdasq_beta_update),
    byrow = TRUE
  )

  return(psi_beta_update_s)
}

#######################################################################################

# Update psi_gamma (for all m = 1, 2, ..., M, and k = 1, 2, ..., K)
update_psi_gamma_HHS_MVN_cov2 <- function(lambdasq_gamma_update)
{
  scale_psi_gamma <- 1 + 1 / lambdasq_gamma_update

  psi_gamma_update_s <- matrix(
    rinvgamma(length(scale_psi_gamma),
              shape = 1,
              scale = as.vector(t(scale_psi_gamma))),
    nrow = nrow(lambdasq_gamma_update),
    byrow = TRUE
  )
  return(psi_gamma_update_s)
}



#######################################################################################

# Update psi_delta (for all j , m, and k = 1, 2, ..., K)
update_psi_delta_HHS_MVN_cov2 <- function(lambdasq_delta_update)
{
  scale_psi_delta <- 1 + 1 / lambdasq_delta_update

  psi_delta_update_s <- matrix(
    rinvgamma(length(scale_psi_delta),
              shape = 1,
              scale = as.vector(t(scale_psi_delta))),
    nrow = nrow(lambdasq_delta_update),
    byrow = TRUE
  )

  return(psi_delta_update_s)
}

#######################################################################################

# Update xi_beta (for j = 1, 2, ..., J)
update_xi_beta_HHS_MVN_cov <- function(J, tausq_beta_update)
{
  # Use more informative prior (xi_beta ~ IG(15, 3))
  shape_xi_beta <- ((J/2) + 15)
  scale_xi_beta <- (3 + sum(1 / tausq_beta_update))

  xi_beta_update_s <- rinvgamma(1, shape = shape_xi_beta, scale = scale_xi_beta) # IG distribution

  return(xi_beta_update_s)
}



#######################################################################################

# Update xi_gamma (for m = 1, 2, ..., M)
update_xi_gamma_HHS_MVN_cov <- function(M, tausq_gamma_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_gamma <- ((M/2) + 15)
  scale_xi_gamma <- (3 + sum(1 / tausq_gamma_update))

  xi_gamma_update_s <- rinvgamma(1, shape = shape_xi_gamma, scale = scale_xi_gamma) # IG distribution

   return(xi_gamma_update_s)
}



#######################################################################################

# Update xi_delta (for j and m)
update_xi_delta_HHS_MVN_cov <- function(J, M, tausq_delta_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_delta <- ((J*M/2) + 15)
  scale_xi_delta<- (3 + sum(1 / tausq_delta_update))

  xi_delta_update_s <- rinvgamma(1, shape = shape_xi_delta, scale = scale_xi_delta) # IG distribution

  return(xi_delta_update_s)
}

##################################################################################


#' Update parameters for HHS-MVN-cov model
#' @export
fit_HHS_MVN_cov_modify3 <- function(niter = 6000, burn_in = 1000, thin = 5,
                                    n, K, Y, W, WTW, n_all_par, J, M, O,
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
                                    Sigma_init = diag(K),
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
    theta_update_s <- update_theta_HHS_MVN_cov_modify(Y, K, W, WTW, n_all_par,
                                                      J, M, O,
                                                      sigmasq_alpha = sigmasq_alpha,
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
    Sigma_update_s <- update_Sigma_HHS_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s


    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_HHS_MVN_cov2(beta_update = beta_update[s, , ],
                                                                psi_beta_update = psi_beta_update[(s-1), , ],
                                                                tausq_beta_update = tausq_beta_update[(s-1), ])
    lambdasq_beta_update[s, , ] <- lambdasq_beta_update_s


    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_HHS_MVN_cov2(gamma_update = gamma_update[s, , ],
                                                                  psi_gamma_update = psi_gamma_update[(s-1), , ],
                                                                  tausq_gamma_update = tausq_gamma_update[(s-1), ])
    lambdasq_gamma_update[s, , ] <- lambdasq_gamma_update_s

    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_HHS_MVN_cov2(J, M, K,
                                                                  delta_update = delta_update[s, , ],
                                                                  lambdasq_beta_update = lambdasq_beta_update[s, , ],
                                                                  lambdasq_gamma_update = lambdasq_gamma_update[s, , ],
                                                                  psi_delta_update = psi_delta_update[(s-1), , ],
                                                                  tausq_delta_update = tausq_delta_update[(s-1), ])
    lambdasq_delta_update[s, , ] <- lambdasq_delta_update_s


    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_HHS_MVN_cov2(J, K,
                                                          beta_update = beta_update[s, , ],
                                                          lambdasq_beta_update = lambdasq_beta_update[s, , ],
                                                          xi_beta_update = xi_beta_update[(s-1)])
    tausq_beta_update[s, ] <- tausq_beta_update_s


    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_HHS_MVN_cov2(M, K,
                                                            gamma_update = gamma_update[s, , ],
                                                            lambdasq_gamma_update = lambdasq_gamma_update[s, , ],
                                                            xi_gamma_update = xi_gamma_update[(s-1)])
    tausq_gamma_update[s, ] <- tausq_gamma_update_s


    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_HHS_MVN_cov2(J, M, K,
                                                            delta_update = delta_update[s, , ],
                                                            lambdasq_delta_update = lambdasq_delta_update[s, , ],
                                                            lambdasq_beta_update = lambdasq_beta_update[s, , ],
                                                            lambdasq_gamma_update = lambdasq_gamma_update[s, , ],
                                                            xi_delta_update = xi_delta_update[(s-1)])
    tausq_delta_update[s, ] <- tausq_delta_update_s


    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_HHS_MVN_cov2(lambdasq_beta_update = lambdasq_beta_update[s, , ])
    psi_beta_update[s, , ] <- psi_beta_update_s


    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_HHS_MVN_cov2(lambdasq_gamma_update = lambdasq_gamma_update[s, , ])
    psi_gamma_update[s, , ] <- psi_gamma_update_s


    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_HHS_MVN_cov2(lambdasq_delta_update = lambdasq_delta_update[s, , ])
    psi_delta_update[s, , ] <- psi_delta_update_s


    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_HHS_MVN_cov(J, tausq_beta_update = tausq_beta_update[s, ])
    xi_beta_update[s] <- xi_beta_update_s


    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_HHS_MVN_cov(M, tausq_gamma_update = tausq_gamma_update[s, ])
    xi_gamma_update[s] <- xi_gamma_update_s


    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_HHS_MVN_cov(J, M, tausq_delta_update = tausq_delta_update[s, ])
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

