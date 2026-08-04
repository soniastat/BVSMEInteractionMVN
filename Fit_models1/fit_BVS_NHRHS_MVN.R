
# Bayesian variable selection with non-hierarchical regularized Horseshoe priors
# for multivariate normal response (BVS-NHRHS-MVN)



##############################################################
# Matrix
#' @importFrom Matrix bdiag Diagonal chol2inv

# MASS
#' @importFrom MASS mvrnorm

# mvtnorm
#' @importFrom mvtnorm dmvnorm

# actuar
#' @importFrom actuar rinvgamma

# invwishart
#' @importFrom LaplacesDemon rinvwishart

# folded normal
#' @importFrom greybox rfnorm

NULL





####################################################################################

# Update theta
update_theta_BVS_NHRHS_MVN <- function(Y, K, W, n_all_par,
                         J, M,
                         c,
                         lambdasq_beta_update, tausq_beta_update,
                         lambdasq_gamma_update, tausq_gamma_update,
                         lambdasq_delta_update, tausq_delta_update,
                         Sigma_update
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

        sigmasq_delta_update[jm, k] <- (lambdasq_delta_update_tilde *
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

    V_theta_list[[k]] <- as.matrix(bdiag(mat_var_alpha_0k, var_beta_k, var_gamma_k,
                                         var_delta_k))
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
update_Sigma_BVS_NHRHS_MVN <- function(Y, n, W, theta_update, Psi_0, nu_0)
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
update_lambdasq_beta_BVS_NHRHS_MVN <- function(c, J, K, beta_update, psi_beta_update, tausq_beta_update,
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
update_lambdasq_gamma_BVS_NHRHS_MVN <- function(c, M, K, gamma_update, psi_gamma_update,
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
update_lambdasq_delta_BVS_NHRHS_MVN <- function(c, J, M, K, delta_update, psi_delta_update,
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
update_tausq_beta_BVS_NHRHS_MVN <- function(c, J, K, beta_update,
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
update_tausq_gamma_BVS_NHRHS_MVN <- function(c, M, K, gamma_update,
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
update_tausq_delta_BVS_NHRHS_MVN <- function(c, J, M, K, delta_update, xi_delta_update,
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
update_psi_beta_BVS_NHRHS_MVN <- function(J, K, lambdasq_beta_update)
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
update_psi_gamma_BVS_NHRHS_MVN <- function(M, K, lambdasq_gamma_update)
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
update_psi_delta_BVS_NHRHS_MVN <- function(J, M, K, lambdasq_delta_update)
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
update_xi_beta_BVS_NHRHS_MVN <- function(J, tausq_beta_update)
{
  # Use more informative prior (xi_beta ~ IG(15, 3))
  shape_xi_beta <- ((J/2) + 15)
  scale_xi_beta <- (3 + sum(1 / tausq_beta_update))

  xi_beta_update_s <- rinvgamma(1, shape = shape_xi_beta, scale = scale_xi_beta) # IG distribution
  return(xi_beta_update_s)
}



#######################################################################################

# Update xi_gamma
update_xi_gamma_BVS_NHRHS_MVN <- function(M, tausq_gamma_update)
{
  # Use more informative prior (xi_gamma ~ IG(15, 3))
  shape_xi_gamma <- ((M/2) + 15)
  scale_xi_gamma <- (3 + sum(1 / tausq_gamma_update))

  xi_gamma_update_s <- rinvgamma(1, shape = shape_xi_gamma, scale = scale_xi_gamma) # IG distribution
  return(xi_gamma_update_s)
}





#######################################################################################

# Update xi_delta
update_xi_delta_BVS_NHRHS_MVN <- function(J, M, tausq_delta_update)
{
  # Use more informative prior (xi_delta ~ IG(15, 3))
  shape_xi_delta <- ((J*M/2) + 15)
  scale_xi_delta <- (3 + sum(1 / tausq_delta_update))

  xi_delta_update_s <- rinvgamma(1, shape = shape_xi_delta, scale = scale_xi_delta) # IG distribution
  return(xi_delta_update_s)
}



##################################################################################


#' Update parameters for BVS-NHRHS-MVN model
#' @export
fit_BVS_NHRHS_MVN <- function(niter = 20, burn_in = 2, thin = 1,
                            n, K, Y, W, n_all_par, J, M,
                            c = 2.5,
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
                            accept_tausq_delta_init = rep(1, times = J*M))
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



  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_BVS_NHRHS_MVN(Y, K, W, n_all_par,
                                   J, M,
                                   c,
                                   lambdasq_beta_update = lambdasq_beta_update[(s-1), , ],
                                   tausq_beta_update = tausq_beta_update[(s-1), ],
                                   lambdasq_gamma_update = lambdasq_gamma_update[(s-1), , ],
                                   tausq_gamma_update = tausq_gamma_update[(s-1), ],
                                   lambdasq_delta_update = lambdasq_delta_update[(s-1), , ],
                                   tausq_delta_update = tausq_delta_update[(s-1), ],
                                   Sigma_update = Sigma_update[(s-1), , ]
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
    Sigma_update_s <- update_Sigma_BVS_NHRHS_MVN(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s


    # lambdasq_beta_update
    lambdasq_beta_update_s <- update_lambdasq_beta_BVS_NHRHS_MVN(c, J, K,
                                                   beta_update = beta_update[s, , ],
                                                   psi_beta_update = psi_beta_update[(s-1), , ],
                                                   tausq_beta_update = tausq_beta_update[(s-1), ],
                                                   lambdasq_beta_update_s_1 = lambdasq_beta_update[(s-1),,],
                                                   omegasq_beta_lambdasq)
    lambdasq_beta_update[s, , ] <- lambdasq_beta_update_s$lambdasq_beta_update
    accept_lambdasq_beta_update[s, ,] <- lambdasq_beta_update_s$accept_lambdasq_beta_update


    # lambdasq_gamma_update
    lambdasq_gamma_update_s <- update_lambdasq_gamma_BVS_NHRHS_MVN(c, M, K,
                                                     gamma_update = gamma_update[s, , ],
                                                     psi_gamma_update = psi_gamma_update[(s-1), , ],
                                                     tausq_gamma_update = tausq_gamma_update[(s-1), ],
                                                     lambdasq_gamma_update_s_1 = lambdasq_gamma_update[(s-1),,],
                                                     omegasq_gamma_lambdasq)
    lambdasq_gamma_update[s, , ] <- lambdasq_gamma_update_s$lambdasq_gamma_update
    accept_lambdasq_gamma_update[s, ,] <- lambdasq_gamma_update_s$accept_lambdasq_gamma_update

    # lambdasq_delta_update
    lambdasq_delta_update_s <- update_lambdasq_delta_BVS_NHRHS_MVN(c, J, M, K,
                                                     delta_update = delta_update[s, , ],
                                                     psi_delta_update = psi_delta_update[(s-1), , ],
                                                     tausq_delta_update = tausq_delta_update[(s-1), ],
                                                     lambdasq_delta_update_s_1 = lambdasq_delta_update[(s-1),,],
                                                     omegasq_delta_lambdasq)
    lambdasq_delta_update[s, , ] <- lambdasq_delta_update_s$lambdasq_delta_update
    accept_lambdasq_delta_update[s, ,] <- lambdasq_delta_update_s$accept_lambdasq_delta_update



    # Update tausq_beta
    tausq_beta_update_s <- update_tausq_beta_BVS_NHRHS_MVN(c, J, K,
                                             beta_update = beta_update[s, , ],
                                             lambdasq_beta_update = lambdasq_beta_update[s, , ],
                                             xi_beta_update = xi_beta_update[(s-1)],
                                             tausq_beta_update_s_1 = tausq_beta_update[(s-1), ],
                                             omegasq_beta_tausq)
    tausq_beta_update[s, ] <- tausq_beta_update_s$tausq_beta_update
    accept_tausq_beta_update[s,] <- tausq_beta_update_s$accept_tausq_beta_update


    # Update tausq_gamma
    tausq_gamma_update_s <- update_tausq_gamma_BVS_NHRHS_MVN(c, M, K,
                                               gamma_update = gamma_update[s, , ],
                                               lambdasq_gamma_update = lambdasq_gamma_update[s, , ],
                                               xi_gamma_update = xi_gamma_update[(s-1)],
                                               tausq_gamma_update_s_1 = tausq_gamma_update[(s-1), ],
                                               omegasq_gamma_tausq)
    tausq_gamma_update[s, ] <- tausq_gamma_update_s$tausq_gamma_update
    accept_tausq_gamma_update[s,] <- tausq_gamma_update_s$accept_tausq_gamma_update


    # Update tausq_delta
    tausq_delta_update_s <- update_tausq_delta_BVS_NHRHS_MVN(c, J, M, K,
                                               delta_update = delta_update[s, , ],
                                               xi_delta_update = xi_delta_update[(s-1)],
                                               lambdasq_delta_update = lambdasq_delta_update[s, , ],
                                               tausq_delta_update_s_1 = tausq_delta_update[(s-1),],
                                               omegasq_delta_tausq)
    tausq_delta_update[s, ] <- tausq_delta_update_s$tausq_delta_update
    accept_tausq_delta_update[s,] <- tausq_delta_update_s$accept_tausq_delta_update


    # Update psi_beta
    psi_beta_update_s <- update_psi_beta_BVS_NHRHS_MVN(J, K,
                                         lambdasq_beta_update = lambdasq_beta_update[s, , ])
    psi_beta_update[s, , ] <- psi_beta_update_s


    # Update psi_gamma
    psi_gamma_update_s <- update_psi_gamma_BVS_NHRHS_MVN(M, K,
                                           lambdasq_gamma_update = lambdasq_gamma_update[s, , ])
    psi_gamma_update[s, , ] <- psi_gamma_update_s


    # Update psi_delta
    psi_delta_update_s <- update_psi_delta_BVS_NHRHS_MVN(J, M, K,
                                           lambdasq_delta_update = lambdasq_delta_update[s, , ])
    psi_delta_update[s, , ] <- psi_delta_update_s



    # Update xi_beta
    xi_beta_update_s <- update_xi_beta_BVS_NHRHS_MVN(J, tausq_beta_update = tausq_beta_update[s, ])
    xi_beta_update[s] <- xi_beta_update_s


    # Update xi_gamma
    xi_gamma_update_s <- update_xi_gamma_BVS_NHRHS_MVN(M, tausq_gamma_update = tausq_gamma_update[s, ])
    xi_gamma_update[s] <- xi_gamma_update_s


    # Update xi_delta
    xi_delta_update_s <- update_xi_delta_BVS_NHRHS_MVN(J, M, tausq_delta_update = tausq_delta_update[s, ])
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
