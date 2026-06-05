

# Normal prior model for multivariate normal response (NP-MVN-cov)



##############################################################
# Matrix
#' @importFrom Matrix bdiag Diagonal chol2inv

# MASS
#' @importFrom MASS mvrnorm

# mvtnorm
#' @importFrom mvtnorm dmvnorm


# invwishart
#' @importFrom LaplacesDemon rinvwishart


NULL


####################################################################################

# Update theta
update_theta_NP_MVN_cov <- function(Y, K, W, n_all_par, J, M, O,
                                    sigmasq_varphi, Sigma_update)
{

  # Build V_theta (same for all k)
  var_alpha_0 <- 100 # dim 1*1
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)
  mat_var_varphi <- diag(vec_sigmasq_varphi)

  V_theta <- as.matrix(bdiag(var_alpha_0, var_beta, var_gamma,
                             var_delta, mat_var_varphi))

  # prior precision
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


  return(theta_update_s)
}


#######################################################################################

# Update Sigma
update_Sigma_NP_MVN_cov <- function(Y, n, W, theta_update, Psi_0, nu_0)
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

#' Update parameters for NP-MVN-cov model
#' @export
fit_NP_MVN_cov <- function(niter = 20, burn_in = 2, thin = 1,
                            n, K, Y, W, n_all_par, J, M, O,
                            theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                            Sigma_init = matrix(0.5, nrow = K, ncol = K),
                            nu_0, Psi_0,
                            sigmasq_varphi = 10)
{
  theta_update <- array(NA, dim = c(niter, n_all_par, K))
  Sigma_update <- array(NA, dim = c(niter, K, K))

  # Initialize
  theta_update[1, , ] <- theta_init
  Sigma_update[1, , ] <- Sigma_init

  for(s in 2:niter)
  {
    if (s %% 50 == 0) cat("Iteration:", s, "\n")
    # theta_update
    theta_update_s <- update_theta_NP_MVN_cov(Y, K, W, n_all_par, J, M, O,
                                              sigmasq_varphi = sigmasq_varphi,
                                              Sigma_update = Sigma_update[(s-1), , ]
    )
    theta_update[s, , ] <- theta_update_s


    # Sigma_update
    Sigma_update_s <- update_Sigma_NP_MVN_cov(Y, n, W, theta_update[s, , ], Psi_0, nu_0)
    Sigma_update[s, , ] <- Sigma_update_s
  }

  # Apply burn-in and thinning
  indices_to_save <- seq(from = burn_in + 1, to = niter, by = thin)

  # Save only the thinned samples
  results <- list(
    theta_update = theta_update[indices_to_save, , ],
    Sigma_update = Sigma_update[indices_to_save, , ]
  )
  return(results)
}


#######################################################
