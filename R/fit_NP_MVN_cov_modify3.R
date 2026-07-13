

# Normal prior model for multivariate normal response (NP-MVN-cov)



##############################################################
# Matrix
#' @importFrom Matrix bdiag Diagonal chol2inv


# invwishart
#' @importFrom LaplacesDemon rinvwishart


NULL


####################################################################################

# Update theta
update_theta_NP_MVN_cov_modify <- function(Y, K, W, WTW, n_all_par, J, M, O,
                                           sigmasq_alpha = 100,
                                           sigmasq_varphi, Sigma_update)
{

  # Build V_theta (same for all k)
  var_beta <- diag(rep(1, J)) # dim J*J
  var_gamma <- diag(rep(1, M)) # dim M*M
  var_delta <- diag(rep(1, J*M)) # dim JM*JM

  vec_sigmasq_varphi <- rep(sigmasq_varphi, times = O)

  # Inverse using blocks separately
  V_theta_inv <- as.matrix(
    bdiag(
      1/sigmasq_alpha,
      diag(1/diag(var_beta)),
      diag(1/diag(var_gamma)),
      diag(1/diag(var_delta)),
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
fit_NP_MVN_cov_modify3 <- function(niter = 6000, burn_in = 1000, thin = 5,
                                   n, K, Y, W, WTW, n_all_par, J, M, O,
                                   sigmasq_alpha = 100,
                                   theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
                                   Sigma_init = diag(K),
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
    theta_update_s <- update_theta_NP_MVN_cov_modify(Y, K, W, WTW, n_all_par, J, M, O,
                                                     sigmasq_alpha = 100,
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







