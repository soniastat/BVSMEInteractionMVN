
# devtools::document()
# Convergence check for the proposed Horseshoe model
###############################################################

library(LaplacesDemon) # IW

library(coda) # convergence check
###############################################################

# J = 20, M = 5
  set.seed(23444)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 20, M = 5)

###################################

Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
W <- cbind(1, x, z, u)
K <- dim(Y)[2] # Number of groups
n <- dim(Y)[1]
J <- dim(x)[2]
M <- dim(z)[2]
n_all_par <- dim(W)[2] # number of theta_k terms 1 + J + M + n_E + n_G

Psi_0 <- diag(K)          # K x K
nu_0 <- K + 2 # nu_0 should be greater that K+1
Sigma_init <- rinvwishart(nu_0, Psi_0)
###########################################
# BVSHHS_SI_MVN Model
###########################################

start_time1 <- Sys.time()
set.seed(3567)  # Set seed for reproducibility
# Run MCMC
res_all_par_chain1 <- fit_BVSHHS_SI_MVN(
  niter = 11000, burn_in = 1000, thin = 5,
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
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0
)
end_time1 <- Sys.time()
end_time1 - start_time1



start_time2 <- Sys.time()
set.seed(67832)  # Set seed for reproducibility
# Run MCMC
res_all_par_chain2 <- fit_BVSHHS_SI_MVN(
  niter = 11000, burn_in = 1000, thin = 5,
  n, K, Y, W, n_all_par, J, M,
  theta_init = matrix(0.1, nrow = n_all_par, ncol = K),
  lambdasq_beta_init = rep(0.3, J),
  tausq_beta_init = 10,
  lambdasq_gamma_init = rep(0.45, M),
  tausq_gamma_init = 1,
  lambdasq_delta_init = rep(0.35, J*M),
  tausq_delta_init = .1,
  psi_beta_init = rep(0.55, J),
  psi_gamma_init = rep(0.5, M),
  psi_delta_init = rep(0.5, J*M),
  xi_beta_init = 2, xi_gamma_init = 1,
  xi_delta_init = 4,
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0
)
end_time2 <- Sys.time()
end_time2 - start_time2


start_time3 <- Sys.time()
set.seed(2432)  # Set seed for reproducibility
# Run MCMC
res_all_par_chain3 <- fit_BVSHHS_SI_MVN(
  niter = 11000, burn_in = 1000, thin = 5,
  n, K, Y, W, n_all_par, J, M,
  theta_init = matrix(-0.6, nrow = n_all_par, ncol = K),
  lambdasq_beta_init = rep(0.3, J),
  tausq_beta_init = 3,
  lambdasq_gamma_init = rep(0.5, M),
  tausq_gamma_init = .5,
  lambdasq_delta_init = rep(0.25, J*M),
  tausq_delta_init = 2,
  psi_beta_init = rep(0.5, J),
  psi_gamma_init = rep(0.5, M),
  psi_delta_init = rep(0.54, J*M),
  xi_beta_init = 3, xi_gamma_init = 1,
  xi_delta_init = 2,
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0
)
end_time3 <- Sys.time()
end_time3 - start_time3


####################################################################################

# alpha0
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, 1, ], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, 1, ], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, 1, ], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)


# All convergence indices are good for alpha_0
#####################################################################################

# beta
# array dim = niter, n_all_par, K
# So here all beta for K = 1
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 1][, 2:(1 + J)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 1][, 2:(1 + J)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 1][, 2:(1 + J)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)

##########
# here all beta for K = 2
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 2][, 2:(1 + J)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 2][, 2:(1 + J)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 2][, 2:(1 + J)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)

##########
# here all beta for K = 3
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 3][, 2:(1 + J)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 3][, 2:(1 + J)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 3][, 2:(1 + J)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)

##########
# here all beta for K = 4
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 4][, 2:(1 + J)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 4][, 2:(1 + J)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 4][, 2:(1 + J)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)

##########
# here all beta for K = 5
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 5][, 2:(1 + J)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 5][, 2:(1 + J)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 5][, 2:(1 + J)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
############################################
# All convergence indices are good for beta
#####################################################################################

# gamma
# here all gamma for K = 1
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 1][, (1 + J + 1):(1 + J + M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 1][, (1 + J + 1):(1 + J + M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 1][, (1 + J + 1):(1 + J + M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
###########################
# here all gamma for K = 2
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 2][, (1 + J + 1):(1 + J + M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 2][, (1 + J + 1):(1 + J + M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 2][, (1 + J + 1):(1 + J + M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
###############################
# here all gamma for K = 3
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 3][, (1 + J + 1):(1 + J + M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 3][, (1 + J + 1):(1 + J + M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 3][, (1 + J + 1):(1 + J + M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
#############################
# here all gamma for K = 4
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 4][, (1 + J + 1):(1 + J + M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 4][, (1 + J + 1):(1 + J + M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 4][, (1 + J + 1):(1 + J + M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
###########################
# here all gamma for K = 5
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 5][, (1 + J + 1):(1 + J + M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 5][, (1 + J + 1):(1 + J + M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 5][, (1 + J + 1):(1 + J + M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
######################################################################
# All converges indices are good for gamma
#####################################################################################

# delta
# here all delta for K = 1
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 1][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 1][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 1][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
#################################
# here all delta for K = 2
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 2][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 2][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 2][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)

#################################
# here all delta for K = 3
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 3][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 3][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 3][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
#################################
# here all delta for K = 4
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 4][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 4][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 4][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
#################################
# here all delta for K = 5
chain1 <- mcmc(as.matrix(res_all_par_chain1$theta_update[, , 5][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain2 <- mcmc(as.matrix(res_all_par_chain2$theta_update[, , 5][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
chain3 <- mcmc(as.matrix(res_all_par_chain3$theta_update[, , 5][, (1 + J + M + 1):(1 + J + M + J*M)], type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for delta
#####################################################################################

# lamdasq beta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$lambdasq_beta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$lambdasq_beta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$lambdasq_beta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for lambdasq_beta
#####################################################################################

# lamdasq gamma
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$lambdasq_gamma_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$lambdasq_gamma_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$lambdasq_gamma_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for lambdasq_gamma
#####################################################################################

# lamdasq delta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$lambdasq_delta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$lambdasq_delta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$lambdasq_delta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for lambdasq_delta
#####################################################################################

# tausq beta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$tausq_beta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$tausq_beta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$tausq_beta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for tausq_beta
#####################################################################################

# tausq gamma
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$tausq_gamma_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$tausq_gamma_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$tausq_gamma_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for tausq_gamma
#####################################################################################

# tausq delta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$tausq_delta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$tausq_delta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$tausq_delta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for tausq_delta
#####################################################################################

# psi beta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$psi_beta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$psi_beta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$psi_beta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for psi_beta
#####################################################################################

# psi gamma
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$psi_gamma_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$psi_gamma_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$psi_gamma_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for psi_gamma
#####################################################################################

# psi delta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$psi_delta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$psi_delta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$psi_delta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for psi_delta
#####################################################################################

# xi beta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$xi_beta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$xi_beta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$xi_beta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for xi_beta
#####################################################################################

# xi gamma
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$xi_gamma_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$xi_gamma_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$xi_gamma_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for xi_gamma
#####################################################################################

# xi delta
chain1 <- mcmc(as.matrix(log(res_all_par_chain1$xi_delta_update), type = "l"))
chain2 <- mcmc(as.matrix(log(res_all_par_chain2$xi_delta_update), type = "l"))
chain3 <- mcmc(as.matrix(log(res_all_par_chain3$xi_delta_update), type = "l"))
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for xi_delta
#####################################################################################
# Sigma convergence test
#####################################################################################

# check symmetry & PD (all passed)
apply(res_all_par_chain1$Sigma_update, 1, function(S) {
  isSymmetric(S) && all(eigen(S, symmetric = TRUE)$values > 0)
})

###########
# What you want to see
# ✔ Flat “hairy caterpillars”
# ✔ No drift
# ✔ No sudden jumps to zero
#
# Red flags
# Downward drift → over-shrinkage
# Upward drift → θ not identified
# Eigenvalues ≈ 0 → numerical instability

Sigma_store <- res_all_par_chain1$Sigma_update
nsave <- dim(Sigma_store)[1]

eigvals <- matrix(NA, nsave, K)
for (s in 1:nsave) {
  eigvals[s, ] <- eigen(Sigma_store[s, , ], symmetric = TRUE)$values
}

matplot(eigvals, type = "l", lty = 1,
        xlab = "Iteration", ylab = "Eigenvalues",
        main = "Trace plots of eigenvalues of Sigma")


###########
# Should mix smoothly
# If flat at zero → over-shrinkage
# If exploding → weak prior / collinearity
sd_store <- matrix(NA, nsave, K)
for (k in 1:K) {
  sd_store[, k] <- sqrt(Sigma_store[, k, k])
}

matplot(sd_store, type = "l", lty = 1,
        xlab = "Iteration", ylab = "SD",
        main = "Trace plots of marginal SDs")
#################
cor_store <- array(NA, dim = c(nsave, K, K))
for (s in 1:nsave) {
  cor_store[s, , ] <- cov2cor(Sigma_store[s, , ])
}

plot(cor_store[, 1, 2], type = "l",
     ylab = "Corr(1,2)", xlab = "Iteration",
     main = "Trace plot of correlation (1,2)")
abline(h = 0, lty = 2)

##########
logdet <- sapply(1:nsave, function(s)
  determinant(Sigma_store[s, , ], logarithm = TRUE)$modulus
)

plot(logdet, type = "l",
     xlab = "Iteration", ylab = "log|Sigma|",
     main = "Log-determinant trace")
############

# Multiple chains
Sigma_store_chain1 <- res_all_par_chain1$Sigma_update
nsave_chain1 <- dim(Sigma_store_chain1)[1]
eigvals_chain1 <- matrix(NA, nsave_chain1, K)
for (s in 1:nsave_chain1) {
  eigvals_chain1[s, ] <- eigen(Sigma_store_chain1[s, , ], symmetric = TRUE)$values
}
sd_store_chain1 <- matrix(NA, nsave_chain1, K)
for (k in 1:K) {
  sd_store_chain1[, k] <- sqrt(Sigma_store_chain1[, k, k])
}


Sigma_store_chain2 <- res_all_par_chain2$Sigma_update
nsave_chain2 <- dim(Sigma_store_chain2)[1]
eigvals_chain2 <- matrix(NA, nsave_chain2, K)
for (s in 1:nsave_chain2) {
  eigvals_chain2[s, ] <- eigen(Sigma_store_chain2[s, , ], symmetric = TRUE)$values
}
sd_store_chain2 <- matrix(NA, nsave_chain2, K)
for (k in 1:K) {
  sd_store_chain2[, k] <- sqrt(Sigma_store_chain2[, k, k])
}


Sigma_store_chain3 <- res_all_par_chain3$Sigma_update
nsave_chain3 <- dim(Sigma_store_chain3)[1]
eigvals_chain3 <- matrix(NA, nsave_chain3, K)
for (s in 1:nsave_chain3) {
  eigvals_chain3[s, ] <- eigen(Sigma_store_chain3[s, , ], symmetric = TRUE)$values
}
sd_store_chain3 <- matrix(NA, nsave_chain3, K)
for (k in 1:K) {
  sd_store_chain3[, k] <- sqrt(Sigma_store_chain3[, k, k])
}

chain1 <- mcmc(sd_store_chain1)
chain2 <- mcmc(sd_store_chain2)
chain3 <- mcmc(sd_store_chain3)
traceplot(chain1)

chain_list <- mcmc.list(chain1, chain2, chain3)
traceplot(chain_list)

# Geweke diagnostics
geweke.diag(chain_list)

# Gelman-Rubin diagnostics
gelman.diag(chain_list , autoburnin = FALSE)
gelman.diag(chain_list , autoburnin = FALSE)$psrf[, 1]

# Autocorrelation plot
autocorr.plot(chain1 , auto.layout = FALSE)

# Effective sample sizes
effectiveSize(chain1)
##########################################
# All convergence indices are good for Sigma
################################################
# All are  good
