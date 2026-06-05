

# Testing codes for Bayesian variable selection with hierarchical regularized
# Horseshoe priors for multivariate normal response (BVS-HRHS-MVN-cov)

rm(list = ls())
gc()

###############################

# Run all the functions from the package
devtools::load_all()

library(doParallel) # for parallelization
library(foreach) # for parallelization
library(LaplacesDemon) # IW

library(doRNG) # used for instead of set.seed

library(ggplot2) # plot

#######################################################

# real data
options(echo=TRUE)
options(stringsAsFactors = FALSE)


library(tidyverse)
library(readxl)


load("exposome_data_analysis/exposome.RData")
load("exposome_data_analysis/modifiers_covariates.rda")
load("exposome_data_analysis/exposures.rda")

names(covariates)


names(phenotype)

covariates_post <- codebook %>%
  filter(domain=="Covariates" & period == "Postnatal") %>%
  pull(variable_name) %>%
  as.character()
covariates_post <- covariates_post[covariates_post %in%
                                     c("hs_c_height_None", "hs_c_weight_None")]
cov_base_post <- covariates[,covariates_post]


# Combine responses/covariates
raw_Y <- cbind(phenotype[, 5:6], cov_base_post)
Y <- scale(raw_Y)

Z <- cbind(Organochlorines_Postnatal, Metals_Postnatal) # only postnatal

J <- dim(X)[2]
M <- dim(Z)[2]
O <- dim(D)[2]
K <- dim(Y)[2]
n <- dim(Y)[1]

# Modifier-environment interaction terms
U <- matrix(0, n, J * M)
col_idx <- 1
for (j in 1:J) {
  for (m in 1:M) {
    U[, col_idx] <- X[, j] * Z[, m]
    col_idx <- col_idx + 1
  }
}

# modifiers = X, covariates = D

W <- cbind(1, X, Z, U, D)
n_all_par <- dim(W)[2]


######################################################
# Fit the BVS-HRHS-MVN-cov model

start_time_BVS_HRHS_MVN_cov <- Sys.time()
res_all_par_BVS_HRHS_MVN_cov <- fit_BVS_HRHS_MVN_cov(
  niter = 6000, burn_in = 1000, thin = 5,
  n, K, Y, W, n_all_par, J, M, O,
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
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  omegasq_beta_lambdasq = 0.25, omegasq_beta_tausq = 0.25,
  omegasq_gamma_lambdasq = 0.25, omegasq_gamma_tausq = 0.25,
  omegasq_delta_lambdasq = 0.25, omegasq_delta_tausq = 0.25,
  accept_lambdasq_beta_init = matrix(1, nrow = J, ncol = K),
  accept_tausq_beta_init = rep(1, times = J),
  accept_lambdasq_gamma_init = matrix(1, nrow = M, ncol = K),
  accept_tausq_gamma_init = rep(1, times = M),
  accept_lambdasq_delta_init = matrix(1, nrow = J*M, ncol = K),
  accept_tausq_delta_init = rep(1, times = M),
  sigmasq_varphi = 10
)
end_time_BVS_HRHS_MVN_cov <- Sys.time()
end_time_BVS_HRHS_MVN_cov - start_time_BVS_HRHS_MVN_cov

# save(res_all_par_BVS_HRHS_MVN_cov, file = "exposome_data_analysis/res_n1301_O6_M17_J44_K5_6000ite_1000burn_5thin_BVSHRHSMVNcov_postexpo.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M17_J44_K5_6000ite_1000burn_5thin_BVSHRHSMVNcov_postexpo.rda")


#######################################################

start_time_BVS_HRHS_MVN_cov_chain2 <- Sys.time()
set.seed(453567)
res_all_par_BVS_HRHS_MVN_cov_chain2 <- fit_BVS_HRHS_MVN_cov(
  niter = 6000, burn_in = 1000, thin = 5,
  n = n, K = K, Y = Y, W = W, n_all_par = n_all_par, J = J, M = M, O,
  c = 2.5,
  theta_init = matrix(-0.5, nrow = n_all_par, ncol = K),
  lambdasq_beta_init = matrix(0.15, nrow = J, ncol = K),
  tausq_beta_init = rep(5, J),
  lambdasq_gamma_init = matrix(0.15, nrow = M, ncol = K),
  tausq_gamma_init = rep(5, M),
  lambdasq_delta_init = matrix(0.15, nrow = J*M, ncol = K),
  tausq_delta_init = rep(5, J*M),
  psi_beta_init = matrix(0.15, nrow = J, ncol = K),
  psi_gamma_init = matrix(0.15, nrow = M, ncol = K),
  psi_delta_init = matrix(0.15, nrow = J*M, ncol = K),
  xi_beta_init = 5, xi_gamma_init = 5,
  xi_delta_init = 5,
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  omegasq_beta_lambdasq = 0.25, omegasq_beta_tausq = 0.25,
  omegasq_gamma_lambdasq = 0.25, omegasq_gamma_tausq = 0.25,
  omegasq_delta_lambdasq = 0.25, omegasq_delta_tausq = 0.25,
  accept_lambdasq_beta_init = matrix(5, nrow = J, ncol = K),
  accept_tausq_beta_init = rep(5, times = J),
  accept_lambdasq_gamma_init = matrix(5, nrow = M, ncol = K),
  accept_tausq_gamma_init = rep(5, times = M),
  accept_lambdasq_delta_init = matrix(5, nrow = J*M, ncol = K),
  accept_tausq_delta_init = rep(5, times = M),
  sigmasq_varphi = 10
)
end_time_BVS_HRHS_MVN_cov_chain2 <- Sys.time()
end_time_BVS_HRHS_MVN_cov_chain2 - start_time_BVS_HRHS_MVN_cov_chain2


# save(res_all_par_BVS_HRHS_MVN_cov_chain2, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K5_6000ite_1000burn_5thin_BVSHRHSMVNcov_pregpostexpo_chain2.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K5_6000ite_1000burn_5thin_BVSHRHSMVNcov_pregpostexpo_chain2.rda")

#######################################################

start_time_BVS_HRHS_MVN_cov_chain3 <- Sys.time()
set.seed(98765453)
res_all_par_BVS_HRHS_MVN_cov_chain3 <- fit_BVS_HRHS_MVN_cov(
  niter = 6000, burn_in = 1000, thin = 5,
  n = n, K = K, Y = Y, W = W, n_all_par = n_all_par, J = J, M = M, O,
  c = 2.5,
  theta_init = matrix(5, nrow = n_all_par, ncol = K),
  lambdasq_beta_init = matrix(2.5, nrow = J, ncol = K),
  tausq_beta_init = rep(.1, J),
  lambdasq_gamma_init = matrix(2.5, nrow = M, ncol = K),
  tausq_gamma_init = rep(.1, M),
  lambdasq_delta_init = matrix(2.5, nrow = J*M, ncol = K),
  tausq_delta_init = rep(.1, J*M),
  psi_beta_init = matrix(2.5, nrow = J, ncol = K),
  psi_gamma_init = matrix(2.5, nrow = M, ncol = K),
  psi_delta_init = matrix(2.5, nrow = J*M, ncol = K),
  xi_beta_init = .1, xi_gamma_init = .1,
  xi_delta_init = .1,
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  omegasq_beta_lambdasq = 0.25, omegasq_beta_tausq = 0.25,
  omegasq_gamma_lambdasq = 0.25, omegasq_gamma_tausq = 0.25,
  omegasq_delta_lambdasq = 0.25, omegasq_delta_tausq = 0.25,
  accept_lambdasq_beta_init = matrix(.1, nrow = J, ncol = K),
  accept_tausq_beta_init = rep(.1, times = J),
  accept_lambdasq_gamma_init = matrix(1, nrow = M, ncol = K),
  accept_tausq_gamma_init = rep(.1, times = M),
  accept_lambdasq_delta_init = matrix(.1, nrow = J*M, ncol = K),
  accept_tausq_delta_init = rep(.1, times = M),
  sigmasq_varphi = 10
)
end_time_BVS_HRHS_MVN_cov_chain3 <- Sys.time()
end_time_BVS_HRHS_MVN_cov_chain3 - start_time_BVS_HRHS_MVN_cov_chain3


# save(res_all_par_BVS_HRHS_MVN_cov_chain3, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_BVSHRHSMVNcov_pregpostexpo_chain3.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_BVSHRHSMVNcov_pregpostexpo_chain3.rda")


#################################################################################


