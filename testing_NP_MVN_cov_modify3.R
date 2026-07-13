
# Testing codes for normal priors model for multivariate normal response (NP-MVN-cov)

rm(list = ls())
gc()

###############################

# Run all the functions from the package
devtools::load_all()

library(LaplacesDemon) # IW
library(ggplot2) # plot

#######################################################

# real data
options(echo=TRUE)
options(stringsAsFactors = FALSE)


library(tidyverse)
library(readxl)


load("exposome.RData")
load("modifiers_covariates.rda")
load("exposures.rda")
load("cov_base_post.rda")


raw_Y <- cbind(
  phenotype[,5:6],
  cov_base_post
)

Y0 <- scale(raw_Y)

Z <- cbind(
  Organochlorines_Postnatal,
  Metals_Postnatal
)


Sigma <- cov(Y0)
K <- dim(Y0)[2]
n <- nrow(Y0)
J <- ncol(X)
M <- ncol(Z)
O <- ncol(D)

############################
# Simulate one dataset
############################

set.seed(23444)

sim_data <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D,
                              Sigma=Sigma)

Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
d <- sim_data$d

W <- cbind(1, x, z, u, d)
WTW <- t(W) %*% W

n_all_par <- ncol(W)

nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)



######################################################

# Fit the NP-MVN-cov model
start_time_NP_MVN_cov <- Sys.time()
res_all_par_NP_MVN_cov <- fit_NP_MVN_cov_modify3(
  niter = 6000, burn_in = 1000, thin = 5,
  n=n, K=K, Y=Y, W=W, WTW=WTW, n_all_par=n_all_par,
  J=J, M=M, O=O,
  sigmasq_alpha = 100,
  theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  sigmasq_varphi = 10)
end_time_NP_MVN_cov <- Sys.time()
end_time_NP_MVN_cov - start_time_NP_MVN_cov


# save(res_all_par_NP_MVN_cov, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K5_6000ite_1000burn_5thin_NPMVN_modify3.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K5_6000ite_1000burn_5thin_NPMVN_modify3.rda")

#######################################################

start_time_NP_MVN_cov_chain2 <- Sys.time()
set.seed(453567)
res_all_par_NP_MVN_cov_chain2 <- fit_NP_MVN_cov(
  niter = 6000, burn_in = 1000, thin = 5,
  n=n, K=K, Y=Y, W=W, WTW=WTW, n_all_par=n_all_par,
  J=J, M=M, O=O,
  sigmasq_alpha = 100,
  theta_init = matrix(-5, nrow = n_all_par, ncol = K),
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  sigmasq_varphi = 10)
end_time_NP_MVN_cov_chain2 <- Sys.time()
end_time_NP_MVN_cov_chain2 - start_time_NP_MVN_cov_chain2


# save(res_all_par_NP_MVN_cov_chain2, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVN_modify3_chain2.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVN_modify3_chain2.rda")


#################################################################################


start_time_NP_MVN_cov_chain3 <- Sys.time()
set.seed(98765453)
res_all_par_NP_MVN_cov_chain3 <- fit_NP_MVN_cov(
  niter = 6000, burn_in = 1000, thin = 5,
  n=n, K=K, Y=Y, W=W, WTW=WTW, n_all_par=n_all_par,
  J=J, M=M, O=O,
  sigmasq_alpha = 100,
  theta_init = matrix(5, nrow = n_all_par, ncol = K),
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  sigmasq_varphi = 10)
end_time_NP_MVN_cov_chain3 <- Sys.time()
end_time_NP_MVN_cov_chain3 - start_time_NP_MVN_cov_chain3


# save(res_all_par_NP_MVN_cov_chain3, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVN_modify3_chain3.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVN_modify3_chain3.rda")


#################################################################################






