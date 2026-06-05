
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

# Fit the NP-MVN-cov model
start_time_NP_MVN_cov <- Sys.time()
res_all_par_NP_MVN_cov <- fit_NP_MVN_cov(
        niter = 6000, burn_in = 1000, thin = 5,
        n, K, Y, W, n_all_par, J, M, O,
        theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
        Sigma_init = Sigma_init,
        nu_0 = nu_0, Psi_0 = Psi_0,
        sigmasq_varphi = 10
      )
end_time_NP_MVN_cov <- Sys.time()
end_time_NP_MVN_cov - start_time_NP_MVN_cov


# save(res_all_par_NP_MVN_cov, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K5_6000ite_1000burn_5thin_NPMVNcov_pregpostexpo.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K5_6000ite_1000burn_5thin_NPMVNcov_pregpostexpo.rda")

#######################################################

start_time_NP_MVN_cov_chain2 <- Sys.time()
set.seed(453567)
res_all_par_NP_MVN_cov_chain2 <- fit_NP_MVN_cov(
  niter = 6000, burn_in = 1000, thin = 5,
  n, K, Y, W, n_all_par, J, M, O,
  theta_init = matrix(0.15, nrow = n_all_par, ncol = K),
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  sigmasq_varphi = 10
)
end_time_NP_MVN_cov_chain2 <- Sys.time()
end_time_NP_MVN_cov_chain2 - start_time_NP_MVN_cov_chain2


# save(res_all_par_NP_MVN_cov_chain2, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVNcov_pregpostexpo_chain2.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVNcov_pregpostexpo_chain2.rda")


#################################################################################


start_time_NP_MVN_cov_chain3 <- Sys.time()
set.seed(98765453)
res_all_par_NP_MVN_cov_chain3 <- fit_NP_MVN_cov(
  niter = 6000, burn_in = 1000, thin = 5,
  n, K, Y, W, n_all_par, J, M, O,
  theta_init = matrix(2.5, nrow = n_all_par, ncol = K),
  Sigma_init = Sigma_init,
  nu_0 = nu_0, Psi_0 = Psi_0,
  sigmasq_varphi = 10
)
end_time_NP_MVN_cov_chain3 <- Sys.time()
end_time_NP_MVN_cov_chain3 - start_time_NP_MVN_cov_chain3


# save(res_all_par_NP_MVN_cov_chain3, file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVNcov_pregpostexpo_chain3.rda")
# load(file = "exposome_data_analysis/res_n1301_O6_M34_J44_K4_6000ite_1000burn_5thin_NPMVNcov_pregpostexpo_chain3.rda")


#################################################################################






