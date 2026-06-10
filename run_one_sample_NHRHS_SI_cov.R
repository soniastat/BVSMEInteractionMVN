
# Testing codes for BVS-NHRHS-SI-MVN-cov

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

sample <- as.numeric(args[1])

cat("Running sample:", sample, "\n")


suppressPackageStartupMessages({
  library(BVSMEInteractionMVN)
  library(LaplacesDemon)
  library(tidyverse)
  library(readxl)
})

############################
# Load data
############################

load("exposome.RData")
load("modifiers_covariates.rda")
load("exposures.rda")

covariates_post <- codebook %>%
  filter(domain=="Covariates" & period=="Postnatal") %>%
  pull(variable_name) %>%
  as.character()

covariates_post <- covariates_post[
  covariates_post %in%
    c("hs_c_height_None","hs_c_weight_None")
]

cov_base_post <- covariates[,covariates_post]

raw_Y <- cbind(
  phenotype[,5:6],
  cov_base_post
)

Y0 <- scale(raw_Y)

Z <- cbind(
  Organochlorines_Postnatal,
  Metals_Postnatal
)

K <- 5
n <- nrow(Y0)

J <- ncol(X)
M <- ncol(Z)
O <- ncol(D)

############################
# Simulate one dataset
############################

set.seed(23444 + sample)
if(is.na(sample)){
  stop("Sample number missing")
}

sim_data <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D)

Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
d <- sim_data$d

W <- cbind(1, x, z, u, d)

n_all_par <- ncol(W)

nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)

############################
# Fit model
############################

set.seed(765878 + sample)

theta_update_save <- tryCatch({

  res <- fit_BVS_NHRHS_SI_MVN_cov(
    niter = 6000, burn_in = 1000, thin = 5,
    n=n, K=K, Y=Y, W=W, n_all_par=n_all_par,
    J=J, M=M, O=O,
    c = 2.5,
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
    nu_0 = nu_0, Psi_0 = Psi_0,
    omegasq_beta_lambdasq = 0.25, omegasq_beta_tausq = 0.25,
    omegasq_gamma_lambdasq = 0.25, omegasq_gamma_tausq = 0.25,
    omegasq_delta_lambdasq = 0.25, omegasq_delta_tausq = 0.25,
    accept_lambdasq_beta_init = rep(1, times = J),
    accept_tausq_beta_init = 1,
    accept_lambdasq_gamma_init = rep(1, times = M),
    accept_tausq_gamma_init = 1,
    accept_lambdasq_delta_init = rep(1, times = J*M),
    accept_tausq_delta_init = 1,
    sigmasq_varphi = 10)

  res$theta_update

}, error=function(e){
  list(
    error=TRUE,
    message=e$message,
    sample=sample
  )

})

dir.create("results", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("results/NHRHS_SI_cov_sample_", sample, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished sample", sample, "\n")






# files <- list.files(
#   "results",
#   pattern = "\\.rds$",
#   full.names = TRUE
# )
#
# files <- files[!grepl("theta_update_list", files)]
#
#
#
# theta_update_list <- lapply(files, readRDS)
# length(theta_update_list)
#
#
# saveRDS(
#   theta_update_list,
#   file = "results/theta_update_NHRHS_SI_cov_list.rds"
# )
#
#
# # in terminal, run these
#
# # cd /Users/sonia/Desktop/Packages/BVSMEInteractionMVN
#
# # seq 1 100 | xargs -n 1 -P 5 Rscript run_one_sample_NHRHS_SI_cov.R
#
#
#
