
# Testing codes for BVS-HHS-MVN-cov

sample_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat("Running sample:", sample_id, "\n")

###################################################

source("/projects/soniast@colostate.edu/alpine_sim_HHS.R")



library(LaplacesDemon)
library(tidyverse)
library(readxl)


############################
# Load data
############################

load("/projects/soniast@colostate.edu/exposome.RData")
load("/projects/soniast@colostate.edu/modifiers_covariates.rda")
load("/projects/soniast@colostate.edu/exposures.rda")
load("/projects/soniast@colostate.edu/cov_base_post.rda")



raw_Y <- cbind(
  phenotype[,5:6],
  cov_base_post
)

Y0 <- scale(raw_Y)

Z <- cbind(
  Organochlorines_Postnatal,
  Metals_Postnatal
)

K <- dim(Y0)[2]
n <- nrow(Y0)
J <- ncol(X)
M <- ncol(Z)
O <- ncol(D)

############################
# Simulate one dataset
############################

set.seed(23444 + sample_id)

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

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_HHS_MVN_cov(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W,
    n_all_par = n_all_par, J = J, M = M, O = O,
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
    sigmasq_varphi = 10
  )

  res$theta_update

}, error=function(e){
  list(
    error=TRUE,
    message=e$message,
    sample_id=sample_id
  )

})




dir.create("/projects/soniast@colostate.edu/results_HHS", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_HHS/HHS_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished sample", sample_id, "\n")






