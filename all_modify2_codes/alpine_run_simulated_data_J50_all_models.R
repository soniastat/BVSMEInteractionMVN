
# Testing codes

sample_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat("Running sample:", sample_id, "\n")

###################################################

source("/projects/soniast@colostate.edu/alpine_sim_simulated_data_all_models.R")



library(LaplacesDemon)
library(Matrix)
library(greybox)




############################
# Simulate one dataset
############################

K <- 5
n <- 300
J <- 50
M <- 5
O <- 2

set.seed(23444 + sample_id)
sim_data <- Sim_data_BMVS(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O
)



Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
d <- sim_data$d
W <- sim_data$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)




############################
# Fit HHS model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_HHS_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
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




dir.create("/projects/soniast@colostate.edu/results_HHS_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_HHS_simdata/J50_HHS_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished HHS sample", sample_id, "\n")

########################










############################
# Fit HHSSI model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_HHS_SI_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
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




dir.create("/projects/soniast@colostate.edu/results_HHSSI_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_HHSSI_simdata/J50_HHSSI_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished hhssi sample", sample_id, "\n")





############################
# Fit HRHS model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_HRHS_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
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

  res$theta_update

}, error=function(e){
  list(
    error=TRUE,
    message=e$message,
    sample_id=sample_id
  )

})




dir.create("/projects/soniast@colostate.edu/results_HRHS_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_HRHS_simdata/J50_HRHS_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished hrhs sample", sample_id, "\n")







############################
# Fit HRHSSI model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_HRHS_SI_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
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
    sample_id=sample_id
  )

})




dir.create("/projects/soniast@colostate.edu/results_HRHSSI_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_HRHSSI_simdata/J50_HRHSSI_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished hrhssi sample", sample_id, "\n")







############################
# Fit NHHS model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_NHHS_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par, J = J,
    M = M, O = O,
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




dir.create("/projects/soniast@colostate.edu/results_NHHS_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_NHHS_simdata/J50_NHHS_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished nhhs sample", sample_id, "\n")



####################################################




############################
# Fit NHHSSI model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_NHHS_SI_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW=WTW, n_all_par = n_all_par,
    J = J, M = M, O = O,
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




dir.create("/projects/soniast@colostate.edu/results_NHHSSI_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_NHHSSI_simdata/J50_NHHSSI_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished nhhssi sample", sample_id, "\n")



#############################################################



############################
# Fit NHRHS model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_NHRHS_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    n_all_par = n_all_par,
    J = J, M = M, O,
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

  res$theta_update

}, error=function(e){
  list(
    error=TRUE,
    message=e$message,
    sample_id=sample_id
  )

})




dir.create("/projects/soniast@colostate.edu/results_NHRHS_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_NHRHS_simdata/J50_NHRHS_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished nhrhs sample", sample_id, "\n")





############################
# Fit NHRHSSI model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_BVS_NHRHS_SI_MVN_cov_modify2(
    niter = 6000, burn_in = 1000, thin = 5,
    n=n, K=K, Y=Y, W=W, WTW = WTW, n_all_par=n_all_par,
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
    sample_id=sample_id
  )

})




dir.create("/projects/soniast@colostate.edu/results_NHRHSSI_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_NHRHSSI_simdata/J50_NHRHSSI_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished nhrhssi sample", sample_id, "\n")


######################################################


############################
# Fit NP model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_NP_MVN_cov_modify(
    niter = 6000, burn_in = 1000, thin = 5,
    n=n, K=K, Y=Y, W=W, WTW=WTW, n_all_par=n_all_par,
    J=J, M=M, O=O,
    theta_init = matrix(0.5, nrow = n_all_par, ncol = K),
    Sigma_init = Sigma_init,
    nu_0 = nu_0, Psi_0 = Psi_0,
    sigmasq_varphi = 10)

  res$theta_update

}, error=function(e){
  list(
    error=TRUE,
    message=e$message,
    sample_id=sample_id
  )

})




dir.create("/projects/soniast@colostate.edu/results_NP_simdata", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_NP_simdata/J50_NP_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished np sample", sample_id, "\n")


