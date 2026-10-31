
rm(list = ls())
gc()


########################
library(LaplacesDemon)
library(Matrix)


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



############################################################################
# β, γ, δ, φ ∼ N(0, 1), and α ∼ Unif (−3.5, −2.5), all predictors simulated
############################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5


# First data
set.seed(23444)
sim_data_normal01_simpred1 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)


Y <- sim_data_normal01_simpred1$Y
x <- sim_data_normal01_simpred1$x
z <- sim_data_normal01_simpred1$z
u <- sim_data_normal01_simpred1$u
d <- sim_data_normal01_simpred1$d
W <- sim_data_normal01_simpred1$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(456557)
res_normal01_simpred_hhssi1 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_simpred_hhssi1, file="practice_data/res_normal01_simpred_hhssi1.rda")
# load(file="practice_data/res_normal01_simpred_hhssi1.rda")


start_time <- Sys.time()
set.seed(456557)
res_normal01_simpred_nhhssi1 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time

# save(res_normal01_simpred_nhhssi1, file="practice_data/res_normal01_simpred_nhhssi1.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi1.rda")


# Second data
# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
set.seed(654586798)
sim_data_normal01_simpred2 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)


Y <- sim_data_normal01_simpred2$Y
x <- sim_data_normal01_simpred2$x
z <- sim_data_normal01_simpred2$z
u <- sim_data_normal01_simpred2$u
d <- sim_data_normal01_simpred2$d
W <- sim_data_normal01_simpred2$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(543565778)
res_normal01_simpred_hhssi2 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_simpred_hhssi2, file="practice_data/res_normal01_simpred_hhssi2.rda")
# load(file="practice_data/res_normal01_simpred_hhssi2.rda")


start_time <- Sys.time()
set.seed(543565778)
res_normal01_simpred_nhhssi2 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time

# save(res_normal01_simpred_nhhssi2, file="practice_data/res_normal01_simpred_nhhssi2.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi2.rda")




# Third data
# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
set.seed(23315644)
sim_data_normal01_simpred3 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)


Y <- sim_data_normal01_simpred3$Y
x <- sim_data_normal01_simpred3$x
z <- sim_data_normal01_simpred3$z
u <- sim_data_normal01_simpred3$u
d <- sim_data_normal01_simpred3$d
W <- sim_data_normal01_simpred3$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(546578557)
res_normal01_simpred_hhssi3 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_simpred_hhssi3, file="practice_data/res_normal01_simpred_hhssi3.rda")
# load(file="practice_data/res_normal01_simpred_hhssi3.rda")


start_time <- Sys.time()
set.seed(546578557)
res_normal01_simpred_nhhssi3 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time

# save(res_normal01_simpred_nhhssi3, file="practice_data/res_normal01_simpred_nhhssi3.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi3.rda")






######################################################################
# β, γ, δ, φ ∼ N(0, 1), and α ∼ Unif (−3.5, −2.5), all predictors real
########################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
# First sample
set.seed(23444)
sim_data_normal01_realpred1 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D,
                              Sigma=Sigma)



Y <- sim_data_normal01_realpred1$Y
x <- sim_data_normal01_realpred1$x
z <- sim_data_normal01_realpred1$z
u <- sim_data_normal01_realpred1$u
d <- sim_data_normal01_realpred1$d
W <- sim_data_normal01_realpred1$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(456557)
res_normal01_realpred_hhssi1 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_realpred_hhssi1, file="practice_data/res_normal01_realpred_hhssi1.rda")
# load(file="practice_data/res_normal01_realpred_hhssi1.rda")


start_time <- Sys.time()
set.seed(456557)
res_normal01_realpred_nhhssi1 <- fit_NHHSSI_MVN_cov_modify3(
    niter = 2000, burn_in = 500, thin = 5,
    n=n, K=K, Y=Y, W=W,  WTW = WTW,
    n_all_par = n_all_par, J = J, M = M, O = O,
    sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_realpred_nhhssi1, file="practice_data/res_normal01_realpred_nhhssi1.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi1.rda")







# Second sample
set.seed(654586798)
sim_data_normal01_realpred2 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                O=O, x=X, z=Z, d=D,
                                                Sigma=Sigma)



Y <- sim_data_normal01_realpred2$Y
x <- sim_data_normal01_realpred2$x
z <- sim_data_normal01_realpred2$z
u <- sim_data_normal01_realpred2$u
d <- sim_data_normal01_realpred2$d
W <- sim_data_normal01_realpred2$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(543565778)
res_normal01_realpred_hhssi2 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_realpred_hhssi2, file="practice_data/res_normal01_realpred_hhssi2.rda")
# load(file="practice_data/res_normal01_realpred_hhssi2.rda")


start_time <- Sys.time()
set.seed(543565778)
res_normal01_realpred_nhhssi2 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_realpred_nhhssi2, file="practice_data/res_normal01_realpred_nhhssi2.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi2.rda")




# Third sample
set.seed(23315644)
sim_data_normal01_realpred3 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                O=O, x=X, z=Z, d=D,
                                                Sigma=Sigma)



Y <- sim_data_normal01_realpred3$Y
x <- sim_data_normal01_realpred3$x
z <- sim_data_normal01_realpred3$z
u <- sim_data_normal01_realpred3$u
d <- sim_data_normal01_realpred3$d
W <- sim_data_normal01_realpred3$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(546578557)
res_normal01_realpred_hhssi3 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_realpred_hhssi3, file="practice_data/res_normal01_realpred_hhssi3.rda")
# load(file="practice_data/res_normal01_realpred_hhssi3.rda")


start_time <- Sys.time()
set.seed(546578557)
res_normal01_realpred_nhhssi3 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_normal01_realpred_nhhssi3, file="practice_data/res_normal01_realpred_nhhssi3.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi3.rda")






########################################################
# α, β, γ, δ, φ ∼ unif(2,3), all predictors simulated
#######################################################


# First sample
set.seed(23444)
sim_data_unif23_simpred1 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)



Y <- sim_data_unif23_simpred1$Y
x <- sim_data_unif23_simpred1$x
z <- sim_data_unif23_simpred1$z
u <- sim_data_unif23_simpred1$u
d <- sim_data_unif23_simpred1$d
W <- sim_data_unif23_simpred1$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(456557)
res_unif23_simpred_hhssi1 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_simpred_hhssi1, file="practice_data/res_unif23_simpred_hhssi1.rda")
# load(file="practice_data/res_unif23_simpred_hhssi1.rda")


start_time <- Sys.time()
set.seed(456557)
res_unif23_simpred_nhhssi1 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_simpred_nhhssi1, file="practice_data/res_unif23_simpred_nhhssi1.rda")
# load(file="practice_data/res_unif23_simpred_nhhssi1.rda")



# Second sample
set.seed(654586798)
sim_data_unif23_simpred2 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)



Y <- sim_data_unif23_simpred2$Y
x <- sim_data_unif23_simpred2$x
z <- sim_data_unif23_simpred2$z
u <- sim_data_unif23_simpred2$u
d <- sim_data_unif23_simpred2$d
W <- sim_data_unif23_simpred2$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(543565778)
res_unif23_simpred_hhssi2 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_simpred_hhssi2, file="practice_data/res_unif23_simpred_hhssi2.rda")
# load(file="practice_data/res_unif23_simpred_hhssi2.rda")


start_time <- Sys.time()
set.seed(543565778)
res_unif23_simpred_nhhssi2 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_simpred_nhhssi2, file="practice_data/res_unif23_simpred_nhhssi2.rda")
# load(file="practice_data/res_unif23_simpred_nhhssi2.rda")




# Third sample
set.seed(23315644)
sim_data_unif23_simpred3 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)



Y <- sim_data_unif23_simpred3$Y
x <- sim_data_unif23_simpred3$x
z <- sim_data_unif23_simpred3$z
u <- sim_data_unif23_simpred3$u
d <- sim_data_unif23_simpred3$d
W <- sim_data_unif23_simpred3$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(546578557)
res_unif23_simpred_hhssi3 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_simpred_hhssi3, file="practice_data/res_unif23_simpred_hhssi3.rda")
# load(file="practice_data/res_unif23_simpred_hhssi3.rda")


start_time <- Sys.time()
set.seed(546578557)
res_unif23_simpred_nhhssi3 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_simpred_nhhssi3, file="practice_data/res_unif23_simpred_nhhssi3.rda")
# load(file="practice_data/res_unif23_simpred_nhhssi3.rda")







#################################################
# α, β, γ, δ, φ ∼ unif(2,3), all predictors real
##################################################

# First sample
set.seed(23444)
sim_data_unif23_realpred1 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                O=O, x=X, z=Z, d=D,
                                                Sigma=Sigma)



Y <- sim_data_unif23_realpred1$Y
x <- sim_data_unif23_realpred1$x
z <- sim_data_unif23_realpred1$z
u <- sim_data_unif23_realpred1$u
d <- sim_data_unif23_realpred1$d
W <- sim_data_unif23_realpred1$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)



start_time <- Sys.time()
set.seed(456557)
res_unif23_realpred_hhssi1 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_realpred_hhssi1, file="practice_data/res_unif23_realpred_hhssi1.rda")
# load(file="practice_data/res_unif23_realpred_hhssi1.rda")



start_time <- Sys.time()
set.seed(456557)
res_unif23_realpred_nhhssi1 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_realpred_nhhssi1, file="practice_data/res_unif23_realpred_nhhssi1.rda")
# load(file="practice_data/res_unif23_realpred_nhhssi1.rda")





# Second sample
set.seed(654586798)
sim_data_unif23_realpred2 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                              O=O, x=X, z=Z, d=D,
                                              Sigma=Sigma)



Y <- sim_data_unif23_realpred2$Y
x <- sim_data_unif23_realpred2$x
z <- sim_data_unif23_realpred2$z
u <- sim_data_unif23_realpred2$u
d <- sim_data_unif23_realpred2$d
W <- sim_data_unif23_realpred2$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)



start_time <- Sys.time()
set.seed(543565778)
res_unif23_realpred_hhssi2 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_realpred_hhssi2, file="practice_data/res_unif23_realpred_hhssi2.rda")
# load(file="practice_data/res_unif23_realpred_hhssi2.rda")



start_time <- Sys.time()
set.seed(543565778)
res_unif23_realpred_nhhssi2 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_realpred_nhhssi2, file="practice_data/res_unif23_realpred_nhhssi2.rda")
# load(file="practice_data/res_unif23_realpred_nhhssi2.rda")





# Third sample
set.seed(23315644)
sim_data_unif23_realpred3 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                              O=O, x=X, z=Z, d=D,
                                              Sigma=Sigma)



Y <- sim_data_unif23_realpred3$Y
x <- sim_data_unif23_realpred3$x
z <- sim_data_unif23_realpred3$z
u <- sim_data_unif23_realpred3$u
d <- sim_data_unif23_realpred3$d
W <- sim_data_unif23_realpred3$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)



start_time <- Sys.time()
set.seed(546578557)
res_unif23_realpred_hhssi3 <- fit_HHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n = n, K = K, Y = Y, W = W, WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_realpred_hhssi3, file="practice_data/res_unif23_realpred_hhssi3.rda")
# load(file="practice_data/res_unif23_realpred_hhssi3.rda")



start_time <- Sys.time()
set.seed(546578557)
res_unif23_realpred_nhhssi3 <- fit_NHHSSI_MVN_cov_modify3(
  niter = 2000, burn_in = 500, thin = 5,
  n=n, K=K, Y=Y, W=W,  WTW = WTW,
  n_all_par = n_all_par, J = J, M = M, O = O,
  sigmasq_alpha = 100,
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
end_time <- Sys.time()
end_time - start_time


# save(res_unif23_realpred_nhhssi3, file="practice_data/res_unif23_realpred_nhhssi3.rda")
# load(file="practice_data/res_unif23_realpred_nhhssi3.rda")







###############################################################################
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors simulated (First sample)
###############################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
# First sample
set.seed(23444)
sim_data_normal01_simpred1 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)

Y <- sim_data_normal01_simpred1$Y
x <- sim_data_normal01_simpred1$x
z <- sim_data_normal01_simpred1$z
u <- sim_data_normal01_simpred1$u
d <- sim_data_normal01_simpred1$d
W <- sim_data_normal01_simpred1$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_simpred_hhssi1.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi1.rda")

sim_data_normal01_simpred1$delta_true



dim(res_normal01_simpred_hhssi1$theta_update)

# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
   matplot(res_normal01_simpred_hhssi1$theta_update[, 1, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = expression(alpha),
    main = paste("Outcome", k)
  )

  matlines(res_normal01_simpred_nhhssi1$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_simpred_hhssi1$theta_update[, 1, ])
range(res_normal01_simpred_nhhssi1$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
hhssi_alpha <- res_normal01_simpred_hhssi1$theta_update[, 1, k]
nhhssi_alpha <- res_normal01_simpred_nhhssi1$theta_update[, 1, k]

plot(hhssi_alpha, nhhssi_alpha,
    pch = 19,
    col = "blue",
    xlab = "HHS-SI",
    ylab = "NHHS-SI",
    main = paste("alpha: Outcome", k),
    xlim = range(c(hhssi_alpha, nhhssi_alpha)),
    ylim = range(c(hhssi_alpha, nhhssi_alpha)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
}
par(mfrow = c(1, 1))



# bias alpha
sim_data_normal01_simpred1$alpha_true
res_normal01_simpred_hhssi1$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_simpred_hhssi1$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_simpred_nhhssi1$theta_update[, 1, ], 2, mean)
# Bias for each outcome
alpha_bias_hhssi <- (alpha_est_hhssi - sim_data_normal01_simpred1$alpha_true)
alpha_bias_nhhssi <- (alpha_est_nhhssi - sim_data_normal01_simpred1$alpha_true)


plot(alpha_bias_hhssi,
     type = "b",
     xlab = "Outcome",
     ylab = "bias of alpha",
     lwd = 1,
     col = "red",
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
lines(alpha_bias_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(alpha_bias_hhssi, alpha_bias_nhhssi,
     lwd = 1,
     xlim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)),
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
abline(a=0, b=1, col = "red")



# active beta

plot(res_normal01_simpred_hhssi1$theta_update[, 2:3, 1],
     res_normal01_simpred_nhhssi1$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi1$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi1$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, 2:3, 2],
       res_normal01_simpred_nhhssi1$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, 2:3, 3],
       res_normal01_simpred_nhhssi1$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, 2:3, 4],
       res_normal01_simpred_nhhssi1$theta_update[, 2:3, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))




range(res_normal01_simpred_hhssi1$theta_update[, 2:3, 2])
range(res_normal01_simpred_nhhssi1$theta_update[, 2:3, 2])

range(res_normal01_simpred_hhssi1$theta_update[, 2:3, 3])
range(res_normal01_simpred_nhhssi1$theta_update[, 2:3, 3])

range(res_normal01_simpred_hhssi1$theta_update[, 2:3, 4])
range(res_normal01_simpred_nhhssi1$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, 2:3, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))




# inactive beta
range(res_normal01_simpred_hhssi1$theta_update[, 4:45, 1])
range(res_normal01_simpred_nhhssi1$theta_update[, 4:45, 1])

plot(res_normal01_simpred_hhssi1$theta_update[, 4:45, 1],
     res_normal01_simpred_nhhssi1$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi1$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi1$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, 4:45, 2],
       res_normal01_simpred_nhhssi1$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, 4:45, 3],
       res_normal01_simpred_nhhssi1$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, 4:45, 4],
       res_normal01_simpred_nhhssi1$theta_update[, 4:45, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



range(res_normal01_simpred_hhssi1$theta_update[, 4:45, 2])
range(res_normal01_simpred_nhhssi1$theta_update[, 4:45, 2])

range(res_normal01_simpred_hhssi1$theta_update[, 4:45, 3])
range(res_normal01_simpred_nhhssi1$theta_update[, 4:45, 3])

range(res_normal01_simpred_hhssi1$theta_update[, 4:45, 4])
range(res_normal01_simpred_nhhssi1$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, 4:45, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))






# bias beta
dim(sim_data_normal01_simpred1$beta_true)
dim(res_normal01_simpred_hhssi1$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_simpred_hhssi1$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_simpred_nhhssi1$theta_update[, 2:45, ], c(2,3), mean)
# Bias for each outcome
beta_bias_hhssi <- (beta_est_hhssi - sim_data_normal01_simpred1$beta_true)
beta_bias_nhhssi <- (beta_est_nhhssi - sim_data_normal01_simpred1$beta_true)

dim(beta_bias_hhssi)
dim(beta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:44,
    beta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "modifiers",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(beta_bias_hhssi,
                   beta_bias_nhhssi)))
  lines(1:44,
        beta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  beta_bias_hhssi[, k],
  beta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(beta_bias_hhssi, beta_bias_nhhssi)),
  ylim = range(c(beta_bias_hhssi, beta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))





mean(beta_bias_hhssi>beta_bias_nhhssi)



beta_var_hhssi <- (res_normal01_simpred_hhssi1$lambdasq_beta_update *
                     res_normal01_simpred_hhssi1$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_simpred_nhhssi1$lambdasq_beta_update *
                      res_normal01_simpred_nhhssi1$tausq_beta_update)

mean(beta_var_hhssi>beta_var_nhhssi)
mean(beta_var_hhssi)
mean(beta_var_nhhssi)
mean(beta_bias_hhssi)
mean(beta_bias_nhhssi)

mean(beta_var_hhssi[1:2] > beta_var_nhhssi[1:2])
mean(beta_var_hhssi[1:2])
mean(beta_var_nhhssi[1:2])
mean(beta_bias_hhssi[1:2])
mean(beta_bias_nhhssi[1:2])

mean(beta_var_hhssi[-(1:2)] > beta_var_nhhssi[-(1:2)])
mean(beta_var_hhssi[-(1:2)])
mean(beta_var_nhhssi[-(1:2)])
mean(beta_bias_hhssi[-(1:2)])
mean(beta_bias_nhhssi[-(1:2)])



plot(beta_var_hhssi, beta_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[1:2], beta_var_nhhssi[1:2],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[-(1:2)], beta_var_nhhssi[-(1:2)],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")












# active gamma
range(res_normal01_simpred_hhssi1$theta_update[, 46:48, 1])
range(res_normal01_simpred_nhhssi1$theta_update[, 46:48, 1])


plot(res_normal01_simpred_hhssi1$theta_update[, 46:48, 1],
     res_normal01_simpred_nhhssi1$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi1$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi1$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, 46:48, 2],
       res_normal01_simpred_nhhssi1$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, 46:48, 3],
       res_normal01_simpred_nhhssi1$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, 46:48, 4],
       res_normal01_simpred_nhhssi1$theta_update[, 46:48, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_simpred_hhssi1$theta_update[, 46:48, 2])
range(res_normal01_simpred_nhhssi1$theta_update[, 46:48, 2])

range(res_normal01_simpred_hhssi1$theta_update[, 46:48, 3])
range(res_normal01_simpred_nhhssi1$theta_update[, 46:48, 3])

range(res_normal01_simpred_hhssi1$theta_update[, 46:48, 4])
range(res_normal01_simpred_nhhssi1$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, 46:48, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive gamma
range(res_normal01_simpred_hhssi1$theta_update[, 49:62, 1])
range(res_normal01_simpred_nhhssi1$theta_update[, 49:62, 1])


plot(res_normal01_simpred_hhssi1$theta_update[, 49:62, 1],
     res_normal01_simpred_nhhssi1$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi1$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi1$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, 49:62, 2],
       res_normal01_simpred_nhhssi1$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, 49:62, 3],
       res_normal01_simpred_nhhssi1$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, 49:62, 4],
       res_normal01_simpred_nhhssi1$theta_update[, 49:62, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi1$theta_update[, 49:62, 2])
range(res_normal01_simpred_nhhssi1$theta_update[, 49:62, 2])

range(res_normal01_simpred_hhssi1$theta_update[, 49:62, 3])
range(res_normal01_simpred_nhhssi1$theta_update[, 49:62, 3])

range(res_normal01_simpred_hhssi1$theta_update[, 49:62, 4])
range(res_normal01_simpred_nhhssi1$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, 49:62, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))






# bias gamma
dim(sim_data_normal01_simpred1$gamma_true)
dim(res_normal01_simpred_hhssi1$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_simpred_hhssi1$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_simpred_nhhssi1$theta_update[, 46:62, ], c(2,3), mean)
# Bias for each outcome
gamma_bias_hhssi <- (gamma_est_hhssi - sim_data_normal01_simpred1$gamma_true)
gamma_bias_nhhssi <- (gamma_est_nhhssi - sim_data_normal01_simpred1$gamma_true)

dim(gamma_bias_hhssi)
dim(gamma_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:17,
    gamma_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "exposures",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(gamma_bias_hhssi,
                   gamma_bias_nhhssi)))
  lines(1:17,
        gamma_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  gamma_bias_hhssi[, k],
  gamma_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)),
  ylim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))



gamma_var_hhssi <- (res_normal01_simpred_hhssi1$lambdasq_gamma_update *
                      res_normal01_simpred_hhssi1$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_simpred_nhhssi1$lambdasq_gamma_update *
                       res_normal01_simpred_nhhssi1$tausq_gamma_update)

mean(gamma_var_hhssi>gamma_var_nhhssi)
mean(gamma_var_hhssi)
mean(gamma_var_nhhssi)
mean(gamma_bias_hhssi)
mean(gamma_bias_nhhssi)

mean(gamma_var_hhssi[1:2] > gamma_var_nhhssi[1:2])
mean(gamma_var_hhssi[1:2])
mean(gamma_var_nhhssi[1:2])
mean(gamma_bias_hhssi[1:2])
mean(gamma_bias_nhhssi[1:2])

mean(gamma_var_hhssi[-(1:2)] > gamma_var_nhhssi[-(1:2)])
mean(gamma_var_hhssi[-(1:2)])
mean(gamma_var_nhhssi[-(1:2)])
mean(gamma_bias_hhssi[-(1:2)])
mean(gamma_bias_nhhssi[-(1:2)])



plot(gamma_var_hhssi, gamma_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[1:3], gamma_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[-(1:3)], gamma_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")










# active delta
range(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 1],
     res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 2],
       res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 3],
       res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 4],
       res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_simpred_hhssi1$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_simpred_nhhssi1$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, c(63:65,80:82), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive delta

range(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 1],
     res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 2],
       res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 3],
       res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 4],
       res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_simpred_hhssi1$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_simpred_nhhssi1$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, c(66:79,83:810), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# bias delta
dim(sim_data_normal01_simpred1$delta_true)
dim(res_normal01_simpred_hhssi1$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_simpred_hhssi1$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_simpred_nhhssi1$theta_update[, 63:810, ], c(2,3), mean)
# Bias for each outcome
delta_bias_hhssi <- (delta_est_hhssi - sim_data_normal01_simpred1$delta_true)
delta_bias_nhhssi <- (delta_est_nhhssi - sim_data_normal01_simpred1$delta_true)

dim(delta_bias_hhssi)
dim(delta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:748,
    delta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "interactions",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(delta_bias_hhssi,
                   delta_bias_nhhssi)))
  lines(1:748,
        delta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  delta_bias_hhssi[, k],
  delta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(delta_bias_hhssi, delta_bias_nhhssi)),
  ylim = range(c(delta_bias_hhssi, delta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






# active varphi
range(res_normal01_simpred_hhssi1$theta_update[, 811:812, 1])
range(res_normal01_simpred_nhhssi1$theta_update[, 811:812, 1])

plot(res_normal01_simpred_hhssi1$theta_update[, 811:812, 1],
     res_normal01_simpred_nhhssi1$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi1$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi1$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, 811:812, 2],
       res_normal01_simpred_nhhssi1$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, 811:812, 3],
       res_normal01_simpred_nhhssi1$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, 811:812, 4],
       res_normal01_simpred_nhhssi1$theta_update[, 811:812, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi1$theta_update[, 811:812, 2])
range(res_normal01_simpred_nhhssi1$theta_update[, 811:812, 2])

range(res_normal01_simpred_hhssi1$theta_update[, 811:812, 3])
range(res_normal01_simpred_nhhssi1$theta_update[, 811:812, 3])

range(res_normal01_simpred_hhssi1$theta_update[, 811:812, 4])
range(res_normal01_simpred_nhhssi1$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, 811:812, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive varphi
range(res_normal01_simpred_hhssi1$theta_update[, 813:816, 1])
range(res_normal01_simpred_nhhssi1$theta_update[, 813:816, 1])

plot(res_normal01_simpred_hhssi1$theta_update[, 813:816, 1],
     res_normal01_simpred_nhhssi1$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi1$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_simpred_hhssi1$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi1$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_simpred_hhssi1$theta_update[, 813:816, 2],
       res_normal01_simpred_nhhssi1$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$theta_update[, 813:816, 3],
       res_normal01_simpred_nhhssi1$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$theta_update[, 813:816, 4],
       res_normal01_simpred_nhhssi1$theta_update[, 813:816, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_simpred_hhssi1$theta_update[, 813:816, 2])
range(res_normal01_simpred_nhhssi1$theta_update[, 813:816, 2])

range(res_normal01_simpred_hhssi1$theta_update[, 813:816, 3])
range(res_normal01_simpred_nhhssi1$theta_update[, 813:816, 3])

range(res_normal01_simpred_hhssi1$theta_update[, 813:816, 4])
range(res_normal01_simpred_nhhssi1$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi1$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi1$theta_update[, 813:816, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))





# bias varphi
dim(sim_data_normal01_simpred1$varphi_true)
dim(res_normal01_simpred_hhssi1$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_simpred_hhssi1$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_simpred_nhhssi1$theta_update[, 811:816, ], c(2,3), mean)
# Bias for each outcome
varphi_bias_hhssi <- (varphi_est_hhssi - sim_data_normal01_simpred1$varphi_true)
varphi_bias_nhhssi <- (varphi_est_nhhssi - sim_data_normal01_simpred1$varphi_true)

dim(varphi_bias_hhssi)
dim(varphi_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:6,
    varphi_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "covariates",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(varphi_bias_hhssi,
                   varphi_bias_nhhssi)))
  lines(1:6,
        varphi_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  varphi_bias_hhssi[, k],
  varphi_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)),
  ylim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






varphi_var_hhssi <- apply(res_normal01_simpred_hhssi1$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_simpred_nhhssi1$theta_update[, 811:816, ], c(2,3), var)

mean(varphi_var_hhssi>varphi_var_nhhssi)
mean(varphi_var_hhssi)
mean(varphi_var_nhhssi)
mean(varphi_bias_hhssi)
mean(varphi_bias_nhhssi)

mean(varphi_var_hhssi[1:2] > varphi_var_nhhssi[1:2])
mean(varphi_var_hhssi[1:2])
mean(varphi_var_nhhssi[1:2])
mean(varphi_bias_hhssi[1:2])
mean(varphi_bias_nhhssi[1:2])

mean(varphi_var_hhssi[-(1:2)] > varphi_var_nhhssi[-(1:2)])
mean(varphi_var_hhssi[-(1:2)])
mean(varphi_var_nhhssi[-(1:2)])
mean(varphi_bias_hhssi[-(1:2)])
mean(varphi_bias_nhhssi[-(1:2)])



plot(varphi_var_hhssi, varphi_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[1:3], varphi_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[-(1:3)], varphi_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")











dim(sim_data_normal01_simpred1$gamma_true)
dim(res_normal01_simpred_hhssi1$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_simpred_hhssi1$lambdasq_delta_update)
dim(res_normal01_simpred_nhhssi1$lambdasq_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi1$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi1$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi1$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi1$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_simpred_hhssi1$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_simpred_nhhssi1$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_simpred_hhssi1$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_simpred_nhhssi1$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi1$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_simpred_nhhssi1$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_simpred_hhssi1$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_simpred_nhhssi1$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_simpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_simpred_nhhssi1$lambdasq_delta_update, 2, mean)
plot(lambdasq_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


# active deltas
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[c(1:3, 18:20)])
hist(lambdasq_delta_nhhssi[c(1:3, 18:20)])
range(lambdasq_delta_hhssi[c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[c(1:3, 18:20)])) # active deltas
hist(log(lambdasq_delta_nhhssi[c(1:3, 18:20)])) # active deltas

# inactive deltas
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[-c(1:3, 18:20)]) # inactive deltas
hist(lambdasq_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(lambdasq_delta_nhhssi[-c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)])) # inactive deltas
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)])) # inactive deltas


# active delta
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive delta
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")





# Note: the active lambdasq_delta values are lower for hhssi model than nhhssi, and range of values is (0.1485711, 2.5748392)
# the active lambdasq_delta values are higher for nhhssi model, and range of values is (58819.5, 403810.4)

# Note: whereas, the inactive lambdasq_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (1.246500e-01, 1.844546e+10)
# the inactive lambdasq_delta values are lower for nhhssi model, and range of values is (2.776089 234.788481)

# Overall, we can say that lambdasq_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(lambdasq_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(lambdasq_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
lambdasq_delta_hhssi_new <- lambdasq_delta_hhssi[-remove_hhssi]
lambdasq_delta_nhhssi_new <- lambdasq_delta_nhhssi[-remove_nhhssi]
# Plot
plot(lambdasq_delta_hhssi_new, lambdasq_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
lambdasq_delta_hhssi[remove_hhssi]
lambdasq_delta_nhhssi[remove_nhhssi]




hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", lambda[delta]^2)),
     xlab = expression(lambda[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]), log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]))))
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

mean(log(lambdasq_delta_hhssi))
mean(log(lambdasq_delta_nhhssi))

mean(lambdasq_delta_hhssi)
mean(lambdasq_delta_nhhssi)

# 1:2 active betas
plot(apply(res_normal01_simpred_hhssi1$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update, 2, mean)
plot(lambdasq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[beta]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(lambdasq_beta_hhssi, lambdasq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)),
     ylim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_beta looks similar for hhssi and nhhssi models





# 1:3 active gammas
plot(apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update, 2, mean)
plot(lambdasq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[gamma]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)),
     ylim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_gamma looks similar for hhssi and nhhssi models




names(res_normal01_simpred_hhssi1)
names(res_normal01_simpred_nhhssi1)





# tausq beta
mean(res_normal01_simpred_hhssi1$tausq_beta_update)
mean(res_normal01_simpred_nhhssi1$tausq_beta_update)
abs(mean(res_normal01_simpred_hhssi1$tausq_beta_update) -
      mean(res_normal01_simpred_nhhssi1$tausq_beta_update))

range(res_normal01_simpred_hhssi1$tausq_beta_update)
range(res_normal01_simpred_nhhssi1$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_simpred_hhssi1$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_simpred_nhhssi1$tausq_beta_update
plot(tausq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[beta]^2),
     lwd = 1,
     col = "red")
lines(tausq_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_beta_hhssi, tausq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)),
     ylim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
abline(a=0, b=1, col = "red")




hist(tausq_beta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
hist(tausq_beta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))


d1 <- density(tausq_beta_hhssi)
d2 <- density(tausq_beta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)




# tausq_gamma
mean(res_normal01_simpred_hhssi1$tausq_gamma_update)
mean(res_normal01_simpred_nhhssi1$tausq_gamma_update)
abs(mean(res_normal01_simpred_hhssi1$tausq_gamma_update) -
      mean(res_normal01_simpred_nhhssi1$tausq_gamma_update))

range(res_normal01_simpred_hhssi1$tausq_gamma_update)
range(res_normal01_simpred_nhhssi1$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_simpred_hhssi1$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_simpred_nhhssi1$tausq_gamma_update
plot(tausq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[gamma]^2),
     lwd = 1,
     col = "red")
lines(tausq_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_gamma_hhssi, tausq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)),
     ylim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")






# tausq_delta
mean(res_normal01_simpred_hhssi1$tausq_delta_update)
mean(res_normal01_simpred_nhhssi1$tausq_delta_update)
abs(mean(res_normal01_simpred_hhssi1$tausq_delta_update) -
      mean(res_normal01_simpred_nhhssi1$tausq_delta_update))

range(res_normal01_simpred_hhssi1$tausq_delta_update)
range(res_normal01_simpred_nhhssi1$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_simpred_hhssi1$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi1$tausq_delta_update
plot(tausq_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[delta]^2),
     lwd = 1,
     col = "red",
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
lines(tausq_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

hist(tausq_delta_hhssi, breaks = 50)
hist(tausq_delta_nhhssi, breaks = 50)

plot(tausq_delta_hhssi, tausq_delta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
abline(a=0, b=1, col = "red")



hist(tausq_delta_hhssi)
hist(tausq_delta_nhhssi)

hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

# log x axis
hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     log = "x")
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))



d1 <- density(tausq_delta_hhssi)
d2 <- density(tausq_delta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)


# Note: tausq_beta and tausq_gamma are same for the hhssi and nhhssi models.
# However, tausq_delta from hhssi and nhhssi models are far different.
# tausq_delta from hhssi are higher than that of nhhssi model.






# xi beta
mean(res_normal01_simpred_hhssi1$xi_beta_update)
mean(res_normal01_simpred_nhhssi1$xi_beta_update)
abs(mean(res_normal01_simpred_hhssi1$xi_beta_update) -
      mean(res_normal01_simpred_nhhssi1$xi_beta_update))

range(res_normal01_simpred_hhssi1$xi_beta_update)
range(res_normal01_simpred_nhhssi1$xi_beta_update)

xi_beta_hhssi <- res_normal01_simpred_hhssi1$xi_beta_update
xi_beta_nhhssi <- res_normal01_simpred_nhhssi1$xi_beta_update
plot(xi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[beta]),
     lwd = 1,
     col = "red")
lines(xi_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_beta_hhssi, xi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_beta_hhssi, xi_beta_nhhssi)),
     ylim = range(c(xi_beta_hhssi, xi_beta_nhhssi)))
abline(a=0, b=1, col = "red")





# xi_gamma
mean(res_normal01_simpred_hhssi1$xi_gamma_update)
mean(res_normal01_simpred_nhhssi1$xi_gamma_update)
abs(mean(res_normal01_simpred_hhssi1$xi_gamma_update) -
      mean(res_normal01_simpred_nhhssi1$xi_gamma_update))

range(res_normal01_simpred_hhssi1$xi_gamma_update)
range(res_normal01_simpred_nhhssi1$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_simpred_hhssi1$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_simpred_nhhssi1$xi_gamma_update
plot(xi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[gamma]),
     lwd = 1,
     col = "red")
lines(xi_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_gamma_hhssi, xi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)),
     ylim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")




# xi_delta
mean(res_normal01_simpred_hhssi1$xi_delta_update)
mean(res_normal01_simpred_nhhssi1$xi_delta_update)
abs(mean(res_normal01_simpred_hhssi1$xi_delta_update) -
      mean(res_normal01_simpred_nhhssi1$xi_delta_update))

range(res_normal01_simpred_hhssi1$xi_delta_update)
range(res_normal01_simpred_nhhssi1$xi_delta_update)

xi_delta_hhssi <- res_normal01_simpred_hhssi1$xi_delta_update
xi_delta_nhhssi <- res_normal01_simpred_nhhssi1$xi_delta_update
plot(xi_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[delta]),
     lwd = 1,
     col = "red",
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
lines(xi_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_delta_hhssi, xi_delta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_delta_hhssi, xi_delta_nhhssi)),
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
abline(a=0, b=1, col = "red")

hist(xi_delta_hhssi)
hist(xi_delta_nhhssi)



# Note: xi_beta and xi_gamma are same for both the hhssi and nhhssi models.
# However, xi_delta from hhssi and nhhssi models are far different.
# xi_delta from nhhssi are much higher than that of hhssi model.

##########################################################

# psi_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_simpred_nhhssi1$psi_delta_update)
dim(res_normal01_simpred_hhssi1$psi_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi1$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi1$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi1$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi1$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_simpred_hhssi1$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_simpred_nhhssi1$psi_delta_update, 2, mean)

# active
plot(psi_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive
plot(psi_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")



# active deltas
plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])
hist(psi_delta_hhssi[c(1:3, 18:20)])

range(psi_delta_hhssi[c(1:3, 18:20)])

# inactive deltas
plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])
hist(psi_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(psi_delta_nhhssi[-c(1:3, 18:20)])


# Note: the active psi_delta values are lower for hhssi model than nhhssi, and range of values is (0.1485711, 2.5748392)
# the active psi_delta values are higher for nhhssi model, and range of values is (58819.5, 403810.4)

# Note: whereas, the inactive psi_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (1.246500e-01, 1.844546e+10)
# the inactive psi_delta values are lower for nhhssi model, and range of values is (2.776089 234.788481)

# Overall, we can say that psi_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(psi_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(psi_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
psi_delta_hhssi_new <- psi_delta_hhssi[-remove_hhssi]
psi_delta_nhhssi_new <- psi_delta_nhhssi[-remove_nhhssi]
# Plot
plot(psi_delta_hhssi_new, psi_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
psi_delta_hhssi[remove_hhssi]
psi_delta_nhhssi[remove_nhhssi]




hist(log(psi_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", psi[delta]^2)),
     xlab = expression(psi[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(psi_delta_hhssi), log(psi_delta_nhhssi))))
hist(log(psi_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))




# 1:2 active betas
plot(apply(res_normal01_simpred_hhssi1$psi_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi1$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_simpred_hhssi1$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_simpred_nhhssi1$psi_beta_update, 2, mean)
plot(psi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[beta])),
     lwd = 1,
     col = "red")
lines(psi_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_beta_hhssi, psi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(psi_beta_hhssi, psi_beta_nhhssi)),
     ylim = range(c(psi_beta_hhssi, psi_beta_nhhssi)))
abline(a=0, b=1, col = "red")


# Note: psi_beta looks similar for hhssi and nhhssi models


# 1:3 active gammas
plot(apply(res_normal01_simpred_hhssi1$psi_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi1$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_simpred_hhssi1$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi1$psi_gamma_update, 2, mean)
plot(psi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[gamma])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
lines(psi_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_gamma_hhssi, psi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)),
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: psi_gamma looks similar for hhssi and nhhssi models


names(res_normal01_simpred_hhssi1)
names(res_normal01_simpred_nhhssi1)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_simpred_hhssi1$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_simpred_nhhssi1$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_simpred_hhssi1$Sigma_update[, , 1], type = "l")
matplot(res_normal01_simpred_nhhssi1$Sigma_update[, , 1], type = "l")

matplot(res_normal01_simpred_hhssi1$Sigma_update[, , 2], type = "l")
matplot(res_normal01_simpred_nhhssi1$Sigma_update[, , 2], type = "l")

matplot(res_normal01_simpred_hhssi1$Sigma_update[, , 3], type = "l")
matplot(res_normal01_simpred_nhhssi1$Sigma_update[, , 3], type = "l")

matplot(res_normal01_simpred_hhssi1$Sigma_update[, , 4], type = "l")
matplot(res_normal01_simpred_nhhssi1$Sigma_update[, , 4], type = "l")



plot(res_normal01_simpred_hhssi1$Sigma_update[, , 1],
     res_normal01_simpred_nhhssi1$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi1$Sigma_update,
                    res_normal01_simpred_hhssi1$Sigma_update)),
     ylim = range(c(res_normal01_simpred_hhssi1$Sigma_update,
                    res_normal01_simpred_hhssi1$Sigma_update)))
# k = 2
points(res_normal01_simpred_hhssi1$Sigma_update[, , 2],
       res_normal01_simpred_nhhssi1$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi1$Sigma_update[, , 3],
       res_normal01_simpred_nhhssi1$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi1$Sigma_update[, , 4],
       res_normal01_simpred_nhhssi1$Sigma_update[, , 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



# Note: Sigmasq are same for both hhssi and nhhssi models





###########

lambdasq_beta_hhssi
lambdasq_gamma_hhssi
lambdasq_delta_hhssi
tausq_delta_hhssi


dim(res_normal01_simpred_hhssi1$lambdasq_beta_update)
dim(res_normal01_simpred_hhssi1$lambdasq_gamma_update)
dim(res_normal01_simpred_hhssi1$lambdasq_delta_update)
length(res_normal01_simpred_hhssi1$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_simpred_hhssi1$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_simpred_hhssi1$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_simpred_hhssi1$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_simpred_hhssi1$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_simpred_nhhssi1$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi1$tausq_delta_update


n_iter <- dim(beta_lam_hhssi)[1]

# Construct beta x gamma component
beta_gamma_lam_hhssi <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_hhssi[t, col_idx] <- beta_lam_hhssi[t, j] * gamma_lam_hhssi[t, m]
      col_idx <- col_idx + 1
    }
  }
}
# prior variance
delta_var_hhssi <- delta_lam_hhssi * beta_gamma_lam_hhssi * tausq_delta_hhssi
# prior variance
delta_var_nhhssi <- delta_lam_nhhssi * tausq_delta_nhhssi



beta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi1$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_simpred_hhssi1$tausq_delta_update,
     type = "l",
     xlab = "MCMC iteration",
     ylab = expression(tau[delta]^2),
     main = "HHS-SI: global delta shrinkage"
)




hist(log10(beta_gamma_lam_hhssi),
     breaks = 50,
     main = expression("HHS-SI: " * log[10](lambda[beta]^2 * lambda[gamma]^2)),
     xlab = expression(log[10](lambda[beta]^2 * lambda[gamma]^2)))
abline(v = 0, lty = 2, lwd = 2)





delta_var_mean_hhssi <- apply(delta_var_hhssi, 2, mean)
beta_gamma_lam_mean_hhssi <- apply(beta_gamma_lam_hhssi, 2, mean)
# plot(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi,
#   log = "xy",
#   pch = 19,
#   xlab = "Delta prior var of hhssi from hhssi",
#   ylab = "Delta prior var of nhhssi from hhssi",
#   xlim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)),
#   ylim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)))
# abline(0, 1, lty = 2, lwd = 2, col = "red")


delta_var_mean_nhhssi <- apply(delta_var_nhhssi, 2, mean)
plot(delta_var_mean_hhssi, delta_var_mean_nhhssi,
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)),
     ylim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

which(delta_var_mean_hhssi[-c(1:3, 18:20)]>delta_var_mean_nhhssi[-c(1:3, 18:20)])

plot(delta_var_mean_nhhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")






delta_var_mean_hhssi[c(1:3, 18:20)] > delta_var_mean_nhhssi[c(1:3, 18:20)]
delta_var_mean_hhssi[-c(1:3, 18:20)] > delta_var_mean_nhhssi[-c(1:3, 18:20)]


ratio_hhssi_nhhssi <- delta_var_hhssi / delta_var_nhhssi
ratio_mean <- apply(ratio_hhssi_nhhssi, 2, mean)
hist(log10(ratio_mean), breaks = 50,
     main = "HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[delta[HHS-SI]] / V[delta[NHHS-SI]])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


active_delta_var_mean_hhssi <- delta_var_mean_hhssi[c(1:3, 18:20)]
active_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[c(1:3, 18:20)]
inactive_delta_var_mean_hhssi <- delta_var_mean_hhssi[-c(1:3, 18:20)]
inactive_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[-c(1:3, 18:20)]
active_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,c(1:3, 18:20)] /
                                          delta_var_nhhssi[,c(1:3, 18:20)])
inactive_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,-c(1:3, 18:20)] /
                                            delta_var_nhhssi[,-c(1:3, 18:20)])
active_delta_var_ratio_mean_hhssi_nhhssi <- apply(active_delta_var_ratio_hhssi_nhhssi, 2, mean)
inactive_delta_var_ratio_mean_hhssi_nhhssi <- apply(inactive_delta_var_ratio_hhssi_nhhssi, 2, mean)

hist(log10(active_delta_var_ratio_mean_hhssi_nhhssi), breaks = 5,
     main = "active HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[activedeltaHHS-SI] / V[activedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


hist(log10(inactive_delta_var_ratio_mean_hhssi_nhhssi), breaks = 50,
     main = "inactive HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[inactivedeltaHHS-SI] / V[inactivedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")









plot(beta_gamma_lam_mean_hhssi,
     pch = 19,
     log = "y",
     xlab = "Interaction index",
     ylab = expression(E[posterior](lambda[beta]^2 * lambda[gamma]^2)),
     main = "HHS-SI beta-gamma multiplier")
abline(h = 1, lty = 2, lwd = 2, col = "red")



















###############################################################################
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors real (First sample)
###############################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
# First sample
set.seed(23444)
sim_data_normal01_realpred1 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                O=O, x=X, z=Z, d=D,
                                                Sigma=Sigma)

Y <- sim_data_normal01_realpred1$Y
x <- sim_data_normal01_realpred1$x
z <- sim_data_normal01_realpred1$z
u <- sim_data_normal01_realpred1$u
d <- sim_data_normal01_realpred1$d
W <- sim_data_normal01_realpred1$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_realpred_hhssi1.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi1.rda")

sim_data_normal01_realpred1$delta_true
dim(res_normal01_realpred_hhssi1$theta_update)






# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(res_normal01_realpred_hhssi1$theta_update[, 1, k],
          type = "l",
          lty = 1,
          lwd = 2,
          xlab = "MCMC iteration",
          ylab = expression(alpha),
          main = paste("Outcome", k)
  )

  matlines(res_normal01_realpred_nhhssi1$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_realpred_hhssi1$theta_update[, 1, ])
range(res_normal01_realpred_nhhssi1$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
  hhssi_alpha <- res_normal01_realpred_hhssi1$theta_update[, 1, k]
  nhhssi_alpha <- res_normal01_realpred_nhhssi1$theta_update[, 1, k]

  plot(hhssi_alpha, nhhssi_alpha,
       pch = 19,
       col = "blue",
       xlab = "HHS-SI",
       ylab = "NHHS-SI",
       main = paste("alpha: Outcome", k),
       xlim = range(c(hhssi_alpha, nhhssi_alpha)),
       ylim = range(c(hhssi_alpha, nhhssi_alpha)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
}
par(mfrow = c(1, 1))



# bias alpha
sim_data_normal01_realpred$alpha_true
res_normal01_realpred_hhssi1$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_realpred_hhssi1$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_realpred_nhhssi1$theta_update[, 1, ], 2, mean)
# Bias for each outcome
alpha_bias_hhssi <- (alpha_est_hhssi - sim_data_normal01_realpred$alpha_true)
alpha_bias_nhhssi <- (alpha_est_nhhssi - sim_data_normal01_realpred$alpha_true)


plot(alpha_bias_hhssi,
     type = "b",
     xlab = "Outcome",
     ylab = "bias of alpha",
     lwd = 1,
     col = "red",
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
lines(alpha_bias_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(alpha_bias_hhssi, alpha_bias_nhhssi,
     lwd = 1,
     xlim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)),
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
abline(a=0, b=1, col = "red")



# active beta

plot(res_normal01_realpred_hhssi1$theta_update[, 2:3, 1],
     res_normal01_realpred_nhhssi1$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi1$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi1$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, 2:3, 2],
       res_normal01_realpred_nhhssi1$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, 2:3, 3],
       res_normal01_realpred_nhhssi1$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, 2:3, 4],
       res_normal01_realpred_nhhssi1$theta_update[, 2:3, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))




range(res_normal01_realpred_hhssi1$theta_update[, 2:3, 2])
range(res_normal01_realpred_nhhssi1$theta_update[, 2:3, 2])

range(res_normal01_realpred_hhssi1$theta_update[, 2:3, 3])
range(res_normal01_realpred_nhhssi1$theta_update[, 2:3, 3])

range(res_normal01_realpred_hhssi1$theta_update[, 2:3, 4])
range(res_normal01_realpred_nhhssi1$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, 2:3, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))




# inactive beta
range(res_normal01_realpred_hhssi1$theta_update[, 4:45, 1])
range(res_normal01_realpred_nhhssi1$theta_update[, 4:45, 1])

plot(res_normal01_realpred_hhssi1$theta_update[, 4:45, 1],
     res_normal01_realpred_nhhssi1$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi1$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi1$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, 4:45, 2],
       res_normal01_realpred_nhhssi1$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, 4:45, 3],
       res_normal01_realpred_nhhssi1$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, 4:45, 4],
       res_normal01_realpred_nhhssi1$theta_update[, 4:45, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



range(res_normal01_realpred_hhssi1$theta_update[, 4:45, 2])
range(res_normal01_realpred_nhhssi1$theta_update[, 4:45, 2])

range(res_normal01_realpred_hhssi1$theta_update[, 4:45, 3])
range(res_normal01_realpred_nhhssi1$theta_update[, 4:45, 3])

range(res_normal01_realpred_hhssi1$theta_update[, 4:45, 4])
range(res_normal01_realpred_nhhssi1$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, 4:45, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))






# bias beta
dim(sim_data_normal01_realpred$beta_true)
dim(res_normal01_realpred_hhssi1$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_realpred_hhssi1$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_realpred_nhhssi1$theta_update[, 2:45, ], c(2,3), mean)
# Bias for each outcome
beta_bias_hhssi <- (beta_est_hhssi - sim_data_normal01_realpred$beta_true)
beta_bias_nhhssi <- (beta_est_nhhssi - sim_data_normal01_realpred$beta_true)

dim(beta_bias_hhssi)
dim(beta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:44,
    beta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "modifiers",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(beta_bias_hhssi,
                   beta_bias_nhhssi)))
  lines(1:44,
        beta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  beta_bias_hhssi[, k],
  beta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(beta_bias_hhssi, beta_bias_nhhssi)),
  ylim = range(c(beta_bias_hhssi, beta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))





mean(beta_bias_hhssi>beta_bias_nhhssi)



beta_var_hhssi <- (res_normal01_realpred_hhssi1$lambdasq_beta_update *
                     res_normal01_realpred_hhssi1$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_realpred_nhhssi1$lambdasq_beta_update *
                      res_normal01_realpred_nhhssi1$tausq_beta_update)

mean(beta_var_hhssi>beta_var_nhhssi)
mean(beta_var_hhssi)
mean(beta_var_nhhssi)
mean(beta_bias_hhssi)
mean(beta_bias_nhhssi)

mean(beta_var_hhssi[1:2] > beta_var_nhhssi[1:2])
mean(beta_var_hhssi[1:2])
mean(beta_var_nhhssi[1:2])
mean(beta_bias_hhssi[1:2])
mean(beta_bias_nhhssi[1:2])

mean(beta_var_hhssi[-(1:2)] > beta_var_nhhssi[-(1:2)])
mean(beta_var_hhssi[-(1:2)])
mean(beta_var_nhhssi[-(1:2)])
mean(beta_bias_hhssi[-(1:2)])
mean(beta_bias_nhhssi[-(1:2)])



plot(beta_var_hhssi, beta_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[1:2], beta_var_nhhssi[1:2],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[-(1:2)], beta_var_nhhssi[-(1:2)],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")












# active gamma
range(res_normal01_realpred_hhssi1$theta_update[, 46:48, 1])
range(res_normal01_realpred_nhhssi1$theta_update[, 46:48, 1])


plot(res_normal01_realpred_hhssi1$theta_update[, 46:48, 1],
     res_normal01_realpred_nhhssi1$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi1$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi1$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, 46:48, 2],
       res_normal01_realpred_nhhssi1$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, 46:48, 3],
       res_normal01_realpred_nhhssi1$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, 46:48, 4],
       res_normal01_realpred_nhhssi1$theta_update[, 46:48, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_realpred_hhssi1$theta_update[, 46:48, 2])
range(res_normal01_realpred_nhhssi1$theta_update[, 46:48, 2])

range(res_normal01_realpred_hhssi1$theta_update[, 46:48, 3])
range(res_normal01_realpred_nhhssi1$theta_update[, 46:48, 3])

range(res_normal01_realpred_hhssi1$theta_update[, 46:48, 4])
range(res_normal01_realpred_nhhssi1$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, 46:48, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive gamma
range(res_normal01_realpred_hhssi1$theta_update[, 49:62, 1])
range(res_normal01_realpred_nhhssi1$theta_update[, 49:62, 1])


plot(res_normal01_realpred_hhssi1$theta_update[, 49:62, 1],
     res_normal01_realpred_nhhssi1$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi1$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi1$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, 49:62, 2],
       res_normal01_realpred_nhhssi1$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, 49:62, 3],
       res_normal01_realpred_nhhssi1$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, 49:62, 4],
       res_normal01_realpred_nhhssi1$theta_update[, 49:62, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi1$theta_update[, 49:62, 2])
range(res_normal01_realpred_nhhssi1$theta_update[, 49:62, 2])

range(res_normal01_realpred_hhssi1$theta_update[, 49:62, 3])
range(res_normal01_realpred_nhhssi1$theta_update[, 49:62, 3])

range(res_normal01_realpred_hhssi1$theta_update[, 49:62, 4])
range(res_normal01_realpred_nhhssi1$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, 49:62, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))






# bias gamma
dim(sim_data_normal01_realpred$gamma_true)
dim(res_normal01_realpred_hhssi1$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_realpred_hhssi1$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_realpred_nhhssi1$theta_update[, 46:62, ], c(2,3), mean)
# Bias for each outcome
gamma_bias_hhssi <- (gamma_est_hhssi - sim_data_normal01_realpred$gamma_true)
gamma_bias_nhhssi <- (gamma_est_nhhssi - sim_data_normal01_realpred$gamma_true)

dim(gamma_bias_hhssi)
dim(gamma_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:17,
    gamma_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "exposures",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(gamma_bias_hhssi,
                   gamma_bias_nhhssi)))
  lines(1:17,
        gamma_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  gamma_bias_hhssi[, k],
  gamma_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)),
  ylim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))



gamma_var_hhssi <- (res_normal01_realpred_hhssi1$lambdasq_gamma_update *
                      res_normal01_realpred_hhssi1$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_realpred_nhhssi1$lambdasq_gamma_update *
                       res_normal01_realpred_nhhssi1$tausq_gamma_update)

mean(gamma_var_hhssi>gamma_var_nhhssi)
mean(gamma_var_hhssi)
mean(gamma_var_nhhssi)
mean(gamma_bias_hhssi)
mean(gamma_bias_nhhssi)

mean(gamma_var_hhssi[1:2] > gamma_var_nhhssi[1:2])
mean(gamma_var_hhssi[1:2])
mean(gamma_var_nhhssi[1:2])
mean(gamma_bias_hhssi[1:2])
mean(gamma_bias_nhhssi[1:2])

mean(gamma_var_hhssi[-(1:2)] > gamma_var_nhhssi[-(1:2)])
mean(gamma_var_hhssi[-(1:2)])
mean(gamma_var_nhhssi[-(1:2)])
mean(gamma_bias_hhssi[-(1:2)])
mean(gamma_bias_nhhssi[-(1:2)])



plot(gamma_var_hhssi, gamma_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[1:3], gamma_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[-(1:3)], gamma_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")










# active delta
range(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 1],
     res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 2],
       res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 3],
       res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 4],
       res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_realpred_hhssi1$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_realpred_nhhssi1$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, c(63:65,80:82), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive delta

range(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 1],
     res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 2],
       res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 3],
       res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 4],
       res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_realpred_hhssi1$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_realpred_nhhssi1$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, c(66:79,83:810), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# bias delta
dim(sim_data_normal01_realpred$delta_true)
dim(res_normal01_realpred_hhssi1$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_realpred_hhssi1$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_realpred_nhhssi1$theta_update[, 63:810, ], c(2,3), mean)
# Bias for each outcome
delta_bias_hhssi <- (delta_est_hhssi - sim_data_normal01_realpred$delta_true)
delta_bias_nhhssi <- (delta_est_nhhssi - sim_data_normal01_realpred$delta_true)

dim(delta_bias_hhssi)
dim(delta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:748,
    delta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "interactions",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(delta_bias_hhssi,
                   delta_bias_nhhssi)))
  lines(1:748,
        delta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  delta_bias_hhssi[, k],
  delta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(delta_bias_hhssi, delta_bias_nhhssi)),
  ylim = range(c(delta_bias_hhssi, delta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






# active varphi
range(res_normal01_realpred_hhssi1$theta_update[, 811:812, 1])
range(res_normal01_realpred_nhhssi1$theta_update[, 811:812, 1])

plot(res_normal01_realpred_hhssi1$theta_update[, 811:812, 1],
     res_normal01_realpred_nhhssi1$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi1$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi1$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, 811:812, 2],
       res_normal01_realpred_nhhssi1$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, 811:812, 3],
       res_normal01_realpred_nhhssi1$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, 811:812, 4],
       res_normal01_realpred_nhhssi1$theta_update[, 811:812, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi1$theta_update[, 811:812, 2])
range(res_normal01_realpred_nhhssi1$theta_update[, 811:812, 2])

range(res_normal01_realpred_hhssi1$theta_update[, 811:812, 3])
range(res_normal01_realpred_nhhssi1$theta_update[, 811:812, 3])

range(res_normal01_realpred_hhssi1$theta_update[, 811:812, 4])
range(res_normal01_realpred_nhhssi1$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, 811:812, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive varphi
range(res_normal01_realpred_hhssi1$theta_update[, 813:816, 1])
range(res_normal01_realpred_nhhssi1$theta_update[, 813:816, 1])

plot(res_normal01_realpred_hhssi1$theta_update[, 813:816, 1],
     res_normal01_realpred_nhhssi1$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi1$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_realpred_hhssi1$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi1$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_realpred_hhssi1$theta_update[, 813:816, 2],
       res_normal01_realpred_nhhssi1$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$theta_update[, 813:816, 3],
       res_normal01_realpred_nhhssi1$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$theta_update[, 813:816, 4],
       res_normal01_realpred_nhhssi1$theta_update[, 813:816, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_realpred_hhssi1$theta_update[, 813:816, 2])
range(res_normal01_realpred_nhhssi1$theta_update[, 813:816, 2])

range(res_normal01_realpred_hhssi1$theta_update[, 813:816, 3])
range(res_normal01_realpred_nhhssi1$theta_update[, 813:816, 3])

range(res_normal01_realpred_hhssi1$theta_update[, 813:816, 4])
range(res_normal01_realpred_nhhssi1$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi1$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi1$theta_update[, 813:816, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))





# bias varphi
dim(sim_data_normal01_realpred$varphi_true)
dim(res_normal01_realpred_hhssi1$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_realpred_hhssi1$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_realpred_nhhssi1$theta_update[, 811:816, ], c(2,3), mean)
# Bias for each outcome
varphi_bias_hhssi <- (varphi_est_hhssi - sim_data_normal01_realpred$varphi_true)
varphi_bias_nhhssi <- (varphi_est_nhhssi - sim_data_normal01_realpred$varphi_true)

dim(varphi_bias_hhssi)
dim(varphi_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:6,
    varphi_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "covariates",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(varphi_bias_hhssi,
                   varphi_bias_nhhssi)))
  lines(1:6,
        varphi_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  varphi_bias_hhssi[, k],
  varphi_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)),
  ylim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






varphi_var_hhssi <- apply(res_normal01_realpred_hhssi1$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_realpred_nhhssi1$theta_update[, 811:816, ], c(2,3), var)

mean(varphi_var_hhssi>varphi_var_nhhssi)
mean(varphi_var_hhssi)
mean(varphi_var_nhhssi)
mean(varphi_bias_hhssi)
mean(varphi_bias_nhhssi)

mean(varphi_var_hhssi[1:2] > varphi_var_nhhssi[1:2])
mean(varphi_var_hhssi[1:2])
mean(varphi_var_nhhssi[1:2])
mean(varphi_bias_hhssi[1:2])
mean(varphi_bias_nhhssi[1:2])

mean(varphi_var_hhssi[-(1:2)] > varphi_var_nhhssi[-(1:2)])
mean(varphi_var_hhssi[-(1:2)])
mean(varphi_var_nhhssi[-(1:2)])
mean(varphi_bias_hhssi[-(1:2)])
mean(varphi_bias_nhhssi[-(1:2)])



plot(varphi_var_hhssi, varphi_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[1:3], varphi_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[-(1:3)], varphi_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")











dim(sim_data_normal01_realpred$gamma_true)
dim(res_normal01_realpred_hhssi1$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_realpred_hhssi1$lambdasq_delta_update)
dim(res_normal01_realpred_nhhssi1$lambdasq_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi1$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi1$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi1$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi1$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_realpred_hhssi1$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_realpred_nhhssi1$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_realpred_hhssi1$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_realpred_nhhssi1$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi1$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_realpred_nhhssi1$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_realpred_hhssi1$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_realpred_nhhssi1$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_realpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_realpred_nhhssi1$lambdasq_delta_update, 2, mean)
plot(lambdasq_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


# active deltas
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[c(1:3, 18:20)])
hist(lambdasq_delta_nhhssi[c(1:3, 18:20)])
range(lambdasq_delta_hhssi[c(1:3, 18:20)])
range(lambdasq_delta_nhhssi[c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[c(1:3, 18:20)])) # active deltas
hist(log(lambdasq_delta_nhhssi[c(1:3, 18:20)])) # active deltas

# inactive deltas
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[-c(1:3, 18:20)]) # inactive deltas
hist(lambdasq_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(lambdasq_delta_hhssi[-c(1:3, 18:20)])
range(lambdasq_delta_nhhssi[-c(1:3, 18:20)])


hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)])) # inactive deltas
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)])) # inactive deltas


# active delta
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive delta
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")





# Note: the active lambdasq_delta values are lower for hhssi model than nhhssi, and range of values is  (3.790582 1546.643649)
# the active lambdasq_delta values are higher for nhhssi model, and range of values is (154557.2 2087763.1)

# Note: whereas, the inactive lambdasq_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (2.723865e-01 1.864608e+07)
# the inactive lambdasq_delta values are lower for nhhssi model, and range of values is (1.771442 956.855391)

# Overall, we can say that lambdasq_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(lambdasq_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(lambdasq_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
lambdasq_delta_hhssi_new <- lambdasq_delta_hhssi[-remove_hhssi]
lambdasq_delta_nhhssi_new <- lambdasq_delta_nhhssi[-remove_nhhssi]
# Plot
plot(lambdasq_delta_hhssi_new, lambdasq_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
lambdasq_delta_hhssi[remove_hhssi]
lambdasq_delta_nhhssi[remove_nhhssi]




hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", lambda[delta]^2)),
     xlab = expression(lambda[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]), log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]))))
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

mean(log(lambdasq_delta_hhssi))
mean(log(lambdasq_delta_nhhssi))

mean(lambdasq_delta_hhssi)
mean(lambdasq_delta_nhhssi)

# 1:2 active betas
plot(apply(res_normal01_realpred_hhssi1$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi1$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_realpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_realpred_nhhssi1$lambdasq_beta_update, 2, mean)
plot(lambdasq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[beta]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(lambdasq_beta_hhssi, lambdasq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)),
     ylim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_beta looks larger for nhhssi than hhssi model





# 1:3 active gammas
plot(apply(res_normal01_realpred_hhssi1$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi1$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_realpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi1$lambdasq_gamma_update, 2, mean)
plot(lambdasq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[gamma]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)),
     ylim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_gamma looks similar for hhssi and nhhssi models




names(res_normal01_realpred_hhssi1)
names(res_normal01_realpred_nhhssi1)





# tausq beta
mean(res_normal01_realpred_hhssi1$tausq_beta_update)
mean(res_normal01_realpred_nhhssi1$tausq_beta_update)
abs(mean(res_normal01_realpred_hhssi1$tausq_beta_update) -
      mean(res_normal01_realpred_nhhssi1$tausq_beta_update))

range(res_normal01_realpred_hhssi1$tausq_beta_update)
range(res_normal01_realpred_nhhssi1$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_realpred_hhssi1$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_realpred_nhhssi1$tausq_beta_update
plot(tausq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[beta]^2),
     lwd = 1,
     col = "red")
lines(tausq_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_beta_hhssi, tausq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)),
     ylim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
abline(a=0, b=1, col = "red")




hist(tausq_beta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
hist(tausq_beta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))


d1 <- density(tausq_beta_hhssi)
d2 <- density(tausq_beta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)




# tausq_gamma
mean(res_normal01_realpred_hhssi1$tausq_gamma_update)
mean(res_normal01_realpred_nhhssi1$tausq_gamma_update)
abs(mean(res_normal01_realpred_hhssi1$tausq_gamma_update) -
      mean(res_normal01_realpred_nhhssi1$tausq_gamma_update))

range(res_normal01_realpred_hhssi1$tausq_gamma_update)
range(res_normal01_realpred_nhhssi1$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_realpred_hhssi1$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_realpred_nhhssi1$tausq_gamma_update
plot(tausq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[gamma]^2),
     lwd = 1,
     col = "red")
lines(tausq_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_gamma_hhssi, tausq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)),
     ylim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")






# tausq_delta
mean(res_normal01_realpred_hhssi1$tausq_delta_update)
mean(res_normal01_realpred_nhhssi1$tausq_delta_update)
abs(mean(res_normal01_realpred_hhssi1$tausq_delta_update) -
      mean(res_normal01_realpred_nhhssi1$tausq_delta_update))

range(res_normal01_realpred_hhssi1$tausq_delta_update)
range(res_normal01_realpred_nhhssi1$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_realpred_hhssi1$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi1$tausq_delta_update
plot(tausq_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[delta]^2),
     lwd = 1,
     col = "red",
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
lines(tausq_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(tausq_delta_hhssi, tausq_delta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
abline(a=0, b=1, col = "red")



hist(tausq_delta_hhssi)
hist(tausq_delta_nhhssi)

hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

# log x axis
hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     log = "x")
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))



d1 <- density(tausq_delta_hhssi)
d2 <- density(tausq_delta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)


# Note: tausq_beta is larger for hhssi than nhhssi
# tausq_gamma is same for the hhssi and nhhssi models.
# tausq_delta from hhssi and nhhssi models differs.
# tausq_delta from hhssi are higher than that of nhhssi model.






# xi beta
mean(res_normal01_realpred_hhssi1$xi_beta_update)
mean(res_normal01_realpred_nhhssi1$xi_beta_update)
abs(mean(res_normal01_realpred_hhssi1$xi_beta_update) -
      mean(res_normal01_realpred_nhhssi1$xi_beta_update))

range(res_normal01_realpred_hhssi1$xi_beta_update)
range(res_normal01_realpred_nhhssi1$xi_beta_update)

xi_beta_hhssi <- res_normal01_realpred_hhssi1$xi_beta_update
xi_beta_nhhssi <- res_normal01_realpred_nhhssi1$xi_beta_update
plot(xi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[beta]),
     lwd = 1,
     col = "red")
lines(xi_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_beta_hhssi, xi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_beta_hhssi, xi_beta_nhhssi)),
     ylim = range(c(xi_beta_hhssi, xi_beta_nhhssi)))
abline(a=0, b=1, col = "red")





# xi_gamma
mean(res_normal01_realpred_hhssi1$xi_gamma_update)
mean(res_normal01_realpred_nhhssi1$xi_gamma_update)
abs(mean(res_normal01_realpred_hhssi1$xi_gamma_update) -
      mean(res_normal01_realpred_nhhssi1$xi_gamma_update))

range(res_normal01_realpred_hhssi1$xi_gamma_update)
range(res_normal01_realpred_nhhssi1$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_realpred_hhssi1$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_realpred_nhhssi1$xi_gamma_update
plot(xi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[gamma]),
     lwd = 1,
     col = "red")
lines(xi_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_gamma_hhssi, xi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)),
     ylim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")




# xi_delta
mean(res_normal01_realpred_hhssi1$xi_delta_update)
mean(res_normal01_realpred_nhhssi1$xi_delta_update)
abs(mean(res_normal01_realpred_hhssi1$xi_delta_update) -
      mean(res_normal01_realpred_nhhssi1$xi_delta_update))

range(res_normal01_realpred_hhssi1$xi_delta_update)
range(res_normal01_realpred_nhhssi1$xi_delta_update)

xi_delta_hhssi <- res_normal01_realpred_hhssi1$xi_delta_update
xi_delta_nhhssi <- res_normal01_realpred_nhhssi1$xi_delta_update
plot(xi_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[delta]),
     lwd = 1,
     col = "red",
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
lines(xi_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_delta_hhssi, xi_delta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_delta_hhssi, xi_delta_nhhssi)),
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
abline(a=0, b=1, col = "red")

hist(xi_delta_hhssi)
hist(xi_delta_nhhssi)



# Note: xi_beta for nhhssi is larger than hhssi
# xi_gamma similar for both the hhssi and nhhssi models.
# However, xi_delta from hhssi and nhhssi models are far different.
# xi_delta from nhhssi are much higher than that of hhssi model.

##########################################################

# psi_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_realpred_nhhssi1$psi_delta_update)
dim(res_normal01_realpred_hhssi1$psi_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi1$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi1$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi1$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi1$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_realpred_hhssi1$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_realpred_nhhssi1$psi_delta_update, 2, mean)

# active
plot(psi_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive
plot(psi_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")



# active deltas
plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])
hist(psi_delta_hhssi[c(1:3, 18:20)])
range(psi_delta_hhssi[c(1:3, 18:20)])
range(psi_delta_nhhssi[c(1:3, 18:20)])


# inactive deltas
plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])
hist(psi_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(psi_delta_hhssi[-c(1:3, 18:20)])
range(psi_delta_nhhssi[-c(1:3, 18:20)])


# Note: the active psi_delta values are larger for hhssi model than nhhssi, and range of values is (4.404875 61.890418)
# the active psi_delta values are lower for nhhssi model, and range of values is (4.088125 58.551846)

# Note: whereas, the inactive psi_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (9.028653e+00 2.638961e+07)
# the inactive psi_delta values are lower for nhhssi model, and range of values is (1.640430e+01 3.505838e+08)

# Overall, we can say that psi_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(psi_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(psi_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
psi_delta_hhssi_new <- psi_delta_hhssi[-remove_hhssi]
psi_delta_nhhssi_new <- psi_delta_nhhssi[-remove_nhhssi]
# Plot
plot(psi_delta_hhssi_new, psi_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
psi_delta_hhssi[remove_hhssi]
psi_delta_nhhssi[remove_nhhssi]




hist(log(psi_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", psi[delta])),
     xlab = expression(psi[delta]),
     ylab = "Frequency",
     xlim = range(c(log(psi_delta_hhssi), log(psi_delta_nhhssi))))
hist(log(psi_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))




# 1:2 active betas
plot(apply(res_normal01_realpred_hhssi1$psi_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi1$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_realpred_hhssi1$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_realpred_nhhssi1$psi_beta_update, 2, mean)
plot(psi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[beta])),
     lwd = 1,
     col = "red")
lines(psi_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_beta_hhssi, psi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(psi_beta_hhssi, psi_beta_nhhssi)),
     ylim = range(c(psi_beta_hhssi, psi_beta_nhhssi)))
abline(a=0, b=1, col = "red")


# Note: psi_beta looks similar for hhssi and nhhssi models


# 1:3 active gammas
plot(apply(res_normal01_realpred_hhssi1$psi_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi1$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_realpred_hhssi1$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi1$psi_gamma_update, 2, mean)
plot(psi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[gamma])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
lines(psi_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_gamma_hhssi, psi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)),
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: psi_gamma looks similar for hhssi and nhhssi models


names(res_normal01_realpred_hhssi1)
names(res_normal01_realpred_nhhssi1)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_realpred_hhssi1$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_realpred_nhhssi1$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
    Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_realpred_hhssi1$Sigma_update[, , 1], type = "l")
matplot(res_normal01_realpred_nhhssi1$Sigma_update[, , 1], type = "l")

matplot(res_normal01_realpred_hhssi1$Sigma_update[, , 2], type = "l")
matplot(res_normal01_realpred_nhhssi1$Sigma_update[, , 2], type = "l")

matplot(res_normal01_realpred_hhssi1$Sigma_update[, , 3], type = "l")
matplot(res_normal01_realpred_nhhssi1$Sigma_update[, , 3], type = "l")

matplot(res_normal01_realpred_hhssi1$Sigma_update[, , 4], type = "l")
matplot(res_normal01_realpred_nhhssi1$Sigma_update[, , 4], type = "l")



plot(res_normal01_realpred_hhssi1$Sigma_update[, , 1],
     res_normal01_realpred_nhhssi1$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi1$Sigma_update,
                    res_normal01_realpred_hhssi1$Sigma_update)),
     ylim = range(c(res_normal01_realpred_hhssi1$Sigma_update,
                    res_normal01_realpred_hhssi1$Sigma_update)))
# k = 2
points(res_normal01_realpred_hhssi1$Sigma_update[, , 2],
       res_normal01_realpred_nhhssi1$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi1$Sigma_update[, , 3],
       res_normal01_realpred_nhhssi1$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi1$Sigma_update[, , 4],
       res_normal01_realpred_nhhssi1$Sigma_update[, , 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



# Note: Sigmasq are same for both hhssi and nhhssi models





###########

lambdasq_beta_hhssi
lambdasq_gamma_hhssi
lambdasq_delta_hhssi
tausq_delta_hhssi


dim(res_normal01_realpred_hhssi1$lambdasq_beta_update)
dim(res_normal01_realpred_hhssi1$lambdasq_gamma_update)
dim(res_normal01_realpred_hhssi1$lambdasq_delta_update)
length(res_normal01_realpred_hhssi1$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_realpred_hhssi1$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_realpred_hhssi1$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_realpred_hhssi1$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_realpred_hhssi1$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_realpred_nhhssi1$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi1$tausq_delta_update


n_iter <- dim(beta_lam_hhssi)[1]

# Construct beta x gamma component
beta_gamma_lam_hhssi <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_hhssi[t, col_idx] <- beta_lam_hhssi[t, j] * gamma_lam_hhssi[t, m]
      col_idx <- col_idx + 1
    }
  }
}
# prior variance
delta_var_hhssi <- delta_lam_hhssi * beta_gamma_lam_hhssi * tausq_delta_hhssi
# prior variance
delta_var_nhhssi <- delta_lam_nhhssi * tausq_delta_nhhssi



beta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi1$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi1$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi1$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_realpred_hhssi1$tausq_delta_update,
     type = "l",
     xlab = "MCMC iteration",
     ylab = expression(tau[delta]^2),
     main = "HHS-SI: global delta shrinkage"
)




hist(log10(beta_gamma_lam_hhssi),
     breaks = 50,
     main = expression("HHS-SI: " * log[10](lambda[beta]^2 * lambda[gamma]^2)),
     xlab = expression(log[10](lambda[beta]^2 * lambda[gamma]^2)))
abline(v = 0, lty = 2, lwd = 2)





delta_var_mean_hhssi <- apply(delta_var_hhssi, 2, mean)
beta_gamma_lam_mean_hhssi <- apply(beta_gamma_lam_hhssi, 2, mean)
# plot(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi,
#   log = "xy",
#   pch = 19,
#   xlab = "Delta prior var of hhssi from hhssi",
#   ylab = "Delta prior var of nhhssi from hhssi",
#   xlim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)),
#   ylim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)))
# abline(0, 1, lty = 2, lwd = 2, col = "red")


delta_var_mean_nhhssi <- apply(delta_var_nhhssi, 2, mean)
plot(delta_var_mean_hhssi, delta_var_mean_nhhssi,
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)),
     ylim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_nhhssi,
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)



plot(delta_var_mean_nhhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")






delta_var_mean_hhssi[c(1:3, 18:20)] > delta_var_mean_nhhssi[c(1:3, 18:20)]
delta_var_mean_hhssi[-c(1:3, 18:20)] > delta_var_mean_nhhssi[-c(1:3, 18:20)]


ratio_hhssi_nhhssi <- delta_var_hhssi / delta_var_nhhssi
ratio_mean <- apply(ratio_hhssi_nhhssi, 2, mean)
hist(log10(ratio_mean), breaks = 50,
     main = "HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[delta,HHS-SI] / V[delta,NHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


active_delta_var_mean_hhssi <- delta_var_mean_hhssi[c(1:3, 18:20)]
active_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[c(1:3, 18:20)]
inactive_delta_var_mean_hhssi <- delta_var_mean_hhssi[-c(1:3, 18:20)]
inactive_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[-c(1:3, 18:20)]
active_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,c(1:3, 18:20)] /
                                          delta_var_nhhssi[,c(1:3, 18:20)])
inactive_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,-c(1:3, 18:20)] /
                                            delta_var_nhhssi[,-c(1:3, 18:20)])
active_delta_var_ratio_mean_hhssi_nhhssi <- apply(active_delta_var_ratio_hhssi_nhhssi, 2, mean)
inactive_delta_var_ratio_mean_hhssi_nhhssi <- apply(inactive_delta_var_ratio_hhssi_nhhssi, 2, mean)

hist(log10(active_delta_var_ratio_mean_hhssi_nhhssi), breaks = 5,
     main = "active HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[activedeltaHHS-SI] / V[activedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


hist(log10(inactive_delta_var_ratio_mean_hhssi_nhhssi), breaks = 50,
     main = "inactive HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[inactivedeltaHHS-SI] / V[inactivedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")









plot(beta_gamma_lam_mean_hhssi,
     pch = 19,
     log = "y",
     xlab = "Interaction index",
     ylab = expression(E[posterior](lambda[beta]^2 * lambda[gamma]^2)),
     main = "HHS-SI beta-gamma multiplier")
abline(h = 1, lty = 2, lwd = 2, col = "red")









###########################################################################################
###########################################################################################


###########################################################################################
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors simulated (Second sample)
###########################################################################################

# Second sample
set.seed(654586798)
sim_data_normal01_simpred2 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)

Y <- sim_data_normal01_simpred2$Y
x <- sim_data_normal01_simpred2$x
z <- sim_data_normal01_simpred2$z
u <- sim_data_normal01_simpred2$u
d <- sim_data_normal01_simpred2$d
W <- sim_data_normal01_simpred2$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_simpred_hhssi2.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi2.rda")

sim_data_normal01_simpred2$delta_true



dim(res_normal01_simpred_hhssi2$theta_update)

# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(res_normal01_simpred_hhssi2$theta_update[, 1, k],
          type = "l",
          lty = 1,
          lwd = 2,
          xlab = "MCMC iteration",
          ylab = expression(alpha),
          main = paste("Outcome", k)
  )

  matlines(res_normal01_simpred_nhhssi2$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_simpred_hhssi2$theta_update[, 1, ])
range(res_normal01_simpred_nhhssi2$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
  hhssi_alpha <- res_normal01_simpred_hhssi2$theta_update[, 1, k]
  nhhssi_alpha <- res_normal01_simpred_nhhssi2$theta_update[, 1, k]

  plot(hhssi_alpha, nhhssi_alpha,
       pch = 19,
       col = "blue",
       xlab = "HHS-SI",
       ylab = "NHHS-SI",
       main = paste("alpha: Outcome", k),
       xlim = range(c(hhssi_alpha, nhhssi_alpha)),
       ylim = range(c(hhssi_alpha, nhhssi_alpha)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
}
par(mfrow = c(1, 1))



# bias alpha
sim_data_normal01_simpred2$alpha_true
res_normal01_simpred_hhssi2$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_simpred_hhssi2$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_simpred_nhhssi2$theta_update[, 1, ], 2, mean)
# Bias for each outcome
alpha_bias_hhssi <- (alpha_est_hhssi - sim_data_normal01_simpred2$alpha_true)
alpha_bias_nhhssi <- (alpha_est_nhhssi - sim_data_normal01_simpred2$alpha_true)


plot(alpha_bias_hhssi,
     type = "b",
     xlab = "Outcome",
     ylab = "bias of alpha",
     lwd = 1,
     col = "red",
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
lines(alpha_bias_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(alpha_bias_hhssi, alpha_bias_nhhssi,
     lwd = 1,
     xlim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)),
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
abline(a=0, b=1, col = "red")



# active beta

plot(res_normal01_simpred_hhssi2$theta_update[, 2:3, 1],
     res_normal01_simpred_nhhssi2$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi2$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi2$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, 2:3, 2],
       res_normal01_simpred_nhhssi2$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, 2:3, 3],
       res_normal01_simpred_nhhssi2$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, 2:3, 4],
       res_normal01_simpred_nhhssi2$theta_update[, 2:3, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))




range(res_normal01_simpred_hhssi2$theta_update[, 2:3, 2])
range(res_normal01_simpred_nhhssi2$theta_update[, 2:3, 2])

range(res_normal01_simpred_hhssi2$theta_update[, 2:3, 3])
range(res_normal01_simpred_nhhssi2$theta_update[, 2:3, 3])

range(res_normal01_simpred_hhssi2$theta_update[, 2:3, 4])
range(res_normal01_simpred_nhhssi2$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, 2:3, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))




# inactive beta
range(res_normal01_simpred_hhssi2$theta_update[, 4:45, 1])
range(res_normal01_simpred_nhhssi2$theta_update[, 4:45, 1])

plot(res_normal01_simpred_hhssi2$theta_update[, 4:45, 1],
     res_normal01_simpred_nhhssi2$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi2$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi2$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, 4:45, 2],
       res_normal01_simpred_nhhssi2$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, 4:45, 3],
       res_normal01_simpred_nhhssi2$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, 4:45, 4],
       res_normal01_simpred_nhhssi2$theta_update[, 4:45, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



range(res_normal01_simpred_hhssi2$theta_update[, 4:45, 2])
range(res_normal01_simpred_nhhssi2$theta_update[, 4:45, 2])

range(res_normal01_simpred_hhssi2$theta_update[, 4:45, 3])
range(res_normal01_simpred_nhhssi2$theta_update[, 4:45, 3])

range(res_normal01_simpred_hhssi2$theta_update[, 4:45, 4])
range(res_normal01_simpred_nhhssi2$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, 4:45, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))






# bias beta
dim(sim_data_normal01_simpred2$beta_true)
dim(res_normal01_simpred_hhssi2$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_simpred_hhssi2$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_simpred_nhhssi2$theta_update[, 2:45, ], c(2,3), mean)
# Bias for each outcome
beta_bias_hhssi <- (beta_est_hhssi - sim_data_normal01_simpred2$beta_true)
beta_bias_nhhssi <- (beta_est_nhhssi - sim_data_normal01_simpred2$beta_true)

dim(beta_bias_hhssi)
dim(beta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:44,
    beta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "modifiers",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(beta_bias_hhssi,
                   beta_bias_nhhssi)))
  lines(1:44,
        beta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  beta_bias_hhssi[, k],
  beta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(beta_bias_hhssi, beta_bias_nhhssi)),
  ylim = range(c(beta_bias_hhssi, beta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))





mean(beta_bias_hhssi>beta_bias_nhhssi)



beta_var_hhssi <- (res_normal01_simpred_hhssi2$lambdasq_beta_update *
                     res_normal01_simpred_hhssi2$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_simpred_nhhssi2$lambdasq_beta_update *
                      res_normal01_simpred_nhhssi2$tausq_beta_update)

mean(beta_var_hhssi>beta_var_nhhssi)
mean(beta_var_hhssi)
mean(beta_var_nhhssi)
mean(beta_bias_hhssi)
mean(beta_bias_nhhssi)

mean(beta_var_hhssi[1:2] > beta_var_nhhssi[1:2])
mean(beta_var_hhssi[1:2])
mean(beta_var_nhhssi[1:2])
mean(beta_bias_hhssi[1:2])
mean(beta_bias_nhhssi[1:2])

mean(beta_var_hhssi[-(1:2)] > beta_var_nhhssi[-(1:2)])
mean(beta_var_hhssi[-(1:2)])
mean(beta_var_nhhssi[-(1:2)])
mean(beta_bias_hhssi[-(1:2)])
mean(beta_bias_nhhssi[-(1:2)])



plot(beta_var_hhssi, beta_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[1:2], beta_var_nhhssi[1:2],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[-(1:2)], beta_var_nhhssi[-(1:2)],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")












# active gamma
range(res_normal01_simpred_hhssi2$theta_update[, 46:48, 1])
range(res_normal01_simpred_nhhssi2$theta_update[, 46:48, 1])


plot(res_normal01_simpred_hhssi2$theta_update[, 46:48, 1],
     res_normal01_simpred_nhhssi2$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi2$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi2$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, 46:48, 2],
       res_normal01_simpred_nhhssi2$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, 46:48, 3],
       res_normal01_simpred_nhhssi2$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, 46:48, 4],
       res_normal01_simpred_nhhssi2$theta_update[, 46:48, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_simpred_hhssi2$theta_update[, 46:48, 2])
range(res_normal01_simpred_nhhssi2$theta_update[, 46:48, 2])

range(res_normal01_simpred_hhssi2$theta_update[, 46:48, 3])
range(res_normal01_simpred_nhhssi2$theta_update[, 46:48, 3])

range(res_normal01_simpred_hhssi2$theta_update[, 46:48, 4])
range(res_normal01_simpred_nhhssi2$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, 46:48, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive gamma
range(res_normal01_simpred_hhssi2$theta_update[, 49:62, 1])
range(res_normal01_simpred_nhhssi2$theta_update[, 49:62, 1])


plot(res_normal01_simpred_hhssi2$theta_update[, 49:62, 1],
     res_normal01_simpred_nhhssi2$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi2$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi2$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, 49:62, 2],
       res_normal01_simpred_nhhssi2$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, 49:62, 3],
       res_normal01_simpred_nhhssi2$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, 49:62, 4],
       res_normal01_simpred_nhhssi2$theta_update[, 49:62, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi2$theta_update[, 49:62, 2])
range(res_normal01_simpred_nhhssi2$theta_update[, 49:62, 2])

range(res_normal01_simpred_hhssi2$theta_update[, 49:62, 3])
range(res_normal01_simpred_nhhssi2$theta_update[, 49:62, 3])

range(res_normal01_simpred_hhssi2$theta_update[, 49:62, 4])
range(res_normal01_simpred_nhhssi2$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, 49:62, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))






# bias gamma
dim(sim_data_normal01_simpred2$gamma_true)
dim(res_normal01_simpred_hhssi2$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_simpred_hhssi2$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_simpred_nhhssi2$theta_update[, 46:62, ], c(2,3), mean)
# Bias for each outcome
gamma_bias_hhssi <- (gamma_est_hhssi - sim_data_normal01_simpred2$gamma_true)
gamma_bias_nhhssi <- (gamma_est_nhhssi - sim_data_normal01_simpred2$gamma_true)

dim(gamma_bias_hhssi)
dim(gamma_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:17,
    gamma_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "exposures",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(gamma_bias_hhssi,
                   gamma_bias_nhhssi)))
  lines(1:17,
        gamma_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  gamma_bias_hhssi[, k],
  gamma_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)),
  ylim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))



gamma_var_hhssi <- (res_normal01_simpred_hhssi2$lambdasq_gamma_update *
                      res_normal01_simpred_hhssi2$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_simpred_nhhssi2$lambdasq_gamma_update *
                       res_normal01_simpred_nhhssi2$tausq_gamma_update)

mean(gamma_var_hhssi>gamma_var_nhhssi)
mean(gamma_var_hhssi)
mean(gamma_var_nhhssi)
mean(gamma_bias_hhssi)
mean(gamma_bias_nhhssi)

mean(gamma_var_hhssi[1:2] > gamma_var_nhhssi[1:2])
mean(gamma_var_hhssi[1:2])
mean(gamma_var_nhhssi[1:2])
mean(gamma_bias_hhssi[1:2])
mean(gamma_bias_nhhssi[1:2])

mean(gamma_var_hhssi[-(1:2)] > gamma_var_nhhssi[-(1:2)])
mean(gamma_var_hhssi[-(1:2)])
mean(gamma_var_nhhssi[-(1:2)])
mean(gamma_bias_hhssi[-(1:2)])
mean(gamma_bias_nhhssi[-(1:2)])



plot(gamma_var_hhssi, gamma_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[1:3], gamma_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[-(1:3)], gamma_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")










# active delta
range(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 1],
     res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 2],
       res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 3],
       res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 4],
       res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_simpred_hhssi2$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_simpred_nhhssi2$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, c(63:65,80:82), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive delta

range(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 1],
     res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 2],
       res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 3],
       res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 4],
       res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_simpred_hhssi2$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_simpred_nhhssi2$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, c(66:79,83:810), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# bias delta
dim(sim_data_normal01_simpred2$delta_true)
dim(res_normal01_simpred_hhssi2$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_simpred_hhssi2$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_simpred_nhhssi2$theta_update[, 63:810, ], c(2,3), mean)
# Bias for each outcome
delta_bias_hhssi <- (delta_est_hhssi - sim_data_normal01_simpred2$delta_true)
delta_bias_nhhssi <- (delta_est_nhhssi - sim_data_normal01_simpred2$delta_true)

dim(delta_bias_hhssi)
dim(delta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:748,
    delta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "interactions",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(delta_bias_hhssi,
                   delta_bias_nhhssi)))
  lines(1:748,
        delta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  delta_bias_hhssi[, k],
  delta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(delta_bias_hhssi, delta_bias_nhhssi)),
  ylim = range(c(delta_bias_hhssi, delta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






# active varphi
range(res_normal01_simpred_hhssi2$theta_update[, 811:812, 1])
range(res_normal01_simpred_nhhssi2$theta_update[, 811:812, 1])

plot(res_normal01_simpred_hhssi2$theta_update[, 811:812, 1],
     res_normal01_simpred_nhhssi2$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi2$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi2$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, 811:812, 2],
       res_normal01_simpred_nhhssi2$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, 811:812, 3],
       res_normal01_simpred_nhhssi2$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, 811:812, 4],
       res_normal01_simpred_nhhssi2$theta_update[, 811:812, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi2$theta_update[, 811:812, 2])
range(res_normal01_simpred_nhhssi2$theta_update[, 811:812, 2])

range(res_normal01_simpred_hhssi2$theta_update[, 811:812, 3])
range(res_normal01_simpred_nhhssi2$theta_update[, 811:812, 3])

range(res_normal01_simpred_hhssi2$theta_update[, 811:812, 4])
range(res_normal01_simpred_nhhssi2$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, 811:812, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive varphi
range(res_normal01_simpred_hhssi2$theta_update[, 813:816, 1])
range(res_normal01_simpred_nhhssi2$theta_update[, 813:816, 1])

plot(res_normal01_simpred_hhssi2$theta_update[, 813:816, 1],
     res_normal01_simpred_nhhssi2$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi2$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_simpred_hhssi2$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi2$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_simpred_hhssi2$theta_update[, 813:816, 2],
       res_normal01_simpred_nhhssi2$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$theta_update[, 813:816, 3],
       res_normal01_simpred_nhhssi2$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$theta_update[, 813:816, 4],
       res_normal01_simpred_nhhssi2$theta_update[, 813:816, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_simpred_hhssi2$theta_update[, 813:816, 2])
range(res_normal01_simpred_nhhssi2$theta_update[, 813:816, 2])

range(res_normal01_simpred_hhssi2$theta_update[, 813:816, 3])
range(res_normal01_simpred_nhhssi2$theta_update[, 813:816, 3])

range(res_normal01_simpred_hhssi2$theta_update[, 813:816, 4])
range(res_normal01_simpred_nhhssi2$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi2$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi2$theta_update[, 813:816, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))





# bias varphi
dim(sim_data_normal01_simpred2$varphi_true)
dim(res_normal01_simpred_hhssi2$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_simpred_hhssi2$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_simpred_nhhssi2$theta_update[, 811:816, ], c(2,3), mean)
# Bias for each outcome
varphi_bias_hhssi <- (varphi_est_hhssi - sim_data_normal01_simpred2$varphi_true)
varphi_bias_nhhssi <- (varphi_est_nhhssi - sim_data_normal01_simpred2$varphi_true)

dim(varphi_bias_hhssi)
dim(varphi_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:6,
    varphi_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "covariates",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(varphi_bias_hhssi,
                   varphi_bias_nhhssi)))
  lines(1:6,
        varphi_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  varphi_bias_hhssi[, k],
  varphi_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)),
  ylim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






varphi_var_hhssi <- apply(res_normal01_simpred_hhssi2$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_simpred_nhhssi2$theta_update[, 811:816, ], c(2,3), var)

mean(varphi_var_hhssi>varphi_var_nhhssi)
mean(varphi_var_hhssi)
mean(varphi_var_nhhssi)
mean(varphi_bias_hhssi)
mean(varphi_bias_nhhssi)

mean(varphi_var_hhssi[1:2] > varphi_var_nhhssi[1:2])
mean(varphi_var_hhssi[1:2])
mean(varphi_var_nhhssi[1:2])
mean(varphi_bias_hhssi[1:2])
mean(varphi_bias_nhhssi[1:2])

mean(varphi_var_hhssi[-(1:2)] > varphi_var_nhhssi[-(1:2)])
mean(varphi_var_hhssi[-(1:2)])
mean(varphi_var_nhhssi[-(1:2)])
mean(varphi_bias_hhssi[-(1:2)])
mean(varphi_bias_nhhssi[-(1:2)])



plot(varphi_var_hhssi, varphi_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[1:3], varphi_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[-(1:3)], varphi_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")











dim(sim_data_normal01_simpred1$gamma_true)
dim(res_normal01_simpred_hhssi2$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_simpred_hhssi2$lambdasq_delta_update)
dim(res_normal01_simpred_nhhssi2$lambdasq_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi2$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi2$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi2$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi2$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_simpred_hhssi2$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_simpred_nhhssi2$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_simpred_hhssi2$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_simpred_nhhssi2$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi2$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_simpred_nhhssi2$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_simpred_hhssi2$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_simpred_nhhssi2$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_simpred_hhssi2$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_simpred_nhhssi2$lambdasq_delta_update, 2, mean)
plot(lambdasq_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


# active deltas
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[c(1:3, 18:20)])
hist(lambdasq_delta_nhhssi[c(1:3, 18:20)])
range(lambdasq_delta_hhssi[c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[c(1:3, 18:20)])) # active deltas
hist(log(lambdasq_delta_nhhssi[c(1:3, 18:20)])) # active deltas

# inactive deltas
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[-c(1:3, 18:20)]) # inactive deltas
hist(lambdasq_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(lambdasq_delta_nhhssi[-c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)])) # inactive deltas
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)])) # inactive deltas


# active delta
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive delta
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")





# Note: the active lambdasq_delta values are lower for hhssi model than nhhssi, and range of values is (0.1485711, 2.5748392)
# the active lambdasq_delta values are higher for nhhssi model, and range of values is (58819.5, 403810.4)

# Note: whereas, the inactive lambdasq_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (1.246500e-01, 1.844546e+10)
# the inactive lambdasq_delta values are lower for nhhssi model, and range of values is (2.776089 234.788481)

# Overall, we can say that lambdasq_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(lambdasq_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(lambdasq_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
lambdasq_delta_hhssi_new <- lambdasq_delta_hhssi[-remove_hhssi]
lambdasq_delta_nhhssi_new <- lambdasq_delta_nhhssi[-remove_nhhssi]
# Plot
plot(lambdasq_delta_hhssi_new, lambdasq_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
lambdasq_delta_hhssi[remove_hhssi]
lambdasq_delta_nhhssi[remove_nhhssi]




hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", lambda[delta]^2)),
     xlab = expression(lambda[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]), log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]))))
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

mean(log(lambdasq_delta_hhssi))
mean(log(lambdasq_delta_nhhssi))

mean(lambdasq_delta_hhssi)
mean(lambdasq_delta_nhhssi)

# 1:2 active betas
plot(apply(res_normal01_simpred_hhssi2$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi2$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_simpred_hhssi2$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_simpred_nhhssi2$lambdasq_beta_update, 2, mean)
plot(lambdasq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[beta]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(lambdasq_beta_hhssi, lambdasq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)),
     ylim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_beta looks similar for hhssi and nhhssi models





# 1:3 active gammas
plot(apply(res_normal01_simpred_hhssi2$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi2$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_simpred_hhssi2$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi2$lambdasq_gamma_update, 2, mean)
plot(lambdasq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[gamma]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)),
     ylim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_gamma looks similar for hhssi and nhhssi models




names(res_normal01_simpred_hhssi2)
names(res_normal01_simpred_nhhssi2)





# tausq beta
mean(res_normal01_simpred_hhssi2$tausq_beta_update)
mean(res_normal01_simpred_nhhssi2$tausq_beta_update)
abs(mean(res_normal01_simpred_hhssi2$tausq_beta_update) -
      mean(res_normal01_simpred_nhhssi2$tausq_beta_update))

range(res_normal01_simpred_hhssi2$tausq_beta_update)
range(res_normal01_simpred_nhhssi2$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_simpred_hhssi2$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_simpred_nhhssi2$tausq_beta_update
plot(tausq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[beta]^2),
     lwd = 1,
     col = "red")
lines(tausq_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_beta_hhssi, tausq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)),
     ylim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
abline(a=0, b=1, col = "red")




hist(tausq_beta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
hist(tausq_beta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))


d1 <- density(tausq_beta_hhssi)
d2 <- density(tausq_beta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)




# tausq_gamma
mean(res_normal01_simpred_hhssi2$tausq_gamma_update)
mean(res_normal01_simpred_nhhssi2$tausq_gamma_update)
abs(mean(res_normal01_simpred_hhssi2$tausq_gamma_update) -
      mean(res_normal01_simpred_nhhssi2$tausq_gamma_update))

range(res_normal01_simpred_hhssi2$tausq_gamma_update)
range(res_normal01_simpred_nhhssi2$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_simpred_hhssi2$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_simpred_nhhssi2$tausq_gamma_update
plot(tausq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[gamma]^2),
     lwd = 1,
     col = "red")
lines(tausq_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_gamma_hhssi, tausq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)),
     ylim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")






# tausq_delta
mean(res_normal01_simpred_hhssi2$tausq_delta_update)
mean(res_normal01_simpred_nhhssi2$tausq_delta_update)
abs(mean(res_normal01_simpred_hhssi2$tausq_delta_update) -
      mean(res_normal01_simpred_nhhssi2$tausq_delta_update))

range(res_normal01_simpred_hhssi2$tausq_delta_update)
range(res_normal01_simpred_nhhssi2$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_simpred_hhssi2$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi2$tausq_delta_update
plot(tausq_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[delta]^2),
     lwd = 1,
     col = "red",
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
lines(tausq_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

hist(tausq_delta_hhssi, breaks = 50)
hist(tausq_delta_nhhssi, breaks = 50)

plot(tausq_delta_hhssi, tausq_delta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
abline(a=0, b=1, col = "red")



hist(tausq_delta_hhssi)
hist(tausq_delta_nhhssi)

hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

# log x axis
hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     log = "x")
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))



d1 <- density(tausq_delta_hhssi)
d2 <- density(tausq_delta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)


# Note: tausq_beta and tausq_gamma are same for the hhssi and nhhssi models.
# However, tausq_delta from hhssi and nhhssi models are far different.
# tausq_delta from hhssi are higher than that of nhhssi model.






# xi beta
mean(res_normal01_simpred_hhssi2$xi_beta_update)
mean(res_normal01_simpred_nhhssi2$xi_beta_update)
abs(mean(res_normal01_simpred_hhssi2$xi_beta_update) -
      mean(res_normal01_simpred_nhhssi2$xi_beta_update))

range(res_normal01_simpred_hhssi2$xi_beta_update)
range(res_normal01_simpred_nhhssi2$xi_beta_update)

xi_beta_hhssi <- res_normal01_simpred_hhssi2$xi_beta_update
xi_beta_nhhssi <- res_normal01_simpred_nhhssi2$xi_beta_update
plot(xi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[beta]),
     lwd = 1,
     col = "red")
lines(xi_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_beta_hhssi, xi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_beta_hhssi, xi_beta_nhhssi)),
     ylim = range(c(xi_beta_hhssi, xi_beta_nhhssi)))
abline(a=0, b=1, col = "red")





# xi_gamma
mean(res_normal01_simpred_hhssi2$xi_gamma_update)
mean(res_normal01_simpred_nhhssi2$xi_gamma_update)
abs(mean(res_normal01_simpred_hhssi2$xi_gamma_update) -
      mean(res_normal01_simpred_nhhssi2$xi_gamma_update))

range(res_normal01_simpred_hhssi2$xi_gamma_update)
range(res_normal01_simpred_nhhssi2$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_simpred_hhssi2$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_simpred_nhhssi2$xi_gamma_update
plot(xi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[gamma]),
     lwd = 1,
     col = "red")
lines(xi_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_gamma_hhssi, xi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)),
     ylim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")




# xi_delta
mean(res_normal01_simpred_hhssi2$xi_delta_update)
mean(res_normal01_simpred_nhhssi2$xi_delta_update)
abs(mean(res_normal01_simpred_hhssi2$xi_delta_update) -
      mean(res_normal01_simpred_nhhssi2$xi_delta_update))

range(res_normal01_simpred_hhssi2$xi_delta_update)
range(res_normal01_simpred_nhhssi2$xi_delta_update)

xi_delta_hhssi <- res_normal01_simpred_hhssi2$xi_delta_update
xi_delta_nhhssi <- res_normal01_simpred_nhhssi2$xi_delta_update
plot(xi_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[delta]),
     lwd = 1,
     col = "red",
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
lines(xi_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_delta_hhssi, xi_delta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_delta_hhssi, xi_delta_nhhssi)),
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
abline(a=0, b=1, col = "red")

hist(xi_delta_hhssi)
hist(xi_delta_nhhssi)



# Note: xi_beta and xi_gamma are same for both the hhssi and nhhssi models.
# However, xi_delta from hhssi and nhhssi models are far different.
# xi_delta from nhhssi are much higher than that of hhssi model.

##########################################################

# psi_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_simpred_nhhssi2$psi_delta_update)
dim(res_normal01_simpred_hhssi2$psi_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi2$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi2$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi2$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi2$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_simpred_hhssi2$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_simpred_nhhssi2$psi_delta_update, 2, mean)

# active
plot(psi_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
lines(psi_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive
plot(psi_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])))
lines(psi_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")



# active deltas
plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])
hist(psi_delta_hhssi[c(1:3, 18:20)])

range(psi_delta_hhssi[c(1:3, 18:20)])

# inactive deltas
plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])
hist(psi_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(psi_delta_nhhssi[-c(1:3, 18:20)])


# Note: the active psi_delta values are lower for hhssi model than nhhssi, and range of values is (0.1485711, 2.5748392)
# the active psi_delta values are higher for nhhssi model, and range of values is (58819.5, 403810.4)

# Note: whereas, the inactive psi_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (1.246500e-01, 1.844546e+10)
# the inactive psi_delta values are lower for nhhssi model, and range of values is (2.776089 234.788481)

# Overall, we can say that psi_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(psi_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(psi_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
psi_delta_hhssi_new <- psi_delta_hhssi[-remove_hhssi]
psi_delta_nhhssi_new <- psi_delta_nhhssi[-remove_nhhssi]
# Plot
plot(psi_delta_hhssi_new, psi_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
psi_delta_hhssi[remove_hhssi]
psi_delta_nhhssi[remove_nhhssi]




hist(log(psi_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", psi[delta]^2)),
     xlab = expression(psi[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(psi_delta_hhssi), log(psi_delta_nhhssi))))
hist(log(psi_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))




# 1:2 active betas
plot(apply(res_normal01_simpred_hhssi2$psi_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi2$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_simpred_hhssi2$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_simpred_nhhssi2$psi_beta_update, 2, mean)
plot(psi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[beta])),
     lwd = 1,
     col = "red")
lines(psi_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_beta_hhssi, psi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(psi_beta_hhssi, psi_beta_nhhssi)),
     ylim = range(c(psi_beta_hhssi, psi_beta_nhhssi)))
abline(a=0, b=1, col = "red")


# Note: psi_beta looks similar for hhssi and nhhssi models


# 1:3 active gammas
plot(apply(res_normal01_simpred_hhssi2$psi_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi2$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_simpred_hhssi2$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi2$psi_gamma_update, 2, mean)
plot(psi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[gamma])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
lines(psi_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_gamma_hhssi, psi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)),
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: psi_gamma looks similar for hhssi and nhhssi models


names(res_normal01_simpred_hhssi2)
names(res_normal01_simpred_nhhssi2)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_simpred_hhssi2$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_simpred_nhhssi2$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
    Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_simpred_hhssi2$Sigma_update[, , 1], type = "l")
matplot(res_normal01_simpred_nhhssi2$Sigma_update[, , 1], type = "l")

matplot(res_normal01_simpred_hhssi2$Sigma_update[, , 2], type = "l")
matplot(res_normal01_simpred_nhhssi2$Sigma_update[, , 2], type = "l")

matplot(res_normal01_simpred_hhssi2$Sigma_update[, , 3], type = "l")
matplot(res_normal01_simpred_nhhssi2$Sigma_update[, , 3], type = "l")

matplot(res_normal01_simpred_hhssi2$Sigma_update[, , 4], type = "l")
matplot(res_normal01_simpred_nhhssi2$Sigma_update[, , 4], type = "l")



plot(res_normal01_simpred_hhssi2$Sigma_update[, , 1],
     res_normal01_simpred_nhhssi2$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi2$Sigma_update,
                    res_normal01_simpred_hhssi2$Sigma_update)),
     ylim = range(c(res_normal01_simpred_hhssi2$Sigma_update,
                    res_normal01_simpred_hhssi2$Sigma_update)))
# k = 2
points(res_normal01_simpred_hhssi2$Sigma_update[, , 2],
       res_normal01_simpred_nhhssi2$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi2$Sigma_update[, , 3],
       res_normal01_simpred_nhhssi2$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi2$Sigma_update[, , 4],
       res_normal01_simpred_nhhssi2$Sigma_update[, , 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



# Note: Sigmasq are same for both hhssi and nhhssi models





###########

lambdasq_beta_hhssi
lambdasq_gamma_hhssi
lambdasq_delta_hhssi
tausq_delta_hhssi


dim(res_normal01_simpred_hhssi2$lambdasq_beta_update)
dim(res_normal01_simpred_hhssi2$lambdasq_gamma_update)
dim(res_normal01_simpred_hhssi2$lambdasq_delta_update)
length(res_normal01_simpred_hhssi2$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_simpred_hhssi2$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_simpred_hhssi2$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_simpred_hhssi2$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_simpred_hhssi2$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_simpred_nhhssi2$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi2$tausq_delta_update


n_iter <- dim(beta_lam_hhssi)[1]

# Construct beta x gamma component
beta_gamma_lam_hhssi <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_hhssi[t, col_idx] <- beta_lam_hhssi[t, j] * gamma_lam_hhssi[t, m]
      col_idx <- col_idx + 1
    }
  }
}
# prior variance
delta_var_hhssi <- delta_lam_hhssi * beta_gamma_lam_hhssi * tausq_delta_hhssi
# prior variance
delta_var_nhhssi <- delta_lam_nhhssi * tausq_delta_nhhssi



beta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi2$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi2$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi2$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_simpred_hhssi2$tausq_delta_update,
     type = "l",
     xlab = "MCMC iteration",
     ylab = expression(tau[delta]^2),
     main = "HHS-SI: global delta shrinkage"
)




hist(log10(beta_gamma_lam_hhssi),
     breaks = 50,
     main = expression("HHS-SI: " * log[10](lambda[beta]^2 * lambda[gamma]^2)),
     xlab = expression(log[10](lambda[beta]^2 * lambda[gamma]^2)))
abline(v = 0, lty = 2, lwd = 2)





delta_var_mean_hhssi <- apply(delta_var_hhssi, 2, mean)
beta_gamma_lam_mean_hhssi <- apply(beta_gamma_lam_hhssi, 2, mean)
# plot(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi,
#   log = "xy",
#   pch = 19,
#   xlab = "Delta prior var of hhssi from hhssi",
#   ylab = "Delta prior var of nhhssi from hhssi",
#   xlim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)),
#   ylim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)))
# abline(0, 1, lty = 2, lwd = 2, col = "red")


delta_var_mean_nhhssi <- apply(delta_var_nhhssi, 2, mean)
plot(delta_var_mean_hhssi, delta_var_mean_nhhssi,
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)),
     ylim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

which(delta_var_mean_hhssi[-c(1:3, 18:20)]>delta_var_mean_nhhssi[-c(1:3, 18:20)])

plot(delta_var_mean_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red",
     ylim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])))
lines(delta_var_mean_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")






delta_var_mean_hhssi[c(1:3, 18:20)] > delta_var_mean_nhhssi[c(1:3, 18:20)]
delta_var_mean_hhssi[-c(1:3, 18:20)] > delta_var_mean_nhhssi[-c(1:3, 18:20)]


ratio_hhssi_nhhssi <- delta_var_hhssi / delta_var_nhhssi
ratio_mean <- apply(ratio_hhssi_nhhssi, 2, mean)
hist(log10(ratio_mean), breaks = 50,
     main = "HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[delta,HHS-SI] / V[delta,NHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


active_delta_var_mean_hhssi <- delta_var_mean_hhssi[c(1:3, 18:20)]
active_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[c(1:3, 18:20)]
inactive_delta_var_mean_hhssi <- delta_var_mean_hhssi[-c(1:3, 18:20)]
inactive_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[-c(1:3, 18:20)]
active_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,c(1:3, 18:20)] /
                                          delta_var_nhhssi[,c(1:3, 18:20)])
inactive_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,-c(1:3, 18:20)] /
                                            delta_var_nhhssi[,-c(1:3, 18:20)])
active_delta_var_ratio_mean_hhssi_nhhssi <- apply(active_delta_var_ratio_hhssi_nhhssi, 2, mean)
inactive_delta_var_ratio_mean_hhssi_nhhssi <- apply(inactive_delta_var_ratio_hhssi_nhhssi, 2, mean)

hist(log10(active_delta_var_ratio_mean_hhssi_nhhssi), breaks = 5,
     main = "active HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[activedeltaHHS-SI] / V[activedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


hist(log10(inactive_delta_var_ratio_mean_hhssi_nhhssi), breaks = 50,
     main = "inactive HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[inactivedeltaHHS-SI] / V[inactivedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")









plot(beta_gamma_lam_mean_hhssi,
     pch = 19,
     log = "y",
     xlab = "Interaction index",
     ylab = expression(E[posterior](lambda[beta]^2 * lambda[gamma]^2)),
     main = "HHS-SI beta-gamma multiplier")
abline(h = 1, lty = 2, lwd = 2, col = "red")



















###############################################################################
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors real (Second sample)
###############################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
# Second sample
set.seed(654586798)
sim_data_normal01_realpred2 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                O=O, x=X, z=Z, d=D,
                                                Sigma=Sigma)

Y <- sim_data_normal01_realpred2$Y
x <- sim_data_normal01_realpred2$x
z <- sim_data_normal01_realpred2$z
u <- sim_data_normal01_realpred2$u
d <- sim_data_normal01_realpred2$d
W <- sim_data_normal01_realpred2$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_realpred_hhssi2.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi2.rda")

sim_data_normal01_realpred2$delta_true
dim(res_normal01_realpred_hhssi2$theta_update)






# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(res_normal01_realpred_hhssi2$theta_update[, 1, k],
          type = "l",
          lty = 1,
          lwd = 2,
          xlab = "MCMC iteration",
          ylab = expression(alpha),
          main = paste("Outcome", k)
  )

  matlines(res_normal01_realpred_nhhssi2$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_realpred_hhssi2$theta_update[, 1, ])
range(res_normal01_realpred_nhhssi2$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
  hhssi_alpha <- res_normal01_realpred_hhssi2$theta_update[, 1, k]
  nhhssi_alpha <- res_normal01_realpred_nhhssi2$theta_update[, 1, k]

  plot(hhssi_alpha, nhhssi_alpha,
       pch = 19,
       col = "blue",
       xlab = "HHS-SI",
       ylab = "NHHS-SI",
       main = paste("alpha: Outcome", k),
       xlim = range(c(hhssi_alpha, nhhssi_alpha)),
       ylim = range(c(hhssi_alpha, nhhssi_alpha)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
}
par(mfrow = c(1, 1))



# bias alpha
sim_data_normal01_realpred2$alpha_true
res_normal01_realpred_hhssi2$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_realpred_hhssi2$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_realpred_nhhssi2$theta_update[, 1, ], 2, mean)
# Bias for each outcome
alpha_bias_hhssi <- (alpha_est_hhssi - sim_data_normal01_realpred2$alpha_true)
alpha_bias_nhhssi <- (alpha_est_nhhssi - sim_data_normal01_realpred2$alpha_true)


plot(alpha_bias_hhssi,
     type = "b",
     xlab = "Outcome",
     ylab = "bias of alpha",
     lwd = 1,
     col = "red",
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
lines(alpha_bias_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(alpha_bias_hhssi, alpha_bias_nhhssi,
     lwd = 1,
     xlim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)),
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
abline(a=0, b=1, col = "red")



# active beta

plot(res_normal01_realpred_hhssi2$theta_update[, 2:3, 1],
     res_normal01_realpred_nhhssi2$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi2$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi2$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, 2:3, 2],
       res_normal01_realpred_nhhssi2$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, 2:3, 3],
       res_normal01_realpred_nhhssi2$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, 2:3, 4],
       res_normal01_realpred_nhhssi2$theta_update[, 2:3, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))




range(res_normal01_realpred_hhssi2$theta_update[, 2:3, 2])
range(res_normal01_realpred_nhhssi2$theta_update[, 2:3, 2])

range(res_normal01_realpred_hhssi2$theta_update[, 2:3, 3])
range(res_normal01_realpred_nhhssi2$theta_update[, 2:3, 3])

range(res_normal01_realpred_hhssi2$theta_update[, 2:3, 4])
range(res_normal01_realpred_nhhssi2$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, 2:3, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))




# inactive beta
range(res_normal01_realpred_hhssi2$theta_update[, 4:45, 1])
range(res_normal01_realpred_nhhssi2$theta_update[, 4:45, 1])

plot(res_normal01_realpred_hhssi2$theta_update[, 4:45, 1],
     res_normal01_realpred_nhhssi2$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi2$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi2$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, 4:45, 2],
       res_normal01_realpred_nhhssi2$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, 4:45, 3],
       res_normal01_realpred_nhhssi2$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, 4:45, 4],
       res_normal01_realpred_nhhssi2$theta_update[, 4:45, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



range(res_normal01_realpred_hhssi2$theta_update[, 4:45, 2])
range(res_normal01_realpred_nhhssi2$theta_update[, 4:45, 2])

range(res_normal01_realpred_hhssi2$theta_update[, 4:45, 3])
range(res_normal01_realpred_nhhssi2$theta_update[, 4:45, 3])

range(res_normal01_realpred_hhssi2$theta_update[, 4:45, 4])
range(res_normal01_realpred_nhhssi2$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, 4:45, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))






# bias beta
dim(sim_data_normal01_realpred2$beta_true)
dim(res_normal01_realpred_hhssi2$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_realpred_hhssi2$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_realpred_nhhssi2$theta_update[, 2:45, ], c(2,3), mean)
# Bias for each outcome
beta_bias_hhssi <- (beta_est_hhssi - sim_data_normal01_realpred2$beta_true)
beta_bias_nhhssi <- (beta_est_nhhssi - sim_data_normal01_realpred2$beta_true)

dim(beta_bias_hhssi)
dim(beta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:44,
    beta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "modifiers",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(beta_bias_hhssi,
                   beta_bias_nhhssi)))
  lines(1:44,
        beta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  beta_bias_hhssi[, k],
  beta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(beta_bias_hhssi, beta_bias_nhhssi)),
  ylim = range(c(beta_bias_hhssi, beta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))





mean(beta_bias_hhssi>beta_bias_nhhssi)



beta_var_hhssi <- (res_normal01_realpred_hhssi2$lambdasq_beta_update *
                     res_normal01_realpred_hhssi2$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_realpred_nhhssi2$lambdasq_beta_update *
                      res_normal01_realpred_nhhssi2$tausq_beta_update)

mean(beta_var_hhssi>beta_var_nhhssi)
mean(beta_var_hhssi)
mean(beta_var_nhhssi)
mean(beta_bias_hhssi)
mean(beta_bias_nhhssi)

mean(beta_var_hhssi[1:2] > beta_var_nhhssi[1:2])
mean(beta_var_hhssi[1:2])
mean(beta_var_nhhssi[1:2])
mean(beta_bias_hhssi[1:2])
mean(beta_bias_nhhssi[1:2])

mean(beta_var_hhssi[-(1:2)] > beta_var_nhhssi[-(1:2)])
mean(beta_var_hhssi[-(1:2)])
mean(beta_var_nhhssi[-(1:2)])
mean(beta_bias_hhssi[-(1:2)])
mean(beta_bias_nhhssi[-(1:2)])



plot(beta_var_hhssi, beta_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[1:2], beta_var_nhhssi[1:2],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[-(1:2)], beta_var_nhhssi[-(1:2)],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")












# active gamma
range(res_normal01_realpred_hhssi2$theta_update[, 46:48, 1])
range(res_normal01_realpred_nhhssi2$theta_update[, 46:48, 1])


plot(res_normal01_realpred_hhssi2$theta_update[, 46:48, 1],
     res_normal01_realpred_nhhssi2$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi2$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi2$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, 46:48, 2],
       res_normal01_realpred_nhhssi2$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, 46:48, 3],
       res_normal01_realpred_nhhssi2$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, 46:48, 4],
       res_normal01_realpred_nhhssi2$theta_update[, 46:48, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_realpred_hhssi2$theta_update[, 46:48, 2])
range(res_normal01_realpred_nhhssi2$theta_update[, 46:48, 2])

range(res_normal01_realpred_hhssi2$theta_update[, 46:48, 3])
range(res_normal01_realpred_nhhssi2$theta_update[, 46:48, 3])

range(res_normal01_realpred_hhssi2$theta_update[, 46:48, 4])
range(res_normal01_realpred_nhhssi2$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, 46:48, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive gamma
range(res_normal01_realpred_hhssi2$theta_update[, 49:62, 1])
range(res_normal01_realpred_nhhssi2$theta_update[, 49:62, 1])


plot(res_normal01_realpred_hhssi2$theta_update[, 49:62, 1],
     res_normal01_realpred_nhhssi2$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi2$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi2$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, 49:62, 2],
       res_normal01_realpred_nhhssi2$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, 49:62, 3],
       res_normal01_realpred_nhhssi2$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, 49:62, 4],
       res_normal01_realpred_nhhssi2$theta_update[, 49:62, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi2$theta_update[, 49:62, 2])
range(res_normal01_realpred_nhhssi2$theta_update[, 49:62, 2])

range(res_normal01_realpred_hhssi2$theta_update[, 49:62, 3])
range(res_normal01_realpred_nhhssi2$theta_update[, 49:62, 3])

range(res_normal01_realpred_hhssi2$theta_update[, 49:62, 4])
range(res_normal01_realpred_nhhssi2$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, 49:62, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))






# bias gamma
dim(sim_data_normal01_realpred2$gamma_true)
dim(res_normal01_realpred_hhssi2$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_realpred_hhssi2$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_realpred_nhhssi2$theta_update[, 46:62, ], c(2,3), mean)
# Bias for each outcome
gamma_bias_hhssi <- (gamma_est_hhssi - sim_data_normal01_realpred2$gamma_true)
gamma_bias_nhhssi <- (gamma_est_nhhssi - sim_data_normal01_realpred2$gamma_true)

dim(gamma_bias_hhssi)
dim(gamma_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:17,
    gamma_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "exposures",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(gamma_bias_hhssi,
                   gamma_bias_nhhssi)))
  lines(1:17,
        gamma_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  gamma_bias_hhssi[, k],
  gamma_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)),
  ylim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))



gamma_var_hhssi <- (res_normal01_realpred_hhssi2$lambdasq_gamma_update *
                      res_normal01_realpred_hhssi2$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_realpred_nhhssi2$lambdasq_gamma_update *
                       res_normal01_realpred_nhhssi2$tausq_gamma_update)

mean(gamma_var_hhssi>gamma_var_nhhssi)
mean(gamma_var_hhssi)
mean(gamma_var_nhhssi)
mean(gamma_bias_hhssi)
mean(gamma_bias_nhhssi)

mean(gamma_var_hhssi[1:2] > gamma_var_nhhssi[1:2])
mean(gamma_var_hhssi[1:2])
mean(gamma_var_nhhssi[1:2])
mean(gamma_bias_hhssi[1:2])
mean(gamma_bias_nhhssi[1:2])

mean(gamma_var_hhssi[-(1:2)] > gamma_var_nhhssi[-(1:2)])
mean(gamma_var_hhssi[-(1:2)])
mean(gamma_var_nhhssi[-(1:2)])
mean(gamma_bias_hhssi[-(1:2)])
mean(gamma_bias_nhhssi[-(1:2)])



plot(gamma_var_hhssi, gamma_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[1:3], gamma_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[-(1:3)], gamma_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")










# active delta
range(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 1],
     res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 2],
       res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 3],
       res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 4],
       res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_realpred_hhssi2$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_realpred_nhhssi2$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, c(63:65,80:82), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive delta

range(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 1],
     res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 2],
       res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 3],
       res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 4],
       res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_realpred_hhssi2$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_realpred_nhhssi2$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, c(66:79,83:810), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# bias delta
dim(sim_data_normal01_realpred2$delta_true)
dim(res_normal01_realpred_hhssi2$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_realpred_hhssi2$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_realpred_nhhssi2$theta_update[, 63:810, ], c(2,3), mean)
# Bias for each outcome
delta_bias_hhssi <- (delta_est_hhssi - sim_data_normal01_realpred2$delta_true)
delta_bias_nhhssi <- (delta_est_nhhssi - sim_data_normal01_realpred2$delta_true)

dim(delta_bias_hhssi)
dim(delta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:748,
    delta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "interactions",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(delta_bias_hhssi,
                   delta_bias_nhhssi)))
  lines(1:748,
        delta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  delta_bias_hhssi[, k],
  delta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(delta_bias_hhssi, delta_bias_nhhssi)),
  ylim = range(c(delta_bias_hhssi, delta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






# active varphi
range(res_normal01_realpred_hhssi2$theta_update[, 811:812, 1])
range(res_normal01_realpred_nhhssi2$theta_update[, 811:812, 1])

plot(res_normal01_realpred_hhssi2$theta_update[, 811:812, 1],
     res_normal01_realpred_nhhssi2$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi2$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi2$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, 811:812, 2],
       res_normal01_realpred_nhhssi2$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, 811:812, 3],
       res_normal01_realpred_nhhssi2$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, 811:812, 4],
       res_normal01_realpred_nhhssi2$theta_update[, 811:812, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi2$theta_update[, 811:812, 2])
range(res_normal01_realpred_nhhssi2$theta_update[, 811:812, 2])

range(res_normal01_realpred_hhssi2$theta_update[, 811:812, 3])
range(res_normal01_realpred_nhhssi2$theta_update[, 811:812, 3])

range(res_normal01_realpred_hhssi2$theta_update[, 811:812, 4])
range(res_normal01_realpred_nhhssi2$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, 811:812, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive varphi
range(res_normal01_realpred_hhssi2$theta_update[, 813:816, 1])
range(res_normal01_realpred_nhhssi2$theta_update[, 813:816, 1])

plot(res_normal01_realpred_hhssi2$theta_update[, 813:816, 1],
     res_normal01_realpred_nhhssi2$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi2$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_realpred_hhssi2$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi2$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_realpred_hhssi2$theta_update[, 813:816, 2],
       res_normal01_realpred_nhhssi2$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$theta_update[, 813:816, 3],
       res_normal01_realpred_nhhssi2$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$theta_update[, 813:816, 4],
       res_normal01_realpred_nhhssi2$theta_update[, 813:816, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_realpred_hhssi2$theta_update[, 813:816, 2])
range(res_normal01_realpred_nhhssi2$theta_update[, 813:816, 2])

range(res_normal01_realpred_hhssi2$theta_update[, 813:816, 3])
range(res_normal01_realpred_nhhssi2$theta_update[, 813:816, 3])

range(res_normal01_realpred_hhssi2$theta_update[, 813:816, 4])
range(res_normal01_realpred_nhhssi2$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi2$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi2$theta_update[, 813:816, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))





# bias varphi
dim(sim_data_normal01_realpred2$varphi_true)
dim(res_normal01_realpred_hhssi2$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_realpred_hhssi2$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_realpred_nhhssi2$theta_update[, 811:816, ], c(2,3), mean)
# Bias for each outcome
varphi_bias_hhssi <- (varphi_est_hhssi - sim_data_normal01_realpred2$varphi_true)
varphi_bias_nhhssi <- (varphi_est_nhhssi - sim_data_normal01_realpred2$varphi_true)

dim(varphi_bias_hhssi)
dim(varphi_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:6,
    varphi_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "covariates",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(varphi_bias_hhssi,
                   varphi_bias_nhhssi)))
  lines(1:6,
        varphi_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  varphi_bias_hhssi[, k],
  varphi_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)),
  ylim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






varphi_var_hhssi <- apply(res_normal01_realpred_hhssi2$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_realpred_nhhssi2$theta_update[, 811:816, ], c(2,3), var)

mean(varphi_var_hhssi>varphi_var_nhhssi)
mean(varphi_var_hhssi)
mean(varphi_var_nhhssi)
mean(varphi_bias_hhssi)
mean(varphi_bias_nhhssi)

mean(varphi_var_hhssi[1:2] > varphi_var_nhhssi[1:2])
mean(varphi_var_hhssi[1:2])
mean(varphi_var_nhhssi[1:2])
mean(varphi_bias_hhssi[1:2])
mean(varphi_bias_nhhssi[1:2])

mean(varphi_var_hhssi[-(1:2)] > varphi_var_nhhssi[-(1:2)])
mean(varphi_var_hhssi[-(1:2)])
mean(varphi_var_nhhssi[-(1:2)])
mean(varphi_bias_hhssi[-(1:2)])
mean(varphi_bias_nhhssi[-(1:2)])



plot(varphi_var_hhssi, varphi_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[1:3], varphi_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[-(1:3)], varphi_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")











dim(sim_data_normal01_realpred2$gamma_true)
dim(res_normal01_realpred_hhssi2$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_realpred_hhssi2$lambdasq_delta_update)
dim(res_normal01_realpred_nhhssi2$lambdasq_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi2$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi2$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi2$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi2$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_realpred_hhssi2$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_realpred_nhhssi2$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_realpred_hhssi2$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_realpred_nhhssi2$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi2$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_realpred_nhhssi2$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_realpred_hhssi2$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_realpred_nhhssi2$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_realpred_hhssi2$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_realpred_nhhssi2$lambdasq_delta_update, 2, mean)
plot(lambdasq_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


# active deltas
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[c(1:3, 18:20)])
hist(lambdasq_delta_nhhssi[c(1:3, 18:20)])
range(lambdasq_delta_hhssi[c(1:3, 18:20)])
range(lambdasq_delta_nhhssi[c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[c(1:3, 18:20)])) # active deltas
hist(log(lambdasq_delta_nhhssi[c(1:3, 18:20)])) # active deltas

# inactive deltas
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[-c(1:3, 18:20)]) # inactive deltas
hist(lambdasq_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(lambdasq_delta_hhssi[-c(1:3, 18:20)])
range(lambdasq_delta_nhhssi[-c(1:3, 18:20)])


hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)])) # inactive deltas
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)])) # inactive deltas


# active delta
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive delta
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")





# Note: the active lambdasq_delta values are lower for hhssi model than nhhssi, and range of values is  (3.790582 1546.643649)
# the active lambdasq_delta values are higher for nhhssi model, and range of values is (154557.2 2087763.1)

# Note: whereas, the inactive lambdasq_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (2.723865e-01 1.864608e+07)
# the inactive lambdasq_delta values are lower for nhhssi model, and range of values is (1.771442 956.855391)

# Overall, we can say that lambdasq_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(lambdasq_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(lambdasq_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
lambdasq_delta_hhssi_new <- lambdasq_delta_hhssi[-remove_hhssi]
lambdasq_delta_nhhssi_new <- lambdasq_delta_nhhssi[-remove_nhhssi]
# Plot
plot(lambdasq_delta_hhssi_new, lambdasq_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
lambdasq_delta_hhssi[remove_hhssi]
lambdasq_delta_nhhssi[remove_nhhssi]




hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", lambda[delta]^2)),
     xlab = expression(lambda[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]), log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]))))
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

mean(log(lambdasq_delta_hhssi))
mean(log(lambdasq_delta_nhhssi))

mean(lambdasq_delta_hhssi)
mean(lambdasq_delta_nhhssi)

# 1:2 active betas
plot(apply(res_normal01_realpred_hhssi2$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi2$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_realpred_hhssi2$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_realpred_nhhssi2$lambdasq_beta_update, 2, mean)
plot(lambdasq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[beta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
lines(lambdasq_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(lambdasq_beta_hhssi, lambdasq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)),
     ylim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_beta looks larger for nhhssi than hhssi model





# 1:3 active gammas
plot(apply(res_normal01_realpred_hhssi2$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi2$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_realpred_hhssi2$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi2$lambdasq_gamma_update, 2, mean)
plot(lambdasq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[gamma]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)),
     ylim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_gamma looks similar for hhssi and nhhssi models




names(res_normal01_realpred_hhssi2)
names(res_normal01_realpred_nhhssi2)





# tausq beta
mean(res_normal01_realpred_hhssi2$tausq_beta_update)
mean(res_normal01_realpred_nhhssi2$tausq_beta_update)
abs(mean(res_normal01_realpred_hhssi2$tausq_beta_update) -
      mean(res_normal01_realpred_nhhssi2$tausq_beta_update))

range(res_normal01_realpred_hhssi2$tausq_beta_update)
range(res_normal01_realpred_nhhssi2$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_realpred_hhssi2$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_realpred_nhhssi2$tausq_beta_update
plot(tausq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[beta]^2),
     lwd = 1,
     col = "red")
lines(tausq_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_beta_hhssi, tausq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)),
     ylim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
abline(a=0, b=1, col = "red")




hist(tausq_beta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
hist(tausq_beta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))


d1 <- density(tausq_beta_hhssi)
d2 <- density(tausq_beta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)




# tausq_gamma
mean(res_normal01_realpred_hhssi2$tausq_gamma_update)
mean(res_normal01_realpred_nhhssi2$tausq_gamma_update)
abs(mean(res_normal01_realpred_hhssi2$tausq_gamma_update) -
      mean(res_normal01_realpred_nhhssi2$tausq_gamma_update))

range(res_normal01_realpred_hhssi2$tausq_gamma_update)
range(res_normal01_realpred_nhhssi2$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_realpred_hhssi2$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_realpred_nhhssi2$tausq_gamma_update
plot(tausq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[gamma]^2),
     lwd = 1,
     col = "red")
lines(tausq_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_gamma_hhssi, tausq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)),
     ylim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")






# tausq_delta
mean(res_normal01_realpred_hhssi2$tausq_delta_update)
mean(res_normal01_realpred_nhhssi2$tausq_delta_update)
abs(mean(res_normal01_realpred_hhssi2$tausq_delta_update) -
      mean(res_normal01_realpred_nhhssi2$tausq_delta_update))

range(res_normal01_realpred_hhssi2$tausq_delta_update)
range(res_normal01_realpred_nhhssi2$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_realpred_hhssi2$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi2$tausq_delta_update
plot(tausq_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[delta]^2),
     lwd = 1,
     col = "red",
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
lines(tausq_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(tausq_delta_hhssi, tausq_delta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
abline(a=0, b=1, col = "red")



hist(tausq_delta_hhssi)
hist(tausq_delta_nhhssi)

hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

# log x axis
hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     log = "x")
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))



d1 <- density(tausq_delta_hhssi)
d2 <- density(tausq_delta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)


# Note: tausq_beta is larger for hhssi than nhhssi
# tausq_gamma is same for the hhssi and nhhssi models.
# tausq_delta from hhssi and nhhssi models differs.
# tausq_delta from hhssi are higher than that of nhhssi model.






# xi beta
mean(res_normal01_realpred_hhssi2$xi_beta_update)
mean(res_normal01_realpred_nhhssi2$xi_beta_update)
abs(mean(res_normal01_realpred_hhssi2$xi_beta_update) -
      mean(res_normal01_realpred_nhhssi2$xi_beta_update))

range(res_normal01_realpred_hhssi2$xi_beta_update)
range(res_normal01_realpred_nhhssi2$xi_beta_update)

xi_beta_hhssi <- res_normal01_realpred_hhssi2$xi_beta_update
xi_beta_nhhssi <- res_normal01_realpred_nhhssi2$xi_beta_update
plot(xi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[beta]),
     lwd = 1,
     col = "red",
     ylim=range(c(xi_beta_hhssi, xi_beta_nhhssi)))
lines(xi_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_beta_hhssi, xi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_beta_hhssi, xi_beta_nhhssi)),
     ylim = range(c(xi_beta_hhssi, xi_beta_nhhssi)))
abline(a=0, b=1, col = "red")





# xi_gamma
mean(res_normal01_realpred_hhssi2$xi_gamma_update)
mean(res_normal01_realpred_nhhssi2$xi_gamma_update)
abs(mean(res_normal01_realpred_hhssi2$xi_gamma_update) -
      mean(res_normal01_realpred_nhhssi2$xi_gamma_update))

range(res_normal01_realpred_hhssi2$xi_gamma_update)
range(res_normal01_realpred_nhhssi2$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_realpred_hhssi2$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_realpred_nhhssi2$xi_gamma_update
plot(xi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[gamma]),
     lwd = 1,
     col = "red")
lines(xi_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_gamma_hhssi, xi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)),
     ylim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")




# xi_delta
mean(res_normal01_realpred_hhssi2$xi_delta_update)
mean(res_normal01_realpred_nhhssi2$xi_delta_update)
abs(mean(res_normal01_realpred_hhssi2$xi_delta_update) -
      mean(res_normal01_realpred_nhhssi2$xi_delta_update))

range(res_normal01_realpred_hhssi2$xi_delta_update)
range(res_normal01_realpred_nhhssi2$xi_delta_update)

xi_delta_hhssi <- res_normal01_realpred_hhssi2$xi_delta_update
xi_delta_nhhssi <- res_normal01_realpred_nhhssi2$xi_delta_update
plot(xi_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[delta]),
     lwd = 1,
     col = "red",
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
lines(xi_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_delta_hhssi, xi_delta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_delta_hhssi, xi_delta_nhhssi)),
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
abline(a=0, b=1, col = "red")

hist(xi_delta_hhssi)
hist(xi_delta_nhhssi)



# Note: xi_beta for nhhssi is larger than hhssi
# xi_gamma similar for both the hhssi and nhhssi models.
# However, xi_delta from hhssi and nhhssi models are far different.
# xi_delta from nhhssi are much higher than that of hhssi model.

##########################################################

# psi_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_realpred_nhhssi2$psi_delta_update)
dim(res_normal01_realpred_hhssi2$psi_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi2$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi2$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi2$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi2$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_realpred_hhssi2$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_realpred_nhhssi2$psi_delta_update, 2, mean)

# active
plot(psi_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red",
     ylim=range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
lines(psi_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive
plot(psi_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")



# active deltas
plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])
hist(psi_delta_hhssi[c(1:3, 18:20)])
range(psi_delta_hhssi[c(1:3, 18:20)])
range(psi_delta_nhhssi[c(1:3, 18:20)])


# inactive deltas
plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])
hist(psi_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(psi_delta_hhssi[-c(1:3, 18:20)])
range(psi_delta_nhhssi[-c(1:3, 18:20)])


# Note: the active psi_delta values are larger for hhssi model than nhhssi, and range of values is (4.404875 61.890418)
# the active psi_delta values are lower for nhhssi model, and range of values is (4.088125 58.551846)

# Note: whereas, the inactive psi_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (9.028653e+00 2.638961e+07)
# the inactive psi_delta values are lower for nhhssi model, and range of values is (1.640430e+01 3.505838e+08)

# Overall, we can say that psi_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(psi_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(psi_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
psi_delta_hhssi_new <- psi_delta_hhssi[-remove_hhssi]
psi_delta_nhhssi_new <- psi_delta_nhhssi[-remove_nhhssi]
# Plot
plot(psi_delta_hhssi_new, psi_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
psi_delta_hhssi[remove_hhssi]
psi_delta_nhhssi[remove_nhhssi]




hist(log(psi_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", psi[delta])),
     xlab = expression(psi[delta]),
     ylab = "Frequency",
     xlim = range(c(log(psi_delta_hhssi), log(psi_delta_nhhssi))))
hist(log(psi_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))




# 1:2 active betas
plot(apply(res_normal01_realpred_hhssi2$psi_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi2$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_realpred_hhssi2$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_realpred_nhhssi2$psi_beta_update, 2, mean)
plot(psi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[beta])),
     lwd = 1,
     col = "red",
     ylim=range(c(psi_beta_hhssi, psi_beta_nhhssi)))
lines(psi_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_beta_hhssi, psi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(psi_beta_hhssi, psi_beta_nhhssi)),
     ylim = range(c(psi_beta_hhssi, psi_beta_nhhssi)))
abline(a=0, b=1, col = "red")


# Note: psi_beta looks similar for hhssi and nhhssi models


# 1:3 active gammas
plot(apply(res_normal01_realpred_hhssi2$psi_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi2$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_realpred_hhssi2$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi2$psi_gamma_update, 2, mean)
plot(psi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[gamma])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
lines(psi_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_gamma_hhssi, psi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)),
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: psi_gamma looks similar for hhssi and nhhssi models


names(res_normal01_realpred_hhssi2)
names(res_normal01_realpred_nhhssi2)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_realpred_hhssi2$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_realpred_nhhssi2$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
    Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_realpred_hhssi2$Sigma_update[, , 1], type = "l")
matplot(res_normal01_realpred_nhhssi2$Sigma_update[, , 1], type = "l")

matplot(res_normal01_realpred_hhssi2$Sigma_update[, , 2], type = "l")
matplot(res_normal01_realpred_nhhssi2$Sigma_update[, , 2], type = "l")

matplot(res_normal01_realpred_hhssi2$Sigma_update[, , 3], type = "l")
matplot(res_normal01_realpred_nhhssi2$Sigma_update[, , 3], type = "l")

matplot(res_normal01_realpred_hhssi2$Sigma_update[, , 4], type = "l")
matplot(res_normal01_realpred_nhhssi2$Sigma_update[, , 4], type = "l")



plot(res_normal01_realpred_hhssi2$Sigma_update[, , 1],
     res_normal01_realpred_nhhssi2$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi2$Sigma_update,
                    res_normal01_realpred_hhssi2$Sigma_update)),
     ylim = range(c(res_normal01_realpred_hhssi2$Sigma_update,
                    res_normal01_realpred_hhssi2$Sigma_update)))
# k = 2
points(res_normal01_realpred_hhssi2$Sigma_update[, , 2],
       res_normal01_realpred_nhhssi2$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi2$Sigma_update[, , 3],
       res_normal01_realpred_nhhssi2$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi2$Sigma_update[, , 4],
       res_normal01_realpred_nhhssi2$Sigma_update[, , 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



# Note: Sigmasq are same for both hhssi and nhhssi models





###########

lambdasq_beta_hhssi
lambdasq_gamma_hhssi
lambdasq_delta_hhssi
tausq_delta_hhssi


dim(res_normal01_realpred_hhssi2$lambdasq_beta_update)
dim(res_normal01_realpred_hhssi2$lambdasq_gamma_update)
dim(res_normal01_realpred_hhssi2$lambdasq_delta_update)
length(res_normal01_realpred_hhssi2$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_realpred_hhssi2$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_realpred_hhssi2$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_realpred_hhssi2$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_realpred_hhssi2$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_realpred_nhhssi2$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi2$tausq_delta_update


n_iter <- dim(beta_lam_hhssi)[1]

# Construct beta x gamma component
beta_gamma_lam_hhssi <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_hhssi[t, col_idx] <- beta_lam_hhssi[t, j] * gamma_lam_hhssi[t, m]
      col_idx <- col_idx + 1
    }
  }
}
# prior variance
delta_var_hhssi <- delta_lam_hhssi * beta_gamma_lam_hhssi * tausq_delta_hhssi
# prior variance
delta_var_nhhssi <- delta_lam_nhhssi * tausq_delta_nhhssi



beta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi2$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi2$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi2$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_realpred_hhssi2$tausq_delta_update,
     type = "l",
     xlab = "MCMC iteration",
     ylab = expression(tau[delta]^2),
     main = "HHS-SI: global delta shrinkage"
)




hist(log10(beta_gamma_lam_hhssi),
     breaks = 50,
     main = expression("HHS-SI: " * log[10](lambda[beta]^2 * lambda[gamma]^2)),
     xlab = expression(log[10](lambda[beta]^2 * lambda[gamma]^2)))
abline(v = 0, lty = 2, lwd = 2)





delta_var_mean_hhssi <- apply(delta_var_hhssi, 2, mean)
beta_gamma_lam_mean_hhssi <- apply(beta_gamma_lam_hhssi, 2, mean)
# plot(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi,
#   log = "xy",
#   pch = 19,
#   xlab = "Delta prior var of hhssi from hhssi",
#   ylab = "Delta prior var of nhhssi from hhssi",
#   xlim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)),
#   ylim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)))
# abline(0, 1, lty = 2, lwd = 2, col = "red")


delta_var_mean_nhhssi <- apply(delta_var_nhhssi, 2, mean)
plot(delta_var_mean_hhssi, delta_var_mean_nhhssi,
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)),
     ylim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_nhhssi,
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)



plot(delta_var_mean_nhhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")






delta_var_mean_hhssi[c(1:3, 18:20)] > delta_var_mean_nhhssi[c(1:3, 18:20)]
delta_var_mean_hhssi[-c(1:3, 18:20)] > delta_var_mean_nhhssi[-c(1:3, 18:20)]


ratio_hhssi_nhhssi <- delta_var_hhssi / delta_var_nhhssi
ratio_mean <- apply(ratio_hhssi_nhhssi, 2, mean)
hist(log10(ratio_mean), breaks = 50,
     main = "HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[delta,HHS-SI] / V[delta,NHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


active_delta_var_mean_hhssi <- delta_var_mean_hhssi[c(1:3, 18:20)]
active_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[c(1:3, 18:20)]
inactive_delta_var_mean_hhssi <- delta_var_mean_hhssi[-c(1:3, 18:20)]
inactive_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[-c(1:3, 18:20)]
active_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,c(1:3, 18:20)] /
                                          delta_var_nhhssi[,c(1:3, 18:20)])
inactive_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,-c(1:3, 18:20)] /
                                            delta_var_nhhssi[,-c(1:3, 18:20)])
active_delta_var_ratio_mean_hhssi_nhhssi <- apply(active_delta_var_ratio_hhssi_nhhssi, 2, mean)
inactive_delta_var_ratio_mean_hhssi_nhhssi <- apply(inactive_delta_var_ratio_hhssi_nhhssi, 2, mean)

hist(log10(active_delta_var_ratio_mean_hhssi_nhhssi), breaks = 5,
     main = "active HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[activedeltaHHS-SI] / V[activedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


hist(log10(inactive_delta_var_ratio_mean_hhssi_nhhssi), breaks = 50,
     main = "inactive HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[inactivedeltaHHS-SI] / V[inactivedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")









plot(beta_gamma_lam_mean_hhssi,
     pch = 19,
     log = "y",
     xlab = "Interaction index",
     ylab = expression(E[posterior](lambda[beta]^2 * lambda[gamma]^2)),
     main = "HHS-SI beta-gamma multiplier")
abline(h = 1, lty = 2, lwd = 2, col = "red")




















#####################################################################################
#######################################################################################

###########################################################################################
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors simulated (Third sample)
###########################################################################################

# Third sample
set.seed(23315644)
sim_data_normal01_simpred3 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)

Y <- sim_data_normal01_simpred3$Y
x <- sim_data_normal01_simpred3$x
z <- sim_data_normal01_simpred3$z
u <- sim_data_normal01_simpred3$u
d <- sim_data_normal01_simpred3$d
W <- sim_data_normal01_simpred3$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_simpred_hhssi3.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi3.rda")

sim_data_normal01_simpred3$delta_true



dim(res_normal01_simpred_hhssi3$theta_update)

# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(res_normal01_simpred_hhssi3$theta_update[, 1, k],
          type = "l",
          lty = 1,
          lwd = 2,
          xlab = "MCMC iteration",
          ylab = expression(alpha),
          main = paste("Outcome", k)
  )

  matlines(res_normal01_simpred_nhhssi3$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_simpred_hhssi3$theta_update[, 1, ])
range(res_normal01_simpred_nhhssi3$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
  hhssi_alpha <- res_normal01_simpred_hhssi3$theta_update[, 1, k]
  nhhssi_alpha <- res_normal01_simpred_nhhssi3$theta_update[, 1, k]

  plot(hhssi_alpha, nhhssi_alpha,
       pch = 19,
       col = "blue",
       xlab = "HHS-SI",
       ylab = "NHHS-SI",
       main = paste("alpha: Outcome", k),
       xlim = range(c(hhssi_alpha, nhhssi_alpha)),
       ylim = range(c(hhssi_alpha, nhhssi_alpha)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
}
par(mfrow = c(1, 1))



# bias alpha
sim_data_normal01_simpred3$alpha_true
res_normal01_simpred_hhssi3$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_simpred_hhssi3$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_simpred_nhhssi3$theta_update[, 1, ], 2, mean)
# Bias for each outcome
alpha_bias_hhssi <- (alpha_est_hhssi - sim_data_normal01_simpred3$alpha_true)
alpha_bias_nhhssi <- (alpha_est_nhhssi - sim_data_normal01_simpred3$alpha_true)


plot(alpha_bias_hhssi,
     type = "b",
     xlab = "Outcome",
     ylab = "bias of alpha",
     lwd = 1,
     col = "red",
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
lines(alpha_bias_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(alpha_bias_hhssi, alpha_bias_nhhssi,
     lwd = 1,
     xlim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)),
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
abline(a=0, b=1, col = "red")



# active beta

plot(res_normal01_simpred_hhssi3$theta_update[, 2:3, 1],
     res_normal01_simpred_nhhssi3$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi3$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi3$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, 2:3, 2],
       res_normal01_simpred_nhhssi3$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, 2:3, 3],
       res_normal01_simpred_nhhssi3$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, 2:3, 4],
       res_normal01_simpred_nhhssi3$theta_update[, 2:3, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))




range(res_normal01_simpred_hhssi3$theta_update[, 2:3, 2])
range(res_normal01_simpred_nhhssi3$theta_update[, 2:3, 2])

range(res_normal01_simpred_hhssi3$theta_update[, 2:3, 3])
range(res_normal01_simpred_nhhssi3$theta_update[, 2:3, 3])

range(res_normal01_simpred_hhssi3$theta_update[, 2:3, 4])
range(res_normal01_simpred_nhhssi3$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, 2:3, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))




# inactive beta
range(res_normal01_simpred_hhssi3$theta_update[, 4:45, 1])
range(res_normal01_simpred_nhhssi3$theta_update[, 4:45, 1])

plot(res_normal01_simpred_hhssi3$theta_update[, 4:45, 1],
     res_normal01_simpred_nhhssi3$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi3$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi3$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, 4:45, 2],
       res_normal01_simpred_nhhssi3$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, 4:45, 3],
       res_normal01_simpred_nhhssi3$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, 4:45, 4],
       res_normal01_simpred_nhhssi3$theta_update[, 4:45, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



range(res_normal01_simpred_hhssi3$theta_update[, 4:45, 2])
range(res_normal01_simpred_nhhssi3$theta_update[, 4:45, 2])

range(res_normal01_simpred_hhssi3$theta_update[, 4:45, 3])
range(res_normal01_simpred_nhhssi3$theta_update[, 4:45, 3])

range(res_normal01_simpred_hhssi3$theta_update[, 4:45, 4])
range(res_normal01_simpred_nhhssi3$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, 4:45, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))






# bias beta
dim(sim_data_normal01_simpred3$beta_true)
dim(res_normal01_simpred_hhssi3$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_simpred_hhssi3$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_simpred_nhhssi3$theta_update[, 2:45, ], c(2,3), mean)
# Bias for each outcome
beta_bias_hhssi <- (beta_est_hhssi - sim_data_normal01_simpred3$beta_true)
beta_bias_nhhssi <- (beta_est_nhhssi - sim_data_normal01_simpred3$beta_true)

dim(beta_bias_hhssi)
dim(beta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:44,
    beta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "modifiers",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(beta_bias_hhssi,
                   beta_bias_nhhssi)))
  lines(1:44,
        beta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  beta_bias_hhssi[, k],
  beta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(beta_bias_hhssi, beta_bias_nhhssi)),
  ylim = range(c(beta_bias_hhssi, beta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))





mean(beta_bias_hhssi>beta_bias_nhhssi)



beta_var_hhssi <- (res_normal01_simpred_hhssi3$lambdasq_beta_update *
                     res_normal01_simpred_hhssi3$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_simpred_nhhssi3$lambdasq_beta_update *
                      res_normal01_simpred_nhhssi3$tausq_beta_update)

mean(beta_var_hhssi>beta_var_nhhssi)
mean(beta_var_hhssi)
mean(beta_var_nhhssi)
mean(beta_bias_hhssi)
mean(beta_bias_nhhssi)

mean(beta_var_hhssi[1:2] > beta_var_nhhssi[1:2])
mean(beta_var_hhssi[1:2])
mean(beta_var_nhhssi[1:2])
mean(beta_bias_hhssi[1:2])
mean(beta_bias_nhhssi[1:2])

mean(beta_var_hhssi[-(1:2)] > beta_var_nhhssi[-(1:2)])
mean(beta_var_hhssi[-(1:2)])
mean(beta_var_nhhssi[-(1:2)])
mean(beta_bias_hhssi[-(1:2)])
mean(beta_bias_nhhssi[-(1:2)])



plot(beta_var_hhssi, beta_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[1:2], beta_var_nhhssi[1:2],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[-(1:2)], beta_var_nhhssi[-(1:2)],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")












# active gamma
range(res_normal01_simpred_hhssi3$theta_update[, 46:48, 1])
range(res_normal01_simpred_nhhssi3$theta_update[, 46:48, 1])


plot(res_normal01_simpred_hhssi3$theta_update[, 46:48, 1],
     res_normal01_simpred_nhhssi3$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi3$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi3$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, 46:48, 2],
       res_normal01_simpred_nhhssi3$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, 46:48, 3],
       res_normal01_simpred_nhhssi3$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, 46:48, 4],
       res_normal01_simpred_nhhssi3$theta_update[, 46:48, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_simpred_hhssi3$theta_update[, 46:48, 2])
range(res_normal01_simpred_nhhssi3$theta_update[, 46:48, 2])

range(res_normal01_simpred_hhssi3$theta_update[, 46:48, 3])
range(res_normal01_simpred_nhhssi3$theta_update[, 46:48, 3])

range(res_normal01_simpred_hhssi3$theta_update[, 46:48, 4])
range(res_normal01_simpred_nhhssi3$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, 46:48, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive gamma
range(res_normal01_simpred_hhssi3$theta_update[, 49:62, 1])
range(res_normal01_simpred_nhhssi3$theta_update[, 49:62, 1])


plot(res_normal01_simpred_hhssi3$theta_update[, 49:62, 1],
     res_normal01_simpred_nhhssi3$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi3$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi3$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, 49:62, 2],
       res_normal01_simpred_nhhssi3$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, 49:62, 3],
       res_normal01_simpred_nhhssi3$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, 49:62, 4],
       res_normal01_simpred_nhhssi3$theta_update[, 49:62, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi3$theta_update[, 49:62, 2])
range(res_normal01_simpred_nhhssi3$theta_update[, 49:62, 2])

range(res_normal01_simpred_hhssi3$theta_update[, 49:62, 3])
range(res_normal01_simpred_nhhssi3$theta_update[, 49:62, 3])

range(res_normal01_simpred_hhssi3$theta_update[, 49:62, 4])
range(res_normal01_simpred_nhhssi3$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, 49:62, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))






# bias gamma
dim(sim_data_normal01_simpred3$gamma_true)
dim(res_normal01_simpred_hhssi3$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_simpred_hhssi3$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_simpred_nhhssi3$theta_update[, 46:62, ], c(2,3), mean)
# Bias for each outcome
gamma_bias_hhssi <- (gamma_est_hhssi - sim_data_normal01_simpred3$gamma_true)
gamma_bias_nhhssi <- (gamma_est_nhhssi - sim_data_normal01_simpred3$gamma_true)

dim(gamma_bias_hhssi)
dim(gamma_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:17,
    gamma_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "exposures",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(gamma_bias_hhssi,
                   gamma_bias_nhhssi)))
  lines(1:17,
        gamma_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  gamma_bias_hhssi[, k],
  gamma_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)),
  ylim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))



gamma_var_hhssi <- (res_normal01_simpred_hhssi3$lambdasq_gamma_update *
                      res_normal01_simpred_hhssi3$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_simpred_nhhssi3$lambdasq_gamma_update *
                       res_normal01_simpred_nhhssi3$tausq_gamma_update)

mean(gamma_var_hhssi>gamma_var_nhhssi)
mean(gamma_var_hhssi)
mean(gamma_var_nhhssi)
mean(gamma_bias_hhssi)
mean(gamma_bias_nhhssi)

mean(gamma_var_hhssi[1:2] > gamma_var_nhhssi[1:2])
mean(gamma_var_hhssi[1:2])
mean(gamma_var_nhhssi[1:2])
mean(gamma_bias_hhssi[1:2])
mean(gamma_bias_nhhssi[1:2])

mean(gamma_var_hhssi[-(1:2)] > gamma_var_nhhssi[-(1:2)])
mean(gamma_var_hhssi[-(1:2)])
mean(gamma_var_nhhssi[-(1:2)])
mean(gamma_bias_hhssi[-(1:2)])
mean(gamma_bias_nhhssi[-(1:2)])



plot(gamma_var_hhssi, gamma_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[1:3], gamma_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[-(1:3)], gamma_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")










# active delta
range(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 1],
     res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 2],
       res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 3],
       res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 4],
       res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_simpred_hhssi3$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_simpred_nhhssi3$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, c(63:65,80:82), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive delta

range(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 1],
     res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 2],
       res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 3],
       res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 4],
       res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_simpred_hhssi3$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_simpred_nhhssi3$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, c(66:79,83:810), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# bias delta
dim(sim_data_normal01_simpred3$delta_true)
dim(res_normal01_simpred_hhssi3$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_simpred_hhssi3$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_simpred_nhhssi3$theta_update[, 63:810, ], c(2,3), mean)
# Bias for each outcome
delta_bias_hhssi <- (delta_est_hhssi - sim_data_normal01_simpred3$delta_true)
delta_bias_nhhssi <- (delta_est_nhhssi - sim_data_normal01_simpred3$delta_true)

dim(delta_bias_hhssi)
dim(delta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:748,
    delta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "interactions",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(delta_bias_hhssi,
                   delta_bias_nhhssi)))
  lines(1:748,
        delta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  delta_bias_hhssi[, k],
  delta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(delta_bias_hhssi, delta_bias_nhhssi)),
  ylim = range(c(delta_bias_hhssi, delta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






# active varphi
range(res_normal01_simpred_hhssi3$theta_update[, 811:812, 1])
range(res_normal01_simpred_nhhssi3$theta_update[, 811:812, 1])

plot(res_normal01_simpred_hhssi3$theta_update[, 811:812, 1],
     res_normal01_simpred_nhhssi3$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi3$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi3$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, 811:812, 2],
       res_normal01_simpred_nhhssi3$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, 811:812, 3],
       res_normal01_simpred_nhhssi3$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, 811:812, 4],
       res_normal01_simpred_nhhssi3$theta_update[, 811:812, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_simpred_hhssi3$theta_update[, 811:812, 2])
range(res_normal01_simpred_nhhssi3$theta_update[, 811:812, 2])

range(res_normal01_simpred_hhssi3$theta_update[, 811:812, 3])
range(res_normal01_simpred_nhhssi3$theta_update[, 811:812, 3])

range(res_normal01_simpred_hhssi3$theta_update[, 811:812, 4])
range(res_normal01_simpred_nhhssi3$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, 811:812, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive varphi
range(res_normal01_simpred_hhssi3$theta_update[, 813:816, 1])
range(res_normal01_simpred_nhhssi3$theta_update[, 813:816, 1])

plot(res_normal01_simpred_hhssi3$theta_update[, 813:816, 1],
     res_normal01_simpred_nhhssi3$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi3$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_simpred_hhssi3$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi3$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_simpred_hhssi3$theta_update[, 813:816, 2],
       res_normal01_simpred_nhhssi3$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$theta_update[, 813:816, 3],
       res_normal01_simpred_nhhssi3$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$theta_update[, 813:816, 4],
       res_normal01_simpred_nhhssi3$theta_update[, 813:816, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_simpred_hhssi3$theta_update[, 813:816, 2])
range(res_normal01_simpred_nhhssi3$theta_update[, 813:816, 2])

range(res_normal01_simpred_hhssi3$theta_update[, 813:816, 3])
range(res_normal01_simpred_nhhssi3$theta_update[, 813:816, 3])

range(res_normal01_simpred_hhssi3$theta_update[, 813:816, 4])
range(res_normal01_simpred_nhhssi3$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi3$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi3$theta_update[, 813:816, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))





# bias varphi
dim(sim_data_normal01_simpred3$varphi_true)
dim(res_normal01_simpred_hhssi3$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_simpred_hhssi3$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_simpred_nhhssi3$theta_update[, 811:816, ], c(2,3), mean)
# Bias for each outcome
varphi_bias_hhssi <- (varphi_est_hhssi - sim_data_normal01_simpred3$varphi_true)
varphi_bias_nhhssi <- (varphi_est_nhhssi - sim_data_normal01_simpred3$varphi_true)

dim(varphi_bias_hhssi)
dim(varphi_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:6,
    varphi_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "covariates",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(varphi_bias_hhssi,
                   varphi_bias_nhhssi)))
  lines(1:6,
        varphi_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  varphi_bias_hhssi[, k],
  varphi_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)),
  ylim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






varphi_var_hhssi <- apply(res_normal01_simpred_hhssi3$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_simpred_nhhssi3$theta_update[, 811:816, ], c(2,3), var)

mean(varphi_var_hhssi>varphi_var_nhhssi)
mean(varphi_var_hhssi)
mean(varphi_var_nhhssi)
mean(varphi_bias_hhssi)
mean(varphi_bias_nhhssi)

mean(varphi_var_hhssi[1:2] > varphi_var_nhhssi[1:2])
mean(varphi_var_hhssi[1:2])
mean(varphi_var_nhhssi[1:2])
mean(varphi_bias_hhssi[1:2])
mean(varphi_bias_nhhssi[1:2])

mean(varphi_var_hhssi[-(1:2)] > varphi_var_nhhssi[-(1:2)])
mean(varphi_var_hhssi[-(1:2)])
mean(varphi_var_nhhssi[-(1:2)])
mean(varphi_bias_hhssi[-(1:2)])
mean(varphi_bias_nhhssi[-(1:2)])



plot(varphi_var_hhssi, varphi_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[1:3], varphi_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[-(1:3)], varphi_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")











dim(sim_data_normal01_simpred3$gamma_true)
dim(res_normal01_simpred_hhssi3$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_simpred_hhssi3$lambdasq_delta_update)
dim(res_normal01_simpred_nhhssi3$lambdasq_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi3$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi3$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi3$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi3$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_simpred_hhssi3$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_simpred_nhhssi3$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_simpred_hhssi3$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_simpred_nhhssi3$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi3$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_simpred_nhhssi3$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_simpred_hhssi3$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_simpred_nhhssi3$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_simpred_hhssi3$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_simpred_nhhssi3$lambdasq_delta_update, 2, mean)
plot(lambdasq_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


# active deltas
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[c(1:3, 18:20)])
hist(lambdasq_delta_nhhssi[c(1:3, 18:20)])
range(lambdasq_delta_hhssi[c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[c(1:3, 18:20)])) # active deltas
hist(log(lambdasq_delta_nhhssi[c(1:3, 18:20)])) # active deltas

# inactive deltas
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[-c(1:3, 18:20)]) # inactive deltas
hist(lambdasq_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(lambdasq_delta_nhhssi[-c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)])) # inactive deltas
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)])) # inactive deltas


# active delta
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive delta
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")





# Note: the active lambdasq_delta values are lower for hhssi model than nhhssi, and range of values is (0.1485711, 2.5748392)
# the active lambdasq_delta values are higher for nhhssi model, and range of values is (58819.5, 403810.4)

# Note: whereas, the inactive lambdasq_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (1.246500e-01, 1.844546e+10)
# the inactive lambdasq_delta values are lower for nhhssi model, and range of values is (2.776089 234.788481)

# Overall, we can say that lambdasq_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(lambdasq_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(lambdasq_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
lambdasq_delta_hhssi_new <- lambdasq_delta_hhssi[-remove_hhssi]
lambdasq_delta_nhhssi_new <- lambdasq_delta_nhhssi[-remove_nhhssi]
# Plot
plot(lambdasq_delta_hhssi_new, lambdasq_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
lambdasq_delta_hhssi[remove_hhssi]
lambdasq_delta_nhhssi[remove_nhhssi]




hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", lambda[delta]^2)),
     xlab = expression(lambda[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]), log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]))))
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

mean(log(lambdasq_delta_hhssi))
mean(log(lambdasq_delta_nhhssi))

mean(lambdasq_delta_hhssi)
mean(lambdasq_delta_nhhssi)

# 1:2 active betas
plot(apply(res_normal01_simpred_hhssi3$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi3$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_simpred_hhssi3$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_simpred_nhhssi3$lambdasq_beta_update, 2, mean)
plot(lambdasq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[beta]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(lambdasq_beta_hhssi, lambdasq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)),
     ylim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_beta looks similar for hhssi and nhhssi models





# 1:3 active gammas
plot(apply(res_normal01_simpred_hhssi3$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi3$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_simpred_hhssi3$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi3$lambdasq_gamma_update, 2, mean)
plot(lambdasq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[gamma]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)),
     ylim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_gamma looks similar for hhssi and nhhssi models




names(res_normal01_simpred_hhssi3)
names(res_normal01_simpred_nhhssi3)





# tausq beta
mean(res_normal01_simpred_hhssi3$tausq_beta_update)
mean(res_normal01_simpred_nhhssi3$tausq_beta_update)
abs(mean(res_normal01_simpred_hhssi3$tausq_beta_update) -
      mean(res_normal01_simpred_nhhssi3$tausq_beta_update))

range(res_normal01_simpred_hhssi3$tausq_beta_update)
range(res_normal01_simpred_nhhssi3$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_simpred_hhssi3$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_simpred_nhhssi3$tausq_beta_update
plot(tausq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[beta]^2),
     lwd = 1,
     col = "red")
lines(tausq_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_beta_hhssi, tausq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)),
     ylim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
abline(a=0, b=1, col = "red")




hist(tausq_beta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
hist(tausq_beta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))


d1 <- density(tausq_beta_hhssi)
d2 <- density(tausq_beta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)




# tausq_gamma
mean(res_normal01_simpred_hhssi3$tausq_gamma_update)
mean(res_normal01_simpred_nhhssi3$tausq_gamma_update)
abs(mean(res_normal01_simpred_hhssi3$tausq_gamma_update) -
      mean(res_normal01_simpred_nhhssi3$tausq_gamma_update))

range(res_normal01_simpred_hhssi3$tausq_gamma_update)
range(res_normal01_simpred_nhhssi3$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_simpred_hhssi3$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_simpred_nhhssi3$tausq_gamma_update
plot(tausq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[gamma]^2),
     lwd = 1,
     col = "red")
lines(tausq_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_gamma_hhssi, tausq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)),
     ylim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")






# tausq_delta
mean(res_normal01_simpred_hhssi3$tausq_delta_update)
mean(res_normal01_simpred_nhhssi3$tausq_delta_update)
abs(mean(res_normal01_simpred_hhssi3$tausq_delta_update) -
      mean(res_normal01_simpred_nhhssi3$tausq_delta_update))

range(res_normal01_simpred_hhssi3$tausq_delta_update)
range(res_normal01_simpred_nhhssi3$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_simpred_hhssi3$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi3$tausq_delta_update
plot(tausq_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[delta]^2),
     lwd = 1,
     col = "red",
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
lines(tausq_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

hist(tausq_delta_hhssi, breaks = 50)
hist(tausq_delta_nhhssi, breaks = 50)

plot(tausq_delta_hhssi, tausq_delta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
abline(a=0, b=1, col = "red")



hist(tausq_delta_hhssi)
hist(tausq_delta_nhhssi)

hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

# log x axis
hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     log = "x")
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))



d1 <- density(tausq_delta_hhssi)
d2 <- density(tausq_delta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)


# Note: tausq_beta and tausq_gamma are same for the hhssi and nhhssi models.
# However, tausq_delta from hhssi and nhhssi models are far different.
# tausq_delta from hhssi are higher than that of nhhssi model.






# xi beta
mean(res_normal01_simpred_hhssi3$xi_beta_update)
mean(res_normal01_simpred_nhhssi3$xi_beta_update)
abs(mean(res_normal01_simpred_hhssi3$xi_beta_update) -
      mean(res_normal01_simpred_nhhssi3$xi_beta_update))

range(res_normal01_simpred_hhssi3$xi_beta_update)
range(res_normal01_simpred_nhhssi3$xi_beta_update)

xi_beta_hhssi <- res_normal01_simpred_hhssi3$xi_beta_update
xi_beta_nhhssi <- res_normal01_simpred_nhhssi3$xi_beta_update
plot(xi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[beta]),
     lwd = 1,
     col = "red")
lines(xi_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_beta_hhssi, xi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_beta_hhssi, xi_beta_nhhssi)),
     ylim = range(c(xi_beta_hhssi, xi_beta_nhhssi)))
abline(a=0, b=1, col = "red")





# xi_gamma
mean(res_normal01_simpred_hhssi3$xi_gamma_update)
mean(res_normal01_simpred_nhhssi3$xi_gamma_update)
abs(mean(res_normal01_simpred_hhssi3$xi_gamma_update) -
      mean(res_normal01_simpred_nhhssi3$xi_gamma_update))

range(res_normal01_simpred_hhssi3$xi_gamma_update)
range(res_normal01_simpred_nhhssi3$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_simpred_hhssi3$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_simpred_nhhssi3$xi_gamma_update
plot(xi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[gamma]),
     lwd = 1,
     col = "red")
lines(xi_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_gamma_hhssi, xi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)),
     ylim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")




# xi_delta
mean(res_normal01_simpred_hhssi3$xi_delta_update)
mean(res_normal01_simpred_nhhssi3$xi_delta_update)
abs(mean(res_normal01_simpred_hhssi3$xi_delta_update) -
      mean(res_normal01_simpred_nhhssi3$xi_delta_update))

range(res_normal01_simpred_hhssi3$xi_delta_update)
range(res_normal01_simpred_nhhssi3$xi_delta_update)

xi_delta_hhssi <- res_normal01_simpred_hhssi3$xi_delta_update
xi_delta_nhhssi <- res_normal01_simpred_nhhssi3$xi_delta_update
plot(xi_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[delta]),
     lwd = 1,
     col = "red",
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
lines(xi_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_delta_hhssi, xi_delta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_delta_hhssi, xi_delta_nhhssi)),
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
abline(a=0, b=1, col = "red")

hist(xi_delta_hhssi)
hist(xi_delta_nhhssi)



# Note: xi_beta and xi_gamma are same for both the hhssi and nhhssi models.
# However, xi_delta from hhssi and nhhssi models are far different.
# xi_delta from nhhssi are much higher than that of hhssi model.

##########################################################

# psi_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_simpred_nhhssi3$psi_delta_update)
dim(res_normal01_simpred_hhssi3$psi_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi3$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi3$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi3$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi3$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_simpred_hhssi3$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_simpred_nhhssi3$psi_delta_update, 2, mean)

# active
plot(psi_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red",
     ylim=range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
lines(psi_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive
plot(psi_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")



# active deltas
plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])
hist(psi_delta_hhssi[c(1:3, 18:20)])

range(psi_delta_hhssi[c(1:3, 18:20)])

# inactive deltas
plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])
hist(psi_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(psi_delta_nhhssi[-c(1:3, 18:20)])


# Note: the active psi_delta values are lower for hhssi model than nhhssi, and range of values is (0.1485711, 2.5748392)
# the active psi_delta values are higher for nhhssi model, and range of values is (58819.5, 403810.4)

# Note: whereas, the inactive psi_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (1.246500e-01, 1.844546e+10)
# the inactive psi_delta values are lower for nhhssi model, and range of values is (2.776089 234.788481)

# Overall, we can say that psi_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(psi_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(psi_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
psi_delta_hhssi_new <- psi_delta_hhssi[-remove_hhssi]
psi_delta_nhhssi_new <- psi_delta_nhhssi[-remove_nhhssi]
# Plot
plot(psi_delta_hhssi_new, psi_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
psi_delta_hhssi[remove_hhssi]
psi_delta_nhhssi[remove_nhhssi]




hist(log(psi_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", psi[delta]^2)),
     xlab = expression(psi[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(psi_delta_hhssi), log(psi_delta_nhhssi))))
hist(log(psi_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))




# 1:2 active betas
plot(apply(res_normal01_simpred_hhssi3$psi_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi3$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_simpred_hhssi3$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_simpred_nhhssi3$psi_beta_update, 2, mean)
plot(psi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[beta])),
     lwd = 1,
     col = "red")
lines(psi_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_beta_hhssi, psi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(psi_beta_hhssi, psi_beta_nhhssi)),
     ylim = range(c(psi_beta_hhssi, psi_beta_nhhssi)))
abline(a=0, b=1, col = "red")


# Note: psi_beta looks similar for hhssi and nhhssi models


# 1:3 active gammas
plot(apply(res_normal01_simpred_hhssi3$psi_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi3$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_simpred_hhssi3$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi3$psi_gamma_update, 2, mean)
plot(psi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[gamma])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
lines(psi_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_gamma_hhssi, psi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)),
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: psi_gamma looks similar for hhssi and nhhssi models


names(res_normal01_simpred_hhssi3)
names(res_normal01_simpred_nhhssi3)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_simpred_hhssi3$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_simpred_nhhssi3$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
    Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_simpred_hhssi3$Sigma_update[, , 1], type = "l")
matplot(res_normal01_simpred_nhhssi3$Sigma_update[, , 1], type = "l")

matplot(res_normal01_simpred_hhssi3$Sigma_update[, , 2], type = "l")
matplot(res_normal01_simpred_nhhssi3$Sigma_update[, , 2], type = "l")

matplot(res_normal01_simpred_hhssi3$Sigma_update[, , 3], type = "l")
matplot(res_normal01_simpred_nhhssi3$Sigma_update[, , 3], type = "l")

matplot(res_normal01_simpred_hhssi3$Sigma_update[, , 4], type = "l")
matplot(res_normal01_simpred_nhhssi3$Sigma_update[, , 4], type = "l")



plot(res_normal01_simpred_hhssi3$Sigma_update[, , 1],
     res_normal01_simpred_nhhssi3$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi3$Sigma_update,
                    res_normal01_simpred_hhssi3$Sigma_update)),
     ylim = range(c(res_normal01_simpred_hhssi3$Sigma_update,
                    res_normal01_simpred_hhssi3$Sigma_update)))
# k = 2
points(res_normal01_simpred_hhssi3$Sigma_update[, , 2],
       res_normal01_simpred_nhhssi3$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi3$Sigma_update[, , 3],
       res_normal01_simpred_nhhssi3$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi3$Sigma_update[, , 4],
       res_normal01_simpred_nhhssi3$Sigma_update[, , 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



# Note: Sigmasq are same for both hhssi and nhhssi models





###########

lambdasq_beta_hhssi
lambdasq_gamma_hhssi
lambdasq_delta_hhssi
tausq_delta_hhssi


dim(res_normal01_simpred_hhssi3$lambdasq_beta_update)
dim(res_normal01_simpred_hhssi3$lambdasq_gamma_update)
dim(res_normal01_simpred_hhssi3$lambdasq_delta_update)
length(res_normal01_simpred_hhssi3$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_simpred_hhssi3$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_simpred_hhssi3$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_simpred_hhssi3$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_simpred_hhssi3$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_simpred_nhhssi3$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi3$tausq_delta_update


n_iter <- dim(beta_lam_hhssi)[1]

# Construct beta x gamma component
beta_gamma_lam_hhssi <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_hhssi[t, col_idx] <- beta_lam_hhssi[t, j] * gamma_lam_hhssi[t, m]
      col_idx <- col_idx + 1
    }
  }
}
# prior variance
delta_var_hhssi <- delta_lam_hhssi * beta_gamma_lam_hhssi * tausq_delta_hhssi
# prior variance
delta_var_nhhssi <- delta_lam_nhhssi * tausq_delta_nhhssi



beta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi3$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi3$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi3$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_simpred_hhssi3$tausq_delta_update,
     type = "l",
     xlab = "MCMC iteration",
     ylab = expression(tau[delta]^2),
     main = "HHS-SI: global delta shrinkage"
)




hist(log10(beta_gamma_lam_hhssi),
     breaks = 50,
     main = expression("HHS-SI: " * log[10](lambda[beta]^2 * lambda[gamma]^2)),
     xlab = expression(log[10](lambda[beta]^2 * lambda[gamma]^2)))
abline(v = 0, lty = 2, lwd = 2)





delta_var_mean_hhssi <- apply(delta_var_hhssi, 2, mean)
beta_gamma_lam_mean_hhssi <- apply(beta_gamma_lam_hhssi, 2, mean)
# plot(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi,
#   log = "xy",
#   pch = 19,
#   xlab = "Delta prior var of hhssi from hhssi",
#   ylab = "Delta prior var of nhhssi from hhssi",
#   xlim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)),
#   ylim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)))
# abline(0, 1, lty = 2, lwd = 2, col = "red")


delta_var_mean_nhhssi <- apply(delta_var_nhhssi, 2, mean)
plot(delta_var_mean_hhssi, delta_var_mean_nhhssi,
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)),
     ylim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

which(delta_var_mean_hhssi[-c(1:3, 18:20)]>delta_var_mean_nhhssi[-c(1:3, 18:20)])

plot(delta_var_mean_nhhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")






delta_var_mean_hhssi[c(1:3, 18:20)] > delta_var_mean_nhhssi[c(1:3, 18:20)]
delta_var_mean_hhssi[-c(1:3, 18:20)] > delta_var_mean_nhhssi[-c(1:3, 18:20)]


ratio_hhssi_nhhssi <- delta_var_hhssi / delta_var_nhhssi
ratio_mean <- apply(ratio_hhssi_nhhssi, 2, mean)
hist(log10(ratio_mean), breaks = 50,
     main = "HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[delta,HHS-SI] / V[delta,NHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


active_delta_var_mean_hhssi <- delta_var_mean_hhssi[c(1:3, 18:20)]
active_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[c(1:3, 18:20)]
inactive_delta_var_mean_hhssi <- delta_var_mean_hhssi[-c(1:3, 18:20)]
inactive_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[-c(1:3, 18:20)]
active_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,c(1:3, 18:20)] /
                                          delta_var_nhhssi[,c(1:3, 18:20)])
inactive_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,-c(1:3, 18:20)] /
                                            delta_var_nhhssi[,-c(1:3, 18:20)])
active_delta_var_ratio_mean_hhssi_nhhssi <- apply(active_delta_var_ratio_hhssi_nhhssi, 2, mean)
inactive_delta_var_ratio_mean_hhssi_nhhssi <- apply(inactive_delta_var_ratio_hhssi_nhhssi, 2, mean)

hist(log10(active_delta_var_ratio_mean_hhssi_nhhssi), breaks = 5,
     main = "active HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[activedeltaHHS-SI] / V[activedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


hist(log10(inactive_delta_var_ratio_mean_hhssi_nhhssi), breaks = 50,
     main = "inactive HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[inactivedeltaHHS-SI] / V[inactivedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")









plot(beta_gamma_lam_mean_hhssi,
     pch = 19,
     log = "y",
     xlab = "Interaction index",
     ylab = expression(E[posterior](lambda[beta]^2 * lambda[gamma]^2)),
     main = "HHS-SI beta-gamma multiplier")
abline(h = 1, lty = 2, lwd = 2, col = "red")



















###############################################################################
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors real (Third sample)
###############################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
# Third sample
set.seed(23315644)
sim_data_normal01_realpred3 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                 O=O, x=X, z=Z, d=D,
                                                 Sigma=Sigma)

Y <- sim_data_normal01_realpred3$Y
x <- sim_data_normal01_realpred3$x
z <- sim_data_normal01_realpred3$z
u <- sim_data_normal01_realpred3$u
d <- sim_data_normal01_realpred3$d
W <- sim_data_normal01_realpred3$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_realpred_hhssi3.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi3.rda")

sim_data_normal01_realpred3$delta_true
dim(res_normal01_realpred_hhssi3$theta_update)






# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(res_normal01_realpred_hhssi3$theta_update[, 1, k],
          type = "l",
          lty = 1,
          lwd = 2,
          xlab = "MCMC iteration",
          ylab = expression(alpha),
          main = paste("Outcome", k)
  )

  matlines(res_normal01_realpred_nhhssi3$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_realpred_hhssi3$theta_update[, 1, ])
range(res_normal01_realpred_nhhssi3$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
  hhssi_alpha <- res_normal01_realpred_hhssi3$theta_update[, 1, k]
  nhhssi_alpha <- res_normal01_realpred_nhhssi3$theta_update[, 1, k]

  plot(hhssi_alpha, nhhssi_alpha,
       pch = 19,
       col = "blue",
       xlab = "HHS-SI",
       ylab = "NHHS-SI",
       main = paste("alpha: Outcome", k),
       xlim = range(c(hhssi_alpha, nhhssi_alpha)),
       ylim = range(c(hhssi_alpha, nhhssi_alpha)))
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
}
par(mfrow = c(1, 1))



# bias alpha
sim_data_normal01_realpred3$alpha_true
res_normal01_realpred_hhssi3$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_realpred_hhssi3$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_realpred_nhhssi3$theta_update[, 1, ], 2, mean)
# Bias for each outcome
alpha_bias_hhssi <- (alpha_est_hhssi - sim_data_normal01_realpred3$alpha_true)
alpha_bias_nhhssi <- (alpha_est_nhhssi - sim_data_normal01_realpred3$alpha_true)


plot(alpha_bias_hhssi,
     type = "b",
     xlab = "Outcome",
     ylab = "bias of alpha",
     lwd = 1,
     col = "red",
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
lines(alpha_bias_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(alpha_bias_hhssi, alpha_bias_nhhssi,
     lwd = 1,
     xlim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)),
     ylim = range(c(alpha_bias_hhssi, alpha_bias_nhhssi)))
abline(a=0, b=1, col = "red")



# active beta

plot(res_normal01_realpred_hhssi3$theta_update[, 2:3, 1],
     res_normal01_realpred_nhhssi3$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi3$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi3$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, 2:3, 2],
       res_normal01_realpred_nhhssi3$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, 2:3, 3],
       res_normal01_realpred_nhhssi3$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, 2:3, 4],
       res_normal01_realpred_nhhssi3$theta_update[, 2:3, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))




range(res_normal01_realpred_hhssi3$theta_update[, 2:3, 2])
range(res_normal01_realpred_nhhssi3$theta_update[, 2:3, 2])

range(res_normal01_realpred_hhssi3$theta_update[, 2:3, 3])
range(res_normal01_realpred_nhhssi3$theta_update[, 2:3, 3])

range(res_normal01_realpred_hhssi3$theta_update[, 2:3, 4])
range(res_normal01_realpred_nhhssi3$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, 2:3, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))




# inactive beta
range(res_normal01_realpred_hhssi3$theta_update[, 4:45, 1])
range(res_normal01_realpred_nhhssi3$theta_update[, 4:45, 1])

plot(res_normal01_realpred_hhssi3$theta_update[, 4:45, 1],
     res_normal01_realpred_nhhssi3$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi3$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi3$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, 4:45, 2],
       res_normal01_realpred_nhhssi3$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, 4:45, 3],
       res_normal01_realpred_nhhssi3$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, 4:45, 4],
       res_normal01_realpred_nhhssi3$theta_update[, 4:45, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



range(res_normal01_realpred_hhssi3$theta_update[, 4:45, 2])
range(res_normal01_realpred_nhhssi3$theta_update[, 4:45, 2])

range(res_normal01_realpred_hhssi3$theta_update[, 4:45, 3])
range(res_normal01_realpred_nhhssi3$theta_update[, 4:45, 3])

range(res_normal01_realpred_hhssi3$theta_update[, 4:45, 4])
range(res_normal01_realpred_nhhssi3$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, 4:45, k],
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
par(mfrow = c(1, 1))






# bias beta
dim(sim_data_normal01_realpred3$beta_true)
dim(res_normal01_realpred_hhssi3$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_realpred_hhssi3$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_realpred_nhhssi3$theta_update[, 2:45, ], c(2,3), mean)
# Bias for each outcome
beta_bias_hhssi <- (beta_est_hhssi - sim_data_normal01_realpred3$beta_true)
beta_bias_nhhssi <- (beta_est_nhhssi - sim_data_normal01_realpred3$beta_true)

dim(beta_bias_hhssi)
dim(beta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:44,
    beta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "modifiers",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(beta_bias_hhssi,
                   beta_bias_nhhssi)))
  lines(1:44,
        beta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  beta_bias_hhssi[, k],
  beta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(beta_bias_hhssi, beta_bias_nhhssi)),
  ylim = range(c(beta_bias_hhssi, beta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))





mean(beta_bias_hhssi>beta_bias_nhhssi)



beta_var_hhssi <- (res_normal01_realpred_hhssi3$lambdasq_beta_update *
                     res_normal01_realpred_hhssi3$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_realpred_nhhssi3$lambdasq_beta_update *
                      res_normal01_realpred_nhhssi3$tausq_beta_update)

mean(beta_var_hhssi>beta_var_nhhssi)
mean(beta_var_hhssi)
mean(beta_var_nhhssi)
mean(beta_bias_hhssi)
mean(beta_bias_nhhssi)

mean(beta_var_hhssi[1:2] > beta_var_nhhssi[1:2])
mean(beta_var_hhssi[1:2])
mean(beta_var_nhhssi[1:2])
mean(beta_bias_hhssi[1:2])
mean(beta_bias_nhhssi[1:2])

mean(beta_var_hhssi[-(1:2)] > beta_var_nhhssi[-(1:2)])
mean(beta_var_hhssi[-(1:2)])
mean(beta_var_nhhssi[-(1:2)])
mean(beta_bias_hhssi[-(1:2)])
mean(beta_bias_nhhssi[-(1:2)])



plot(beta_var_hhssi, beta_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[1:2], beta_var_nhhssi[1:2],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(beta_var_hhssi[-(1:2)], beta_var_nhhssi[-(1:2)],
     pch = 19,
     xlab = "HHS-SI beta var",
     ylab = "NHHS-SI beta var",
     xlim = range(c(beta_var_hhssi, beta_var_nhhssi)),
     ylim = range(c(beta_var_hhssi, beta_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")












# active gamma
range(res_normal01_realpred_hhssi3$theta_update[, 46:48, 1])
range(res_normal01_realpred_nhhssi3$theta_update[, 46:48, 1])


plot(res_normal01_realpred_hhssi3$theta_update[, 46:48, 1],
     res_normal01_realpred_nhhssi3$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi3$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi3$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, 46:48, 2],
       res_normal01_realpred_nhhssi3$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, 46:48, 3],
       res_normal01_realpred_nhhssi3$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, 46:48, 4],
       res_normal01_realpred_nhhssi3$theta_update[, 46:48, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_realpred_hhssi3$theta_update[, 46:48, 2])
range(res_normal01_realpred_nhhssi3$theta_update[, 46:48, 2])

range(res_normal01_realpred_hhssi3$theta_update[, 46:48, 3])
range(res_normal01_realpred_nhhssi3$theta_update[, 46:48, 3])

range(res_normal01_realpred_hhssi3$theta_update[, 46:48, 4])
range(res_normal01_realpred_nhhssi3$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, 46:48, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive gamma
range(res_normal01_realpred_hhssi3$theta_update[, 49:62, 1])
range(res_normal01_realpred_nhhssi3$theta_update[, 49:62, 1])


plot(res_normal01_realpred_hhssi3$theta_update[, 49:62, 1],
     res_normal01_realpred_nhhssi3$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi3$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi3$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, 49:62, 2],
       res_normal01_realpred_nhhssi3$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, 49:62, 3],
       res_normal01_realpred_nhhssi3$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, 49:62, 4],
       res_normal01_realpred_nhhssi3$theta_update[, 49:62, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi3$theta_update[, 49:62, 2])
range(res_normal01_realpred_nhhssi3$theta_update[, 49:62, 2])

range(res_normal01_realpred_hhssi3$theta_update[, 49:62, 3])
range(res_normal01_realpred_nhhssi3$theta_update[, 49:62, 3])

range(res_normal01_realpred_hhssi3$theta_update[, 49:62, 4])
range(res_normal01_realpred_nhhssi3$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, 49:62, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))






# bias gamma
dim(sim_data_normal01_realpred3$gamma_true)
dim(res_normal01_realpred_hhssi3$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_realpred_hhssi3$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_realpred_nhhssi3$theta_update[, 46:62, ], c(2,3), mean)
# Bias for each outcome
gamma_bias_hhssi <- (gamma_est_hhssi - sim_data_normal01_realpred3$gamma_true)
gamma_bias_nhhssi <- (gamma_est_nhhssi - sim_data_normal01_realpred3$gamma_true)

dim(gamma_bias_hhssi)
dim(gamma_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:17,
    gamma_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "exposures",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(gamma_bias_hhssi,
                   gamma_bias_nhhssi)))
  lines(1:17,
        gamma_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  gamma_bias_hhssi[, k],
  gamma_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)),
  ylim = range(c(gamma_bias_hhssi, gamma_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))



gamma_var_hhssi <- (res_normal01_realpred_hhssi3$lambdasq_gamma_update *
                      res_normal01_realpred_hhssi3$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_realpred_nhhssi3$lambdasq_gamma_update *
                       res_normal01_realpred_nhhssi3$tausq_gamma_update)

mean(gamma_var_hhssi>gamma_var_nhhssi)
mean(gamma_var_hhssi)
mean(gamma_var_nhhssi)
mean(gamma_bias_hhssi)
mean(gamma_bias_nhhssi)

mean(gamma_var_hhssi[1:2] > gamma_var_nhhssi[1:2])
mean(gamma_var_hhssi[1:2])
mean(gamma_var_nhhssi[1:2])
mean(gamma_bias_hhssi[1:2])
mean(gamma_bias_nhhssi[1:2])

mean(gamma_var_hhssi[-(1:2)] > gamma_var_nhhssi[-(1:2)])
mean(gamma_var_hhssi[-(1:2)])
mean(gamma_var_nhhssi[-(1:2)])
mean(gamma_bias_hhssi[-(1:2)])
mean(gamma_bias_nhhssi[-(1:2)])



plot(gamma_var_hhssi, gamma_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[1:3], gamma_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(gamma_var_hhssi[-(1:3)], gamma_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI gamma var",
     ylab = "NHHS-SI gamma var",
     xlim = range(c(gamma_var_hhssi, gamma_var_nhhssi)),
     ylim = range(c(gamma_var_hhssi, gamma_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")










# active delta
range(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 1],
     res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 2],
       res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 3],
       res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 4],
       res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_realpred_hhssi3$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_realpred_nhhssi3$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, c(63:65,80:82), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive delta

range(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 1],
     res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 2],
       res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 3],
       res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 4],
       res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_realpred_hhssi3$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_realpred_nhhssi3$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, c(66:79,83:810), k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# bias delta
dim(sim_data_normal01_realpred3$delta_true)
dim(res_normal01_realpred_hhssi3$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_realpred_hhssi3$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_realpred_nhhssi3$theta_update[, 63:810, ], c(2,3), mean)
# Bias for each outcome
delta_bias_hhssi <- (delta_est_hhssi - sim_data_normal01_realpred3$delta_true)
delta_bias_nhhssi <- (delta_est_nhhssi - sim_data_normal01_realpred3$delta_true)

dim(delta_bias_hhssi)
dim(delta_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:748,
    delta_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "interactions",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(delta_bias_hhssi,
                   delta_bias_nhhssi)))
  lines(1:748,
        delta_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  delta_bias_hhssi[, k],
  delta_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(delta_bias_hhssi, delta_bias_nhhssi)),
  ylim = range(c(delta_bias_hhssi, delta_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






# active varphi
range(res_normal01_realpred_hhssi3$theta_update[, 811:812, 1])
range(res_normal01_realpred_nhhssi3$theta_update[, 811:812, 1])

plot(res_normal01_realpred_hhssi3$theta_update[, 811:812, 1],
     res_normal01_realpred_nhhssi3$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi3$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi3$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, 811:812, 2],
       res_normal01_realpred_nhhssi3$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, 811:812, 3],
       res_normal01_realpred_nhhssi3$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, 811:812, 4],
       res_normal01_realpred_nhhssi3$theta_update[, 811:812, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))


range(res_normal01_realpred_hhssi3$theta_update[, 811:812, 2])
range(res_normal01_realpred_nhhssi3$theta_update[, 811:812, 2])

range(res_normal01_realpred_hhssi3$theta_update[, 811:812, 3])
range(res_normal01_realpred_nhhssi3$theta_update[, 811:812, 3])

range(res_normal01_realpred_hhssi3$theta_update[, 811:812, 4])
range(res_normal01_realpred_nhhssi3$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, 811:812, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))




# inactive varphi
range(res_normal01_realpred_hhssi3$theta_update[, 813:816, 1])
range(res_normal01_realpred_nhhssi3$theta_update[, 813:816, 1])

plot(res_normal01_realpred_hhssi3$theta_update[, 813:816, 1],
     res_normal01_realpred_nhhssi3$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi3$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_realpred_hhssi3$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi3$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_realpred_hhssi3$theta_update[, 813:816, 2],
       res_normal01_realpred_nhhssi3$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$theta_update[, 813:816, 3],
       res_normal01_realpred_nhhssi3$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$theta_update[, 813:816, 4],
       res_normal01_realpred_nhhssi3$theta_update[, 813:816, 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))

range(res_normal01_realpred_hhssi3$theta_update[, 813:816, 2])
range(res_normal01_realpred_nhhssi3$theta_update[, 813:816, 2])

range(res_normal01_realpred_hhssi3$theta_update[, 813:816, 3])
range(res_normal01_realpred_nhhssi3$theta_update[, 813:816, 3])

range(res_normal01_realpred_hhssi3$theta_update[, 813:816, 4])
range(res_normal01_realpred_nhhssi3$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi3$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi3$theta_update[, 813:816, k],
    lty = 2,
    lwd = 2)
  legend(
    "topright",
    legend = c("HHS-SI", "NHHS-SI"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n")
}
par(mfrow = c(1, 1))





# bias varphi
dim(sim_data_normal01_realpred3$varphi_true)
dim(res_normal01_realpred_hhssi3$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_realpred_hhssi3$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_realpred_nhhssi3$theta_update[, 811:816, ], c(2,3), mean)
# Bias for each outcome
varphi_bias_hhssi <- (varphi_est_hhssi - sim_data_normal01_realpred3$varphi_true)
varphi_bias_nhhssi <- (varphi_est_nhhssi - sim_data_normal01_realpred3$varphi_true)

dim(varphi_bias_hhssi)
dim(varphi_bias_nhhssi)


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    1:6,
    varphi_bias_hhssi[, k],
    type = "b",
    pch = 19,
    col = "blue",
    xlab = "covariates",
    ylab = "Bias",
    main = paste("Outcome", k),
    ylim = range(c(varphi_bias_hhssi,
                   varphi_bias_nhhssi)))
  lines(1:6,
        varphi_bias_nhhssi[, k],
        type = "b",
        pch = 17,
        col = "red")

  abline(h = 0, lty = 2, lwd = 2)
  legend("topright",
         legend = c("HHS-SI", "NHHS-SI"),
         col = c("blue", "red"),
         pch = c(19, 17),
         lty = 1,
         bty = "n")
}


par(mfrow = c(2, 2))
for (k in 1:4) {plot(
  varphi_bias_hhssi[, k],
  varphi_bias_nhhssi[, k],
  pch = 19,
  xlab = "HHS-SI Bias",
  ylab = "NHHS-SI Bias",
  main = paste("Outcome", k),
  xlim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)),
  ylim = range(c(varphi_bias_hhssi, varphi_bias_nhhssi)))

  abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
  abline(h = 0, v = 0, lty = 3)
}
par(mfrow = c(1, 1))






varphi_var_hhssi <- apply(res_normal01_realpred_hhssi3$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_realpred_nhhssi3$theta_update[, 811:816, ], c(2,3), var)

mean(varphi_var_hhssi>varphi_var_nhhssi)
mean(varphi_var_hhssi)
mean(varphi_var_nhhssi)
mean(varphi_bias_hhssi)
mean(varphi_bias_nhhssi)

mean(varphi_var_hhssi[1:2] > varphi_var_nhhssi[1:2])
mean(varphi_var_hhssi[1:2])
mean(varphi_var_nhhssi[1:2])
mean(varphi_bias_hhssi[1:2])
mean(varphi_bias_nhhssi[1:2])

mean(varphi_var_hhssi[-(1:2)] > varphi_var_nhhssi[-(1:2)])
mean(varphi_var_hhssi[-(1:2)])
mean(varphi_var_nhhssi[-(1:2)])
mean(varphi_bias_hhssi[-(1:2)])
mean(varphi_bias_nhhssi[-(1:2)])



plot(varphi_var_hhssi, varphi_var_nhhssi,
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[1:3], varphi_var_nhhssi[1:3],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(varphi_var_hhssi[-(1:3)], varphi_var_nhhssi[-(1:3)],
     pch = 19,
     xlab = "HHS-SI varphi var",
     ylab = "NHHS-SI varphi var",
     xlim = range(c(varphi_var_hhssi, varphi_var_nhhssi)),
     ylim = range(c(varphi_var_hhssi, varphi_var_nhhssi)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")











dim(sim_data_normal01_realpred3$gamma_true)
dim(res_normal01_realpred_hhssi3$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_realpred_hhssi3$lambdasq_delta_update)
dim(res_normal01_realpred_nhhssi3$lambdasq_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi3$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi3$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi3$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi3$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_realpred_hhssi3$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_realpred_nhhssi3$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_realpred_hhssi3$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_realpred_nhhssi3$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi3$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_realpred_nhhssi3$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_realpred_hhssi3$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_realpred_nhhssi3$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_realpred_hhssi3$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_realpred_nhhssi3$lambdasq_delta_update, 2, mean)
plot(lambdasq_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[delta]^2)),
     lwd = 1,
     col = "red",
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)],
                    lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
lines(lambdasq_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


# active deltas
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[c(1:3, 18:20)])
hist(lambdasq_delta_nhhssi[c(1:3, 18:20)])
range(lambdasq_delta_hhssi[c(1:3, 18:20)])
range(lambdasq_delta_nhhssi[c(1:3, 18:20)])

hist(log(lambdasq_delta_hhssi[c(1:3, 18:20)])) # active deltas
hist(log(lambdasq_delta_nhhssi[c(1:3, 18:20)])) # active deltas

# inactive deltas
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])
hist(lambdasq_delta_hhssi[-c(1:3, 18:20)]) # inactive deltas
hist(lambdasq_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(lambdasq_delta_hhssi[-c(1:3, 18:20)])
range(lambdasq_delta_nhhssi[-c(1:3, 18:20)])


hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)])) # inactive deltas
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)])) # inactive deltas


# active delta
plot(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[c(1:3, 18:20)], lambdasq_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive delta
plot(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_hhssi[-c(1:3, 18:20)], lambdasq_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")





# Note: the active lambdasq_delta values are lower for hhssi model than nhhssi, and range of values is  (3.790582 1546.643649)
# the active lambdasq_delta values are higher for nhhssi model, and range of values is (154557.2 2087763.1)

# Note: whereas, the inactive lambdasq_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (2.723865e-01 1.864608e+07)
# the inactive lambdasq_delta values are lower for nhhssi model, and range of values is (1.771442 956.855391)

# Overall, we can say that lambdasq_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(lambdasq_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(lambdasq_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
lambdasq_delta_hhssi_new <- lambdasq_delta_hhssi[-remove_hhssi]
lambdasq_delta_nhhssi_new <- lambdasq_delta_nhhssi[-remove_nhhssi]
# Plot
plot(lambdasq_delta_hhssi_new, lambdasq_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
lambdasq_delta_hhssi[remove_hhssi]
lambdasq_delta_nhhssi[remove_nhhssi]




hist(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", lambda[delta]^2)),
     xlab = expression(lambda[delta]^2),
     ylab = "Frequency",
     xlim = range(c(log(lambdasq_delta_hhssi[-c(1:3, 18:20)]), log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]))))
hist(log(lambdasq_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

mean(log(lambdasq_delta_hhssi))
mean(log(lambdasq_delta_nhhssi))

mean(lambdasq_delta_hhssi)
mean(lambdasq_delta_nhhssi)

# 1:2 active betas
plot(apply(res_normal01_realpred_hhssi3$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi3$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_realpred_hhssi3$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_realpred_nhhssi3$lambdasq_beta_update, 2, mean)
plot(lambdasq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[beta]^2)),
     lwd = 1,
     col = "red",
     ylim=range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
lines(lambdasq_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(lambdasq_beta_hhssi, lambdasq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)),
     ylim = range(c(lambdasq_beta_hhssi, lambdasq_beta_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_beta looks larger for nhhssi than hhssi model





# 1:3 active gammas
plot(apply(res_normal01_realpred_hhssi3$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi3$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_realpred_hhssi3$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi3$lambdasq_gamma_update, 2, mean)
plot(lambdasq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(lambda[gamma]^2)),
     lwd = 1,
     col = "red")
lines(lambdasq_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)),
     ylim = range(c(lambdasq_gamma_hhssi, lambdasq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: lambdasq_gamma looks similar for hhssi and nhhssi models




names(res_normal01_realpred_hhssi3)
names(res_normal01_realpred_nhhssi3)





# tausq beta
mean(res_normal01_realpred_hhssi3$tausq_beta_update)
mean(res_normal01_realpred_nhhssi3$tausq_beta_update)
abs(mean(res_normal01_realpred_hhssi3$tausq_beta_update) -
      mean(res_normal01_realpred_nhhssi3$tausq_beta_update))

range(res_normal01_realpred_hhssi3$tausq_beta_update)
range(res_normal01_realpred_nhhssi3$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_realpred_hhssi3$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_realpred_nhhssi3$tausq_beta_update
plot(tausq_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[beta]^2),
     lwd = 1,
     col = "red")
lines(tausq_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_beta_hhssi, tausq_beta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)),
     ylim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
abline(a=0, b=1, col = "red")




hist(tausq_beta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_beta_hhssi, tausq_beta_nhhssi)))
hist(tausq_beta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))


d1 <- density(tausq_beta_hhssi)
d2 <- density(tausq_beta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[beta]^2)),
     xlab = expression(tau[beta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)




# tausq_gamma
mean(res_normal01_realpred_hhssi3$tausq_gamma_update)
mean(res_normal01_realpred_nhhssi3$tausq_gamma_update)
abs(mean(res_normal01_realpred_hhssi3$tausq_gamma_update) -
      mean(res_normal01_realpred_nhhssi3$tausq_gamma_update))

range(res_normal01_realpred_hhssi3$tausq_gamma_update)
range(res_normal01_realpred_nhhssi3$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_realpred_hhssi3$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_realpred_nhhssi3$tausq_gamma_update
plot(tausq_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[gamma]^2),
     lwd = 1,
     col = "red")
lines(tausq_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(tausq_gamma_hhssi, tausq_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)),
     ylim = range(c(tausq_gamma_hhssi, tausq_gamma_nhhssi)))
abline(a=0, b=1, col = "red")






# tausq_delta
mean(res_normal01_realpred_hhssi3$tausq_delta_update)
mean(res_normal01_realpred_nhhssi3$tausq_delta_update)
abs(mean(res_normal01_realpred_hhssi3$tausq_delta_update) -
      mean(res_normal01_realpred_nhhssi3$tausq_delta_update))

range(res_normal01_realpred_hhssi3$tausq_delta_update)
range(res_normal01_realpred_nhhssi3$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_realpred_hhssi3$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi3$tausq_delta_update
plot(tausq_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(tau[delta]^2),
     lwd = 1,
     col = "red",
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
lines(tausq_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(tausq_delta_hhssi, tausq_delta_nhhssi,
     lwd = 1,
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     ylim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
abline(a=0, b=1, col = "red")



hist(tausq_delta_hhssi)
hist(tausq_delta_nhhssi)

hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)))
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))

# log x axis
hist(tausq_delta_hhssi,
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylab = "Frequency",
     xlim = range(c(tausq_delta_hhssi, tausq_delta_nhhssi)),
     log = "x")
hist(tausq_delta_nhhssi,
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))



d1 <- density(tausq_delta_hhssi)
d2 <- density(tausq_delta_nhhssi)
plot(d1,
     col = "red",
     lwd = 2,
     main = expression(paste("Density of ", tau[delta]^2)),
     xlab = expression(tau[delta]^2),
     ylim = range(c(d1$y, d2$y)))
lines(d2,
      col = "blue",
      lwd = 2)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lwd = 2)


# Note: tausq_beta is larger for hhssi than nhhssi
# tausq_gamma is same for the hhssi and nhhssi models.
# tausq_delta from hhssi and nhhssi models differs.
# tausq_delta from hhssi are higher than that of nhhssi model.






# xi beta
mean(res_normal01_realpred_hhssi3$xi_beta_update)
mean(res_normal01_realpred_nhhssi3$xi_beta_update)
abs(mean(res_normal01_realpred_hhssi3$xi_beta_update) -
      mean(res_normal01_realpred_nhhssi3$xi_beta_update))

range(res_normal01_realpred_hhssi3$xi_beta_update)
range(res_normal01_realpred_nhhssi3$xi_beta_update)

xi_beta_hhssi <- res_normal01_realpred_hhssi3$xi_beta_update
xi_beta_nhhssi <- res_normal01_realpred_nhhssi3$xi_beta_update
plot(xi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[beta]),
     lwd = 1,
     col = "red",
     ylim=range(c(xi_beta_hhssi,xi_beta_nhhssi)))
lines(xi_beta_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_beta_hhssi, xi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_beta_hhssi, xi_beta_nhhssi)),
     ylim = range(c(xi_beta_hhssi, xi_beta_nhhssi)))
abline(a=0, b=1, col = "red")





# xi_gamma
mean(res_normal01_realpred_hhssi3$xi_gamma_update)
mean(res_normal01_realpred_nhhssi3$xi_gamma_update)
abs(mean(res_normal01_realpred_hhssi3$xi_gamma_update) -
      mean(res_normal01_realpred_nhhssi3$xi_gamma_update))

range(res_normal01_realpred_hhssi3$xi_gamma_update)
range(res_normal01_realpred_nhhssi3$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_realpred_hhssi3$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_realpred_nhhssi3$xi_gamma_update
plot(xi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[gamma]),
     lwd = 1,
     col = "red")
lines(xi_gamma_nhhssi,
      type = "p",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_gamma_hhssi, xi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)),
     ylim = range(c(xi_gamma_hhssi, xi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")




# xi_delta
mean(res_normal01_realpred_hhssi3$xi_delta_update)
mean(res_normal01_realpred_nhhssi3$xi_delta_update)
abs(mean(res_normal01_realpred_hhssi3$xi_delta_update) -
      mean(res_normal01_realpred_nhhssi3$xi_delta_update))

range(res_normal01_realpred_hhssi3$xi_delta_update)
range(res_normal01_realpred_nhhssi3$xi_delta_update)

xi_delta_hhssi <- res_normal01_realpred_hhssi3$xi_delta_update
xi_delta_nhhssi <- res_normal01_realpred_nhhssi3$xi_delta_update
plot(xi_delta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(xi[delta]),
     lwd = 1,
     col = "red",
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
lines(xi_delta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(xi_delta_hhssi, xi_delta_nhhssi,
     lwd = 1,
     xlim = range(c(xi_delta_hhssi, xi_delta_nhhssi)),
     ylim = range(c(xi_delta_hhssi, xi_delta_nhhssi)))
abline(a=0, b=1, col = "red")

hist(xi_delta_hhssi)
hist(xi_delta_nhhssi)



# Note: xi_beta for nhhssi is larger than hhssi
# xi_gamma similar for both the hhssi and nhhssi models.
# However, xi_delta from hhssi and nhhssi models are far different.
# xi_delta from nhhssi are much higher than that of hhssi model.

##########################################################

# psi_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_realpred_nhhssi3$psi_delta_update)
dim(res_normal01_realpred_hhssi3$psi_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi3$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi3$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi3$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi3$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_realpred_hhssi3$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_realpred_nhhssi3$psi_delta_update, 2, mean)

# active
plot(psi_delta_hhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")


# inactive
plot(psi_delta_hhssi[-c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[delta])),
     lwd = 1,
     col = "red")
lines(psi_delta_nhhssi[-c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)],
     lwd = 1,
     xlim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])))
abline(a=0, b=1, col = "red")



# active deltas
plot(psi_delta_hhssi[c(1:3, 18:20)], psi_delta_nhhssi[c(1:3, 18:20)])
hist(psi_delta_hhssi[c(1:3, 18:20)])
range(psi_delta_hhssi[c(1:3, 18:20)])
range(psi_delta_nhhssi[c(1:3, 18:20)])


# inactive deltas
plot(psi_delta_hhssi[-c(1:3, 18:20)], psi_delta_nhhssi[-c(1:3, 18:20)])
hist(psi_delta_nhhssi[-c(1:3, 18:20)]) # inactive deltas
range(psi_delta_hhssi[-c(1:3, 18:20)])
range(psi_delta_nhhssi[-c(1:3, 18:20)])


# Note: the active psi_delta values are larger for hhssi model than nhhssi, and range of values is (4.404875 61.890418)
# the active psi_delta values are lower for nhhssi model, and range of values is (4.088125 58.551846)

# Note: whereas, the inactive psi_delta values are higher (there are some lower too) for hhssi model than nhhssi, and range of values is (9.028653e+00 2.638961e+07)
# the inactive psi_delta values are lower for nhhssi model, and range of values is (1.640430e+01 3.505838e+08)

# Overall, we can say that psi_delta values behave differently for hhssi and nhhssi models

# Find the indices of the 6 highest values in each vector
remove_hhssi <- order(psi_delta_hhssi, decreasing = TRUE)[1:6]
remove_nhhssi <- order(psi_delta_nhhssi, decreasing = TRUE)[1:6]
# Remove those values
psi_delta_hhssi_new <- psi_delta_hhssi[-remove_hhssi]
psi_delta_nhhssi_new <- psi_delta_nhhssi[-remove_nhhssi]
# Plot
plot(psi_delta_hhssi_new, psi_delta_nhhssi_new,
     xlab = "HHS-SI",
     ylab = "NHHS-SI")
psi_delta_hhssi[remove_hhssi]
psi_delta_nhhssi[remove_nhhssi]




hist(log(psi_delta_hhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     main = expression(paste("Distribution of log of ", psi[delta])),
     xlab = expression(psi[delta]),
     ylab = "Frequency",
     xlim = range(c(log(psi_delta_hhssi), log(psi_delta_nhhssi))))
hist(log(psi_delta_nhhssi[-c(1:3, 18:20)]),
     breaks = 30,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     add = TRUE)
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)),
       border = c("red", "blue"))




# 1:2 active betas
plot(apply(res_normal01_realpred_hhssi3$psi_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi3$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_realpred_hhssi3$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_realpred_nhhssi3$psi_beta_update, 2, mean)
plot(psi_beta_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[beta])),
     lwd = 1,
     col = "red",
     ylim=range(c(psi_beta_hhssi,psi_beta_nhhssi)))
lines(psi_beta_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_beta_hhssi, psi_beta_nhhssi,
     lwd = 1,
     xlim = range(c(psi_beta_hhssi, psi_beta_nhhssi)),
     ylim = range(c(psi_beta_hhssi, psi_beta_nhhssi)))
abline(a=0, b=1, col = "red")


# Note: psi_beta looks similar for hhssi and nhhssi models


# 1:3 active gammas
plot(apply(res_normal01_realpred_hhssi3$psi_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi3$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_realpred_hhssi3$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi3$psi_gamma_update, 2, mean)
plot(psi_gamma_hhssi,
     type = "b",
     xlab = "Parameter",
     ylab = expression(mean(psi[gamma])),
     lwd = 1,
     col = "red",
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
lines(psi_gamma_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)

plot(psi_gamma_hhssi, psi_gamma_nhhssi,
     lwd = 1,
     xlim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)),
     ylim = range(c(psi_gamma_hhssi, psi_gamma_nhhssi)))
abline(a=0, b=1, col = "red")

# Note: psi_gamma looks similar for hhssi and nhhssi models


names(res_normal01_realpred_hhssi3)
names(res_normal01_realpred_nhhssi3)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_realpred_hhssi3$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_realpred_nhhssi3$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
    Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_realpred_hhssi3$Sigma_update[, , 1], type = "l")
matplot(res_normal01_realpred_nhhssi3$Sigma_update[, , 1], type = "l")

matplot(res_normal01_realpred_hhssi3$Sigma_update[, , 2], type = "l")
matplot(res_normal01_realpred_nhhssi3$Sigma_update[, , 2], type = "l")

matplot(res_normal01_realpred_hhssi3$Sigma_update[, , 3], type = "l")
matplot(res_normal01_realpred_nhhssi3$Sigma_update[, , 3], type = "l")

matplot(res_normal01_realpred_hhssi3$Sigma_update[, , 4], type = "l")
matplot(res_normal01_realpred_nhhssi3$Sigma_update[, , 4], type = "l")



plot(res_normal01_realpred_hhssi3$Sigma_update[, , 1],
     res_normal01_realpred_nhhssi3$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi3$Sigma_update,
                    res_normal01_realpred_hhssi3$Sigma_update)),
     ylim = range(c(res_normal01_realpred_hhssi3$Sigma_update,
                    res_normal01_realpred_hhssi3$Sigma_update)))
# k = 2
points(res_normal01_realpred_hhssi3$Sigma_update[, , 2],
       res_normal01_realpred_nhhssi3$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi3$Sigma_update[, , 3],
       res_normal01_realpred_nhhssi3$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi3$Sigma_update[, , 4],
       res_normal01_realpred_nhhssi3$Sigma_update[, , 4],
       pch = 19,
       col = "salmon")
# Reference line
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")
# Legend
legend("topleft",
       legend = c("k = 1", "k = 2", "k = 3", "k = 4", "y = x"),
       col = c("yellow", "lightblue", "green", "salmon", "red"),
       pch = c(19, 19, 19, 19, NA),
       lty = c(NA, NA, NA, NA, 2),
       lwd = c(NA, NA, NA, NA, 2))



# Note: Sigmasq are same for both hhssi and nhhssi models





###########

lambdasq_beta_hhssi
lambdasq_gamma_hhssi
lambdasq_delta_hhssi
tausq_delta_hhssi


dim(res_normal01_realpred_hhssi3$lambdasq_beta_update)
dim(res_normal01_realpred_hhssi3$lambdasq_gamma_update)
dim(res_normal01_realpred_hhssi3$lambdasq_delta_update)
length(res_normal01_realpred_hhssi3$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_realpred_hhssi3$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_realpred_hhssi3$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_realpred_hhssi3$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_realpred_hhssi3$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_realpred_nhhssi3$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi3$tausq_delta_update


n_iter <- dim(beta_lam_hhssi)[1]

# Construct beta x gamma component
beta_gamma_lam_hhssi <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_hhssi[t, col_idx] <- beta_lam_hhssi[t, j] * gamma_lam_hhssi[t, m]
      col_idx <- col_idx + 1
    }
  }
}
# prior variance
delta_var_hhssi <- delta_lam_hhssi * beta_gamma_lam_hhssi * tausq_delta_hhssi
# prior variance
delta_var_nhhssi <- delta_lam_nhhssi * tausq_delta_nhhssi



beta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi3$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi3$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi3$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_realpred_hhssi3$tausq_delta_update,
     type = "l",
     xlab = "MCMC iteration",
     ylab = expression(tau[delta]^2),
     main = "HHS-SI: global delta shrinkage"
)




hist(log10(beta_gamma_lam_hhssi),
     breaks = 50,
     main = expression("HHS-SI: " * log[10](lambda[beta]^2 * lambda[gamma]^2)),
     xlab = expression(log[10](lambda[beta]^2 * lambda[gamma]^2)))
abline(v = 0, lty = 2, lwd = 2)





delta_var_mean_hhssi <- apply(delta_var_hhssi, 2, mean)
beta_gamma_lam_mean_hhssi <- apply(beta_gamma_lam_hhssi, 2, mean)
# plot(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi,
#   log = "xy",
#   pch = 19,
#   xlab = "Delta prior var of hhssi from hhssi",
#   ylab = "Delta prior var of nhhssi from hhssi",
#   xlim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)),
#   ylim = range(c(delta_var_mean_hhssi, beta_gamma_lam_mean_hhssi)))
# abline(0, 1, lty = 2, lwd = 2, col = "red")


delta_var_mean_nhhssi <- apply(delta_var_nhhssi, 2, mean)
plot(delta_var_mean_hhssi, delta_var_mean_nhhssi,
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)),
     ylim = range(c(delta_var_mean_hhssi, delta_var_mean_nhhssi)))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_nhhssi,
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi,
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)



plot(delta_var_mean_nhhssi[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red")
lines(delta_var_mean_nhhssi[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


plot(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[-c(1:3, 18:20)], delta_var_mean_nhhssi[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)],
     log = "xy",
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_hhssi[c(1:3, 18:20)], delta_var_mean_nhhssi[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")






delta_var_mean_hhssi[c(1:3, 18:20)] > delta_var_mean_nhhssi[c(1:3, 18:20)]
delta_var_mean_hhssi[-c(1:3, 18:20)] > delta_var_mean_nhhssi[-c(1:3, 18:20)]


ratio_hhssi_nhhssi <- delta_var_hhssi / delta_var_nhhssi
ratio_mean <- apply(ratio_hhssi_nhhssi, 2, mean)
hist(log10(ratio_mean), breaks = 50,
     main = "HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[delta,HHS-SI] / V[delta,NHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


active_delta_var_mean_hhssi <- delta_var_mean_hhssi[c(1:3, 18:20)]
active_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[c(1:3, 18:20)]
inactive_delta_var_mean_hhssi <- delta_var_mean_hhssi[-c(1:3, 18:20)]
inactive_delta_var_mean_nhhssi <- delta_var_mean_nhhssi[-c(1:3, 18:20)]
active_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,c(1:3, 18:20)] /
                                          delta_var_nhhssi[,c(1:3, 18:20)])
inactive_delta_var_ratio_hhssi_nhhssi <- (delta_var_hhssi[,-c(1:3, 18:20)] /
                                            delta_var_nhhssi[,-c(1:3, 18:20)])
active_delta_var_ratio_mean_hhssi_nhhssi <- apply(active_delta_var_ratio_hhssi_nhhssi, 2, mean)
inactive_delta_var_ratio_mean_hhssi_nhhssi <- apply(inactive_delta_var_ratio_hhssi_nhhssi, 2, mean)

hist(log10(active_delta_var_ratio_mean_hhssi_nhhssi), breaks = 5,
     main = "active HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[activedeltaHHS-SI] / V[activedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")


hist(log10(inactive_delta_var_ratio_mean_hhssi_nhhssi), breaks = 50,
     main = "inactive HHS-SI / NHHS-SI delta variance ratio",
     xlab = expression(log[10](V[inactivedeltaHHS-SI] / V[inactivedeltaNHHS-SI])))
abline(v = 0, lty = 2, lwd = 2, col = "red")









plot(beta_gamma_lam_mean_hhssi,
     pch = 19,
     log = "y",
     xlab = "Interaction index",
     ylab = expression(E[posterior](lambda[beta]^2 * lambda[gamma]^2)),
     main = "HHS-SI beta-gamma multiplier")
abline(h = 1, lty = 2, lwd = 2, col = "red")













###############################################################################################
# Combined plot
###############################################################################################


# First sample
set.seed(23444)
sim_data_normal01_simpred1 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n ,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)
# load(file="practice_data/res_normal01_simpred_hhssi1.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi1.rda")
set.seed(23444)
sim_data_normal01_realpred1 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                 O=O, x=X, z=Z, d=D,
                                                 Sigma=Sigma)
# load(file="practice_data/res_normal01_realpred_hhssi1.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi1.rda")





# Second sample
set.seed(654586798)
sim_data_normal01_simpred2 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n ,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)
# load(file="practice_data/res_normal01_simpred_hhssi2.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi2.rda")
set.seed(654586798)
sim_data_normal01_realpred2 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                 O=O, x=X, z=Z, d=D,
                                                 Sigma=Sigma)
# load(file="practice_data/res_normal01_realpred_hhssi2.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi2.rda")





# Third sample
set.seed(23315644)
sim_data_normal01_simpred3 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n ,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)
# load(file="practice_data/res_normal01_simpred_hhssi3.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi3.rda")
set.seed(23315644)
sim_data_normal01_realpred3 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                 O=O, x=X, z=Z, d=D,
                                                 Sigma=Sigma)
# load(file="practice_data/res_normal01_realpred_hhssi3.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi3.rda")







#fs sim
beta_lam_simpred_hhssi1  <- res_normal01_simpred_hhssi1$lambdasq_beta_update
gamma_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_gamma_update
delta_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_delta_update
tausq_delta_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update

delta_lam_simpred_nhhssi1  <- res_normal01_simpred_nhhssi1$lambdasq_delta_update
tausq_delta_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update

#fs real
beta_lam_realpred_hhssi1  <- res_normal01_realpred_hhssi1$lambdasq_beta_update
gamma_lam_realpred_hhssi1 <- res_normal01_realpred_hhssi1$lambdasq_gamma_update
delta_lam_realpred_hhssi1 <- res_normal01_realpred_hhssi1$lambdasq_delta_update
tausq_delta_realpred_hhssi1 <- res_normal01_realpred_hhssi1$tausq_delta_update

delta_lam_realpred_nhhssi1  <- res_normal01_realpred_nhhssi1$lambdasq_delta_update
tausq_delta_realpred_nhhssi1 <- res_normal01_realpred_nhhssi1$tausq_delta_update


#ss sim
beta_lam_simpred_hhssi2  <- res_normal01_simpred_hhssi2$lambdasq_beta_update
gamma_lam_simpred_hhssi2 <- res_normal01_simpred_hhssi2$lambdasq_gamma_update
delta_lam_simpred_hhssi2 <- res_normal01_simpred_hhssi2$lambdasq_delta_update
tausq_delta_simpred_hhssi2 <- res_normal01_simpred_hhssi2$tausq_delta_update

delta_lam_simpred_nhhssi2  <- res_normal01_simpred_nhhssi2$lambdasq_delta_update
tausq_delta_simpred_nhhssi2 <- res_normal01_simpred_nhhssi2$tausq_delta_update

#ss real
beta_lam_realpred_hhssi2  <- res_normal01_realpred_hhssi2$lambdasq_beta_update
gamma_lam_realpred_hhssi2 <- res_normal01_realpred_hhssi2$lambdasq_gamma_update
delta_lam_realpred_hhssi2 <- res_normal01_realpred_hhssi2$lambdasq_delta_update
tausq_delta_realpred_hhssi2 <- res_normal01_realpred_hhssi2$tausq_delta_update

delta_lam_realpred_nhhssi2  <- res_normal01_realpred_nhhssi2$lambdasq_delta_update
tausq_delta_realpred_nhhssi2 <- res_normal01_realpred_nhhssi2$tausq_delta_update


#ts sim
beta_lam_simpred_hhssi3  <- res_normal01_simpred_hhssi3$lambdasq_beta_update
gamma_lam_simpred_hhssi3 <- res_normal01_simpred_hhssi3$lambdasq_gamma_update
delta_lam_simpred_hhssi3 <- res_normal01_simpred_hhssi3$lambdasq_delta_update
tausq_delta_simpred_hhssi3 <- res_normal01_simpred_hhssi3$tausq_delta_update

delta_lam_simpred_nhhssi3  <- res_normal01_simpred_nhhssi3$lambdasq_delta_update
tausq_delta_simpred_nhhssi3 <- res_normal01_simpred_nhhssi3$tausq_delta_update

#ts real
beta_lam_realpred_hhssi3  <- res_normal01_realpred_hhssi3$lambdasq_beta_update
gamma_lam_realpred_hhssi3 <- res_normal01_realpred_hhssi3$lambdasq_gamma_update
delta_lam_realpred_hhssi3 <- res_normal01_realpred_hhssi3$lambdasq_delta_update
tausq_delta_realpred_hhssi3 <- res_normal01_realpred_hhssi3$tausq_delta_update

delta_lam_realpred_nhhssi3  <- res_normal01_realpred_nhhssi3$lambdasq_delta_update
tausq_delta_realpred_nhhssi3 <- res_normal01_realpred_nhhssi3$tausq_delta_update




n_iter <- dim(beta_lam_realpred_hhssi3)[1]



####################################################################
# for normal01_simpred_hhssi1, normal01_simpred_nhhssi1

# Construct beta x gamma component
beta_gamma_lam_simpred_hhssi1 <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_simpred_hhssi1[t, col_idx] <- (beta_lam_simpred_hhssi1[t, j] *
                                                      gamma_lam_simpred_hhssi1[t, m])
      col_idx <- col_idx + 1
    }
  }
}
beta_gamma_lam_mean_simpred_hhssi1 <- apply(beta_gamma_lam_simpred_hhssi1, 2, mean)

# prior variance from hhssi
delta_var_simpred_hhssi1 <- (delta_lam_simpred_hhssi1 *
                               beta_gamma_lam_simpred_hhssi1 *
                               tausq_delta_simpred_hhssi1)
delta_var_mean_simpred_hhssi1 <- apply(delta_var_simpred_hhssi1, 2, mean)
dim(delta_var_simpred_hhssi1)
# prior variance from nhhssi
delta_var_simpred_nhhssi1 <- (delta_lam_simpred_nhhssi1 *
                                tausq_delta_simpred_nhhssi1)
delta_var_mean_simpred_nhhssi1 <- apply(delta_var_simpred_nhhssi1, 2, mean)


# prior var estimates from hhssi vs nhhssi
plot(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1),
     pch = 19,
     xlab = "log HHS-SI delta var mean",
     ylab = "log NHHS-SI delta var mean",
     main ="log transformed",
     xlim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))),
     ylim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)),
     ylim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
plot(delta_var_mean_simpred_hhssi1[idx],
     delta_var_mean_simpred_nhhssi1[idx],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_simpred_hhssi1[idx],
                    delta_var_mean_simpred_nhhssi1[idx])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[idx],
                    delta_var_mean_simpred_nhhssi1[idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")

plot(delta_var_mean_simpred_hhssi1[!idx],
     delta_var_mean_simpred_nhhssi1[!idx],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                    delta_var_mean_simpred_nhhssi1[!idx])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                    delta_var_mean_simpred_nhhssi1[!idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




# for active delta
plot(delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     main="for corresponding active deltas",
     xlim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")



# for inactive delta
plot(delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     main="for corresponding inactive deltas",
     xlim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")



inactive_simpred_hhssi1 <- delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)]
inactive_simpred_nhhssi1 <- delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
inxinc <- inactive_simpred_nhhssi1 > inactive_simpred_hhssi1
plot(inactive_simpred_hhssi1[inxinc],
     inactive_simpred_nhhssi1[inxinc],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(inactive_simpred_hhssi1[inxinc],
                    inactive_simpred_hhssi1[inxinc])),
     ylim = range(c(inactive_simpred_hhssi1[inxinc],
                    inactive_simpred_hhssi1[inxinc])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")

plot(inactive_simpred_hhssi1[!inxinc],
     inactive_simpred_nhhssi1[!inxinc],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(inactive_simpred_hhssi1[!inxinc],
                    inactive_simpred_hhssi1[!inxinc])),
     ylim = range(c(inactive_simpred_hhssi1[!inxinc],
                    inactive_simpred_hhssi1[!inxinc])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



# so the prior variance of lambdasq which are bigger from nhhssi model than hhssi model,
# comes from the corresponding inactive delta values


# other variance terms
# lambdasq_beta
lambdsq_beta_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdsq_beta_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update, 2, mean)
plot(lambdsq_beta_mean_simpred_hhssi1,
     lambdsq_beta_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq beta mean",
     ylab = "NHHS-SI lambdasq beta mean",
     main="overall",
     xlim = range(c(lambdsq_beta_mean_simpred_hhssi1,
                    lambdsq_beta_mean_simpred_nhhssi1)),
     ylim = range(c(lambdsq_beta_mean_simpred_hhssi1,
                    lambdsq_beta_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_beta_mean_simpred_hhssi1[1:2],
     lambdsq_beta_mean_simpred_nhhssi1[1:2],
     pch = 19,
     xlab = "HHS-SI lambdasq beta mean",
     ylab = "NHHS-SI lambdasq beta mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_beta_mean_simpred_hhssi1[1:2],
                    lambdsq_beta_mean_simpred_nhhssi1[1:2])),
     ylim = range(c(lambdsq_beta_mean_simpred_hhssi1[1:2],
                    lambdsq_beta_mean_simpred_nhhssi1[1:2])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_beta_mean_simpred_hhssi1[-(1:2)],
     lambdsq_beta_mean_simpred_nhhssi1[-(1:2)],
     pch = 19,
     xlab = "HHS-SI lambdasq beta mean",
     ylab = "NHHS-SI lambdasq beta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_beta_mean_simpred_hhssi1[-(1:2)],
                    lambdsq_beta_mean_simpred_nhhssi1[-(1:2)])),
     ylim = range(c(lambdsq_beta_mean_simpred_hhssi1[-(1:2)],
                    lambdsq_beta_mean_simpred_nhhssi1[-(1:2)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




# lambdasq_gamma
lambdsq_gamma_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdsq_gamma_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update, 2, mean)
plot(lambdsq_gamma_mean_simpred_hhssi1,
     lambdsq_gamma_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq gamma mean",
     ylab = "NHHS-SI lambdasq gamma mean",
     main="overall",
     xlim = range(c(lambdsq_gamma_mean_simpred_hhssi1,
                    lambdsq_gamma_mean_simpred_nhhssi1)),
     ylim = range(c(lambdsq_gamma_mean_simpred_hhssi1,
                    lambdsq_gamma_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_gamma_mean_simpred_hhssi1[1:3],
     lambdsq_gamma_mean_simpred_nhhssi1[1:3],
     pch = 19,
     xlab = "HHS-SI lambdasq gamma mean",
     ylab = "NHHS-SI lambdasq gamma mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_gamma_mean_simpred_hhssi1[1:3],
                    lambdsq_gamma_mean_simpred_nhhssi1[1:3])),
     ylim = range(c(lambdsq_gamma_mean_simpred_hhssi1[1:3],
                    lambdsq_gamma_mean_simpred_nhhssi1[1:3])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_gamma_mean_simpred_hhssi1[-(1:3)],
     lambdsq_gamma_mean_simpred_nhhssi1[-(1:3)],
     pch = 19,
     xlab = "HHS-SI lambdasq gamma mean",
     ylab = "NHHS-SI lambdasq gamma mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_gamma_mean_simpred_hhssi1[-(1:3)],
                    lambdsq_gamma_mean_simpred_nhhssi1[-(1:3)])),
     ylim = range(c(lambdsq_gamma_mean_simpred_hhssi1[-(1:3)],
                    lambdsq_gamma_mean_simpred_nhhssi1[-(1:3)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")





# lambdasq_delta
lambdsq_delta_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdsq_delta_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_delta_update, 2, mean)
plot(lambdsq_delta_mean_simpred_hhssi1,
     lambdsq_delta_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="overall",
     xlim = range(c(lambdsq_delta_mean_simpred_hhssi1,
                    lambdsq_delta_mean_simpred_nhhssi1)),
     ylim = range(c(lambdsq_delta_mean_simpred_hhssi1,
                    lambdsq_delta_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
     lambdsq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(lambdsq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
     lambdsq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))])),
     ylim = range(c(lambdsq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




# tausq_beta
tausq_beta_mean_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_beta_update
tausq_beta_mean_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_beta_update
plot(tausq_beta_mean_simpred_hhssi1,
     tausq_beta_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI tausq beta",
     ylab = "NHHS-SI tausq beta",
     main="overall",
     xlim = range(c(tausq_beta_mean_simpred_hhssi1,
                    tausq_beta_mean_simpred_nhhssi1)),
     ylim = range(c(tausq_beta_mean_simpred_hhssi1,
                    tausq_beta_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_beta_mean_simpred_hhssi1[1:2],
     tausq_beta_mean_simpred_nhhssi1[1:2],
     pch = 19,
     xlab = "HHS-SI tausq beta",
     ylab = "NHHS-SI tausq beta",
     main="for corresponding active parameters",
     xlim = range(c(tausq_beta_mean_simpred_hhssi1[1:2],
                    tausq_beta_mean_simpred_nhhssi1[1:2])),
     ylim = range(c(tausq_beta_mean_simpred_hhssi1[1:2],
                    tausq_beta_mean_simpred_nhhssi1[1:2])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_beta_mean_simpred_hhssi1[-(1:2)],
     tausq_beta_mean_simpred_nhhssi1[-(1:2)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(tausq_beta_mean_simpred_hhssi1[-(1:2)],
                    tausq_beta_mean_simpred_nhhssi1[-(1:2)])),
     ylim = range(c(tausq_beta_mean_simpred_hhssi1[-(1:2)],
                    tausq_beta_mean_simpred_nhhssi1[-(1:2)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



# tausq_gamma
tausq_gamma_mean_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_gamma_update
tausq_gamma_mean_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_gamma_update
plot(tausq_gamma_mean_simpred_hhssi1,
     tausq_gamma_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI tausq gamma",
     ylab = "NHHS-SI tausq gamma",
     main="overall",
     xlim = range(c(tausq_gamma_mean_simpred_hhssi1,
                    tausq_gamma_mean_simpred_nhhssi1)),
     ylim = range(c(tausq_gamma_mean_simpred_hhssi1,
                    tausq_gamma_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_gamma_mean_simpred_hhssi1[1:3],
     tausq_gamma_mean_simpred_nhhssi1[1:3],
     pch = 19,
     xlab = "HHS-SI tausq gamma",
     ylab = "NHHS-SI tausq gamma",
     main="for corresponding active parameters",
     xlim = range(c(tausq_gamma_mean_simpred_hhssi1[1:3],
                    tausq_gamma_mean_simpred_nhhssi1[1:3])),
     ylim = range(c(tausq_gamma_mean_simpred_hhssi1[1:3],
                    tausq_gamma_mean_simpred_nhhssi1[1:3])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_gamma_mean_simpred_hhssi1[-(1:3)],
     tausq_gamma_mean_simpred_nhhssi1[-(1:3)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(tausq_gamma_mean_simpred_hhssi1[-(1:3)],
                    tausq_gamma_mean_simpred_nhhssi1[-(1:3)])),
     ylim = range(c(tausq_gamma_mean_simpred_hhssi1[-(1:3)],
                    tausq_gamma_mean_simpred_nhhssi1[-(1:3)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")







# tausq_delta
tausq_delta_mean_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update
tausq_delta_mean_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update
plot(tausq_delta_mean_simpred_hhssi1,
     tausq_delta_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI tausq delta",
     ylab = "NHHS-SI tausq delta",
     main="overall",
     xlim = range(c(tausq_delta_mean_simpred_hhssi1,
                    tausq_delta_mean_simpred_nhhssi1)),
     ylim = range(c(tausq_delta_mean_simpred_hhssi1,
                    tausq_delta_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
     tausq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI tausq delta",
     ylab = "NHHS-SI tausq delta",
     main="for corresponding active parameters",
     xlim = range(c(tausq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
                    tausq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(tausq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
                    tausq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
     tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(tausq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
                    tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))])),
     ylim = range(c(tausq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
                    tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")






# lambdasq_delta*tausq_delta
lambdsq_tausq_delta_mean_simpred_hhssi1 <- apply((res_normal01_simpred_hhssi1$lambdasq_delta_update*
                                                 res_normal01_simpred_hhssi1$tausq_delta_update),
                                                 2, mean)
lambdsq_tausq_delta_mean_simpred_nhhssi1 <- apply((res_normal01_simpred_nhhssi1$lambdasq_delta_update*
                                                     res_normal01_simpred_nhhssi1$tausq_delta_update),
                                                  2, mean)
plot(lambdsq_tausq_delta_mean_simpred_hhssi1,
     lambdsq_tausq_delta_mean_simpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq*tausq delta mean",
     ylab = "NHHS-SI lambdasq*tausq delta mean",
     main="overall",
     xlim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1,
                    lambdsq_tausq_delta_mean_simpred_nhhssi1)),
     ylim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1,
                    lambdsq_tausq_delta_mean_simpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_tausq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
     lambdsq_tausq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_tausq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_tausq_delta_mean_simpred_nhhssi1[c(1:3, 18:20)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_tausq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
     lambdsq_tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))])),
     ylim = range(c(lambdsq_tausq_delta_mean_simpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_tausq_delta_mean_simpred_nhhssi1[-(c(1:3, 18:20))])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")








#########################
# for normal01_realpred_hhssi1, normal01_realpred_nhhssi1

# Construct beta x gamma component
beta_gamma_lam_realpred_hhssi1 <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_realpred_hhssi1[t, col_idx] <- (beta_lam_realpred_hhssi1[t, j] *
                                                      gamma_lam_realpred_hhssi1[t, m])
      col_idx <- col_idx + 1
    }
  }
}
beta_gamma_lam_mean_realpred_hhssi1 <- apply(beta_gamma_lam_realpred_hhssi1, 2, mean)

# prior variance from hhssi
delta_var_realpred_hhssi1 <- (delta_lam_realpred_hhssi1 *
                               beta_gamma_lam_realpred_hhssi1 *
                               tausq_delta_realpred_hhssi1)
delta_var_mean_realpred_hhssi1 <- apply(delta_var_realpred_hhssi1, 2, mean)
dim(delta_var_realpred_hhssi1)
# prior variance from nhhssi
delta_var_realpred_nhhssi1 <- (delta_lam_realpred_nhhssi1 *
                                tausq_delta_realpred_nhhssi1)
delta_var_mean_realpred_nhhssi1 <- apply(delta_var_realpred_nhhssi1, 2, mean)


# prior var estimates from hhssi vs nhhssi
plot(log(delta_var_mean_realpred_hhssi1), log(delta_var_mean_realpred_nhhssi1),
     pch = 19,
     xlab = "log HHS-SI delta var mean",
     ylab = "log NHHS-SI delta var mean",
     main ="log transformed",
     xlim = range(c(log(delta_var_mean_realpred_hhssi1), log(delta_var_mean_realpred_nhhssi1))),
     ylim = range(c(log(delta_var_mean_realpred_hhssi1), log(delta_var_mean_realpred_nhhssi1))))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(delta_var_mean_realpred_hhssi1, delta_var_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_realpred_hhssi1, delta_var_mean_realpred_nhhssi1)),
     ylim = range(c(delta_var_mean_realpred_hhssi1, delta_var_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




idx <- delta_var_mean_realpred_nhhssi1 > delta_var_mean_realpred_hhssi1
plot(delta_var_mean_realpred_hhssi1[idx],
     delta_var_mean_realpred_nhhssi1[idx],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_realpred_hhssi1[idx],
                    delta_var_mean_realpred_nhhssi1[idx])),
     ylim = range(c(delta_var_mean_realpred_hhssi1[idx],
                    delta_var_mean_realpred_nhhssi1[idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


# for active delta
plot(delta_var_mean_realpred_hhssi1[c(1:3, 18:20)],
     type = "b",
     xlab = "Parameter",
     lwd = 1,
     col = "red",
     main="for corresponding active deltas",
     ylim=range(c(delta_var_mean_realpred_hhssi1[c(1:3, 18:20)],delta_var_mean_realpred_nhhssi1[c(1:3, 18:20)])))
lines(delta_var_mean_realpred_nhhssi1[c(1:3, 18:20)],
      type = "b",
      lwd = 1,
      lty = 2,
      col = "blue")
legend("topright",
       legend = c("HHS-SI", "NHHS-SI"),
       col = c("red", "blue"),
       lty = c(1, 2),
       pch = c(1, 1),
       lwd = 1)


# for inactive delta
plot(delta_var_mean_realpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_realpred_nhhssi1[-c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     main="for corresponding inactive deltas",
     xlim = range(c(delta_var_mean_realpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_realpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_realpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_realpred_nhhssi1[-c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")

inactive_realpred_hhssi1 <- delta_var_mean_realpred_hhssi1[-c(1:3, 18:20)]
inactive_realpred_nhhssi1 <- delta_var_mean_realpred_nhhssi1[-c(1:3, 18:20)]
inxinc <- delta_var_mean_realpred_nhhssi1[-c(1:3, 18:20)] > delta_var_mean_realpred_hhssi1[-c(1:3, 18:20)]
plot(inactive_realpred_hhssi1[inxinc],
     inactive_realpred_nhhssi1[inxinc],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(inactive_realpred_hhssi1[inxinc],
                    inactive_realpred_hhssi1[inxinc])),
     ylim = range(c(inactive_realpred_hhssi1[inxinc],
                    inactive_realpred_hhssi1[inxinc])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



# so the prior variance of lambdasq which are bigger from nhhssi model than hhssi model,
# comes from the corresponding inactive delta values


# other variance terms
# lambdasq_beta
lambdsq_beta_mean_realpred_hhssi1 <- apply(res_normal01_realpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdsq_beta_mean_realpred_nhhssi1 <- apply(res_normal01_realpred_nhhssi1$lambdasq_beta_update, 2, mean)
plot(lambdsq_beta_mean_realpred_hhssi1,
     lambdsq_beta_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq beta mean",
     ylab = "NHHS-SI lambdasq beta mean",
     main="overall",
     xlim = range(c(lambdsq_beta_mean_realpred_hhssi1,
                    lambdsq_beta_mean_realpred_nhhssi1)),
     ylim = range(c(lambdsq_beta_mean_realpred_hhssi1,
                    lambdsq_beta_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_beta_mean_realpred_hhssi1[1:2],
     lambdsq_beta_mean_realpred_nhhssi1[1:2],
     pch = 19,
     xlab = "HHS-SI lambdasq beta mean",
     ylab = "NHHS-SI lambdasq beta mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_beta_mean_realpred_hhssi1[1:2],
                    lambdsq_beta_mean_realpred_nhhssi1[1:2])),
     ylim = range(c(lambdsq_beta_mean_realpred_hhssi1[1:2],
                    lambdsq_beta_mean_realpred_nhhssi1[1:2])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_beta_mean_realpred_hhssi1[-(1:2)],
     lambdsq_beta_mean_realpred_nhhssi1[-(1:2)],
     pch = 19,
     xlab = "HHS-SI lambdasq beta mean",
     ylab = "NHHS-SI lambdasq beta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_beta_mean_realpred_hhssi1[-(1:2)],
                    lambdsq_beta_mean_realpred_nhhssi1[-(1:2)])),
     ylim = range(c(lambdsq_beta_mean_realpred_hhssi1[-(1:2)],
                    lambdsq_beta_mean_realpred_nhhssi1[-(1:2)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




# lambdasq_gamma
lambdsq_gamma_mean_realpred_hhssi1 <- apply(res_normal01_realpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdsq_gamma_mean_realpred_nhhssi1 <- apply(res_normal01_realpred_nhhssi1$lambdasq_gamma_update, 2, mean)
plot(lambdsq_gamma_mean_realpred_hhssi1,
     lambdsq_gamma_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq gamma mean",
     ylab = "NHHS-SI lambdasq gamma mean",
     main="overall",
     xlim = range(c(lambdsq_gamma_mean_realpred_hhssi1,
                    lambdsq_gamma_mean_realpred_nhhssi1)),
     ylim = range(c(lambdsq_gamma_mean_realpred_hhssi1,
                    lambdsq_gamma_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_gamma_mean_realpred_hhssi1[1:3],
     lambdsq_gamma_mean_realpred_nhhssi1[1:3],
     pch = 19,
     xlab = "HHS-SI lambdasq gamma mean",
     ylab = "NHHS-SI lambdasq gamma mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_gamma_mean_realpred_hhssi1[1:3],
                    lambdsq_gamma_mean_realpred_nhhssi1[1:3])),
     ylim = range(c(lambdsq_gamma_mean_realpred_hhssi1[1:3],
                    lambdsq_gamma_mean_realpred_nhhssi1[1:3])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_gamma_mean_realpred_hhssi1[-(1:3)],
     lambdsq_gamma_mean_realpred_nhhssi1[-(1:3)],
     pch = 19,
     xlab = "HHS-SI lambdasq gamma mean",
     ylab = "NHHS-SI lambdasq gamma mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_gamma_mean_realpred_hhssi1[-(1:3)],
                    lambdsq_gamma_mean_realpred_nhhssi1[-(1:3)])),
     ylim = range(c(lambdsq_gamma_mean_realpred_hhssi1[-(1:3)],
                    lambdsq_gamma_mean_realpred_nhhssi1[-(1:3)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")





# lambdasq_delta
lambdsq_delta_mean_realpred_hhssi1 <- apply(res_normal01_realpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdsq_delta_mean_realpred_nhhssi1 <- apply(res_normal01_realpred_nhhssi1$lambdasq_delta_update, 2, mean)
plot(lambdsq_delta_mean_realpred_hhssi1,
     lambdsq_delta_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="overall",
     xlim = range(c(lambdsq_delta_mean_realpred_hhssi1,
                    lambdsq_delta_mean_realpred_nhhssi1)),
     ylim = range(c(lambdsq_delta_mean_realpred_hhssi1,
                    lambdsq_delta_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
     lambdsq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(lambdsq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
     lambdsq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))])),
     ylim = range(c(lambdsq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




# tausq_beta
tausq_beta_mean_realpred_hhssi1 <- res_normal01_realpred_hhssi1$tausq_beta_update
tausq_beta_mean_realpred_nhhssi1 <- res_normal01_realpred_nhhssi1$tausq_beta_update
plot(tausq_beta_mean_realpred_hhssi1,
     tausq_beta_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI tausq beta",
     ylab = "NHHS-SI tausq beta",
     main="overall",
     xlim = range(c(tausq_beta_mean_realpred_hhssi1,
                    tausq_beta_mean_realpred_nhhssi1)),
     ylim = range(c(tausq_beta_mean_realpred_hhssi1,
                    tausq_beta_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_beta_mean_realpred_hhssi1[1:2],
     tausq_beta_mean_realpred_nhhssi1[1:2],
     pch = 19,
     xlab = "HHS-SI tausq beta",
     ylab = "NHHS-SI tausq beta",
     main="for corresponding active parameters",
     xlim = range(c(tausq_beta_mean_realpred_hhssi1[1:2],
                    tausq_beta_mean_realpred_nhhssi1[1:2])),
     ylim = range(c(tausq_beta_mean_realpred_hhssi1[1:2],
                    tausq_beta_mean_realpred_nhhssi1[1:2])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_beta_mean_realpred_hhssi1[-(1:2)],
     tausq_beta_mean_realpred_nhhssi1[-(1:2)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(tausq_beta_mean_realpred_hhssi1[-(1:2)],
                    tausq_beta_mean_realpred_nhhssi1[-(1:2)])),
     ylim = range(c(tausq_beta_mean_realpred_hhssi1[-(1:2)],
                    tausq_beta_mean_realpred_nhhssi1[-(1:2)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



# tausq_gamma
tausq_gamma_mean_realpred_hhssi1 <- res_normal01_realpred_hhssi1$tausq_gamma_update
tausq_gamma_mean_realpred_nhhssi1 <- res_normal01_realpred_nhhssi1$tausq_gamma_update
plot(tausq_gamma_mean_realpred_hhssi1,
     tausq_gamma_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI tausq gamma",
     ylab = "NHHS-SI tausq gamma",
     main="overall",
     xlim = range(c(tausq_gamma_mean_realpred_hhssi1,
                    tausq_gamma_mean_realpred_nhhssi1)),
     ylim = range(c(tausq_gamma_mean_realpred_hhssi1,
                    tausq_gamma_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_gamma_mean_realpred_hhssi1[1:3],
     tausq_gamma_mean_realpred_nhhssi1[1:3],
     pch = 19,
     xlab = "HHS-SI tausq gamma",
     ylab = "NHHS-SI tausq gamma",
     main="for corresponding active parameters",
     xlim = range(c(tausq_gamma_mean_realpred_hhssi1[1:3],
                    tausq_gamma_mean_realpred_nhhssi1[1:3])),
     ylim = range(c(tausq_gamma_mean_realpred_hhssi1[1:3],
                    tausq_gamma_mean_realpred_nhhssi1[1:3])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_gamma_mean_realpred_hhssi1[-(1:3)],
     tausq_gamma_mean_realpred_nhhssi1[-(1:3)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(tausq_gamma_mean_realpred_hhssi1[-(1:3)],
                    tausq_gamma_mean_realpred_nhhssi1[-(1:3)])),
     ylim = range(c(tausq_gamma_mean_realpred_hhssi1[-(1:3)],
                    tausq_gamma_mean_realpred_nhhssi1[-(1:3)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")







# tausq_delta
tausq_delta_mean_realpred_hhssi1 <- res_normal01_realpred_hhssi1$tausq_delta_update
tausq_delta_mean_realpred_nhhssi1 <- res_normal01_realpred_nhhssi1$tausq_delta_update
plot(tausq_delta_mean_realpred_hhssi1,
     tausq_delta_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI tausq delta",
     ylab = "NHHS-SI tausq delta",
     main="overall",
     xlim = range(c(tausq_delta_mean_realpred_hhssi1,
                    tausq_delta_mean_realpred_nhhssi1)),
     ylim = range(c(tausq_delta_mean_realpred_hhssi1,
                    tausq_delta_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
     tausq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI tausq delta",
     ylab = "NHHS-SI tausq delta",
     main="for corresponding active parameters",
     xlim = range(c(tausq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
                    tausq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(tausq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
                    tausq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(tausq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
     tausq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(tausq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
                    tausq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))])),
     ylim = range(c(tausq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
                    tausq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")






# lambdasq_delta*tausq_delta
lambdsq_tausq_delta_mean_realpred_hhssi1 <- apply((res_normal01_realpred_hhssi1$lambdasq_delta_update*
                                                    res_normal01_realpred_hhssi1$tausq_delta_update),
                                                 2, mean)
lambdsq_tausq_delta_mean_realpred_nhhssi1 <- apply((res_normal01_realpred_nhhssi1$lambdasq_delta_update*
                                                     res_normal01_realpred_nhhssi1$tausq_delta_update),
                                                  2, mean)
plot(lambdsq_tausq_delta_mean_realpred_hhssi1,
     lambdsq_tausq_delta_mean_realpred_nhhssi1,
     pch = 19,
     xlab = "HHS-SI lambdasq*tausq delta mean",
     ylab = "NHHS-SI lambdasq*tausq delta mean",
     main="overall",
     xlim = range(c(lambdsq_tausq_delta_mean_realpred_hhssi1,
                    lambdsq_tausq_delta_mean_realpred_nhhssi1)),
     ylim = range(c(lambdsq_tausq_delta_mean_realpred_hhssi1,
                    lambdsq_tausq_delta_mean_realpred_nhhssi1)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_tausq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
     lambdsq_tausq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding active parameters",
     xlim = range(c(lambdsq_tausq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_tausq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(lambdsq_tausq_delta_mean_realpred_hhssi1[c(1:3, 18:20)],
                    lambdsq_tausq_delta_mean_realpred_nhhssi1[c(1:3, 18:20)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")


plot(lambdsq_tausq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
     lambdsq_tausq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))],
     pch = 19,
     xlab = "HHS-SI lambdasq delta mean",
     ylab = "NHHS-SI lambdasq delta mean",
     main="for corresponding inactive parameters",
     xlim = range(c(lambdsq_tausq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_tausq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))])),
     ylim = range(c(lambdsq_tausq_delta_mean_realpred_hhssi1[-(c(1:3, 18:20))],
                    lambdsq_tausq_delta_mean_realpred_nhhssi1[-(c(1:3, 18:20))])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



#######################################################
# combined plot from three samples
#####################################################

####################################################################
# Function to construct beta x gamma component

get_beta_gamma <- function(beta_lam, gamma_lam, J, M) {

  n_iter <- nrow(beta_lam)

  beta_gamma_lam <- matrix(NA, nrow = n_iter, ncol = J * M)

  for (t in 1:n_iter) {
    col_idx <- 1

    for (j in 1:J) {
      for (m in 1:M) {

        beta_gamma_lam[t, col_idx] <-
          beta_lam[t, j] * gamma_lam[t, m]

        col_idx <- col_idx + 1
      }
    }
  }

  return(beta_gamma_lam)
}



# Simulated predictors
beta_gamma_simpred_hhssi1 <-
  get_beta_gamma(beta_lam_simpred_hhssi1,
                 gamma_lam_simpred_hhssi1, J, M)

beta_gamma_simpred_hhssi2 <-
  get_beta_gamma(beta_lam_simpred_hhssi2,
                 gamma_lam_simpred_hhssi2, J, M)

beta_gamma_simpred_hhssi3 <-
  get_beta_gamma(beta_lam_simpred_hhssi3,
                 gamma_lam_simpred_hhssi3, J, M)


# Real predictors
beta_gamma_realpred_hhssi1 <-
  get_beta_gamma(beta_lam_realpred_hhssi1,
                 gamma_lam_realpred_hhssi1, J, M)

beta_gamma_realpred_hhssi2 <-
  get_beta_gamma(beta_lam_realpred_hhssi2,
                 gamma_lam_realpred_hhssi2, J, M)

beta_gamma_realpred_hhssi3 <-
  get_beta_gamma(beta_lam_realpred_hhssi3,
                 gamma_lam_realpred_hhssi3, J, M)



delta_var_simpred_hhssi1 <-
  delta_lam_simpred_hhssi1 *
  beta_gamma_simpred_hhssi1 *
  tausq_delta_simpred_hhssi1

delta_var_simpred_hhssi2 <-
  delta_lam_simpred_hhssi2 *
  beta_gamma_simpred_hhssi2 *
  tausq_delta_simpred_hhssi2

delta_var_simpred_hhssi3 <-
  delta_lam_simpred_hhssi3 *
  beta_gamma_simpred_hhssi3 *
  tausq_delta_simpred_hhssi3


delta_var_simpred_nhhssi1 <-
  delta_lam_simpred_nhhssi1 *
  tausq_delta_simpred_nhhssi1

delta_var_simpred_nhhssi2 <-
  delta_lam_simpred_nhhssi2 *
  tausq_delta_simpred_nhhssi2

delta_var_simpred_nhhssi3 <-
  delta_lam_simpred_nhhssi3 *
  tausq_delta_simpred_nhhssi3


delta_var_realpred_hhssi1 <-
  delta_lam_realpred_hhssi1 *
  beta_gamma_realpred_hhssi1 *
  tausq_delta_realpred_hhssi1

delta_var_realpred_hhssi2 <-
  delta_lam_realpred_hhssi2 *
  beta_gamma_realpred_hhssi2 *
  tausq_delta_realpred_hhssi2

delta_var_realpred_hhssi3 <-
  delta_lam_realpred_hhssi3 *
  beta_gamma_realpred_hhssi3 *
  tausq_delta_realpred_hhssi3


delta_var_realpred_nhhssi1 <-
  delta_lam_realpred_nhhssi1 *
  tausq_delta_realpred_nhhssi1

delta_var_realpred_nhhssi2 <-
  delta_lam_realpred_nhhssi2 *
  tausq_delta_realpred_nhhssi2

delta_var_realpred_nhhssi3 <-
  delta_lam_realpred_nhhssi3 *
  tausq_delta_realpred_nhhssi3


# Simulated predictors
delta_var_mean_simpred_hhssi1 <- apply(delta_var_simpred_hhssi1, 2, mean)
delta_var_mean_simpred_hhssi2 <- apply(delta_var_simpred_hhssi2, 2, mean)
delta_var_mean_simpred_hhssi3 <- apply(delta_var_simpred_hhssi3, 2, mean)

delta_var_mean_simpred_nhhssi1 <- apply(delta_var_simpred_nhhssi1, 2, mean)
delta_var_mean_simpred_nhhssi2 <- apply(delta_var_simpred_nhhssi2, 2, mean)
delta_var_mean_simpred_nhhssi3 <- apply(delta_var_simpred_nhhssi3, 2, mean)


# Average across the 3 samples
delta_var_mean_simpred_hhssi_avg <-
  rowMeans(cbind(delta_var_mean_simpred_hhssi1,
                 delta_var_mean_simpred_hhssi2,
                 delta_var_mean_simpred_hhssi3))

delta_var_mean_simpred_nhhssi_avg <-
  rowMeans(cbind(delta_var_mean_simpred_nhhssi1,
                 delta_var_mean_simpred_nhhssi2,
                 delta_var_mean_simpred_nhhssi3))


delta_var_mean_realpred_hhssi1 <- apply(delta_var_realpred_hhssi1, 2, mean)
delta_var_mean_realpred_hhssi2 <- apply(delta_var_realpred_hhssi2, 2, mean)
delta_var_mean_realpred_hhssi3 <- apply(delta_var_realpred_hhssi3, 2, mean)

delta_var_mean_realpred_nhhssi1 <- apply(delta_var_realpred_nhhssi1, 2, mean)
delta_var_mean_realpred_nhhssi2 <- apply(delta_var_realpred_nhhssi2, 2, mean)
delta_var_mean_realpred_nhhssi3 <- apply(delta_var_realpred_nhhssi3, 2, mean)


delta_var_mean_realpred_hhssi_avg <-
  rowMeans(cbind(delta_var_mean_realpred_hhssi1,
                 delta_var_mean_realpred_hhssi2,
                 delta_var_mean_realpred_hhssi3))

delta_var_mean_realpred_nhhssi_avg <-
  rowMeans(cbind(delta_var_mean_realpred_nhhssi1,
                 delta_var_mean_realpred_nhhssi2,
                 delta_var_mean_realpred_nhhssi3))


delta_var_mean_simpred_hhssi_avg
delta_var_mean_simpred_nhhssi_avg

delta_var_mean_realpred_hhssi_avg
delta_var_mean_realpred_nhhssi_avg




plot(log(delta_var_mean_simpred_hhssi_avg),
     log(delta_var_mean_simpred_nhhssi_avg),
     pch = 19,
     xlab = "log HHS-SI delta variance",
     ylab = "log NHHS-SI delta variance",
     main = "Average across 3 samples",
     xlim = range(c(log(delta_var_mean_simpred_hhssi_avg),
                    log(delta_var_mean_simpred_nhhssi_avg))),
     ylim = range(c(log(delta_var_mean_simpred_hhssi_avg),
                    log(delta_var_mean_simpred_nhhssi_avg))))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(delta_var_mean_simpred_hhssi_avg,
     delta_var_mean_simpred_nhhssi_avg,
     pch = 19,
     xlab = "log HHS-SI delta variance",
     ylab = "log NHHS-SI delta variance",
     main = "Average across 3 samples",
     xlim = range(c(delta_var_mean_simpred_hhssi_avg,
                   delta_var_mean_simpred_nhhssi_avg)),
     ylim = range(c(delta_var_mean_simpred_hhssi_avg,
                    delta_var_mean_simpred_nhhssi_avg)))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")





idx <- delta_var_mean_simpred_hhssi_avg < delta_var_mean_simpred_nhhssi_avg
plot(delta_var_mean_simpred_hhssi_avg[idx],
     delta_var_mean_simpred_nhhssi_avg[idx],
     pch = 19,
     xlab = "log HHS-SI delta variance",
     ylab = "log NHHS-SI delta variance",
     main = "Average across 3 samples",
     xlim = range(c(delta_var_mean_simpred_hhssi_avg[idx],
                    delta_var_mean_simpred_nhhssi_avg[idx])),
     ylim = range(c(delta_var_mean_simpred_hhssi_avg[idx],
                    delta_var_mean_simpred_nhhssi_avg[idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(delta_var_mean_simpred_hhssi_avg[!idx],
     delta_var_mean_simpred_nhhssi_avg[!idx],
     pch = 19,
     xlab = "log HHS-SI delta variance",
     ylab = "log NHHS-SI delta variance",
     main = "Average across 3 samples",
     xlim = range(c(delta_var_mean_simpred_hhssi_avg[!idx],
                    delta_var_mean_simpred_nhhssi_avg[!idx])),
     ylim = range(c(delta_var_mean_simpred_hhssi_avg[!idx],
                    delta_var_mean_simpred_nhhssi_avg[!idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(delta_var_mean_simpred_hhssi_avg[c(1:3,18:20)],
     delta_var_mean_simpred_nhhssi_avg[c(1:3,18:20)],
     pch = 19,
     xlab = "log HHS-SI delta variance",
     ylab = "log NHHS-SI delta variance",
     main = "Average across 3 samples",
     xlim = range(c(delta_var_mean_simpred_hhssi_avg[c(1:3,18:20)],
                    delta_var_mean_simpred_nhhssi_avg[c(1:3,18:20)])),
     ylim = range(c(delta_var_mean_simpred_hhssi_avg[c(1:3,18:20)],
                    delta_var_mean_simpred_nhhssi_avg[c(1:3,18:20)])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



idav <- which(delta_var_mean_simpred_nhhssi_avg>delta_var_mean_simpred_hhssi_avg)
id1 <- which(delta_var_mean_simpred_nhhssi1>delta_var_mean_simpred_hhssi1)
id2 <- which(delta_var_mean_simpred_nhhssi2>delta_var_mean_simpred_hhssi2)
id3 <- which(delta_var_mean_simpred_nhhssi3>delta_var_mean_simpred_hhssi3)

length(id1)
length(id2)
length(id3)
length(idav)

id_match <- intersect(id1, id2)

length(id_match)
id_match

which(!id_match)
####################################################################################
# delta prior variance estimates
####################################################################
# for normal01_simpred_hhssi1, normal01_simpred_nhhssi1

#fs sim
beta_lam_simpred_hhssi1  <- res_normal01_simpred_hhssi1$lambdasq_beta_update
gamma_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_gamma_update
delta_lam_simpred_hhssi1 <- res_normal01_simpred_hhssi1$lambdasq_delta_update
tausq_delta_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update

delta_lam_simpred_nhhssi1  <- res_normal01_simpred_nhhssi1$lambdasq_delta_update
tausq_delta_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update

n_iter <- dim(beta_lam_simpred_hhssi1)[1]

# Construct beta x gamma component
beta_gamma_lam_simpred_hhssi1 <- matrix(NA, nrow = n_iter, ncol = J * M)
for (t in 1:n_iter)
{
  col_idx <- 1
  for (j in 1:J) {
    for (m in 1:M) {
      beta_gamma_lam_simpred_hhssi1[t, col_idx] <- (beta_lam_simpred_hhssi1[t, j] *
                                                      gamma_lam_simpred_hhssi1[t, m])
      col_idx <- col_idx + 1
    }
  }
}
beta_gamma_lam_mean_simpred_hhssi1 <- apply(beta_gamma_lam_simpred_hhssi1, 2, mean)

# prior variance from hhssi
delta_var_simpred_hhssi1 <- (delta_lam_simpred_hhssi1 *
                               beta_gamma_lam_simpred_hhssi1 *
                               tausq_delta_simpred_hhssi1)
delta_var_mean_simpred_hhssi1 <- apply(delta_var_simpred_hhssi1, 2, mean)
dim(delta_var_simpred_hhssi1)
# prior variance from nhhssi
delta_var_simpred_nhhssi1 <- (delta_lam_simpred_nhhssi1 *
                                tausq_delta_simpred_nhhssi1)
delta_var_mean_simpred_nhhssi1 <- apply(delta_var_simpred_nhhssi1, 2, mean)


# prior var estimates from hhssi vs nhhssi
idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
point_col <- ifelse(idx, "blue", "red")

par(mfrow=c(1,1))
plot(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1),
     pch = 19,
     col = point_col,
     xlab = "log HHS-SI delta var mean",
     ylab = "log NHHS-SI delta var mean",
     main ="log transformed",
     xlim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))),
     ylim = range(c(log(delta_var_mean_simpred_hhssi1), log(delta_var_mean_simpred_nhhssi1))))
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




plot(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1,
     pch = 19,
     col = point_col,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)),
     ylim = range(c(delta_var_mean_simpred_hhssi1, delta_var_mean_simpred_nhhssi1)))
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
plot(delta_var_mean_simpred_hhssi1[idx],
     delta_var_mean_simpred_nhhssi1[idx],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_simpred_hhssi1[idx],
                    delta_var_mean_simpred_nhhssi1[idx])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[idx],
                    delta_var_mean_simpred_nhhssi1[idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(delta_var_mean_simpred_hhssi1[!idx],
     delta_var_mean_simpred_nhhssi1[!idx],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     xlim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                    delta_var_mean_simpred_nhhssi1[!idx])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[!idx],
                    delta_var_mean_simpred_nhhssi1[!idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




# for active delta
plot(delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     main="for corresponding active deltas",
     xlim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")



# for inactive delta
idx_inct <- delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)] > delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)]
point_col_inct <- ifelse(idx_inct, "blue", "red")
plot(delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)],
     pch = 19,
     col= point_col_inct,
     xlab = "HHS-SI delta var mean",
     ylab = "NHHS-SI delta var mean",
     main="for corresponding inactive deltas",
     xlim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)])))
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(0, 1, lty = 2, lwd = 2, col = "red")


#######################################################################################
# sd of prior variance estimate
sqrt_delta_var_mean_simpred_hhssi1 <- sqrt(delta_var_mean_simpred_hhssi1)
sqrt_delta_var_mean_simpred_nhhssi1 <- sqrt(delta_var_mean_simpred_nhhssi1)

plot(sqrt_delta_var_mean_simpred_hhssi1, sqrt_delta_var_mean_simpred_nhhssi1,
     pch = 19,
     col = point_col,
     xlab = "HHS-SI sqrt of delta var",
     ylab = "NHHS-SI sqrt of delta var",
     xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1, sqrt_delta_var_mean_simpred_nhhssi1)),
     ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1, sqrt_delta_var_mean_simpred_nhhssi1)))
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




idx <- sqrt_delta_var_mean_simpred_nhhssi1 > sqrt_delta_var_mean_simpred_hhssi1
plot(sqrt_delta_var_mean_simpred_hhssi1[idx],
     sqrt_delta_var_mean_simpred_nhhssi1[idx],
     pch = 19,
     xlab = "HHS-SI sqrt of delta var",
     ylab = "NHHS-SI sqrt of delta var",
     xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[idx],
                    sqrt_delta_var_mean_simpred_nhhssi1[idx])),
     ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[idx],
                    sqrt_delta_var_mean_simpred_nhhssi1[idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")



plot(sqrt_delta_var_mean_simpred_hhssi1[!idx],
     sqrt_delta_var_mean_simpred_nhhssi1[!idx],
     pch = 19,
     xlab = "HHS-SI sqrt of delta var",
     ylab = "NHHS-SI sqrt of delta var",
     xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[!idx],
                    sqrt_delta_var_mean_simpred_nhhssi1[!idx])),
     ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[!idx],
                    sqrt_delta_var_mean_simpred_nhhssi1[!idx])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")




# for active delta
plot(sqrt_delta_var_mean_simpred_hhssi1[c(1:3, 18:20)],
     sqrt_delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)],
     pch = 19,
     xlab = "HHS-SI sqrt of delta var",
     ylab = "NHHS-SI sqrt of delta var",
     main="for corresponding active deltas",
     xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], sqrt_delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[c(1:3, 18:20)], sqrt_delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)])))
abline(0, 1, lty = 2, lwd = 2, col = "red")



# for inactive delta
plot(sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)],
     sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)],
     pch = 19,
     col = point_col_inct,
     xlab = "HHS-SI sqrt of delta var",
     ylab = "NHHS-SI sqrt of delta var",
     main="for corresponding inactive deltas",
     xlim = range(c(sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)], sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)])))
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(0, 1, lty = 2, lwd = 2, col = "red")



sqrt_inactive_simpred_hhssi1 <- sqrt_delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)]
sqrt_inactive_simpred_nhhssi1 <- sqrt_delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
inxinc <- sqrt_inactive_simpred_nhhssi1 > sqrt_inactive_simpred_hhssi1
plot(sqrt_inactive_simpred_hhssi1[inxinc],
     sqrt_inactive_simpred_nhhssi1[inxinc],
     pch = 19,
     xlab = "HHS-SI sqrt of delta var",
     ylab = "NHHS-SI sqrt of delta var",
     xlim = range(c(sqrt_inactive_simpred_hhssi1[inxinc],
                    sqrt_inactive_simpred_hhssi1[inxinc])),
     ylim = range(c(sqrt_inactive_simpred_hhssi1[inxinc],
                    sqrt_inactive_simpred_hhssi1[inxinc])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")

plot(sqrt_inactive_simpred_hhssi1[!inxinc],
     sqrt_inactive_simpred_nhhssi1[!inxinc],
     pch = 19,
     xlab = "HHS-SI sqrt of delta var",
     ylab = "NHHS-SI sqrt of delta var",
     xlim = range(c(sqrt_inactive_simpred_hhssi1[!inxinc],
                    sqrt_inactive_simpred_hhssi1[!inxinc])),
     ylim = range(c(sqrt_inactive_simpred_hhssi1[!inxinc],
                    sqrt_inactive_simpred_hhssi1[!inxinc])))
abline(a = 0, b = 1, lty = 2, lwd = 2, col = "red")





#############################################################################################
# Posterior sd estimate of delta
delta_sd_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 63:810, ], c(2,3), sd)
delta_sd_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 63:810, ], c(2,3), sd)

delta_est_var_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 63:810, ], c(2,3), var)
delta_est_var_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 63:810, ], c(2,3), var)



# Posterior mean estimate of delta
delta_est_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 63:810, ], c(2,3), mean)
delta_est_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 63:810, ], c(2,3), mean)

# Bias for each outcome
delta_bias_normal01_simpred_hhssi1 <- (delta_est_normal01_simpred_hhssi1 - sim_data_normal01_simpred1$delta_true)
delta_bias_normal01_simpred_nhhssi1 <- (delta_est_normal01_simpred_nhhssi1 - sim_data_normal01_simpred1$delta_true)



mean(delta_est_var_normal01_simpred_hhssi1)
mean(delta_est_var_normal01_simpred_nhhssi1)

mean(delta_bias_normal01_simpred_hhssi1)
mean(delta_bias_normal01_simpred_nhhssi1)

##########################################################################################


# sd of Delta Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_normal01_simpred_hhssi1[, k],
    delta_sd_normal01_simpred_nhhssi1[, k],
    col =point_col,
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Overall sd: Outcome", k),
    xlim = range(c(delta_sd_normal01_simpred_hhssi1[,k],
                   delta_sd_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(delta_sd_normal01_simpred_hhssi1[,k],
                   delta_sd_normal01_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")



# active delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_normal01_simpred_hhssi1[c(1:3, 18:20), k],
    delta_sd_normal01_simpred_nhhssi1[c(1:3, 18:20), k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("active: Outcome", k),
    xlim = range(c(delta_sd_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_sd_normal01_simpred_nhhssi1[c(1:3, 18:20),k])),
    ylim = range(c(delta_sd_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_sd_normal01_simpred_nhhssi1[c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20), k],
    delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20), k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    col = point_col_inct,
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])),
    ylim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
}






idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
# idx
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_normal01_simpred_hhssi1[idx, k],
    delta_sd_normal01_simpred_nhhssi1[idx, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_sd_normal01_simpred_hhssi1[idx,k],
                   delta_sd_normal01_simpred_nhhssi1[idx,k])),
    ylim = range(c(delta_sd_normal01_simpred_hhssi1[idx,k],
                   delta_sd_normal01_simpred_nhhssi1[idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


# !idx
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_normal01_simpred_hhssi1[!idx, k],
    delta_sd_normal01_simpred_nhhssi1[!idx, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_sd_normal01_simpred_hhssi1[!idx,k],
                   delta_sd_normal01_simpred_nhhssi1[!idx,k])),
    ylim = range(c(delta_sd_normal01_simpred_hhssi1[!idx,k],
                   delta_sd_normal01_simpred_nhhssi1[!idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





idx_act <- delta_var_mean_simpred_hhssi1[c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)]
# no active delta with nhhssi var bigger than hhssi

# inactive delta
idx_inact <- delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20), ][idx_inact, k],
    delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20), ][idx_inact, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Outcome", k),
    xlim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20),][idx_inact, k],
                   delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20),][idx_inact, k])),
    ylim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20),][idx_inact, k],
                   delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20),][idx_inact, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20), ][!idx_inact, k],
    delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20), ][!idx_inact, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Outcome", k),
    xlim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20),][!idx_inact, k],
                   delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20),][!idx_inact, k])),
    ylim = range(c(delta_sd_normal01_simpred_hhssi1[-c(1:3, 18:20),][!idx_inact, k],
                   delta_sd_normal01_simpred_nhhssi1[-c(1:3, 18:20),][!idx_inact, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





##########################################################################################

# variance of Delta Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_var_normal01_simpred_hhssi1[, k],
    delta_est_var_normal01_simpred_nhhssi1[, k],
    col =point_col,
    xlab = "HHSSI delta var",
    ylab = "NHHSSI delta var",
    main = paste("Overall sd: Outcome", k),
    xlim = range(c(delta_est_var_normal01_simpred_hhssi1[,k],
                   delta_est_var_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(delta_est_var_normal01_simpred_hhssi1[,k],
                   delta_est_var_normal01_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")



# active delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_var_normal01_simpred_hhssi1[c(1:3, 18:20), k],
    delta_est_var_normal01_simpred_nhhssi1[c(1:3, 18:20), k],
    xlab = "HHSSI delta var",
    ylab = "NHHSSI delta var",
    main = paste("active: Outcome", k),
    xlim = range(c(delta_est_var_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_est_var_normal01_simpred_nhhssi1[c(1:3, 18:20),k])),
    ylim = range(c(delta_est_var_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_est_var_normal01_simpred_nhhssi1[c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_var_normal01_simpred_hhssi1[-c(1:3, 18:20), k],
    delta_est_var_normal01_simpred_nhhssi1[-c(1:3, 18:20), k],
    xlab = "HHSSI delta var",
    ylab = "NHHSSI delta var",
    col = point_col_inct,
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_est_var_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_est_var_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])),
    ylim = range(c(delta_est_var_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_est_var_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
}





##########################################################################################


# Delta Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_normal01_simpred_hhssi1[, k],
    delta_est_normal01_simpred_nhhssi1[, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    col = point_col,
    main = paste("Overall estimates: Outcome", k),
    xlim = range(c(delta_est_normal01_simpred_hhssi1[,k],
                   delta_est_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(delta_est_normal01_simpred_hhssi1[,k],
                   delta_est_normal01_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
}

# active delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_normal01_simpred_hhssi1[c(1:3, 18:20), k],
    delta_est_normal01_simpred_nhhssi1[c(1:3, 18:20), k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("active: Outcome", k),
    xlim = range(c(delta_est_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_est_normal01_simpred_nhhssi1[c(1:3, 18:20),k])),
    ylim = range(c(delta_est_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_est_normal01_simpred_nhhssi1[c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20), k],
    delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20), k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    col = point_col_inct,
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])),
    ylim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])))
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0,1, lty = 2, lwd = 1, col = "red")
}






idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_normal01_simpred_hhssi1[idx, k],
    delta_est_normal01_simpred_nhhssi1[idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_est_normal01_simpred_hhssi1[idx,k],
                   delta_est_normal01_simpred_nhhssi1[idx,k])),
    ylim = range(c(delta_est_normal01_simpred_hhssi1[idx,k],
                   delta_est_normal01_simpred_nhhssi1[idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_normal01_simpred_hhssi1[!idx, k],
    delta_est_normal01_simpred_nhhssi1[!idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_est_normal01_simpred_hhssi1[!idx,k],
                   delta_est_normal01_simpred_nhhssi1[!idx,k])),
    ylim = range(c(delta_est_normal01_simpred_hhssi1[!idx,k],
                   delta_est_normal01_simpred_nhhssi1[!idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





idx <- delta_var_mean_simpred_hhssi1[c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)]
# no active delta with nhhssi var bigger than hhssi

# inactive delta
idx <- delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20), ][idx, k],
    delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20), ][idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])),
    ylim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20), ][!idx, k],
    delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20), ][!idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])),
    ylim = range(c(delta_est_normal01_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_est_normal01_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





##########################################################################################
# Delta Bias
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_normal01_simpred_hhssi1[, k],
    delta_bias_normal01_simpred_nhhssi1[, k],
    col = point_col,
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Overall bias: Outcome", k),
    xlim = range(c(delta_bias_normal01_simpred_hhssi1[,k],
                   delta_bias_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(delta_bias_normal01_simpred_hhssi1[,k],
                   delta_bias_normal01_simpred_nhhssi1[,k])))
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




# active delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_normal01_simpred_hhssi1[c(1:3, 18:20), k],
    delta_bias_normal01_simpred_nhhssi1[c(1:3, 18:20), k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("active: Outcome", k),
    xlim = range(c(delta_bias_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_bias_normal01_simpred_nhhssi1[c(1:3, 18:20),k])),
    ylim = range(c(delta_bias_normal01_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_bias_normal01_simpred_nhhssi1[c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20), k],
    delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20), k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    col = point_col_inct,
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])),
    ylim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
  legend("topleft",
         legend = c(
           "NHHSSI variance > HHSSI variance",
           "NHHSSI variance <= HHSSI variance"),
         col = c("blue", "red"),
         pch = 19,
         bty = "n")
}






idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_normal01_simpred_hhssi1[idx, k],
    delta_bias_normal01_simpred_nhhssi1[idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_bias_normal01_simpred_hhssi1[idx,k],
                   delta_bias_normal01_simpred_nhhssi1[idx,k])),
    ylim = range(c(delta_bias_normal01_simpred_hhssi1[idx,k],
                   delta_bias_normal01_simpred_nhhssi1[idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_normal01_simpred_hhssi1[!idx, k],
    delta_bias_normal01_simpred_nhhssi1[!idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_bias_normal01_simpred_hhssi1[!idx,k],
                   delta_bias_normal01_simpred_nhhssi1[!idx,k])),
    ylim = range(c(delta_bias_normal01_simpred_hhssi1[!idx,k],
                   delta_bias_normal01_simpred_nhhssi1[!idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}






idx <- delta_var_mean_simpred_hhssi1[c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)]
# no active delta with nhhssi var bigger than hhssi

# inactive delta
idx <- delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20), ][idx, k],
    delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20), ][idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Outcome", k),
    xlim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])),
    ylim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20), ][!idx, k],
    delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20), ][!idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Outcome", k),
    xlim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])),
    ylim = range(c(delta_bias_normal01_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_bias_normal01_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}



######################################################################################

idx <- delta_var_mean_simpred_hhssi1 < delta_var_mean_simpred_nhhssi1
# delta_idx <- which(!idx)
delta_idx <- which(idx)
# delta_idx <- 1:748
beta_idx <- ceiling(delta_idx / M)
gamma_idx <- ((delta_idx - 1) %% M) + 1
data.frame(
  delta_index = delta_idx,
  beta_index = beta_idx,
  gamma_index = gamma_idx
)



length(beta_idx)
length(gamma_idx)

unique_beta_idx <- unique(beta_idx) # they are all inactive beta for NHHSSI>HHSSI var
unique_gamma_idx <- unique(gamma_idx) # they are all inactive gamma for NHHSSI>HHSSI var


length(unique_beta_idx)
length(unique_gamma_idx)
###########################

# # this diag_dat for inactive with nhhssi smaller hhssi only
# # Keep only terms with inactive beta and gamma indices
# keep <- !(beta_idx %in% c(1, 2) |
#             gamma_idx %in% c(1, 2, 3))
# delta_idx <- delta_idx[keep]
# beta_idx  <- beta_idx[keep]
# gamma_idx <- gamma_idx[keep]
# diag_dat <- data.frame(
#   delta = delta_idx,
#   beta = beta_idx,
#   gamma = gamma_idx
# )


diag_dat <- data.frame(
  delta = delta_idx,
  beta  = beta_idx,
  gamma = gamma_idx
)
beta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$theta_update[, 2:45, ]
gamma_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$theta_update[, 46:62, ]
delta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$theta_update[, 63:810, ]


beta_mean_hhssi <- apply(beta_normal01_simpred_hhssi1, c(2, 3), mean)
gamma_mean_hhssi <- apply(gamma_normal01_simpred_hhssi1, c(2, 3), mean)
delta_mean_hhssi <- apply(delta_normal01_simpred_hhssi1, c(2, 3), mean)

diag_dat$delta_est1 <- delta_mean_hhssi[delta_idx, 1]
diag_dat$delta_est2 <- delta_mean_hhssi[delta_idx, 2]
diag_dat$delta_est3 <- delta_mean_hhssi[delta_idx, 3]
diag_dat$delta_est4 <- delta_mean_hhssi[delta_idx, 4]

diag_dat$beta_est1 <- beta_mean_hhssi[beta_idx, 1]
diag_dat$beta_est2 <- beta_mean_hhssi[beta_idx, 2]
diag_dat$beta_est3 <- beta_mean_hhssi[beta_idx, 3]
diag_dat$beta_est4 <- beta_mean_hhssi[beta_idx, 4]


diag_dat$gamma_est1 <- gamma_mean_hhssi[gamma_idx, 1]
diag_dat$gamma_est2 <- gamma_mean_hhssi[gamma_idx, 2]
diag_dat$gamma_est3 <- gamma_mean_hhssi[gamma_idx, 3]
diag_dat$gamma_est4 <- gamma_mean_hhssi[gamma_idx, 4]

dim(diag_dat)
head(diag_dat)

library(dplyr)
library(tidyr)

diag_long <- diag_dat %>%
  pivot_longer(
    cols = matches("^(delta|beta|gamma)_est[1-4]$"),
    names_to = c("parameter", "outcome"),
    names_pattern = "(delta|beta|gamma)_est([1-4])",
    values_to = "estimate"
  )

diag_long <- data.frame(
  delta = rep(diag_dat$delta, 4),
  beta  = rep(diag_dat$beta, 4),
  gamma = rep(diag_dat$gamma, 4),
  outcome = rep(1:4, each = nrow(diag_dat)),
  delta_est = c(
    diag_dat$delta_est1,
    diag_dat$delta_est2,
    diag_dat$delta_est3,
    diag_dat$delta_est4
  ),
  beta_est = c(
    diag_dat$beta_est1,
    diag_dat$beta_est2,
    diag_dat$beta_est3,
    diag_dat$beta_est4
  ),
  gamma_est = c(
    diag_dat$gamma_est1,
    diag_dat$gamma_est2,
    diag_dat$gamma_est3,
    diag_dat$gamma_est4
  )
)


diag_long$beta_gamma_product <-
  diag_long$beta_est * diag_long$gamma_est
par(mfrow=c(1,1))
par(mfrow = c(2, 2))
for (k in 1:4) {
  tmp <- subset(diag_long, outcome == k)
  plot(
    tmp$beta_est * tmp$gamma_est,
    tmp$delta_est,
    xlab = expression(hat(beta) %*% hat(gamma)),
    ylab = expression(hat(delta)),
    main = paste("Outcome", k),
    pch = 16,
    cex = 0.6
  )
abline(h = 0, v = 0, lty = 2)
}


##########################################################################################

# Posterior sd estimate of beta
beta_sd_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 2:45, ], c(2,3), sd)
beta_sd_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 2:45, ], c(2,3), sd)

# Posterior mean estimate of beta
beta_est_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 2:45, ], c(2,3), mean)
beta_est_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 2:45, ], c(2,3), mean)

# Bias for each outcome
beta_bias_normal01_simpred_hhssi1 <- (beta_est_normal01_simpred_hhssi1 - sim_data_normal01_simpred1$beta_true)
beta_bias_normal01_simpred_nhhssi1 <- (beta_est_normal01_simpred_nhhssi1 - sim_data_normal01_simpred1$beta_true)

# prior var estimate of beta
beta_var_normal01_simpred_hhssi1 <- apply((res_normal01_simpred_hhssi1$lambdasq_beta_update*
                                            res_normal01_simpred_hhssi1$tausq_beta_update), 2, mean)
beta_var_normal01_simpred_nhhssi1 <- apply((res_normal01_simpred_nhhssi1$lambdasq_beta_update*
                                              res_normal01_simpred_nhhssi1$tausq_beta_update), 2, mean)


# beta Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_normal01_simpred_hhssi1[, k],
    beta_est_normal01_simpred_nhhssi1[, k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("Overall estimates: Outcome", k),
    xlim = range(c(beta_est_normal01_simpred_hhssi1[,k],
                   beta_est_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(beta_est_normal01_simpred_hhssi1[,k],
                   beta_est_normal01_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_normal01_simpred_hhssi1[1:2, k],
    beta_est_normal01_simpred_nhhssi1[1:2, k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("active: Outcome", k),
    xlim = range(c(beta_est_normal01_simpred_hhssi1[1:2,k],
                   beta_est_normal01_simpred_nhhssi1[1:2,k])),
    ylim = range(c(beta_est_normal01_simpred_hhssi1[1:2,k],
                   beta_est_normal01_simpred_nhhssi1[1:2,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_normal01_simpred_hhssi1[-c(1:2), k],
    beta_est_normal01_simpred_nhhssi1[-c(1:2), k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(beta_est_normal01_simpred_hhssi1[-c(1:2),k],
                   beta_est_normal01_simpred_nhhssi1[-c(1:2),k])),
    ylim = range(c(beta_est_normal01_simpred_hhssi1[-c(1:2),k],
                   beta_est_normal01_simpred_nhhssi1[-c(1:2),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}







# all the indices are all inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_normal01_simpred_hhssi1[unique_beta_idx, k],
    beta_est_normal01_simpred_nhhssi1[unique_beta_idx, k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("Outcome", k),
    xlim = range(c(beta_est_normal01_simpred_hhssi1[unique_beta_idx,k],
                   beta_est_normal01_simpred_nhhssi1[unique_beta_idx,k])),
    ylim = range(c(beta_est_normal01_simpred_hhssi1[unique_beta_idx,k],
                   beta_est_normal01_simpred_nhhssi1[unique_beta_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




##########################################################################################
# beta Bias
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_normal01_simpred_hhssi1[, k],
    beta_bias_normal01_simpred_nhhssi1[, k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("Overall bias: Outcome", k),
    xlim = range(c(beta_bias_normal01_simpred_hhssi1[,k],
                   beta_bias_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(beta_bias_normal01_simpred_hhssi1[,k],
                   beta_bias_normal01_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_normal01_simpred_hhssi1[c(1:2), k],
    beta_bias_normal01_simpred_nhhssi1[c(1:2), k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("active: Outcome", k),
    xlim = range(c(beta_bias_normal01_simpred_hhssi1[c(1:2),k],
                   beta_bias_normal01_simpred_nhhssi1[c(1:2),k])),
    ylim = range(c(beta_bias_normal01_simpred_hhssi1[c(1:2),k],
                   beta_bias_normal01_simpred_nhhssi1[c(1:2),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_normal01_simpred_hhssi1[-c(1:2), k],
    beta_bias_normal01_simpred_nhhssi1[-c(1:2), k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("inactive: Outcome", k),
    xlim = range(c(beta_bias_normal01_simpred_hhssi1[-c(1:2),k],
                   beta_bias_normal01_simpred_nhhssi1[-c(1:2),k])),
    ylim = range(c(beta_bias_normal01_simpred_hhssi1[-c(1:2),k],
                   beta_bias_normal01_simpred_nhhssi1[-c(1:2),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}



# all the indices are all inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_normal01_simpred_hhssi1[unique_beta_idx, k],
    beta_bias_normal01_simpred_nhhssi1[unique_beta_idx, k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(beta_bias_normal01_simpred_hhssi1[unique_beta_idx,k],
                   beta_bias_normal01_simpred_nhhssi1[unique_beta_idx,k])),
    ylim = range(c(beta_bias_normal01_simpred_hhssi1[unique_beta_idx,k],
                   beta_bias_normal01_simpred_nhhssi1[unique_beta_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




######################################################################################

# beta prior variance estimate
# overall
par(mfrow = c(1,1))
plot(beta_var_normal01_simpred_hhssi1,
    beta_var_normal01_simpred_nhhssi1,
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_normal01_simpred_hhssi1,
                   beta_var_normal01_simpred_nhhssi1)),
    ylim = range(c(beta_var_normal01_simpred_hhssi1,
                   beta_var_normal01_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


# active beta
plot(beta_var_normal01_simpred_hhssi1[c(1:2)],
    beta_var_normal01_simpred_nhhssi1[c(1:2)],
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_normal01_simpred_hhssi1[c(1:2)],
                   beta_var_normal01_simpred_nhhssi1[c(1:2)])),
    ylim = range(c(beta_var_normal01_simpred_hhssi1[c(1:2)],
                   beta_var_normal01_simpred_nhhssi1[c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive beta
plot(beta_var_normal01_simpred_hhssi1[-c(1:2)],
    beta_var_normal01_simpred_nhhssi1[-c(1:2)],
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_normal01_simpred_hhssi1[-c(1:2)],
                   beta_var_normal01_simpred_nhhssi1[-c(1:2)])),
    ylim = range(c(beta_var_normal01_simpred_hhssi1[-c(1:2)],
                   beta_var_normal01_simpred_nhhssi1[-c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive beta
plot(beta_var_normal01_simpred_hhssi1[unique_beta_idx],
    beta_var_normal01_simpred_nhhssi1[unique_beta_idx],
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_normal01_simpred_hhssi1[unique_beta_idx],
                   beta_var_normal01_simpred_nhhssi1[unique_beta_idx])),
    ylim = range(c(beta_var_normal01_simpred_hhssi1[unique_beta_idx],
                   beta_var_normal01_simpred_nhhssi1[unique_beta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")





######################################################################################
# Posterior sd estimate of gamma
gamma_sd_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 46:62, ], c(2,3), sd)
gamma_sd_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 46:62, ], c(2,3), sd)

# Posterior mean estimate of gamma
gamma_est_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$theta_update[, 46:62, ], c(2,3), mean)

# Bias for each outcome
gamma_bias_normal01_simpred_hhssi1 <- (gamma_est_normal01_simpred_hhssi1 - sim_data_normal01_simpred1$gamma_true)
gamma_bias_normal01_simpred_nhhssi1 <- (gamma_est_normal01_simpred_nhhssi1 - sim_data_normal01_simpred1$gamma_true)

# prior var estimate of gamma
gamma_var_normal01_simpred_hhssi1 <- apply((res_normal01_simpred_hhssi1$lambdasq_gamma_update*
                                             res_normal01_simpred_hhssi1$tausq_gamma_update), 2, mean)
gamma_var_normal01_simpred_nhhssi1 <- apply((res_normal01_simpred_nhhssi1$lambdasq_gamma_update*
                                              res_normal01_simpred_nhhssi1$tausq_gamma_update), 2, mean)


# gamma Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_normal01_simpred_hhssi1[, k],
    gamma_est_normal01_simpred_nhhssi1[, k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("Overall estimates: Outcome", k),
    xlim = range(c(gamma_est_normal01_simpred_hhssi1[,k],
                   gamma_est_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(gamma_est_normal01_simpred_hhssi1[,k],
                   gamma_est_normal01_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_normal01_simpred_hhssi1[1:3, k],
    gamma_est_normal01_simpred_nhhssi1[1:3, k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("active: Outcome", k),
    xlim = range(c(gamma_est_normal01_simpred_hhssi1[1:3,k],
                   gamma_est_normal01_simpred_nhhssi1[1:3,k])),
    ylim = range(c(gamma_est_normal01_simpred_hhssi1[1:3,k],
                   gamma_est_normal01_simpred_nhhssi1[1:3,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_normal01_simpred_hhssi1[-c(1:3), k],
    gamma_est_normal01_simpred_nhhssi1[-c(1:3), k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(gamma_est_normal01_simpred_hhssi1[-c(1:3),k],
                   gamma_est_normal01_simpred_nhhssi1[-c(1:3),k])),
    ylim = range(c(gamma_est_normal01_simpred_hhssi1[-c(1:3),k],
                   gamma_est_normal01_simpred_nhhssi1[-c(1:3),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}







# overall (infact all are inactive gamma)
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_normal01_simpred_hhssi1[unique_gamma_idx, k],
    gamma_est_normal01_simpred_nhhssi1[unique_gamma_idx, k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("Outcome", k),
    xlim = range(c(gamma_est_normal01_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_est_normal01_simpred_nhhssi1[unique_gamma_idx,k])),
    ylim = range(c(gamma_est_normal01_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_est_normal01_simpred_nhhssi1[unique_gamma_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




##########################################################################################
# gamma Bias
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_normal01_simpred_hhssi1[, k],
    gamma_bias_normal01_simpred_nhhssi1[, k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("Overall bias: Outcome", k),
    xlim = range(c(gamma_bias_normal01_simpred_hhssi1[,k],
                   gamma_bias_normal01_simpred_nhhssi1[,k])),
    ylim = range(c(gamma_bias_normal01_simpred_hhssi1[,k],
                   gamma_bias_normal01_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_normal01_simpred_hhssi1[c(1:3), k],
    gamma_bias_normal01_simpred_nhhssi1[c(1:3), k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("active: Outcome", k),
    xlim = range(c(gamma_bias_normal01_simpred_hhssi1[c(1:3),k],
                   gamma_bias_normal01_simpred_nhhssi1[c(1:3),k])),
    ylim = range(c(gamma_bias_normal01_simpred_hhssi1[c(1:3),k],
                   gamma_bias_normal01_simpred_nhhssi1[c(1:3),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_normal01_simpred_hhssi1[-c(1:3), k],
    gamma_bias_normal01_simpred_nhhssi1[-c(1:3), k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("inactive: Outcome", k),
    xlim = range(c(gamma_bias_normal01_simpred_hhssi1[-c(1:3),k],
                   gamma_bias_normal01_simpred_nhhssi1[-c(1:3),k])),
    ylim = range(c(gamma_bias_normal01_simpred_hhssi1[-c(1:3),k],
                   gamma_bias_normal01_simpred_nhhssi1[-c(1:3),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}



# overall (infact all are inactive gamma)
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_normal01_simpred_hhssi1[unique_gamma_idx, k],
    gamma_bias_normal01_simpred_nhhssi1[unique_gamma_idx, k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(gamma_bias_normal01_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_bias_normal01_simpred_nhhssi1[unique_gamma_idx,k])),
    ylim = range(c(gamma_bias_normal01_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_bias_normal01_simpred_nhhssi1[unique_gamma_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




######################################################################################
# gamma prior variance estimate
# overall
par(mfrow = c(1,1))
plot(gamma_var_normal01_simpred_hhssi1,
     gamma_var_normal01_simpred_nhhssi1,
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_normal01_simpred_hhssi1,
                    gamma_var_normal01_simpred_nhhssi1)),
     ylim = range(c(gamma_var_normal01_simpred_hhssi1,
                    gamma_var_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active gamma
plot(gamma_var_normal01_simpred_hhssi1[c(1:3)],
     gamma_var_normal01_simpred_nhhssi1[c(1:3)],
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_normal01_simpred_hhssi1[c(1:3)],
                    gamma_var_normal01_simpred_nhhssi1[c(1:3)])),
     ylim = range(c(gamma_var_normal01_simpred_hhssi1[c(1:3)],
                    gamma_var_normal01_simpred_nhhssi1[c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive gamma
plot(gamma_var_normal01_simpred_hhssi1[-c(1:3)],
     gamma_var_normal01_simpred_nhhssi1[-c(1:3)],
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_normal01_simpred_hhssi1[-c(1:3)],
                    gamma_var_normal01_simpred_nhhssi1[-c(1:3)])),
     ylim = range(c(gamma_var_normal01_simpred_hhssi1[-c(1:3)],
                    gamma_var_normal01_simpred_nhhssi1[-c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# overall (infact all are inactive gamma)
plot(gamma_var_normal01_simpred_hhssi1[unique_gamma_idx],
     gamma_var_normal01_simpred_nhhssi1[unique_gamma_idx],
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_normal01_simpred_hhssi1[unique_gamma_idx],
                    gamma_var_normal01_simpred_nhhssi1[unique_gamma_idx])),
     ylim = range(c(gamma_var_normal01_simpred_hhssi1[unique_gamma_idx],
                    gamma_var_normal01_simpred_nhhssi1[unique_gamma_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")

######################################################################################

# lambdasq terms
lambdasq_beta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdasq_beta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update, 2, mean)

lambdasq_gamma_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update, 2, mean)

lambdasq_delta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdasq_delta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_delta_update, 2, mean)

######################################################################################

# lambdasq_beta estimate
# overall
par(mfrow = c(1,1))
plot(lambdasq_beta_normal01_simpred_hhssi1,
     lambdasq_beta_normal01_simpred_nhhssi1,
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1,
                    lambdasq_beta_normal01_simpred_nhhssi1)),
     ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1,
                    lambdasq_beta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active beta
plot(lambdasq_beta_normal01_simpred_hhssi1[c(1:2)],
     lambdasq_beta_normal01_simpred_nhhssi1[c(1:2)],
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1[c(1:2)],
                    lambdasq_beta_normal01_simpred_nhhssi1[c(1:2)])),
     ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1[c(1:2)],
                    lambdasq_beta_normal01_simpred_nhhssi1[c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive beta
plot(lambdasq_beta_normal01_simpred_hhssi1[-c(1:2)],
     lambdasq_beta_normal01_simpred_nhhssi1[-c(1:2)],
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1[-c(1:2)],
                    lambdasq_beta_normal01_simpred_nhhssi1[-c(1:2)])),
     ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1[-c(1:2)],
                    lambdasq_beta_normal01_simpred_nhhssi1[-c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive beta
plot(lambdasq_beta_normal01_simpred_hhssi1[unique_beta_idx],
     lambdasq_beta_normal01_simpred_nhhssi1[unique_beta_idx],
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_normal01_simpred_hhssi1[unique_beta_idx],
                    lambdasq_beta_normal01_simpred_nhhssi1[unique_beta_idx])),
     ylim = range(c(lambdasq_beta_normal01_simpred_hhssi1[unique_beta_idx],
                    lambdasq_beta_normal01_simpred_nhhssi1[unique_beta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# lambdasq_gamma estimate
# overall
par(mfrow = c(1,1))
plot(lambdasq_gamma_normal01_simpred_hhssi1,
     lambdasq_gamma_normal01_simpred_nhhssi1,
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1,
                    lambdasq_gamma_normal01_simpred_nhhssi1)),
     ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1,
                    lambdasq_gamma_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active gamma
plot(lambdasq_gamma_normal01_simpred_hhssi1[c(1:3)],
     lambdasq_gamma_normal01_simpred_nhhssi1[c(1:3)],
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[c(1:3)],
                    lambdasq_gamma_normal01_simpred_nhhssi1[c(1:3)])),
     ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[c(1:3)],
                    lambdasq_gamma_normal01_simpred_nhhssi1[c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive gamma
plot(lambdasq_gamma_normal01_simpred_hhssi1[-c(1:3)],
     lambdasq_gamma_normal01_simpred_nhhssi1[-c(1:3)],
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[-c(1:3)],
                    lambdasq_gamma_normal01_simpred_nhhssi1[-c(1:3)])),
     ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[-c(1:3)],
                    lambdasq_gamma_normal01_simpred_nhhssi1[-c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive gamma
plot(lambdasq_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
     lambdasq_gamma_normal01_simpred_nhhssi1[unique_gamma_idx],
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                    lambdasq_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])),
     ylim = range(c(lambdasq_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                    lambdasq_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# lambdasq_delta estimate
# overall
par(mfrow = c(1,1))
plot(lambdasq_delta_normal01_simpred_hhssi1,
     lambdasq_delta_normal01_simpred_nhhssi1,
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1,
                    lambdasq_delta_normal01_simpred_nhhssi1)),
     ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1,
                    lambdasq_delta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active delta
plot(lambdasq_delta_normal01_simpred_hhssi1[c(1:3, 18:20)],
     lambdasq_delta_normal01_simpred_nhhssi1[c(1:3, 18:20)],
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1[c(1:3, 18:20)],
                    lambdasq_delta_normal01_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1[c(1:3, 18:20)],
                    lambdasq_delta_normal01_simpred_nhhssi1[c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive delta
plot(lambdasq_delta_normal01_simpred_hhssi1[-c(1:3, 18:20)],
     lambdasq_delta_normal01_simpred_nhhssi1[-c(1:3, 18:20)],
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1[-c(1:3, 18:20)],
                    lambdasq_delta_normal01_simpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1[-c(1:3, 18:20)],
                    lambdasq_delta_normal01_simpred_nhhssi1[-c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive delta
plot(lambdasq_delta_normal01_simpred_hhssi1[delta_idx],
     lambdasq_delta_normal01_simpred_nhhssi1[delta_idx],
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_normal01_simpred_hhssi1[delta_idx],
                    lambdasq_delta_normal01_simpred_nhhssi1[delta_idx])),
     ylim = range(c(lambdasq_delta_normal01_simpred_hhssi1[delta_idx],
                    lambdasq_delta_normal01_simpred_nhhssi1[delta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")









######################################################################################

# psi terms
psi_beta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$psi_beta_update, 2, mean)
psi_beta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$psi_beta_update, 2, mean)

psi_gamma_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$psi_gamma_update, 2, mean)
psi_gamma_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$psi_gamma_update, 2, mean)

psi_delta_normal01_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$psi_delta_update, 2, mean)
psi_delta_normal01_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$psi_delta_update, 2, mean)

######################################################################################

# psi_beta estimate
# overall
par(mfrow = c(1,1))
plot(psi_beta_normal01_simpred_hhssi1,
     psi_beta_normal01_simpred_nhhssi1,
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_normal01_simpred_hhssi1,
                    psi_beta_normal01_simpred_nhhssi1)),
     ylim = range(c(psi_beta_normal01_simpred_hhssi1,
                    psi_beta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active beta
plot(psi_beta_normal01_simpred_hhssi1[c(1:2)],
     psi_beta_normal01_simpred_nhhssi1[c(1:2)],
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_normal01_simpred_hhssi1[c(1:2)],
                    psi_beta_normal01_simpred_nhhssi1[c(1:2)])),
     ylim = range(c(psi_beta_normal01_simpred_hhssi1[c(1:2)],
                    psi_beta_normal01_simpred_nhhssi1[c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive beta
plot(psi_beta_normal01_simpred_hhssi1[-c(1:2)],
     psi_beta_normal01_simpred_nhhssi1[-c(1:2)],
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_normal01_simpred_hhssi1[-c(1:2)],
                    psi_beta_normal01_simpred_nhhssi1[-c(1:2)])),
     ylim = range(c(psi_beta_normal01_simpred_hhssi1[-c(1:2)],
                    psi_beta_normal01_simpred_nhhssi1[-c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive beta
plot(psi_beta_normal01_simpred_hhssi1[unique_beta_idx],
     psi_beta_normal01_simpred_nhhssi1[unique_beta_idx],
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_normal01_simpred_hhssi1[unique_beta_idx],
                    psi_beta_normal01_simpred_nhhssi1[unique_beta_idx])),
     ylim = range(c(psi_beta_normal01_simpred_hhssi1[unique_beta_idx],
                    psi_beta_normal01_simpred_nhhssi1[unique_beta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# psi_gamma estimate
# overall
par(mfrow = c(1,1))
plot(psi_gamma_normal01_simpred_hhssi1,
     psi_gamma_normal01_simpred_nhhssi1,
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_normal01_simpred_hhssi1,
                    psi_gamma_normal01_simpred_nhhssi1)),
     ylim = range(c(psi_gamma_normal01_simpred_hhssi1,
                    psi_gamma_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active gamma
plot(psi_gamma_normal01_simpred_hhssi1[c(1:3)],
     psi_gamma_normal01_simpred_nhhssi1[c(1:3)],
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_normal01_simpred_hhssi1[c(1:3)],
                    psi_gamma_normal01_simpred_nhhssi1[c(1:3)])),
     ylim = range(c(psi_gamma_normal01_simpred_hhssi1[c(1:3)],
                    psi_gamma_normal01_simpred_nhhssi1[c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive gamma
plot(psi_gamma_normal01_simpred_hhssi1[-c(1:3)],
     psi_gamma_normal01_simpred_nhhssi1[-c(1:3)],
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_normal01_simpred_hhssi1[-c(1:3)],
                    psi_gamma_normal01_simpred_nhhssi1[-c(1:3)])),
     ylim = range(c(psi_gamma_normal01_simpred_hhssi1[-c(1:3)],
                    psi_gamma_normal01_simpred_nhhssi1[-c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive gamma
plot(psi_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
     psi_gamma_normal01_simpred_nhhssi1[unique_gamma_idx],
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                    psi_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])),
     ylim = range(c(psi_gamma_normal01_simpred_hhssi1[unique_gamma_idx],
                    psi_gamma_normal01_simpred_nhhssi1[unique_gamma_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# psi_delta estimate
# overall
par(mfrow = c(1,1))
plot(psi_delta_normal01_simpred_hhssi1,
     psi_delta_normal01_simpred_nhhssi1,
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_normal01_simpred_hhssi1,
                    psi_delta_normal01_simpred_nhhssi1)),
     ylim = range(c(psi_delta_normal01_simpred_hhssi1,
                    psi_delta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active delta
plot(psi_delta_normal01_simpred_hhssi1[c(1:3, 18:20)],
     psi_delta_normal01_simpred_nhhssi1[c(1:3, 18:20)],
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_normal01_simpred_hhssi1[c(1:3, 18:20)],
                    psi_delta_normal01_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_normal01_simpred_hhssi1[c(1:3, 18:20)],
                    psi_delta_normal01_simpred_nhhssi1[c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive delta
plot(psi_delta_normal01_simpred_hhssi1[-c(1:3, 18:20)],
     psi_delta_normal01_simpred_nhhssi1[-c(1:3, 18:20)],
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_normal01_simpred_hhssi1[-c(1:3, 18:20)],
                    psi_delta_normal01_simpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_normal01_simpred_hhssi1[-c(1:3, 18:20)],
                    psi_delta_normal01_simpred_nhhssi1[-c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive delta
plot(psi_delta_normal01_simpred_hhssi1[delta_idx],
     psi_delta_normal01_simpred_nhhssi1[delta_idx],
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_normal01_simpred_hhssi1[delta_idx],
                    psi_delta_normal01_simpred_nhhssi1[delta_idx])),
     ylim = range(c(psi_delta_normal01_simpred_hhssi1[delta_idx],
                    psi_delta_normal01_simpred_nhhssi1[delta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")


######################################################################################

# tausq terms
tausq_beta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_beta_update
tausq_beta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_beta_update

tausq_gamma_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_gamma_update
tausq_gamma_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_gamma_update

tausq_delta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$tausq_delta_update
tausq_delta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$tausq_delta_update

######################################################################################

# tausq_beta estimate
plot(tausq_beta_normal01_simpred_hhssi1,
     tausq_beta_normal01_simpred_nhhssi1,
     xlab = "HHSSI tausq beta",
     ylab = "NHHSSI tausq beta",
     xlim = range(c(tausq_beta_normal01_simpred_hhssi1,
                    tausq_beta_normal01_simpred_nhhssi1)),
     ylim = range(c(tausq_beta_normal01_simpred_hhssi1,
                    tausq_beta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# tausq_gamma estimate
plot(tausq_gamma_normal01_simpred_hhssi1,
     tausq_gamma_normal01_simpred_nhhssi1,
     xlab = "HHSSI tausq gamma",
     ylab = "NHHSSI tausq gamma",
     xlim = range(c(tausq_gamma_normal01_simpred_hhssi1,
                    tausq_gamma_normal01_simpred_nhhssi1)),
     ylim = range(c(tausq_gamma_normal01_simpred_hhssi1,
                    tausq_gamma_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")



# tausq_delta estimate
plot(tausq_delta_normal01_simpred_hhssi1,
     tausq_delta_normal01_simpred_nhhssi1,
     xlab = "HHSSI tausq delta",
     ylab = "NHHSSI tausq delta",
     xlim = range(c(tausq_delta_normal01_simpred_hhssi1,
                    tausq_delta_normal01_simpred_nhhssi1)),
     ylim = range(c(tausq_delta_normal01_simpred_hhssi1,
                    tausq_delta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")





######################################################################################

# xi terms
xi_beta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$xi_beta_update
xi_beta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$xi_beta_update

xi_gamma_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$xi_gamma_update
xi_gamma_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$xi_gamma_update

xi_delta_normal01_simpred_hhssi1 <- res_normal01_simpred_hhssi1$xi_delta_update
xi_delta_normal01_simpred_nhhssi1 <- res_normal01_simpred_nhhssi1$xi_delta_update

######################################################################################

# xi_beta estimate
plot(xi_beta_normal01_simpred_hhssi1,
     xi_beta_normal01_simpred_nhhssi1,
     xlab = "HHSSI xi beta",
     ylab = "NHHSSI xi beta",
     xlim = range(c(xi_beta_normal01_simpred_hhssi1,
                    xi_beta_normal01_simpred_nhhssi1)),
     ylim = range(c(xi_beta_normal01_simpred_hhssi1,
                    xi_beta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# xi_gamma estimate
plot(xi_gamma_normal01_simpred_hhssi1,
     xi_gamma_normal01_simpred_nhhssi1,
     xlab = "HHSSI xi gamma",
     ylab = "NHHSSI xi gamma",
     xlim = range(c(xi_gamma_normal01_simpred_hhssi1,
                    xi_gamma_normal01_simpred_nhhssi1)),
     ylim = range(c(xi_gamma_normal01_simpred_hhssi1,
                    xi_gamma_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")



# xi_delta estimate
plot(xi_delta_normal01_simpred_hhssi1,
     xi_delta_normal01_simpred_nhhssi1,
     xlab = "HHSSI xi delta",
     ylab = "NHHSSI xi delta",
     xlim = range(c(xi_delta_normal01_simpred_hhssi1,
                    xi_delta_normal01_simpred_nhhssi1)),
     ylim = range(c(xi_delta_normal01_simpred_hhssi1,
                    xi_delta_normal01_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


#############################################################################################
# prior var og=f beta, gamma
beta_var_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_beta_update *
                                        res_normal01_simpred_hhssi1$tausq_beta_update, 2, mean)
beta_var_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_beta_update *
                                        res_normal01_simpred_nhhssi1$tausq_beta_update, 2, mean)


gamma_var_mean_simpred_hhssi1 <- apply(res_normal01_simpred_hhssi1$lambdasq_gamma_update *
                                        res_normal01_simpred_hhssi1$tausq_gamma_update, 2, mean)
gamma_var_mean_simpred_nhhssi1 <- apply(res_normal01_simpred_nhhssi1$lambdasq_gamma_update *
                                        res_normal01_simpred_nhhssi1$tausq_gamma_update, 2, mean)


point_col_beta <- ifelse(seq_along(beta_var_mean_simpred_hhssi1) %in% 1:2, "blue", "red")
plot(
  log(beta_var_mean_simpred_hhssi1),
  log(beta_var_mean_simpred_nhhssi1),
  pch = 19,
  col = point_col_beta,
  main= "log transformed",
  xlab = "HHSSI beta prior var",
  ylab = "NHHSSI beta prior var",
  xlim=range(c(log(beta_var_mean_simpred_hhssi1), log(beta_var_mean_simpred_nhhssi1))),
  ylim=range(c(log(beta_var_mean_simpred_hhssi1), log(beta_var_mean_simpred_nhhssi1)))
)
# Legend
legend("topleft",
       legend = c(
         "active beta",
         "inactive"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(0, 1, lty = 2, col="red")


plot(
  beta_var_mean_simpred_hhssi1[1:2],
  beta_var_mean_simpred_nhhssi1[1:2],
  pch = 19,
  col = "blue",
  xlab = "HHSSI beta prior var",
  ylab = "NHHSSI beta prior var",
  xlim=range(c(beta_var_mean_simpred_hhssi1[1:2], beta_var_mean_simpred_nhhssi1[1:2])),
  ylim=range(c(beta_var_mean_simpred_hhssi1[1:2], beta_var_mean_simpred_nhhssi1[1:2]))
)
abline(0, 1, lty = 2, col="red")


plot(
  beta_var_mean_simpred_hhssi1[-(1:2)],
  beta_var_mean_simpred_nhhssi1[-(1:2)],
  pch = 19,
  col = "red",
  xlab = "HHSSI beta prior var",
  ylab = "NHHSSI beta prior var",
  xlim=range(c(beta_var_mean_simpred_hhssi1[-(1:2)], beta_var_mean_simpred_nhhssi1[-(1:2)])),
  ylim=range(c(beta_var_mean_simpred_hhssi1[-(1:2)], beta_var_mean_simpred_nhhssi1[-(1:2)]))
)
abline(0, 1, lty = 2, col="red")






point_col_gamma <- ifelse(seq_along(gamma_var_mean_simpred_hhssi1) %in% 1:3, "blue", "red")
plot(
  log(gamma_var_mean_simpred_hhssi1),
  log(gamma_var_mean_simpred_nhhssi1),
  pch = 19,
  col = point_col_gamma,
  main= "log transformed",
  xlab = "HHSSI gamma prior var",
  ylab = "NHHSSI gamma prior var",
  xlim=range(c(log(gamma_var_mean_simpred_hhssi1), log(gamma_var_mean_simpred_nhhssi1))),
  ylim=range(c(log(gamma_var_mean_simpred_hhssi1), log(gamma_var_mean_simpred_nhhssi1)))
)
# Legend
legend("topleft",
       legend = c(
         "active gamma",
         "inactive"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(0, 1, lty = 2, col="red")


plot(
  gamma_var_mean_simpred_hhssi1[1:3],
  gamma_var_mean_simpred_nhhssi1[1:3],
  pch = 19,
  col = "blue",
  xlab = "HHSSI gamma prior var",
  ylab = "NHHSSI gamma prior var",
  xlim=range(c(gamma_var_mean_simpred_hhssi1[1:3], gamma_var_mean_simpred_nhhssi1[1:3])),
  ylim=range(c(gamma_var_mean_simpred_hhssi1[1:3], gamma_var_mean_simpred_nhhssi1[1:3]))
)
abline(0, 1, lty = 2, col="red")


plot(
  gamma_var_mean_simpred_hhssi1[-(1:3)],
  gamma_var_mean_simpred_nhhssi1[-(1:3)],
  pch = 19,
  col = "red",
  xlab = "HHSSI gamma prior var",
  ylab = "NHHSSI gamma prior var",
  xlim=range(c(gamma_var_mean_simpred_hhssi1[-(1:3)], gamma_var_mean_simpred_nhhssi1[-(1:3)])),
  ylim=range(c(gamma_var_mean_simpred_hhssi1[-(1:3)], gamma_var_mean_simpred_nhhssi1[-(1:3)]))
)
abline(0, 1, lty = 2, col="red")

##############################################################################################

# ratio of prior variance estimate from NHHSSI and HHSSI vs lambdasq_beta*lambdasq_gamma

beta_gamma_lam_mean_simpred_hhssi1
delta_tau_simpred_hhssi1 <- (delta_lam_simpred_hhssi1 *
                                tausq_delta_simpred_hhssi1)
delta_tau_mean_simpred_hhssi1 <- apply(delta_tau_simpred_hhssi1, 2, mean)


plot(
  log(beta_gamma_lam_mean_simpred_hhssi1),
  log(delta_tau_mean_simpred_hhssi1),
  pch = 19,
  main= "log transformed",
  ylab = expression(hat(lambda)[delta]^2 * hat(tau)[delta]^2),
  xlab = expression(hat(lambda)[beta]^2 * hat(lambda)[gamma]^2),
  xlim=range(c(log(beta_gamma_lam_mean_simpred_hhssi1), log(delta_tau_mean_simpred_hhssi1))),
  ylim=range(c(log(beta_gamma_lam_mean_simpred_hhssi1), log(delta_tau_mean_simpred_hhssi1)))
)
abline(0, 1, lty = 2, col="red")


# Identify whether NHHSSI variance is greater than HHSSI variance
idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
# Colors
point_col <- ifelse(idx, "blue", "red")
# Common limits
x_vals <- log(beta_gamma_lam_mean_simpred_hhssi1)
y_vals <- log(delta_tau_mean_simpred_hhssi1)
plot(x_vals, y_vals,
  pch = 19,
  col = point_col,
  main = "HHSSI prior variance structure",
  ylab = expression(log(hat(lambda)[delta]^2 * hat(tau)[delta]^2)),
  xlab = expression(log(hat(lambda)[beta]^2 * hat(lambda)[gamma]^2)),
  xlim = range(c(x_vals, y_vals)),
  ylim = range(c(x_vals, y_vals))
)
# 45-degree reference line
abline(0, 1, lty = 2, col = "black")
# Legend
legend("topright",
  legend = c(
    "NHHSSI variance > HHSSI variance",
    "NHHSSI variance <= HHSSI variance"),
  col = c("blue", "red"),
  pch = 19,
  bty = "n"
)



#################################################

# Beta-gamma shrinkage term
bg <- beta_gamma_lam_mean_simpred_hhssi1
# Log transform
log_bg <- log(bg)
# Rescale to a reasonable point-size range
cex_vals <- 0.5 +
  2.5 * (log_bg - min(log_bg)) /
  (max(log_bg) - min(log_bg))
# Prior variances
x_vals <- log(delta_var_mean_simpred_hhssi1)
y_vals <- log(delta_var_mean_simpred_nhhssi1)
idx <- delta_var_mean_simpred_nhhssi1 > delta_var_mean_simpred_hhssi1
# Colors
point_col <- ifelse(idx, "blue", "red")
# Plot
plot(x_vals, y_vals,
  pch = 21,
  cex = cex_vals,
  bg = point_col,
  col = "black",
  main = "Log-transformed prior variance",
  ylab = "log(NHHSSI prior variance)",
  xlab = "log(HHSSI prior variance)",
  xlim = range(c(x_vals, y_vals)),
  ylim = range(c(x_vals, y_vals))
)
# 45-degree reference line
abline(
  0, 1,
  lty = 2,
  col = "black"
)
# Legend for variance comparison
legend("topleft",
  legend = c("NHHSSI variance > HHSSI variance",
    "NHHSSI variance <= HHSSI variance"
  ),
  pt.bg = c("blue", "red"),
  col = "black",
  pch = 21,
  bty = "n"
)





# Colors
point_col <- ifelse(idx, "blue", "red")
# Common limits
x_vals <- log(delta_tau_mean_simpred_hhssi1)
y_vals <- log(delta_var_mean_simpred_nhhssi1)
plot(x_vals, y_vals,
     pch = 19,
     col = point_col,
     main = "HHSSI prior variance structure",
     ylab = expression(NHHSSI_log(hat(lambda)[delta]^2 * hat(tau)[delta]^2)),
     xlab = expression(HHSSI_log(hat(lambda)[delta]^2 * hat(tau)[delta]^2)),
     xlim = range(c(x_vals, y_vals)),
     ylim = range(c(x_vals, y_vals))
)
# 45-degree reference line
abline(0, 1, lty = 2, col = "black")
# Legend
legend("topright",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n"
)


dat <- round(cbind(log(delta_tau_mean_simpred_hhssi1),
            log(delta_var_mean_simpred_hhssi1),
            log(delta_var_mean_simpred_nhhssi1),
            log(beta_gamma_lam_mean_simpred_hhssi1)), 5)
dat_sorted <- dat[order(dat[, 3], decreasing = TRUE), ]

dat_sorted




##########################################################

active_idx <- c(1:3, 18:20)

J <- 44
M <- 17

active_info <- data.frame(
  index = active_idx,
  j = ((active_idx - 1) %/% M) + 1,
  m = ((active_idx - 1) %% M) + 1,
  hhssi_var = delta_var_mean_simpred_hhssi1[active_idx],
  nhhssi_var = delta_var_mean_simpred_nhhssi1[active_idx],
  bg = beta_gamma_lam_mean_simpred_hhssi1[active_idx]
)
active_info



dim(delta_est_normal01_simpred_hhssi1)
delta_est_normal01_simpred_nhhssi1
delta_true <- sim_data_normal01_simpred1$delta_true

rmse_hhssi <- sqrt(rowMeans((delta_est_normal01_simpred_hhssi1 - delta_true)^2))
rmse_nhhssi <- sqrt(rowMeans((delta_est_normal01_simpred_nhhssi1 - delta_true)^2 ))
active_rmse <- data.frame(
  index = active_idx,
  j = ((active_idx - 1) %/% M) + 1,
  m = ((active_idx - 1) %% M) + 1,
  HHSSI = rmse_hhssi[active_idx],
  NHHSSI = rmse_nhhssi[active_idx]
)
active_rmse




bias_hhssi <- rowMeans(delta_est_normal01_simpred_hhssi1 - delta_true)
bias_nhhssi <- rowMeans(delta_est_normal01_simpred_nhhssi1 - delta_true)
active_bias <- data.frame(
  index = active_idx,
  j = ((active_idx - 1) %/% M) + 1,
  m = ((active_idx - 1) %% M) + 1,
  HHSSI = bias_hhssi[active_idx],
  NHHSSI = bias_nhhssi[active_idx]
)
active_bias




active_absbias <- data.frame(
  index = active_idx,
  j = ((active_idx - 1) %/% M) + 1,
  m = ((active_idx - 1) %% M) + 1,
  HHSSI = abs(bias_hhssi[active_idx]),
  NHHSSI = abs(bias_nhhssi[active_idx])
)
active_absbias


active_results <- data.frame(
  index = active_idx,
  j = ((active_idx - 1) %/% M) + 1,
  m = ((active_idx - 1) %% M) + 1,

  log_bg = log(beta_gamma_lam_mean_simpred_hhssi1[active_idx]),

  log_var_hhssi =
    log(delta_var_mean_simpred_hhssi1[active_idx]),

  log_var_nhhssi =
    log(delta_var_mean_simpred_nhhssi1[active_idx]),

  RMSE_HHSSI = rmse_hhssi[active_idx],
  RMSE_NHHSSI = rmse_nhhssi[active_idx],

  Bias_HHSSI = bias_hhssi[active_idx],
  Bias_NHHSSI = bias_nhhssi[active_idx],

  AbsBias_HHSSI = abs(bias_hhssi[active_idx]),
  AbsBias_NHHSSI = abs(bias_nhhssi[active_idx])
)

round(active_results, 4)




plot(active_rmse$HHSSI, active_rmse$NHHSSI,
  pch = 19,
  xlab = "HHSSI RMSE",
  ylab = "NHHSSI RMSE",
  main = "RMSE for active interaction coefficients",
  xlim = range(c(active_rmse$HHSSI, active_rmse$NHHSSI)),
  ylim = range(c(active_rmse$HHSSI, active_rmse$NHHSSI))
)
abline(0, 1, lty = 2)
text(active_rmse$HHSSI, active_rmse$NHHSSI,
  labels = paste0("(", active_rmse$j, ",", active_rmse$m, ")"),
  pos = 4,
  cex = 0.8)


plot(active_bias$HHSSI, active_bias$NHHSSI,
  pch = 19,
  xlab = "HHSSI bias",
  ylab = "NHHSSI bias",
  main = "Bias for active interaction coefficients")
abline(h = 0, v = 0, lty = 2)
abline(0, 1, lty = 2)
text(active_bias$HHSSI, active_bias$NHHSSI,
  labels = paste0("(", active_bias$j, ",", active_bias$m, ")"),
  pos = 4,
  cex = 0.8)




point_col <- ifelse(idx, "blue", "red")
plot(log(delta_var_mean_simpred_hhssi1), log(rmse_hhssi),
    pch = 19,
    cex = 0.7,
    col = point_col,
    xlab = "log(HHSSI prior variance)",
    ylab = "log(HHSSI RMSE)",
    main = "HHSSI prior variance vs RMSE")
points(log(delta_var_mean_simpred_hhssi1[active_idx]),
    log(rmse_hhssi[active_idx]),
    pch = 19,
    cex = 1.5,
    col="red")
text(log(delta_var_mean_simpred_hhssi1[active_idx]),
    log(rmse_hhssi[active_idx]),
    labels = active_idx,
    pos = 4,
    cex = 0.8)
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(0,1)




plot(log(delta_var_mean_simpred_nhhssi1),
  log(rmse_nhhssi),
  pch = 19,
  cex = 0.7,
  col = point_col,
  xlab = "log(NHHSSI prior variance)",
  ylab = "log(NHHSSI RMSE)",
  main = "NHHSSI prior variance vs RMSE")
points(log(delta_var_mean_simpred_nhhssi1[active_idx]),
  log(rmse_nhhssi[active_idx]),
  pch = 19,
  cex = 1.5,
  col="red")
text(log(delta_var_mean_simpred_nhhssi1[active_idx]),
  log(rmse_nhhssi[active_idx]),
  labels = active_idx,
  pos = 4,
  cex = 0.8)
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(0,1)




plot(log(delta_var_mean_simpred_hhssi1),
     bias_hhssi,
     pch = 19,
     cex = 0.7,
     col = point_col,
     xlab = "log(HHSSI prior variance)",
     ylab = "HHSSI bias",
     main = "HHSSI prior variance vs bias")
points(log(delta_var_mean_simpred_hhssi1[active_idx]),
       bias_hhssi[active_idx],
       pch = 19,
       cex = 1.5,
       col="red")
text(log(delta_var_mean_simpred_hhssi1[active_idx]),
     bias_hhssi[active_idx],
     labels = active_idx,
     pos = 4,
     cex = 0.8)
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(h=0)





plot(log(delta_var_mean_simpred_nhhssi1),
     bias_nhhssi,
     pch = 19,
     cex = 0.7,
     col = point_col,
     xlab = "log(NHHSSI prior variance)",
     ylab = "NHHSSI bias",
     main = "NHHSSI prior variance vs bias")
points(log(delta_var_mean_simpred_nhhssi1[active_idx]),
       bias_nhhssi[active_idx],
       pch = 19,
       cex = 1.5,
       col="red")
text(log(delta_var_mean_simpred_nhhssi1[active_idx]),
     bias_nhhssi[active_idx],
     labels = active_idx,
     pos = 4,
     cex = 0.8)
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n")
abline(h=0)




point_col <- ifelse(idx, "blue", "red")
plot(log(rmse_hhssi), log(rmse_nhhssi),
     pch = 19,
     cex = 0.7,
     col = point_col,
     xlab = "log(HHSSI RMSE)",
     ylab = "log(NHHSSI RMSE)",
     main = "RMSE of HHSSI vs NHHSSI")
points(log(rmse_hhssi)[active_idx], log(rmse_nhhssi)[active_idx],
       pch = 19,
       cex = 1.5,
       col="red")
text(log(rmse_hhssi[active_idx]),
     log(rmse_nhhssi[active_idx]),
     labels = active_idx,
     pos = 4,
     cex = 0.8)
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n"
)
abline(0,1)





point_col <- ifelse(idx, "blue", "red")
plot(bias_hhssi, bias_nhhssi,
     pch = 19,
     cex = 0.7,
     col = point_col,
     xlab = "HHSSI bias",
     ylab = "NHHSSI bias",
     main = "Bias of HHSSI vs NHHSSI")
points(bias_hhssi[active_idx], bias_nhhssi[active_idx],
       pch = 19,
       cex = 1.5,
       col="red")
text(bias_hhssi[active_idx],
     bias_nhhssi[active_idx],
     labels = active_idx,
     pos = 4,
     cex = 0.8)
legend("topleft",
       legend = c(
         "NHHSSI variance > HHSSI variance",
         "NHHSSI variance <= HHSSI variance"),
       col = c("blue", "red"),
       pch = 19,
       bty = "n"
)
abline(h=0)



##########################################################################

delta_true <- sim_data_normal01_simpred1$delta_true

rmse_hhssi1 <- sqrt(rowMeans((delta_est_normal01_simpred_hhssi1 - delta_true)^2))
rmse_nhhssi1 <- sqrt(rowMeans((delta_est_normal01_simpred_nhhssi1 - delta_true)^2 ))

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$RMSE_HHSSI <- (rmse_hhssi1)
df$RMSE_NHHSSI <- (rmse_nhhssi1)
df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
  modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "RMSE") +

  labs(
    title = "RMSE of Modifier–Exposure Interaction Estimates (with simulated predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))


##########################################################################

delta_est_normal01_realpred_hhssi1 <- apply(res_normal01_realpred_hhssi1$theta_update[, 63:810, ], c(2,3), mean)
delta_est_normal01_realpred_nhhssi1 <- apply(res_normal01_realpred_nhhssi1$theta_update[, 63:810, ], c(2,3), mean)

delta_true <- sim_data_normal01_realpred1$delta_true

rmse_hhssi1 <- sqrt(rowMeans((delta_est_normal01_realpred_hhssi1 - delta_true)^2))
rmse_nhhssi1 <- sqrt(rowMeans((delta_est_normal01_realpred_nhhssi1 - delta_true)^2 ))

df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$RMSE_HHSSI <- rmse_hhssi1
df$RMSE_NHHSSI <- rmse_nhhssi1
df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            RMSE = df$RMSE_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, RMSE = df$RMSE_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = RMSE)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "RMSE") +

  labs(
    title = "RMSE of Modifier–Exposure Interaction Estimates (with real predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))


######################################################

# prior var heatmap

# with real predictors
delta_var_mean_simpred_hhssi1
delta_var_mean_simpred_nhhssi1

library(ggplot2)
df <- expand.grid(exposure = 1:M, modifier = 1:J)
df$prior_var_HHSSI <- log(delta_var_mean_simpred_hhssi1)
df$prior_var_NHHSSI <- log(delta_var_mean_simpred_nhhssi1)
df_long <- rbind(data.frame(exposure = df$exposure, modifier = df$modifier,
                            Prior_var = df$prior_var_HHSSI,
                            Method = "HHSSI"),
                 data.frame(exposure = df$exposure, modifier = df$modifier, Prior_var = df$prior_var_NHHSSI,
                            Method = "NHHSSI"))
# Identify active interactions
df_long$Active <- with(df_long,
                       modifier %in% 1:2 & exposure %in% 1:3)
ggplot(df_long, aes(x = exposure, y = modifier, fill = Prior_var)) +
  geom_tile() +
  # Outline active interaction tiles
  geom_tile(
    data = subset(df_long, Active),
    fill = NA,
    color = "black",
    linewidth = 0.8) +
  facet_wrap(~Method) +

  scale_x_continuous(breaks = 1:M) +
  scale_y_continuous(breaks = seq(1, J)) +

  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Prior variance",
    limits = range(c(log(delta_var_mean_simpred_hhssi1),
                     log(delta_var_mean_simpred_nhhssi1)))) +

  labs(
    title = "RMSE of Modifier–Exposure Interaction Estimates (with 100 datasets of real predictors)",
    subtitle = "(Black borders indicate active interactions)",
    x = "Exposure index",
    y = "Modifier index") +

  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))






