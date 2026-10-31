
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



############################
# β, γ, δ, φ ∼ N(0, 1), and α ∼ Unif (−3.5, −2.5), all predictors simulated
############################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
set.seed(23444)
sim_data_normal01_simpred <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)



Y <- sim_data_normal01_simpred$Y
x <- sim_data_normal01_simpred$x
z <- sim_data_normal01_simpred$z
u <- sim_data_normal01_simpred$u
d <- sim_data_normal01_simpred$d
W <- sim_data_normal01_simpred$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(456557)
res_normal01_simpred_hhssi <- fit_HHSSI_MVN_cov_modify3(
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


# save(res_normal01_simpred_hhssi, file="practice_data/res_normal01_simpred_hhssi.rda")
# load(file="practice_data/res_normal01_simpred_hhssi.rda")


start_time <- Sys.time()
set.seed(456557)
res_normal01_simpred_nhhssi <- fit_NHHSSI_MVN_cov_modify3(
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

# save(res_normal01_simpred_nhhssi, file="practice_data/res_normal01_simpred_nhhssi.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi.rda")

############################
# β, γ, δ, φ ∼ N(0, 1), and α ∼ Unif (−3.5, −2.5), all predictors real
############################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
set.seed(23444)
sim_data_normal01_realpred <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D,
                              Sigma=Sigma)



Y <- sim_data_normal01_realpred$Y
x <- sim_data_normal01_realpred$x
z <- sim_data_normal01_realpred$z
u <- sim_data_normal01_realpred$u
d <- sim_data_normal01_realpred$d
W <- sim_data_normal01_realpred$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(456557)
res_normal01_realpred_hhssi <- fit_HHSSI_MVN_cov_modify3(
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


# save(res_normal01_realpred_hhssi, file="practice_data/res_normal01_realpred_hhssi.rda")
# load(file="practice_data/res_normal01_realpred_hhssi.rda")


start_time <- Sys.time()
set.seed(456557)
res_normal01_realpred_nhhssi <- fit_NHHSSI_MVN_cov_modify3(
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


# save(res_normal01_realpred_nhhssi, file="practice_data/res_normal01_realpred_nhhssi.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi.rda")









########################################################
# α, β, γ, δ, φ ∼ unif(2,3), all predictors simulated
#######################################################

set.seed(23444)
sim_data_unif23_simpred <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)



Y <- sim_data_unif23_simpred$Y
x <- sim_data_unif23_simpred$x
z <- sim_data_unif23_simpred$z
u <- sim_data_unif23_simpred$u
d <- sim_data_unif23_simpred$d
W <- sim_data_unif23_simpred$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)


start_time <- Sys.time()
set.seed(456557)
res_unif23_simpred_hhssi <- fit_HHSSI_MVN_cov_modify3(
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


# save(res_unif23_simpred_hhssi, file="practice_data/res_unif23_simpred_hhssi.rda")
# load(file="practice_data/res_unif23_simpred_hhssi.rda")


start_time <- Sys.time()
set.seed(456557)
res_unif23_simpred_nhhssi <- fit_NHHSSI_MVN_cov_modify3(
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


# save(res_unif23_simpred_nhhssi, file="practice_data/res_unif23_simpred_nhhssi.rda")
# load(file="practice_data/res_unif23_simpred_nhhssi.rda")


############################
# α, β, γ, δ, φ ∼ unif(2,3), all predictors real
############################

set.seed(23444)
sim_data_unif23_realpred <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                O=O, x=X, z=Z, d=D,
                                                Sigma=Sigma)



Y <- sim_data_unif23_realpred$Y
x <- sim_data_unif23_realpred$x
z <- sim_data_unif23_realpred$z
u <- sim_data_unif23_realpred$u
d <- sim_data_unif23_realpred$d
W <- sim_data_unif23_realpred$W

WTW <- t(W) %*% W

n_all_par <- ncol(W)
nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)



start_time <- Sys.time()
set.seed(456557)
res_unif23_realpred_hhssi <- fit_HHSSI_MVN_cov_modify3(
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


# save(res_unif23_realpred_hhssi, file="practice_data/res_unif23_realpred_hhssi.rda")
# load(file="practice_data/res_unif23_realpred_hhssi.rda")



start_time <- Sys.time()
set.seed(456557)
res_unif23_realpred_nhhssi <- fit_NHHSSI_MVN_cov_modify3(
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


# save(res_unif23_realpred_nhhssi, file="practice_data/res_unif23_realpred_nhhssi.rda")
# load(file="practice_data/res_unif23_realpred_nhhssi.rda")











###############################################################################
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors simulated
###############################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
set.seed(23444)
sim_data_normal01_simpred <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)

Y <- sim_data_normal01_simpred$Y
x <- sim_data_normal01_simpred$x
z <- sim_data_normal01_simpred$z
u <- sim_data_normal01_simpred$u
d <- sim_data_normal01_simpred$d
W <- sim_data_normal01_simpred$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_simpred_hhssi.rda")
# load(file="practice_data/res_normal01_simpred_nhhssi.rda")

sim_data_normal01_simpred$delta_true



dim(res_normal01_simpred_hhssi$theta_update)

# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
   matplot(res_normal01_simpred_hhssi$theta_update[, 1, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = expression(alpha),
    main = paste("Outcome", k)
  )

  matlines(res_normal01_simpred_nhhssi$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_simpred_hhssi$theta_update[, 1, ])
range(res_normal01_simpred_nhhssi$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
hhssi_alpha <- res_normal01_simpred_hhssi$theta_update[, 1, k]
nhhssi_alpha <- res_normal01_simpred_nhhssi$theta_update[, 1, k]

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
sim_data_normal01_simpred$alpha_true
res_normal01_simpred_hhssi$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_simpred_hhssi$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_simpred_nhhssi$theta_update[, 1, ], 2, mean)
# Bias for each outcome
alpha_bias_hhssi <- (alpha_est_hhssi - sim_data_normal01_simpred$alpha_true)
alpha_bias_nhhssi <- (alpha_est_nhhssi - sim_data_normal01_simpred$alpha_true)


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

plot(res_normal01_simpred_hhssi$theta_update[, 2:3, 1],
     res_normal01_simpred_nhhssi$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, 2:3, ],
                    res_normal01_simpred_hhssi$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, 2:3, 2],
       res_normal01_simpred_nhhssi$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, 2:3, 3],
       res_normal01_simpred_nhhssi$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, 2:3, 4],
       res_normal01_simpred_nhhssi$theta_update[, 2:3, 4],
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




range(res_normal01_simpred_hhssi$theta_update[, 2:3, 2])
range(res_normal01_simpred_nhhssi$theta_update[, 2:3, 2])

range(res_normal01_simpred_hhssi$theta_update[, 2:3, 3])
range(res_normal01_simpred_nhhssi$theta_update[, 2:3, 3])

range(res_normal01_simpred_hhssi$theta_update[, 2:3, 4])
range(res_normal01_simpred_nhhssi$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi$theta_update[, 2:3, k],
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
range(res_normal01_simpred_hhssi$theta_update[, 4:45, 1])
range(res_normal01_simpred_nhhssi$theta_update[, 4:45, 1])

plot(res_normal01_simpred_hhssi$theta_update[, 4:45, 1],
     res_normal01_simpred_nhhssi$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, 4:45, ],
                    res_normal01_simpred_hhssi$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, 4:45, 2],
       res_normal01_simpred_nhhssi$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, 4:45, 3],
       res_normal01_simpred_nhhssi$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, 4:45, 4],
       res_normal01_simpred_nhhssi$theta_update[, 4:45, 4],
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



range(res_normal01_simpred_hhssi$theta_update[, 4:45, 2])
range(res_normal01_simpred_nhhssi$theta_update[, 4:45, 2])

range(res_normal01_simpred_hhssi$theta_update[, 4:45, 3])
range(res_normal01_simpred_nhhssi$theta_update[, 4:45, 3])

range(res_normal01_simpred_hhssi$theta_update[, 4:45, 4])
range(res_normal01_simpred_nhhssi$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_simpred_nhhssi$theta_update[, 4:45, k],
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
dim(sim_data_normal01_simpred$beta_true)
dim(res_normal01_simpred_hhssi$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_simpred_hhssi$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_simpred_nhhssi$theta_update[, 2:45, ], c(2,3), mean)
# Bias for each outcome
beta_bias_hhssi <- (beta_est_hhssi - sim_data_normal01_simpred$beta_true)
beta_bias_nhhssi <- (beta_est_nhhssi - sim_data_normal01_simpred$beta_true)

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



beta_var_hhssi <- (res_normal01_simpred_hhssi$lambdasq_beta_update *
                     res_normal01_simpred_hhssi$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_simpred_nhhssi$lambdasq_beta_update *
                      res_normal01_simpred_nhhssi$tausq_beta_update)

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
range(res_normal01_simpred_hhssi$theta_update[, 46:48, 1])
range(res_normal01_simpred_nhhssi$theta_update[, 46:48, 1])


plot(res_normal01_simpred_hhssi$theta_update[, 46:48, 1],
     res_normal01_simpred_nhhssi$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, 46:48, ],
                    res_normal01_simpred_hhssi$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, 46:48, 2],
       res_normal01_simpred_nhhssi$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, 46:48, 3],
       res_normal01_simpred_nhhssi$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, 46:48, 4],
       res_normal01_simpred_nhhssi$theta_update[, 46:48, 4],
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

range(res_normal01_simpred_hhssi$theta_update[, 46:48, 2])
range(res_normal01_simpred_nhhssi$theta_update[, 46:48, 2])

range(res_normal01_simpred_hhssi$theta_update[, 46:48, 3])
range(res_normal01_simpred_nhhssi$theta_update[, 46:48, 3])

range(res_normal01_simpred_hhssi$theta_update[, 46:48, 4])
range(res_normal01_simpred_nhhssi$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi$theta_update[, 46:48, k],
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
range(res_normal01_simpred_hhssi$theta_update[, 49:62, 1])
range(res_normal01_simpred_nhhssi$theta_update[, 49:62, 1])


plot(res_normal01_simpred_hhssi$theta_update[, 49:62, 1],
     res_normal01_simpred_nhhssi$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, 49:62, ],
                    res_normal01_simpred_hhssi$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, 49:62, 2],
       res_normal01_simpred_nhhssi$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, 49:62, 3],
       res_normal01_simpred_nhhssi$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, 49:62, 4],
       res_normal01_simpred_nhhssi$theta_update[, 49:62, 4],
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


range(res_normal01_simpred_hhssi$theta_update[, 49:62, 2])
range(res_normal01_simpred_nhhssi$theta_update[, 49:62, 2])

range(res_normal01_simpred_hhssi$theta_update[, 49:62, 3])
range(res_normal01_simpred_nhhssi$theta_update[, 49:62, 3])

range(res_normal01_simpred_hhssi$theta_update[, 49:62, 4])
range(res_normal01_simpred_nhhssi$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi$theta_update[, 49:62, k],
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
dim(sim_data_normal01_simpred$gamma_true)
dim(res_normal01_simpred_hhssi$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_simpred_hhssi$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_simpred_nhhssi$theta_update[, 46:62, ], c(2,3), mean)
# Bias for each outcome
gamma_bias_hhssi <- (gamma_est_hhssi - sim_data_normal01_simpred$gamma_true)
gamma_bias_nhhssi <- (gamma_est_nhhssi - sim_data_normal01_simpred$gamma_true)

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



gamma_var_hhssi <- (res_normal01_simpred_hhssi$lambdasq_gamma_update *
                      res_normal01_simpred_hhssi$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_simpred_nhhssi$lambdasq_gamma_update *
                       res_normal01_simpred_nhhssi$tausq_gamma_update)

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
range(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 1],
     res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), ],
                    res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 2],
       res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 3],
       res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 4],
       res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 4],
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


range(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_simpred_hhssi$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_simpred_nhhssi$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi$theta_update[, c(63:65,80:82), k],
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

range(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 1],
     res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), ],
                    res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 2],
       res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 3],
       res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 4],
       res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 4],
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


range(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_simpred_hhssi$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_simpred_nhhssi$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi$theta_update[, c(66:79,83:810), k],
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
dim(sim_data_normal01_simpred$delta_true)
dim(res_normal01_simpred_hhssi$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_simpred_hhssi$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_simpred_nhhssi$theta_update[, 63:810, ], c(2,3), mean)
# Bias for each outcome
delta_bias_hhssi <- (delta_est_hhssi - sim_data_normal01_simpred$delta_true)
delta_bias_nhhssi <- (delta_est_nhhssi - sim_data_normal01_simpred$delta_true)

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
range(res_normal01_simpred_hhssi$theta_update[, 811:812, 1])
range(res_normal01_simpred_nhhssi$theta_update[, 811:812, 1])

plot(res_normal01_simpred_hhssi$theta_update[, 811:812, 1],
     res_normal01_simpred_nhhssi$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, 811:812, ],
                    res_normal01_simpred_hhssi$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, 811:812, 2],
       res_normal01_simpred_nhhssi$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, 811:812, 3],
       res_normal01_simpred_nhhssi$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, 811:812, 4],
       res_normal01_simpred_nhhssi$theta_update[, 811:812, 4],
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


range(res_normal01_simpred_hhssi$theta_update[, 811:812, 2])
range(res_normal01_simpred_nhhssi$theta_update[, 811:812, 2])

range(res_normal01_simpred_hhssi$theta_update[, 811:812, 3])
range(res_normal01_simpred_nhhssi$theta_update[, 811:812, 3])

range(res_normal01_simpred_hhssi$theta_update[, 811:812, 4])
range(res_normal01_simpred_nhhssi$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi$theta_update[, 811:812, k],
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
range(res_normal01_simpred_hhssi$theta_update[, 813:816, 1])
range(res_normal01_simpred_nhhssi$theta_update[, 813:816, 1])

plot(res_normal01_simpred_hhssi$theta_update[, 813:816, 1],
     res_normal01_simpred_nhhssi$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_simpred_hhssi$theta_update[, 813:816, ],
                    res_normal01_simpred_hhssi$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_simpred_hhssi$theta_update[, 813:816, 2],
       res_normal01_simpred_nhhssi$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$theta_update[, 813:816, 3],
       res_normal01_simpred_nhhssi$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$theta_update[, 813:816, 4],
       res_normal01_simpred_nhhssi$theta_update[, 813:816, 4],
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

range(res_normal01_simpred_hhssi$theta_update[, 813:816, 2])
range(res_normal01_simpred_nhhssi$theta_update[, 813:816, 2])

range(res_normal01_simpred_hhssi$theta_update[, 813:816, 3])
range(res_normal01_simpred_nhhssi$theta_update[, 813:816, 3])

range(res_normal01_simpred_hhssi$theta_update[, 813:816, 4])
range(res_normal01_simpred_nhhssi$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_simpred_hhssi$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_simpred_nhhssi$theta_update[, 813:816, k],
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
dim(sim_data_normal01_simpred$varphi_true)
dim(res_normal01_simpred_hhssi$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_simpred_hhssi$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_simpred_nhhssi$theta_update[, 811:816, ], c(2,3), mean)
# Bias for each outcome
varphi_bias_hhssi <- (varphi_est_hhssi - sim_data_normal01_simpred$varphi_true)
varphi_bias_nhhssi <- (varphi_est_nhhssi - sim_data_normal01_simpred$varphi_true)

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






varphi_var_hhssi <- apply(res_normal01_simpred_hhssi$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_simpred_nhhssi$theta_update[, 811:816, ], c(2,3), var)

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











dim(sim_data_normal01_simpred$gamma_true)
dim(res_normal01_simpred_hhssi$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_simpred_hhssi$lambdasq_delta_update)
dim(res_normal01_simpred_nhhssi$lambdasq_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_simpred_hhssi$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_simpred_nhhssi$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_simpred_hhssi$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_simpred_nhhssi$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_simpred_nhhssi$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_simpred_hhssi$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_simpred_nhhssi$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_simpred_hhssi$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_simpred_nhhssi$lambdasq_delta_update, 2, mean)
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
plot(apply(res_normal01_simpred_hhssi$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_simpred_hhssi$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_simpred_nhhssi$lambdasq_beta_update, 2, mean)
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
plot(apply(res_normal01_simpred_hhssi$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_simpred_hhssi$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi$lambdasq_gamma_update, 2, mean)
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




names(res_normal01_simpred_hhssi)
names(res_normal01_simpred_nhhssi)





# tausq beta
mean(res_normal01_simpred_hhssi$tausq_beta_update)
mean(res_normal01_simpred_nhhssi$tausq_beta_update)
abs(mean(res_normal01_simpred_hhssi$tausq_beta_update) -
      mean(res_normal01_simpred_nhhssi$tausq_beta_update))

range(res_normal01_simpred_hhssi$tausq_beta_update)
range(res_normal01_simpred_nhhssi$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_simpred_hhssi$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_simpred_nhhssi$tausq_beta_update
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
mean(res_normal01_simpred_hhssi$tausq_gamma_update)
mean(res_normal01_simpred_nhhssi$tausq_gamma_update)
abs(mean(res_normal01_simpred_hhssi$tausq_gamma_update) -
      mean(res_normal01_simpred_nhhssi$tausq_gamma_update))

range(res_normal01_simpred_hhssi$tausq_gamma_update)
range(res_normal01_simpred_nhhssi$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_simpred_hhssi$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_simpred_nhhssi$tausq_gamma_update
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
mean(res_normal01_simpred_hhssi$tausq_delta_update)
mean(res_normal01_simpred_nhhssi$tausq_delta_update)
abs(mean(res_normal01_simpred_hhssi$tausq_delta_update) -
      mean(res_normal01_simpred_nhhssi$tausq_delta_update))

range(res_normal01_simpred_hhssi$tausq_delta_update)
range(res_normal01_simpred_nhhssi$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_simpred_hhssi$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi$tausq_delta_update
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
mean(res_normal01_simpred_hhssi$xi_beta_update)
mean(res_normal01_simpred_nhhssi$xi_beta_update)
abs(mean(res_normal01_simpred_hhssi$xi_beta_update) -
      mean(res_normal01_simpred_nhhssi$xi_beta_update))

range(res_normal01_simpred_hhssi$xi_beta_update)
range(res_normal01_simpred_nhhssi$xi_beta_update)

xi_beta_hhssi <- res_normal01_simpred_hhssi$xi_beta_update
xi_beta_nhhssi <- res_normal01_simpred_nhhssi$xi_beta_update
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
mean(res_normal01_simpred_hhssi$xi_gamma_update)
mean(res_normal01_simpred_nhhssi$xi_gamma_update)
abs(mean(res_normal01_simpred_hhssi$xi_gamma_update) -
      mean(res_normal01_simpred_nhhssi$xi_gamma_update))

range(res_normal01_simpred_hhssi$xi_gamma_update)
range(res_normal01_simpred_nhhssi$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_simpred_hhssi$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_simpred_nhhssi$xi_gamma_update
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
mean(res_normal01_simpred_hhssi$xi_delta_update)
mean(res_normal01_simpred_nhhssi$xi_delta_update)
abs(mean(res_normal01_simpred_hhssi$xi_delta_update) -
      mean(res_normal01_simpred_nhhssi$xi_delta_update))

range(res_normal01_simpred_hhssi$xi_delta_update)
range(res_normal01_simpred_nhhssi$xi_delta_update)

xi_delta_hhssi <- res_normal01_simpred_hhssi$xi_delta_update
xi_delta_nhhssi <- res_normal01_simpred_nhhssi$xi_delta_update
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
dim(res_normal01_simpred_nhhssi$psi_delta_update)
dim(res_normal01_simpred_hhssi$psi_delta_update)

# active deltas
matplot(res_normal01_simpred_hhssi$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_simpred_hhssi$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_simpred_nhhssi$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_simpred_hhssi$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_simpred_nhhssi$psi_delta_update, 2, mean)

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
plot(apply(res_normal01_simpred_hhssi$psi_beta_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_simpred_hhssi$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_simpred_nhhssi$psi_beta_update, 2, mean)
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
plot(apply(res_normal01_simpred_hhssi$psi_gamma_update, 2, mean))
plot(apply(res_normal01_simpred_nhhssi$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_simpred_hhssi$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_simpred_nhhssi$psi_gamma_update, 2, mean)
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


names(res_normal01_simpred_hhssi)
names(res_normal01_simpred_nhhssi)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_simpred_hhssi$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_simpred_nhhssi$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_simpred_hhssi$Sigma_update[, , 1], type = "l")
matplot(res_normal01_simpred_nhhssi$Sigma_update[, , 1], type = "l")

matplot(res_normal01_simpred_hhssi$Sigma_update[, , 2], type = "l")
matplot(res_normal01_simpred_nhhssi$Sigma_update[, , 2], type = "l")

matplot(res_normal01_simpred_hhssi$Sigma_update[, , 3], type = "l")
matplot(res_normal01_simpred_nhhssi$Sigma_update[, , 3], type = "l")

matplot(res_normal01_simpred_hhssi$Sigma_update[, , 4], type = "l")
matplot(res_normal01_simpred_nhhssi$Sigma_update[, , 4], type = "l")



plot(res_normal01_simpred_hhssi$Sigma_update[, , 1],
     res_normal01_simpred_nhhssi$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_simpred_hhssi$Sigma_update,
                    res_normal01_simpred_hhssi$Sigma_update)),
     ylim = range(c(res_normal01_simpred_hhssi$Sigma_update,
                    res_normal01_simpred_hhssi$Sigma_update)))
# k = 2
points(res_normal01_simpred_hhssi$Sigma_update[, , 2],
       res_normal01_simpred_nhhssi$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_simpred_hhssi$Sigma_update[, , 3],
       res_normal01_simpred_nhhssi$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_simpred_hhssi$Sigma_update[, , 4],
       res_normal01_simpred_nhhssi$Sigma_update[, , 4],
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


dim(res_normal01_simpred_hhssi$lambdasq_beta_update)
dim(res_normal01_simpred_hhssi$lambdasq_gamma_update)
dim(res_normal01_simpred_hhssi$lambdasq_delta_update)
length(res_normal01_simpred_hhssi$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_simpred_hhssi$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_simpred_hhssi$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_simpred_hhssi$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_simpred_hhssi$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_simpred_nhhssi$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_simpred_nhhssi$tausq_delta_update


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



beta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_simpred_hhssi$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_simpred_hhssi$tausq_delta_update,
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
# β, γ, δ, φ ∼ normal(0,1), alpha ~ unif(-3.5, -2.5), all predictors real
###############################################################################

# file alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5
set.seed(23444)
sim_data_normal01_realpred <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                O=O, x=X, z=Z, d=D,
                                                Sigma=Sigma)

Y <- sim_data_normal01_realpred$Y
x <- sim_data_normal01_realpred$x
z <- sim_data_normal01_realpred$z
u <- sim_data_normal01_realpred$u
d <- sim_data_normal01_realpred$d
W <- sim_data_normal01_realpred$W
WTW <- t(W) %*% W
n_all_par <- ncol(W)

# load(file="practice_data/res_normal01_realpred_hhssi.rda")
# load(file="practice_data/res_normal01_realpred_nhhssi.rda")

sim_data_normal01_realpred$delta_true
dim(res_normal01_realpred_hhssi$theta_update)






# alpha
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(res_normal01_realpred_hhssi$theta_update[, 1, k],
          type = "l",
          lty = 1,
          lwd = 2,
          xlab = "MCMC iteration",
          ylab = expression(alpha),
          main = paste("Outcome", k)
  )

  matlines(res_normal01_realpred_nhhssi$theta_update[, 1, k], lty = 2, lwd = 2)

  legend("topright", legend = c("HHS-SI", "NHHS-SI"), lty = c(1, 2), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))

range(res_normal01_realpred_hhssi$theta_update[, 1, ])
range(res_normal01_realpred_nhhssi$theta_update[, 1, ])

par(mfrow = c(2, 2))
for (k in 1:4) {
  hhssi_alpha <- res_normal01_realpred_hhssi$theta_update[, 1, k]
  nhhssi_alpha <- res_normal01_realpred_nhhssi$theta_update[, 1, k]

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
res_normal01_realpred_hhssi$theta_update[, 1, ]

# Posterior mean estimate of alpha
alpha_est_hhssi <- apply(res_normal01_realpred_hhssi$theta_update[, 1, ], 2, mean)
alpha_est_nhhssi <- apply(res_normal01_realpred_nhhssi$theta_update[, 1, ], 2, mean)
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

plot(res_normal01_realpred_hhssi$theta_update[, 2:3, 1],
     res_normal01_realpred_nhhssi$theta_update[, 2:3, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi$theta_update[, 2:3, ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, 2:3, ],
                    res_normal01_realpred_hhssi$theta_update[, 2:3, ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, 2:3, 2],
       res_normal01_realpred_nhhssi$theta_update[, 2:3, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, 2:3, 3],
       res_normal01_realpred_nhhssi$theta_update[, 2:3, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, 2:3, 4],
       res_normal01_realpred_nhhssi$theta_update[, 2:3, 4],
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




range(res_normal01_realpred_hhssi$theta_update[, 2:3, 2])
range(res_normal01_realpred_nhhssi$theta_update[, 2:3, 2])

range(res_normal01_realpred_hhssi$theta_update[, 2:3, 3])
range(res_normal01_realpred_nhhssi$theta_update[, 2:3, 3])

range(res_normal01_realpred_hhssi$theta_update[, 2:3, 4])
range(res_normal01_realpred_nhhssi$theta_update[, 2:3, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, 2:3, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active beta",
    main = paste("Active beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi$theta_update[, 2:3, k],
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
range(res_normal01_realpred_hhssi$theta_update[, 4:45, 1])
range(res_normal01_realpred_nhhssi$theta_update[, 4:45, 1])

plot(res_normal01_realpred_hhssi$theta_update[, 4:45, 1],
     res_normal01_realpred_nhhssi$theta_update[, 4:45, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi$theta_update[, 4:45, ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, 4:45, ],
                    res_normal01_realpred_hhssi$theta_update[, 4:45, ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, 4:45, 2],
       res_normal01_realpred_nhhssi$theta_update[, 4:45, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, 4:45, 3],
       res_normal01_realpred_nhhssi$theta_update[, 4:45, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, 4:45, 4],
       res_normal01_realpred_nhhssi$theta_update[, 4:45, 4],
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



range(res_normal01_realpred_hhssi$theta_update[, 4:45, 2])
range(res_normal01_realpred_nhhssi$theta_update[, 4:45, 2])

range(res_normal01_realpred_hhssi$theta_update[, 4:45, 3])
range(res_normal01_realpred_nhhssi$theta_update[, 4:45, 3])

range(res_normal01_realpred_hhssi$theta_update[, 4:45, 4])
range(res_normal01_realpred_nhhssi$theta_update[, 4:45, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, 4:45, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive beta",
    main = paste("Inactive beta traceplots: k =", k)
  )

  matlines(
    res_normal01_realpred_nhhssi$theta_update[, 4:45, k],
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
dim(res_normal01_realpred_hhssi$theta_update[, 2:3, ])

# Posterior mean estimate of beta
beta_est_hhssi <- apply(res_normal01_realpred_hhssi$theta_update[, 2:45, ], c(2,3), mean)
beta_est_nhhssi <- apply(res_normal01_realpred_nhhssi$theta_update[, 2:45, ], c(2,3), mean)
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



beta_var_hhssi <- (res_normal01_realpred_hhssi$lambdasq_beta_update *
                     res_normal01_realpred_hhssi$tausq_beta_update)
beta_var_nhhssi <- (res_normal01_realpred_nhhssi$lambdasq_beta_update *
                      res_normal01_realpred_nhhssi$tausq_beta_update)

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
range(res_normal01_realpred_hhssi$theta_update[, 46:48, 1])
range(res_normal01_realpred_nhhssi$theta_update[, 46:48, 1])


plot(res_normal01_realpred_hhssi$theta_update[, 46:48, 1],
     res_normal01_realpred_nhhssi$theta_update[, 46:48, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi$theta_update[, 46:48, ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, 46:48, ],
                    res_normal01_realpred_hhssi$theta_update[, 46:48, ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, 46:48, 2],
       res_normal01_realpred_nhhssi$theta_update[, 46:48, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, 46:48, 3],
       res_normal01_realpred_nhhssi$theta_update[, 46:48, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, 46:48, 4],
       res_normal01_realpred_nhhssi$theta_update[, 46:48, 4],
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

range(res_normal01_realpred_hhssi$theta_update[, 46:48, 2])
range(res_normal01_realpred_nhhssi$theta_update[, 46:48, 2])

range(res_normal01_realpred_hhssi$theta_update[, 46:48, 3])
range(res_normal01_realpred_nhhssi$theta_update[, 46:48, 3])

range(res_normal01_realpred_hhssi$theta_update[, 46:48, 4])
range(res_normal01_realpred_nhhssi$theta_update[, 46:48, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, 46:48, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active gamma",
    main = paste("active gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi$theta_update[, 46:48, k],
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
range(res_normal01_realpred_hhssi$theta_update[, 49:62, 1])
range(res_normal01_realpred_nhhssi$theta_update[, 49:62, 1])


plot(res_normal01_realpred_hhssi$theta_update[, 49:62, 1],
     res_normal01_realpred_nhhssi$theta_update[, 49:62, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi$theta_update[, 49:62, ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, 49:62, ],
                    res_normal01_realpred_hhssi$theta_update[, 49:62, ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, 49:62, 2],
       res_normal01_realpred_nhhssi$theta_update[, 49:62, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, 49:62, 3],
       res_normal01_realpred_nhhssi$theta_update[, 49:62, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, 49:62, 4],
       res_normal01_realpred_nhhssi$theta_update[, 49:62, 4],
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


range(res_normal01_realpred_hhssi$theta_update[, 49:62, 2])
range(res_normal01_realpred_nhhssi$theta_update[, 49:62, 2])

range(res_normal01_realpred_hhssi$theta_update[, 49:62, 3])
range(res_normal01_realpred_nhhssi$theta_update[, 49:62, 3])

range(res_normal01_realpred_hhssi$theta_update[, 49:62, 4])
range(res_normal01_realpred_nhhssi$theta_update[, 49:62, 4])




par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, 49:62, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive gamma",
    main = paste("inactive gamma traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi$theta_update[, 49:62, k],
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
dim(res_normal01_realpred_hhssi$theta_update[, 46:62, ])

# Posterior mean estimate of gamma
gamma_est_hhssi <- apply(res_normal01_realpred_hhssi$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_nhhssi <- apply(res_normal01_realpred_nhhssi$theta_update[, 46:62, ], c(2,3), mean)
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



gamma_var_hhssi <- (res_normal01_realpred_hhssi$lambdasq_gamma_update *
                      res_normal01_realpred_hhssi$tausq_gamma_update)
gamma_var_nhhssi <- (res_normal01_realpred_nhhssi$lambdasq_gamma_update *
                       res_normal01_realpred_nhhssi$tausq_gamma_update)

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
range(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 1])
range(res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 1])

plot(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 1],
     res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), ],
                    res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 2],
       res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 3],
       res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 4],
       res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 4],
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


range(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 2])
range(res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 2])

range(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 3])
range(res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 3])

range(res_normal01_realpred_hhssi$theta_update[, c(63:65, 80:82), 4])
range(res_normal01_realpred_nhhssi$theta_update[, c(63:65, 80:82), 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, c(63:65,80:82), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active delta",
    main = paste("active delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi$theta_update[, c(63:65,80:82), k],
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

range(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 1])
range(res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 1])

plot(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 1],
     res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), ],
                    res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 2],
       res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 3],
       res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 4],
       res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 4],
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


range(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 2])
range(res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 2])

range(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 3])
range(res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 3])

range(res_normal01_realpred_hhssi$theta_update[, c(66:79, 83:810), 4])
range(res_normal01_realpred_nhhssi$theta_update[, c(66:79, 83:810), 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, c(66:79,83:810), k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive delta",
    main = paste("inactive delta traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi$theta_update[, c(66:79,83:810), k],
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
dim(res_normal01_realpred_hhssi$theta_update[, 63:810, ])

# Posterior mean estimate of delta
delta_est_hhssi <- apply(res_normal01_realpred_hhssi$theta_update[, 63:810, ], c(2,3), mean)
delta_est_nhhssi <- apply(res_normal01_realpred_nhhssi$theta_update[, 63:810, ], c(2,3), mean)
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
range(res_normal01_realpred_hhssi$theta_update[, 811:812, 1])
range(res_normal01_realpred_nhhssi$theta_update[, 811:812, 1])

plot(res_normal01_realpred_hhssi$theta_update[, 811:812, 1],
     res_normal01_realpred_nhhssi$theta_update[, 811:812, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi$theta_update[, 811:812, ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, 811:812, ],
                    res_normal01_realpred_hhssi$theta_update[, 811:812, ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, 811:812, 2],
       res_normal01_realpred_nhhssi$theta_update[, 811:812, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, 811:812, 3],
       res_normal01_realpred_nhhssi$theta_update[, 811:812, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, 811:812, 4],
       res_normal01_realpred_nhhssi$theta_update[, 811:812, 4],
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


range(res_normal01_realpred_hhssi$theta_update[, 811:812, 2])
range(res_normal01_realpred_nhhssi$theta_update[, 811:812, 2])

range(res_normal01_realpred_hhssi$theta_update[, 811:812, 3])
range(res_normal01_realpred_nhhssi$theta_update[, 811:812, 3])

range(res_normal01_realpred_hhssi$theta_update[, 811:812, 4])
range(res_normal01_realpred_nhhssi$theta_update[, 811:812, 4])



par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, 811:812, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "active varphi",
    main = paste("active varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi$theta_update[, 811:812, k],
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
range(res_normal01_realpred_hhssi$theta_update[, 813:816, 1])
range(res_normal01_realpred_nhhssi$theta_update[, 813:816, 1])

plot(res_normal01_realpred_hhssi$theta_update[, 813:816, 1],
     res_normal01_realpred_nhhssi$theta_update[, 813:816, 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi$theta_update[, 813:816, ])),
     ylim = range(c(res_normal01_realpred_hhssi$theta_update[, 813:816, ],
                    res_normal01_realpred_hhssi$theta_update[, 813:816, ])))
# k = 2
points(res_normal01_realpred_hhssi$theta_update[, 813:816, 2],
       res_normal01_realpred_nhhssi$theta_update[, 813:816, 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$theta_update[, 813:816, 3],
       res_normal01_realpred_nhhssi$theta_update[, 813:816, 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$theta_update[, 813:816, 4],
       res_normal01_realpred_nhhssi$theta_update[, 813:816, 4],
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

range(res_normal01_realpred_hhssi$theta_update[, 813:816, 2])
range(res_normal01_realpred_nhhssi$theta_update[, 813:816, 2])

range(res_normal01_realpred_hhssi$theta_update[, 813:816, 3])
range(res_normal01_realpred_nhhssi$theta_update[, 813:816, 3])

range(res_normal01_realpred_hhssi$theta_update[, 813:816, 4])
range(res_normal01_realpred_nhhssi$theta_update[, 813:816, 4])


par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(
    res_normal01_realpred_hhssi$theta_update[, 813:816, k],
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "MCMC iteration",
    ylab = "inactive varphi",
    main = paste("inactive varphi traceplots: k =", k))
  matlines(
    res_normal01_realpred_nhhssi$theta_update[, 813:816, k],
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
dim(res_normal01_realpred_hhssi$theta_update[, 811:816, ])

# Posterior mean estimate of varphi
varphi_est_hhssi <- apply(res_normal01_realpred_hhssi$theta_update[, 811:816, ], c(2,3), mean)
varphi_est_nhhssi <- apply(res_normal01_realpred_nhhssi$theta_update[, 811:816, ], c(2,3), mean)
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






varphi_var_hhssi <- apply(res_normal01_realpred_hhssi$theta_update[, 811:816, ], c(2,3), var)
varphi_var_nhhssi <- apply(res_normal01_realpred_nhhssi$theta_update[, 811:816, ], c(2,3), var)

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
dim(res_normal01_realpred_hhssi$theta_update)





# lambdasq_delta, 1:3 and 18:20 are active deltas
dim(res_normal01_realpred_hhssi$lambdasq_delta_update)
dim(res_normal01_realpred_nhhssi$lambdasq_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi$lambdasq_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi$lambdasq_delta_update[,-c(1:3, 18:20)], type = "l")


# active deltas
matplot(res_normal01_realpred_hhssi$lambdasq_delta_update[,1:3], type = "l")
matplot(res_normal01_realpred_nhhssi$lambdasq_delta_update[,1:3], type = "l")

matplot(res_normal01_realpred_hhssi$lambdasq_delta_update[,18:20], type = "l")
matplot(res_normal01_realpred_nhhssi$lambdasq_delta_update[,18:20], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi$lambdasq_delta_update[,4:6], type = "l")
matplot(res_normal01_realpred_nhhssi$lambdasq_delta_update[,4:6], type = "l")

matplot(res_normal01_realpred_hhssi$lambdasq_delta_update[,21:23], type = "l")
matplot(res_normal01_realpred_nhhssi$lambdasq_delta_update[,21:23], type = "l")


lambdasq_delta_hhssi <- apply(res_normal01_realpred_hhssi$lambdasq_delta_update, 2, mean)
lambdasq_delta_nhhssi <- apply(res_normal01_realpred_nhhssi$lambdasq_delta_update, 2, mean)
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
plot(apply(res_normal01_realpred_hhssi$lambdasq_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi$lambdasq_beta_update, 2, mean))

lambdasq_beta_hhssi <- apply(res_normal01_realpred_hhssi$lambdasq_beta_update, 2, mean)
lambdasq_beta_nhhssi <- apply(res_normal01_realpred_nhhssi$lambdasq_beta_update, 2, mean)
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
plot(apply(res_normal01_realpred_hhssi$lambdasq_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi$lambdasq_gamma_update, 2, mean))

lambdasq_gamma_hhssi <- apply(res_normal01_realpred_hhssi$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi$lambdasq_gamma_update, 2, mean)
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




names(res_normal01_realpred_hhssi)
names(res_normal01_realpred_nhhssi)





# tausq beta
mean(res_normal01_realpred_hhssi$tausq_beta_update)
mean(res_normal01_realpred_nhhssi$tausq_beta_update)
abs(mean(res_normal01_realpred_hhssi$tausq_beta_update) -
      mean(res_normal01_realpred_nhhssi$tausq_beta_update))

range(res_normal01_realpred_hhssi$tausq_beta_update)
range(res_normal01_realpred_nhhssi$tausq_beta_update)

tausq_beta_hhssi <- res_normal01_realpred_hhssi$tausq_beta_update
tausq_beta_nhhssi <- res_normal01_realpred_nhhssi$tausq_beta_update
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
mean(res_normal01_realpred_hhssi$tausq_gamma_update)
mean(res_normal01_realpred_nhhssi$tausq_gamma_update)
abs(mean(res_normal01_realpred_hhssi$tausq_gamma_update) -
      mean(res_normal01_realpred_nhhssi$tausq_gamma_update))

range(res_normal01_realpred_hhssi$tausq_gamma_update)
range(res_normal01_realpred_nhhssi$tausq_gamma_update)

tausq_gamma_hhssi <- res_normal01_realpred_hhssi$tausq_gamma_update
tausq_gamma_nhhssi <- res_normal01_realpred_nhhssi$tausq_gamma_update
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
mean(res_normal01_realpred_hhssi$tausq_delta_update)
mean(res_normal01_realpred_nhhssi$tausq_delta_update)
abs(mean(res_normal01_realpred_hhssi$tausq_delta_update) -
      mean(res_normal01_realpred_nhhssi$tausq_delta_update))

range(res_normal01_realpred_hhssi$tausq_delta_update)
range(res_normal01_realpred_nhhssi$tausq_delta_update)

tausq_delta_hhssi <- res_normal01_realpred_hhssi$tausq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi$tausq_delta_update
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
mean(res_normal01_realpred_hhssi$xi_beta_update)
mean(res_normal01_realpred_nhhssi$xi_beta_update)
abs(mean(res_normal01_realpred_hhssi$xi_beta_update) -
      mean(res_normal01_realpred_nhhssi$xi_beta_update))

range(res_normal01_realpred_hhssi$xi_beta_update)
range(res_normal01_realpred_nhhssi$xi_beta_update)

xi_beta_hhssi <- res_normal01_realpred_hhssi$xi_beta_update
xi_beta_nhhssi <- res_normal01_realpred_nhhssi$xi_beta_update
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
mean(res_normal01_realpred_hhssi$xi_gamma_update)
mean(res_normal01_realpred_nhhssi$xi_gamma_update)
abs(mean(res_normal01_realpred_hhssi$xi_gamma_update) -
      mean(res_normal01_realpred_nhhssi$xi_gamma_update))

range(res_normal01_realpred_hhssi$xi_gamma_update)
range(res_normal01_realpred_nhhssi$xi_gamma_update)

xi_gamma_hhssi <- res_normal01_realpred_hhssi$xi_gamma_update
xi_gamma_nhhssi <- res_normal01_realpred_nhhssi$xi_gamma_update
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
mean(res_normal01_realpred_hhssi$xi_delta_update)
mean(res_normal01_realpred_nhhssi$xi_delta_update)
abs(mean(res_normal01_realpred_hhssi$xi_delta_update) -
      mean(res_normal01_realpred_nhhssi$xi_delta_update))

range(res_normal01_realpred_hhssi$xi_delta_update)
range(res_normal01_realpred_nhhssi$xi_delta_update)

xi_delta_hhssi <- res_normal01_realpred_hhssi$xi_delta_update
xi_delta_nhhssi <- res_normal01_realpred_nhhssi$xi_delta_update
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
dim(res_normal01_realpred_nhhssi$psi_delta_update)
dim(res_normal01_realpred_hhssi$psi_delta_update)

# active deltas
matplot(res_normal01_realpred_hhssi$psi_delta_update[,c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi$psi_delta_update[,c(1:3, 18:20)], type = "l")

# inactive deltas
matplot(res_normal01_realpred_hhssi$psi_delta_update[,-c(1:3, 18:20)], type = "l")
matplot(res_normal01_realpred_nhhssi$psi_delta_update[,-c(1:3, 18:20)], type = "l")


psi_delta_hhssi <- apply(res_normal01_realpred_hhssi$psi_delta_update, 2, mean)
psi_delta_nhhssi <- apply(res_normal01_realpred_nhhssi$psi_delta_update, 2, mean)

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
plot(apply(res_normal01_realpred_hhssi$psi_beta_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi$psi_beta_update, 2, mean))

psi_beta_hhssi <- apply(res_normal01_realpred_hhssi$psi_beta_update, 2, mean)
psi_beta_nhhssi <- apply(res_normal01_realpred_nhhssi$psi_beta_update, 2, mean)
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
plot(apply(res_normal01_realpred_hhssi$psi_gamma_update, 2, mean))
plot(apply(res_normal01_realpred_nhhssi$psi_gamma_update, 2, mean))

psi_gamma_hhssi <- apply(res_normal01_realpred_hhssi$psi_gamma_update, 2, mean)
psi_gamma_nhhssi <- apply(res_normal01_realpred_nhhssi$psi_gamma_update, 2, mean)
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


names(res_normal01_realpred_hhssi)
names(res_normal01_realpred_nhhssi)



# Sigma
# Posterior mean Sigma for HHS-SI
Sigma_mean_hhssi <- apply(res_normal01_realpred_hhssi$Sigma_update,
                          c(2, 3), mean)

# Posterior mean Sigma for NHHS-SI
Sigma_mean_nhhssi <- apply(res_normal01_realpred_nhhssi$Sigma_update,
                           c(2, 3), mean)

(Sigma_mean_hhssi-
    Sigma_mean_nhhssi)

range(Sigma_mean_hhssi)
range(Sigma_mean_nhhssi)

matplot(Sigma_mean_hhssi, type = "l")
matplot(Sigma_mean_nhhssi, type = "l")

matplot(res_normal01_realpred_hhssi$Sigma_update[, , 1], type = "l")
matplot(res_normal01_realpred_nhhssi$Sigma_update[, , 1], type = "l")

matplot(res_normal01_realpred_hhssi$Sigma_update[, , 2], type = "l")
matplot(res_normal01_realpred_nhhssi$Sigma_update[, , 2], type = "l")

matplot(res_normal01_realpred_hhssi$Sigma_update[, , 3], type = "l")
matplot(res_normal01_realpred_nhhssi$Sigma_update[, , 3], type = "l")

matplot(res_normal01_realpred_hhssi$Sigma_update[, , 4], type = "l")
matplot(res_normal01_realpred_nhhssi$Sigma_update[, , 4], type = "l")



plot(res_normal01_realpred_hhssi$Sigma_update[, , 1],
     res_normal01_realpred_nhhssi$Sigma_update[, , 1],
     pch = 19,
     xlab = "HHS-SI",
     ylab = "NHHS-SI",
     col = "yellow",
     xlim = range(c(res_normal01_realpred_hhssi$Sigma_update,
                    res_normal01_realpred_hhssi$Sigma_update)),
     ylim = range(c(res_normal01_realpred_hhssi$Sigma_update,
                    res_normal01_realpred_hhssi$Sigma_update)))
# k = 2
points(res_normal01_realpred_hhssi$Sigma_update[, , 2],
       res_normal01_realpred_nhhssi$Sigma_update[, , 2],
       pch = 19,
       col = "lightblue")
# k = 3
points(res_normal01_realpred_hhssi$Sigma_update[, , 3],
       res_normal01_realpred_nhhssi$Sigma_update[, , 3],
       pch = 19,
       col = "green")
points(res_normal01_realpred_hhssi$Sigma_update[, , 4],
       res_normal01_realpred_nhhssi$Sigma_update[, , 4],
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


dim(res_normal01_realpred_hhssi$lambdasq_beta_update)
dim(res_normal01_realpred_hhssi$lambdasq_gamma_update)
dim(res_normal01_realpred_hhssi$lambdasq_delta_update)
length(res_normal01_realpred_hhssi$tausq_delta_update)



beta_lam_hhssi  <- res_normal01_realpred_hhssi$lambdasq_beta_update
gamma_lam_hhssi <- res_normal01_realpred_hhssi$lambdasq_gamma_update
delta_lam_hhssi <- res_normal01_realpred_hhssi$lambdasq_delta_update
tausq_delta_hhssi <- res_normal01_realpred_hhssi$tausq_delta_update


delta_lam_nhhssi  <- res_normal01_realpred_nhhssi$lambdasq_delta_update
tausq_delta_nhhssi <- res_normal01_realpred_nhhssi$tausq_delta_update


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



beta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi$lambdasq_beta_update, 2, mean)
plot(beta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Beta parameter",
     ylab = expression(paste("Posterior mean ", lambda[beta]^2)),
     main = "HHS-SI: beta local shrinkage (lambdasq beta)")


gamma_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi$lambdasq_gamma_update, 2, mean)
plot(gamma_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[gamma]^2)),
     main = "HHS-SI: gamma local shrinkage (lambdasq gamma)")


delta_lambda_mean_hhssi <- apply(res_normal01_realpred_hhssi$lambdasq_delta_update, 2, mean)
plot(delta_lambda_mean_hhssi,
     pch = 19,
     xlab = "Gamma parameter",
     ylab = expression(paste("Posterior mean ", lambda[delta]^2)),
     main = "HHS-SI: delta local shrinkage (lambdasq delta)")
hist(delta_lambda_mean_hhssi, breaks = 1000)



plot(res_normal01_realpred_hhssi$tausq_delta_update,
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





