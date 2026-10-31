
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


###############################################################################################
# Combined plot
###############################################################################################


# First sample
set.seed(23444)
sim_data_unif23_simpred1 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n ,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)
# load(file="practice_data/res_unif23_simpred_hhssi1.rda")
# load(file="practice_data/res_unif23_simpred_nhhssi1.rda")
set.seed(23444)
sim_data_unif23_realpred1 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                 O=O, x=X, z=Z, d=D,
                                                 Sigma=Sigma)
# load(file="practice_data/res_unif23_realpred_hhssi1.rda")
# load(file="practice_data/res_unif23_realpred_nhhssi1.rda")





# Second sample
set.seed(654586798)
sim_data_unif23_simpred2 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n ,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)
# load(file="practice_data/res_unif23_simpred_hhssi2.rda")
# load(file="practice_data/res_unif23_simpred_nhhssi2.rda")
set.seed(654586798)
sim_data_unif23_realpred2 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                 O=O, x=X, z=Z, d=D,
                                                 Sigma=Sigma)
# load(file="practice_data/res_unif23_realpred_hhssi2.rda")
# load(file="practice_data/res_unif23_realpred_nhhssi2.rda")





# Third sample
set.seed(23315644)
sim_data_unif23_simpred3 <- Sim_data_BVS_real_sim_var(
  K = K,
  n = n ,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)
# load(file="practice_data/res_unif23_simpred_hhssi3.rda")
# load(file="practice_data/res_unif23_simpred_nhhssi3.rda")
set.seed(23315644)
sim_data_unif23_realpred3 <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                                                 O=O, x=X, z=Z, d=D,
                                                 Sigma=Sigma)
# load(file="practice_data/res_unif23_realpred_hhssi3.rda")
# load(file="practice_data/res_unif23_realpred_nhhssi3.rda")







#fs sim
beta_lam_simpred_hhssi1  <- res_unif23_simpred_hhssi1$lambdasq_beta_update
gamma_lam_simpred_hhssi1 <- res_unif23_simpred_hhssi1$lambdasq_gamma_update
delta_lam_simpred_hhssi1 <- res_unif23_simpred_hhssi1$lambdasq_delta_update
tausq_delta_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_delta_update

delta_lam_simpred_nhhssi1  <- res_unif23_simpred_nhhssi1$lambdasq_delta_update
tausq_delta_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_delta_update

#fs real
beta_lam_realpred_hhssi1  <- res_unif23_realpred_hhssi1$lambdasq_beta_update
gamma_lam_realpred_hhssi1 <- res_unif23_realpred_hhssi1$lambdasq_gamma_update
delta_lam_realpred_hhssi1 <- res_unif23_realpred_hhssi1$lambdasq_delta_update
tausq_delta_realpred_hhssi1 <- res_unif23_realpred_hhssi1$tausq_delta_update

delta_lam_realpred_nhhssi1  <- res_unif23_realpred_nhhssi1$lambdasq_delta_update
tausq_delta_realpred_nhhssi1 <- res_unif23_realpred_nhhssi1$tausq_delta_update


#ss sim
beta_lam_simpred_hhssi2  <- res_unif23_simpred_hhssi2$lambdasq_beta_update
gamma_lam_simpred_hhssi2 <- res_unif23_simpred_hhssi2$lambdasq_gamma_update
delta_lam_simpred_hhssi2 <- res_unif23_simpred_hhssi2$lambdasq_delta_update
tausq_delta_simpred_hhssi2 <- res_unif23_simpred_hhssi2$tausq_delta_update

delta_lam_simpred_nhhssi2  <- res_unif23_simpred_nhhssi2$lambdasq_delta_update
tausq_delta_simpred_nhhssi2 <- res_unif23_simpred_nhhssi2$tausq_delta_update

#ss real
beta_lam_realpred_hhssi2  <- res_unif23_realpred_hhssi2$lambdasq_beta_update
gamma_lam_realpred_hhssi2 <- res_unif23_realpred_hhssi2$lambdasq_gamma_update
delta_lam_realpred_hhssi2 <- res_unif23_realpred_hhssi2$lambdasq_delta_update
tausq_delta_realpred_hhssi2 <- res_unif23_realpred_hhssi2$tausq_delta_update

delta_lam_realpred_nhhssi2  <- res_unif23_realpred_nhhssi2$lambdasq_delta_update
tausq_delta_realpred_nhhssi2 <- res_unif23_realpred_nhhssi2$tausq_delta_update


#ts sim
beta_lam_simpred_hhssi3  <- res_unif23_simpred_hhssi3$lambdasq_beta_update
gamma_lam_simpred_hhssi3 <- res_unif23_simpred_hhssi3$lambdasq_gamma_update
delta_lam_simpred_hhssi3 <- res_unif23_simpred_hhssi3$lambdasq_delta_update
tausq_delta_simpred_hhssi3 <- res_unif23_simpred_hhssi3$tausq_delta_update

delta_lam_simpred_nhhssi3  <- res_unif23_simpred_nhhssi3$lambdasq_delta_update
tausq_delta_simpred_nhhssi3 <- res_unif23_simpred_nhhssi3$tausq_delta_update

#ts real
beta_lam_realpred_hhssi3  <- res_unif23_realpred_hhssi3$lambdasq_beta_update
gamma_lam_realpred_hhssi3 <- res_unif23_realpred_hhssi3$lambdasq_gamma_update
delta_lam_realpred_hhssi3 <- res_unif23_realpred_hhssi3$lambdasq_delta_update
tausq_delta_realpred_hhssi3 <- res_unif23_realpred_hhssi3$tausq_delta_update

delta_lam_realpred_nhhssi3  <- res_unif23_realpred_nhhssi3$lambdasq_delta_update
tausq_delta_realpred_nhhssi3 <- res_unif23_realpred_nhhssi3$tausq_delta_update




n_iter <- dim(beta_lam_realpred_hhssi3)[1]



####################################################################
# for unif23_simpred_hhssi1, unif23_simpred_nhhssi1

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
lambdsq_beta_mean_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdsq_beta_mean_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$lambdasq_beta_update, 2, mean)
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
lambdsq_gamma_mean_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdsq_gamma_mean_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$lambdasq_gamma_update, 2, mean)
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
lambdsq_delta_mean_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdsq_delta_mean_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$lambdasq_delta_update, 2, mean)
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
tausq_beta_mean_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_beta_update
tausq_beta_mean_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_beta_update
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
tausq_gamma_mean_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_gamma_update
tausq_gamma_mean_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_gamma_update
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
tausq_delta_mean_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_delta_update
tausq_delta_mean_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_delta_update
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
lambdsq_tausq_delta_mean_simpred_hhssi1 <- apply((res_unif23_simpred_hhssi1$lambdasq_delta_update*
                                                 res_unif23_simpred_hhssi1$tausq_delta_update),
                                                 2, mean)
lambdsq_tausq_delta_mean_simpred_nhhssi1 <- apply((res_unif23_simpred_nhhssi1$lambdasq_delta_update*
                                                     res_unif23_simpred_nhhssi1$tausq_delta_update),
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
# for unif23_realpred_hhssi1, unif23_realpred_nhhssi1

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
lambdsq_beta_mean_realpred_hhssi1 <- apply(res_unif23_realpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdsq_beta_mean_realpred_nhhssi1 <- apply(res_unif23_realpred_nhhssi1$lambdasq_beta_update, 2, mean)
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
lambdsq_gamma_mean_realpred_hhssi1 <- apply(res_unif23_realpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdsq_gamma_mean_realpred_nhhssi1 <- apply(res_unif23_realpred_nhhssi1$lambdasq_gamma_update, 2, mean)
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
lambdsq_delta_mean_realpred_hhssi1 <- apply(res_unif23_realpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdsq_delta_mean_realpred_nhhssi1 <- apply(res_unif23_realpred_nhhssi1$lambdasq_delta_update, 2, mean)
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
tausq_beta_mean_realpred_hhssi1 <- res_unif23_realpred_hhssi1$tausq_beta_update
tausq_beta_mean_realpred_nhhssi1 <- res_unif23_realpred_nhhssi1$tausq_beta_update
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
tausq_gamma_mean_realpred_hhssi1 <- res_unif23_realpred_hhssi1$tausq_gamma_update
tausq_gamma_mean_realpred_nhhssi1 <- res_unif23_realpred_nhhssi1$tausq_gamma_update
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
tausq_delta_mean_realpred_hhssi1 <- res_unif23_realpred_hhssi1$tausq_delta_update
tausq_delta_mean_realpred_nhhssi1 <- res_unif23_realpred_nhhssi1$tausq_delta_update
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
lambdsq_tausq_delta_mean_realpred_hhssi1 <- apply((res_unif23_realpred_hhssi1$lambdasq_delta_update*
                                                    res_unif23_realpred_hhssi1$tausq_delta_update),
                                                 2, mean)
lambdsq_tausq_delta_mean_realpred_nhhssi1 <- apply((res_unif23_realpred_nhhssi1$lambdasq_delta_update*
                                                     res_unif23_realpred_nhhssi1$tausq_delta_update),
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
# for unif23_simpred_hhssi1, unif23_simpred_nhhssi1

#fs sim
beta_lam_simpred_hhssi1  <- res_unif23_simpred_hhssi1$lambdasq_beta_update
gamma_lam_simpred_hhssi1 <- res_unif23_simpred_hhssi1$lambdasq_gamma_update
delta_lam_simpred_hhssi1 <- res_unif23_simpred_hhssi1$lambdasq_delta_update
tausq_delta_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_delta_update

delta_lam_simpred_nhhssi1  <- res_unif23_simpred_nhhssi1$lambdasq_delta_update
tausq_delta_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_delta_update


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
delta_sd_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$theta_update[, 63:810, ], c(2,3), sd)
delta_sd_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$theta_update[, 63:810, ], c(2,3), sd)



# Posterior mean estimate of delta
delta_est_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$theta_update[, 63:810, ], c(2,3), mean)
delta_est_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$theta_update[, 63:810, ], c(2,3), mean)

# Bias for each outcome
delta_bias_unif23_simpred_hhssi1 <- (delta_est_unif23_simpred_hhssi1 - sim_data_unif23_simpred1$delta_true)
delta_bias_unif23_simpred_nhhssi1 <- (delta_est_unif23_simpred_nhhssi1 - sim_data_unif23_simpred1$delta_true)




##########################################################################################
# Delta Estimates of sd
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_unif23_simpred_hhssi1[, k],
    delta_sd_unif23_simpred_nhhssi1[, k],
    col =point_col,
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Overall sd: Outcome", k),
    xlim = range(c(delta_sd_unif23_simpred_hhssi1[,k],
                   delta_sd_unif23_simpred_nhhssi1[,k])),
    ylim = range(c(delta_sd_unif23_simpred_hhssi1[,k],
                   delta_sd_unif23_simpred_nhhssi1[,k])))
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
    delta_sd_unif23_simpred_hhssi1[c(1:3, 18:20), k],
    delta_sd_unif23_simpred_nhhssi1[c(1:3, 18:20), k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("active: Outcome", k),
    xlim = range(c(delta_sd_unif23_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_sd_unif23_simpred_nhhssi1[c(1:3, 18:20),k])),
    ylim = range(c(delta_sd_unif23_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_sd_unif23_simpred_nhhssi1[c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20), k],
    delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20), k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    col = point_col_inct,
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20),k])),
    ylim = range(c(delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20),k])))
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
    delta_sd_unif23_simpred_hhssi1[idx, k],
    delta_sd_unif23_simpred_nhhssi1[idx, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_sd_unif23_simpred_hhssi1[idx,k],
                   delta_sd_unif23_simpred_nhhssi1[idx,k])),
    ylim = range(c(delta_sd_unif23_simpred_hhssi1[idx,k],
                   delta_sd_unif23_simpred_nhhssi1[idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


# !idx
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_unif23_simpred_hhssi1[!idx, k],
    delta_sd_unif23_simpred_nhhssi1[!idx, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_sd_unif23_simpred_hhssi1[!idx,k],
                   delta_sd_unif23_simpred_nhhssi1[!idx,k])),
    ylim = range(c(delta_sd_unif23_simpred_hhssi1[!idx,k],
                   delta_sd_unif23_simpred_nhhssi1[!idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





idx_act <- delta_var_mean_simpred_hhssi1[c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)]
# no active delta with nhhssi var bigger than hhssi

# inactive delta
idx_inact <- delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20), ][idx_inact, k],
    delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20), ][idx_inact, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Outcome", k),
    xlim = range(c(delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20),][idx_inact, k],
                   delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20),][idx_inact, k])),
    ylim = range(c(delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20),][idx_inact, k],
                   delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20),][idx_inact, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20), ][!idx_inact, k],
    delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20), ][!idx_inact, k],
    xlab = "HHSSI delta sd",
    ylab = "NHHSSI delta sd",
    main = paste("Outcome", k),
    xlim = range(c(delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20),][!idx_inact, k],
                   delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20),][!idx_inact, k])),
    ylim = range(c(delta_sd_unif23_simpred_hhssi1[-c(1:3, 18:20),][!idx_inact, k],
                   delta_sd_unif23_simpred_nhhssi1[-c(1:3, 18:20),][!idx_inact, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





##########################################################################################
# Delta Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_unif23_simpred_hhssi1[, k],
    delta_est_unif23_simpred_nhhssi1[, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    col = point_col,
    main = paste("Overall estimates: Outcome", k),
    xlim = range(c(delta_est_unif23_simpred_hhssi1[,k],
                   delta_est_unif23_simpred_nhhssi1[,k])),
    ylim = range(c(delta_est_unif23_simpred_hhssi1[,k],
                   delta_est_unif23_simpred_nhhssi1[,k])))
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
    delta_est_unif23_simpred_hhssi1[c(1:3, 18:20), k],
    delta_est_unif23_simpred_nhhssi1[c(1:3, 18:20), k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("active: Outcome", k),
    xlim = range(c(delta_est_unif23_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_est_unif23_simpred_nhhssi1[c(1:3, 18:20),k])),
    ylim = range(c(delta_est_unif23_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_est_unif23_simpred_nhhssi1[c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20), k],
    delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20), k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    col = point_col_inct,
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20),k])),
    ylim = range(c(delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20),k])))
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
    delta_est_unif23_simpred_hhssi1[idx, k],
    delta_est_unif23_simpred_nhhssi1[idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_est_unif23_simpred_hhssi1[idx,k],
                   delta_est_unif23_simpred_nhhssi1[idx,k])),
    ylim = range(c(delta_est_unif23_simpred_hhssi1[idx,k],
                   delta_est_unif23_simpred_nhhssi1[idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_unif23_simpred_hhssi1[!idx, k],
    delta_est_unif23_simpred_nhhssi1[!idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_est_unif23_simpred_hhssi1[!idx,k],
                   delta_est_unif23_simpred_nhhssi1[!idx,k])),
    ylim = range(c(delta_est_unif23_simpred_hhssi1[!idx,k],
                   delta_est_unif23_simpred_nhhssi1[!idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





idx <- delta_var_mean_simpred_hhssi1[c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)]
# no active delta with nhhssi var bigger than hhssi

# inactive delta
idx <- delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20), ][idx, k],
    delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20), ][idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])),
    ylim = range(c(delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20), ][!idx, k],
    delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20), ][!idx, k],
    xlab = "HHSSI delta estimates",
    ylab = "NHHSSI delta estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])),
    ylim = range(c(delta_est_unif23_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_est_unif23_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}





##########################################################################################
# Delta Bias
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_unif23_simpred_hhssi1[, k],
    delta_bias_unif23_simpred_nhhssi1[, k],
    col = point_col,
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Overall bias: Outcome", k),
    xlim = range(c(delta_bias_unif23_simpred_hhssi1[,k],
                   delta_bias_unif23_simpred_nhhssi1[,k])),
    ylim = range(c(delta_bias_unif23_simpred_hhssi1[,k],
                   delta_bias_unif23_simpred_nhhssi1[,k])))
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
    delta_bias_unif23_simpred_hhssi1[c(1:3, 18:20), k],
    delta_bias_unif23_simpred_nhhssi1[c(1:3, 18:20), k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("active: Outcome", k),
    xlim = range(c(delta_bias_unif23_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_bias_unif23_simpred_nhhssi1[c(1:3, 18:20),k])),
    ylim = range(c(delta_bias_unif23_simpred_hhssi1[c(1:3, 18:20),k],
                   delta_bias_unif23_simpred_nhhssi1[c(1:3, 18:20),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive delta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20), k],
    delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20), k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    col = point_col_inct,
    main = paste("inactive: Outcome", k),
    xlim = range(c(delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20),k])),
    ylim = range(c(delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20),k],
                   delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20),k])))
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
    delta_bias_unif23_simpred_hhssi1[idx, k],
    delta_bias_unif23_simpred_nhhssi1[idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_bias_unif23_simpred_hhssi1[idx,k],
                   delta_bias_unif23_simpred_nhhssi1[idx,k])),
    ylim = range(c(delta_bias_unif23_simpred_hhssi1[idx,k],
                   delta_bias_unif23_simpred_nhhssi1[idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_unif23_simpred_hhssi1[!idx, k],
    delta_bias_unif23_simpred_nhhssi1[!idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(delta_bias_unif23_simpred_hhssi1[!idx,k],
                   delta_bias_unif23_simpred_nhhssi1[!idx,k])),
    ylim = range(c(delta_bias_unif23_simpred_hhssi1[!idx,k],
                   delta_bias_unif23_simpred_nhhssi1[!idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}






idx <- delta_var_mean_simpred_hhssi1[c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[c(1:3, 18:20)]
# no active delta with nhhssi var bigger than hhssi

# inactive delta
idx <- delta_var_mean_simpred_hhssi1[-c(1:3, 18:20)] < delta_var_mean_simpred_nhhssi1[-c(1:3, 18:20)]
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20), ][idx, k],
    delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20), ][idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Outcome", k),
    xlim = range(c(delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])),
    ylim = range(c(delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20),][idx, k],
                   delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20),][idx, k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}


par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20), ][!idx, k],
    delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20), ][!idx, k],
    xlab = "HHSSI delta bias",
    ylab = "NHHSSI delta bias",
    main = paste("Outcome", k),
    xlim = range(c(delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])),
    ylim = range(c(delta_bias_unif23_simpred_hhssi1[-c(1:3, 18:20),][!idx, k],
                   delta_bias_unif23_simpred_nhhssi1[-c(1:3, 18:20),][!idx, k])))
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

unique_beta_idx <- unique(beta_idx) # they are all inactive beta
unique_gamma_idx <- unique(gamma_idx) # they are all inactive gamma

length(unique_beta_idx)
length(unique_gamma_idx)
###########################
diag_dat <- data.frame(
  delta = delta_idx,
  beta  = beta_idx,
  gamma = gamma_idx
)
beta_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$theta_update[, 2:45, ]
gamma_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$theta_update[, 46:62, ]
delta_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$theta_update[, 63:810, ]


beta_mean_hhssi <- apply(beta_unif23_simpred_hhssi1, c(2, 3), mean)
gamma_mean_hhssi <- apply(gamma_unif23_simpred_hhssi1, c(2, 3), mean)
delta_mean_hhssi <- apply(delta_unif23_simpred_hhssi1, c(2, 3), mean)

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
beta_sd_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$theta_update[, 2:45, ], c(2,3), sd)
beta_sd_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$theta_update[, 2:45, ], c(2,3), sd)

# Posterior mean estimate of beta
beta_est_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$theta_update[, 2:45, ], c(2,3), mean)
beta_est_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$theta_update[, 2:45, ], c(2,3), mean)

# Bias for each outcome
beta_bias_unif23_simpred_hhssi1 <- (beta_est_unif23_simpred_hhssi1 - sim_data_unif23_simpred1$beta_true)
beta_bias_unif23_simpred_nhhssi1 <- (beta_est_unif23_simpred_nhhssi1 - sim_data_unif23_simpred1$beta_true)

# prior var estimate of beta
beta_var_unif23_simpred_hhssi1 <- apply((res_unif23_simpred_hhssi1$lambdasq_beta_update*
                                            res_unif23_simpred_hhssi1$tausq_beta_update), 2, mean)
beta_var_unif23_simpred_nhhssi1 <- apply((res_unif23_simpred_nhhssi1$lambdasq_beta_update*
                                              res_unif23_simpred_nhhssi1$tausq_beta_update), 2, mean)


# beta Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_unif23_simpred_hhssi1[, k],
    beta_est_unif23_simpred_nhhssi1[, k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("Overall estimates: Outcome", k),
    xlim = range(c(beta_est_unif23_simpred_hhssi1[,k],
                   beta_est_unif23_simpred_nhhssi1[,k])),
    ylim = range(c(beta_est_unif23_simpred_hhssi1[,k],
                   beta_est_unif23_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_unif23_simpred_hhssi1[1:2, k],
    beta_est_unif23_simpred_nhhssi1[1:2, k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("active: Outcome", k),
    xlim = range(c(beta_est_unif23_simpred_hhssi1[1:2,k],
                   beta_est_unif23_simpred_nhhssi1[1:2,k])),
    ylim = range(c(beta_est_unif23_simpred_hhssi1[1:2,k],
                   beta_est_unif23_simpred_nhhssi1[1:2,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_unif23_simpred_hhssi1[-c(1:2), k],
    beta_est_unif23_simpred_nhhssi1[-c(1:2), k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(beta_est_unif23_simpred_hhssi1[-c(1:2),k],
                   beta_est_unif23_simpred_nhhssi1[-c(1:2),k])),
    ylim = range(c(beta_est_unif23_simpred_hhssi1[-c(1:2),k],
                   beta_est_unif23_simpred_nhhssi1[-c(1:2),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}







# all the indices are all inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_est_unif23_simpred_hhssi1[unique_beta_idx, k],
    beta_est_unif23_simpred_nhhssi1[unique_beta_idx, k],
    xlab = "HHSSI beta estimates",
    ylab = "NHHSSI beta estimates",
    main = paste("Outcome", k),
    xlim = range(c(beta_est_unif23_simpred_hhssi1[unique_beta_idx,k],
                   beta_est_unif23_simpred_nhhssi1[unique_beta_idx,k])),
    ylim = range(c(beta_est_unif23_simpred_hhssi1[unique_beta_idx,k],
                   beta_est_unif23_simpred_nhhssi1[unique_beta_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




##########################################################################################
# beta Bias
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_unif23_simpred_hhssi1[, k],
    beta_bias_unif23_simpred_nhhssi1[, k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("Overall bias: Outcome", k),
    xlim = range(c(beta_bias_unif23_simpred_hhssi1[,k],
                   beta_bias_unif23_simpred_nhhssi1[,k])),
    ylim = range(c(beta_bias_unif23_simpred_hhssi1[,k],
                   beta_bias_unif23_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_unif23_simpred_hhssi1[c(1:2), k],
    beta_bias_unif23_simpred_nhhssi1[c(1:2), k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("active: Outcome", k),
    xlim = range(c(beta_bias_unif23_simpred_hhssi1[c(1:2),k],
                   beta_bias_unif23_simpred_nhhssi1[c(1:2),k])),
    ylim = range(c(beta_bias_unif23_simpred_hhssi1[c(1:2),k],
                   beta_bias_unif23_simpred_nhhssi1[c(1:2),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_unif23_simpred_hhssi1[-c(1:2), k],
    beta_bias_unif23_simpred_nhhssi1[-c(1:2), k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("inactive: Outcome", k),
    xlim = range(c(beta_bias_unif23_simpred_hhssi1[-c(1:2),k],
                   beta_bias_unif23_simpred_nhhssi1[-c(1:2),k])),
    ylim = range(c(beta_bias_unif23_simpred_hhssi1[-c(1:2),k],
                   beta_bias_unif23_simpred_nhhssi1[-c(1:2),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}



# all the indices are all inactive beta
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    beta_bias_unif23_simpred_hhssi1[unique_beta_idx, k],
    beta_bias_unif23_simpred_nhhssi1[unique_beta_idx, k],
    xlab = "HHSSI beta bias",
    ylab = "NHHSSI beta bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(beta_bias_unif23_simpred_hhssi1[unique_beta_idx,k],
                   beta_bias_unif23_simpred_nhhssi1[unique_beta_idx,k])),
    ylim = range(c(beta_bias_unif23_simpred_hhssi1[unique_beta_idx,k],
                   beta_bias_unif23_simpred_nhhssi1[unique_beta_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




######################################################################################

# beta prior variance estimate
# overall
par(mfrow = c(1,1))
plot(beta_var_unif23_simpred_hhssi1,
    beta_var_unif23_simpred_nhhssi1,
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_unif23_simpred_hhssi1,
                   beta_var_unif23_simpred_nhhssi1)),
    ylim = range(c(beta_var_unif23_simpred_hhssi1,
                   beta_var_unif23_simpred_nhhssi1)))
  abline(0,1, lty = 2, lwd = 1, col = "red")


# active beta
plot(beta_var_unif23_simpred_hhssi1[c(1:2)],
    beta_var_unif23_simpred_nhhssi1[c(1:2)],
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_unif23_simpred_hhssi1[c(1:2)],
                   beta_var_unif23_simpred_nhhssi1[c(1:2)])),
    ylim = range(c(beta_var_unif23_simpred_hhssi1[c(1:2)],
                   beta_var_unif23_simpred_nhhssi1[c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive beta
plot(beta_var_unif23_simpred_hhssi1[-c(1:2)],
    beta_var_unif23_simpred_nhhssi1[-c(1:2)],
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_unif23_simpred_hhssi1[-c(1:2)],
                   beta_var_unif23_simpred_nhhssi1[-c(1:2)])),
    ylim = range(c(beta_var_unif23_simpred_hhssi1[-c(1:2)],
                   beta_var_unif23_simpred_nhhssi1[-c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive beta
plot(beta_var_unif23_simpred_hhssi1[unique_beta_idx],
    beta_var_unif23_simpred_nhhssi1[unique_beta_idx],
    xlab = "HHSSI beta var",
    ylab = "NHHSSI beta var",
    xlim = range(c(beta_var_unif23_simpred_hhssi1[unique_beta_idx],
                   beta_var_unif23_simpred_nhhssi1[unique_beta_idx])),
    ylim = range(c(beta_var_unif23_simpred_hhssi1[unique_beta_idx],
                   beta_var_unif23_simpred_nhhssi1[unique_beta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")





######################################################################################
# Posterior sd estimate of gamma
gamma_sd_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$theta_update[, 46:62, ], c(2,3), sd)
gamma_sd_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$theta_update[, 46:62, ], c(2,3), sd)

# Posterior mean estimate of gamma
gamma_est_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$theta_update[, 46:62, ], c(2,3), mean)
gamma_est_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$theta_update[, 46:62, ], c(2,3), mean)

# Bias for each outcome
gamma_bias_unif23_simpred_hhssi1 <- (gamma_est_unif23_simpred_hhssi1 - sim_data_unif23_simpred1$gamma_true)
gamma_bias_unif23_simpred_nhhssi1 <- (gamma_est_unif23_simpred_nhhssi1 - sim_data_unif23_simpred1$gamma_true)

# prior var estimate of gamma
gamma_var_unif23_simpred_hhssi1 <- apply((res_unif23_simpred_hhssi1$lambdasq_gamma_update*
                                             res_unif23_simpred_hhssi1$tausq_gamma_update), 2, mean)
gamma_var_unif23_simpred_nhhssi1 <- apply((res_unif23_simpred_nhhssi1$lambdasq_gamma_update*
                                              res_unif23_simpred_nhhssi1$tausq_gamma_update), 2, mean)


# gamma Estimates
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_unif23_simpred_hhssi1[, k],
    gamma_est_unif23_simpred_nhhssi1[, k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("Overall estimates: Outcome", k),
    xlim = range(c(gamma_est_unif23_simpred_hhssi1[,k],
                   gamma_est_unif23_simpred_nhhssi1[,k])),
    ylim = range(c(gamma_est_unif23_simpred_hhssi1[,k],
                   gamma_est_unif23_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_unif23_simpred_hhssi1[1:3, k],
    gamma_est_unif23_simpred_nhhssi1[1:3, k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("active: Outcome", k),
    xlim = range(c(gamma_est_unif23_simpred_hhssi1[1:3,k],
                   gamma_est_unif23_simpred_nhhssi1[1:3,k])),
    ylim = range(c(gamma_est_unif23_simpred_hhssi1[1:3,k],
                   gamma_est_unif23_simpred_nhhssi1[1:3,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_unif23_simpred_hhssi1[-c(1:3), k],
    gamma_est_unif23_simpred_nhhssi1[-c(1:3), k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("inactive: Outcome", k),
    xlim = range(c(gamma_est_unif23_simpred_hhssi1[-c(1:3),k],
                   gamma_est_unif23_simpred_nhhssi1[-c(1:3),k])),
    ylim = range(c(gamma_est_unif23_simpred_hhssi1[-c(1:3),k],
                   gamma_est_unif23_simpred_nhhssi1[-c(1:3),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}







# overall (infact all are inactive gamma)
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_est_unif23_simpred_hhssi1[unique_gamma_idx, k],
    gamma_est_unif23_simpred_nhhssi1[unique_gamma_idx, k],
    xlab = "HHSSI gamma estimates",
    ylab = "NHHSSI gamma estimates",
    main = paste("Outcome", k),
    xlim = range(c(gamma_est_unif23_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_est_unif23_simpred_nhhssi1[unique_gamma_idx,k])),
    ylim = range(c(gamma_est_unif23_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_est_unif23_simpred_nhhssi1[unique_gamma_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




##########################################################################################
# gamma Bias
# overall
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_unif23_simpred_hhssi1[, k],
    gamma_bias_unif23_simpred_nhhssi1[, k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("Overall bias: Outcome", k),
    xlim = range(c(gamma_bias_unif23_simpred_hhssi1[,k],
                   gamma_bias_unif23_simpred_nhhssi1[,k])),
    ylim = range(c(gamma_bias_unif23_simpred_hhssi1[,k],
                   gamma_bias_unif23_simpred_nhhssi1[,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# active gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_unif23_simpred_hhssi1[c(1:3), k],
    gamma_bias_unif23_simpred_nhhssi1[c(1:3), k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("active: Outcome", k),
    xlim = range(c(gamma_bias_unif23_simpred_hhssi1[c(1:3),k],
                   gamma_bias_unif23_simpred_nhhssi1[c(1:3),k])),
    ylim = range(c(gamma_bias_unif23_simpred_hhssi1[c(1:3),k],
                   gamma_bias_unif23_simpred_nhhssi1[c(1:3),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}

# inactive gamma
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_unif23_simpred_hhssi1[-c(1:3), k],
    gamma_bias_unif23_simpred_nhhssi1[-c(1:3), k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("inactive: Outcome", k),
    xlim = range(c(gamma_bias_unif23_simpred_hhssi1[-c(1:3),k],
                   gamma_bias_unif23_simpred_nhhssi1[-c(1:3),k])),
    ylim = range(c(gamma_bias_unif23_simpred_hhssi1[-c(1:3),k],
                   gamma_bias_unif23_simpred_nhhssi1[-c(1:3),k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}



# overall (infact all are inactive gamma)
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(
    gamma_bias_unif23_simpred_hhssi1[unique_gamma_idx, k],
    gamma_bias_unif23_simpred_nhhssi1[unique_gamma_idx, k],
    xlab = "HHSSI gamma bias",
    ylab = "NHHSSI gamma bias",
    main = paste("Overall: Outcome", k),
    xlim = range(c(gamma_bias_unif23_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_bias_unif23_simpred_nhhssi1[unique_gamma_idx,k])),
    ylim = range(c(gamma_bias_unif23_simpred_hhssi1[unique_gamma_idx,k],
                   gamma_bias_unif23_simpred_nhhssi1[unique_gamma_idx,k])))
  abline(0,1, lty = 2, lwd = 1, col = "red")
}




######################################################################################
# gamma prior variance estimate
# overall
par(mfrow = c(1,1))
plot(gamma_var_unif23_simpred_hhssi1,
     gamma_var_unif23_simpred_nhhssi1,
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_unif23_simpred_hhssi1,
                    gamma_var_unif23_simpred_nhhssi1)),
     ylim = range(c(gamma_var_unif23_simpred_hhssi1,
                    gamma_var_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active gamma
plot(gamma_var_unif23_simpred_hhssi1[c(1:3)],
     gamma_var_unif23_simpred_nhhssi1[c(1:3)],
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_unif23_simpred_hhssi1[c(1:3)],
                    gamma_var_unif23_simpred_nhhssi1[c(1:3)])),
     ylim = range(c(gamma_var_unif23_simpred_hhssi1[c(1:3)],
                    gamma_var_unif23_simpred_nhhssi1[c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive gamma
plot(gamma_var_unif23_simpred_hhssi1[-c(1:3)],
     gamma_var_unif23_simpred_nhhssi1[-c(1:3)],
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_unif23_simpred_hhssi1[-c(1:3)],
                    gamma_var_unif23_simpred_nhhssi1[-c(1:3)])),
     ylim = range(c(gamma_var_unif23_simpred_hhssi1[-c(1:3)],
                    gamma_var_unif23_simpred_nhhssi1[-c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# overall (infact all are inactive gamma)
plot(gamma_var_unif23_simpred_hhssi1[unique_gamma_idx],
     gamma_var_unif23_simpred_nhhssi1[unique_gamma_idx],
     xlab = "HHSSI gamma var",
     ylab = "NHHSSI gamma var",
     xlim = range(c(gamma_var_unif23_simpred_hhssi1[unique_gamma_idx],
                    gamma_var_unif23_simpred_nhhssi1[unique_gamma_idx])),
     ylim = range(c(gamma_var_unif23_simpred_hhssi1[unique_gamma_idx],
                    gamma_var_unif23_simpred_nhhssi1[unique_gamma_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")

######################################################################################

# lambdasq terms
lambdasq_beta_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$lambdasq_beta_update, 2, mean)
lambdasq_beta_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$lambdasq_beta_update, 2, mean)

lambdasq_gamma_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$lambdasq_gamma_update, 2, mean)
lambdasq_gamma_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$lambdasq_gamma_update, 2, mean)

lambdasq_delta_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$lambdasq_delta_update, 2, mean)
lambdasq_delta_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$lambdasq_delta_update, 2, mean)

######################################################################################

# lambdasq_beta estimate
# overall
par(mfrow = c(1,1))
plot(lambdasq_beta_unif23_simpred_hhssi1,
     lambdasq_beta_unif23_simpred_nhhssi1,
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_unif23_simpred_hhssi1,
                    lambdasq_beta_unif23_simpred_nhhssi1)),
     ylim = range(c(lambdasq_beta_unif23_simpred_hhssi1,
                    lambdasq_beta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active beta
plot(lambdasq_beta_unif23_simpred_hhssi1[c(1:2)],
     lambdasq_beta_unif23_simpred_nhhssi1[c(1:2)],
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_unif23_simpred_hhssi1[c(1:2)],
                    lambdasq_beta_unif23_simpred_nhhssi1[c(1:2)])),
     ylim = range(c(lambdasq_beta_unif23_simpred_hhssi1[c(1:2)],
                    lambdasq_beta_unif23_simpred_nhhssi1[c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive beta
plot(lambdasq_beta_unif23_simpred_hhssi1[-c(1:2)],
     lambdasq_beta_unif23_simpred_nhhssi1[-c(1:2)],
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_unif23_simpred_hhssi1[-c(1:2)],
                    lambdasq_beta_unif23_simpred_nhhssi1[-c(1:2)])),
     ylim = range(c(lambdasq_beta_unif23_simpred_hhssi1[-c(1:2)],
                    lambdasq_beta_unif23_simpred_nhhssi1[-c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive beta
plot(lambdasq_beta_unif23_simpred_hhssi1[unique_beta_idx],
     lambdasq_beta_unif23_simpred_nhhssi1[unique_beta_idx],
     xlab = "HHSSI lambdasq beta",
     ylab = "NHHSSI lambdasq beta",
     xlim = range(c(lambdasq_beta_unif23_simpred_hhssi1[unique_beta_idx],
                    lambdasq_beta_unif23_simpred_nhhssi1[unique_beta_idx])),
     ylim = range(c(lambdasq_beta_unif23_simpred_hhssi1[unique_beta_idx],
                    lambdasq_beta_unif23_simpred_nhhssi1[unique_beta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# lambdasq_gamma estimate
# overall
par(mfrow = c(1,1))
plot(lambdasq_gamma_unif23_simpred_hhssi1,
     lambdasq_gamma_unif23_simpred_nhhssi1,
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_unif23_simpred_hhssi1,
                    lambdasq_gamma_unif23_simpred_nhhssi1)),
     ylim = range(c(lambdasq_gamma_unif23_simpred_hhssi1,
                    lambdasq_gamma_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active gamma
plot(lambdasq_gamma_unif23_simpred_hhssi1[c(1:3)],
     lambdasq_gamma_unif23_simpred_nhhssi1[c(1:3)],
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_unif23_simpred_hhssi1[c(1:3)],
                    lambdasq_gamma_unif23_simpred_nhhssi1[c(1:3)])),
     ylim = range(c(lambdasq_gamma_unif23_simpred_hhssi1[c(1:3)],
                    lambdasq_gamma_unif23_simpred_nhhssi1[c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive gamma
plot(lambdasq_gamma_unif23_simpred_hhssi1[-c(1:3)],
     lambdasq_gamma_unif23_simpred_nhhssi1[-c(1:3)],
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_unif23_simpred_hhssi1[-c(1:3)],
                    lambdasq_gamma_unif23_simpred_nhhssi1[-c(1:3)])),
     ylim = range(c(lambdasq_gamma_unif23_simpred_hhssi1[-c(1:3)],
                    lambdasq_gamma_unif23_simpred_nhhssi1[-c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive gamma
plot(lambdasq_gamma_unif23_simpred_hhssi1[unique_gamma_idx],
     lambdasq_gamma_unif23_simpred_nhhssi1[unique_gamma_idx],
     xlab = "HHSSI lambdasq gamma",
     ylab = "NHHSSI lambdasq gamma",
     xlim = range(c(lambdasq_gamma_unif23_simpred_hhssi1[unique_gamma_idx],
                    lambdasq_gamma_unif23_simpred_nhhssi1[unique_gamma_idx])),
     ylim = range(c(lambdasq_gamma_unif23_simpred_hhssi1[unique_gamma_idx],
                    lambdasq_gamma_unif23_simpred_nhhssi1[unique_gamma_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# lambdasq_delta estimate
# overall
par(mfrow = c(1,1))
plot(lambdasq_delta_unif23_simpred_hhssi1,
     lambdasq_delta_unif23_simpred_nhhssi1,
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_unif23_simpred_hhssi1,
                    lambdasq_delta_unif23_simpred_nhhssi1)),
     ylim = range(c(lambdasq_delta_unif23_simpred_hhssi1,
                    lambdasq_delta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active delta
plot(lambdasq_delta_unif23_simpred_hhssi1[c(1:3, 18:20)],
     lambdasq_delta_unif23_simpred_nhhssi1[c(1:3, 18:20)],
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_unif23_simpred_hhssi1[c(1:3, 18:20)],
                    lambdasq_delta_unif23_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_unif23_simpred_hhssi1[c(1:3, 18:20)],
                    lambdasq_delta_unif23_simpred_nhhssi1[c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive delta
plot(lambdasq_delta_unif23_simpred_hhssi1[-c(1:3, 18:20)],
     lambdasq_delta_unif23_simpred_nhhssi1[-c(1:3, 18:20)],
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_unif23_simpred_hhssi1[-c(1:3, 18:20)],
                    lambdasq_delta_unif23_simpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(lambdasq_delta_unif23_simpred_hhssi1[-c(1:3, 18:20)],
                    lambdasq_delta_unif23_simpred_nhhssi1[-c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive delta
plot(lambdasq_delta_unif23_simpred_hhssi1[delta_idx],
     lambdasq_delta_unif23_simpred_nhhssi1[delta_idx],
     xlab = "HHSSI lambdasq delta",
     ylab = "NHHSSI lambdasq delta",
     xlim = range(c(lambdasq_delta_unif23_simpred_hhssi1[delta_idx],
                    lambdasq_delta_unif23_simpred_nhhssi1[delta_idx])),
     ylim = range(c(lambdasq_delta_unif23_simpred_hhssi1[delta_idx],
                    lambdasq_delta_unif23_simpred_nhhssi1[delta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")









######################################################################################

# psi terms
psi_beta_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$psi_beta_update, 2, mean)
psi_beta_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$psi_beta_update, 2, mean)

psi_gamma_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$psi_gamma_update, 2, mean)
psi_gamma_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$psi_gamma_update, 2, mean)

psi_delta_unif23_simpred_hhssi1 <- apply(res_unif23_simpred_hhssi1$psi_delta_update, 2, mean)
psi_delta_unif23_simpred_nhhssi1 <- apply(res_unif23_simpred_nhhssi1$psi_delta_update, 2, mean)

######################################################################################

# psi_beta estimate
# overall
par(mfrow = c(1,1))
plot(psi_beta_unif23_simpred_hhssi1,
     psi_beta_unif23_simpred_nhhssi1,
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_unif23_simpred_hhssi1,
                    psi_beta_unif23_simpred_nhhssi1)),
     ylim = range(c(psi_beta_unif23_simpred_hhssi1,
                    psi_beta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active beta
plot(psi_beta_unif23_simpred_hhssi1[c(1:2)],
     psi_beta_unif23_simpred_nhhssi1[c(1:2)],
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_unif23_simpred_hhssi1[c(1:2)],
                    psi_beta_unif23_simpred_nhhssi1[c(1:2)])),
     ylim = range(c(psi_beta_unif23_simpred_hhssi1[c(1:2)],
                    psi_beta_unif23_simpred_nhhssi1[c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive beta
plot(psi_beta_unif23_simpred_hhssi1[-c(1:2)],
     psi_beta_unif23_simpred_nhhssi1[-c(1:2)],
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_unif23_simpred_hhssi1[-c(1:2)],
                    psi_beta_unif23_simpred_nhhssi1[-c(1:2)])),
     ylim = range(c(psi_beta_unif23_simpred_hhssi1[-c(1:2)],
                    psi_beta_unif23_simpred_nhhssi1[-c(1:2)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive beta
plot(psi_beta_unif23_simpred_hhssi1[unique_beta_idx],
     psi_beta_unif23_simpred_nhhssi1[unique_beta_idx],
     xlab = "HHSSI psi beta",
     ylab = "NHHSSI psi beta",
     xlim = range(c(psi_beta_unif23_simpred_hhssi1[unique_beta_idx],
                    psi_beta_unif23_simpred_nhhssi1[unique_beta_idx])),
     ylim = range(c(psi_beta_unif23_simpred_hhssi1[unique_beta_idx],
                    psi_beta_unif23_simpred_nhhssi1[unique_beta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# psi_gamma estimate
# overall
par(mfrow = c(1,1))
plot(psi_gamma_unif23_simpred_hhssi1,
     psi_gamma_unif23_simpred_nhhssi1,
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_unif23_simpred_hhssi1,
                    psi_gamma_unif23_simpred_nhhssi1)),
     ylim = range(c(psi_gamma_unif23_simpred_hhssi1,
                    psi_gamma_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active gamma
plot(psi_gamma_unif23_simpred_hhssi1[c(1:3)],
     psi_gamma_unif23_simpred_nhhssi1[c(1:3)],
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_unif23_simpred_hhssi1[c(1:3)],
                    psi_gamma_unif23_simpred_nhhssi1[c(1:3)])),
     ylim = range(c(psi_gamma_unif23_simpred_hhssi1[c(1:3)],
                    psi_gamma_unif23_simpred_nhhssi1[c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive gamma
plot(psi_gamma_unif23_simpred_hhssi1[-c(1:3)],
     psi_gamma_unif23_simpred_nhhssi1[-c(1:3)],
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_unif23_simpred_hhssi1[-c(1:3)],
                    psi_gamma_unif23_simpred_nhhssi1[-c(1:3)])),
     ylim = range(c(psi_gamma_unif23_simpred_hhssi1[-c(1:3)],
                    psi_gamma_unif23_simpred_nhhssi1[-c(1:3)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive gamma
plot(psi_gamma_unif23_simpred_hhssi1[unique_gamma_idx],
     psi_gamma_unif23_simpred_nhhssi1[unique_gamma_idx],
     xlab = "HHSSI psi gamma",
     ylab = "NHHSSI psi gamma",
     xlim = range(c(psi_gamma_unif23_simpred_hhssi1[unique_gamma_idx],
                    psi_gamma_unif23_simpred_nhhssi1[unique_gamma_idx])),
     ylim = range(c(psi_gamma_unif23_simpred_hhssi1[unique_gamma_idx],
                    psi_gamma_unif23_simpred_nhhssi1[unique_gamma_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")




######################################################################################

# psi_delta estimate
# overall
par(mfrow = c(1,1))
plot(psi_delta_unif23_simpred_hhssi1,
     psi_delta_unif23_simpred_nhhssi1,
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_unif23_simpred_hhssi1,
                    psi_delta_unif23_simpred_nhhssi1)),
     ylim = range(c(psi_delta_unif23_simpred_hhssi1,
                    psi_delta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# active delta
plot(psi_delta_unif23_simpred_hhssi1[c(1:3, 18:20)],
     psi_delta_unif23_simpred_nhhssi1[c(1:3, 18:20)],
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_unif23_simpred_hhssi1[c(1:3, 18:20)],
                    psi_delta_unif23_simpred_nhhssi1[c(1:3, 18:20)])),
     ylim = range(c(psi_delta_unif23_simpred_hhssi1[c(1:3, 18:20)],
                    psi_delta_unif23_simpred_nhhssi1[c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")

# inactive delta
plot(psi_delta_unif23_simpred_hhssi1[-c(1:3, 18:20)],
     psi_delta_unif23_simpred_nhhssi1[-c(1:3, 18:20)],
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_unif23_simpred_hhssi1[-c(1:3, 18:20)],
                    psi_delta_unif23_simpred_nhhssi1[-c(1:3, 18:20)])),
     ylim = range(c(psi_delta_unif23_simpred_hhssi1[-c(1:3, 18:20)],
                    psi_delta_unif23_simpred_nhhssi1[-c(1:3, 18:20)])))
abline(0,1, lty = 2, lwd = 1, col = "red")



# all are infact inactive delta
plot(psi_delta_unif23_simpred_hhssi1[delta_idx],
     psi_delta_unif23_simpred_nhhssi1[delta_idx],
     xlab = "HHSSI psi delta",
     ylab = "NHHSSI psi delta",
     xlim = range(c(psi_delta_unif23_simpred_hhssi1[delta_idx],
                    psi_delta_unif23_simpred_nhhssi1[delta_idx])),
     ylim = range(c(psi_delta_unif23_simpred_hhssi1[delta_idx],
                    psi_delta_unif23_simpred_nhhssi1[delta_idx])))
abline(0,1, lty = 2, lwd = 1, col = "red")


######################################################################################

# tausq terms
tausq_beta_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_beta_update
tausq_beta_unif23_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_beta_update

tausq_gamma_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_gamma_update
tausq_gamma_unif23_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_gamma_update

tausq_delta_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$tausq_delta_update
tausq_delta_unif23_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$tausq_delta_update

######################################################################################

# tausq_beta estimate
plot(tausq_beta_unif23_simpred_hhssi1,
     tausq_beta_unif23_simpred_nhhssi1,
     xlab = "HHSSI tausq beta",
     ylab = "NHHSSI tausq beta",
     xlim = range(c(tausq_beta_unif23_simpred_hhssi1,
                    tausq_beta_unif23_simpred_nhhssi1)),
     ylim = range(c(tausq_beta_unif23_simpred_hhssi1,
                    tausq_beta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# tausq_gamma estimate
plot(tausq_gamma_unif23_simpred_hhssi1,
     tausq_gamma_unif23_simpred_nhhssi1,
     xlab = "HHSSI tausq gamma",
     ylab = "NHHSSI tausq gamma",
     xlim = range(c(tausq_gamma_unif23_simpred_hhssi1,
                    tausq_gamma_unif23_simpred_nhhssi1)),
     ylim = range(c(tausq_gamma_unif23_simpred_hhssi1,
                    tausq_gamma_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")



# tausq_delta estimate
plot(tausq_delta_unif23_simpred_hhssi1,
     tausq_delta_unif23_simpred_nhhssi1,
     xlab = "HHSSI tausq delta",
     ylab = "NHHSSI tausq delta",
     xlim = range(c(tausq_delta_unif23_simpred_hhssi1,
                    tausq_delta_unif23_simpred_nhhssi1)),
     ylim = range(c(tausq_delta_unif23_simpred_hhssi1,
                    tausq_delta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")





######################################################################################

# xi terms
xi_beta_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$xi_beta_update
xi_beta_unif23_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$xi_beta_update

xi_gamma_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$xi_gamma_update
xi_gamma_unif23_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$xi_gamma_update

xi_delta_unif23_simpred_hhssi1 <- res_unif23_simpred_hhssi1$xi_delta_update
xi_delta_unif23_simpred_nhhssi1 <- res_unif23_simpred_nhhssi1$xi_delta_update

######################################################################################

# xi_beta estimate
plot(xi_beta_unif23_simpred_hhssi1,
     xi_beta_unif23_simpred_nhhssi1,
     xlab = "HHSSI xi beta",
     ylab = "NHHSSI xi beta",
     xlim = range(c(xi_beta_unif23_simpred_hhssi1,
                    xi_beta_unif23_simpred_nhhssi1)),
     ylim = range(c(xi_beta_unif23_simpred_hhssi1,
                    xi_beta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")


# xi_gamma estimate
plot(xi_gamma_unif23_simpred_hhssi1,
     xi_gamma_unif23_simpred_nhhssi1,
     xlab = "HHSSI xi gamma",
     ylab = "NHHSSI xi gamma",
     xlim = range(c(xi_gamma_unif23_simpred_hhssi1,
                    xi_gamma_unif23_simpred_nhhssi1)),
     ylim = range(c(xi_gamma_unif23_simpred_hhssi1,
                    xi_gamma_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")



# xi_delta estimate
plot(xi_delta_unif23_simpred_hhssi1,
     xi_delta_unif23_simpred_nhhssi1,
     xlab = "HHSSI xi delta",
     ylab = "NHHSSI xi delta",
     xlim = range(c(xi_delta_unif23_simpred_hhssi1,
                    xi_delta_unif23_simpred_nhhssi1)),
     ylim = range(c(xi_delta_unif23_simpred_hhssi1,
                    xi_delta_unif23_simpred_nhhssi1)))
abline(0,1, lty = 2, lwd = 1, col = "red")



# bias of beta, gamma, delta on the same plot

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



dim(delta_est_unif23_simpred_hhssi1)
delta_est_unif23_simpred_nhhssi1
delta_true <- sim_data_unif23_simpred1$delta_true

rmse_hhssi <- sqrt(rowMeans((delta_est_unif23_simpred_hhssi1 - delta_true)^2))
rmse_nhhssi <- sqrt(rowMeans((delta_est_unif23_simpred_nhhssi1 - delta_true)^2 ))
active_rmse <- data.frame(
  index = active_idx,
  j = ((active_idx - 1) %/% M) + 1,
  m = ((active_idx - 1) %% M) + 1,
  HHSSI = rmse_hhssi[active_idx],
  NHHSSI = rmse_nhhssi[active_idx]
)
active_rmse




bias_hhssi <- rowMeans(delta_est_unif23_simpred_hhssi1 - delta_true)
bias_nhhssi <- rowMeans(delta_est_unif23_simpred_nhhssi1 - delta_true)
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



