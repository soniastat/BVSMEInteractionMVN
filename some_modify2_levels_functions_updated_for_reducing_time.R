
# some practice codes to reduce time

##########################################################################
# HHS model
#########################################################################

# Update lambdasq_beta (for all j = 1, 2, ..., J, and k = 1, 2, ..., K)
update_lambdasq_beta_BVS_HHS_MVN_cov <- function(J, K, beta_update, psi_beta_update, tausq_beta_update)
{
  lambdasq_beta_update_s <- matrix(NA, nrow = J, ncol = K)

  for(j in 1:J)
  {
    for(k in 1:K)
    {
      scale_lambdasq_beta <- ((1/psi_beta_update[j, k]) +
                                ((beta_update[j, k])^2 / (2 * tausq_beta_update[j])))
      lambdasq_beta_update_s[j, k] <- rinvgamma(1, shape = 1, scale = scale_lambdasq_beta) # IG distribution
    }
  }
  return(lambdasq_beta_update_s)
}

#######
update_lambdasq_beta_BVS_HHS_MVN_cov2 <- function(beta_update,
                                                  psi_beta_update,
                                                  tausq_beta_update)
{
  scale_lambdasq_beta <-
    1 / psi_beta_update +
    beta_update^2 / (2 * tausq_beta_update)

  matrix(
    rinvgamma(length(scale_lambdasq_beta),
              shape = 1,
              scale = c(t(scale_lambdasq_beta))),
    nrow = nrow(beta_update),
    byrow = TRUE
  )
}



J <- 5
K <- 4
beta_update <- matrix(rnorm(J*K), nrow = J, ncol = K)
psi_beta_update <- matrix(rgamma(J*K, 2, 3), nrow = J, ncol = K)
tausq_beta_update <- rgamma(J, 2, 3)


library(LaplacesDemon)

scale1 <- matrix(NA, J, K)
for(j in 1:J)
  for(k in 1:K)
    scale1[j,k] <-
  1/psi_beta_update[j,k] +
  beta_update[j,k]^2/(2*tausq_beta_update[j])

scale2 <-
  1/psi_beta_update +
  beta_update^2/(2*tausq_beta_update)

all.equal(scale1, scale2)




set.seed(1)
x1 <- rinvgamma(20, shape = 1, scale = 2)

set.seed(1)
x2 <- replicate(20, rinvgamma(1, shape = 1, scale = 2))

identical(x1, x2)
all.equal(x1, x2)

cbind(x1, x2)





set.seed(788)
x1 <- update_lambdasq_beta_BVS_HHS_MVN_cov(
  J, K,
  beta_update,
  psi_beta_update,
  tausq_beta_update
)

set.seed(788)
x2 <- update_lambdasq_beta_BVS_HHS_MVN_cov2(
  beta_update,
  psi_beta_update,
  tausq_beta_update
)

identical(x1, x2)
all.equal(x1, x2)


#######################################################################################

# Update lambdasq_gamma (for all m = 1, 2, ..., M, and k = 1, 2, ..., K)
update_lambdasq_gamma_BVS_HHS_MVN_cov <- function(M, K, gamma_update,
                                                  psi_gamma_update,
                                                  tausq_gamma_update)
{
  lambdasq_gamma_update_s <- matrix(NA, nrow = M, ncol = K)

  for(m in 1:M)
  {
    for(k in 1:K)
    {
      scale_lambdasq_gamma <- ((1/psi_gamma_update[m, k]) +
                                 ((gamma_update[m, k])^2 / (2 * tausq_gamma_update[m])))
      lambdasq_gamma_update_s[m, k] <- rinvgamma(1, shape = 1, scale = scale_lambdasq_gamma) # IG distribution
    }
  }
  return(lambdasq_gamma_update_s)
}


#######
update_lambdasq_gamma_BVS_HHS_MVN_cov2 <- function(gamma_update,
                                                  psi_gamma_update,
                                                  tausq_gamma_update)
{
  scale_lambdasq_gamma <-
    1 / psi_gamma_update +
    gamma_update^2 / (2 * tausq_gamma_update)

  matrix(
    rinvgamma(length(scale_lambdasq_gamma),
              shape = 1,
              scale = c(t(scale_lambdasq_gamma))),
    nrow = nrow(gamma_update),
    byrow = TRUE
  )
}



M <- 6
K <- 4
gamma_update <- matrix(rnorm(M*K), nrow = M, ncol = K)
psi_gamma_update <- matrix(rgamma(M*K, 2, 3), nrow = M, ncol = K)
tausq_gamma_update <- rgamma(M, 2, 3)


set.seed(788)
x1 <- update_lambdasq_gamma_BVS_HHS_MVN_cov(
  M, K,
  gamma_update,
  psi_gamma_update,
  tausq_gamma_update
)

set.seed(788)
x2 <- update_lambdasq_gamma_BVS_HHS_MVN_cov2(
  gamma_update,
  psi_gamma_update,
  tausq_gamma_update
)

identical(x1, x2)
all.equal(x1, x2)




#####################################################################################

# Update lambdasq_delta (for all j and m, and k = 1, 2, ..., K)
update_lambdasq_delta_BVS_HHS_MVN_cov <- function(J, M, K, delta_update, lambdasq_beta_update,
                                                  lambdasq_gamma_update, psi_delta_update,
                                                  tausq_delta_update)
{
  lambdasq_delta_update_s <- matrix(NA, nrow = J*M, ncol = K)

  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction
      for(k in 1:K)
      {
        scale_lambdasq_delta <- ((1/psi_delta_update[jm, k]) +
                                   ((delta_update[jm, k])^2 / (2 * lambdasq_beta_update[j, k]
                                                               * lambdasq_gamma_update[m, k] * tausq_delta_update[jm])))

        lambdasq_delta_update_s[jm, k] <- rinvgamma(1, shape = 1, scale = scale_lambdasq_delta) # IG distribution
      }
    }
  }
  return(lambdasq_delta_update_s)
}

#######

update_lambdasq_delta_BVS_HHS_MVN_cov2 <- function(J, M, K,
                                                   delta_update,
                                                   lambdasq_beta_update,
                                                   lambdasq_gamma_update,
                                                   psi_delta_update,
                                                   tausq_delta_update)
{
  beta_part <- lambdasq_beta_update[rep(1:J, each = M),
                                    ,
                                    drop = FALSE]

  gamma_part <- lambdasq_gamma_update[rep(1:M, times = J),
                                      ,
                                      drop = FALSE]

  scale_lambdasq_delta <-
    1 / psi_delta_update +
    delta_update^2 /
    (2 * beta_part * gamma_part *
       tausq_delta_update)

  matrix(
    rinvgamma(length(scale_lambdasq_delta),
              shape = 1,
              scale = as.vector(t(scale_lambdasq_delta))),
    nrow = J*M,
    byrow = TRUE
  )
}


delta_update <- matrix(rnorm(J*M*K), nrow = J*M, ncol = K)
psi_delta_update <- matrix(rgamma(J*M*K, 2, 3), nrow = J*M, ncol = K)
tausq_delta_update <- rgamma(J*M, 2, 3)
lambdasq_beta_update <- matrix(rgamma(J*K, 2, 3), nrow = J, ncol = K)
lambdasq_gamma_update <- matrix(rgamma(M*K, 2, 3), nrow = M, ncol = K)



set.seed(788)
x1 <- update_lambdasq_delta_BVS_HHS_MVN_cov(
  J, M, K, delta_update, lambdasq_beta_update,
  lambdasq_gamma_update, psi_delta_update,
  tausq_delta_update
)

set.seed(788)
x2 <- update_lambdasq_delta_BVS_HHS_MVN_cov2(
  J, M, K,
  delta_update,
  lambdasq_beta_update,
  lambdasq_gamma_update,
  psi_delta_update,
  tausq_delta_update
)

identical(x1, x2)
all.equal(x1, x2)


#############

library(microbenchmark)


set.seed(123)

benchmark_result <- microbenchmark(

  Original = update_lambdasq_delta_BVS_HHS_MVN_cov(
    J, M, K,
    delta_update,
    lambdasq_beta_update,
    lambdasq_gamma_update,
    psi_delta_update,
    tausq_delta_update
  ),

  Vectorized = update_lambdasq_delta_BVS_HHS_MVN_cov2(
    J, M, K,
    delta_update,
    lambdasq_beta_update,
    lambdasq_gamma_update,
    psi_delta_update,
    tausq_delta_update
  ),

  times = 1000,
  unit = "microseconds"
)

benchmark_result

##############################################

# Update tausq_delta (for all j, and m )
update_tausq_delta_BVS_HHS_MVN_cov <- function(J, M, K, delta_update, lambdasq_delta_update,
                                               lambdasq_beta_update, lambdasq_gamma_update, xi_delta_update)
{
  tausq_delta_update_s <- rep(NA, J * M)

  shape_tausq_delta <- (K + 1) / 2
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction

      sum_term_tausq_delta <- 0
      for (k in 1:K)
      {
        sum_term_tausq_delta <- sum_term_tausq_delta + ((delta_update[jm, k])^2 / (2 *
                                                                                     lambdasq_delta_update[jm, k] * lambdasq_beta_update[j, k] * lambdasq_gamma_update[m, k]))
      }

      scale_tausq_delta <- ((1/xi_delta_update) + sum_term_tausq_delta)
      tausq_delta_update_s[jm] <- rinvgamma(1, shape = shape_tausq_delta, scale = scale_tausq_delta) # IG distribution
    }
  }
  return(tausq_delta_update_s)
}



###
update_tausq_delta_BVS_HHS_MVN_cov2 <- function(J, M, K,
                                                delta_update,
                                                lambdasq_delta_update,
                                                lambdasq_beta_update,
                                                lambdasq_gamma_update,
                                                xi_delta_update)
{
  shape_tausq_delta <- (K + 1) / 2

  ## Expand beta and gamma according to jm = (j-1)M + m
  beta_part  <- lambdasq_beta_update[rep(seq_len(J), each = M),
                                     , drop = FALSE]

  gamma_part <- lambdasq_gamma_update[rep(seq_len(M), times = J),
                                      , drop = FALSE]

  ## Compute all K terms simultaneously
  sum_term_tausq_delta <- rowSums(
    delta_update^2 /
      (2 * lambdasq_delta_update *
         beta_part *
         gamma_part)
  )

  scale_tausq_delta <- 1 / xi_delta_update + sum_term_tausq_delta

  rinvgamma(J * M,
            shape = shape_tausq_delta,
            scale = scale_tausq_delta)
}




lambdasq_delta_update <- matrix(rgamma(J*M*K, 2, 3), nrow = J*M, ncol = K)
xi_delta_update <- 1.2


set.seed(788)
x1 <- update_tausq_delta_BVS_HHS_MVN_cov(
  J, M, K, delta_update,
  lambdasq_delta_update,
  lambdasq_beta_update,
  lambdasq_gamma_update,
  xi_delta_update
)

set.seed(788)
x2 <- update_tausq_delta_BVS_HHS_MVN_cov2(
  J, M, K,
  delta_update,
  lambdasq_delta_update,
  lambdasq_beta_update,
  lambdasq_gamma_update,
  xi_delta_update
)

identical(x1, x2)
all.equal(x1, x2)




######################################################################

# Update psi_beta (for all j = 1, 2, ..., J, and k = 1, 2, ..., K)
update_psi_beta_BVS_HHS_MVN_cov <- function(J, K, lambdasq_beta_update)
{
  psi_beta_update_s <- matrix(NA, nrow = J, ncol = K)

  for(j in 1:J)
  {
    for(k in 1:K)
    {
      scale_psi_beta <- (1 + (1/lambdasq_beta_update[j, k]))
      psi_beta_update_s[j, k] <- rinvgamma(1, shape = 1, scale = scale_psi_beta) # IG distribution
    }
  }
  return(psi_beta_update_s)
}

#########
update_psi_beta_BVS_HHS_MVN_cov2 <- function(lambdasq_beta_update)
{
  scale_psi_beta <- 1 + 1 / lambdasq_beta_update

  matrix(
    rinvgamma(length(scale_psi_beta),
              shape = 1,
              scale = as.vector(t(scale_psi_beta))),
    nrow = nrow(lambdasq_beta_update),
    byrow = TRUE
  )
}




set.seed(123)
x1 <- update_psi_beta_BVS_HHS_MVN_cov(
  J, K,
  lambdasq_beta_update
)

set.seed(123)
x2 <- update_psi_beta_BVS_HHS_MVN_cov2(
  lambdasq_beta_update
)

identical(x1, x2)
all.equal(x1, x2)



#######################################################################################

# Update psi_gamma (for all m = 1, 2, ..., M, and k = 1, 2, ..., K)
update_psi_gamma_BVS_HHS_MVN_cov <- function(M, K, lambdasq_gamma_update)
{
  psi_gamma_update_s <- matrix(NA, nrow = M, ncol = K)

  for(m in 1:M)
  {
    for(k in 1:K)
    {
      scale_psi_gamma <- (1 + (1/lambdasq_gamma_update[m, k]))
      psi_gamma_update_s[m, k] <- rinvgamma(1, shape = 1, scale = scale_psi_gamma) # IG distribution
    }
  }
  return(psi_gamma_update_s)
}

######
update_psi_gamma_BVS_HHS_MVN_cov2 <- function(lambdasq_gamma_update)
{
  scale_psi_gamma <- 1 + 1 / lambdasq_gamma_update

  matrix(
    rinvgamma(length(scale_psi_gamma),
              shape = 1,
              scale = as.vector(t(scale_psi_gamma))),
    nrow = nrow(lambdasq_gamma_update),
    byrow = TRUE
  )
}




set.seed(123)
x1 <- update_psi_gamma_BVS_HHS_MVN_cov(
  M, K,
  lambdasq_gamma_update
)

set.seed(123)
x2 <- update_psi_gamma_BVS_HHS_MVN_cov2(
  lambdasq_gamma_update
)

identical(x1, x2)
all.equal(x1, x2)




#######################################################################################

# Update psi_delta (for all j , m, and k = 1, 2, ..., K)
update_psi_delta_BVS_HHS_MVN_cov <- function(J, M, K, lambdasq_delta_update)
{
  psi_delta_update_s <- matrix(NA, nrow = J*M, ncol = K)

  for(jm in 1:(J*M))
  {
    for(k in 1:K)
    {
      scale_psi_delta <- (1 + (1/lambdasq_delta_update[jm, k]))
      psi_delta_update_s[jm, k] <- rinvgamma(1, shape = 1, scale = scale_psi_delta) # IG distribution
    }
  }
  return(psi_delta_update_s)
}

###
update_psi_delta_BVS_HHS_MVN_cov2 <- function(lambdasq_delta_update)
{
  scale_psi_delta <- 1 + 1 / lambdasq_delta_update

  matrix(
    rinvgamma(length(scale_psi_delta),
              shape = 1,
              scale = as.vector(t(scale_psi_delta))),
    nrow = nrow(lambdasq_delta_update),
    byrow = TRUE
  )
}



set.seed(123)
x1 <- update_psi_delta_BVS_HHS_MVN_cov(
  J, M, K,
  lambdasq_delta_update)

set.seed(123)
x2 <- update_psi_delta_BVS_HHS_MVN_cov2(
  lambdasq_delta_update)

identical(x1, x2)
all.equal(x1, x2)

################################################################
# Time check
#################################################################
load("exposome.RData")
load("modifiers_covariates.rda")
load("exposures.rda")
load("cov_base_post.rda")


raw_Y <- cbind(
  phenotype[,5:6],
  cov_base_post)

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



set.seed(6352)
sim_data <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D,
                              Sigma=Sigma)

Y <- sim_data$Y
x <- sim_data$x
z <- sim_data$z
u <- sim_data$u
d <- sim_data$d

W <- cbind(1, x, z, u, d)
WTW = t(W) %*% W
n_all_par <- ncol(W)

library(LaplacesDemon)

nu_0 <- K + 2
Psi_0 <- diag(K)
Sigma_init <- rinvwishart(nu_0, Psi_0)

library(Matrix)

benchmark_result <- microbenchmark(
  res_modify2 <- fit_BVS_HHS_MVN_cov_modify2(
    niter = 5, burn_in = 2, thin = 1,
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
    sigmasq_varphi = 10),

  res_modify3 = fit_HHS_MVN_cov_modify3(
    niter = 5, burn_in = 2, thin = 1,
    n = n, K = K, Y = Y, W = W, WTW = WTW,
    sigmasq_alpha = 100,
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
    sigmasq_varphi = 10),

  times = 5,
  unit = "microseconds"
)

benchmark_result



#################################################################
# NHHS model
#################################################################

###############################################################
# HHSSI model
#################################################################

J <- 5
M <- 6
K <- 4

lambdasq_beta_update <- rgamma(J, 2, 3)
lambdasq_gamma_update <- rgamma(M, 2, 3)
lambdasq_delta_update <- rgamma(J*M, 2, 3)
tausq_delta_update <- 1.2

sigmasq_delta_update <- rep(NA, J * M)
for (j in 1:J)
{
  for (m in 1:M)
  {
    jm <- (j - 1) * M + m   # row index for (j,m)
    sigmasq_delta_update[jm] <- (lambdasq_delta_update[jm] *
                                   lambdasq_beta_update[j] *
                                   lambdasq_gamma_update[m] *
                                   tausq_delta_update)
  }
}
######
beta_rep  <- rep(lambdasq_beta_update, each = M)
gamma_rep <- rep(lambdasq_gamma_update, times = J)

sigmasq_delta_update2 <-
  lambdasq_delta_update * beta_rep * gamma_rep * tausq_delta_update


identical(sigmasq_delta_update, sigmasq_delta_update2)
all.equal(sigmasq_delta_update, sigmasq_delta_update2)

cbind(sigmasq_delta_update, sigmasq_delta_update2)

#################################################################

update_lambdasq_delta_HHSSI_MVN_cov <- function(J, M, K, delta_update, lambdasq_beta_update,
                                                lambdasq_gamma_update, psi_delta_update, tausq_delta_update)
{
  lambdasq_delta_update_s <- rep(NA, J*M)
  shape_lambdasq_delta <- (K+1)/2

  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction
      scale_lambdasq_delta <- ((1/psi_delta_update[jm]) + (sum((delta_update[jm, ])^2) /
                                                             (2 * lambdasq_beta_update[j] * lambdasq_gamma_update[m] * tausq_delta_update)))

      lambdasq_delta_update_s[jm] <- rinvgamma(1, shape = shape_lambdasq_delta, scale = scale_lambdasq_delta) # IG distribution
    }
  }
  return(lambdasq_delta_update_s)
}

#####
update_lambdasq_delta_HHSSI_MVN_cov2 <- function(
    J, M, K, delta_update,
    lambdasq_beta_update,
    lambdasq_gamma_update,
    psi_delta_update,
    tausq_delta_update)
{
  shape_lambdasq_delta <- (K + 1) / 2

  # Sum of squares for each (j,m)
  delta_ss <- rowSums(delta_update^2)

  # Repeat beta and gamma to match the ordering of delta_update
  beta_rep  <- rep(lambdasq_beta_update, each = M)
  gamma_rep <- rep(lambdasq_gamma_update, times = J)

  scale_lambdasq_delta <-
    1 / psi_delta_update +
    delta_ss / (2 * beta_rep * gamma_rep * tausq_delta_update)

  rinvgamma(
    J * M,
    shape = shape_lambdasq_delta,
    scale = scale_lambdasq_delta
  )
}


library(LaplacesDemon)
psi_delta_update <- rgamma(J*M, 2, 3)
delta_update <- matrix(rnorm(J*M*K), nrow=J*M, ncol = K)

set.seed(123)
x1 <- update_lambdasq_delta_HHSSI_MVN_cov(J, M, K, delta_update,
                                          lambdasq_beta_update,
                                          lambdasq_gamma_update,
                                          psi_delta_update,
                                          tausq_delta_update)

set.seed(123)
x2 <- update_lambdasq_delta_HHSSI_MVN_cov2(J, M, K, delta_update,
                                           lambdasq_beta_update,
                                           lambdasq_gamma_update,
                                           psi_delta_update,
                                           tausq_delta_update)

identical(x1, x2)
all.equal(x1, x2)

##################################################################

update_tausq_delta_HHSSI_MVN_cov <- function(J, M, K, delta_update, lambdasq_delta_update,
                                             lambdasq_beta_update, lambdasq_gamma_update, xi_delta_update)
{
  sum_term_tausq_delta <- 0
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction

      for (k in 1:K)
      {
        sum_term_tausq_delta <- sum_term_tausq_delta + ((delta_update[jm, k])^2 / (2 *
                                                                                     lambdasq_delta_update[jm] * lambdasq_beta_update[j] * lambdasq_gamma_update[m]))
      }
    }
  }
  shape_tausq_delta <- (J*M*K + 1) / 2
  scale_tausq_delta <- ((1/xi_delta_update) + sum_term_tausq_delta)
  tausq_delta_update_s <- rinvgamma(1, shape = shape_tausq_delta, scale = scale_tausq_delta) # IG distribution

  return(tausq_delta_update_s)
}

#########
update_tausq_delta_HHSSI_MVN_cov2 <- function(
    J, M, K,
    delta_update,
    lambdasq_delta_update,
    lambdasq_beta_update,
    lambdasq_gamma_update,
    xi_delta_update)
{
  # Repeat beta and gamma to match rows of delta_update
  beta_rep  <- rep(lambdasq_beta_update, each = M)
  gamma_rep <- rep(lambdasq_gamma_update, times = J)

  # Sum of squares for each (j,m)
  delta_ss <- rowSums(delta_update * delta_update)

  sum_term_tausq_delta <- sum(
    delta_ss /
      (2 * lambdasq_delta_update * beta_rep * gamma_rep)
  )

  shape_tausq_delta <- (J * M * K + 1) / 2
  scale_tausq_delta <- 1 / xi_delta_update + sum_term_tausq_delta

  rinvgamma(
    1,
    shape = shape_tausq_delta,
    scale = scale_tausq_delta
  )
}


xi_delta_update <- 1.2

set.seed(123)
x1 <- update_tausq_delta_HHSSI_MVN_cov(J, M, K,
                                          delta_update,
                                          lambdasq_delta_update,
                                          lambdasq_beta_update,
                                          lambdasq_gamma_update,
                                          xi_delta_update)

set.seed(123)
x2 <- update_tausq_delta_HHSSI_MVN_cov2(J, M, K,
                                           delta_update,
                                           lambdasq_delta_update,
                                           lambdasq_beta_update,
                                           lambdasq_gamma_update,
                                           xi_delta_update)

identical(x1, x2)
all.equal(x1, x2)


beta_rep  <- rep(lambdasq_beta_update, each = M)
gamma_rep <- rep(lambdasq_gamma_update, times = J)
delta_ss <- rowSums(delta_update * delta_update)
sum_loop <- 0
for (j in 1:J)
  for (m in 1:M) {
    jm <- (j - 1) * M + m
    for (k in 1:K)
      sum_loop <- sum_loop +
        delta_update[jm, k]^2 /
        (2 * lambdasq_delta_update[jm] *
           lambdasq_beta_update[j] *
           lambdasq_gamma_update[m])
  }

sum_vec <- sum(delta_ss /
                 (2 * lambdasq_delta_update *
                    beta_rep *
                    gamma_rep))

sprintf("%.17f", sum_loop)
sprintf("%.17f", sum_vec)
sum_loop - sum_vec

###############################################################

update_lambdasq_delta_NHHSSI_MVN_cov <- function(J, M, K, delta_update,
                                                 psi_delta_update,
                                                 tausq_delta_update)
{
  lambdasq_delta_update_s <- rep(NA, J*M)
  shape_lambdasq_delta <- (K+1)/2

  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction
      scale_lambdasq_delta <- ((1/psi_delta_update[jm]) + (sum((delta_update[jm, ])^2) /
                                                             (2 * tausq_delta_update)))

      lambdasq_delta_update_s[jm] <- rinvgamma(1, shape = shape_lambdasq_delta, scale = scale_lambdasq_delta) # IG distribution
    }
  }
  return(lambdasq_delta_update_s)
}

###
update_lambdasq_delta_NHHSSI_MVN_cov2 <- function(
    J, M, K, delta_update,
    lambdasq_beta_update,
    lambdasq_gamma_update,
    psi_delta_update,
    tausq_delta_update)
{
  shape_lambdasq_delta <- (K + 1) / 2

  # Sum of squares for each (j,m)
  delta_ss <- rowSums(delta_update^2)
  scale_lambdasq_delta <- 1 / psi_delta_update + delta_ss / (2 * tausq_delta_update)

  lambdasq_delta_update_s <- rinvgamma(
    J * M,
    shape = shape_lambdasq_delta,
    scale = scale_lambdasq_delta
  )
}




set.seed(123)
x1 <- update_lambdasq_delta_NHHSSI_MVN_cov(J, M, K, delta_update,
                                           psi_delta_update,
                                           tausq_delta_update)

set.seed(123)
x2 <- update_lambdasq_delta_NHHSSI_MVN_cov2(J, M, K, delta_update,
                                            lambdasq_beta_update,
                                            lambdasq_gamma_update,
                                            psi_delta_update,
                                            tausq_delta_update)

identical(x1, x2)
all.equal(x1, x2)

##########################################

update_tausq_delta_NHHSSI_MVN_cov <- function(J, M, K, delta_update,
                                              lambdasq_delta_update,
                                              xi_delta_update)
{
  sum_term_tausq_delta <- 0
  for (j in 1:J)
  {
    for (m in 1:M)
    {
      jm <- (j - 1) * M + m  # correct order based on E matrix construction

      for (k in 1:K)
      {
        sum_term_tausq_delta <- sum_term_tausq_delta + ((delta_update[jm, k])^2 / (2 *
                                                                                     lambdasq_delta_update[jm]))
      }
    }
  }
  shape_tausq_delta <- (J*M*K + 1) / 2
  scale_tausq_delta <- ((1/xi_delta_update) + sum_term_tausq_delta)
  tausq_delta_update_s <- rinvgamma(1, shape = shape_tausq_delta, scale = scale_tausq_delta) # IG distribution

  return(tausq_delta_update_s)
}
#######
update_tausq_delta_NHHSSI_MVN_cov2 <- function(J, M, K, delta_update,
                                               lambdasq_delta_update,
                                               xi_delta_update)
{
  # Sum of squares for each (j,m)
  delta_ss <- rowSums(delta_update * delta_update)

  sum_term_tausq_delta <- sum(delta_ss / (2 * lambdasq_delta_update))

  shape_tausq_delta <- (J * M * K + 1) / 2
  scale_tausq_delta <- 1 / xi_delta_update + sum_term_tausq_delta

  tausq_delta_update_s <- rinvgamma(
    1,
    shape = shape_tausq_delta,
    scale = scale_tausq_delta
  )
  return(tausq_delta_update_s)
}



set.seed(123)
x1 <- update_tausq_delta_NHHSSI_MVN_cov(J, M, K, delta_update,
                                        lambdasq_delta_update,
                                        xi_delta_update)

set.seed(123)
x2 <- update_tausq_delta_NHHSSI_MVN_cov2(J, M, K, delta_update,
                                         lambdasq_delta_update,
                                         xi_delta_update)

identical(x1, x2)
all.equal(x1, x2)

sprintf("%.17f", x1)
sprintf("%.17f", x2)





