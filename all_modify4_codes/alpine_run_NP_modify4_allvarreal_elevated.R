
# Testing codes for NP-MVN-cov

sample_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat("Running sample:", sample_id, "\n")

###################################################

source("/projects/soniast@colostate.edu/alpine_sim_NP_modify4.R")



library(LaplacesDemon)
library(Matrix)



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


Sigma <- cov(Y0)
K <- dim(Y0)[2]
n <- nrow(Y0)
J <- ncol(X)
M <- ncol(Z)
O <- ncol(D)

############################
# Simulate one dataset
############################

set.seed(23444 + sample_id)

sim_data <- Sim_data_BVS_real_elevated(K=K, n=n, J=J, M=M,
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





############################
# Fit model
############################

set.seed(765878 + sample_id)

theta_update_save <- tryCatch({

  res <- fit_NP_MVN_cov_modify3(
    niter = 6000, burn_in = 1000, thin = 5,
    n=n, K=K, Y=Y, W=W, WTW=WTW, n_all_par=n_all_par,
    J=J, M=M, O=O,
    sigmasq_alpha = 100,
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




dir.create("/projects/soniast@colostate.edu/results_NP_modify4_allvarreal_elevated", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_NP_modify4_allvarreal_elevated/modify4_allvarreal_elevated_NP_cov_sample_", sample_id, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished sample", sample_id, "\n")






