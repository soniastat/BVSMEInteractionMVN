
# Testing codes for BVS-NP-MVN-cov

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
load("cov_base_post.rda")

# covariates_post <- codebook %>%
#   filter(domain=="Covariates" & period=="Postnatal") %>%
#   pull(variable_name) %>%
#   as.character()
#
# covariates_post <- covariates_post[
#   covariates_post %in%
#     c("hs_c_height_None","hs_c_weight_None")
# ]
#
# cov_base_post <- covariates[,covariates_post]
# save(cov_base_post,
#      file="cov_base_post.rda")



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

  res <- fit_NP_MVN_cov(
    niter = 2, burn_in = 0, thin = 1,
    n=n, K=K, Y=Y, W=W, n_all_par=n_all_par,
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
    sample=sample
  )

})

dir.create("results", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("results/NP_cov_sample_", sample, ".rds")

saveRDS(theta_update_save, file=outfile)

cat("Finished sample", sample, "\n")






files <- list.files(
  "results",
  pattern = "\\.rds$",
  full.names = TRUE
)

files <- files[!grepl("theta_update_list", files)]



theta_update_list <- lapply(files, readRDS)
length(theta_update_list)


saveRDS(
  theta_update_list,
  file = "results/theta_update_NP_cov_list.rds"
)


# in terminal, run these

# cd /Users/sonia/Desktop/Packages/BVSMEInteractionMVN

# seq 1 5 | xargs -n 1 -P 5 Rscript run_one_sample_NP_cov.R



