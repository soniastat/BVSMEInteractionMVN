
# Testing codes

sample_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat("Running sample:", sample_id, "\n")

###################################################

source("/projects/soniast@colostate.edu/alpine_HHSSI_function_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.R")



library(LaplacesDemon)
library(Matrix)




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
sim_data <- Sim_data_BVS_weak_real_sim_var(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  Sigma = Sigma
)




dir.create("/projects/soniast@colostate.edu/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak_sample_", sample_id, ".rds")

saveRDS(sim_data, file=outfile)

cat("Finished datasimulate sample", sample_id, "\n")



#############################################################


