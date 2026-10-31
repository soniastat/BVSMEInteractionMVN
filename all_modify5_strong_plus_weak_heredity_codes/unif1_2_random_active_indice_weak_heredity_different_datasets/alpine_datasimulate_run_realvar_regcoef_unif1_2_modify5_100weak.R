
# Testing codes

sample_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat("Running sample:", sample_id, "\n")

###################################################

source("/projects/soniast@colostate.edu/alpine_HHSSI_function_regcoef_unif1_2_modify5_weak.R")



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
sim_data <- Sim_data_BVS_real(
  K = K,
  n = n,
  J = J,
  M = M,
  O = O,
  x=X, z=Z, d=D,
  Sigma = Sigma,
  weak_per = 1
)




dir.create("/projects/soniast@colostate.edu/results_datasimulate_realvar_regcoef_unif1_2_modify5_100weak", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_datasimulate_realvar_regcoef_unif1_2_modify5_100weak/results_datasimulate_realvar_regcoef_unif1_2_modify5_100weak_sample_", sample_id, ".rds")

saveRDS(sim_data, file=outfile)

cat("Finished datasimulate sample", sample_id, "\n")



#############################################################


