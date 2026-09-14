

# Simulated data for all var real

sample_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat("Running sample:", sample_id, "\n")

###################################################

source("/projects/soniast@colostate.edu/alpine_sim_HHS_modify4.R")



library(LaplacesDemon)



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
# Simulate one dataset all var real
############################

set.seed(23444 + sample_id)

sim_data <- Sim_data_BVS_real_elevated(K=K, n=n, J=J, M=M,
                                       O=O, x=X, z=Z, d=D, Sigma=Sigma)


dir.create("/projects/soniast@colostate.edu/results_sample_data_allvarreal_modify4_elevated", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_sample_data_allvarreal_modify4_elevated/modify4_elevated_sample_data_allvarreal_", sample_id, ".rds")

saveRDS(sim_data, file=outfile)

cat("Finished sample", sample_id, "\n")







