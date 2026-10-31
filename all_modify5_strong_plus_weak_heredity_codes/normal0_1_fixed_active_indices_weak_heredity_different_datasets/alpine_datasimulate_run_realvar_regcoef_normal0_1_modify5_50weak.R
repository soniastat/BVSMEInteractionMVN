
# Testing codes

sample_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat("Running sample:", sample_id, "\n")

###################################################

source("/projects/soniast@colostate.edu/alpine_HHSSI_function_regcoef_normal0_1_modify5_weak.R")



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

# Fix active weak-heredity interaction rows across all samples
grid <- expand.grid(m = 1:M, j = 1:J)
weak_only_idx <- with(grid, (j - 1) * M + m)[xor(grid$m %in% 1:3, grid$j %in% 1:2)]
weak_per <- 0.5
set.seed(87665467)
fixed_weak_rows <- sample(weak_only_idx, round(weak_per * length(weak_only_idx)))


set.seed(23444 + sample_id)
sim_data <- Sim_data_BVS_real(K=K, n=n, J=J, M=M,
                              O=O, x=X, z=Z, d=D,
                              Sigma=Sigma,
                              weak_per=0.5,
                              active_weak_idx=fixed_weak_rows)




dir.create("/projects/soniast@colostate.edu/results_datasimulate_realvar_regcoef_normal0_1_modify5_50weak", recursive=TRUE, showWarnings=FALSE)

outfile <- paste0("/projects/soniast@colostate.edu/results_datasimulate_realvar_regcoef_normal0_1_modify5_50weak/results_datasimulate_realvar_regcoef_normal0_1_modify5_50weak_sample_", sample_id, ".rds")

saveRDS(sim_data, file=outfile)

cat("Finished datasimulate sample", sample_id, "\n")



#############################################################


