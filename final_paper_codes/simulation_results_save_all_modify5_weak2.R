
rm(list = ls())
gc()

######################################################

# strong + 0% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_0weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_modify5_0weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_0weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_0weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_0weak.rda")


files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_0weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_0weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_0weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_0weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_0weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_0weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_0weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_modify5_0weak.rda")

######################################################

# strong + 10% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_10weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_modify5_10weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_10weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_10weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_10weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_10weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_10weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_10weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_10weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_10weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_10weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_10weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_modify5_10weak.rda")


#############


# strong + 25% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_25weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_modify5_25weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_25weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_25weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_25weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_25weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_25weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_25weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_25weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_25weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_25weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_25weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_modify5_25weak.rda")


########

# strong + 50% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_50weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_modify5_50weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_50weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_50weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_50weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_50weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_50weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_50weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_50weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_50weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_50weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_50weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_modify5_50weak.rda")




########

# strong + 100% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_100weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_modify5_100weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5_100weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_100weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_modify5_100weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_100weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5_100weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_100weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5_100weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_100weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5_100weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_modify5_100weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_modify5_100weak.rda")








######################################################
# unif(1, 2)
######################################################

# strong + 0% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_0weak",
                   paste0("results_datasimulate_realvar_regcoef_unif1_2_modify5_0weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_0weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_0weak",
                   paste0("results_HHS_realvar_regcoef_unif1_2_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_0weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_0weak",
                   paste0("results_HHSSI_realvar_regcoef_unif1_2_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_0weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_0weak",
                   paste0("results_NHHS_realvar_regcoef_unif1_2_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_0weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_0weak",
                   paste0("results_NHHSSI_realvar_regcoef_unif1_2_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_0weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_0weak",
                   paste0("results_NP_realvar_regcoef_unif1_2_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_unif1_2_modify5_0weak.rda")

######################################################

# strong + 10% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_10weak",
                   paste0("results_datasimulate_realvar_regcoef_unif1_2_modify5_10weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_10weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_10weak",
                   paste0("results_HHS_realvar_regcoef_unif1_2_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_10weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_10weak",
                   paste0("results_HHSSI_realvar_regcoef_unif1_2_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_10weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_10weak",
                   paste0("results_NHHS_realvar_regcoef_unif1_2_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_10weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_10weak",
                   paste0("results_NHHSSI_realvar_regcoef_unif1_2_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_10weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_10weak",
                   paste0("results_NP_realvar_regcoef_unif1_2_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_unif1_2_modify5_10weak.rda")


#############


# strong + 25% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_25weak",
                   paste0("results_datasimulate_realvar_regcoef_unif1_2_modify5_25weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_25weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_25weak",
                   paste0("results_HHS_realvar_regcoef_unif1_2_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_25weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_25weak",
                   paste0("results_HHSSI_realvar_regcoef_unif1_2_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_25weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_25weak",
                   paste0("results_NHHS_realvar_regcoef_unif1_2_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_25weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_25weak",
                   paste0("results_NHHSSI_realvar_regcoef_unif1_2_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_25weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_25weak",
                   paste0("results_NP_realvar_regcoef_unif1_2_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_unif1_2_modify5_25weak.rda")


########

# strong + 50% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_50weak",
                   paste0("results_datasimulate_realvar_regcoef_unif1_2_modify5_50weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_50weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_50weak",
                   paste0("results_HHS_realvar_regcoef_unif1_2_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_50weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_50weak",
                   paste0("results_HHSSI_realvar_regcoef_unif1_2_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_50weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_50weak",
                   paste0("results_NHHS_realvar_regcoef_unif1_2_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_50weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_50weak",
                   paste0("results_NHHSSI_realvar_regcoef_unif1_2_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_50weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_50weak",
                   paste0("results_NP_realvar_regcoef_unif1_2_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_unif1_2_modify5_50weak.rda")




########

# strong + 100% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_100weak",
                   paste0("results_datasimulate_realvar_regcoef_unif1_2_modify5_100weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5_100weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_100weak",
                   paste0("results_HHS_realvar_regcoef_unif1_2_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_unif1_2_modify5_100weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_100weak",
                   paste0("results_HHSSI_realvar_regcoef_unif1_2_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5_100weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_100weak",
                   paste0("results_NHHS_realvar_regcoef_unif1_2_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5_100weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_100weak",
                   paste0("results_NHHSSI_realvar_regcoef_unif1_2_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5_100weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_unif1_2_modify5_100weak",
                   paste0("results_NP_realvar_regcoef_unif1_2_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_unif1_2_modify5_100weak.rda")











################################################################
# beta, gamma, delta, varphi ~ N(0,1), alpha ~ unif(-3.5, -2.5)
################################################################

# strong + 0% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_0weak.rda")

######################################################

# strong + 10% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_10weak.rda")


#############


# strong + 25% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_25weak.rda")


########

# strong + 50% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_50weak.rda")




########

# strong + 100% of only Weak heredity
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")


files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
length(theta_update_HHS_list)
# 100
save(theta_update_HHS_list, file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")




files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")



files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
length(theta_update_NHHS_list)
# 100
save(theta_update_NHHS_list, file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")




files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")




files <- file.path("practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak",
                   paste0("results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak_sample_", 1:100, ".rds"))
theta_update_NP_list <- lapply(files, readRDS)
length(theta_update_NP_list)
# 100
save(theta_update_NP_list, file = "practice_data/results_NP_realvar_regcoef_normal0_1_alpha_3.5_2.5_modify5_100weak.rda")




