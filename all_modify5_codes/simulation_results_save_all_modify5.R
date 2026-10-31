

##################################################
# HHS model

# save 100 datest as a list
files <- file.path("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")


# save 100 datest as a list
files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")


# save 100 datest as a list
files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")





# save 100 datest as a list
files <- file.path("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")


# save 100 datest as a list
files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_HHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_HHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHS_simvar_regcoef_normal0_1_modify5",
                   paste0("results_HHS_simvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_simvar_regcoef_normal0_1_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHS_realvar_regcoef_normal0_1_modify5",
                   paste0("results_HHS_realvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_realvar_regcoef_normal0_1_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHS_simvar_regcoef_unif1_2_modify5",
                   paste0("results_HHS_simvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_simvar_regcoef_unif1_2_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHS_realvar_regcoef_unif1_2_modify5",
                   paste0("results_HHS_realvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_realvar_regcoef_unif1_2_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHS_simvar_regcoef_unif2_3_modify5",
                   paste0("results_HHS_simvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_simvar_regcoef_unif2_3_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHS_realvar_regcoef_unif2_3_modify5",
                   paste0("results_HHS_realvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_HHS_list <- lapply(files, readRDS)
save(theta_update_HHS_list,
     file = "practice_data/results_HHS_realvar_regcoef_unif2_3_modify5.rda")

##################################################


##################################################
# HHSSI model

# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")





# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_HHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_simvar_regcoef_normal0_1_modify5",
                   paste0("results_HHSSI_simvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_simvar_regcoef_normal0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_simvar_regcoef_unif1_2_modify5",
                   paste0("results_HHSSI_simvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_simvar_regcoef_unif1_2_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5",
                   paste0("results_HHSSI_realvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_realvar_regcoef_unif1_2_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_simvar_regcoef_unif2_3_modify5",
                   paste0("results_HHSSI_simvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_simvar_regcoef_unif2_3_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_HHSSI_realvar_regcoef_unif2_3_modify5",
                   paste0("results_HHSSI_realvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
save(theta_update_HHSSI_list,
     file = "practice_data/results_HHSSI_realvar_regcoef_unif2_3_modify5.rda")

##################################################






##################################################
# NHHS model

# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_NHHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5",
                   paste0("results_NHHS_realvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_realvar_regcoef_normal0_1_modify5.rda")





# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_simvar_regcoef_normal0_1_modify5",
                   paste0("results_NHHS_simvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_simvar_regcoef_normal0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_simvar_regcoef_unif1_2_modify5",
                   paste0("results_NHHS_simvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_simvar_regcoef_unif1_2_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5",
                   paste0("results_NHHS_realvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_realvar_regcoef_unif1_2_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_simvar_regcoef_unif2_3_modify5",
                   paste0("results_NHHS_simvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_simvar_regcoef_unif2_3_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHS_realvar_regcoef_unif2_3_modify5",
                   paste0("results_NHHS_realvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_NHHS_list <- lapply(files, readRDS)
save(theta_update_NHHS_list,
     file = "practice_data/results_NHHS_realvar_regcoef_unif2_3_modify5.rda")

##################################################





##################################################
# NHHSSI model

# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_delta_notdepends_betagamma_modify5",
                   paste0("results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif0_1_delta_notdepends_betagamma_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_modify5",
                   paste0("results_NHHSSI_simvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_simvar_regcoef_normal0_1_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_unif1_2_modify5",
                   paste0("results_NHHSSI_simvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_simvar_regcoef_unif1_2_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5",
                   paste0("results_NHHSSI_realvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_realvar_regcoef_unif1_2_modify5.rda")




# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_unif2_3_modify5",
                   paste0("results_NHHSSI_simvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_simvar_regcoef_unif2_3_modify5.rda")



# save 100 datest as a list
files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_unif2_3_modify5",
                   paste0("results_NHHSSI_realvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
save(theta_update_NHHSSI_list,
     file = "practice_data/results_NHHSSI_realvar_regcoef_unif2_3_modify5.rda")

##################################################









######################################################
# Simulated data
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
               paste0("results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5.rda")



# Simulated data
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_depends_betagamma_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_delta_notdepends_betagamma_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_datasimulate_simvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif0_1_modify5",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_unif0_1_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif0_1_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_modify5.rda")





# Simulated data
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_normal0_1_modify5",
                   paste0("results_datasimulate_simvar_regcoef_normal0_1_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_normal0_1_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5",
                   paste0("results_datasimulate_realvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_unif1_2_modify5.rda")



# Simulated data
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_unif1_2_modify5",
                   paste0("results_datasimulate_simvar_regcoef_unif1_2_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_unif1_2_modify5.rda")





# Simulated data
files <- file.path("practice_data/results_datasimulate_realvar_regcoef_unif2_3_modify5",
                   paste0("results_datasimulate_realvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_unif2_3_modify5.rda")




# Simulated data
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_unif2_3_modify5",
                   paste0("results_datasimulate_simvar_regcoef_unif2_3_modify5_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_unif2_3_modify5.rda")



######################################################

# Weak heredity
files <- file.path("practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak",
                   paste0("results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")



files <- file.path("practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak",
                   paste0("results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak_sample_", 1:100, ".rds"))
sim_data_list <- lapply(files, readRDS)
length(sim_data_list)
# 100
save(sim_data_list, file = "practice_data/results_datasimulate_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")




files <- file.path("practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak",
                   paste0("results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")





files <- file.path("practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak",
                   paste0("results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak_sample_", 1:100, ".rds"))
theta_update_HHSSI_list <- lapply(files, readRDS)
length(theta_update_HHSSI_list)
# 100
save(theta_update_HHSSI_list, file = "practice_data/results_HHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")





files <- file.path("practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak",
                   paste0("results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_simvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")





files <- file.path("practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak",
                   paste0("results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak_sample_", 1:100, ".rds"))
theta_update_NHHSSI_list <- lapply(files, readRDS)
length(theta_update_NHHSSI_list)
# 100
save(theta_update_NHHSSI_list, file = "practice_data/results_NHHSSI_realvar_regcoef_normal0_1_alpha_unif_3.5_2.5_modify5_weak.rda")



