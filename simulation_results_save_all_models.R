
# Simulated data
files <- file.path("results_sample_data",
               paste0("sample_data_", 1:100, ".rds"))

sim_data_list <- lapply(files, readRDS)

length(sim_data_list)
# 100

save(sim_data_list,
  file = "results/sample_data_list_Sigma.rda"
)




##################################################
# # HHS model
#
# # save 100 datest as a list
# files <- file.path("results_HHS",
#                    paste0("modify2_HHS_cov_sample_", 1:100, ".rds"))
#
# theta_update_HHS_list <- lapply(files, readRDS)
#
# save(theta_update_HHS_SI_list,
#      file = "results/theta_update_HHS_list.rda")


##################################################
# HHS model

# save 100 datest as a list
files <- file.path("results_HHS_modify2",
                   paste0("modify2_HHS_cov_sample_", 1:100, ".rds"))

theta_update_HHS_list <- lapply(files, readRDS)

save(theta_update_HHS_list,
     file = "results/theta_update_HHS_list_Sigma.rda")

##################################################
# HHS_SI model

# # save 100 datest as a list
# files <- file.path("results_HHS_SI",
#            paste0("modify2_HHS_SI_cov_sample_", 1:100, ".rds"))
#
# theta_update_HHS_SI_list <- lapply(files, readRDS)
#
# save(theta_update_HHS_SI_list,
#      file = "results/theta_update_HHS_SI_list.rda")

##################################################
# HHS_SI model

# save 100 datest as a list
files <- file.path("results_HHS_SI_modify2",
                   paste0("modify2_HHS_SI_cov_sample_", 1:100, ".rds"))

theta_update_HHS_SI_list <- lapply(files, readRDS)

save(theta_update_HHS_SI_list,
     file = "results/theta_update_HHS_SI_list_Sigma.rda")

###################################
# # HRHS model
#
# # save 100 datest as a list
# files <- file.path(
#   "results_HRHS",
#   paste0("modify2_HRHS_cov_sample_", 1:100, ".rds")
# )
#
# theta_update_HRHS_list <- lapply(files, readRDS)
#
# save(theta_update_HRHS_list,
#      file = "results/theta_update_HRHS_list.rda")

###################################
# HRHS model

# save 100 datest as a list
files <- file.path(
  "results_HRHS_modify2",
  paste0("modify2_HRHS_cov_sample_", 1:100, ".rds")
)

theta_update_HRHS_list <- lapply(files, readRDS)

save(theta_update_HRHS_list,
     file = "results/theta_update_HRHS_list_Sigma.rda")

###################
# # HRHS_SI model
#
# # save 100 datest as a list
# files <- file.path(
#   "results_HRHS_SI",
#   paste0("modify2_HRHS_SI_cov_sample_", 1:100, ".rds")
# )
#
# theta_update_HRHS_SI_list <- lapply(files, readRDS)
#
# save(theta_update_HRHS_SI_list,
#      file = "results/theta_update_HRHS_SI_list.rda")


###################
# HRHS_SI model

# save 100 datest as a list
files <- file.path(
  "results_HRHS_SI_modify2",
  paste0("modify2_HRHS_SI_cov_sample_", 1:100, ".rds")
)

theta_update_HRHS_SI_list <- lapply(files, readRDS)

save(theta_update_HRHS_SI_list,
     file = "results/theta_update_HRHS_SI_list_Sigma.rda")


######################################
# # NHHS model
#
# # save 100 datest as a list
# files <- file.path(
#   "results_NHHS",
#   paste0("modify2_NHHS_cov_sample_", 1:100, ".rds")
# )
#
# theta_update_NHHS_list <- lapply(files, readRDS)
#
# save(theta_update_NHHS_list,
#      file = "results/theta_update_NHHS_list.rda")

######################################
# NHHS model

# save 100 datest as a list
files <- file.path(
  "results_NHHS_modify2",
  paste0("modify2_NHHS_cov_sample_", 1:100, ".rds")
)

theta_update_NHHS_list <- lapply(files, readRDS)

save(theta_update_NHHS_list,
     file = "results/theta_update_NHHS_list_Sigma.rda")



###################
# # NHHS_SI model
#
# # save 100 datest as a list
# files <- file.path(
#   "results_NHHS_SI",
#   paste0("modify2_NHHS_SI_cov_sample_", 1:100, ".rds")
# )
#
# theta_update_NHHS_SI_list <- lapply(files, readRDS)
#
# save(theta_update_NHHS_SI_list,
#      file = "results/theta_update_NHHS_SI_list.rda")


###################
# NHHS_SI model

# save 100 datest as a list
files <- file.path(
  "results_NHHS_SI_modify2",
  paste0("modify2_NHHS_SI_cov_sample_", 1:100, ".rds")
)

theta_update_NHHS_SI_list <- lapply(files, readRDS)

save(theta_update_NHHS_SI_list,
     file = "results/theta_update_NHHS_SI_list_Sigma.rda")

#############################
# # NHRHS model
#
# # save 100 datest as a list
# files <- file.path(
#   "results_NHRHS",
#   paste0("modify2_NHRHS_cov_sample_", 1:100, ".rds")
# )
#
# theta_update_NHRHS_list <- lapply(files, readRDS)
#
# save(theta_update_NHRHS_list,
#      file = "results/theta_update_NHRHS_list.rda")

#############################
# NHRHS model

# save 100 datest as a list
files <- file.path(
  "results_NHRHS_modify2",
  paste0("modify2_NHRHS_cov_sample_", 1:100, ".rds")
)

theta_update_NHRHS_list <- lapply(files, readRDS)

save(theta_update_NHRHS_list,
     file = "results/theta_update_NHRHS_list_Sigma.rda")


####################
# # NHRHS_SI model
#
# # save 100 datest as a list
# files <- file.path(
#   "results_NHRHS_SI",
#   paste0("modify2_NHRHS_SI_cov_sample_", 1:100, ".rds")
# )
#
# theta_update_NHRHS_SI_list <- lapply(files, readRDS)
#
# save(theta_update_NHRHS_SI_list,
#      file = "results/theta_update_NHRHS_SI_list.rda")

####################
# NHRHS_SI model

# save 100 datest as a list
files <- file.path(
  "results_NHRHS_SI_modify2",
  paste0("modify2_NHRHS_SI_cov_sample_", 1:100, ".rds")
)

theta_update_NHRHS_SI_list <- lapply(files, readRDS)

save(theta_update_NHRHS_SI_list,
     file = "results/theta_update_NHRHS_SI_list_Sigma.rda")

####################################
# # NP model
#
# # save 100 datest as a list
# files <- file.path(
#   "results_NP",
#   paste0("modify_NP_cov_sample_", 1:100, ".rds")
# )
#
# theta_update_NP_list <- lapply(files, readRDS)
#
# save(theta_update_NP_list,
#      file = "results/theta_update_NP_list.rda")



####################################
# NP model

# save 100 datest as a list
files <- file.path(
  "results_NP_modify",
  paste0("modify_NP_cov_sample_", 1:100, ".rds")
)

theta_update_NP_list <- lapply(files, readRDS)

save(theta_update_NP_list,
     file = "results/theta_update_NP_list_Sigma.rda")




