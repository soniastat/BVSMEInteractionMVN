##################################################################################
# Different number of covariates
##################################################################################

# J = 10, M = 5
sim_data_list1 <- list()
for(sample in 1:100)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 10, M = 5)
  sim_data_list1[[sample]] <- sim_data
}


sim_data$alpha_true
sim_data$beta_true
sim_data$gamma_true
sim_data$delta_true
sim_data$x
###############################################

# J = 15, M = 5
sim_data_list1 <- list()
for(sample in 1:100)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 15, M = 5)
  sim_data_list1[[sample]] <- sim_data
}

###############################################

# J = 20, M = 5
sim_data_list1 <- list()
for(sample in 1:100)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 20, M = 5)
  sim_data_list1[[sample]] <- sim_data
}

###############################################

# J = 25, M = 5
sim_data_list1 <- list()
for(sample in 1:100)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 25, M = 5)
  sim_data_list1[[sample]] <- sim_data
}

###############################################

# J = 50, M = 5
sim_data_list1 <- list()
for(sample in 1:100)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 50, M = 5)
  sim_data_list1[[sample]] <- sim_data
}


###############################################

# J = 75, M = 5
sim_data_list1 <- list()
for(sample in 1:100)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 75, M = 5)
  sim_data_list1[[sample]] <- sim_data
}

###############################################

# J = 100, M = 5
sim_data_list1 <- list()
for(sample in 1:100)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 100, M = 5)
  sim_data_list1[[sample]] <- sim_data
}

#####################################################################################

# J = 10, M = 20
sim_data_list1 <- list()
for(sample in 1:10)
{
  print(sample)
  set.seed(23444 + sample)
  sim_data <- Sim_data_BVS(K = 5, n = 300, J = 20, M = 20)
  sim_data_list1[[sample]] <- sim_data
}

#####################################################################################

