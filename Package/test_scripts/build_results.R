### Script for building on many datasets with local pc

### Initial setup

important_vars <- list()

source("~/Desktop/Research/package/helper_scripts/initialSetup.R")

# Adding helper functions
source(paste0(path_start,"Desktop/Research/package/helper_scripts/sourcingScript.R"))

# Storing Results
results <- matrix(nrow=0,ncol = 5)
colnames(results) <- c("Network","Target","Missing","Added","Same_V_Structures")

timeDiffs <- sapply(nets,function(net){
  cat("Working on network",net,"\n")
  start_time <- Sys.time()

  # Store results for current network
  current_network <- matrix(nrow=0,ncol = 5)
  colnames(current_network) <- c("Network","Target","Missing","Added","Same_V_Structures")

  source(paste0(path_start,"Desktop/Research/package/helper_scripts/netSetup.R"),local = TRUE)

  source(paste0(path_start,"Desktop/Research/package/helper_scripts/saveDAGphotos.R"),local = TRUE)
  # Move to next network if there are more than 50 nodes (only want small networks for now)
  if (ncol(true_dag) > 50)
    return(NA)

  # Get all results for one target
  if (!file.exists("./one_target")) dir.create("./one_target")
  setwd("./one_target")
  lpc_builds(true_dag,true_cpdag,1,net)

  # Get all results for two targets
  setwd("..")
  if (!file.exists("./two_targets")) dir.create("./two_targets")
  setwd("./two_targets")
  lpc_builds(true_dag,true_cpdag,2,net)


  # Get all results for three targets
  setwd("..")
  if (!file.exists("./three_targets")) dir.create("./three_targets")
  setwd("./three_targets")
  lpc_builds(true_dag,true_cpdag,3,net)

  # Get all results for four targets
  setwd("..")
  if (!file.exists("./four_targets")) dir.create("./four_targets")
  setwd("./four_targets")
  lpc_builds(true_dag,true_cpdag,4,net)

  # Get all results for five targets
  setwd("..")
  if (!file.exists("./five_targets")) dir.create("./five_targets")
  setwd("./five_targets")
  lpc_builds(true_dag,true_cpdag,5,net)


  setwd("..")

  write.table(current_network,paste0(net,"_results.txt"))

  setwd("..")

  end_time <- Sys.time()
  return(end_time-start_time)

})

setwd("..")
write.table(results,file = "finalResults.txt")
