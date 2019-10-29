######################################################################################################

# Author: Stephen Smith
# Date: 10/22/2019
# Second attempt at the build_results file

######################################################################################################

##### Variable Creation #####

# Storing Results
result_cols <- c("Network","Target","Number_of_Targets","Time_to_Build","Missing","Added","Wrong_Direction","SHD","Same_V_Structures")
results <- matrix(nrow=0,ncol = length(result_cols))
colnames(results) <- result_cols

###### Workhorse Function #####

br <- function(net,vars){
  cat("Working on network",net,"\n")
  
  # Store results for current network
  vars$current_network <- matrix(nrow=0,ncol = length(result_cols))
  colnames(vars$current_network) <- result_cols
  
  vars <- NetSetup(net,vars)
  if (nrow(vars$true_dag)>50) return()
    
  vars <- SaveDAGPhotos(net,vars)
  
  for (num_sep_nodes in 1:5){
    # Building the result file
    if (!dir.exists(paste0("./",num_sep_nodes," target"))) dir.create(paste0("./",num_sep_nodes," target"))
    setwd(paste0("./",num_sep_nodes," target"))
    
    vars <- lpc_builds(vars,net,num_sep_nodes)
    
    setwd('..')
  }
  
  write.table(vars$current_network,paste0(net," results.txt"))
  
  return(vars)
}



##### Script #####

var_list <- list()
# Storing Results
var_list$results <- matrix(nrow=0,ncol = 5)
colnames(var_list$results) <- c("Network","Target","Missing","Added","Same_V_Structures")

### Define paths for storing
package_loc <- '~/Desktop/Research/package/'
var_list$ps2 <- '/home/stephen/'
var_list$result_dir <- '~/Desktop/Research/Results/'
var_list$simulation_dir <- '/home/stephen/Desktop/Research/Simulations/'
var_list$rds_dir <- paste0(var_list$ps2,"Desktop/Research/projects/bn_data_generation/networks/rds/")

### Source Functions  
source('~/Desktop/Research/package/test_scripts/helper_files.R')

### Load Libraries

LoadLibraries(package_loc)

### Generate Data Grid
var_list$data_grid <- GenerateDataGrid() 

### Names of the Networks
var_list$net_names <- GetNetworkNames(var_list$ps2)

### Build results for each network

all_var_lists <- sapply(var_list$net_names,br,var_list)

















