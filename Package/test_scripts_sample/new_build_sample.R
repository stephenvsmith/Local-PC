######################################################################################################

# Author: Stephen Smith
# Date: 10/22/2019
# Second attempt at the build_results file for sample version

######################################################################################################

##### Variable Creation #####

# Storing Results
result_cols <- c("Network","Target","Number_of_Targets","Time_to_Build",
                                "Undirected Missing","Directed Missing","Total Missing",
                                "Undirected Added","Directed Added","Total Added",
                                "Wrong_Direction","Directed_Undirected","Undirected_Directed",
                                "SHD","tp","fp","fn","Same_V_Structures")
results <- matrix(nrow=0,ncol = length(result_cols))
colnames(results) <- result_cols

###### Workhorse Function #####

br <- function(net,vars){
  cat("Working on network",net,"\n")
  max_target <- 5
  # Store results for current network
  vars$current_network <- matrix(nrow=0,ncol = length(result_cols))
  colnames(vars$current_network) <- result_cols
  
  vars <- NetSetup(net,vars)
  if (nrow(vars$true_dag)>45 | nrow(vars$true_dag)<15) {
    max_target <- 1
  } 
  
  vars <- SaveDAGPhotos(net,vars)
  
  for (num_sep_nodes in 1:max_target){
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
var_list$results <- results

### Define paths for storing
package_loc <- '~/Desktop/Research/Package/'
var_list$ps2 <- '/home/stephen/'
var_list$result_dir <- '~/Desktop/Research/Results/sample/'
var_list$simulation_dir <- '/home/stephen/Desktop/Research/Simulations/'
var_list$rds_dir <- paste0(var_list$ps2,"Desktop/Research/Networks/rds/")

### Source Functions  
source('~/Desktop/Research/Package/test_scripts_sample/helper_files_sample.R')

### Load Libraries

LoadLibraries(package_loc)

### Generate Data Grid
var_list$data_grid <- GenerateDataGrid() 

### Names of the Networks
var_list$net_names <- GetNetworkNames(var_list$ps2)
completed <- list.files("~/Desktop/Research/Results/sample")
var_list$net_names <- setdiff(var_list$net_names,completed)


### Build results for each network

all_var_lists <- sapply(var_list$net_names,br,var_list)

















