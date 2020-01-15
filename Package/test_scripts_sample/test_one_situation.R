#### File for testing individual situations

network <- "link"
target <- 20

#### Get simulation information
setwd(paste0("~/Desktop/Research/Simulations/",network))
dir <- list.files()
for (d in dir){
  if (file.exists(paste0("./",d,"/trueDAG.txt"))){
    true_dag <- read.table(paste0("./",d,"/trueDAG.txt"))
    data <- read.table(paste0("./",d,"/data1.txt"))
  }
}
rds_dir <- paste0("~/Desktop/Research/Networks/rds/")
file <- paste0(rds_dir,network,".rds")

# Save names of nodes (variables)
node_names <- names(readRDS(file))
rownames(true_dag) <- node_names
colnames(true_dag) <- node_names
colnames(data) <- node_names


res <- local_pc2(true_dag = true_dag,target = target,lmax=3,verbose = FALSE,verbose_small = FALSE,pop = TRUE)


######################################################################################################################

vars <- list()
# Storing Results
vars$results <- matrix(nrow=0,ncol = 5)
colnames(vars$results) <- c("Network","Target","Missing","Added","Same_V_Structures")

### Define paths for storing
package_loc <- '~/Desktop/Research/Package/'
vars$ps2 <- '/home/stephen/'
vars$result_dir <- '~/Desktop/Research/Results/population'
vars$simulation_dir <- '/home/stephen/Desktop/Research/Simulations/'
vars$rds_dir <- paste0(vars$ps2,"Desktop/Research/Networks/rds/")

### Source Functions  
source('~/Desktop/Research/Package/test_scripts/helper_files.R')

### Load Libraries

LoadLibraries(package_loc)

### Generate Data Grid
vars$data_grid <- GenerateDataGrid() 

### Names of the Networks
vars$net_names <- GetNetworkNames(vars$ps2)

cat("Working on network",network,"\n")
# Store results for current network
vars$current_network <- matrix(nrow=0,ncol = length(result_cols))
colnames(vars$current_network) <- result_cols

vars <- NetSetup(network,vars)

vars <- SaveDAGPhotos(network,vars)


lpc_builds <- function(vars,net,num_sep_nodes=1,target) {
  
  setwd("~/Desktop")
  
  # Create directory to store results for this test
  if (!file.exists(paste0(net,"; target=",paste(target,collapse = ","))))
    dir.create(paste0(net,"; target=",paste(target,collapse = ",")))
  setwd(paste0("./",net,"; target=",paste(target,collapse = ",")))
  # for ( i in 1:nrow(vars$true_dag)){
  #   cat(i,":",vars$node_names[i],"\n")
  # }
  vars <- test_bn(net,target,vars)
  metrics <- local_pc_dist(target,vars)
  
  vars$results <- rbind(results,c(net,paste(target,collapse = ","),length(target),vars$time,metrics))
  vars$current_network <- rbind(vars$current_network,c(net,paste(target,collapse = ","),length(target),vars$time,metrics))
  
  write.table(vars$results[nrow(vars$results),],"measurements.txt")
  setwd("..")
  
  return(vars)
}
tmp <- lpc_builds(vars,network,target = target)
vars$node_names
