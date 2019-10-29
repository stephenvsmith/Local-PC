### Script for building on many datasets with local pc
computer_name <- Sys.info()["nodename"]
my_laptop <- computer_name == "Stephens-MacBook-Pro.local"

library(bnlearn)
if (!my_laptop){
  setwd("/home/stephen/Dropbox/Academics/Research/Desktop/Research/package")
  install("localpc")
}
library(localpc)

path_start <- ifelse(my_laptop,"/Users/stephensmith/","/home/stephen/")
ps2 <- ifelse(my_laptop,"/Users/stephensmith/","/home/stephen/Dropbox/Academics/Research/")

# Grid for generating data
data.grid = data.frame(network = "asia",
                       data.type = "continuous",
                       n.dat = 1,
                       n.obs = 1000,
                       c.ratio = 0,
                       max.in.degree = Inf,
                       lb = 0,  # lower bound of coefficients
                       ub = 10,  # upper bound of coefficients
                       low = 0,  # lower bound of variances if continuous, of number of levels if discrete
                       high = 5,  # upper bound of variances if continuous, of number of levels if discrete
                       scale = TRUE,
                       stringsAsFactors = FALSE)

# This function removes all variables that aren't used
simplify_dag <- function(mat){
  zeros <- rep(0,nrow(mat))
  rm_ind <- c()
  for (i in 1:nrow(mat)){
    if (all(mat[i,]==zeros) & all(mat[,i]==zeros))
      rm_ind <- c(rm_ind,i)
  }
  return(rm_ind)
}

# Loading Function to calculate HD and SHD
source(paste0(ps2,"Desktop/Research/package/test_scripts/local_pc_dist.R"))
# Loading Data Generation Function
source(paste0(ps2,"Desktop/Research/projects/code/data_gen_R.R"))

# Obtain the names of all the networks
file_vec <- list.files(paste0(ps2,"Desktop/Research/projects/bn_data_generation/networks/rds/"))
nets <- sapply(file_vec, function(x) sub(".rds","",x))

# Function to visualize true DAG in neighborhood against local PC DAG
source(paste0(ps2,"Desktop/Research/package/test_scripts/test_bn.R"))

# Storing Results
results <- matrix(nrow=0,ncol = 6)
colnames(results) <- c("Network","Observations","Target","HD","SHD","Same_V_Structures")

# Set up result file
if (!dir.exists(paste0(path_start,"Dropbox/Academics/Research/Results/Sample")))
  dir.create( paste0(path_start,"Dropbox/Academics/Research/Results/Sample") )

for (net in nets){
  
  if (net == "andes" | net == "alarm") next
  
  # Store results for current network
  current_network <- matrix(nrow=0,ncol = 6)
  colnames(current_network) <- c("Network","Observations","Target","HD","SHD","Same_V_Structures")
  # Change Working Directory to Store Results
  setwd(paste0(path_start,"Dropbox/Academics/Research/Results/Sample"))
  # Get the true DAG
  data.grid$network <- net
  
  # Set up File Storage for Data and Results
  if (!file.exists(paste0(path_start,"Dropbox/Academics/Research/Data/",net)))
    dir.create(paste0(path_start,"Dropbox/Academics/Research/Data/",net))
  
  # Generate Data and Get True DAG
  gdg <- generate.data.grid(data.grid,out.dir=paste0(path_start,"Dropbox/Academics/Research/Data/",net),verbose=FALSE,path.start=ps2)
  dir.name <- paste0(net,"; ","n = 1000; c = 0")
  true_dag <- as.matrix(read.table(paste0(path_start,"Dropbox/Academics/Research/Data/",net,"/",dir.name,"/trueDAG.txt")))
  data <- as.matrix(read.table(paste0(path_start,"Dropbox/Academics/Research/Data/",net,"/",dir.name,"/data1.txt")))
  
  # Create a File for the Network
  if (!file.exists(paste0("./",net))) dir.create(paste0("./",net))
  setwd(paste0("./",net))
  
  # Save File of True DAG
  file <- paste0(ps2,"Desktop/Research/projects/bn_data_generation/networks/rds/",net,".rds")
  names <- names(readRDS(file))
  t_dag <- empty.graph(names)
  rownames(true_dag) <- names
  colnames(true_dag) <- names
  amat(t_dag) <- true_dag
  png(file = "true_dag.png")
  par(mfrow=c(1,1))
  bnlearn::graphviz.plot(t_dag,main = "DAG")
  dev.off()
  
  # Save File of True CPDAG
  tcp_dag <- cpdag(t_dag)
  png(file = "true_cp_dag.png")
  par(mfrow=c(1,1))
  bnlearn::graphviz.plot(t_dag,main = "CPDAG")
  dev.off()
  
  # Get all results for one target
  if (!file.exists("./one_target")) dir.create("./one_target")
  setwd("./one_target")
  for (i in 1:nrow(true_dag)) {
    target <- i
    if (!file.exists(paste0(net,"; target=",target))) dir.create(paste0(net,"; target=",target))
    setwd(paste0("./",net,"; target=",target))
    temp <- test_bn(net,target,true_dag=true_dag,data = data,path_start=path_start)
    cat(paste(temp,collapse = '\n'),file = "p_vals.txt")
    res <- local_pc_dist(target,true_dag=true_dag,data = data)
    results <- rbind(results,c(net,data.grid$n.obs,target,res))
    current_network <- rbind(current_network,c(net,data.grid$n.obs,target,res))
    write.table(results[nrow(results),],"measurements.txt")
    
    setwd("..")
  }
  
  # Get all results for two targets
  setwd("..")
  if (!file.exists("./two_targets")) dir.create("./two_targets")
  setwd("./two_targets")
  nodes <- 1:nrow(true_dag)
  target_sets <- combn(nodes,2)
  num <- min(30,ncol(target_sets))
  target_sets <- target_sets[,sample(1:ncol(target_sets),num,replace = FALSE)]
  tmp <- apply(target_sets,2,function(target){
    if (!file.exists(paste0(net,"; target=",paste(target,collapse = ",")))) 
      dir.create(paste0(net,"; target=",paste(target,collapse = ",")))
    
    setwd(paste0("./",net,"; target=",paste(target,collapse = ",")))
    temp <- test_bn(net,target,true_dag,path_start=path_start,data = data)
    cat(paste(temp,collapse = '\n'),file = "p_vals.txt")
    
    res <- local_pc_dist(target,true_dag,data = data)
    results <<- rbind(results,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    current_network <<- rbind(current_network,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    
    write.table(results[nrow(results),],"measurements.txt")
    setwd("..")
  })
  
  
  # Get all results for three targets
  setwd("..")
  if (!file.exists("./three_targets")) dir.create("./three_targets")
  setwd("./three_targets")
  nodes <- 1:nrow(true_dag)
  target_sets <- combn(nodes,3)
  num <- min(30,ncol(target_sets))
  target_sets <- target_sets[,sample(1:ncol(target_sets),num,replace = FALSE)]
  tmp <- apply(target_sets,2,function(target){
    if (!file.exists(paste0(net,"; target=",paste(target,collapse = ",")))) 
      dir.create(paste0(net,"; target=",paste(target,collapse = ",")))
    
    setwd(paste0("./",net,"; target=",paste(target,collapse = ",")))
    temp <- test_bn(net,target,true_dag,path_start=path_start,data = data)
    cat(paste(temp,collapse = '\n'),file = "p_vals.txt")
    res <- local_pc_dist(target,true_dag,data = data)
    results <<- rbind(results,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    current_network <<- rbind(current_network,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    
    write.table(results[nrow(results),],"measurements.txt")
    setwd("..")
  })
  
  # Get all results for four targets
  setwd("..")
  if (!file.exists("./four_targets")) dir.create("./four_targets")
  setwd("./four_targets")
  nodes <- 1:nrow(true_dag)
  target_sets <- combn(nodes,4)
  num <- min(30,ncol(target_sets))
  target_sets <- target_sets[,sample(1:ncol(target_sets),num,replace = FALSE)]
  tmp <- apply(target_sets,2,function(target){
    if (!file.exists(paste0(net,"; target=",paste(target,collapse = ",")))) 
      dir.create(paste0(net,"; target=",paste(target,collapse = ",")))
    
    setwd(paste0("./",net,"; target=",paste(target,collapse = ",")))
    temp <- test_bn(net,target,true_dag,path_start=path_start,data = data)
    cat(paste(temp,collapse = '\n'),file = "p_vals.txt")
    res <- local_pc_dist(target,true_dag,data = data)
    results <<- rbind(results,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    current_network <<- rbind(current_network,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    
    write.table(results[nrow(results),],"measurements.txt")
    setwd("..")
  })
  
  # Get all results for five targets
  setwd("..")
  if (!file.exists("./five_targets")) dir.create("./five_targets")
  setwd("./five_targets")
  nodes <- 1:nrow(true_dag)
  target_sets <- combn(nodes,5)
  num <- min(30,ncol(target_sets))
  target_sets <- target_sets[,sample(1:ncol(target_sets),num,replace = FALSE)]
  tmp <- apply(target_sets,2,function(target){
    if (!file.exists(paste0(net,"; target=",paste(target,collapse = ",")))) 
      dir.create(paste0(net,"; target=",paste(target,collapse = ",")))
    
    setwd(paste0("./",net,"; target=",paste(target,collapse = ",")))
    temp <- test_bn(net,target,true_dag,path_start=path_start,data = data)
    cat(paste(temp,collapse = '\n'),file = "p_vals.txt")
    res <- local_pc_dist(target,true_dag,data = data)
    results <<- rbind(results,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    current_network <<- rbind(current_network,c(net,data.grid$n.obs,paste(target,collapse = ","),res))
    
    write.table(results[nrow(results),],"measurements.txt")
    setwd("..")
  })
  
  
  
  setwd("..")
  
  write.table(current_network,paste0(net,"_results.txt"))
  
}

setwd("..")
write.table(results,file = "finalResults.txt")
