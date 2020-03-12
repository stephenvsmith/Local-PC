##############################################################################################################################

# Author: Stephen Smith
# Date: 1/15/2020
# Most Recent Update: 2/20/2020
# Description: In this file, we will set up simulations for the Local PC Algorithm (After adjusting it for local FCI)

##############################################################################################################################

set.seed(1000)
library(tidyverse,quietly = TRUE)
library(bnlearn,quietly = TRUE)
library(pcalg,quietly = TRUE)
library(parallel,quietly = TRUE)

# Setup -------------------------------------------------------------------

vars <- list() # List for storing important variables
alpha_vals <- c(0.01,0.001,0.0001)
hoffman <- FALSE


# population or sample version
vars$pop <- FALSE 
if (vars$pop) vars$lmax <- 4

pop_string <- ifelse(vars$pop,"population","sample")
vars$data.grid <- data.frame(network = "asia",
                       data.type = "continuous",
                       n.dat = 1,
                       n.obs = 1000,
                       c.ratio = 0,
                       max.in.degree = Inf,
                       lb = 0.1,  # lower bound of coefficients
                       ub = 5,  # upper bound of coefficients
                       low = 0.1,  # lower bound of variances if continuous, of number of levels if discrete
                       high = 5,  # upper bound of variances if continuous, of number of levels if discrete
                       scale = TRUE,
                       stringsAsFactors = FALSE)
if (vars$pop){
  vars$alpha_vals <- "population"
} else {
  vars$alpha_vals <- alpha_vals
}

cat("The significance levels for conditional independence tests are\n",paste(vars$alpha_vals,collapse = '\n'),"\n",sep = "")
vars$networks_to_skip <- c("munin4","pigs","pathfinder")

# Define paths for storing
vars$home_dir <- '/home/stephen'
#vars$home_dir <- '/Users/stephensmith'
vars$package_loc <- paste0(vars$home_dir,'/Dropbox/Academics/Research/Code/Packages/LocalPC')
vars$data_gen_file <- paste0(vars$home_dir,'/Dropbox/Academics/Research/Code/Scripts/Current Script/data_gen.R')
vars$result_dir <- paste0(vars$home_dir,'/Research_Results/Current/Results/',pop_string,'/')
if (!vars$pop & !hoffman) {
  vars$result_dirs <- sapply(vars$alpha_vals,function(a){
    result_dir <- paste0(vars$result_dir,'alpha=',a)
    if (!dir.exists(vars$result_dir)){
      dir.create(vars$result_dir)
    }
    return(result_dir)
  })
  names(vars$result_dirs) <- paste0("alpha=",vars$alpha_vals)
}
vars$simulation_dir <- paste0(vars$home_dir,'/Research_Results/Current/Simulations/')
vars$rds_dir <- paste0(vars$home_dir,'/Dropbox/Academics/Research/Code/Simulation/Networks/rds')

if (hoffman){
  source('/u/scratch/s/stephens/Hoffman/hoffman_dirs.R')
  vars <- hoffman_dirs(vars)
}

source(vars$data_gen_file)

vars$metric_names <- c("Algorithm","Network","Network Size","Max Sep Set",
  "Target","Neighborhood Size","Ground Truth","Number of Tests",
  "Time to Build","n","alpha","Neighborhood Edges",
  "Undirected Missing","Directed Missing","Total Missing",
  "Undirected Added","Directed Added","Total Added",
  "Wrong_Direction","Directed_Undirected","Undirected_Directed","My TP",
  "SHD","tp","fp","fn","Same_V_Structures")

devtools::install(vars$package_loc,quiet = TRUE)
library(localpc,quietly = TRUE)

numCores <- min(detectCores() - 1,4)
numCores <- max(1,numCores)
cat("We will be using",numCores,"cores.\n\n")


# Grab Networks and Network Info -------------------------------------
create_specs <- function(vars){

  # Grab all networks and network names/sizes
  vars <- get_networks(vars)

  # Grab True DAGs and CPDAGs
  vars <- get_dags_cpdags(vars)

  # Get Targets and Their Neighbors
  vars <- get_network_targets(vars)
  vars <- set_lmax_tracker(vars)
  vars <- get_vertex_sets(vars)
  
  # Remove unwanted networks
  vars$net_names <- setdiff(vars$net_names,vars$networks_to_skip) 

  return(vars)
}

# Create list with networks and network info (size,names)
get_networks <- function(vars){
  setwd(vars$rds_dir)
  network_files <- list.files()
  vars$networks <- lapply(network_files,readRDS)
  
  vars$net_names <- stringr::str_replace(network_files,".rds","")
  names(vars$networks) <- vars$net_names
  
  vars$num_nodes <- lapply(vars$networks,function(x) length(names(x)))
  names(vars$num_nodes) <- vars$net_names
  return(vars)
}

# Create list with DAGs and CPDAGs
get_dags_cpdags <- function(vars){
  vars$true_dags <- lapply(vars$networks,amat)
  names(vars$true_dags) <- vars$net_names
  
  vars$true_cpdags <- lapply(vars$networks,cpdag)
  vars$true_cpdags_amat <- lapply(vars$true_cpdags,amat)
  names(vars$true_cpdags) <- vars$net_names
  names(vars$true_cpdags_amat) <- vars$net_names
  return(vars)
}

# Generate the targets for each network
get_network_targets <- function(vars){
  vars$targets <- lapply(vars$num_nodes,get_num_targets)
  names(vars$targets) <- vars$net_names
  
  return(vars)
}

get_num_targets <- function(p){
  # minimize the number of potential targets
  set.seed(10000)
  if (p>50 & p<500){
    targets <- sample(1:p,50)
  } else if (p > 500){
    targets <- sample(1:p,10)
  } else {
    targets <- 1:p
  }
  return(targets)
}

# Set the tracker for lmax counting
set_lmax_tracker <- function(vars){
  lmax_tracker <- lapply(vars$net_names,function(net){
    res <- lapply(vars$targets[[net]],function(t){
      return(NULL)
    })
    names(res) <- as.character(vars$targets[[net]])
    return(res)
  })
  names(lmax_tracker) <- vars$net_names
  vars$lmax_tracker <- lmax_tracker
  return(vars)
}

# Create Initial Graphs for each network and target node
# First get vertex sets for each target containing the neighborhood of target
get_vertex_sets <- function(vars){
  vars$vertex_sets <- lapply(vars$net_names,function(net){
    targets <- vars$targets[[net]]
    true_dag <- vars$true_dags[[net]]
    vertex_sets <- lapply(targets,function(t){
      neighbors <- find_neighbors(t,true_dag)
      return(c(t,neighbors))
    })
    names(vertex_sets) <- as.character(targets)
    return(vertex_sets)
  })
  names(vars$vertex_sets) <- vars$net_names
  
  return(vars)
}

find_neighbors <- function(t,true_dag){
  # Parents of node t
  parents <- which(true_dag[,t]==1)
  # Children of node t
  children <- which(true_dag[t,]==1)
  # Spouses share the same children
  spouses <- c()
  for (c in children){ # loop through all of t's children
    # Find all parents of current child
    potential <- which(true_dag[,c]==1)
    # Add all parents that are not t to spouses
    spouses <- c(spouses,potential[potential!=t])
  }
  return(unique(c(parents,children,spouses)))
}

# Build Network Simulations -----------------------------------------------

create_sims <- function(vars){
  cat("Generating Datasets\n\n")
  vars$n_vals <- lapply(vars$net_names,function(net){
    vars <- network_information(net,vars)
    # checks whether or not we already have simulations built
    vars <- check_sims_created(vars)
    n_vals <- sort(unique(c(ceiling(vars$p/4),floor(vars$p/2),25,50,100,1000)),decreasing = TRUE)
    n_vals <- n_vals[n_vals >= 25]
    
    # Generate Data
    if (vars$sims_not_created){
      for (n in n_vals){
        sim_net_dir <- paste0(vars$simulation_dir,net)
        vars$data.grid$n.obs <- n
        gdg <- generate.data.grid(vars$data.grid,
                                  out.dir=sim_net_dir,
                                  verbose=FALSE,
                                  path.start=vars$home_dir,
                                  var_list=vars)
      }
    }
    return(n_vals)
  })
  names(vars$n_vals) <- vars$net_names
  
  return(vars)
  
}

# Check whether or not simulations have been created so they do not have to be created again
check_sims_created <- function(vars){
  files <- list.files(vars$simulation_dir)
  vars$sims_not_created <- length(files) == 0 | !(vars$net_current %in% files)
  return(vars)
}

# Track current important network information
network_information <- function(net,vars){
  p <- vars$num_nodes[[net]]
  vars$node_names <- names(vars$networks[[net]])
  vars$data.grid$network <- net
  vars$net_current <- net
  vars$p <- p
  return(vars)
}

# Main Functions ----------------------------------------------------------

simulate_local_pc <- function(vars){
  
  cat("We are beginning our simulations at",date(),'\n\n')

  # Generate Simulation Data in Output Directory
  vars <- create_specs(vars)

  if (!vars$pop){
    vars <- create_sims(vars)
  }

  # Run various simulation for different tolerance levels
  all_results <- lapply(vars$alpha_vals,function(alpha){
    completed_tol_level <- alpha %in% is_alpha_completed(vars)
    if (!completed_tol_level){
      vars$alpha <- alpha
      go_to_dir(vars$result_dirs[paste0("alpha=",alpha)])
      # Run algorithms on all datasets for given tolerance level
      local_pc_results <- mclapply(vars$net_names,run_pc,vars,mc.cores = numCores,mc.preschedule = FALSE)
      go_to_dir(vars$result_dirs[paste0("alpha=",alpha)])
      saveRDS(local_pc_results,paste0("results_alpha_",alpha,".rds"))
      return(local_pc_results)
    } else {
      cat("We have already completed tolerance level",alpha,"and we are moving on.\n")
      return(get_alpha_results(vars,alpha))
    }
    
  })

  setwd(vars$result_dir)
  saveRDS(all_results,"AllResults.rds")
  
  cat("We are ending our simulations at",date(),'\n\n')

  return(local_pc_results)

}

is_alpha_completed <- function(vars){
  current_wd <- getwd()
  
  # First, find alpha values that have been started
  alpha_existence <- sapply(vars$alpha_vals,function(alpha){
    existence <- dir.exists(vars$result_dirs[paste0("alpha=",alpha)])
    return(existence)
  })
  # This is the vector of alpha values that have a folder
  remaining_alphas <- vars$alpha_vals[alpha_existence]
  
  if (length(remaining_alphas)==0) {
    return(c())
  } else {
    # Second, find the networks within each alpha value that have been completed
    completed_nets <- lapply(remaining_alphas,function(alpha){
      setwd(vars$result_dirs[paste0("alpha=",alpha)])
      net_folders <- list.files()
      nets_to_complete <- setdiff(vars$net_names,vars$networks_to_skip)
      completed_nets_logical <- sapply(nets_to_complete,function(net){
        complete <- file.exists(paste0(net,"/results.rds"))
        return(complete)
      })
      nets_completed <- nets_to_complete[completed_nets_logical]
    })
    
    # Determine if all the networks for this alpha value are completed
    completed_alphas_logical <- sapply(completed_nets,function(nets){
      nets_to_complete <- setdiff(vars$net_names,vars$networks_to_skip)
      res <- all(nets_to_complete %in% nets)
      return(res)
    })
    
    setwd(current_wd)
    
    completed_alphas <- remaining_alphas[completed_alphas_logical]
    
    return(completed_alphas)
  }
}

get_alpha_results <- function(vars,alpha){
  current_wd <- getwd()
  cat("Grabbing results for alpha = ",alpha,"\n")
  setwd(vars$result_dirs[paste0("alpha=",alpha)])
  nets <- list.files()
  nets <- setdiff(nets,vars$networks_to_skip)
  nets <- nets[dir.exists(nets)]
  alpha_results <- lapply(nets,function(net){
    res <- readRDS(paste0(net,"/results.rds"))
    return(res)
  })

  setwd(current_wd)
  return(alpha_results)
}

# The wrapper function for each network's results
run_pc <- function(net,vars){
  
  skip_net <- is_net_completed(vars,net)
  
  if (!skip_net){
    net_start <- Sys.time()
    cat("We are working on network",net,"for tolerance level",vars$alpha,"\n")
    
    # Build directories for the results for each target
    new_location <- build_net_directory(net,targets,vars)
    # puts us in the targets folder for this network
    setwd(new_location)
    
    # Current Network Information
    vars <- network_information(net,vars)
    
    # Run our algorithms on all of our targets
    results <- lapply(vars$targets[[net]],run_pc_target,vars=vars)
    results <- do.call(rbind,results)
    
    saveRDS(results,"../results.rds")
    net_end <- Sys.time()
    diff <- net_end - net_start
    units(diff) <- "mins"
    cat("Network",net,"took",diff,"minutes to complete\n")
    
    go_to_dir(vars$result_dirs[paste0("alpha=",vars$alpha)])
  } else {
    # Grab the results from the rds file
    results <- grab_net_results(vars,net)
    cat(net,"has already been completed for alpha=",vars$alpha,"so we are moving on.\n")
  }
  
  

  return(results)
}

# checks to see if we have completed this network already
is_net_completed <- function(vars,net){
  
  current_wd <- getwd()
  
  setwd(vars$result_dirs[paste0("alpha=",vars$alpha)])
  net_folders <- list.files()
  
  # Check to see if the network folder is created
  folder_created <- net %in% net_folders
  if (folder_created){
    results_finished <- file.exists(paste0(net,"/results.rds"))
    setwd(current_wd)
    return(results_finished)
  } else {
    setwd(current_wd)
    return(folder_created)
  }
}

# Grabs the rds file for the tolerance level and the network
grab_net_results <- function(vars,net) {
  current_wd <- getwd()
  
  setwd(vars$result_dirs[paste0("alpha=",vars$alpha)])
  results <- readRDS(paste0("./",net,"/results.rds"))
  
  setwd(current_wd)
  return(results)
}


# Creates network directory for tolerance level to store results
# Also stores picture of DAG
# returns directory to store info about targets
build_net_directory <- function(net,targets,vars){
  
  go_to_dir(net)

  # Save photo of the true dag
  png(file = "true_dag.png")
  bnlearn::graphviz.plot(vars$networks[[net]],main = "True DAG")
  dev.off()

  if (!dir.exists("targets")){
    dir.create("targets")
  }

  return("targets")
}


# Target Functions --------------------------------------------------------

# Main function run for each target
run_pc_target <- function(t,vars){
  net <- vars$net_current
  q1 <- format(paste("Tolerance:",vars$alpha),width = 20,justify = "centre")
  q2 <- format(paste("Net:",net),width = 25,justify = "centre")
  q3 <- format(paste("Target:",t),width = 15,justify = "centre")
  cat(q1,"|",q2,"|",q3,"\n")
  vars <- target_setup(t,vars)

  if (vars$pop){
    vars$run_pop_target(t,vars)
    result_mat <- data.frame(vars$metrics)
  } else {
    result_mat <- lapply(vars$n_vals[[net]],function(n){
      # Create directory for sample size
      vars$n <- n
      go_to_dir(paste0("n=",n))
      df <- grab_data(net,n,vars)

      # Run global PC
      vars <- run_global_pc(df,vars)
      # Run local PC
      vars <- run_local_pc(df,vars)
      
      rownames(vars$metrics) <- NULL
      return(vars$metrics)
    })
    
    result_mat <- do.call(rbind,result_mat)

    result_mat <- data.frame(result_mat)
    write.table(result_mat,"./target_results.txt")
    setwd('..')
  }

  return(result_mat)
}

target_setup <- function(t,vars){
  # Change working directory to the target's folder
  go_to_dir(paste0("Target=",t))
  
  vars$t <- t
  # Store vertex sets neighbors
  net <- vars$net_current
  vars$V <- vars$vertex_sets[[net]][[as.character(t)]]
  
  # Create the initial graph for the local pc algorithm
  G_start <- get_initial_graph(vars)
  
  # Store True DAG, True CPDAG, and True CPDAG adj. mat.
  vars$G_start <- G_start
  td <- vars$true_dags[[net]]
  vars$td <- td
  vars$tcpdag <- vars$true_cpdags[[net]]
  vars$tcpdag_amat <- vars$true_cpdags_amat[[net]]
  return(vars)
}

go_to_dir <- function(file){
  if (!dir.exists(file)){
    dir.create(file)
  }
  setwd(file)
}

get_initial_graph <- function(vars){
  V <- vars$V
  p <- vars$p
  node_names <- vars$node_names
  g <- matrix(0,nrow = p,ncol = p)
  rownames(g) <- node_names
  colnames(g) <- node_names

  # Create the initial connected graph to pass into the PC algorithm
  for (i in V){
    remaining <- setdiff(V,i)
    for (j in remaining){
      g[i,j] <- 1
    }
  }

  return(g)
}

# Run the simulation for the population version at the current target
run_pop_target <- function(t,vars){
  start <- Sys.time()
  vars$localpc_result <- localpc(true_dag = vars$td,target = t,G = vars$G_start,lmax = vars$lmax,
                                 verbose = FALSE, verbose_small = FALSE)
  end <- Sys.time()
  diff <- end - start
  units(diff) <- "mins"
  vars$diff <- diff
  vars$num_tests <- vars$localpc_result$num_tests
  # After this step, working directory is in the network folder
  vars <- neighborhood_results(vars,vars$localpc_result$G,localpc=TRUE) 
  
  return(vars)
}

# Grab dataframe for network at given sample size
grab_data <- function(net,n,vars){
  sim_net_dir <- paste0(vars$simulation_dir,net)
  dir_name <- paste0(net,"; ","n = ",n,"; c = 0")
  df <- read.table(paste0(sim_net_dir,'/',dir_name,"/data1.txt"))
  return(df)
}

# Run global PC 
run_global_pc <- function(df,vars){
  
  largest_possible_sepset <- get_lmax_from_tracker(vars)
  sink(file = "./pc_tests.txt")
  start <- Sys.time()
  vars$pc.fit <- as(pc(suffStat = list(C = cor(df), n = vars$n),
                       indepTest = gaussCItest, ## indep.test: partial correlations
                       alpha=vars$alpha, labels = vars$node_names,
                       verbose = TRUE,m.max=largest_possible_sepset),"amat")
  end <- Sys.time()
  sink(file = NULL)
  diff <- end - start
  units(diff) <- "mins"
  vars$diff <- diff
  vars$lmax <- get_lmax('./pc_tests.txt')
  vars <- update_lmax_tracker(vars)
  vars$num_tests <- get_pc_test_num('./pc_tests.txt')
  # After this step, working directory is in the network folder
  vars <- neighborhood_results(vars,vars$pc.fit,localpc = FALSE)
  return(vars)
}

# Get largest possible sep set size
get_lmax_from_tracker <- function(vars){
  tracker_location <- paste0(vars$result_dir,"/lmax_tracker.rds")
  if (file.exists(tracker_location)){
    vars$lmax_tracker <- readRDS(tracker_location)
  }
  current_lmax <- vars$lmax_tracker[[vars$net_current]][[as.character(vars$t)]]
  if (is.null(current_lmax)){
    return(5)
  }
  return(current_lmax)
}

# Update the lmax tracker
update_lmax_tracker <- function(vars){
  
  tracker_location <- paste0(vars$result_dir,"/lmax_tracker.rds")
  if (file.exists(tracker_location)){
    vars$lmax_tracker <- readRDS(tracker_location)
  }
  
  current_lmax <- vars$lmax_tracker[[vars$net_current]][[as.character(vars$t)]]
  

  if (is.null(current_lmax)){
    vars$lmax_tracker[[vars$net_current]][[as.character(vars$t)]] <- vars$lmax
  }
  
  tracker_location <- paste0(vars$result_dir,"/lmax_tracker.rds")
  saveRDS(vars$lmax_tracker,tracker_location)
  return(vars)
}

# Run Local PC
run_local_pc <- function(df,vars){
  # Run Local PC
  
  val <- get_lmax_from_tracker(vars)
  vars$lmax <- max(1,val)
  start <- Sys.time()
  vars$localpc_result <- localpc(true_dag = vars$td,data = df,target = vars$t,G = vars$G_start,lmax = vars$lmax,
                                 verbose = FALSE, verbose_small = FALSE,pop = FALSE)
  end <- Sys.time()
  diff <- end - start
  units(diff) <- "mins"
  vars$diff <- diff
  vars$num_tests <- vars$localpc_result$num_tests
  vars <- neighborhood_results(vars,vars$localpc_result$G,localpc=TRUE)
  return(vars)
}

# Get the maximum size of separating sets used in PC algorithm
get_lmax <- function(file){
  pc_file <- read_file(file = file)
  
  sep_sets <- unlist(str_extract_all(pc_file,'S= .* :'))
  counts <- sapply(sep_sets,function(s){
    res <- str_count(s,"[0-9]+")
  })
  
  if (is.list(counts)){
    counts <- unlist(counts)
    if (length(counts)==0){
      counts <- 0
    }
  }
  
  return(max(counts))
}

get_pc_test_num <- function(file){
  pc_file <- read_file(file = file)
  tests <- unlist(str_extract_all(pc_file,'pval = [0-9]+'))
  return(length(tests))
}


# Compile all results about the simulation
neighborhood_results <- function(vars,estimated_amat,localpc=TRUE){
  
  vars$algorithm <- ifelse(localpc,"Local PC","PC")
  
  # Zoom in on estimated and true DAGs (only the target and first-order neighbors)
  vars$nodes_zoom <- vars$node_names[vars$V]

  # Get bnlearn object for estimated subgraph
  vars$estimated_dag_subgraph <- get_bnlearn_graph(estimated_amat,vars)
  
  # Save global PC estimate for comparison with local PC
  if (!localpc){
    vars$estimated_globalpc <- vars$estimated_dag_subgraph
  }
  
  # Get bnlearn object for local ground truth
  true_dag_subgraph <- get_bnlearn_graph(vars$td,vars)
  # CPDAG of subgraph
  subgraph_cpdag <- cpdag(true_dag_subgraph)
  vars$subgraph_cpdag <- subgraph_cpdag 
  
  # Get bnlearn object for global ground truth
  cpdag_subgraph <- get_bnlearn_graph(vars$tcpdag_amat,vars)
  vars$cpdag_subgraph <- cpdag_subgraph

  # Compare results
  vars <- combine_results(vars,localpc)
  
  return(vars)
}

# Return bnlearn object for adjacency matrix given
get_bnlearn_graph <- function(adj.mat,vars){
  bnlearn_graph <- empty.graph(nodes = vars$nodes_zoom)
  if (length(vars$V) > 1){
    amat(bnlearn_graph,check.cycles=FALSE) <- adj.mat[vars$V,vars$V]
  }
  return(bnlearn_graph)
}

# Combine all important metric results
combine_results <- function(vars,localpc){
  
  if (vars$pop) {
    vars$n <- "population"
  }
  
  # Find distance between estimated DAG and Local Ground Truth
  results_est_local <- graph_dist(vars$subgraph_cpdag,vars$estimated_dag_subgraph)
  metrics <- get_metrics_vector(results_est_local,vars,
                                ground_truth="Local",algo=vars$algorithm)
  
  # Distance between estimated DAG and Global Ground Truth
  results_est_global <- graph_dist(vars$cpdag_subgraph,vars$estimated_dag_subgraph)
  metrics <- rbind(metrics,get_metrics_vector(results_est_global,vars,
                                              ground_truth = "Global",
                                              algo = vars$algorithm))
  
  if (localpc){
    results_local_global <- graph_dist(vars$cpdag_subgraph,vars$subgraph_cpdag)
    metrics <- rbind(metrics,get_metrics_vector(results_est_global,vars,
                                                ground_truth = "Global",
                                                algo = "Local Ground Truth"))
  }
  
  colnames(metrics) <- vars$metric_names
  if (localpc){
    vars$metrics <- rbind(vars$metrics,metrics)
  } else {
    vars$metrics <- metrics
  }
  
  # Store results and create photos
  if (vars$pop){
    next_location <- save_results(vars,file = paste0("Results.txt"),vars$subgraph_cpdag,vars$estimated_dag_subgraph)
    setwd(next_location)
  } else {
    if (!localpc){
      next_location <- save_results(vars,file = paste0("PC Results.txt"),vars$cpdag_subgraph,vars$estimated_dag_subgraph)
    } else {
      next_location <- save_results(vars,file = paste0("Local PC Results.txt"),vars$subgraph_cpdag,vars$estimated_dag_subgraph)
    }
    setwd(next_location)
  }
  
  return(vars)
}

# Get metrics vector based on results 
get_metrics_vector <- function(results,vars,ground_truth,algo){
  
  if (algo=='Local Ground Truth'){
    metrics <- c("Algorithm"=algo,"Network"=vars$net_current,
                 "Network Size"=vars$p,"Max Sep Set"=NA,
                 "Target"=vars$t,"Neighborhood Size"=length(vars$V),
                 "Ground Truth"=ground_truth,"Number of Tests"=NA,
                 "Time to Build"=NA,"n"=NA,"alpha"=NA,results)
  } else {
    metrics <- c("Algorithm"=algo,"Network"=vars$net_current,
               "Network Size"=vars$p,"Max Sep Set"=vars$lmax,
               "Target"=vars$t,"Neighborhood Size"=length(vars$V),
               "Ground Truth"=ground_truth,"Number of Tests"=vars$num_tests,
               "Time to Build"=vars$diff,"n"=vars$n,"alpha"=vars$alpha,results)
  }
  return(metrics)
}

save_results <- function(vars,file,truth,estimate){
  
  results_txt_file(vars,file)

  truth_string <- ifelse(vars$algorithm=='PC','Global Ground Truth','Local Ground Truth')
  png(filename = paste(vars$algorithm,"comparison.png"),width = 1200,height = 1200)
  par(mfrow=c(1,2))
  bnlearn::graphviz.plot(estimate,main = vars$algorithm,highlight = list(nodes=vars$node_names[vars$t],col="green",fill="light blue",textCol="white"))
  bnlearn::graphviz.plot(truth,main = truth_string,highlight = list(nodes=vars$node_names[vars$t],col="green",fill="light blue",textCol="white"))
  dev.off()

  if (vars$algorithm=='Local PC' & !vars$pop){
    png(filename = paste("estimates comparison.png"),width = 1200,height = 1200)
    par(mfrow=c(1,2))
    bnlearn::graphviz.plot(estimate,main = vars$algorithm,highlight = list(nodes=vars$node_names[vars$t],col="green",fill="light blue",textCol="white"))
    bnlearn::graphviz.plot(vars$estimated_globalpc,main = "PC",highlight = list(nodes=vars$node_names[vars$t],col="green",fill="light blue",textCol="white"))
    dev.off()
    
    png(filename = paste("ground truth comparison.png"),width = 1200,height = 1200)
    par(mfrow=c(1,2))
    bnlearn::graphviz.plot(truth,main = "Local Ground Truth",highlight = list(nodes=vars$node_names[vars$t],col="green",fill="light blue",textCol="white"))
    bnlearn::graphviz.plot(vars$cpdag_subgraph,main = "Global Ground Truth",highlight = list(nodes=vars$node_names[vars$t],col="green",fill="light blue",textCol="white"))
    dev.off()
  } 

  if (!file.exists("target_in_DAG.png")){
    png(filename = "target_in_DAG.png")
    par(mfrow=c(1,1))
    bnlearn::graphviz.plot(vars$networks[[vars$net_current]],
                           highlight =list(nodes=vars$node_names[vars$t],
                                           col="green",fill="light blue",textCol="white"))
    dev.off()
  }

  next_location <- ifelse(vars$pop | vars$algorithm!='PC','..','.')

  return(next_location)
}

results_txt_file <- function(vars,file){
  if (vars$pop | vars$algorithm=='PC'){
    cols_print <- 1:2
  } else {
    cols_print <- 3:4
  }

  cat("",file = file,append = FALSE)
  for (i in 1:length(vars$metric_names)){
    q1 <- format(colnames(vars$metrics)[i],width = 25,justify = "left")
    q2 <- format(vars$metrics[cols_print[1],i],width = 25,justify = "centre")
    q3 <- format(vars$metrics[cols_print[2],i],width = 25,justify = "right")
    cat(q1,q2,q3,"\n",file = file,append = TRUE)
  }
}

  

########## Distance Metrics ##########

### Here we will define functions to compare the DAGs obtained from the algorithm
### and the true CPDAG

# function that determines whether or not a pair has already been considered (True if it has not been considered, False otherwise)
checkPair <- function(pair_list,i,j){
  if (length(pair_list)==0) return(TRUE)
  for (k in 1:length(pair_list)){
    pair_considered <- pair_list[[k]]
    if (i %in% pair_considered & j %in% pair_considered) return(FALSE)
  }
  # pair has not been considered otherwise
  return(TRUE)
}

# This will compute the number of extraneous edges in our estimated CPDAG
addedEdges <- function(true_cpdag,lpc_cpdag){

  if (is.null(true_cpdag) & is.null(lpc_cpdag)){
    return(rep(0,3))
  }
  total <- 0
  directed <- 0
  undirected <- 0
  pairs_counted <- list()
  for (i in 1:nrow(lpc_cpdag)){
    for (j in 1:ncol(lpc_cpdag)){

      if (lpc_cpdag[i,j] == 1){
        if (true_cpdag[i,j]==0 & true_cpdag[j,i]==0){

          # Case where we have added an edge that doesn't exist at all
          if (checkPair(pairs_counted,i,j)){
            total <- total + 1
            if (lpc_cpdag[j,i]==1){
              undirected <- undirected + 1
            } else {
              directed <- directed + 1
            }

            pairs_counted[[total]] <- c(i,j)
          }

        }
      }

    }
  }
  return(c("total"=total,"undirected"=undirected,"directed"=directed))
}

# This will compute the number of missing edges in the estimated DAG
missingEdges <- function(true_cpdag,lpc_cpdag){
  if (is.null(true_cpdag) & is.null(lpc_cpdag)){
    return(rep(0,3))
  }
  total <- 0
  undirected <- 0
  directed <- 0
  pairs_counted <- list()
  for (i in 1:nrow(true_cpdag)){
    for (j in 1:ncol(true_cpdag)){

      if (true_cpdag[i,j]==1){
        if (lpc_cpdag[i,j]==0 & lpc_cpdag[j,i]==0){

          # Check to see if this pair has been considered
          if (checkPair(pairs_counted,i,j)){
            total <- total + 1
            if (true_cpdag[j,i]==1){
              undirected <- undirected + 1
            } else {
              directed <- directed + 1
            }
            pairs_counted[[total]] <- c(i,j)
          }
        }
      }
    }
  }
  return(c("undirected"=undirected,"directed"=directed,"total"=total))
}

# This function will calculate all the wrongly oriented edges
wrongDirection <- function(true_cpdag,lpc_cpdag) {
  if (is.null(true_cpdag) & is.null(lpc_cpdag)){
    return(0)
  }
  total <- 0
  for (i in 1:nrow(true_cpdag)){
    for (j in 1:ncol(true_cpdag)){

      if (true_cpdag[i,j] == 1 & true_cpdag[j,i] == 0) {
        if (lpc_cpdag[i,j] == 0 & lpc_cpdag[j,i] == 1){
          total <- total + 1
        }
      }

    }
  }

  return(total)
}

# This function will calculate all the edges that are directed when they should be undirected
directed_undirected <- function(true_cpdag,lpc_cpdag) {
  if (is.null(true_cpdag) & is.null(lpc_cpdag)){
    return(0)
  }
  total <- 0
  for (i in 1:nrow(true_cpdag)){
    for (j in 1:ncol(true_cpdag)){

      if (true_cpdag[i,j] == 1 & true_cpdag[j,i] == 1) {
        if (lpc_cpdag[i,j] == 1 & lpc_cpdag[j,i] == 0){
          total <- total + 1
        }
      }

    }
  }

  return(total)
}

# This function will calculate all the edges that are undirected (Estimate) when they should be directed (True)
undirected_directed <- function(true_cpdag,lpc_cpdag) {
  if (is.null(true_cpdag) & is.null(lpc_cpdag)){
    return(0)
  }
  total <- 0
  for (i in 1:nrow(true_cpdag)){
    for (j in 1:ncol(true_cpdag)){

      if (true_cpdag[i,j] == 1 & true_cpdag[j,i] == 0) { # Truth is directed
        if (lpc_cpdag[i,j] == 1 & lpc_cpdag[j,i] == 1){ # Estimated CPDAG is undirected
          total <- total + 1
        }
      }

    }
  }

  return(total)
}

# Returns the number of true positives in the graph
get_tp <- function(true_cpdag,lpc_cpdag){
  if (is.null(true_cpdag) & is.null(lpc_cpdag)){
    return(0)
  }
  total <- 0
  for (i in 1:nrow(true_cpdag)){
    for (j in 1:ncol(true_cpdag)){
      
      if (true_cpdag[i,j] == 1 & true_cpdag[j,i] == 0) { # Truth is directed
        if (lpc_cpdag[i,j] == 1 & lpc_cpdag[j,i] == 0){ # Estimated CPDAG matches
          total <- total + 1
        }
      }
      
      if (true_cpdag[i,j] == 1 & true_cpdag[j,i] == 1) { # Truth is undirected
        if (lpc_cpdag[i,j] == 1 & lpc_cpdag[j,i] == 1){ # Estimated CPDAG matches
          total <- total + 0.5
        }
      } 
      
    }
  }
  
  return(total)
}

get_neighborhood_edge_num <- function(truth){
  num_edges <- 0

  p_neighborhood <- 1:ncol(truth)
  for (i in p_neighborhood){
    remaining <- setdiff(p_neighborhood,i)
    for (j in remaining){
      if (truth[i,j]==1 & truth[j,i]==1){
        num_edges <- num_edges + 0.5
      } else if (truth[i,j]==1 & truth[j,i]==0){
        num_edges <- num_edges + 1
      }
    }
  }
  return(num_edges)
}


# Wrapper function to compute all the distance metrics
graph_dist <- function(truth,estimate){
  vs_dag <- vstructs(truth)
  vs_t_cpdag <- vstructs(estimate)

  same_v_structs <- nrow(vs_dag) == nrow(vs_t_cpdag)
  if (same_v_structs){
    same_v_structs <- all(vs_dag == vs_t_cpdag)
  } else {
    if (!acyclic(estimate)){
      same_v_structs <- 'CYCLES'
    }
  }

  truth_amat <- amat(truth)
  num_edges <- get_neighborhood_edge_num(truth_amat)
  estimate_amat <- amat(estimate)

  added <- addedEdges(truth_amat,estimate_amat)
  missing <- missingEdges(truth_amat,estimate_amat)
  
  my_tp <- get_tp(truth_amat,estimate_amat)

  metrics <- c(num_edges,missing["undirected"],missing["directed"],missing["total"],
               added["undirected"],added["directed"],added["total"],
               wrongDirection(truth_amat,estimate_amat),
               directed_undirected(truth_amat,estimate_amat),
               undirected_directed(truth_amat,estimate_amat),my_tp,
               bnlearn::shd(estimate,truth),unlist(bnlearn::compare(truth,estimate)),
               same_v_structs)

  return(metrics)
}

results <- simulate_local_pc(vars)


