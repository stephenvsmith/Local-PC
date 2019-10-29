### Local PC Test Script

# Loading Data Generation Function
source("/Users/stephensmith/Desktop/Research/projects/code/data_gen_R.R")

# Obtain the names of all the networks
file_vec <- list.files("~/Desktop/Research/projects/bn_data_generation/networks/rds/")
nets <- sapply(file_vec, function(x) sub(".rds","",x))
nets

# Write functions to be used on each dataset

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

# Running local pc for visualization
test_bn <- function(net,target,true_dag,verbose=FALSE) {
  file <- paste0("~/Desktop/Research/projects/bn_data_generation/networks/rds/",net,".rds")
  names <- names(readRDS(file))
  result <- local_pc2(true_dag = true_dag,target = target,verbose = FALSE,verbose_small = FALSE)
  g <- result$G
  ind <- simplify_dag(g)
  if (!is.null(ind)){
    if (length(ind) == nrow(g)) return()
    g <- g[-ind,-ind]
    names1 <- names[-ind]
    rownames(g) <- names1
    colnames(g) <- names1
    true_dag <- true_dag[-ind,-ind]
    rownames(true_dag) <- names1
    colnames(true_dag) <- names1
  } else {
    names1 <- names
  }
  dag <- empty.graph(names1)
  amat(dag) <- g
  t_dag <- empty.graph(names1)
  amat(t_dag) <- true_dag
  if (verbose) cat("Target Node(s): ",paste(names[target],collapse = " | "),"\n")
  if (verbose) bnlearn::graphviz.plot(dag)
  jpeg(file = "comparison.jpeg")
  par(mfrow=c(1,2))
  bnlearn::graphviz.plot(dag,main = "Local PC")
  bnlearn::graphviz.plot(t_dag,main = "Truth")
  dev.off()
  return(g)
}

# Testing for distance
local_pc_dist <- function(target,true_dag){
  # Run Local PC for the target
  result <- local_pc2(true_dag = true_dag,target = target,verbose = FALSE,verbose_small = FALSE)
  g <- result$G

  # Create graph objects for the Local PC DAG and for the true DAG
  names <- as.character(1:ncol(g))
  colnames(g) <- as.character(names)
  rownames(g) <- as.character(names)
  colnames(true_dag) <- as.character(names)
  rownames(true_dag) <- as.character(names)

  # Simplify the DAGs to only contain neighborhoods of the targets
  ind <- simplify_dag(g)
  if (!is.null(ind)){
    if (length(ind) == nrow(g)) return(c(0,0))
    g <- g[-ind,-ind]
    true_dag <- true_dag[-ind,-ind]
  }

  lpc_dag <- empty.graph(names[-ind])
  amat(lpc_dag) <- g

  t_dag <- empty.graph(names[-ind])
  amat(t_dag) <- true_dag

  return(c(bnlearn::hamming(lpc_dag,t_dag),bnlearn::shd(lpc_dag,t_dag)))
}

local_pc_dist(target = 17,true_dag = td[["andes"]])
test_bn("asia",c(1,4),td[["asia"]],verbose = TRUE)
test_bn("andes",17,td[["andes"]])

local_pc2(true_dag = td[["asia"]],target = c(1,4),verbose = TRUE)
res <- local_pc2(true_dag = td[["andes"]],target = 17,verbose = TRUE)

results <- matrix(nrow=0,ncol = 4)
colnames(results) <- c("Network","Target","HD","SHD")
setwd("~/Dropbox/Academics/Research/Results")

data.grid = data.frame(network = "asia",
                       data.type = "continuous",
                       n.dat = 1,
                       n.obs = 1000,
                       c.ratio = 0,
                       max.in.degree = Inf,
                       lb = 0.5,  # lower bound of coefficients
                       ub = 1,  # upper bound of coefficients
                       low = 0,  # lower bound of variances if continuous, of number of levels if discrete
                       high = 1,  # upper bound of variances if continuous, of number of levels if discrete
                       scale = TRUE,
                       stringsAsFactors = FALSE)

for (net in nets){
  setwd("~/Dropbox/Academics/Research/Results")
  # Get the true DAG
  data.grid$network = net
  if (!file.exists(paste0("~/Dropbox/Academics/Research/Data/",net)))
    dir.create(paste0("~/Dropbox/Academics/Research/Data/",net))
  gdg <- generate.data.grid(data.grid,out.dir=paste0("/Users/stephensmith/Dropbox/Academics/Research/Data/",net),verbose=FALSE)
  dir.name <- paste0(net,"; ","n = 1000; c = 0")
  true_dag <- as.matrix(read.table(paste0("~/Dropbox/Academics/Research/Data/",net,"/",dir.name,"/trueDAG.txt")))

  if (!file.exists(paste0("./",net))) dir.create(paste0("./",net))
  setwd(paste0("./",net))

  # Get all results for one target
  if (!file.exists("./one_target")) dir.create("./one_target")
  setwd("./one_target")
  for (i in 1:nrow(true_dag)) {
    target <- i
    if (!file.exists(paste0(net,"; target=",target))) dir.create(paste0(net,"; target=",target))
    setwd(paste0("./",net,"; target=",target))
    temp <- test_bn(net,target,true_dag)
    res <- local_pc_dist(target,true_dag)
    results <- rbind(results,c(net,target,res))
    write.table(results[nrow(results),],"measurements.txt")

    setwd("..")
  }
}
