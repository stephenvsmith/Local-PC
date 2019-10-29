### Function for network results

# Change Working Directory to Store Results
setwd(paste0(path_start,"Dropbox/Academics/Research/Results/Population"))
# Get the true DAG
data.grid$network <- net

# Set up File Storage for Data and Results
if (!file.exists(paste0(path_start,"Dropbox/Academics/Research/Data/",net)))
  dir.create(paste0(path_start,"Dropbox/Academics/Research/Data/",net))

# Generate Data and Get True DAG
gdg <- generate.data.grid(data.grid,
                          out.dir=paste0(path_start,"Dropbox/Academics/Research/Data/",net),
                          verbose=FALSE,path.start=ps2)
dir.name <- paste0(net,"; ","n = 1000; c = 0")
true_dag <- as.matrix(read.table(paste0(path_start,
                                        "Dropbox/Academics/Research/Data/",
                                        net,"/",dir.name,"/trueDAG.txt")))

true_cpdag <- as.matrix(read.table(paste0(path_start,
                                        "Dropbox/Academics/Research/Data/",
                                        net,"/",dir.name,"/trueCPDAG.txt")))
true_cpdag <- matrix(as.numeric(true_cpdag),nrow = nrow(true_cpdag))
