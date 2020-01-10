#### File for testing individual situations

network <- "hepar2"
target <- 51

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


res <- local_pc2(true_dag = true_dag,data = data,target = target,lmax=3,verbose = FALSE,verbose_small = FALSE,pop = FALSE)
