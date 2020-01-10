lpc_builds <- function(true_dag,true_cpdag,num_sep_nodes,net) {
  
  ### Generating all possible target sets
  nodes <- 1:nrow(true_dag)
  target_sets <- combn(nodes,num_sep_nodes)

  # 30 maximum different target sets are tested
  num <- min(30,ncol(target_sets))
  target_sets <- target_sets[,sample(1:ncol(target_sets),num,replace = FALSE)]
  if (!is.matrix(target_sets)){
    target_sets <- matrix(target_sets,ncol = length(target_sets))
  }
  tmp <- apply(target_sets,2,function(target,net){
    # Create the necessary files
    if (!file.exists(paste0(net,"; target=",paste(target,collapse = ","))))
      dir.create(paste0(net,"; target=",paste(target,collapse = ",")))
    setwd(paste0("./",net,"; target=",paste(target,collapse = ",")))

    temp <- test_bn(net,target,true_dag,path_start=path_start)
    res <- local_pc_dist(target,true_dag,true_cpdag)

    results <<- rbind(results,c(net,paste(target,collapse = ","),res))
    current_network <<- rbind(current_network,c(net,paste(target,collapse = ","),res))

    write.table(results[nrow(results),],"measurements.txt")
    setwd("..")
  },net)

  return()
}
