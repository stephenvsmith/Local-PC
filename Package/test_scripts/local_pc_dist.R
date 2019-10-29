# Testing for distance

### NEED TO CHANGE THE DISTANCE METRICS TO BE MORE APPROPRIATE
source("~/Desktop/Research/package/helper_scripts/distance/missingEdges.R")
source("~/Desktop/Research/package/helper_scripts/distance/addedEdges.R")

local_pc_dist <- function(target,true_dag,true_cpdag,data=NULL){
  # Run Local PC for the target
  if (is.null(data)){
    result <- local_pc2(true_dag = true_dag,target = target,verbose = FALSE,verbose_small = FALSE,pop=TRUE)
  } else {
    result <- local_pc2(data = data,target = target,true_dag = true_dag,verbose = FALSE,verbose_small = FALSE,pop = FALSE)
  }
  g <- result$G

  # Create graph objects for the Local PC DAG and for the true DAG
  names <- as.character(1:ncol(g))
  colnames(g) <- as.character(names)
  rownames(g) <- as.character(names)
  colnames(true_dag) <- as.character(names)
  rownames(true_dag) <- as.character(names)
  colnames(true_cpdag) <- as.character(names)
  rownames(true_cpdag) <- as.character(names)

  # Simplify the DAGs to only contain neighborhoods of the targets
  ind <- simplify_dag(g,target)
  if (!is.null(ind)){
    if (length(ind) == nrow(g)) return(c(0,0))
    g <- g[-ind,-ind]
    true_dag <- true_dag[-ind,-ind]
    true_cpdag <- true_cpdag[-ind,-ind]
  }
  lpc_dag <- empty.graph(names[-ind])
  amat(lpc_dag) <- g

  t_cpdag <- empty.graph(names[-ind])
  amat(t_cpdag) <- true_cpdag

  vs_dag <- vstructs(lpc_dag)
  vs_t_cpdag <- vstructs(t_cpdag)
  same_v_structs <- nrow(vs_dag) == nrow(vs_t_cpdag)
  if (same_v_structs){
    same_v_structs <- all(vs_dag == vs_t_cpdag)
  }

  return(c(missingEdges(true_cpdag,g),missingEdges(true_cpdag,g),same_v_structs))
}
