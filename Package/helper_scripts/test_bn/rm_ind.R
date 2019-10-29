# rm_ind

rm_ind <- function(result,true_dag,target,path_start){
  g <- result$G
  ind <- simplify_dag(g,target)

  cat("Number of Vertices:",ncol(true_dag),"\n",
      file = paste0(path_start,"Dropbox/Academics/Research/Build Notes/test_bn_notes.txt"),append = TRUE)
  cat("Target:",paste(target,sep = ","),"\n",
      file = paste0(path_start,"Dropbox/Academics/Research/Build Notes/test_bn_notes.txt"),append = TRUE)
  cat("Indices to remove:",paste(ind,sep = ", "),"\n",
      file = paste0(path_start,"Dropbox/Academics/Research/Build Notes/test_bn_notes.txt"),append = TRUE)

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

  cat("Building the simplified DAGs for comparison\n",
      file = paste0(path_start,"Dropbox/Academics/Research/Build Notes/test_bn_notes.txt"),append = TRUE)
  dag <- empty.graph(names1)
  amat(dag) <- g

  t_dag <- empty.graph(names1)
  amat(t_dag) <- true_dag

  return(list("dag"=dag,"t_dag"=t_dag))
}


