# Running local pc for visualization
test_bn <- function(net,target,true_dag=NULL,dag_data=NULL,verbose=FALSE,path_start="~/") {
  cat("Target: ",target,"\n")
  # Setup
  source(paste0(path_start,"Desktop/Research/package/helper_scripts/test_bn/test_bn_setup.R"),local = TRUE)

  # Simplifying DAGs
  source(paste0(path_start,"Desktop/Research/package/helper_scripts/test_bn/rm_ind.R"),local = TRUE)
  imp_vars <- rm_ind(result,true_dag,target,path_start)

  if (verbose) cat("Target Node(s): ",paste(names[target],collapse = " | "),"\n")
  if (verbose) bnlearn::graphviz.plot(imp_vars$dag)
  png(filename = "comparison.png",width = 1200,height = 1200)
  par(mfrow=c(1,2))
  bnlearn::graphviz.plot(imp_vars$dag,main = "Local PC",highlight = list(nodes=names[target],col="green",fill="blue",textCol="white"))
  bnlearn::graphviz.plot(imp_vars$t_dag,main = "Truth",highlight = list(nodes=names[target],col="green",fill="blue",textCol="white"))
  dev.off()
  png(filename = "OriginalDAG.png",width = 1200,height = 1200)
  bnlearn::graphviz.plot(t_dag_original,main = "True DAG",highlight = list(nodes=names[target],col="green",fill="blue",textCol="white"))
  dev.off()
  return(result$p_vals)
}
