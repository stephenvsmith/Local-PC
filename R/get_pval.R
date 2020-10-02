get_pval <- function(i,j,true_dag,names,k=c()){
  e <- bnlearn::empty.graph(names)
  bnlearn::amat(e) <- true_dag
  if (length(k)==0){
    pval <- as.numeric(bnlearn::dsep(e,names[i+1],names[j+1]))
  } else {
    pval <- as.numeric(bnlearn::dsep(e,names[i+1],names[j+1],names[k+1]))
  }
  
  return(pval)
}