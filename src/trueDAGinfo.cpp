#include <Rcpp.h>
using namespace Rcpp;

/*
 * This function helps us to get the neighbors from the true DAG.
 * Input: target node i, number of nodes p, and the true DAG
 * Returns: Vector of neighbors 
 */
NumericVector get_neighbors_from_dag(int i,int p,NumericMatrix true_dag){
  NumericVector neighbors;
  NumericVector parents;
  NumericVector children;
  
  for (int j = 0;j<p;++j){
    if (true_dag(j,i)==1){
      parents.push_back(j);
      //Rcout << "Call from get_neighbors_from_dag. Node " << j << " is a parent.\n";
    } else if (true_dag(i,j)==1){
      children.push_back(j);
      //Rcout << "Call from get_neighbors_from_dag. Node " << j << " is a child.\n";
    }
  }
  
  NumericVector potential_spouses;
  int current_val;
  for (NumericVector::iterator it = children.begin(); it != children.end(); ++it){
    //Rcout << "Call from get_neighbors_from_dag. We are evaluating the following child: " << *it << std::endl;
    for (int j = 0; j<p; ++j){
      current_val = true_dag(j,*it);
      if (current_val == 1 & i != j){
        potential_spouses.push_back(j);
        //Rcout << "Call from get_neighbors_from_dag. Node " << j << " is a potential spouse.\n";
      }
    }
  }
  
  neighbors = union_(parents,children);
  neighbors = union_(neighbors,potential_spouses);
  
  return neighbors;
}