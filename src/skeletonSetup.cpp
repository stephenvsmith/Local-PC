#include <Rcpp.h>
#include "trueDAGinfo.h"
using namespace Rcpp;

/*
 * This function sets up the nested lists that will hold separating sets
 */
List create_conditioning_sets_cpp(int p){
  List S(p);
  for (int i=0;i<p;++i){
    List sublist = List(p);
    for (int j=0;j<p;++j){
      sublist[j] = NA_REAL;
    }
    S[i] = sublist;
  }
  return S;
}

// [[Rcpp::export]]
NumericMatrix get_initial_graph(int target,int p,NumericMatrix true_dag){
  
  NumericMatrix C_tilde(p);
  
  NumericVector neighbors = get_neighbors_from_dag(target,p,true_dag);
  
  int node1;
  int node2;
  NumericVector::iterator it1;
  NumericVector::iterator it2;
  
  for (it1 = neighbors.begin(); it1 != neighbors.end(); ++it1){
    it2 = it1+1;
    node1 = *it1;
    node2 = *it2;
    C_tilde(node1,node2) = 1;
    C_tilde(node2,node1) = 1;
  }
  
  return C_tilde;
  
}

/*
 * The following function sets up the basic data structures for the skeleton algorithm
 */

List pc_pop_skeleton_setup_cpp(NumericMatrix true_dag,StringVector names,int lmax,bool verbose){
  // Number of nodes
  int p = 0;
  p = true_dag.nrow();
  
  if (verbose){
    Rcout << "There are " << p << " nodes in the DAG.\n";
  }
  
  NumericMatrix C_tilde(p);
  std::fill(C_tilde.begin(), C_tilde.end(), 1);
  C_tilde.fill_diag(0);

  if (verbose){
    Rcout << "Our starting matrix is " << C_tilde.nrow() << "x" << C_tilde.ncol() << ".\n";
  }
  
  List S = create_conditioning_sets_cpp(p);
  
  std::vector<double> p_vals;
  
  return List::create(
    _["p"] = p,
    _["C_tilde"]=C_tilde,
    _["true_dag"]=true_dag,
    _["names"]=names,
    _["lmax"]=lmax,
    _["S"]=S,
    _["verbose"]=verbose,
    _["p_vals"]=p_vals);
  
}

/*
 * The following function sets up the basic data structures for the skeleton algorithm
 */

List pc_sample_skeleton_setup_cpp(NumericMatrix true_dag,int target,StringVector names,int lmax,bool verbose){
  // Number of nodes
  int p = 0;
  p = true_dag.nrow();
  
  if (verbose){
    Rcout << "There are " << p << " nodes in the DAG.\n";
  }
  
  NumericMatrix C_tilde = get_initial_graph(target,p,true_dag);
  /*
   NumericMatrix C_tilde(p,p);
   std::fill(C_tilde.begin(), C_tilde.end(), 1);
   C_tilde.fill_diag(0);
   */
  if (verbose){
    Rcout << "Our starting matrix is " << C_tilde.nrow() << "x" << C_tilde.ncol() << ".\n";
  }
  
  List S = create_conditioning_sets_cpp(p);
  
  std::vector<double> p_vals;
  
  return List::create(
    _["p"] = p,
    _["C_tilde"]=C_tilde,
    _["true_dag"]=true_dag,
    _["names"]=names,
    _["lmax"]=lmax,
    _["S"]=S,
    _["verbose"]=verbose,
    _["p_vals"]=p_vals);
  
}
