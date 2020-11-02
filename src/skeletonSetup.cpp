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

/*
 * This function sets up the nested lists that will hold separating sets
 * It will do this efficiently by only setting up lists for the neighbors of interest to the algorithm
 */
// [[Rcpp::export]]
void create_conditioning_sets_efficient_cpp(List &S,NumericVector &neighbors){
  int l = neighbors.length();
  String i_char;
  String j_char;
  for (int i=0;i<l;++i){
    List sublist = List(); // Creating the sublist for neighbor i
    for (int j=0;j<l;++j){
      if (j != i){
        j_char = String((char) neighbors(j));
        sublist[j_char] = NA_REAL;
      }
    }
    i_char = String((char) neighbors(i));
    S[i_char] = sublist;
  }
}

/*
 * This function sets up the nested lists that will hold separating sets
 * It will do this efficiently by only setting up lists for the neighbors of interest to the algorithm
 */
// [[Rcpp::export]]
List create_conditioning_sets_efficient_cpp2(NumericVector &neighbors){
  int l = neighbors.length();
  String i_char;
  String j_char;
  List S(0);
  for (int i=0;i<l;++i){
    List sublist = List(); // Creating the sublist for neighbor i
    for (int j=0;j<l;++j){
      if (j != i){
        j_char = String((char) neighbors(j));
        sublist[j_char] = NA_REAL;
      }
    }
    i_char = String((char) neighbors(i));
    S[i_char] = sublist;
  }
  return S;
}

/*
 * This function returns a complete graph for the neighbors of the target node
 */

// [[Rcpp::export]]
NumericMatrix get_initial_graph(int target,int p,NumericMatrix &true_dag){
  
  NumericMatrix C_tilde(p);
  
  // Find the neighborhood of the target node
  NumericVector neighbors = get_neighbors_from_dag(target,p,true_dag);
  neighbors.push_front(target);
  
  int node1;
  int node2;
  NumericVector::iterator it1;
  NumericVector::iterator it2;
  
  for (it1 = neighbors.begin(); it1 != neighbors.end()-1; ++it1){
    node1 = *it1;
    it2 = it1+1;
    while(it2 != neighbors.end()){
      node2 = *it2;
      C_tilde(node1,node2) = 1;
      C_tilde(node2,node1) = 1;
      ++it2;
    }
  }
  
  return C_tilde;
  
}

/*
 * The following function sets up the basic data structures for the skeleton algorithm
 */

List pc_pop_skeleton_setup_cpp(NumericMatrix &true_dag,StringVector &names,const int &lmax,bool &verbose){
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
 * The following function sets up the basic data structures for the skeleton algorithm that is used
 * for the sample version of the algorithm
 */
// [[Rcpp::export]]
List pc_sample_skeleton_setup_cpp(NumericMatrix &true_dag,const int &target,StringVector &names,const int &lmax,bool &verbose){
  // Number of nodes
  int p = 0;
  p = true_dag.nrow();
  
  if (verbose){
    Rcout << "There are " << p << " nodes in the DAG.\n";
  }
  
  // Initial graph that will be modified through the process of the algorithm
  NumericMatrix C_tilde = get_initial_graph(target,p,true_dag);
  
  if (verbose){
    Rcout << "Our starting matrix is " << C_tilde.nrow() << "x" << C_tilde.ncol() << ".\n";
  }
  
  // Create the list that will store 
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
 * The following function sets up the basic data structures for the skeleton algorithm that is used
 * for the sample version of the algorithm
 */
// [[Rcpp::export]]
List pc_sample_skeleton_setup_efficient_cpp(NumericMatrix &true_dag,
                                            const int &target,StringVector &names,
                                            const int &lmax,bool &verbose){
  // Number of nodes
  int p;
  p = true_dag.nrow();
  
  // Find the neighborhood of the target node
  NumericVector neighbors = get_neighbors_from_dag(target,p,true_dag);
  neighbors.push_front(target);
  std::sort(neighbors.begin(),neighbors.end());
  
  int N = neighbors.size();
  
  if (verbose){
    Rcout << "There are " << p << " nodes in the DAG.\n";
    Rcout << "There are " << N << " nodes in the neighborhood.\n";
  }
  
  // Initial graph that will be modified through the process of the algorithm
  NumericMatrix C_tilde(N);
  std::fill(C_tilde.begin(), C_tilde.end(), 1);
  C_tilde.fill_diag(0);
  
  if (verbose){
    Rcout << "Our starting matrix is " << C_tilde.nrow() << "x" << C_tilde.ncol() << ".\n";
  }
  
  // Create the list that will store 
  List S(0);
  create_conditioning_sets_efficient_cpp(S,neighbors);
  
  std::vector<double> p_vals;
  
  return List::create(
    _["p"] = p,
    _["C_tilde"]=C_tilde,
    _["true_dag"]=true_dag,
    _["names"]=names,
    _["neighborhood"]=neighbors,
    _["lmax"]=lmax,
    _["S"]=S,
    _["verbose"]=verbose,
    _["p_vals"]=p_vals);
  
}
