#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "sharedFunctions.h"
#include "deprecatedFunctions.h"
#include "vStructHelpers.h"
using namespace Rcpp;

/*
 * The following function sets up the basic data structures for the skeleton algorithm that is used
 * for the sample version of the algorithm
 */
// [[Rcpp::export]]
List pc_pop_skeleton_setup_cpp(NumericMatrix &true_dag,StringVector &names,const int &lmax,bool &verbose){
  // Number of nodes
  int p = 0;
  p = true_dag.nrow();
  
  if (verbose){
    Rcout << "There are " << p << " nodes in the DAG.\n";
  }
  
  // Initial graph that will be modified through the process of the algorithm
  NumericMatrix C_tilde(p);
  std::fill(C_tilde.begin(), C_tilde.end(), 1);
  C_tilde.fill_diag(0);
  
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
 * This function checks for separation between node i and node j
 * given any set from the matrix kvals (each column is a potential
 * separating set)
 */
void check_separation(const int &l,const int &i,const int &j,
                      const NumericMatrix &kvals,Function get_pval,
                      NumericVector &sep,NumericMatrix true_dag,
                      const StringVector &names,NumericMatrix C,
                      List S,double &pval,bool &verbose){
  int k;
  int kp = kvals.cols();
  bool keep_checking_k; // Tracks to see whether or not to keep check for separating sets
  
  if (l == 0){
    sep = NA_REAL;
    pval = as<double>(get_pval(i,j,true_dag,names));
    //Rcout << "The p-value is " << pval << std::endl;
    if (pval == 1){
      change_S_0(S,i,j);
      change_S_0(S,j,i);
      
      C(i,j) = 0;
      C(j,i) = 0;
    }
  } else {
    k = 0;
    keep_checking_k = true;
    while (keep_checking_k & (k<kp)){
      sep = kvals( _ , k );
      pval = as<double>(get_pval(i,j,true_dag,names,sep));
      if (verbose){
        Rcout << "The p-value is " << pval << std::endl;
      }
      if (pval==1){
        if (verbose){
          Rcout << names(i) << " is separated from " << names(j) << " by node(s):\n";
          print_vector_elements(sep,names);
        }
        change_S(S,i,j,sep);
        change_S(S,j,i,sep);
        C(i,j) = 0;
        C(j,i) = 0;
        keep_checking_k = false;
      }
      ++k;
    }
  }
}

// [[Rcpp::export]]
List pc_pop_get_skeleton_cpp(List var_list){
  int l = -1;
  int lmax = var_list["lmax"]; // Maximum separating set size
  int p = var_list["p"]; // Number of nodes in the network
  bool verbose = var_list["verbose"]; // Whether or not to print diagnostic information
  NumericMatrix C_tilde = var_list["C_tilde"]; // The initial graph
  NumericMatrix C = clone(C_tilde); // Clone of the initial graph that will be modified
  NumericMatrix true_dag = var_list["true_dag"]; // A pxp matrix of the true DAG
  StringVector names = var_list["names"]; // A vector of length p with the names of the nodes
  List S = var_list["S"]; // A p^2 size list containing the separating sets for the different nodes.

  NumericVector exceptions; //= {0,0};
  NumericVector neighbors; 
  NumericVector edges_i;
  NumericVector sep;

  NumericMatrix kvals;
  double pval=0.0;

  Environment myEnv = Environment::namespace_env("LocalPC");
  Function get_pval= myEnv["get_pval"];

  int num_tests=0;

  while (l < lmax){
    l += 1;
    if (verbose){
      Rcout << "The value of l is " << l << std::endl;
    }

    for (int i=0;i<p;++i){
      if (verbose){
        Rcout << "The value of i is " << i << std::endl;
      }
      // Work through potential neighbors with separating set of size l
      edges_i = get_current_edges(i,p,C);
      for (NumericVector::iterator it = edges_i.begin(); it != edges_i.end(); ++it){
        int j = *it;
        if (j != i){

          if (verbose){
            Rcout << "The value of j is " << j << std::endl;
          }
          // Find neighbors of i and j
          neighbors = union_(get_neighbors_from_dag(i,p,true_dag),get_neighbors_from_dag(j,p,true_dag));
          exceptions = NumericVector::create(i,j);
          //Rcout << "Exceptions: " << exceptions(0) << " and " << exceptions(1) << std::endl;

          //print_vector_elements(neighbors,names,"Neighbors before:","=============");

          neighbors = setdiff(neighbors,exceptions);

          //print_vector_elements(neighbors,names,"Neighbors after:","=============");

          if (neighbors.length()>= l){

            if (verbose){
              Rcout << "There are " << neighbors.length() << " neighbor(s).\n";
            }
            kvals = combn_cpp(neighbors,l);

            check_separation(l,i,j,kvals,get_pval,sep,true_dag,names,C,S,pval,verbose);

            if (verbose){
              Rcout << "l: " << l << " | i: " << i << " | j: " << j << " | k: ";
              if (l == 0){
                Rcout << sep;
              } else {
                print_vector_elements(sep,names);
              }
              Rcout << " | p-val: " << pval;
              Rcout << std::endl;
              //print_S_vals(S);
            }
            ++num_tests;
          }
        }
      }

    }

  }
  if (verbose){
    Rcout << "The final C matrix:\n";
    print_matrix(C);
    Rcout << "Conclusion of algorithm.\n";
  }

  return List::create(
    _["C"]=C,
    _["S"]=S,
    _["NumTests"]=num_tests,
    _["verbose"]=verbose
  );
}



// [[Rcpp::export]]
List pc_pop_cpp(NumericMatrix true_dag,StringVector names,int lmax=3,bool verbose=true,bool verbose_small=true){
  
  List var_list = pc_pop_skeleton_setup_cpp(true_dag,names,lmax,verbose);
  
  List final_skeleton_list = pc_pop_get_skeleton_cpp(var_list);
  return get_v_structures(final_skeleton_list);
}

