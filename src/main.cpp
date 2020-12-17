#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "sharedFunctions.h"
#include "skeletonHelpersEfficient.h"
#include "vStructHelpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
List pc_sample_get_skeleton_efficient_cpp(List var_list,arma::mat df,double signif_level=0.05){
  int l = -1;
  int lmax = var_list["lmax"];
  int p = var_list["p"];
  bool verbose = var_list["verbose"];
  NumericMatrix C = var_list["C_tilde"];
  NumericMatrix true_dag = var_list["true_dag"];
  NumericVector neighborhood = var_list["neighborhood"];
  int N = neighborhood.size();
  StringVector names = var_list["names"];
  List S = var_list["S"];
  
  arma::mat R = arma::cor(df);
  int n = df.n_rows;
  
  NumericVector neighbors;
  NumericVector edges_i;
  NumericVector sep;
  
  NumericMatrix kvals;
  double pval=0.0;
  
  int num_tests=0;
  
  while (l < lmax){
    l += 1;
    if (verbose){
      Rcout << "The value of l is " << l << std::endl;
    }
    
    for (int i=0;i<N;++i){
      if (verbose){
        Rcout << "The value of i is " << i << std::endl;
      }
      // Work through potential neighbors with separating set of size l
      // These potential neighbors are those currently connected to node i in the current iteration's estimated graph
      edges_i = get_current_edges(i,N,C);
      for (NumericVector::iterator it = edges_i.begin(); it != edges_i.end(); ++it){
        int j = *it;
        if (j != i){
    
          if (verbose){
            Rcout << "The value of j is " << j << std::endl;
          }
          // Find neighbors of i and j from the true DAG (or they are estimated)
          // These neighbors are using the true node numbers (check documentation for this function)
          neighbors = get_potential_sep(i,j,neighborhood,N,true_dag);
          
          // If there are enough potential neighbors to match the current separating set size, we continue
          if (neighbors.length()>= l){
            if (verbose){
              Rcout << "There are " << neighbors.length() << " neighbor(s).\n";
            }
            kvals = combn_cpp(neighbors,l);
            
            check_separation_sample_efficient(l,i,j,kvals,sep,true_dag,names,neighborhood,C,S,pval,num_tests,R,n,signif_level,verbose);
            
            if (verbose){
              iteration_print(l,i,j,sep,names,pval);
            }
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
    _["neighborhood"]=neighborhood,
    _["p"]=p,
    _["NumTests"]=num_tests,
    _["verbose"]=verbose
  );
}

// [[Rcpp::export]]
List pc_sample_efficient_cpp(NumericMatrix true_dag,arma::mat df,
                             int target,
                             StringVector names,int lmax=3,
                             double signif_level = 0.05,
                             bool verbose=true,bool verbose_small=true){
  
  List var_list = pc_sample_skeleton_setup_efficient_cpp(true_dag,target,names,lmax,verbose); // Well tested
  
  List final_skeleton_list = pc_sample_get_skeleton_efficient_cpp(var_list,df,signif_level);
  return get_v_structures_efficient(final_skeleton_list);
}
