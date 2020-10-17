#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "printFunctions.h"
#include "trueDAGinfo.h"
#include "skeletonSetup.h"
#include "skeletonHelpers.h"
#include "vStructHelpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
List pc_sample_get_skeleton_cpp(List var_list,arma::mat df,double signif_level=0.95){
  int l = -1;
  int lmax = var_list["lmax"];
  int p = var_list["p"];
  bool verbose = var_list["verbose"];
  NumericMatrix C_tilde = var_list["C_tilde"];
  NumericMatrix C = clone(C_tilde);
  NumericMatrix true_dag = var_list["true_dag"];
  StringVector names = var_list["names"];
  List S = var_list["S"];
  
  arma::mat R = arma::cor(df);
  int n = df.n_rows;
  
  NumericVector exceptions; //= {0,0};
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
    
    for (int i=0;i<p;++i){
      if (verbose){
        Rcout << "The value of i is " << i << std::endl;
      }
      // Work through potential neighbors with separating set of size l
      // These potential neighbors are those currently connected to node i in the estimated graph
      edges_i = get_current_edges(i,p,C);
      for (NumericVector::iterator it = edges_i.begin(); it != edges_i.end(); ++it){
        int j = *it;
        if (j != i){
          
          if (verbose){
            Rcout << "The value of j is " << j << std::endl;
          }
          // Find neighbors of i and j from the true DAG (will have to be estimated at some point)
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
            
            check_separation_sample(l,i,j,kvals,sep,true_dag,names,C,S,pval,R,n,signif_level,verbose);
            
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
List pc_sample_cpp(NumericMatrix true_dag,arma::mat df,
                   int target,
                   StringVector names,int lmax=3,
                   double signif_level = 0.95,
                   bool verbose=true,bool verbose_small=true){
  
  List var_list = pc_sample_skeleton_setup_cpp(true_dag,target,names,lmax,verbose);
  
  List final_skeleton_list = pc_sample_get_skeleton_cpp(var_list,df,signif_level);
  return get_v_structures(final_skeleton_list);
}