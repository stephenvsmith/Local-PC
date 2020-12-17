#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "pCorTest.h"
#include "vStructHelpers.h"
#include "sharedFunctions.h"
#include "deprecatedFunctions.h"
using namespace Rcpp;

/*
 * This function returns a complete graph for the neighbors of the target node
 */
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
 * 
 * This function checks whether or not nodes i and j are separated by any of the 
 * sets in the matrix kvals for a given significance level
 * 
 */
void check_separation_sample(const int &l,const int &i,const int &j,
                             const NumericMatrix &kvals,
                             NumericVector &sep,NumericMatrix true_dag,
                             const StringVector &names,NumericMatrix C,
                             List S,double &pval,arma::mat &R,int &n,
                             double &signif_level,bool &verbose){
  int k;
  int kp = kvals.cols();
  bool keep_checking_k;
  List test_result;
  
  arma::uvec sep_arma;
  
  if (l == 0){
    sep = NA_REAL;
    //Rcout << "Size of arma vector when l=0: " << sep_arma.size() << std::endl;
    test_result = condIndTest(R,i,j,sep_arma,n,signif_level);
    pval = test_result["pval"];
    //pval = as<double>(get_pval(i,j,true_dag,names));
    //Rcout << "The p-value is " << pval << std::endl;
    if (test_result["result"]){
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
      sep_arma = as<arma::uvec>(sep);
      test_result = condIndTest(R,i,j,sep_arma,n,signif_level);
      pval = test_result["pval"];
      //pval = as<double>(get_pval(i,j,true_dag,names,sep));
      if (verbose){
        Rcout << "The p-value is " << pval << std::endl;
      }
      if (test_result["result"]){
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

List get_v_structures(List L) {
  
  List S = L["S"];
  NumericMatrix G = L["C"];
  bool verbose = L["verbose"];
  
  int p = G.ncol();
  int j;
  int k;
  
  bool no_neighbors;
  bool j_invalid;
  NumericVector placeholder;
  
  NumericVector i_adj;
  NumericVector j_adj;
  NumericVector j_vals;
  NumericVector k_vals;
  
  List sublist;
  if (verbose){
    Rcout << "Beginning loops to find v-structures.\n";
  }
  for (int i=0;i<p;++i){
    placeholder = G(i,_);
    no_neighbors = (all(placeholder==0)).is_true();
    if (!no_neighbors){
      if (verbose){
        Rcout << "i: " << i << std::endl;
      }
      i_adj = get_adjacent(G,i);
      j_vals = get_nonadjacent(G,i);
      
      for (NumericVector::iterator it=j_vals.begin();it != j_vals.end();++it){
        j = *it;
        // Node j has no children, j is parent to i, or we are repeating an analysis and this j should not be considered
        placeholder = G(j,_);
        j_invalid = (all(placeholder==0)).is_true();
        j_invalid = j_invalid || G(j,i)==1 || j > i;
        if (!j_invalid){
          if (verbose){
            Rcout << "j: " << j << std::endl;
          }
          j_adj = get_adjacent(G,j);
          k_vals = intersect(i_adj,j_adj);
          // If there are no common neighbors, move to next j
          if (k_vals.length()!=0){
            // We loop through all of the common neighbors
            for (NumericVector::iterator it2 = k_vals.begin();it2!=k_vals.end();++it2){
              k = *it2;
              if (verbose){
                Rcout << "k: " << k << std::endl; 
              }
              sublist = S[i];
              if (!check_membership(sublist[j],k)){
                if (verbose){
                  Rcout << "Separation Set: ";
                  print_vector_elements_nonames(sublist[j]);
                  Rcout << " | V-Structure: " << i << "->" << k << "<-" << j << std::endl; 
                }
                G(k,i) = 0;
                G(k,j) = 0;
              }
            }
          }
        }
      }
    }
  }
  
  
  return List::create(
    _["G"]=G,
    _["NumTests"]=L["NumTests"]
  );
}

/*
 * The following function sets up the basic data structures for the skeleton algorithm that is used
 * for the sample version of the algorithm
 */
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

List pc_sample_get_skeleton_cpp(List var_list,arma::mat df,double signif_level=0.05){
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
  
  NumericVector exceptions; // used to prevent us from considering conditional independence of a node with itself
  NumericVector neighbors;
  NumericVector edges_i;
  NumericVector sep;
  
  NumericVector::iterator it;
  
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
      for (it = edges_i.begin(); it != edges_i.end(); ++it){
        int j = *it;
        if (j != i){
          
          if (verbose){
            Rcout << "The value of j is " << j << std::endl;
          }
          // Find neighbors of i and j from the true DAG (or they are estimated)
          neighbors = union_(get_neighbors_from_dag(i,p,true_dag),get_neighbors_from_dag(j,p,true_dag));
          exceptions = NumericVector::create(i,j);
          neighbors = setdiff(neighbors,exceptions);
          
          // If there are enough potential neighbors to match the current separating set size, we continue
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
                   double signif_level = 0.05,
                   bool verbose=true,bool verbose_small=true){
  
  List var_list = pc_sample_skeleton_setup_cpp(true_dag,target,names,lmax,verbose);
  
  List final_skeleton_list = pc_sample_get_skeleton_cpp(var_list,df,signif_level);
  return get_v_structures(final_skeleton_list);
}

