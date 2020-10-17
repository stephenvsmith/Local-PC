#include "printFunctions.h"

using namespace Rcpp;

NumericVector get_adjacent(NumericMatrix M,int i){
  NumericVector final;
  int p = M.ncol();
  
  for (int j=0;j<p;++j){
    if (M(i,j)!=0){
      final.push_back(j);
    }
  }
  
  return final;
}

NumericVector get_nonadjacent(NumericMatrix M,int i){
  NumericVector final;
  int p = M.ncol();
  
  for (int j=0;j<p;++j){
    if (M(i,j)==0 && j!=i){
      final.push_back(j);
    }
  }
  
  return final;
}

bool check_membership(NumericVector x,int i){
  
  NumericVector::iterator it = x.begin();
  int j;
  
  for (;it!=x.end();++it){
    j = *it;
    if (i == j){
      return true;
    }
  }
  
  return false;
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

