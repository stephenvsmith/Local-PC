#include "deprecatedFunctions.h"
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
 * This function helps us to add a separating set for nodes i and j
 */
List change_S(List S,int i,int j,NumericVector sep){
  NumericVector sep_new;
  sep_new = clone(sep);
  // Rcout << "S before:\n";
  // print_S_vals(S);
  List sublist;
  sublist = S[i];
  sublist[j] = sep_new;
  S[i] = sublist;
  // Rcout << "S after:\n";
  // print_S_vals(S);
  
  return S;
}

/*
 * Allows us to change S to indicate that nodes i and j are separated without a separating set
 * Note: this is the inefficient version
 */
List change_S_0(List S,int i,int j){
  
  // Rcout << "S before:\n";
  // print_S_vals(S);
  List sublist0;
  sublist0 = S[i];
  sublist0[j] = -1;
  S[i] = sublist0;
  // Rcout << "S after:\n";
  // print_S_vals(S);
  
  return S;
}

/*
 * This function checks for separation between node i and node j
 * given any set from the matrix kvals (each column is a potential
 * separating set)

void check_separation_efficient(const int &l,const int &i,const int &j,
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
      change_S_0_efficient(S,i,j);
      change_S_0_efficient(S,j,i);
      
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
        change_S_efficient(S,i,j,sep);
        change_S_efficient(S,j,i,sep);
        C(i,j) = 0;
        C(j,i) = 0;
        keep_checking_k = false;
      }
      ++k;
    }
  }
}
*/
