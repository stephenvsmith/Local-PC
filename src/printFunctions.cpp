#include <Rcpp.h>

using namespace Rcpp;

/*
 * 
 * Functions for printing different values 
 * 
 */

void print_vector_elements(NumericVector v,StringVector names, String opening,String closing){
  int l = v.length();
  Rcout << opening.get_cstring();
  for (int i = 0;i<l;++i) {
    Rcout << names(v(i)) << " ";
  }
  Rcout << closing.get_cstring();
}


void print_vector_elements_nonames(NumericVector v,String opening,String closing,String sep){
  int l = v.length();
  Rcout << opening.get_cstring();
  for (int i = 0;i<l;++i) {
    Rcout << v(i) << sep.get_cstring();
  }
  Rcout << closing.get_cstring();
}


void print_matrix(NumericMatrix m){
  int n = m.nrow();
  int p = m.ncol();
  String ending;
  for (int i = 0;i<n;++i){
    for (int j = 0;j<p;++j){
      if (j==p-1){
        ending = "\n";
      } else{
        ending = " ";
      }
      Rcout << m(i,j) << ending.get_cstring();
    }
  }
}

void print_S_vals(List S){
  List sublist;
  for (int i=0;i<S.length();++i){
    sublist = S[i];
    for (int j=0;j<sublist.length();++j){
      Rcout << "S[[" << i << "]][[" << j << "]] = ";
      print_vector_elements_nonames(sublist[j],"",""," ");
      Rcout << " ";
    }
    Rcout << std::endl;
  }
}
