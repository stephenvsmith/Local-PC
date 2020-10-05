#include <Rcpp.h>

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

