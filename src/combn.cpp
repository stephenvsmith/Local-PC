#include <Rcpp.h>
using namespace Rcpp;


int factorial(int n) {
  if (n == 1){
    return n;
  } else {
    return n * factorial(n-1);
  }
}



int choose_cpp(int n,int k){
  int result = 1;
  for (int i = 1;i<=k;++i){
    result *= (n - (k - i)) / i;
  }
  
  return result;
}
/*
int choose_cpp(int n,int k){
  return factorial(n) / (factorial(k)*factorial(n-k));
}
*/
 
NumericVector successor_set(NumericVector x,NumericVector x1){
  int k = x1.length();
  int n = x.length();
  int m;
  NumericVector x2(x1);
  
  int j = 0;
  if (k % 2 != 0){
    ++j;
    m = n;
  }
  
  return x2;
}


void transfer_over(IntegerVector &main_vector,IntegerVector subset,IntegerVector subset_new){
  if (subset.length() != subset_new.length()){
    Rcout << "We have a problem in transfer_over.\n";
  }
  int l = subset.length();
  int index;
  int index_new;
  
  for (int i = 0;i<l;++i){
    index = subset[i];
    index_new = subset_new[i];
    Rcout << "We are exchanging index " << index << ". main_vector: " << main_vector[index];
    Rcout << " | subset_new: " << index_new << std::endl;
    main_vector[index] = index_new;
  }
}


void transfer_to_column(NumericMatrix &result,int j,NumericVector &x){
  if (result.nrow() != x.length()){
    Rcout << "Dimensions are incorrect.\n";
  } else {
    int k = x.length();
    for (int i = 0;i<k;++i){
      result(i,j) = x(i);
    }
    Rcout << result << std::endl; 
  }
}


NumericMatrix combn_cpp(NumericVector x,int k){
  int n = x.length();
  int total_combns = choose_cpp(n,k);
  //int total_combns = choose(n,k);
  //double a_i;
  NumericMatrix result(k,total_combns);
  
  int e = 0;
  int h = k;
  int i;
  IntegerVector j;
  int nmmp1;
  IntegerVector a = seq_len(k)-1;
  Rcout << a[0] << std::endl;
    
  NumericVector r = x[a];
  transfer_to_column(result,0,r);
    
  if (k > 0) {
    i = 1;
    nmmp1 = n - k; // This is the number of different values we can start with
    while (a[0] != nmmp1) {
      if (e < n - h - 1) {
        h = 1;
        e = a[k-1];
        j = NumericVector::create(1);
      } else {
        e = a[k - h - 1];
        h++;
        j = seq_len(h);
      }
      transfer_over(a,j+k-h-1,j+e);
      r = x[a];
      transfer_to_column(result,i,r);
      ++i;
    }
  }
  
  return result;
}



/*** R
factorial(4)

choose_cpp(6,4)
microbenchmark::microbenchmark(choose(6,4),choose_cpp(6,4))

combn(c(3,7,8,9),2)

combn_cpp(c(3,7,8,9),2)

transfer_over(1:6,3:4,2:3)

transfer_to_column(matrix(c(1:16),nrow = 4),1,c(19:22))

x <- rnorm(15)
#microbenchmark::microbenchmark(combn(x,5),combn_cpp(x,5))
*/
