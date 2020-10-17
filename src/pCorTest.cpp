#include "pCorTest.h"

using namespace Rcpp;

// Here, we are treating C as the correlation matrix
// [[Rcpp::export]]
double get_partial_correlation(arma::mat C,int i,int j,arma::uvec k){
  double pc;
  int k_size = k.size();
  
  if (k_size==0){
    pc = C(i,j);
  } else if (k_size==1){
    pc = (C(i,j)-C(i,k(0))*C(j,k(0))) / sqrt( (1-pow(C(i,k(0)),2))*(1-pow(C(j,k(0)),2)));
  } else {
    arma::uvec indices(k_size+2);
    indices(0) = i;
    indices(1) = j;
    for (int l=0;l<k_size;++l){
      indices(l+2) = k(l);
    }
    arma::mat Cinv = arma::pinv(C(indices,indices));
    pc = -Cinv(0, 1)/sqrt(Cinv(0, 0) * Cinv(1, 1));
  }
  
  return pc;
}

double log_part(double r){
  return log1p((2*r)/(1-r));
}

// [[Rcpp::export]]
double fisherZ(double pc,int n,int k_size){
  return sqrt(n - k_size - 3) * 0.5 * log_part(pc);
}

// [[Rcpp::export]]
List condIndTest(arma::mat C,int i,int j,arma::uvec k,int n,double signif_level){
  double pc = get_partial_correlation(C,i,j,k);
  double pc_transformed = fisherZ(pc,n,k.size());
  bool lower = pc_transformed < 0;
  
  double cutoff = R::qnorm((1+signif_level)/2,0.0,1.0,true,false);
  //Rcpp::Rcout << "Value = " << pc_transformed << " | Cutoff = " << cutoff << std::endl;
  return List::create(
    _["result"]=abs(pc_transformed) <= cutoff,
    _["pval"]=R::pnorm(pc_transformed,0.0,1.0,lower,false)                           
  );
}
