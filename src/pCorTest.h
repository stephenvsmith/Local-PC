#ifndef PCORTEST_H
#define PCORTEST_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double get_partial_correlation(arma::mat C,int i,int j,arma::uvec k);

double log_part(double r);

double fisherZ(double pc,int n,int k_size);

List condIndTest(arma::mat &C,const int &i,const int &j,const arma::uvec &k,const int &n,const double &signif_level=0.05);


#endif