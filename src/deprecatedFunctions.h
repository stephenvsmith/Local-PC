#ifndef DEPRECATEDFUNCTIONS_H
#define DEPRECATEDFUNCTIONS_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

List create_conditioning_sets_cpp(int p);

List change_S(List S,int i,int j,NumericVector sep);

List change_S_0(List S,int i,int j);

void check_separation_efficient(const int &l,const int &i,const int &j,
                                const NumericMatrix &kvals,Function get_pval,
                                NumericVector &sep,NumericMatrix true_dag,
                                const StringVector &names,NumericMatrix C,
                                List S,double &pval,bool &verbose);

#endif