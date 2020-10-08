#ifndef SKELETONHELPERS_H
#define SKELETONHELPERS_H

using namespace Rcpp;

NumericVector get_current_edges(int i,int p,NumericMatrix graph);

NumericMatrix combn_cpp(NumericVector x,int l);

List change_S(List S,int i,int j,NumericVector sep);

List change_S_0(List S,int i,int j);

void check_separation(const int &l,const int &i,const int &j,
                      const NumericMatrix &kvals,Function get_pval,
                      NumericVector &sep,NumericMatrix true_dag,
                      const StringVector &names,NumericMatrix C,
                      List S,double &pval,bool &verbose);

void check_separation_sample(const int &l,const int &i,const int &j,
                             const NumericMatrix &kvals,
                             NumericVector &sep,NumericMatrix true_dag,
                             const StringVector &names,NumericMatrix C,
                             List S,double &pval,arma::mat &df,int &n,
                             double &signif_level,bool &verbose);

#endif