#ifndef SKELETONHELPERSEFFICIENT_H
#define SKELETONHELPERSEFFICIENT_H

using namespace Rcpp;

NumericVector get_current_edges_efficient(int i,int p,NumericMatrix graph);

List change_S_efficient(List S,int i,int j,NumericVector sep);
List change_S_0_efficient(List S,int i,int j);

void check_separation_efficient(const int &l,const int &i,const int &j,
                      const NumericMatrix &kvals,Function get_pval,
                      NumericVector &sep,NumericMatrix true_dag,
                      const StringVector &names,NumericMatrix C,
                      List S,double &pval,bool &verbose);

void check_separation_sample_efficient(const int &l,const int &i,const int &j,
                             const NumericMatrix &kvals,
                             NumericVector &sep,NumericMatrix true_dag,
                             const StringVector &names,const NumericVector &neighborhood,
                             NumericMatrix C,
                             List S,double &pval,arma::mat &df,int &n,
                             double &signif_level,bool &verbose);

NumericVector get_potential_sep(const int &i,const int &j,const NumericVector &neighborhood,const int &N,const NumericMatrix &true_dag);

#endif