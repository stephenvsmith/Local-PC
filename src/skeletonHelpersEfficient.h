#ifndef SKELETONHELPERSEFFICIENT_H
#define SKELETONHELPERSEFFICIENT_H

using namespace Rcpp;

List create_conditioning_sets_efficient_cpp(int p);

List pc_sample_skeleton_setup_efficient_cpp(NumericMatrix &true_dag,const int &target,StringVector &names,const int &lmax,bool &verbose);

List change_S_efficient(List S,int i,int j,NumericVector sep);

List change_S_0_efficient(List S,int i,int j);

void check_separation_sample_efficient(const int &l,const int &i,const int &j,
                             const NumericMatrix &kvals,
                             NumericVector &sep,NumericMatrix true_dag,
                             const StringVector &names,const NumericVector &neighborhood,
                             NumericMatrix C,
                             List S,double &pval,int &num_tests,arma::mat &df,int &n,
                             double &signif_level,bool &verbose);

NumericVector get_potential_sep(const int &i,const int &j,const NumericVector &neighborhood,const int &p,const NumericMatrix &true_dag);

#endif