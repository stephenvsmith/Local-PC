#ifndef SHAREDFUNCTIONS_H
#define SHAREDFUNCTIONS_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

NumericVector get_current_edges(int i,int p,NumericMatrix graph);

NumericMatrix combn_cpp(NumericVector x,int l);

NumericVector get_neighbors_from_dag(int i,int p,NumericMatrix true_dag);

void print_vector_elements(NumericVector v,StringVector names, String opening="",String closing="");

void print_vector_elements_nonames(NumericVector v,String opening="",String closing="",String sep=" ");

void print_matrix(NumericMatrix m);

void print_S_vals(List S);

void iteration_print(const int &l,const int &i,const int &j,const NumericVector &sep,const StringVector &names,const double &pval);

#endif