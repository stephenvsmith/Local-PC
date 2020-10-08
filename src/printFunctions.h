#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

void print_vector_elements(NumericVector v,StringVector names, String opening="",String closing="");

void print_vector_elements_nonames(NumericVector v,String opening="",String closing="",String sep=" ");

void print_matrix(NumericMatrix m);

void print_S_vals(List S);

#endif