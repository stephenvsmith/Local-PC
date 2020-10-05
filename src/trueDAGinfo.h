#ifndef TRUEDAGINFO_H
#define TRUEDAGINFO_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector get_neighbors_from_dag(int i,int p,NumericMatrix true_dag);

#endif