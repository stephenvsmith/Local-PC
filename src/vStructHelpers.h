#ifndef VSTRUCTHELPERS_H
#define VSTRUCTHELPERS_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector get_adjacent(NumericMatrix M,int i);
NumericVector get_nonadjacent(NumericMatrix M,int i);
bool check_membership(NumericVector x,int i);
List get_v_structures(List L);
List get_v_structures_efficient(List L);

#endif