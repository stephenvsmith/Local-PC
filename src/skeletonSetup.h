#ifndef SKELETONSETUP_H
#define SKELETONSETUP_H

#include <Rcpp.h>
using namespace Rcpp;

List create_conditioning_sets_cpp(int p);

List pc_pop_skeleton_setup_cpp(NumericMatrix true_dag,StringVector names,int lmax,bool verbose);

#endif