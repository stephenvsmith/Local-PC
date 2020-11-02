#ifndef SKELETONSETUP_H
#define SKELETONSETUP_H

#include <Rcpp.h>
using namespace Rcpp;

List create_conditioning_sets_cpp(int p);

List pc_pop_skeleton_setup_cpp(NumericMatrix &true_dag,StringVector &names,const int &lmax,bool &verbose);

List pc_sample_skeleton_setup_cpp(NumericMatrix &true_dag,const int &target,StringVector &names,const int &lmax,bool &verbose);

List pc_sample_skeleton_setup_efficient_cpp(NumericMatrix &true_dag,const int &target,StringVector &names,const int &lmax,bool &verbose);

#endif