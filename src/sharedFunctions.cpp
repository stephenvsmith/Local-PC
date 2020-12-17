#include "sharedFunctions.h"
using namespace Rcpp;

// SKELETON HELPER FUNCTIONS

/*
 * This function returns a vector showing the current graph edges coming from node i in the graph
 * being worked on. Assumes a complete graph
 * Tested: 12/16/20
 */
// [[Rcpp::export]]
NumericVector get_current_edges(int i,int p,NumericMatrix graph){
  NumericVector current_edges;
  for (int j=0;j<p;++j){
    if (graph(i,j)==1){
      current_edges.push_back(j);
    }
  }
  return current_edges;
}

/*
 * This function generates all possible combinations of a vector x of length l and puts them in a matrix
 * It uses the R function combn to produce the results.
 * Tested: 12/16/20
 */
// [[Rcpp::export]]
NumericMatrix combn_cpp(NumericVector x,int l){
  int xlength = x.length();
  NumericMatrix result;
  
  if (l > x.length()){
    Rcout << "There aren't enough neighbors for the current value of l\n";
  }
  
  if (l==0){
    result = NumericMatrix(1,1);
    result(0,0) = NA_REAL;
    return result;
  } else if (xlength==1 & l == 1){
    result = NumericMatrix(1,1);
    result(0,0) = x(0);
    return result;
  } else if (l>=1 & xlength>1) {
    Function f("combn");
    result = f(Named("x")=x,_["m"]=l);
    return result;
  } else {
    Rcout << "The value of l is negative: " << l << std::endl;
    Rcout << "Or the value of l is less than the size of the vector. Length of input vector: " << x.length() << std::endl;
  }
  
  return result;
}


// TRUE DAG HELPER FUNCTION

/*
 * This function has been tested.
 * File: test_skeletonhelpers.R
 */

/*
 * This function helps us to get the neighbors from the true DAG.
 * Input: target node i, number of nodes p, and the true DAG
 * Returns: Vector of neighbors 
 * Reviewed: 12/16/20
 */
// [[Rcpp::export]]
NumericVector get_neighbors_from_dag(int i,int p,NumericMatrix true_dag){
  NumericVector neighbors;
  NumericVector parents;
  NumericVector children;
  
  for (int j = 0;j<p;++j){
    if (true_dag(j,i)==1){ 
      parents.push_back(j);
      //Rcout << "Call from get_neighbors_from_dag. Node " << j << " is a parent.\n";
    } else if (true_dag(i,j)==1){
      children.push_back(j);
      //Rcout << "Call from get_neighbors_from_dag. Node " << j << " is a child.\n";
    }
  }
  
  NumericVector potential_spouses;
  int current_val;
  for (NumericVector::iterator it = children.begin(); it != children.end(); ++it){
    //Rcout << "Call from get_neighbors_from_dag. We are evaluating the following child: " << *it << std::endl;
    for (int j = 0; j<p; ++j){
      current_val = true_dag(j,*it);
      if (current_val == 1 & i != j){
        potential_spouses.push_back(j);
        //Rcout << "Call from get_neighbors_from_dag. Node " << j << " is a potential spouse.\n";
      }
    }
  }
  
  neighbors = union_(parents,children);
  neighbors = union_(neighbors,potential_spouses);
  std::sort(neighbors.begin(),neighbors.end());
  
  return neighbors;
}


// PRINTING FUNCTIONS (DEBUGGING)

/*
 * 
 * Functions for printing different values 
 * 
 */

void print_vector_elements(NumericVector v,StringVector names,String opening,String closing){
  int ln = v.length();
  Rcout << opening.get_cstring();
  for (int i = 0;i<ln;++i) {
    Rcout << names(v(i)) << " ";
  }
  Rcout << closing.get_cstring();
}


void print_vector_elements_nonames(NumericVector v,String opening,String closing,String sep){
  int l = v.length();
  Rcout << opening.get_cstring();
  for (int i = 0;i<l;++i) {
    Rcout << v(i) << sep.get_cstring();
  }
  Rcout << closing.get_cstring();
}


void print_matrix(NumericMatrix m){
  int n = m.nrow();
  int p = m.ncol();
  String ending;
  for (int i = 0;i<n;++i){
    for (int j = 0;j<p;++j){
      if (j==p-1){
        ending = "\n";
      } else{
        ending = " ";
      }
      Rcout << m(i,j) << ending.get_cstring();
    }
  }
}

void print_S_vals(List S){
  List sublist;
  for (int i=0;i<S.length();++i){
    sublist = S[i];
    for (int j=0;j<sublist.length();++j){
      Rcout << "S[[" << i << "]][[" << j << "]] = ";
      print_vector_elements_nonames(sublist[j],"",""," ");
      Rcout << " ";
    }
    Rcout << std::endl;
  }
}

void iteration_print(const int &l,const int &i,const int &j,const NumericVector &sep,const StringVector &names,const double &pval){
  Rcout << "l: " << l << " | i: " << i << " | j: " << j << " | k: ";
  if (l == 0){
    Rcout << sep;
  } else {
    print_vector_elements(sep,names);
  }
  Rcout << " | p-val: " << pval;
  Rcout << std::endl;
}
