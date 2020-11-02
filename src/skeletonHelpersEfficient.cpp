#include "printFunctions.h"
#include "pCorTest.h"
#include "skeletonHelpers.h"
#include "trueDAGinfo.h"
using namespace Rcpp;


/*
 * This function helps us to add a separating set for nodes i and j
 */
// [[Rcpp::export]]
void change_S_efficient(List &S,int i,int j,NumericVector sep){
  NumericVector sep_new;
  sep_new = clone(sep);
  String i_char = String((char) i);
  String j_char = String((char) j);
  List sublist;
  sublist = S[i_char];
  sublist[j_char] = sep_new;
  S[i_char] = sublist;
}

/*
 * Allows us to change S to indicate that nodes i and j are separated without a separating set
 */
// [[Rcpp::export]]
void change_S_0_efficient(List &S,int i,int j){
  String i_char = String((char) i);
  String j_char = String((char) j);
  List sublist0;
  sublist0 = S[i_char];
  sublist0[j_char] = -1;
  S[i_char] = sublist0;
}

// [[Rcpp::export]]
NumericVector get_potential_sep(const int &i,const int &j,const NumericVector &neighborhood,const int &N,const NumericMatrix &true_dag){
  //Rcout << "We are finding the neighbors of " << neighborhood(i) << " and " << neighborhood(j) << std::endl;
  NumericVector potential_seps = union_(get_neighbors_from_dag(neighborhood(i),N,true_dag),
                                        get_neighbors_from_dag(neighborhood(j),N,true_dag));
  NumericVector exceptions = NumericVector::create(neighborhood(i),neighborhood(j));
  //print_vector_elements_nonames(potential_seps);
  potential_seps = setdiff(potential_seps,exceptions);
  //print_vector_elements_nonames(potential_seps);
  return potential_seps;
}

/*
 * 
 * This function checks whether or not nodes i and j are separated by any of the 
 * sets in the matrix kvals for a given significance level
 * 
 */
void check_separation_sample_efficient(const int &l,const int &i,const int &j,
                      const NumericMatrix &kvals,
                      NumericVector &sep,NumericMatrix true_dag,
                      const StringVector &names,const NumericVector &neighborhood,
                      NumericMatrix C, // Passed by reference automatically
                      List S,double &pval,arma::mat &R,int &n, // R is the correlation matrix
                      double &signif_level,bool &verbose){
  int k;
  bool keep_checking_k;
  List test_result;
  
  arma::uvec sep_arma;
  
  if (l == 0){
    sep = NA_REAL;
    //Rcout << "Size of arma vector when l=0: " << sep_arma.size() << std::endl;
    test_result = condIndTest(R,neighborhood(i),neighborhood(j),sep_arma,n,signif_level);
    pval = test_result["pval"];
    //pval = as<double>(get_pval(i,j,true_dag,names));
    //Rcout << "The p-value is " << pval << std::endl;
    if (test_result["result"]){ // We have conditional independence established
      change_S_0_efficient(S,neighborhood(i),neighborhood(j));
      change_S_0_efficient(S,neighborhood(j),neighborhood(i));
      
      C(i,j) = 0;
      C(j,i) = 0;
    }
  } else {
    k = 0;
    keep_checking_k = true;
    while (keep_checking_k & (k<kvals.cols())){
      sep = kvals( _ , k );
      sep_arma = as<arma::uvec>(sep);
      test_result = condIndTest(R,neighborhood(i),neighborhood(j),sep_arma,n,signif_level);
      pval = test_result["pval"];
      //pval = as<double>(get_pval(i,j,true_dag,names,sep));
      if (verbose){
        Rcout << "The p-value is " << pval << std::endl;
      }
      if (test_result["result"]){
        if (verbose){
          Rcout << names(neighborhood(i)) << " is separated from " << names(neighborhood(j)) << " by node(s):\n";
          print_vector_elements(sep,names);
        }
        change_S_efficient(S,neighborhood(i),neighborhood(j),sep);
        change_S_efficient(S,neighborhood(j),neighborhood(i),sep);
        C(i,j) = 0;
        C(j,i) = 0;
        keep_checking_k = false;
        if (verbose){
          print_matrix(kvals);
        }
      }
      ++k;
    }
  }
}
