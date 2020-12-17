#include "sharedFunctions.h"
#include "pCorTest.h"
using namespace Rcpp;

/*
 * This function sets up the nested lists that will hold separating sets
 * It will do this efficiently by only setting up lists for the neighbors of interest to the algorithm
 * Tested: 12/16/20
 */
void create_conditioning_sets_efficient_cpp(List &S,NumericVector &neighbors){
  int l = neighbors.length();
  String i_char;
  String j_char;
  List sublist;
  for (int i=0;i<l;++i){
    sublist = List(); // Creating the sublist for neighbor i
    for (int j=0;j<l;++j){
      if (j != i){
        j_char = String((char) neighbors(j));
        sublist[j_char] = NA_REAL;
      }
    }
    i_char = String((char) neighbors(i));
    S[i_char] = sublist;
  }
}



/*
 * This function sets up the nested lists that will hold separating sets
 * It will do this efficiently by only setting up lists for the neighbors of interest to the algorithm
 * This is the version for testing
 * Tested: 12/16/20
 */
// [[Rcpp::export]]
List create_conditioning_sets_efficient_cpp2(NumericVector &neighbors){
  int l = neighbors.length();
  String i_char;
  String j_char;
  List S(0);
  List sublist;
  for (int i=0;i<l;++i){
    sublist = List(); // Creating the sublist for neighbor i
    for (int j=0;j<l;++j){
      if (j != i){
        j_char = String((char) neighbors(j));
        sublist[j_char] = NA_REAL;
      }
    }
    i_char = String((char) neighbors(i));
    S[i_char] = sublist;
  }
  return S;
}

/*
 * The following function sets up the basic data structures for the skeleton algorithm that is used
 * for the sample version of the algorithm
 * Tested: 12/16/20
 */
// [[Rcpp::export]]
List pc_sample_skeleton_setup_efficient_cpp(NumericMatrix &true_dag,
                                            const int &target,
                                            StringVector &names,
                                            const int &lmax,
                                            bool &verbose){
  int p;
  p = true_dag.nrow();
  
  // Find the neighborhood of the target node
  NumericVector neighbors = get_neighbors_from_dag(target,p,true_dag);
  neighbors.push_front(target);
  std::sort(neighbors.begin(),neighbors.end());
  
  int N = neighbors.size();
  
  if (verbose){
    Rcout << "There are " << p << " nodes in the DAG.\n";
    Rcout << "There are " << N << " nodes in the neighborhood.\n";
  }
  
  // Initial graph that will be modified through the process of the algorithm
  NumericMatrix C_tilde(N);
  std::fill(C_tilde.begin(), C_tilde.end(), 1);
  C_tilde.fill_diag(0);
  
  if (verbose){
    Rcout << "Our starting matrix is " << C_tilde.nrow() << "x" << C_tilde.ncol() << ".\n";
  }
  
  // Create the list that will store 
  List S(0);
  create_conditioning_sets_efficient_cpp(S,neighbors);
  
  std::vector<double> p_vals;
  
  return List::create(
    _["p"] = p,
    _["C_tilde"]=C_tilde,
    _["true_dag"]=true_dag,
    _["names"]=names,
    _["neighborhood"]=neighbors,
    _["lmax"]=lmax,
    _["S"]=S,
    _["verbose"]=verbose,
    _["p_vals"]=p_vals);
  
}

/*
 * This function helps us to add a separating set for nodes i and j
 * Tested: 12/16/20
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
 * Tested: 12/16/20
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

/* 
 * This function returns nodes that may separate node i and node j 
 * i and j are with respect to the neighborhood of interest
 * That is, if i=0 and j=2 and the neighborhood is {4,5,6}, then 
 * i referes to node 4 and j refers to node 6
 * Tested: 12/16/20
 */
// [[Rcpp::export]]
NumericVector get_potential_sep(const int &i,const int &j,const NumericVector &neighborhood,const int &p,const NumericMatrix &true_dag){
  //Rcout << "We are finding the neighbors of " << neighborhood(i) << " and " << neighborhood(j) << std::endl;
  NumericVector potential_seps = union_(get_neighbors_from_dag(neighborhood(i),p,true_dag),
                                        get_neighbors_from_dag(neighborhood(j),p,true_dag));
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
 * Tested: 12/16/20
 */
// [[Rcpp::export]]
void check_separation_sample_efficient(const int &l,const int &i,const int &j,
                      const NumericMatrix &kvals,
                      NumericVector &sep,NumericMatrix true_dag,
                      const StringVector &names,const NumericVector &neighborhood,
                      NumericMatrix C, // Passed by reference automatically
                      List S,double &pval,int &num_tests,arma::mat &R,int &n, // R is the correlation matrix
                      double &signif_level,bool &verbose){
  int k;
  int kp = kvals.cols();
  bool keep_checking_k;
  List test_result;
  
  arma::uvec sep_arma;
  
  if (l == 0){
    sep = NA_REAL;
    //Rcout << "Size of arma vector when l=0: " << sep_arma.size() << std::endl;
    test_result = condIndTest(R,neighborhood(i),neighborhood(j),sep_arma,n,signif_level);
    ++num_tests;
    pval = test_result["pval"];
    //pval = as<double>(get_pval(i,j,true_dag,names));
    if (verbose){
      Rcout << "The p-value is " << pval << std::endl;
    }
    if (test_result["result"]){ // We have conditional independence established
      if (verbose){
        Rcout << names(neighborhood(i)) << " is separated from " << names(neighborhood(j)) << std::endl;
        //Rcout << "i = " << i << " | j = " << j << " | N_i = " << neighborhood(i) << " | N_j = " << neighborhood(j) << std::endl;
      }
      change_S_0_efficient(S,neighborhood(i),neighborhood(j));
      change_S_0_efficient(S,neighborhood(j),neighborhood(i));
      
      C(i,j) = 0;
      C(j,i) = 0;
    }
  } else {
    k = 0;
    keep_checking_k = true;
    while (keep_checking_k & (k<kp)){
      sep = kvals( _ , k );
      sep_arma = as<arma::uvec>(sep);
      test_result = condIndTest(R,neighborhood(i),neighborhood(j),sep_arma,n,signif_level);
      ++num_tests;
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
      } else {
        if (verbose){
          Rcout << names(neighborhood(i)) << " is NOT separated from " << names(neighborhood(j)) << " by node(s):\n";
          print_vector_elements(sep,names);
        }
        
      }
      ++k; // Increment to obtain the next potential separating set
    }
  }
}
