#include "printFunctions.h"
#include "pCorTest.h"
using namespace Rcpp;

/*
 * This function returns a vector showing the current graph edges coming from node i in the graph
 * being worked on
 */
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
 * This function generates all possible combinations of a vector x and puts them in a matrix
 */
NumericMatrix combn_cpp(NumericVector x,int l){
  
  NumericMatrix result;
  if (l==0){
    result = NumericMatrix(1,1);
    result(0,0) = NA_REAL;
    return result;
  } else if (l==1){
    result = NumericMatrix(1,1);
    result(0,0) = x(0);
    return result;
  } else if (l>1) {
    Function f("combn");
    result = f(Named("x")=x,_["m"]=l);
    return result;
  } else {
    Rcout << "The value of l is negative: " << l << std::endl;
  }
  
  return result;
}



/*
 * This function helps us to add a separating set for nodes i and j
 */
List change_S(List S,int i,int j,NumericVector sep){
  NumericVector sep_new;
  sep_new = clone(sep);
  // Rcout << "S before:\n";
  // print_S_vals(S);
  List sublist;
  sublist = S[i];
  sublist[j] = sep_new;
  S[i] = sublist;
  // Rcout << "S after:\n";
  // print_S_vals(S);
  
  return S;
}

/*
 * Allows us to change S to indicate that nodes i and j are separated without a separating set
 */
List change_S_0(List S,int i,int j){
  
  // Rcout << "S before:\n";
  // print_S_vals(S);
  List sublist0;
  sublist0 = S[i];
  sublist0[j] = -1;
  S[i] = sublist0;
  // Rcout << "S after:\n";
  // print_S_vals(S);
  
  return S;
}

void check_separation(const int &l,const int &i,const int &j,
                      const NumericMatrix &kvals,Function get_pval,
                      NumericVector &sep,NumericMatrix true_dag,
                      const StringVector &names,NumericMatrix C,
                      List S,double &pval,bool &verbose){
  int k;
  bool keep_checking_k;
  
  if (l == 0){
    sep = NA_REAL;
    pval = as<double>(get_pval(i,j,true_dag,names));
    //Rcout << "The p-value is " << pval << std::endl;
    if (pval == 1){
      change_S_0(S,i,j);
      change_S_0(S,j,i);
      
      C(i,j) = 0;
      C(j,i) = 0;
    }
  } else {
    k = 0;
    keep_checking_k = true;
    while (keep_checking_k & (k<kvals.cols())){
      sep = kvals( _ , k );
      pval = as<double>(get_pval(i,j,true_dag,names,sep));
      if (verbose){
        Rcout << "The p-value is " << pval << std::endl;
      }
      if (pval==1){
        if (verbose){
          Rcout << names(i) << " is separated from " << names(j) << " by node(s):\n";
          print_vector_elements(sep,names);
        }
        change_S(S,i,j,sep);
        change_S(S,j,i,sep);
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

void check_separation_sample(const int &l,const int &i,const int &j,
                      const NumericMatrix &kvals,
                      NumericVector &sep,NumericMatrix true_dag,
                      const StringVector &names,NumericMatrix C,
                      List S,double &pval,arma::mat &R,int &n,
                      double &signif_level,bool &verbose){
  int k;
  bool keep_checking_k;
  List test_result;
  
  arma::uvec sep_arma;
  
  if (l == 0){
    sep = NA_REAL;
    Rcout << "Size of arma vector when l=0: " << sep_arma.size() << std::endl;
    test_result = condIndTest(R,i,j,sep_arma,n,signif_level);
    pval = test_result["pval"];
    //pval = as<double>(get_pval(i,j,true_dag,names));
    //Rcout << "The p-value is " << pval << std::endl;
    if (test_result["result"]){
      change_S_0(S,i,j);
      change_S_0(S,j,i);
      
      C(i,j) = 0;
      C(j,i) = 0;
    }
  } else {
    k = 0;
    keep_checking_k = true;
    while (keep_checking_k & (k<kvals.cols())){
      sep = kvals( _ , k );
      sep_arma = as<arma::uvec>(sep);
      test_result = condIndTest(R,i,j,sep_arma,n,signif_level);
      pval = test_result["pval"];
      //pval = as<double>(get_pval(i,j,true_dag,names,sep));
      if (verbose){
        Rcout << "The p-value is " << pval << std::endl;
      }
      if (test_result["result"]){
        if (verbose){
          Rcout << names(i) << " is separated from " << names(j) << " by node(s):\n";
          print_vector_elements(sep,names);
        }
        change_S(S,i,j,sep);
        change_S(S,j,i,sep);
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
