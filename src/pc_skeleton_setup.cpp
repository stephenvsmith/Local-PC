#include <Rcpp.h>
using namespace Rcpp;

/*
 * This function sets up the nested lists that will hold separating sets
 */
// [[Rcpp::export]]
List create_conditioning_sets_cpp(int p){
  List S(p);
  for (int i=0;i<p;++i){
    List sublist = List(p);
    for (int j=0;j<p;++j){
      sublist[j] = NA_REAL;
    }
    S[i] = sublist;
  }
  return S;
}

/*
 * The following function sets up the basic data structures for the skeleton algorithm
 */
// [[Rcpp::export]]
List pc_pop_skeleton_setup_cpp(NumericMatrix true_dag,StringVector names,int lmax,bool verbose){
  // Number of nodes
  int p = 0;
  p = true_dag.nrow();

  if (verbose){
    Rcout << "There are " << p << " nodes in the DAG.\n";
  }

  NumericMatrix C_tilde(p,p);
  std::fill(C_tilde.begin(), C_tilde.end(), 1);
  C_tilde.fill_diag(0);
  if (verbose){
    Rcout << "Our starting matrix is " << C_tilde.nrow() << "x" << C_tilde.ncol() << ".\n";
  }

  List S = create_conditioning_sets_cpp(p);

  std::vector<double> p_vals;

  return List::create(
    _["p"] = p,
    _["C_tilde"]=C_tilde,
    _["true_dag"]=true_dag,
    _["names"]=names,
    _["lmax"]=lmax,
    _["S"]=S,
    _["verbose"]=verbose,
    _["p_vals"]=p_vals);

}


/*
 * This function helps us to get the neighbors from the true DAG
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

  return neighbors;
}

/*
 * This function returns a vector showing the current graph edges coming from node i
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
 * This function generates all possible combinations of a vector x
 */
// [[Rcpp::export]]
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

// [[Rcpp::export]]
void print_vector_elements(NumericVector v,StringVector names, String opening="",String closing=""){
  int l = v.length();
  Rcout << opening.get_cstring();
  for (int i = 0;i<l;++i) {
    Rcout << names(v(i)) << " ";
  }
  Rcout << closing.get_cstring();
}

// [[Rcpp::export]]
void print_vector_elements_nonames(NumericVector v,String opening="",String closing="",String sep=" "){
  int l = v.length();
  Rcout << opening.get_cstring();
  for (int i = 0;i<l;++i) {
    Rcout << v(i) << sep.get_cstring();
  }
  Rcout << closing.get_cstring();
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
void print_S_vals(List S){
  List sublist;
  for (int i=0;i<S.length();++i){
    sublist = S[i];
    for (int j=0;j<sublist.length();++j){
      Rcout << "S[[" << i << "]][[" << j << "]] = ";
      print_vector_elements_nonames(sublist[j]);
      Rcout << " ";
    }
    Rcout << std::endl;
  }
}

/*
 * This function helps us to add a separating set for nodes i and j
 */
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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

// [[Rcpp::export]]
void check_separation(const int &l,const int &i,const int &j,const NumericMatrix &kvals,Function get_pval,NumericVector &sep,NumericMatrix true_dag,const StringVector &names,NumericMatrix C,List S,double &pval){
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
      Rcout << "The p-value is " << pval << std::endl;
      if (pval==1){
        Rcout << names(i) << " is separated from " << names(j) << " by node(s):\n";
        print_vector_elements(sep,names);
        change_S(S,i,j,sep);
        change_S(S,j,i,sep);
        C(i,j) = 0;
        C(j,i) = 0;
        keep_checking_k = false;
        print_matrix(kvals);
      }
      ++k;
    }
  }
}


// [[Rcpp::export]]
List pc_pop_get_skeleton_cpp(List var_list){
  int l = -1;
  int lmax = var_list["lmax"];
  int p = var_list["p"];
  bool verbose = var_list["verbose"];
  NumericMatrix C_tilde = var_list["C_tilde"];
  NumericMatrix C = clone(C_tilde);
  NumericMatrix true_dag = var_list["true_dag"];
  StringVector names = var_list["names"];
  List S = var_list["S"];

  NumericVector exceptions; //= {0,0};
  NumericVector neighbors;
  NumericVector edges_i;
  NumericVector sep;

  NumericMatrix kvals;
  double pval=0.0;

  Environment myEnv = Environment::namespace_env("LocalPC");
  Function get_pval= myEnv["get_pval"];

  int num_tests=0;

  while (l < lmax){
    l += 1;
    if (verbose){
      Rcout << "The value of l is " << l << std::endl;
    }

    for (int i=0;i<p;++i){
      if (verbose){
        Rcout << "The value of i is " << i << std::endl;
      }
      // Work through potential neighbors with separating set of size l
      edges_i = get_current_edges(i,p,C);
      for (NumericVector::iterator it = edges_i.begin(); it != edges_i.end(); ++it){
        int j = *it;
        if (j != i){

          if (verbose){
            Rcout << "The value of j is " << j << std::endl;
          }
          // Find neighbors of i and j
          neighbors = union_(get_neighbors_from_dag(i,p,true_dag),get_neighbors_from_dag(j,p,true_dag));
          exceptions = NumericVector::create(i,j);
          //Rcout << "Exceptions: " << exceptions(0) << " and " << exceptions(1) << std::endl;

          //print_vector_elements(neighbors,names,"Neighbors before:","=============");

          neighbors = setdiff(neighbors,exceptions);

          //print_vector_elements(neighbors,names,"Neighbors after:","=============");

          if (neighbors.length()>= l){

            if (verbose){
              Rcout << "There are " << neighbors.length() << " neighbor(s).\n";
            }
            kvals = combn_cpp(neighbors,l);

            check_separation(l,i,j,kvals,get_pval,sep,true_dag,names,C,S,pval);

            Rcout << "l: " << l << " | i: " << i << " | j: " << j << " | k: ";
            print_vector_elements(sep,names);
            Rcout << " | p-val: " << pval;
            Rcout << std::endl;
            //print_S_vals(S);
            ++num_tests;
          }
        }
      }

    }

  }
  Rcout << "The final C matrix:\n";
  print_matrix(C);
  Rcout << "Conclusion of algorithm\n";

  return var_list;
}

// [[Rcpp::export]]
List pc_pop_skeleton_cpp(NumericMatrix true_dag,StringVector names,int lmax=3,bool verbose=true,bool verbose_small=true){

  List var_list = pc_pop_skeleton_setup_cpp(true_dag,names,lmax,verbose);


  List final_skeleton_list = pc_pop_get_skeleton_cpp(var_list);
  return var_list;
}

//[[Rcpp::export]]
NumericVector test(NumericVector final,NumericVector remove){
  return setdiff(final,remove);
}

// // [[Rcpp::export]]
// List pc_skeleton_setup_cpp(DataFrame data,NumericMatrix true_dag,NumericMatrix C_tilde,bool pop, int lmax,bool verbose,bool verbose_small,double tol){
//   // Number of nodes
//   int p = 0;
//   if (pop){
//     p = true_dag.nrow();
//     std::fill(C_tilde.begin(), C_tilde.end(), 1);
//     C_tilde.fill_diag(0);
//   } else {
//     p = data.ncol();
//   }
//
//   return List::create(p,C_tilde);
// }


