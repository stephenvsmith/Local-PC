// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pc_pop_get_skeleton_cpp
List pc_pop_get_skeleton_cpp(List var_list);
RcppExport SEXP _LocalPC_pc_pop_get_skeleton_cpp(SEXP var_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type var_list(var_listSEXP);
    rcpp_result_gen = Rcpp::wrap(pc_pop_get_skeleton_cpp(var_list));
    return rcpp_result_gen;
END_RCPP
}
// pc_pop_cpp
List pc_pop_cpp(NumericMatrix true_dag, StringVector names, int lmax, bool verbose, bool verbose_small);
RcppExport SEXP _LocalPC_pc_pop_cpp(SEXP true_dagSEXP, SEXP namesSEXP, SEXP lmaxSEXP, SEXP verboseSEXP, SEXP verbose_smallSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type true_dag(true_dagSEXP);
    Rcpp::traits::input_parameter< StringVector >::type names(namesSEXP);
    Rcpp::traits::input_parameter< int >::type lmax(lmaxSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose_small(verbose_smallSEXP);
    rcpp_result_gen = Rcpp::wrap(pc_pop_cpp(true_dag, names, lmax, verbose, verbose_small));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LocalPC_pc_pop_get_skeleton_cpp", (DL_FUNC) &_LocalPC_pc_pop_get_skeleton_cpp, 1},
    {"_LocalPC_pc_pop_cpp", (DL_FUNC) &_LocalPC_pc_pop_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_LocalPC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
