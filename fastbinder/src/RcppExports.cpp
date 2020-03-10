// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// I_binder_naive
long I_binder_naive(Rcpp::IntegerMatrix CI, int version, double Const_Binder);
RcppExport SEXP _fastbinder_I_binder_naive(SEXP CISEXP, SEXP versionSEXP, SEXP Const_BinderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type CI(CISEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    Rcpp::traits::input_parameter< double >::type Const_Binder(Const_BinderSEXP);
    rcpp_result_gen = Rcpp::wrap(I_binder_naive(CI, version, Const_Binder));
    return rcpp_result_gen;
END_RCPP
}
// I_binder_opt
long I_binder_opt(Rcpp::IntegerMatrix CI, int version, double Const_Binder);
RcppExport SEXP _fastbinder_I_binder_opt(SEXP CISEXP, SEXP versionSEXP, SEXP Const_BinderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type CI(CISEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    Rcpp::traits::input_parameter< double >::type Const_Binder(Const_BinderSEXP);
    rcpp_result_gen = Rcpp::wrap(I_binder_opt(CI, version, Const_Binder));
    return rcpp_result_gen;
END_RCPP
}
// I_binder_parallel
long I_binder_parallel(Rcpp::IntegerMatrix CI, int version, double Const_Binder);
RcppExport SEXP _fastbinder_I_binder_parallel(SEXP CISEXP, SEXP versionSEXP, SEXP Const_BinderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type CI(CISEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    Rcpp::traits::input_parameter< double >::type Const_Binder(Const_BinderSEXP);
    rcpp_result_gen = Rcpp::wrap(I_binder_parallel(CI, version, Const_Binder));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastbinder_I_binder_naive", (DL_FUNC) &_fastbinder_I_binder_naive, 3},
    {"_fastbinder_I_binder_opt", (DL_FUNC) &_fastbinder_I_binder_opt, 3},
    {"_fastbinder_I_binder_parallel", (DL_FUNC) &_fastbinder_I_binder_parallel, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastbinder(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
