// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bsts_posterior_predict
Rcpp::List bsts_posterior_predict(Rcpp::NumericVector model, arma::rowvec Z, arma::mat G, arma::mat Q, arma::mat R, arma::mat init_state, arma::mat X, arma::mat B, arma::mat gpriors, arma::mat qpriors, arma::vec vpriors, arma::mat U);
RcppExport SEXP _tsforeign_bsts_posterior_predict(SEXP modelSEXP, SEXP ZSEXP, SEXP GSEXP, SEXP QSEXP, SEXP RSEXP, SEXP init_stateSEXP, SEXP XSEXP, SEXP BSEXP, SEXP gpriorsSEXP, SEXP qpriorsSEXP, SEXP vpriorsSEXP, SEXP USEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init_state(init_stateSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gpriors(gpriorsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type qpriors(qpriorsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vpriors(vpriorsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    rcpp_result_gen = Rcpp::wrap(bsts_posterior_predict(model, Z, G, Q, R, init_state, X, B, gpriors, qpriors, vpriors, U));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tsforeign_bsts_posterior_predict", (DL_FUNC) &_tsforeign_bsts_posterior_predict, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_tsforeign(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}