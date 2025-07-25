// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fastSigLm
List fastSigLm(const arma::vec& y, const arma::mat& X);
RcppExport SEXP _FusioMR_fastSigLm(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastSigLm(y, X));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_memo_joint
List gibbs_memo_joint(int niter, NumericVector Gamma_hat_1, NumericVector Gamma_hat_2, NumericVector gamma_hat_1, NumericVector gamma_hat_2, NumericVector s_hat_Gamma_1, NumericVector s_hat_Gamma_2, NumericVector s_hat_gamma_1, NumericVector s_hat_gamma_2);
RcppExport SEXP _FusioMR_gibbs_memo_joint(SEXP niterSEXP, SEXP Gamma_hat_1SEXP, SEXP Gamma_hat_2SEXP, SEXP gamma_hat_1SEXP, SEXP gamma_hat_2SEXP, SEXP s_hat_Gamma_1SEXP, SEXP s_hat_Gamma_2SEXP, SEXP s_hat_gamma_1SEXP, SEXP s_hat_gamma_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Gamma_hat_1(Gamma_hat_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Gamma_hat_2(Gamma_hat_2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma_hat_1(gamma_hat_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma_hat_2(gamma_hat_2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_hat_Gamma_1(s_hat_Gamma_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_hat_Gamma_2(s_hat_Gamma_2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_hat_gamma_1(s_hat_gamma_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_hat_gamma_2(s_hat_gamma_2SEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_memo_joint(niter, Gamma_hat_1, Gamma_hat_2, gamma_hat_1, gamma_hat_2, s_hat_Gamma_1, s_hat_Gamma_2, s_hat_gamma_1, s_hat_gamma_2));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_semo_nohp
List gibbs_semo_nohp(int niter, NumericVector Gamma_hat_1, NumericVector Gamma_hat_2, NumericVector gamma_hat, NumericVector s2_hat_Gamma_1, NumericVector s2_hat_Gamma_2, NumericVector s2_hat_gamma);
RcppExport SEXP _FusioMR_gibbs_semo_nohp(SEXP niterSEXP, SEXP Gamma_hat_1SEXP, SEXP Gamma_hat_2SEXP, SEXP gamma_hatSEXP, SEXP s2_hat_Gamma_1SEXP, SEXP s2_hat_Gamma_2SEXP, SEXP s2_hat_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Gamma_hat_1(Gamma_hat_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Gamma_hat_2(Gamma_hat_2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma_hat(gamma_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2_hat_Gamma_1(s2_hat_Gamma_1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2_hat_Gamma_2(s2_hat_Gamma_2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2_hat_gamma(s2_hat_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_semo_nohp(niter, Gamma_hat_1, Gamma_hat_2, gamma_hat, s2_hat_Gamma_1, s2_hat_Gamma_2, s2_hat_gamma));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_seso_nohp
NumericVector gibbs_seso_nohp(int niter, NumericVector Gamma_hat, NumericVector gamma_hat, NumericVector s2_hat_Gamma, NumericVector s2_hat_gamma);
RcppExport SEXP _FusioMR_gibbs_seso_nohp(SEXP niterSEXP, SEXP Gamma_hatSEXP, SEXP gamma_hatSEXP, SEXP s2_hat_GammaSEXP, SEXP s2_hat_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Gamma_hat(Gamma_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma_hat(gamma_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2_hat_Gamma(s2_hat_GammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2_hat_gamma(s2_hat_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_seso_nohp(niter, Gamma_hat, gamma_hat, s2_hat_Gamma, s2_hat_gamma));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_seso_uhp_only
NumericVector gibbs_seso_uhp_only(int niter, NumericVector Gamma_hat, NumericVector gamma_hat, NumericVector s2_hat_Gamma, NumericVector s2_hat_gamma);
RcppExport SEXP _FusioMR_gibbs_seso_uhp_only(SEXP niterSEXP, SEXP Gamma_hatSEXP, SEXP gamma_hatSEXP, SEXP s2_hat_GammaSEXP, SEXP s2_hat_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Gamma_hat(Gamma_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma_hat(gamma_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2_hat_Gamma(s2_hat_GammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2_hat_gamma(s2_hat_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_seso_uhp_only(niter, Gamma_hat, gamma_hat, s2_hat_Gamma, s2_hat_gamma));
    return rcpp_result_gen;
END_RCPP
}
// my_rinvgamma
NumericVector my_rinvgamma(int n, double shape, double rate);
RcppExport SEXP _FusioMR_my_rinvgamma(SEXP nSEXP, SEXP shapeSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(my_rinvgamma(n, shape, rate));
    return rcpp_result_gen;
END_RCPP
}
// my_rinvwishart
arma::mat my_rinvwishart(double nu, arma::mat S);
RcppExport SEXP _FusioMR_my_rinvwishart(SEXP nuSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(my_rinvwishart(nu, S));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _FusioMR_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// my_rdirichlet
NumericVector my_rdirichlet(int n, NumericVector alpha);
RcppExport SEXP _FusioMR_my_rdirichlet(SEXP nSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(my_rdirichlet(n, alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FusioMR_fastSigLm", (DL_FUNC) &_FusioMR_fastSigLm, 2},
    {"_FusioMR_gibbs_memo_joint", (DL_FUNC) &_FusioMR_gibbs_memo_joint, 9},
    {"_FusioMR_gibbs_semo_nohp", (DL_FUNC) &_FusioMR_gibbs_semo_nohp, 7},
    {"_FusioMR_gibbs_seso_nohp", (DL_FUNC) &_FusioMR_gibbs_seso_nohp, 5},
    {"_FusioMR_gibbs_seso_uhp_only", (DL_FUNC) &_FusioMR_gibbs_seso_uhp_only, 5},
    {"_FusioMR_my_rinvgamma", (DL_FUNC) &_FusioMR_my_rinvgamma, 3},
    {"_FusioMR_my_rinvwishart", (DL_FUNC) &_FusioMR_my_rinvwishart, 2},
    {"_FusioMR_mvrnormArma", (DL_FUNC) &_FusioMR_mvrnormArma, 3},
    {"_FusioMR_my_rdirichlet", (DL_FUNC) &_FusioMR_my_rdirichlet, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_FusioMR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
