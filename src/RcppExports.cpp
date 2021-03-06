// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// loadVars
double loadVars(arma::mat X, arma::uvec covarType, arma::mat Y, arma::mat V, arma::mat D, int nPredictors, int cType, int mSplit, int mBucket, double cp, double argStep);
RcppExport SEXP _fTree_loadVars(SEXP XSEXP, SEXP covarTypeSEXP, SEXP YSEXP, SEXP VSEXP, SEXP DSEXP, SEXP nPredictorsSEXP, SEXP cTypeSEXP, SEXP mSplitSEXP, SEXP mBucketSEXP, SEXP cpSEXP, SEXP argStepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type covarType(covarTypeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type nPredictors(nPredictorsSEXP);
    Rcpp::traits::input_parameter< int >::type cType(cTypeSEXP);
    Rcpp::traits::input_parameter< int >::type mSplit(mSplitSEXP);
    Rcpp::traits::input_parameter< int >::type mBucket(mBucketSEXP);
    Rcpp::traits::input_parameter< double >::type cp(cpSEXP);
    Rcpp::traits::input_parameter< double >::type argStep(argStepSEXP);
    rcpp_result_gen = Rcpp::wrap(loadVars(X, covarType, Y, V, D, nPredictors, cType, mSplit, mBucket, cp, argStep));
    return rcpp_result_gen;
END_RCPP
}
// unloadVars
void unloadVars();
RcppExport SEXP _fTree_unloadVars() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    unloadVars();
    return R_NilValue;
END_RCPP
}
// findOneBestContinuousSplit
List findOneBestContinuousSplit(arma::uvec INDEX, int pIndex);
RcppExport SEXP _fTree_findOneBestContinuousSplit(SEXP INDEXSEXP, SEXP pIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    Rcpp::traits::input_parameter< int >::type pIndex(pIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(findOneBestContinuousSplit(INDEX, pIndex));
    return rcpp_result_gen;
END_RCPP
}
// findOneBestCategoricalSplit
List findOneBestCategoricalSplit(arma::uvec INDEX, int pIndex);
RcppExport SEXP _fTree_findOneBestCategoricalSplit(SEXP INDEXSEXP, SEXP pIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    Rcpp::traits::input_parameter< int >::type pIndex(pIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(findOneBestCategoricalSplit(INDEX, pIndex));
    return rcpp_result_gen;
END_RCPP
}
// findAllBestSplits
List findAllBestSplits(arma::uvec INDEX, double currentGoodness);
RcppExport SEXP _fTree_findAllBestSplits(SEXP INDEXSEXP, SEXP currentGoodnessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    Rcpp::traits::input_parameter< double >::type currentGoodness(currentGoodnessSEXP);
    rcpp_result_gen = Rcpp::wrap(findAllBestSplits(INDEX, currentGoodness));
    return rcpp_result_gen;
END_RCPP
}
// fTreeRPart
List fTreeRPart(arma::uvec INDEX, double currentGoodness, int depth);
RcppExport SEXP _fTree_fTreeRPart(SEXP INDEXSEXP, SEXP currentGoodnessSEXP, SEXP depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    Rcpp::traits::input_parameter< double >::type currentGoodness(currentGoodnessSEXP);
    Rcpp::traits::input_parameter< int >::type depth(depthSEXP);
    rcpp_result_gen = Rcpp::wrap(fTreeRPart(INDEX, currentGoodness, depth));
    return rcpp_result_gen;
END_RCPP
}
// fTreeBootstrap
List fTreeBootstrap(arma::mat BOOTINDEX);
RcppExport SEXP _fTree_fTreeBootstrap(SEXP BOOTINDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type BOOTINDEX(BOOTINDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(fTreeBootstrap(BOOTINDEX));
    return rcpp_result_gen;
END_RCPP
}
// computeGoodness
double computeGoodness(arma::uvec INDEX);
RcppExport SEXP _fTree_computeGoodness(SEXP INDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(computeGoodness(INDEX));
    return rcpp_result_gen;
END_RCPP
}
// sseCost
double sseCost(arma::uvec INDEX);
RcppExport SEXP _fTree_sseCost(SEXP INDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(sseCost(INDEX));
    return rcpp_result_gen;
END_RCPP
}
// mahalanobisCost
double mahalanobisCost(arma::uvec INDEX);
RcppExport SEXP _fTree_mahalanobisCost(SEXP INDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(mahalanobisCost(INDEX));
    return rcpp_result_gen;
END_RCPP
}
// l2Cost
double l2Cost(arma::uvec INDEX);
RcppExport SEXP _fTree_l2Cost(SEXP INDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(l2Cost(INDEX));
    return rcpp_result_gen;
END_RCPP
}
// wssCost
double wssCost(arma::uvec INDEX);
RcppExport SEXP _fTree_wssCost(SEXP INDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(wssCost(INDEX));
    return rcpp_result_gen;
END_RCPP
}
// sql2Cost
double sql2Cost(arma::uvec INDEX);
RcppExport SEXP _fTree_sql2Cost(SEXP INDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(sql2Cost(INDEX));
    return rcpp_result_gen;
END_RCPP
}
// rdsCost
double rdsCost(arma::uvec INDEX);
RcppExport SEXP _fTree_rdsCost(SEXP INDEXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type INDEX(INDEXSEXP);
    rcpp_result_gen = Rcpp::wrap(rdsCost(INDEX));
    return rcpp_result_gen;
END_RCPP
}
// getAllCombinations
arma::mat getAllCombinations(int N);
RcppExport SEXP _fTree_getAllCombinations(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(getAllCombinations(N));
    return rcpp_result_gen;
END_RCPP
}
// combToString
std::string combToString(arma::vec combination, arma::vec categories);
RcppExport SEXP _fTree_combToString(SEXP combinationSEXP, SEXP categoriesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type combination(combinationSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type categories(categoriesSEXP);
    rcpp_result_gen = Rcpp::wrap(combToString(combination, categories));
    return rcpp_result_gen;
END_RCPP
}
// getCombinIndices
arma::uvec getCombinIndices(arma::vec Xcov, arma::vec Combination, arma::vec categories);
RcppExport SEXP _fTree_getCombinIndices(SEXP XcovSEXP, SEXP CombinationSEXP, SEXP categoriesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Xcov(XcovSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Combination(CombinationSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type categories(categoriesSEXP);
    rcpp_result_gen = Rcpp::wrap(getCombinIndices(Xcov, Combination, categories));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fTree_loadVars", (DL_FUNC) &_fTree_loadVars, 11},
    {"_fTree_unloadVars", (DL_FUNC) &_fTree_unloadVars, 0},
    {"_fTree_findOneBestContinuousSplit", (DL_FUNC) &_fTree_findOneBestContinuousSplit, 2},
    {"_fTree_findOneBestCategoricalSplit", (DL_FUNC) &_fTree_findOneBestCategoricalSplit, 2},
    {"_fTree_findAllBestSplits", (DL_FUNC) &_fTree_findAllBestSplits, 2},
    {"_fTree_fTreeRPart", (DL_FUNC) &_fTree_fTreeRPart, 3},
    {"_fTree_fTreeBootstrap", (DL_FUNC) &_fTree_fTreeBootstrap, 1},
    {"_fTree_computeGoodness", (DL_FUNC) &_fTree_computeGoodness, 1},
    {"_fTree_sseCost", (DL_FUNC) &_fTree_sseCost, 1},
    {"_fTree_mahalanobisCost", (DL_FUNC) &_fTree_mahalanobisCost, 1},
    {"_fTree_l2Cost", (DL_FUNC) &_fTree_l2Cost, 1},
    {"_fTree_wssCost", (DL_FUNC) &_fTree_wssCost, 1},
    {"_fTree_sql2Cost", (DL_FUNC) &_fTree_sql2Cost, 1},
    {"_fTree_rdsCost", (DL_FUNC) &_fTree_rdsCost, 1},
    {"_fTree_getAllCombinations", (DL_FUNC) &_fTree_getAllCombinations, 1},
    {"_fTree_combToString", (DL_FUNC) &_fTree_combToString, 2},
    {"_fTree_getCombinIndices", (DL_FUNC) &_fTree_getCombinIndices, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fTree(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
