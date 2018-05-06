// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <iostream>
#include <RcppArmadillo.h>
#include <ctime>
using namespace Rcpp;
using namespace std;
using namespace arma;

// PROTOTYPES: #################################################################

// Loads/Unloads all variables into memory
double loadVars(arma::mat X, arma::mat Y, arma::mat V, arma::mat D, int cType);
void unloadVars();

// These find one best split for continuous and categorical (unorderable) predictors.
List findOneBestContinuousSplit(arma::uvec INDEX, int pIndex, int minBucket);
List findOneBestCategoricalSplit(arma::uvec INDEX, int pIndex);

// Main iterator that finds all best splits for all predictors
List findAllBestSplits(arma::uvec INDEX, double currentGoodness, int k);

// Main function:
List fTreeRPart(arma::uvec INDEX, double currentGoodness, int height);

// Boostrap:
List fTreeBootstrap(arma::mat BOOTINDEX);

// Cost functions:
double computeGoodness(arma::uvec INDEX);
double sseCost(arma::uvec INDEX);
double mahalanobisCost(arma::uvec INDEX);
double l2Cost(arma::uvec INDEX);
double wssCost(arma::uvec INDEX);
double rdsCost(arma::uvec INDEX);
double sql2Cost(arma::uvec INDEX);

// Helpers for categorical predictors:
arma::mat getAllCombinations(int N);
std::string combToString(arma::vec combination, arma::vec categories);
arma::uvec getCombinIndices(arma::vec Xcov, arma::vec Combination, arma::vec categories);

// Global Variables: ###########################################################
arma::mat Xcov;
arma::uvec covType; // This is the new feature for categorical predictors
arma::mat Yout;
arma::mat Vinv;
arma::mat Dist;
int    costType;
int    minSPLIT;
int    minBUCKET;
double minDEVIANCE;
int    kP; // the number of predictors to consider in random forest type of tree.
double ArgStep;

// [[Rcpp::export]]
double loadVars(arma::mat X, arma::uvec covarType, arma::mat Y, arma::mat V, arma::mat D,
                int nPredictors, int cType = 1, int mSplit = 15,
                int mBucket = 10, double cp = 0.01, double argStep = 1){

  Xcov        = X;
  covType     = covarType; // new categorical feature!
  Yout        = Y;
  Vinv        = V;
  Dist        = D;
  costType    = cType;

  minSPLIT    = mSplit;
  minBUCKET   = mBucket;
  kP          = nPredictors;

  ArgStep     = argStep;

  int nP = Y.n_cols ;
  arma::uvec INDEX = conv_to<uvec>::from(linspace(0, nP - 1, nP)) ;
  minDEVIANCE      = cp * computeGoodness(INDEX) ;


  return 0;

}

// [[Rcpp::export]]
void unloadVars(){
    Xcov.clear();
    covType.clear(); // a new capability for categorical predictors.
    Yout.clear();
    Vinv.clear();
    Dist.clear();
    costType = 0;
    minBUCKET = 0;
    minDEVIANCE = 0;
    kP = 0;
}

// MAIN TREE fitting routines. To upgrade for unorderable categorical one must
// write a new function and patch it into "loadVars" and "findAllBestSplits"
// note that "splitvalue" in that case becomes a string.

// [[Rcpp::export]]
List findOneBestContinuousSplit(arma::uvec INDEX, int pIndex){

  arma::vec curCovariate = Xcov.col(pIndex); // copies from global
  curCovariate = curCovariate.elem(INDEX);

  List ret;                    // output list
  ret["splitPoint"] = std::string("NA");  // 1 / 0.0; // if a split is not found it will return infinity.

  arma::vec unq     = unique(curCovariate);

  if(unq.n_elem == 1) {
    return(ret);
  }

  // this is where you will check for unorderable predictors and handle appropriately
  // important modification is that "splitPoint" changes to become a string.
  // the rest remains the same.

  arma::vec delta   = diff(unq)/2;
  arma::vec splits  = unq.subvec(0, unq.n_elem - 2) + delta;
  arma::vec goodnessVec(3); goodnessVec.ones(); // 0 - total, 1 - left; 2 - right;
  goodnessVec = goodnessVec/0.0;

  for(int i = 0; i < splits.n_elem; i++){

    arma::uvec leftInd  = INDEX(arma::find(curCovariate < splits(i)));
    arma::uvec rightInd = INDEX(arma::find(curCovariate >= splits(i)));

    if(leftInd.n_elem >= minBUCKET && rightInd.n_elem >= minBUCKET){

      double leftGoodness  = computeGoodness(leftInd);
      double rightGoodness = computeGoodness(rightInd);

      if((leftGoodness + rightGoodness) < goodnessVec(0)){ // if it improves save it.

        goodnessVec(0) = leftGoodness + rightGoodness;
        goodnessVec(1) = leftGoodness;
        goodnessVec(2) = rightGoodness;

        ret["leftInd"]    = leftInd;
        ret["rightInd"]   = rightInd;
        ret["goodness"]   = goodnessVec;
        ret["splitPoint"] = std::to_string(splits(i));   //  has to convert to string

      }
    }
  }

  return ret;

}

// [[Rcpp::export]]
List findOneBestCategoricalSplit(arma::uvec INDEX, int pIndex){

  arma::vec curCovariate = Xcov.col(pIndex); // copies from global
  curCovariate = curCovariate.elem(INDEX);

  List ret;                    // output list
  ret["splitPoint"] = std::string("NA");  // 1 / 0.0; // if a split is not found it will return infinity.

  arma::vec unq     = unique(curCovariate);

  //cout << unq << endl;

  //cout << "FOBCS pt1" << endl;

  if(unq.n_elem == 1) {
    return(ret);
  }

  //cout << "FOBCS pt2" << endl;

  // this is where you will check for unorderable predictors and handle appropriately
  // important modification is that "splitPoint" changes to become a string.
  // the rest remains the same.

  //cout << "FOBCS pt3" << endl;

  arma::mat splits  = getAllCombinations(unq.n_elem);

  //cout<<splits<<endl;

  //cout << splits <<endl;

  arma::vec goodnessVec(3); goodnessVec.ones(); // 0 - total, 1 - left; 2 - right;
  goodnessVec = goodnessVec/0.0;

  //cout << "FOBCS pt4" << endl;

  for(int i = 0; i < splits.n_rows; i++){

    arma::uvec flags    = getCombinIndices(curCovariate, splits.row(i).t(), unq);
    arma::uvec leftInd  = INDEX(arma::find(flags == 1));
    arma::uvec rightInd = INDEX(arma::find(flags == 0));

    //cout << "FOBCS pt6" << endl;

    if(leftInd.n_elem >= minBUCKET && rightInd.n_elem >= minBUCKET){

      double leftGoodness  = computeGoodness(leftInd);
      double rightGoodness = computeGoodness(rightInd);

      if((leftGoodness + rightGoodness) < goodnessVec(0)){ // if it improves save it.

        //cout << "FOBCS pt7" << endl;

        goodnessVec(0) = leftGoodness + rightGoodness;
        goodnessVec(1) = leftGoodness;
        goodnessVec(2) = rightGoodness;

        ret["leftInd"]    = leftInd;
        ret["rightInd"]   = rightInd;
        ret["goodness"]   = goodnessVec;

        //cout << "FOBCS pt8" << endl;
        //cout << splits.row(i)<<endl;

        ret["splitPoint"] = combToString(splits.row(i).t(), unq);   //  has to convert to string
        //cout << "FOBCS pt9" << endl;

      }
    }
  }

  return ret;

}

// [[Rcpp::export]]
List findAllBestSplits(arma::uvec INDEX, double currentGoodness){

  int nP = Xcov.n_cols;
  arma::uvec currentPredictors = conv_to<uvec>::from(linspace(0, nP - 1, nP)); // All predictors!
  arma::mat goodness(3, nP); goodness.fill(datum::nan); // row 0 = total; row 1 = left; row 2 = right;
  double bestGoodness = currentGoodness;
  List ret;
  List splitData;

  //cout << "FABS pt1" << endl;

  ret["bestPredictor"] = -1;

  if(kP != nP){ // modification for random forest approach:
    std::random_shuffle(currentPredictors.begin(), currentPredictors.end() );
    currentPredictors = currentPredictors.subvec(0, kP - 1);
  }

  // Loop through predictors
  for(int i = 0; i < currentPredictors.size(); i++){

    int cIndex = currentPredictors(i);
    //cout << "FABS pt2" << endl;

    if(covType(cIndex) == 1) { // by convention, we code categorical as 1, continuous as 0.

      // here we will modify for unorderable predictors.
      //cout << "FABS pt2a" << endl;

      splitData = findOneBestCategoricalSplit(INDEX, cIndex);

    } else {
      //cout << "FABS pt2b" << endl;

      // Continuous predictors:
      splitData = findOneBestContinuousSplit(INDEX, cIndex);

    }
    //if(R_finite(splitData["splitPoint"]) == 1){ // if the best split was found start checking.

    std::string splitPoint = splitData["splitPoint"];  // Its a string! enhorabuena! :)
    std::string notFound ("NA");

    if(notFound.compare(splitPoint) != 0){

      goodness.col(cIndex) = as<vec>(splitData["goodness"]) ; // always save goodness for sensitivity

      if(bestGoodness > goodness(0, cIndex)){ // if it improves save it.
        ret = splitData ;
        ret["bestPredictor"] = cIndex;
        bestGoodness = goodness(0, cIndex);
      }

    }

  }

ret["allGoodness"] = goodness ;

return ret;

}

// [[Rcpp::export]]
List fTreeRPart(arma::uvec INDEX, double currentGoodness, int depth){

  List ret ;
  ret["indices"]        = INDEX ;
  ret["nodeGoodness"]   = currentGoodness ;
  ret["isLeaf"]         = 1 ;                  // 1 = yes, 0 = no

  ret["nLeaves"]        = 1 ;                  // visualization parameter.
  ret["depth"]          = depth ;              // visualization parameter.
  ret["midpoint"]       = 1 ;                  // visualization parameter.

  if(INDEX.n_elem > minSPLIT){

    //cout << "Point1 made" << endl;

    List allSplits      = findAllBestSplits(INDEX, currentGoodness) ;

    //cout << "Point1a made" << endl;

    int bestPredictor   = allSplits["bestPredictor"] ;

    if(bestPredictor == -1) { // it happens that the best predictor/best split cannot be found.
      return(ret);
    }

    arma::vec cGoodness   = allSplits["goodness"] ;
    //cout << "Point2 made" << endl;

    if( (currentGoodness - cGoodness(0)) > minDEVIANCE ){   // if it improves continue

      ret["isLeaf"]         = 0 ;
      ret["splitsGoodness"] = allSplits["allGoodness"] ;    // a matrix with the goodness of all parameters
      ret["bestPredictor"]  = allSplits["bestPredictor"] ;  // winning predictor
      ret["splitPoint"]     = allSplits["splitPoint"] ;     // winning split (for categorical: A string of left categories!)

      // Recursion:----------------
      List left             = fTreeRPart(allSplits["leftInd"]  , cGoodness(1), depth + 1) ;
      List right            = fTreeRPart(allSplits["rightInd"] , cGoodness(2), depth + 1) ;

      ret["left"]           = left ;
      ret["right"]          = right ;

      // All visualization parameters from here:
      int nLeftLeaves       = left["nLeaves"] ;
      int nRightLeaves      = right["nLeaves"] ;

      ret["nLeaves"]        = nLeftLeaves + nRightLeaves ;
      ret["midpoint"]       = (2 * nLeftLeaves + nRightLeaves) / 2.0 ;

    }
  }

  return(ret) ;

}

// [[Rcpp::export]]
List fTreeBootstrap(arma::mat BOOTINDEX){

  // it should run in parallel.
  int nBoot = BOOTINDEX.n_rows;
  List ret(nBoot);
  for(int j = 0; j < nBoot; j++){
    arma::uvec cIndex = conv_to<uvec>::from(BOOTINDEX.row(j));
    double cGoodness = computeGoodness(cIndex);
    ret(j) = fTreeRPart(cIndex, cGoodness, 0);
  }
  return(ret);
}

/////////////////////////////// SPLIT FUNCTIONS ////////////////////////////////
// [[Rcpp::export]]
double computeGoodness(arma::uvec INDEX){

  double goodness = 0;

  switch(costType){

    case 1: goodness = sseCost(INDEX);         break;
    case 2: goodness = mahalanobisCost(INDEX); break;
    case 3: goodness = l2Cost(INDEX);          break;
    case 4: goodness = wssCost(INDEX);         break;
    case 5: goodness = rdsCost(INDEX);         break;
    case 6: goodness = sql2Cost(INDEX);        break;
    default: cout << "computeGoodness ERROR: cost type was not initialized properly" << endl;

  }

return goodness;

}

// [[Rcpp::export]]
double sseCost(arma::uvec INDEX){

  arma::mat currentMat = Yout.cols(INDEX);
  arma::vec rowSums(Yout.n_rows);

  for(int i = 0; i < currentMat.n_rows; i++){
    rowSums(i) = ArgStep * accu(square(currentMat.row(i) - arma::mean(currentMat.row(i))));
  }

  rowSums(0)        *= 0.5;
  rowSums(currentMat.n_rows - 1) *= 0.5;

  return accu(rowSums);

}

/*
double mahalanobisCost(arma::uvec INDEX){

  arma::mat currentMat = Yout.cols(INDEX);
  arma::vec rowSums(Yout.n_rows);
  arma::mat cVinv = Vinv;

  for(int i = 0; i < currentMat.n_rows; i++){
    currentMat.row(i) -= arma::mean(currentMat.row(i));
  }

  arma::mat outVec = arma::trans(currentMat)*cVinv*currentMat;

  return trace(outVec);

}
*/

// [[Rcpp::export]]
double mahalanobisCost(arma::uvec INDEX){

  arma::mat currentMat = Yout.cols(INDEX);
  arma::vec outVec(currentMat.n_cols);
  arma::mat cVinv = Vinv;

  arma::vec meanVector = currentMat * arma::ones(currentMat.n_cols);
  meanVector /= (double)currentMat.n_cols;

  for(int i = 0; i < currentMat.n_cols; i++){
    arma::mat X = currentMat.col(i) - meanVector;
    arma::mat Y = X.t() * cVinv;
    outVec(i) = as_scalar(Y * X);
  }

  return sum(outVec) ;

}



// [[Rcpp::export]]
double l2Cost(arma::uvec INDEX){

  arma::mat currentMat = Yout.cols(INDEX);
  arma::vec l2norms(currentMat.n_cols);
  l2norms.zeros();

  for(int i = 0; i < currentMat.n_rows; i++){
    l2norms += arma::trans(square(currentMat.row(i) - arma::mean(currentMat.row(i))));
  }

  return accu(sqrt(l2norms));
}

// [[Rcpp::export]]
double wssCost(arma::uvec INDEX){

  arma::mat currentMat = Dist.submat(INDEX, INDEX);
  arma::vec vones(currentMat.n_cols);

  vones.ones();

  return(as_scalar((arma::trans(vones) * currentMat * vones)));

}

// [[Rcpp::export]]
double sql2Cost(arma::uvec INDEX){

  arma::mat currentMat = Yout.cols(INDEX);
  arma::vec rowSums(Yout.n_rows);

  for(int i = 0; i < currentMat.n_rows; i++){
    rowSums(i) = accu(square(currentMat.row(i) - arma::mean(currentMat.row(i))));
  }

  return accu(rowSums);

}

// [[Rcpp::export]]
double rdsCost(arma::uvec INDEX){

  arma::mat currentMat = Dist.submat(INDEX, INDEX);
  arma::vec vones(currentMat.n_cols);

  vones.ones();

  arma::vec medSum = currentMat * vones;

  return(as_scalar(medSum.min()));

}

////////////////////////////// HELPER FUNCTIONS ////////////////////////////////
// [[Rcpp::export]]
arma::mat getAllCombinations(int N)
{

  vec elems  = arma::linspace<vec>(1, N, N);
  mat outMat(pow(2,N) - 1, N, fill::zeros);

  for (int s = 1; s < (1 << elems.n_elem); ++s)
  {
    for (int e = 0; e < elems.n_elem; ++e)
    {
      if (s & (1 << e))
      {
        outMat(s - 1, e) = 1;
      }
    }
  }
  return outMat;
}

// [[Rcpp::export]]
std::string combToString(arma::vec combination, arma::vec categories){

  std::string ret;

  for(int i = 0; i < combination.n_elem; i++){
    if(combination(i) == 1){
      ret +=  "," + std::to_string(int(categories(i)));
    }
  }

  ret.erase(ret.begin());

  return ret;
}

// [[Rcpp::export]]
arma::uvec getCombinIndices(arma::vec Xcov, arma::vec Combination, arma::vec categories){

  //cout << "GCIND pt1" << endl;

  arma::uvec ret(Xcov.n_elem); ret.zeros();

  for(int i = 0; i < Combination.n_elem; i++){
    if(Combination(i) == 1) {
      ret.elem(arma::find(Xcov == categories(i))).ones();
    }
  }

 // return arma::find(ret == lr);
 return ret ;
}

