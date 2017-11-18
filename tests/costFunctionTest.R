# This script tests cost functions.
# Author Ogy Grujic
# Date: October 16th 2016
library(Rcpp)
library(functInterp)
library(RcppArmadillo)
library(functInterp)
library(Matrix)
library(testthat)

# sourceCpp('./src/main.cpp')
load("~/Google Drive/Functional_Research/TreeS_workingDir/TreeS_workingDir/RealDataSetScaled.Rdata")

plotBoth <- function(.ftreeObj, .nodeLocator = NULL, .index = 1){

  A <- plotFtree(.ftreeObj, .index = .index, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 2.4, .node = .nodeLocator, .ggReturn = TRUE)
  B <- plotNodeData(.ftreeObj, .treeIndex = .index, .node = .nodeLocator, .ggR = TRUE)
  multiPlot(A,B, cols=2)
}


# Setup: Subject to Change:
COVARIATES <- as.matrix(cbind(covariatesCompletion, covariatesGeological))
FUNCTIONS  <- as.matrix(completeSet)
SIGMA      <- as.matrix(nearPD(cov(t(FUNCTIONS)))$mat)
DISTANCE   <- as.matrix(dist(t(FUNCTIONS)))
EMPTY      <- matrix(NA, 2, 2)

context("Testing Cpp functions")

fTree:::loadVars(COVARIATES, FUNCTIONS, solve(SIGMA), DISTANCE, 22, cType = 4,
                 mSplit =10, mBucket=3, cp=0.0001)

# TEST Sse: ####################################################################
test_that("SSE cost function works", {
  cppSSE <- fTree:::sseCost((1:ncol(FUNCTIONS))-1)
  colMean <- rowMeans(FUNCTIONS)
  Rsse <- sum(apply(FUNCTIONS, 2, function(x) trapzc(1, (x-colMean)^2)))
  expect_that(cppSSE, equals(Rsse))
  })

# TEST Mahalanobis: ############################################################
test_that("mahalanobis cost function works", {
  cppMahalanobis <- fTree:::mahalanobisCost(1:ncol(FUNCTIONS)-1)
  Q <- FUNCTIONS - apply(FUNCTIONS, 1, mean)
  Rmahalanobis <- sum(diag(t(Q) %*% solve(SIGMA) %*% Q))
  expect_that(cppMahalanobis, equals(Rmahalanobis))
})

# TEST L2 norm: ################################################################
test_that("L2 cost function works", {
  cppL2 <- fTree:::l2Cost(1:ncol(FUNCTIONS)-1)
  Rl2 <- FUNCTIONS - apply(FUNCTIONS, 1, mean)
  Rl2 <- sqrt(apply(Rl2^2, 2, sum))
  expect_that(cppL2, equals(sum(Rl2)))
})

# Test WSS: ####################################################################
test_that("WSS cost function works", {
  cppWSS <- fTree:::wssCost(1:nrow(DISTANCE)-1)
  expect_that(sum(DISTANCE), equals(cppWSS))
})

# Test RDS: ####################################################################
test_that("RDS cost function works", {
  cppWSS <- fTree:::rdsCost(1:nrow(DISTANCE)-1)
  expect_that(min(apply(DISTANCE,1,sum)), equals(cppWSS))
})

# TEST Goodness computation: ###################################################
test_that("Goodness computation SSE", {
  costType = 1
  fTree:::loadVars(COVARIATES, FUNCTIONS, solve(SIGMA), DISTANCE, 22, costType)
  expect_that(fTree:::computeGoodness(1:30-1), equals(fTree:::sseCost((1:30)-1)))
})

test_that("Goodness computation mahalanobis", {
  costType = 2
  fTree:::loadVars(COVARIATES, FUNCTIONS, solve(SIGMA), DISTANCE, 22, costType)
  expect_that(fTree:::computeGoodness(1:30-1), equals(fTree:::mahalanobisCost((1:30)-1)))
})

test_that("Goodness computation L2norm", {
  costType = 3
  fTree:::loadVars(COVARIATES, FUNCTIONS, solve(SIGMA), DISTANCE, 22, costType)
  expect_that(fTree:::computeGoodness(1:30-1), equals(fTree:::l2Cost((1:30)-1)))
})

test_that("Goodness computation WSS", {
  costType = 4
  fTree:::loadVars(COVARIATES, FUNCTIONS, solve(SIGMA), DISTANCE, 22, costType)
  expect_that(fTree:::computeGoodness(1:30-1), equals(fTree:::wssCost((1:30)-1)))
})

# TEST findOneBestSplit: #######################################################
oneBestList <- fTree:::findOneBestContinuousSplit(0:171, 7L)
str(oneBestList,1)

# TEST findAllBestSplits: ######################################################
cGood <- fTree:::computeGoodness(0:171)
allBestList <- fTree:::findAllBestSplits(0:171, cGood)
str(allBestList, 1)

# TEST recursion: ##############################################################
cGood     <- fTree:::computeGoodness(0:171)
system.time(rpartList <- fTree:::fTreeRPart(0:171, cGood, 0))

# TEST bootstrap: ##############################################################
nBoot = 100
BOOTINDEX <- matrix(sample(ncol(FUNCTIONS), ncol(FUNCTIONS) * nBoot, replace=TRUE), ncol = ncol(FUNCTIONS), nrow=nBoot)
system.time(myBoot <- fTree:::fTreeBootstrap(BOOTINDEX - 1))

# TEST R functions: ############################################################
context("Testing R functions")

expect_that(ftree(.Y = FUNCTIONS), throws_error('Covariates were not provided!'))
expect_that(ftree(.X = COVARIATES, .Y = FUNCTIONS, cost.type = "wss"),
            throws_error("For cost type WSS you MUST provide a distance matrix or distance type (i.e. euclidean). "))
expect_that(ftree(.X = COVARIATES, .D = 'euclidean', cost.type = "wss"),
            throws_error("Functions were not provided! Please provide functions and distance type or just a distance matrix."))
expect_that(ftree(.X = COVARIATES, .D = 34, cost.type = "wss"),
            throws_error(".D is not a matrix nor a string. Exiting!"))
expect_that(ftree(.X = COVARIATES, .Y = FUNCTIONS, nP = 34, tree.type = "randomforest"),
            throws_error("nP cannot be larger than the number of covariates in .X.\nRandom forest gives the best result for: nP = sqrt(ncol(.X))"))

train <- sample(ncol(FUNCTIONS), 100, replace = FALSE)
system.time(rfTest <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "randomforest", cost.type = "sse", nBoot = 100))
system.time(bagTest <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], .D = "minkowski", tree.type = "bagging", cost.type = "sse", nBoot = 100))

singleTestDeuclidean  <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], .D = "euclidean", tree.type = "single", cost.type = "wss")
singleTestSSE         <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "single", cost.type = "sse")
singleTestMahalanobis <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "single", cost.type = "mahalanobis")
singleTestL2          <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "single", cost.type = "l2norm")
singleTestRDS         <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], .D = "euclidean", tree.type = "single", cost.type = "rdist")

# Plotting Test: Single Test WSS------------------------------------------------
plotFtree(singleTestDeuclidean, .horizontal = TRUE, .ylimLow = -0.5, .labSize = 3.5, .node = "0")
plotNodeData(singleTestDeuclidean, .node = "0")

plotFtree(singleTestDeuclidean, .horizontal = TRUE, .ylimLow = -0.5, .labSize = 3.5, .node = "0r")
plotNodeData(singleTestDeuclidean, .node = "0r")

plotFtree(singleTestDeuclidean, .horizontal = TRUE, .ylimLow = -0.5, .labSize = 3.5, .node = "0rl")
plotNodeData(singleTestDeuclidean, .node = "0rl")

plotFtree(singleTestDeuclidean, .horizontal = TRUE, .ylimLow = -0.5, .labSize = 3.5, .node = "0rlr")
plotNodeData(singleTestDeuclidean, .node = "0rlr")

# Plotting Test: Single TEST SSE------------------------------------------------
plotFtree(singleTestSSE, .horizontal = TRUE, .ylimLow = -2.5,  .labSize = 3.5, .node = "0")
plotNodeData(singleTestSSE, .node = "0")

plotFtree(singleTestSSE, .horizontal = TRUE, .ylimLow = -2.5,  .labSize = 3.5, .node = "0l")
plotNodeData(singleTestSSE, .node = "0l")

plotFtree(singleTestSSE, .horizontal = TRUE, .ylimLow = -2.5,  .labSize = 3.5, .node = "0lr")
plotNodeData(singleTestSSE, .node = "0lr")

plotFtree(singleTestSSE, .horizontal = TRUE, .ylimLow = -2.5,  .labSize = 3.5, .node = "0lrl")
plotNodeData(singleTestSSE, .node = "0lrl")

plotFtree(singleTestSSE, .horizontal = TRUE, .ylimLow = -2.5,  .labSize = 3.5, .node = "0ll")
plotNodeData(singleTestSSE, .node = "0ll")

plotFtree(singleTestSSE, .horizontal = TRUE, .ylimLow = -2.5,  .labSize = 3.5, .node = "0llr")
plotNodeData(singleTestSSE, .node = "0llr")

plotFtree(singleTestSSE, .horizontal = TRUE, .ylimLow = -2.5,  .labSize = 3.5, .node = "0llrr")
plotNodeData(singleTestSSE, .node = "0llrr")

# Plotting Test: Single Test Mahalanobis----------------------------------------
plotFtree(singleTestMahalanobis, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0")
plotNodeData(singleTestMahalanobis, .node = "0")

plotFtree(singleTestMahalanobis, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0l")
plotNodeData(singleTestMahalanobis, .node = "0l")

plotFtree(singleTestMahalanobis, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lr")
plotNodeData(singleTestMahalanobis, .node = "0lr")

plotFtree(singleTestMahalanobis, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0ll")
plotNodeData(singleTestMahalanobis, .node = "0ll")

plotFtree(singleTestMahalanobis, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lll")
plotNodeData(singleTestMahalanobis, .node = "0lll")

plotFtree(singleTestMahalanobis, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llll")
plotNodeData(singleTestMahalanobis, .node = "0llll")

# Plotting Test: Single Test L2 Norm--------------------------------------------
plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0")
plotNodeData(singleTestL2, .node = "0")

plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0l")
plotNodeData(singleTestL2, .node = "0l")

plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lr")
plotNodeData(singleTestL2, .node = "0lr")

plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lrr")
plotNodeData(singleTestL2, .node = "0lrr")

plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lrrl")
plotNodeData(singleTestL2, .node = "0lrrl")

plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lrrlr")
plotNodeData(singleTestL2, .node = "0lrrlr")

plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lrrlrr")
plotNodeData(singleTestL2, .node = "0lrrlrr")

plotFtree(singleTestL2, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0lrrlrrr")
plotNodeData(singleTestL2, .node = "0lrrlrrr")

# Plotting Test: Single Test RDS------------------------------------------------
plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0")
plotNodeData(singleTestRDS, .node = "0")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0l")
plotNodeData(singleTestRDS, .node = "0l")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0ll")
plotNodeData(singleTestRDS, .node = "0ll")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llr")
plotNodeData(singleTestRDS, .node = "0llr")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llrr")
plotNodeData(singleTestRDS, .node = "0llrr")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llrrr")
plotNodeData(singleTestRDS, .node = "0llrrr")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llrrl")
plotNodeData(singleTestRDS, .node = "0llrrl")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llrrlr")
plotNodeData(singleTestRDS, .node = "0llrrlr")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llrrlrl")
plotNodeData(singleTestRDS, .node = "0llrrlrl")

plotFtree(singleTestRDS, .horizontal = TRUE, .ylimLow = -1.5, .labSize = 3.5, .node = "0llrrlrlr")
plotNodeData(singleTestRDS, .node = "0llrrlrlr")

# Plotting Test: Random Forest--------------------------------------------------
plotFtree(rfTest, .horizontal = TRUE, .labSize = 3.5, .ylimLow = -1)
plotFtree(rfTest, .index = 2, .horizontal = TRUE, .labSize = 3.5, .ylimLow = -1)
plotFtree(rfTest, .index = 1, .horizontal = TRUE, .labSize = 3.5, .ylimLow = -1)

# Plotting Test: Bagging--------------------------------------------------------
plotFtree(bagTest, .horizontal = TRUE, .labSize = 3.5, .ylimLow = -1)
plotFtree(bagTest, .index = 2, .horizontal = TRUE, .labSize = 3.5, .ylimLow = -1)
plotFtree(bagTest, .index = 1, .horizontal = TRUE, .labSize = 3.5, .ylimLow = -1)

# Predictions TEST--------------------------------------------------------------
singleTestSSE.predict         <- predictFtree(singleTestSSE, as.data.frame(COVARIATES[-train,]))
singleTestMahalanobis.predict <- predictFtree(singleTestMahalanobis, as.data.frame(COVARIATES[-train,]))
singleTestL2.predict          <- predictFtree(singleTestL2, as.data.frame(COVARIATES[-train,]))

rfTest.predict <- predictFtree(rfTest, as.data.frame(COVARIATES[-train,]))
bagTest.predict <- predictFtree(bagTest, as.data.frame(COVARIATES[-train,]))

# Plotting the results:
index = 10
plot(FUNCTIONS[,index], col="black", lty = 1, lwd = 5, type = "l", ylim = c(0,500))
lines(Reduce(c,singleTestSSE.predict[[index]]), col = "red", lwd = 3, lty=2)
legend("topright", c("True", "Forecast"), col=c("black","red"), lty=c(1,2), lwd=c(5,3))


index = 15
par(mfrow=c(1,2))
matplot(FUNCTIONS[,train], type = "l", col="blue", xlab="Time (days)", ylab="Oil Rate (stb/day)")
matplot(t(Reduce(cbind, rfTest.predict[[index]])), col = "red", lwd = 1, lty=1, type="l", ylim = c(0,500), add=TRUE)
lines(FUNCTIONS[,-train][,index], col="black", lty = 1, lwd = 5, type = "l", ylim = c(0,500))
lines(Reduce(c,singleTestSSE.predict[[index]]), col="black", lty = 3, lwd = 5, type = "l", ylim = c(0,500))
legend("topright", c("Train", "True", "RF forecast", "Single Tree"), col=c("blue","black","red", "black"), lty=c(1,1,2,3), lwd=c(1,5,3,3))

matplot(FUNCTIONS[,train], type = "l", col="blue", xlab="Time (days)", ylab="Oil Rate (stb/day)")
matplot(t(Reduce(cbind, bagTest.predict[[index]])), col = "red", lwd = 1, lty=1, type="l", ylim = c(0,500), add=TRUE)
lines(FUNCTIONS[,-train][,index], col="black", lty = 1, lwd = 5, type = "l", ylim = c(0,500))
lines(Reduce(c,singleTestSSE.predict[[index]]), col="black", lty = 3, lwd = 5, type = "l", ylim = c(0,500))
legend("topright", c("Train", "True", "Bagging forecast", "Single Tree"), col=c("blue","black","red", "black"), lty=c(1,1,2,3), lwd=c(1,5,3,3))

# Variography:------------------------------------------------------------------
library(functInterp)

fg <- fstat(NULL, "Raw Data", as.data.frame(scaleInput(COVARIATES[train,23:24])), as.data.frame(FUNCTIONS[,train]))
fg <- estimateDrift("~.", fg)
fg <- fvariogram("~.",
                 fg,
                 Nlags  = 100,
                 LagMax = 1,
                 directions = "omni",
                 useResidual = FALSE,
                 ArgStep = 1)
plotVariogram(fg)

Residual <- as.matrix(
  ldply(predictFtree(rfTest, as.data.frame(COVARIATES[train,])),
        function(x) {colMeans(x, na.rm = TRUE)}, .id = NULL)) - t(FUNCTIONS[,train])

fg2 <- fstat(NULL, "Residual", as.data.frame(COVARIATES[train,23:24]+0.001), as.data.frame(t(Residual)))
fg2 <- estimateDrift("~.", fg2, .type = "OLS")
fg2 <- fvariogram("~.",
                 fg2,
                 Nlags  = 100,
                 LagMax = 1,
                 directions = "all",
                 angTolerance = 40,
                 useResidual = FALSE,
                 ArgStep = 1)
fg2 <- fitVariograms(fg2, model = vgm(3e+05, "Sph", 0.1, 0))
plotVariogram(fg2)
fg2 <- addCovariance(fg2)

fg2.predict <- predictFstat(fg2, as.data.frame(COVARIATES[-train,23:24]), .what = "Residual", .type = "OK")

FinalHybrid <- as.matrix(
  ldply(predictFtree(rfTest, as.data.frame(COVARIATES[-train,])),
        function(x) {colMeans(x, na.rm = TRUE)}, .id = NULL)) - t(fg2.predict$Forecast)

RF.Forecast <- ldply(predictFtree(rfTest, as.data.frame(COVARIATES[-train,])),
                     function(x) {colMeans(x, na.rm = TRUE)}, .id = NULL)

matplot(t(FinalHybrid), type = "l", col="blue")
matplot(t(RF.Forecast), type="l", col="red", add=TRUE)


for(index in 1:70){
M <- cbind(t(FinalHybrid)[,index],
           FUNCTIONS[,-train][,index],
           t(RF.Forecast)[,index])
matplot(M, col=c("blue","black","red"), type = "l", lty=c(2,1,3), lwd=c(4,3,4), xlab="Time (days)", ylab="Oil Rate (stb/day)")
legend("topright", c("True", "Random Forest", "RF + Trace Kriging"), col=c("black","blue","red"), lty=c(1,2,3), lwd=c(4,3,4))
}



# Error computation:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Sensitivity TEST::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
varSensitivity(singleTestDeuclidean)
varSensitivity(singleTestSSE)
varSensitivity(singleTestL2)
varSensitivity(singleTestMahalanobis)
varSensitivity(bagTest)
varSensitivity(singleTestRDS)

# Comparison with DGSA::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(DGSA)
kmClust <- kmeans(t(FUNCTIONS),2)
matplot(FUNCTIONS, type = "l", col=kmClust$cluster)
dgsaSensitivity <- dgsa(kmClust$cluster, COVARIATES, .interactions = TRUE)
plotParetoDGSA(dgsaSensitivity)
plotMatrixDGSA(dgsaSensitivity, .hypothesis = FALSE, order = "hclust")

# variogram:
residuals <- predictFtree(singleTestSSE, .Xnew = as.data.frame(COVARIATES))
residuals <- lapply(residuals, function(x) Reduce(cbind, x))
residuals <- FUNCTIONS - Reduce(cbind, residuals)

finterObj <- fstat(NULL, "residuals", as.data.frame(Reservoir1_allCovariates[,5:6]), as.data.frame(residuals), 1)

finterObj <- fvariogram("~.", finterObj, directions = "omni",
                        Nlags=100, LagMax = 20000, ArgStep = 1)

finterObj <- fitVariograms(finterObj,
                           model = vgm(8e+6, 'Sph', 15000, 0),
                           fitSills = TRUE, forceNugget = FALSE, fitRanges = TRUE)

plotVariogram(finterObj)

# Splitting:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# L2norm, Mahalanobis, SSE
# split: 70-30 then monte carlo.
set.seed(1)
train <- sample(nrow(COVARIATES), 0.7*nrow(COVARIATES), replace=FALSE)

# SSE:
sse.Tree_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "single", cost.type = "sse")
sse.rfor_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "randomforest", cost.type = "sse", nBoot = 1000)
sse.bagg_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "bagging", cost.type = "sse", nBoot = 1000)

sse.Tree_0.7.predict  <- predictFtree(sse.Tree_0.7, as.data.frame(COVARIATES[-train,]))
sse.rfor_0.7.predict  <- predictFtree(sse.rfor_0.7, as.data.frame(COVARIATES[-train,]))
sse.bagg_0.7.predict  <- predictFtree(sse.bagg_0.7, as.data.frame(COVARIATES[-train,]))

# Mahalanobis:
mahalanobis.Tree_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "single", cost.type = "mahalanobis")
mahalanobis.rfor_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "randomforest", cost.type = "mahalanobis", nBoot = 1000)
mahalanobis.bagg_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "bagging", cost.type = "mahalanobis", nBoot = 1000)

mahalanobis.Tree_0.7.predict  <- predictFtree(mahalanobis.Tree_0.7, as.data.frame(COVARIATES[-train,]))
mahalanobis.rfor_0.7.predict  <- predictFtree(mahalanobis.rfor_0.7, as.data.frame(COVARIATES[-train,]))
mahalanobis.bagg_0.7.predict  <- predictFtree(mahalanobis.bagg_0.7, as.data.frame(COVARIATES[-train,]))

# L2Norm:
l2norm.Tree_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "single", cost.type = "l2norm")
l2norm.rfor_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "randomforest", cost.type = "l2norm", nBoot = 1000)
l2norm.bagg_0.7 <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train], tree.type = "bagging", cost.type = "l2norm", nBoot = 1000)

l2norm.Tree_0.7.predict  <- predictFtree(l2norm.Tree_0.7, as.data.frame(COVARIATES[-train,]))
l2norm.rfor_0.7.predict  <- predictFtree(l2norm.rfor_0.7, as.data.frame(COVARIATES[-train,]))
l2norm.bagg_0.7.predict  <- predictFtree(l2norm.bagg_0.7, as.data.frame(COVARIATES[-train,]))


for(index in 1:10){
par(mfrow=c(1,2))

matplot(Reduce(cbind, sse.bagg_0.7.predict[[index]]), col = "red", lwd = 1, lty=1, type="l", ylim = c(0,500),
        ylab="Oil Rate stb/day")
lines(FUNCTIONS[,-train][,index], col="black", lty = 1, lwd = 5, type = "l", ylim = c(0,500))
lines(sse.Tree_0.7.predict[[index]][[1]], col="blue", lty=2, lwd=5)
legend("topright", c("True", "Bagging", "single Tree"), col=c("black","red", "blue"), lty=c(1,1,2), lwd=c(5,3,5))

matplot(Reduce(cbind, sse.rfor_0.7.predict[[index]]), col = "red", lwd = 1, lty=1, type="l", ylim = c(0,500),
        ylab="Oil Rate stb/day")
lines(FUNCTIONS[,-train][,index], col="black", lty = 1, lwd = 5, type = "l", ylim = c(0,500))
lines(sse.Tree_0.7.predict[[index]][[1]], col="blue", lty=2, lwd=5)
legend("topright", c("True", "R.Forest", "single Tree"), col=c("black","red", "blue"), lty=c(1,1,2), lwd=c(5,3,5))
}

sse.Tree_0.7.residuals <- predictFtree(sse.Tree_0.7, .Xnew = as.data.frame(COVARIATES[train,]))
sse.Tree_0.7.residuals <- lapply(sse.Tree_0.7.residuals, function(x) Reduce(cbind, x))
sse.Tree_0.7.residuals <- FUNCTIONS[,train] - Reduce(cbind, sse.Tree_0.7.residuals)

sse.Tree_0.7.finterObj <- fstat(NULL, "residuals", as.data.frame(COVARIATES[train,23:24]), as.data.frame(sse.Tree_0.7.residuals), 1)
sse.Tree_0.7.finterObj <- estimateDrift("~.",sse.Tree_0.7.finterObj)
sse.Tree_0.7.finterObj <- fvariogram("~.", sse.Tree_0.7.finterObj, directions = "omni",
                                    Nlags=100, LagMax = 1, ArgStep = 1)
sse.Tree_0.7.finterObj <- fitVariograms(sse.Tree_0.7.finterObj,
                                       model = vgm(8e+6, 'Sph', 0.25, 0),
                                       fitSills = TRUE, forceNugget = FALSE, fitRanges = TRUE)
plotVariogram(sse.Tree_0.7.finterObj, npSize = TRUE)
sse.Tree_0.7.finterObj <- addCovariance(sse.Tree_0.7.finterObj, 'omni')
sse.Tree_0.7.finterObj.predict <- predictFstat(sse.Tree_0.7.finterObj,
                                               .newCoordinates = as.data.frame(COVARIATES[-train,23:24]),
                                               .what = 'residuals', .type = "UK")

trend <- Reduce(cbind,lapply(sse.Tree_0.7.predict, function(x) Reduce(cbind, x)))
forecast <- trend + sse.Tree_0.7.finterObj.predict$Forecast

for(index in 1:30){
  plot(trend[,index], col="red", type = "l", lty=1, lwd=2, ylim=c(0,500))
  lines(FUNCTIONS[,-train][,index], col="black", lty = 1, lwd = 5, type = "l", ylim = c(0,500))
  lines(forecast[,index], col="blue", lty=2, lwd=5)
  legend("topright", c("True", "Tree", "Hybrid"), col=c("black","red", "blue"), lty=c(1,1,2), lwd=c(5,3,5))
}

# What happens when you load all covariates (57 of them)
st.57covariates <- ftree(.X = Reservoir1_allCovariates[,-c(1,2,22,31,34,40,49,55,56,57)],
                         .Y = FUNCTIONS, tree.type = "single", cost.type = "sse")

plotFtree(st.57covariates, .horizontal = TRUE, .labSize = 2.4)
plotBoth(st.57covariates, "0")
plotBoth(st.57covariates, "0l")
plotBoth(st.57covariates, "0lr")

varSensitivity(st.57covariates)

library(PerformanceAnalytics)
chart.Correlation(Reservoir1_allCovariates[,15:21])

# Vqtz and other petrophysical parameters can be of importance for analyses of
# type: why is well X underperforming etc.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Monte Carlo Analysis
################################################################################

PATHS <- matrix(NA, ncol = 172, nrow = 100)
for(i in 1:100){
  set.seed(i)
  PATHS[i,] <- sample(172, 172, replace = FALSE)
}

COVARIATES <- Reservoir1_allCovariates[,-c(1,2,22,31,34,40,49,55,56,57)]
FUNCTIONS  <- as.matrix(completeSet)

ERRORS.single.mean   <- matrix(NA, ncol = 10, nrow = 100)
ERRORS.single.median <- matrix(NA, ncol = 10, nrow = 100)

ERRORS.rf.mean     <- matrix(NA, ncol = 10, nrow = 100)
ERRORS.rf.median   <- matrix(NA, ncol = 10, nrow = 100)

ERRORS.bag.mean    <- matrix(NA, ncol = 10, nrow = 100)
ERRORS.bag.median  <- matrix(NA, ncol = 10, nrow = 100)

totalSSE <- mean(apply((FUNCTIONS - rowMeans(FUNCTIONS))^2, 2, function(x) trapzc(1,x)))

ErrorList <- vector("list", 10)

library(foreach)
require(doParallel)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 2:10){

  ErrorList[[i]] <- foreach(j=1:100, .packages = c("fTree", "plyr","functInterp")) %dopar% {

    train <- PATHS[j, 1:(i*10)]
    test  <- PATHS[j, 101:172]

    bs.SingleTree   <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train],
                             tree.type = "single", cost.type = "sse", .minSplit = 3, .minBucket = 1)
    bs.RandomForest <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train],
                             tree.type = "randomforest", cost.type = "sse", .minSplit = 3, .minBucket = 1,
                             nBoot = 1000, parallel = FALSE)
    bs.Bagging      <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train],
                             tree.type = "bagging", cost.type = "sse", .minSplit = 3, .minBucket = 1,
                             nBoot = 1000, parallel = FALSE)

    bs.SingleTree.predict    <- Reduce(cbind,lapply(predictFtree(bs.SingleTree,   as.data.frame(COVARIATES[test,])), function(x) t(x)))
    bs.RandomForest.predict  <- as.matrix(ldply(predictFtree(bs.RandomForest, as.data.frame(COVARIATES[test,])),
                                      function(x) colMeans(x, na.rm = TRUE), .id=NULL))
    bs.Bagging.predict       <- as.matrix(ldply(predictFtree(bs.Bagging,      as.data.frame(COVARIATES[test,])),
                                      function(x) {colMeans(x, na.rm = TRUE)}, .id = NULL))

    class(bs.RandomForest.predict) <- "numeric"
    class(bs.Bagging.predict)      <- "numeric"

    # Error quantification
    bs.SingleTree.Errors     <- apply((bs.SingleTree.predict   - t(FUNCTIONS[,test]))^2, 1, function(x) trapzc(1,x)) / totalSSE
    bs.RandomForest.Errors   <- apply((bs.RandomForest.predict - t(FUNCTIONS[,test]))^2, 1, function(x) trapzc(1,x)) / totalSSE
    bs.Bagging.predictErrors <- apply((bs.Bagging.predict      - t(FUNCTIONS[,test]))^2, 1, function(x) trapzc(1,x)) / totalSSE

    index = 20
    matplot(FUNCTIONS, type = "l", col="blue")
    matplot(FUNCTIONS[,train], col="red", type = "l", add=TRUE)
    matplot(t(bs.RandomForest.predict), col="green", type = "l", add=TRUE)
    lines(bs.RandomForest.predict[index,], col="red", lwd=5)
    lines(FUNCTIONS[,test[index]], col="black", lty = 1, lwd = 5)
    lines(bs.SingleTree.predict[,index], col="green", lty = 2, lwd = 5)

    legend("topright", c("True", "Forecast"), col=c("black","red"), lty=c(1,2), lwd=c(5,3))


    return(list(ERRORS.single.mean   = mean(bs.SingleTree.Errors),
                ERRORS.single.median = median(bs.SingleTree.Errors),
                ERRORS.rf.mean       = mean(bs.RandomForest.Errors),
                ERRORS.rf.median     = median(bs.RandomForest.Errors),
                ERRORS.bag.mean      = mean(bs.Bagging.predictErrors),
                ERRORS.bag.median    = median(bs.Bagging.predictErrors)))

  }

}

stopCluster(cl)

# Regression based approach:
cl <- makeCluster(8)
registerDoParallel(cl)

ErrorListRegression <- vector("list", 10)

for(i in 2:10){

  ErrorListRegression[[i]] <- foreach(j=1:100, .packages = c("fTree", "plyr","functInterp")) %dopar% {

    train <- PATHS[j, 1:(i*10)]
    test  <- PATHS[j, 101:172]

    # Train Model:


    # Predict:


    # Error Computation:
    bs.Regression.Errors     <- apply((bs.Regression.predict   - t(FUNCTIONS[,test]))^2, 1, function(x) trapzc(1,x)) / totalSSE

    return(list(ERRORS.regression.mean   = mean(bs.Regression.Errors),
                ERRORS.regression.median = median(bs.Regression.Errors)))

  }
}

stopCluster(cl)


# FPCA:
library(fda)
library(MuFiCokriging)

ErrorListCokriging <- vector("list", 10)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 5:10){

  ErrorListCokriging[[i]] <- foreach(j=1:100,
                                     .packages = c("fTree", "plyr","functInterp","MuFiCokriging", "fda"),
                                     .errorhandling = "pass") %dopar% {

    train <- PATHS[j, 1:(i*10)]
    test  <- PATHS[j, 101:172]

    basis  <- create.bspline.basis(c(1, 600), nbasis = 30);
    production.fd      <- smooth.basis(1:600, as.matrix(FUNCTIONS[,train]), basis)$fd
    production.fd.pca  <- pca.fd(production.fd,   nharm = 2);

    UcoKmle.K2 <- MuFicokm(formula    = list(~1, ~1),
                           MuFidesign = NestedDesign(as.matrix(COVARIATES[train,]), nlevel = 2, indices = list(1:length(train))),
                           response   = list(production.fd.pca$scores[,1], production.fd.pca$scores[,2]),
                           nlevel     = 2)

    UcoKmle.K2.predictions <- predict(UcoKmle.K2, newdata = COVARIATES[test,], "UK")

    MEAN <- eval.fd(1:600, production.fd.pca$meanfd)
    fPCS <- eval.fd(1:600, production.fd.pca$harmonics)
    SCORES <- Reduce(rbind, UcoKmle.K2.predictions$mux)
    RESIDUALS <- fPCS %*% SCORES
    FINAL <- apply(RESIDUALS,2,function(x)  x+MEAN)

    ERROR <- apply((FINAL - FUNCTIONS[,test]) ^ 2, 2, function(x) trapzc(1,x)) / totalSSE

    return(list(ERRORS.cokrig.mean   = mean(ERROR),
                ERRORS.cokrig.median = median(ERROR)))

  }
}

stopCluster(cl)


# Universal CoKriging:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ErrorListUniCokriging <- vector("list", 10)

cl <- makeCluster(8)
registerDoParallel(cl)

ucokCOVARIATES <- cbind(covariatesCompletion, covariatesGeological)

for(i in 5:10){

  ErrorListUniCokriging[[i]] <- foreach(j=1:100,
   .packages = c("fTree", "plyr","functInterp","MuFiCokriging", "fda"),
   .errorhandling = "pass") %dopar% {

     train <- PATHS[j, 1:(i*10)]
     test  <- PATHS[j, 101:172]

     basis  <- create.bspline.basis(c(1, 600), nbasis = 30);
     production.fd      <- smooth.basis(1:600, as.matrix(FUNCTIONS[,train]), basis)$fd
     production.fd.pca  <- pca.fd(production.fd,   nharm = 2);

     UcoKmle.K2 <- MuFicokm(formula    = list(~1, ~1),
                            MuFidesign = NestedDesign(as.matrix(ucokCOVARIATES[train,]), nlevel = 2, indices = list(1:length(train))),
                            response   = list(production.fd.pca$scores[,1], production.fd.pca$scores[,2]),
                            nlevel     = 2)

     UcoKmle.K2.predictions <- predict(UcoKmle.K2, newdata = ucokCOVARIATES[test,], "UK")

     MEAN <- eval.fd(1:600, production.fd.pca$meanfd)
     fPCS <- eval.fd(1:600, production.fd.pca$harmonics)
     SCORES <- Reduce(rbind, UcoKmle.K2.predictions$mux)
     RESIDUALS <- fPCS %*% SCORES
     FINAL <- apply(RESIDUALS,2,function(x)  x+MEAN)

     ERROR <- apply((FINAL - FUNCTIONS[,test]) ^ 2, 2, function(x) trapzc(1,x)) / totalSSE

     return(list(ERRORS.cokrig.mean   = mean(ERROR),
                 ERRORS.cokrig.median = median(ERROR)))

   }
}

stopCluster(cl)

# HYBRID MODEL::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ErrorListHybrid <- vector("list", 10)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 5:10){

  ErrorListHybrid[[i]] <- foreach(j=1:100,
     .packages = c("fTree", "plyr","functInterp","MuFiCokriging", "fda"),
     .errorhandling = "pass") %dopar% {

       train <- PATHS[j, 1:(i*10)]
       test  <- PATHS[j, 101:172]

       # Random Forest component::::::::::::::::::::::::::::::::::::::::::::::::
       system.time(bs.RandomForest <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train],
                                tree.type = "randomforest", cost.type = "mahalanobis", .minSplit = 10,
                                nBoot = 100, parallel = FALSE, .cp = 0.001))
       bs.RandomForest.predict.TRAIN  <- as.matrix(ldply(predictFtree(bs.RandomForest, as.data.frame(COVARIATES[train,])),
                                                   function(x) colMeans(x, na.rm = TRUE), .id=NULL))
       bs.RandomForest.predict.TEST   <- as.matrix(ldply(predictFtree(bs.RandomForest, as.data.frame(COVARIATES[test,])),
                                                         function(x) colMeans(x, na.rm = TRUE), .id=NULL))
      ##########################################################################

       # plotBoth(bs.RandomForest, "0", 500)

       # matplot(FUNCTIONS[,train], type = "l", col="blue")
       # matplot(t(bs.RandomForest.predict.TRAIN), type = "l", col="red", add=TRUE)

       # matplot(FUNCTIONS[,test], type = "l", col="blue")
       # matplot(t(bs.RandomForest.predict.TEST), type = "l", col="red", add=TRUE)

       # A <- predictFtree(bs.RandomForest, as.data.frame(COVARIATES[train[1:2],]))
       # matplot(FUNCTIONS, type = "l", col="blue")
       # matplot(t(A[[1]]), type="l", col="red", add=TRUE)

       RESIDUALS <- FUNCTIONS[,train] - t(bs.RandomForest.predict.TRAIN)

       basis  <- create.bspline.basis(c(1, 600), nbasis = 30);
       residual.fd      <- smooth.basis(1:600, as.matrix(RESIDUALS), basis)$fd
       residual.fd.pca  <- pca.fd(residual.fd,   nharm = 2);

       UcoKmle.K2 <- MuFicokm(formula    = list(~1, ~1),
                              MuFidesign = NestedDesign(as.matrix(COVARIATES[train,]),
                                                        nlevel = 2, indices = list(1:length(train))),
                              response   = list(residual.fd.pca$scores[,1], residual.fd.pca$scores[,2]),
                              nlevel     = 2)

       UcoKmle.K2.predictions <- predict(UcoKmle.K2, newdata = COVARIATES[test,], "SK")

       MEAN <- eval.fd(1:600, residual.fd.pca$meanfd)
       fPCS <- eval.fd(1:600, residual.fd.pca$harmonics)
       SCORES <- Reduce(rbind, UcoKmle.K2.predictions$mux)
       RESIDUALS <- fPCS %*% SCORES
       FINAL <- apply(RESIDUALS,2,function(x)  x+MEAN)

       # matplot(FUNCTIONS[,test], type = "l", col="blue")
       # matplot(FINAL + t(bs.RandomForest.predict.TEST), type = "l", col="red", add=TRUE)

       ERROR <- apply((FINAL + t(bs.RandomForest.predict.TEST) - FUNCTIONS[,test]) ^ 2, 2, function(x) trapzc(1,x)) / totalSSE

       return(list(ERRORS.hybrid.mean   = mean(ERROR),
                   ERRORS.hybrid.median = median(ERROR)))

     }
}

stopCluster(cl)


# Spatial Hybrid::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ErrorListHybridSpatial <- vector("list", 10)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 2:10){

  ErrorListHybridSpatial[[i]] <- foreach(j=1:100,
    .packages = c("fTree", "plyr","functInterp","MuFiCokriging", "fda"),
    .errorhandling = "pass") %dopar% {

      train <- PATHS[j, 1:(i*10)]
      test  <- PATHS[j, 101:172]

      # Random Forest component::::::::::::::::::::::::::::::::::::::::::::::::
      system.time(bs.RandomForest <- ftree(.X = COVARIATES[train,], .Y = FUNCTIONS[,train],
                                           tree.type = "randomforest", cost.type = "sse", .minSplit = 3, .minBucket = 1,
                                           nBoot = 500, parallel = FALSE))
      bs.RandomForest.predict.TRAIN  <- as.matrix(ldply(predictFtree(bs.RandomForest, as.data.frame(COVARIATES[train,])),
                                                        function(x) colMeans(x, na.rm = TRUE), .id=NULL))
      bs.RandomForest.predict.TEST   <- as.matrix(ldply(predictFtree(bs.RandomForest, as.data.frame(COVARIATES[test,])),
                                                        function(x) colMeans(x, na.rm = TRUE), .id=NULL))
      ##########################################################################

      RESIDUALS <- FUNCTIONS[,train] - t(bs.RandomForest.predict.TRAIN)

      basis  <- create.bspline.basis(c(1, 600), nbasis = 30);
      residual.fd      <- smooth.basis(1:600, as.matrix(RESIDUALS), basis)$fd
      residual.fd.pca  <- pca.fd(residual.fd,   nharm = 2);

      UcoKmle.K2 <- MuFicokm(formula    = list(~1, ~1),
                             MuFidesign = NestedDesign(as.matrix(COVARIATES[train,3:4]),
                                                       nlevel = 2, indices = list(1:length(train))),
                             response   = list(residual.fd.pca$scores[,1], residual.fd.pca$scores[,2]),
                             nlevel     = 2)

      UcoKmle.K2.predictions <- predict(UcoKmle.K2, newdata = COVARIATES[test,3:4], "SK")

      MEAN <- eval.fd(1:600, residual.fd.pca$meanfd)
      fPCS <- eval.fd(1:600, residual.fd.pca$harmonics)
      SCORES <- Reduce(rbind, UcoKmle.K2.predictions$mux)
      RESIDUALS <- fPCS %*% SCORES
      FINAL <- apply(RESIDUALS,2,function(x)  x+MEAN)

      ERROR <- apply((FINAL + t(bs.RandomForest.predict.TEST) - FUNCTIONS[,test]) ^ 2, 2, function(x) trapzc(1,x)) / totalSSE

      return(list(ERRORS.spatialhybrid.mean   = mean(ERROR),
                  ERRORS.spatialhybrid.median = median(ERROR)))

    }
}

stopCluster(cl)


# Regression::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
source("../../Summer_Project_2016/Code/Lijing Implemented/basisRegress.R")

fReg <- fstat(NULL, 'functions', as.data.frame(ucokCOVARIATES), as.data.frame(FUNCTIONS))
fReg <- basisRegress(~., fReg, CB = FALSE)

ErrorListRegress <- vector("list", 10)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 1:10){

  ErrorListRegress[[i]] <- foreach(j=1:100,
   .packages = c("fTree", "plyr","functInterp","MuFiCokriging", "fda"),
   .errorhandling = "pass") %dopar% {

     train <- PATHS[j, 1:(i*10)]
     test  <- PATHS[j, 101:172]

     fReg <- fstat(NULL, 'functions', as.data.frame(ucokCOVARIATES[train,]), as.data.frame(FUNCTIONS[,train]))
     fReg <- basisRegress(~., fReg, CB = FALSE)
     FINAL <- predBasisRegress(fReg, .newCoordinates = as.data.frame(ucokCOVARIATES[test,]), CB = FALSE)

     ERROR <- apply((t(FINAL$data$functions$predYhat) - FUNCTIONS[,test]) ^ 2, 2, function(x) trapzc(1,x)) / totalSSE

     return(list(ERRORS.spatialhybrid.mean   = mean(ERROR),
                 ERRORS.spatialhybrid.median = median(ERROR)))

   }
}

stopCluster(cl)


# Hybrid::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ErrorListRegressCOK <- vector("list", 10)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 1:10){

  ErrorListRegressCOK[[i]] <- foreach(j=1:100,
       .packages = c("fTree", "plyr","functInterp","MuFiCokriging", "fda"),
       .errorhandling = "pass") %dopar% {

         train <- PATHS[j, 1:(i*10)]
         test  <- PATHS[j, 101:172]

         fReg <- fstat(NULL, 'functions', as.data.frame(ucokCOVARIATES[train,]), as.data.frame(FUNCTIONS[,train]))
         fReg <- basisRegress(~., fReg, CB = FALSE)
         FINALmean <- predBasisRegress(fReg, .newCoordinates = as.data.frame(ucokCOVARIATES[test,]), CB = FALSE)

         RESIDUALS <- FUNCTIONS[,train] - t(FINALmean$data$functions$yhat)

         basis  <- create.bspline.basis(c(1, 600), nbasis = 30);
         residual.fd      <- smooth.basis(1:600, as.matrix(RESIDUALS), basis)$fd
         residual.fd.pca  <- pca.fd(residual.fd,   nharm = 2);

         UcoKmle.K2 <- MuFicokm(formula    = list(~1, ~1),
                                MuFidesign = NestedDesign(as.matrix(ucokCOVARIATES[train,]),
                                                          nlevel = 2, indices = list(1:length(train))),
                                response   = list(residual.fd.pca$scores[,1], residual.fd.pca$scores[,2]),
                                nlevel     = 2)

         UcoKmle.K2.predictions <- predict(UcoKmle.K2, newdata = ucokCOVARIATES[test,], "SK")

         MEAN <- eval.fd(1:600, residual.fd.pca$meanfd)
         fPCS <- eval.fd(1:600, residual.fd.pca$harmonics)
         SCORES <- Reduce(rbind, UcoKmle.K2.predictions$mux)
         RESIDUALS <- fPCS %*% SCORES
         FINALresid <- apply(RESIDUALS,2,function(x)  x+MEAN)

         # matplot(FUNCTIONS[,test], type = "l", col="blue")
         # matplot(FINALresid + t(FINALmean$data$functions$predYhat), type = "l", col="red", add=TRUE)
         # matplot(FINALresid, type = "l", col="blue")
         #
         # ERRORl <- apply((t(FINALmean$data$functions$predYhat) - FUNCTIONS[,test]) ^ 2, 2, function(x) trapzc(1,x)) / totalSSE

         .index = 12
         plot(FUNCTIONS[,test[.index]], col="blue", main = ERROR[.index])
         lines(t(FINALmean$data$functions$predYhat)[,.index], col="red")
         lines((FINALresid + t(FINALmean$data$functions$predYhat))[,.index], col="black")

         ERROR <- apply((FINALresid + t(FINALmean$data$functions$predYhat) - FUNCTIONS[,test]) ^ 2, 2, function(x) trapzc(1,x)) / totalSSE

         return(list(ERRORS.spatialhybrid2.mean   = mean(ERROR),
                     ERRORS.spatialhybrid2.median = median(ERROR)))

       }
}

stopCluster(cl)





# Plotting:
for(i in 5:10) {
  for(j in 1:100){
    ErrorList[[i]][[j]] <- as.data.frame(ErrorList[[i]][[j]])
  }
}

# Fix for cokriging:
for(i in 5:10) {
  for(j in 1:length(ErrorListCokriging[[i]])){
    if(is.numeric(ErrorListCokriging[[i]][[j]][[1]])){
      ErrorListCokriging[[i]][[j]] <- as.data.frame(ErrorListCokriging[[i]][[j]])
    } else {
      ErrorListCokriging[[i]][[j]] <- data.frame(ERRORS.cokrig.mean = NaN,
                                                 ERRORS.cokrig.median = NaN)
    }
  }
}

# Fix for hybrid:
for(i in 5:10) {
  for(j in 1:length(ErrorListHybrid[[i]])){
    if(is.numeric(ErrorListHybrid[[i]][[j]][[1]])){
      ErrorListHybrid[[i]][[j]] <- as.data.frame(ErrorListHybrid[[i]][[j]])
    } else {
      ErrorListHybrid[[i]][[j]] <- data.frame(ERRORS.hybrid.mean = NaN,
                                                 ERRORS.hybrid.median = NaN)
    }
  }
}

for(i in 2:10) {
  for(j in 1:length(ErrorListHybridSpatial[[i]])){
    if(is.numeric(ErrorListHybridSpatial[[i]][[j]][[1]])){
      ErrorListHybridSpatial[[i]][[j]] <- as.data.frame(ErrorListHybridSpatial[[i]][[j]])
    } else {
      ErrorListHybridSpatial[[i]][[j]] <- data.frame(ERRORS.hybrid.mean = NaN,
                                              ERRORS.hybrid.median = NaN)
    }
  }
}

for(i in 1:10) {
  for(j in 1:length(ErrorListRegress[[i]])){
    if(is.numeric(ErrorListRegress[[i]][[j]][[1]])){
      ErrorListRegress[[i]][[j]] <- as.data.frame(ErrorListRegress[[i]][[j]])
    } else {
      ErrorListRegress[[i]][[j]] <- data.frame(ERRORS.hybrid.mean = NaN,
                                                     ERRORS.hybrid.median = NaN)
    }
  }
}

for(i in 1:10) {
  for(j in 1:length(ErrorListRegressCOK[[i]])){
    if(is.numeric(ErrorListRegressCOK[[i]][[j]][[1]])){
      ErrorListRegressCOK[[i]][[j]] <- as.data.frame(ErrorListRegressCOK[[i]][[j]])
    } else {
      ErrorListRegressCOK[[i]][[j]] <- data.frame(ERRORS.spatialhybrid2.mean = NaN,
                                               ERRORS.spatialhybrid2.median = NaN)
    }
  }
}

for(i in 1:10) {
  # temp <- melt(Reduce(rbind, ErrorList[[i]]))
  # temp$bin <- i*10
  # ErrorList[[i]] <- temp

  # temp <- melt(Reduce(rbind, ErrorListCokriging[[i]]))
  # temp$bin <- i*10
  # ErrorListCokriging[[i]] <- temp

  # temp <- melt(Reduce(rbind, ErrorListHybrid[[i]]))
  # temp$bin <- i*10
  # ErrorListHybrid[[i]] <- temp

  # temp <- melt(Reduce(rbind, ErrorListHybridSpatial[[i]]))
  # temp$bin <- i*10
  # ErrorListHybridSpatial[[i]] <- temp

  # temp <- melt(Reduce(rbind, ErrorListRegress[[i]]))
  # temp$bin <- i*10
  # ErrorListRegress[[i]] <- temp

  temp <- melt(Reduce(rbind, ErrorListRegressCOK[[i]]))
  temp$bin <- i*10
  ErrorListRegressCOK[[i]] <- temp

}



ggDATA <- rbind(Reduce(rbind, ErrorList), Reduce(rbind, ErrorListCokriging),
                Reduce(rbind, ErrorListHybrid), Reduce(rbind, ErrorListRegress), Reduce(rbind, ErrorListRegressCOK))

ggDATA$method <- ggDATA$variable
ggDATA$method <- gsub("ERRORS.", "", ggDATA$method)
ggDATA$type   <- ggDATA$method

ggDATA$method <- gsub("rf.mean", "Random Forest", ggDATA$method)
ggDATA$method <- gsub("rf.median", "Random Forest", ggDATA$method)
ggDATA$method <- gsub("bag.mean", "Bagging", ggDATA$method)
ggDATA$method <- gsub("bag.median", "Bagging", ggDATA$method)
ggDATA$method <- gsub("single.mean", "Single Tree", ggDATA$method)
ggDATA$method <- gsub("single.median", "Single Tree", ggDATA$method)
ggDATA$method <- gsub("cokrig.mean", "Cokriging", ggDATA$method)
ggDATA$method <- gsub("cokrig.median", "Cokriging", ggDATA$method)
ggDATA$method <- gsub("spatialhybrid.mean", "Regression", ggDATA$method)
ggDATA$method <- gsub("spatialhybrid.median", "Regression", ggDATA$method)
ggDATA$method <- gsub("spatialhybrid2.mean", "Reg CoK Hybrid (XY)", ggDATA$method)
ggDATA$method <- gsub("spatialhybrid2.median", "Reg CoK Hybrid (XY)", ggDATA$method)
ggDATA$method <- gsub("hybrid.mean", "RF CoK Hybrid (all)", ggDATA$method)
ggDATA$method <- gsub("hybrid.median", "RF CoK Hybrid (all)", ggDATA$method)


ggDATA$type   <- gsub("single.", "", ggDATA$type)
ggDATA$type   <- gsub("rf.", "", ggDATA$type)
ggDATA$type   <- gsub("bag.", "", ggDATA$type)
ggDATA$type   <- gsub("cokrig.", "", ggDATA$type)
ggDATA$type   <- gsub("spatialhybrid2.", "", ggDATA$type)
ggDATA$type   <- gsub("spatialhybrid.", "", ggDATA$type)
ggDATA$type   <- gsub("hybrid.", "", ggDATA$type)


ggplot(ggDATA) + geom_boxplot(aes(x=method, y=value, colour=method), outlier.shape=NA) + facet_grid(type~bin, scales = "free_y" ) +
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + ylim(c(0,3))+
  xlab("The number of training wells") + ylab("Normalized SSE")


ggDATA <- rbind(Reduce(rbind, ErrorListRegress))
ggDATA$method <- ggDATA$variable
ggDATA$method <- gsub("ERRORS.", "", ggDATA$method)
ggDATA$type   <- ggDATA$method
ggDATA$method <- gsub("spatialhybrid.mean", "Regression", ggDATA$method)
ggDATA$method <- gsub("spatialhybrid.median", "Regression", ggDATA$method)
ggDATA$type   <- gsub("spatialhybrid.", "", ggDATA$type)
ggplot(ggDATA) + geom_boxplot(aes(x=method, y=value, colour=type), outlier.shape=NA) + facet_grid(.~bin, scales = "free_y" ) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(c(0,3))+ggtitle("Regression Performance")+
  xlab("The number of training wells") + ylab("Normalized SSE")










# Forecasting Monte Carlo:::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Questions:
# How many wells do we need in order to start predicting accuratelly?
# Is there a difference in error quantification when you split into pads and
# When you split randomly?
# How to interpret the variogram?

# Random Forest alone
# Single Tree alone
# Bagging alone
# Bagging combined with functional kriging
# Random forest combined with functional kriging
# Ultimately cokriging with well logs as well.


# What needs to be done:
# Compare all distance measurements along with the plots of separation. DONE
# Compare sensitivities achieved with each cost function and tree type
# Conduct a forecasting study with the following parameters:
# Random splitting all wells.
# Random splitting into pads.
# Report SSE & containment within confidence bands.

# It is enough to report that the method works on the same domain of influence. It is still practical
# It is also enough to report SSE error.

# It is highly possible that neural network with two output parameters (scores) would outperform
# This entire machine learning approach..... it is to be seen. we don't have any categorical parameters.



################################################################################
# Update tests: and verify DONE
# Test tree split: DONE
# Test visualizations; upgrade to have height as input parameter. DONE
# Write R wrapper; this will check wheher we are doing random forest or bagging or simply one tree.
# Write Interaction function;
# Read about sensitivity in random forest;
# Parallel processing for random forest and bagging. DONE
# Visualization function(s) DONE

# Visualization function. It needs to view the node, the tree, etc:

# Pruning function:

# Variable importance is outside all this. Since we will be implementing several methods.
# The most basic one will just summarize the tree.
# The second one will run permutation.
# The third one will run growth sequence.

# Advantages of the tree based methods:
# A. Works with any type of input parameter (number, category, function)
# B. It works with functional outputs and images/geomodels....
# C. It provides a nice sensitivity/importance summary
# D. It also provides emulation capabilities. Namely we are able to run
#    sensitivity with proxies and full physics all together.

# Plan: See what you need to visualize the tree, then proceed to adjust the tree generating code
# Develop unit tests for everything!
# Visualizers:
#   a. Visualize the tree (ggplot graphics) simple and straightforward; DONE
#   b. Visualize the data. One visualizer for functional data (from the environment);
#   c. Visualize the data. One visualizer for Distance based data;
#   d. Visualize the data. Combined plot
# Need a function to prune the tree
# Need a function to compute all residuals from the tree
# Need a function to compute sensitivities from the tree (all methods)

# Visualizer is one function that calls subfunctions.
# Also need a data explorer from ggplot. Namely to identify the node that is closest to some point.
# One click updates it all.



# compare:
# regression
# random forest    DONE
# bagging          DONE
# single tree      DONE
# cokriging
# hybrid tree
# hybrid forest
# hybrid bagging


