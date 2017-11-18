library(fTree)
load("~/Google Drive/Functional_Research/PhD_Dissertation/Functional_Interpolation_Multivariate/CleanData.RData")

COVARIATES <- WELL_PARAMETERS[,c(4:7,16,19,20,23,24,25:27,32)]

# Trees for data mining:
tree.A <- ftree(.X = COVARIATES[,-1], .Y = OilRates, tree.type = 'single', cost.type = 'l2square')
plotFtree(tree.A, .labSize = 2.5, .horizontal = T, .ylimLow = -2, .round = 1)
varSensitivity(tree.A)

# Similarly we can fit a bagged functional tree and view its sensitivity:
tree.B <- ftree(.X = WELL_PARAMETERS[,-1], .Y = OilRates, tree.type = 'bagging', nBoot = 200, cost.type = 'sse')
plotFtree(tree.B, .labSize = 2.5, .horizontal = T, .ylimLow = -2, .round = 1, .index = 3)
varSensitivity(tree.B)

# The interface allows for different cost functions:
tree.C <- ftree(.X = WELL_PARAMETERS[,-1], .Y = OilRates, tree.type = 'bagging', nBoot = 200, cost.type = 'sse')
plotFtree(tree.C, .labSize = 2.5, .horizontal = T, .ylimLow = -2, .round = 1, .index = 1)
varSensitivity(tree.C)

# Growing trees with distances:
tree.D <- ftree(.X = WELL_PARAMETERS[,-1], .Y = OilRates, .D = "euclidean", tree.type = 'bagging', nBoot = 200, cost.type = 'rdist')
plotFtree(tree.D, .labSize = 2.5, .horizontal = T, .ylimLow = -2, .round = 1, .index = 1)
varSensitivity(tree.D)

# Growing trees with multivariate functional data: part 1




# Growing trees with multivariate functional data: part 2







