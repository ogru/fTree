#' @useDynLib fTree
#' @importFrom Rcpp sourceCpp
#'
#' @title Fits a regression tree to functional data
#'
#' @param .X - An nxp matrix of covariates (n= # of observations, p=# of covariates)
#' @param .Y - A matrix of functions stacked in columns. Assumes that all functions were evaluated on the same time grid
#' @param .D - Optional distance matrix for wss/rdist cost function
#' @param cost.type  - Cost function type. It can be any of the following: "sse", "mahalanobis", "wss", "l2norm", "rdist" see details.
#' @param tree.type  - What type of tree based predictor you want to fit. Currently supported: single tree, random forest, bagging
#' @param nP - The number of predictors to consider on each attempted split. Active only for tree.type = "randomforest"
#' @param nBoot - The number of trees to consider in bootstrapping
#' @param .minSplit - minimum required number of elements in a node in order to attempt a split.
#' @param .minBucket - minimum number of elements in leaf nodes. Defaults to .minSplit/3.
#' @param .cp - complexity parameter, split is accepted if it provides imporovement that is at least cp*rootGoodness
#' @param verbose - print progres (default = TRUE)
#'
#' @details This code implements various functional and multivariate tree splitting routines. See vignette for a detailed description of
#'          each cost function and for a tutorial on how to use this function.
#'
#' @author Ognjen Grujic (\email{ogyg@stanford.edu} or \email{ognjengr@gmail.com})
#'
#' @export
#'
ftree <- function(.X = NULL, .Y = NULL, .D = NULL, SIGMA_inv = NULL, cost.type = 'sse',
                  tree.type = "single", nP = if(tree.type == "randomforest") round((ncol(.X)/3))
                  else ncol(.X), nBoot = 1000,
                  .minSplit = 20, .minBucket = round(.minSplit/3), .cp = 0.005,
                  verbose = TRUE, parallel = TRUE) {

  require(foreach)
  require(doParallel)

  cost.type <- match.arg(cost.type, c("sse", "mahalanobis", "wss", "l2norm", "rdist"))
  tree.type <- match.arg(tree.type, c("single", "randomforest", "bagging"))

  if(is.null(.X)) stop('Covariates were not provided!')

  .X = as.matrix(.X)

  if(cost.type == "wss" | cost.type == "rdist"){

    if(is.null(.D)) {
      stop('For cost type WSS you MUST provide a distance matrix or distance type (i.e. euclidean). ')
    } else if(is.character(.D)){
      if(is.null(.Y)) {
        stop('Functions were not provided! Please provide functions and distance type or just a distance matrix.')
      } else {
        .D <- as.matrix(dist(as.matrix(t(.Y)), method = .D))
        SIGMA_inv = matrix(NA, 2, 2)
      }
    } else if(is.matrix(.D)){
      SIGMA_inv = matrix(NA,2,2)
    } else {
      stop('.D is not a matrix nor a string. Exiting!')
    }

    if(ncol(.D) != nrow(.X)) stop("The number of observations does not match the dimensions of the distance matrix.")

    .nObservations = ncol(.D)

  } else { # sse, mahalanobis, l2norm

    .Y = as.matrix(.Y)
    .nObservations = ncol(.Y)

    if(.nObservations != nrow(.X)) stop("The number of observations does not match the number of functions. Exiting!")

    if(cost.type == "mahalanobis"){
      if(is.null(SIGMA_inv)){
        SIGMA_inv <- solve(as.matrix(nearPD(cov(t(.Y)))$mat))
        .D = matrix(NA, 2, 2)
      } else {
        if((ncol(.Y) != ncol(SIGMA_inv)) && (ncol(.Y) != nrow(.Y))) {
          stop('nrow/ncol(SIGMA_inv) is not the same as ncol(.Y). Exiting!')
        }
      }
    } else {
      SIGMA_inv = matrix(NA, 2, 2)
      .D <- SIGMA_inv
    }
  }

  if(tree.type != "single") { # bagging, random forest....
    if(nP > ncol(.X)) {
      stop("nP cannot be larger than the number of covariates in .X.
           Random forest gives the best result for: nP = sqrt(ncol(.X))")
    }
    } else {
      nP = ncol(.X)
  }

  out <- list()
  out[['covariates']] <- .X
  out[['functions']]  <- .Y
  out[['costType']]   <- cost.type
  out[['treeType']]   <- tree.type
  out[['minSplit']]   <- .minSplit
  out[['minBucket']]  <- .minBucket
  out[['cP']]         <- .cp

  .cType = switch(cost.type,
                 "sse"         = 1,
                 "mahalanobis" = 2,
                 "l2norm"      = 3,
                 "wss"         = 4,
                 "rdist"       = 5)

  if(tree.type=="single"){

    loadVars(.X, .Y, SIGMA_inv, .D, nP, .cType, .minSplit, .minBucket, .cp)

    .INDEX = seq(1, nrow(.X) - 1)
    .cGood         <- fTree:::computeGoodness(.INDEX)
    out[['trees']][[1]] <- fTree:::fTreeRPart(.INDEX, .cGood, 0)

    unloadVars();

  } else {

    .BOOTINDEX <- matrix(sample(nrow(.X), nrow(.X) * nBoot, replace=TRUE),
                         ncol = nrow(.X), nrow=nBoot)

    if(parallel){

      .no_cores <- detectCores() - 1
      cl <- makeCluster(.no_cores)
      registerDoParallel(cl)

      .bootList   <- matToList(.BOOTINDEX, .no_cores)

      out[['trees']] <- Reduce("c", foreach(i = 1:length(.bootList)) %dopar% {
        loadVars(.X, .Y, SIGMA_inv, .D, nP, .cType, .minSplit, .minBucket, .cp)
        output <- fTreeBootstrap(.bootList[[i]] - 1)
        unloadVars();
        return(output)
      })

      stopCluster(cl)

    } else{ # single core mode.

      loadVars(.X, .Y, SIGMA_inv, .D, nP, .cType, .minSplit, .minBucket, .cp)
      out[['trees']] <- fTree:::fTreeBootstrap(.BOOTINDEX - 1)
      unloadVars();

    }
  }

  # Save the structure:
  class(out) <- c('fTree', "list")
  return(out);

}

#' @title A helper function for parallel computing (not exported)
#'
#' @param .matrix - A matrix to split into list;
#' @param .chunksize - A number of rows in each list element;
#'
#' @details This code simply splits rows of a matrix into a list where each element of the list
#'          contains .chunksize number of rows of the input matrix. If .chunksize > nrow(.matrix) it
#'          returns a list with one element, the .matrix itself.
#'
#' @author Ognjen Grujic (\email{ognjengr@gmail.com})
#'
matToList <- function(.matrix, .chunksize){

  if(.chunksize >= nrow(.matrix)){
    outList = list(indices=.matrix)
  } else {
    outList <- vector('list', ceiling(nrow(.matrix)/.chunksize))
    for(i in 1:length(outList)){
      if(nrow(.matrix) < .chunksize){
        outList[[i]] <- .matrix
      } else {
        outList[[i]] <- .matrix[1:.chunksize,]
        .matrix <- .matrix[-c(1:.chunksize), , drop=FALSE]
      }
    }
  }
  return(outList)

}

#' @title Extracts node data into a data frame
#'
#' @param .tree - an "fTree" object
#' @param .index - index of the tree you wish to process
#'
#' @export
extractNodeData <- function(.tree, .index, .node = NULL){

  if(is.null(.tree)) stop('fTree object was not provided. Exiting!')
  if(class(.tree)[1] != 'fTree') stop("Passed object is not of class 'fTree'. Exiting!")

  .segments = data.frame(xStart    = numeric(),
                         xEnd      = numeric(),
                         yStart    = numeric(),
                         yEnd      = numeric(),
                         predictor = character(),
                         value     = character(),
                         position  = character(),
                         left      = character(),
                         right     = character(),
                         colour    = character())

  .segments <-    data.frame(xStart   = .tree$trees[[.index]]$midpoint,
                                xEnd     = .tree$trees[[.index]]$midpoint,
                                yStart   = -0.5,
                                yEnd     = 0,
                                predictor= NA,
                                value    = NA,
                                position = "0",
                                left     = NA,
                                right    = NA,
                                colour   = "black")

  .extractNode <- function(.NODE, .nLeftLeaves, .position){

    if(.NODE$isLeaf == 0){

      # Color:
      if(!is.null(.node) && (nchar(.position)) >= nchar(.node)){

        cp = substr(.position, 1, nchar(.node))
        lr = substr(.position, nchar(.node) + 1, nchar(.node) + 1)

        if(.node == cp){

          if(nchar(.position) == nchar(.node)){
            colour = c("red","blue")
          } else {
            if(lr == "r"){
              colour = c("blue", "blue")
            } else {
              colour = c("red", "red")
            }
          }
        } else {
          colour = c("black", "black")
        }

      } else {
        colour = c("black", "black")
      }


      # Horizontal Line:
      horizontals = data.frame(xStart    = rep(.NODE$midpoint + .nLeftLeaves, 2),
                               xEnd      = c(.NODE$left$midpoint + .nLeftLeaves,
                                             .NODE$left$nLeaves   + .nLeftLeaves + .NODE$right$midpoint),
                               yStart    = rep(.NODE$depth, 2),
                               yEnd      = rep(.NODE$depth, 2),
                               predictor = c(.NODE$bestPredictor + 1, NA),
                               value     = c(.NODE$splitPoint, NA),
                               position  = rep(.position, 2),
                               left      = NA,
                               right     = NA,
                               colour    = colour,
                               stringsAsFactors = FALSE)

      horizontals$left[1]  <- paste(.NODE$left$indices, collapse=",")
      horizontals$right[1] <- paste(.NODE$right$indices, collapse=",")

      # Vertical Lines:
      verticals  = data.frame(xStart   = c(.NODE$left$midpoint + .nLeftLeaves,
                                           .NODE$left$nLeaves  + .nLeftLeaves + .NODE$right$midpoint),
                              xEnd     = c(.NODE$left$midpoint + .nLeftLeaves,
                                           .NODE$left$nLeaves  + .nLeftLeaves + .NODE$right$midpoint),
                              yStart   = rep(.NODE$depth, 2),
                              yEnd     = rep(.NODE$depth + 1, 2),
                              predictor= NA,
                              value    = NA,
                              position = c(paste(.position, "l", sep=""),
                                           paste(.position, "r", sep="")),
                              left     = NA,
                              right    = NA,
                              colour   = colour)

      .segments <<- rbind(.segments, horizontals, verticals)

      # Continue calling:
      .extractNode(.NODE$left,  .nLeftLeaves                     , paste(.position, "l", sep=""))
      .extractNode(.NODE$right, .nLeftLeaves + .NODE$left$nLeaves, paste(.position, "r", sep=""))
    }
  }

  .extractNode(.tree$trees[[.index]], 0, "0")

  .segments$label <- paste(colnames(.tree$covariates)[.segments$predictor], .segments$value, sep=" > ")
  .segments$label[is.na(.segments$predictor)] <- NA
  .segments$colour = factor(.segments$colour, levels=c("black","red","blue"))
  return(.segments)

}

#' @title Plots a functional regression tree
#'
#' @param .object - An fTree object or .segments object produced with \code{extractNodeData}
#' @param .index - Index of the tree you wish you plot
#' @param .horizontal - Whether to plot tree with horizontal orientation
#' @param .lineWidth - Width of the tree lines
#' @param .depth - Maximum plotting depth (useful for large trees)
#' @param .labSize - size of node labels
#' @param .labAng - Angle of node labels
#' @param .ggReturn - Whether to return an ggplot object or not.
#'
#' @details This function relies on "extractNodeData" function and ggplot plotting package
#'          output of this function are publication ready graphics of fitted regression trees.
#'
#' @author Ognjen Grujic (\email{ogyg@stanford.edu})
#'
#' @export
#'
plotFtree <- function(.object,
                      .index = 1,
                      .horizontal = FALSE,
                      .lineWidth  = 1,
                      .depth      = Inf,
                      .labSize    = 1,
                      .labAng     = 0,
                      .ggReturn   = FALSE,
                      .ylimLow    = 0,
                      .ylimHigh   = NA,
                      .node       = NULL
) {

  if(is.null(.object)) stop("Object was not provided! Exiting.")

  if(class(.object) == "fTree"){
      if(.index > length(.object$trees)){
        stop('Provided .index is larger than the number of trees in the structure.')
      }
     .segments <- extractNodeData(.object, .index, .node)
  } else if(class(.object) == "fTreeSegments"){
     .segments <- .object
  } else {
     stop("Provided object is not valid. Either provide an fTree object
         or segments data frame produced with 'extractNode' function. Exiting!")
  }

  ggBASE <- ggplot(data=subset(.segments, yStart < .depth)) +
            geom_segment(aes(x=xStart, y=yStart, xend=xEnd, yend=yEnd, colour = colour), lwd = .lineWidth)  +
            geom_label(aes(x=xStart, y=yStart - 0.05, label=label),
               vjust="center",
               angle = .labAng, size = .labSize) +
            scale_colour_manual(values=c("black","red","blue")) + guides(colour=FALSE)

  if(.horizontal) {
    ggBASE <- ggBASE + ylim(.ylimLow, .ylimHigh) + coord_flip()
  } else {
    ggBASE <- ggBASE + scale_y_reverse(limits=c(8, -1.5)) + xlim(0, max(.segments$xEnd) + 1)
  }

  if(.ggReturn) {
    return(ggBASE)
  } else {
    print(ggBASE)
  }

}


#' @title Plots node data of a functional regression tree
#'
#' @param ftreeObj - An fTree object
#' @param .treeIndex - Index of the tree you wish you plot
#' @param .node - A code for specific node to plot
#' @param .ggRet - Whether to return ggplot structure
#'
#' @details Depending on the type of the fitted tree this function will either plot raw functions
#'          or low dimensional representation with fPCA or MDS.
#' @author Ognjen Grujic \email{ogyg@stanford.edu}
#'
#' @export
#'
plotNodeData <- function(ftreeObj, .treeIndex = 1, .node = "0", .ggRet = FALSE){

  if(is.null(ftreeObj)) stop("Tree was not provided. Exiting")
  if(is.null(ftreeObj$functions)) stop("functions were not found in ftreeObj. Exiting")

  .getIndices <- function(.NODE, .position){

    if(nchar(.position) == 1){
      if(.NODE$isLeaf == 1){
        return(list(centerInd = .NODE$indices + 1))
      } else {
        return(list(leftInd  = .NODE$left$indices + 1,
                    rightInd = .NODE$right$indices + 1))
      }
    } else {

      if(substr(.position, 2, 2) == "r"){
        .getIndices(.NODE$right, substr(.position,2, nchar(.position)))
      } else {
        .getIndices(.NODE$left, substr(.position,2, nchar(.position)))
      }
    }
  }

  .indices <- .getIndices(ftreeObj$trees[[.treeIndex]], .node)

  ggDATA <- melt(ftreeObj$functions)
  ggBASE <- ggplot(ggDATA) + geom_line(aes(x=Var1, y=value, group = Var2), colour = "black", lwd = 1.3, alpha = 0.7)

  if(length(.indices) > 1){
      ggBASE <- ggBASE + geom_line(data=melt(ftreeObj$functions[,.indices$rightInd]),
                                   aes(x=Var1, y=value, group = Var2, colour = "blue"), lwd = 1.15, alpha = 0.9)
      ggBASE <- ggBASE + geom_line(data=melt(ftreeObj$functions[,.indices$leftInd]),
                         aes(x=Var1, y=value, group = Var2, colour = "red"), lwd = 1.15, alpha = 0.9)
      ggBASE <- ggBASE + scale_colour_manual(values=c("black","blue","red")) + guides(colour = FALSE)
  } else {
      ggBASE <- ggBASE + geom_line(data=melt(ftreeObj$functions[,.indices$centerInd]),
                                   aes(x=Var1, y=value, group = Var2, colour = "orange"), lwd = 1.15, alpha = 0.9)
      ggBASE <- ggBASE + scale_colour_manual(values=c("black","orange"))
  }

  if(.ggRet == FALSE){
    print(ggBASE)
  } else {
    return(ggBASE)
  }
}

#' @title matrix to list
#'
#' @param .matrix  - A matrix you wish to divide into a list
#' @param .chunksize - Number of rows of matrix .matrix to be contained in each
#'                     element of the resulting list.
#'
#' @details This function simply splits provided matrix into N elements each containing
#'          .chunksize of rows of the original matrix. List elements are mutually exclusive
#'
#' @author Ogy Grujic (\email{ogyg@stanford.edu})
#'
#' @export
#'
matToList <- function(.matrix, .chunksize){

  if(.chunksize >= nrow(.matrix)){
    outList = list(indices=.matrix)
  } else {
    outList <- vector('list', ceiling(nrow(.matrix)/.chunksize))
    for(i in 1:length(outList)){
      if(nrow(.matrix) < .chunksize){
        outList[[i]] <- .matrix
      } else {
        outList[[i]] <- .matrix[1:.chunksize,]
        .matrix <- .matrix[-c(1:.chunksize), , drop=FALSE]
      }
    }
  }
  return(outList)

}

#' @title Predict function for ftree
#'
#' @param ftreeObj - An object of class "ftree" as produced by the \code{ftree} function
#' @param .Xnew    - A data frame with a new set of covariates for which to make a prediction(s).
#'
#' @details The code returns a list of length nrow(.Xnew). Each element of the list is a matrix whose
#'          number of rows equals the number of the original argument values (discretizations) of the
#'          functions contained in ftreeObj, while the number of columns is equal to the number of
#'          trees contained in ftreeObj (i.e. for a single tree mode it's a vector).
#'          This function does not work with trees of type "wss" and "rdist".
#'
#' @author Ogy Grujic \email{ognjengr@gmail.com}
#'
#' @export
#'
predictFtree <- function(ftreeObj = NULL, .Xnew = NULL){

  if(is.null(ftreeObj) | class(ftreeObj)[1] != "fTree") stop('Ftree Object is irregular or not provided. Stopping!')
  if(is.null(.Xnew)) stop("New set of covariates were not provided. Stopping!")
  if(!is.data.frame((.Xnew))) stop("Xnew should be a data frame!")

  match.arg(ftreeObj$costType, c("mahalanobis", "sse", "l2norm"))

  # this should change to include a situation where Xnew has more than the number of covariates in ftreeObj.
  # if that makes sense at all.
  if(length(setdiff(colnames(ftreeObj$covariates), colnames(.Xnew))) > 0) {
    stop("Covariates in Xnew don't match covariates in ftreeObj!")
  }

  # Makes a prediction of one single tree.
  .treePred <- function(.node, .xnew){
         if(.node$isLeaf == 1){

           temp <- ftreeObj$functions[,(.node$indices + 1)]

           if(is.null(ncol(temp))){
             return(temp)
           } else {
             return(rowMeans(temp))
           }
           # if(is.null(nrow(ftreeObj$functions[,.node$indices]))) {
           #    return(ftreeObj$functions[,(.node$indices + 1)])
           # } else {
           #    return(rowMeans(ftreeObj$functions[,(.node$indices + 1)]))
           # }

         } else {

           # see if its right or left:
           .predictor <- .node$bestPredictor + 1

           if(length(grep(",", .node$splitPoint)) > 0){ # it's categorical:
              # create a set out of .node$splitPoint
              # check if category present in .xnew is inside the set.
              # if yes then its left, if not then its right.

           } else {
             if(.xnew[.predictor] < as.numeric(.node$splitPoint)){

                .treePred(.node$left, .xnew)

             } else {

                .treePred(.node$right, .xnew)

             }
           }
         }
  }

  # computes predictions for each row of Xnew (if applicable).
  # Xnew should be provided as a data frame even if only one row is provided.
  .predictions <- apply(.Xnew, 1, function(x){
    ldply(ftreeObj$trees, function(y){
      .treePred(y, x)
    })
  })

  # n_trees <- length(ftreeObj$trees)
  # # .predictions <- vector("list", nrow(.Xnew))
  # .predictions <- apply(.Xnew, 1, function(x){
  #                       .preds <- vector("list", n_trees)
  #                       for(j in 1:n_trees){
  #                         .preds[[j]] <- .treePred(ftreeObj$trees[[j]], x)
  #                       }
  #                       return(.preds)
  #                     })


  return(.predictions)

}


#' @title Produces a sensitivity/importance plot
#'
#' @param ftreeObj - An object of class "ftree"
#' @param mode     - Mode you wish to run (plot - just plots, ggRet - returns a ggplot structure, data - returns data)
#'
#' @details Simply computes sensitivities based on cost function reduction (no permutation is performed ever).
#'          For plotting, this code uses ggplot graphics. If you wish to modify any part of the plot set mode="ggRet", the code
#'          will then return a ggplot structure that you can modify. Refer to ggplot documentation for more details.
#'
#' @author Ogy Grujic (\email{ognjengr@gmail.com})
#'
#' @export
#'
varSensitivity <- function(ftreeObj = NULL, mode = "plot"){

  if(is.null(ftreeObj)) stop("ftreeObject was not provided. Exiting!")
  if(class(ftreeObj)[1] != "fTree") stop('Passed object is not of class "fTree". Exiting!')

  mode <- match.arg(mode, c("plot", "ggRet", "data"))

  .importance <- matrix(0, nrow = ncol(ftreeObj$covariates), ncol = length(ftreeObj$trees))
  rownames(.importance)  <- colnames(ftreeObj$covariates)

  .extractSensitivities <- function(.node) {
    if(.node$isLeaf == 1) {
      return
    } else {
      .importance[,i] <<- .importance[,i] + (.node$nodeGoodness - .node$splitsGoodness[1,])*(length(.node$indices)/nrow(ftreeObj$covariates))
    }
  }

  for(i in 1:length(ftreeObj$trees)) {
    .extractSensitivities(ftreeObj$trees[[i]])
  }

  .importance[is.nan(.importance)] <- 0

  .ggDATA <- apply(.importance, 1, mean)

  .ggDATA <- .ggDATA/max(.ggDATA)
  .ggDATA <- rev(sort(.ggDATA))

  .ggDATA <- data.frame(parameter = names(.ggDATA), value=.ggDATA)
  .ggDATA$rank <- factor(1:nrow(.ggDATA))
  .ggDATA$parameter <- factor(.ggDATA$parameter, levels = rev(as.character(.ggDATA$parameter)))

  .ggP <- ggplot(.ggDATA, aes(x= parameter, y=value, fill=rank)) +
    geom_bar(stat = "identity", position = "dodge", lwd = 0.2, colour = "black") +
    coord_flip() + guides(fill=FALSE) + ggtitle('Variable importance plot') + xlab('relative importance')

  if(mode == "plot") {
    print(.ggP)
  } else if(mode == "ggRet") {
    return(.ggP)
  } else {
    return(.ggDATA)
  }

}



