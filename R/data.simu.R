#' Simulates observed data given the posterior, W, the block
#' assignments for each variable and the
#' number of observations, N.
#'
#' @param N sample size
#' @param W_DxK a list of matrices with K x D_{i} dimension
#' @param varIdx.by.block Specifies which variables are in which block. Must be input as a list where each element corresponds to the block number and contains a vector of indicies of the variables that belong to that block.
#' @param sd.noise SD of random Gaussian noise in the observations.
#' @return A list containing the following elements:
#' \item{Y.list}{A list where each element corresponds to a block. Each block contains the observed values of the variables in that block.}
#' \item{X}{Contains the factor values of each observation.}
#' @export

data.simu <- function(N, W_DxK, varIdx.by.block, sd.noise){
  B = length(varIdx.by.block)
  K = ncol(W_DxK)
  D = nrow(W_DxK)
  # X <- matrix(rnorm(N*K), N, K) ## 0 score mean not guarantteed
  X <- matrix(NA, N, K)
  for (k in 1:K){
    X[, k] <- rnorm(N) ## each factor (column) has a mean closer to 0
  }
  colnames(X) <- paste0("K", 1:K)
  X <- scale(X)
  Y <- X %*% t(W_DxK) + sd.noise*matrix(rnorm(N*D), N, D)
  colnames(Y) <- paste0("V", 1:D)
  Y.list <- list()
  for (b in 1:B){
    Y.list[[b]] <- Y[, varIdx.by.block[[b]]]
  }
  return(list(Y.list=Y.list, X=X))
}
