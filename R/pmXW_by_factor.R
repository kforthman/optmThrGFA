#' Compute factor-specific posterior mean of reconstructed data (X*t(W)).
#'
#' @param gfa.obj A GFA object as output by the gfa() function in the package GFA.
#' @export

pmXW_by_factor <- function(gfa.obj){
  comps <- list()
  for (k in 1:gfa.obj$K) {
    comp <- crossprod(gfa.obj$posterior$X[,,k], gfa.obj$posterior$W[,,k])
    comps[[k]] <- comp / gfa.obj$opts$iter.saved
  }
  return(comps)
}
