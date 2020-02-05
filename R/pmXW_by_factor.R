##### Compute factor-specific posterior mean of reconstructed data (X*t(W)) #####
pmXW_by_factor <- function(gfa.obj){
  comps <- list()
  for (k in 1:res$K) {
    comp <- crossprod(res$posterior$X[,,k], res$posterior$W[,,k])
    comps[[k]] <- comp / res$opts$iter.saved
  }
  return(comps)
}
