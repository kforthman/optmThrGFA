#' A function to carry out matching in parallel.
#'
#' @param comps A list of GFA replicates with posterior medians.
#' @param corGrids A list of cor thresholds to test. Cor threshold defines how close two components are required to be, in terms of correlation, in order to match them.
#' @param matchGrids A list of match thresholds to test. Match threshold describes the proportion of sampling chains that need to contain a component in order to include it in the robust components.
#' @export
#'

match_dopar <- function(comps, corGrids, matchGrids) {
  reps <- length(comps)               # number of GFA replicates
  maxK <- max(sapply(comps, length))  # largest K across replicates
  # create a list of vectors
  compStorage <- vector("list", length=reps) #Store the components that can still be used
  for(rep in 1:reps) compStorage[[rep]] <- which(sapply(comps[[rep]], function(x) {sum(abs(x)) > 0}))

  rob <-
    foreach (i = 1:length(corGrids) ) %:%
    foreach (j = 1:length(matchGrids) ) %dopar% {
      tmp <- matchFactors(comps, maxK=maxK, corThr=corGrids[i], matchThr=matchGrids[j])
      print(paste0("Elapsed time for (corThr, matchThr) = (", corGrids[i], ", ", matchGrids[j], ") is ",
                   round(tmp$elapsed.time/60, 1), " minutes"))
      if(is.null(tmp)){stop(paste0("matchFactors() at the threshold (corThr, matchThr) = (", corGrids[i], ", ", matchGrids[j], ") returned NULL."))}
      return(tmp)
    }
  names(rob) <- paste0("corThr_", corGrids)
  for (i in 1:length(corGrids)){
    names(rob[[i]]) <- paste0("matchThr_", matchGrids)
  }
  return(rob)
}
