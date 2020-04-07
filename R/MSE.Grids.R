#' A function to compute MSE of the GFA given robust components thresholds.
#'
#' A function to compute MSE b/w observed and reconstructed data matrices for each replicate for
#' each pair of (corThr, matchThr).
#'
#' @param Ymtx The set of observed data as an NxD matrix.
#' @param maxK The maximal K among the GFA replicates.
#' @param comps A list of GFA replicates with posterior medians.
#' @param corGrids is a vector of cor thresholds to be tested. Cor threshold defines how close two components are required to be, in terms of correlation, in order to match them.
#' @param matchGrids is a vector of match thresholds to be tested. Match threshold describes the proportion of sampling chains that need to contain a component in order to include it in the robust components.
#' @inheritParams matchFactors
#' @export

MSE.Grids <- function(Ymtx, maxK, comps, corGrids, matchGrids){
  R <- length(comps)

  indices <- rep( list(list()), length(corGrids) )
  K.grids <-
    matrix(NA, nrow=length(corGrids), ncol=length(matchGrids),
           dimnames=list(corThr=corGrids, matchThr=matchGrids))
  mse <- list(all=NULL, survived=NULL)
  mse$all <- mse$survived <- array(NA, c(length(corGrids), length(matchGrids), R))

  for (i in 1:length(corGrids)){
    for (j in 1:length(matchGrids)){
      print(paste(Sys.time(),": ",i,",",j))

      tmp <- matchFactors(
        comps = comps,
        maxK = maxK,
        corThr = corGrids[i],
        matchThr = matchGrids[j]
        )
      if(ncol(tmp$indices)>0){
        indices[[i]][[j]] <- tmp$indices
        K.grids[i,j] <- tmp$Krobust # ifelse(tmp$Krobust > max(K.by.rep), K.by.rep, tmp$Krobust)
        if (K.grids[i,j]>0){ # at least 1 factor
          for (r in 1:R){
            ## use all "matched" factors
            indx.use <- tmp$indices[r,]
            yhat <- Reduce("+", comps[[r]][abs(indx.use)])
            mse$all[i,j,r] <- mean((Ymtx-yhat)^2, na.rm=T)
            ## use the "matched" factors surviving corThr
            indx.use <- tmp$indices[r,][tmp$indices[r,] >0]
            yhat <- Reduce("+", comps[[r]][indx.use])
            mse$survived[i,j,r] <- mean((Ymtx-yhat)^2, na.rm=T)
          }
        }
      }
      rm(tmp)
    }
  }
  return(list(K.grids=K.grids, indices=indices, mse=mse))
}
