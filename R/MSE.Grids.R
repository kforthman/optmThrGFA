#' A function to compute MSE of the GFA given robust components thresholds.
#'
#' A function to compute MSE b/w observed and reconstructed data matrices for each replicate for
#' each pair of (corThr, matchThr).
#'
#' @param Ymtx The set of observed data as an NxD matrix.
#' @param comps A list of GFA replicates with posterior medians.
#' @param match.res No description.
#' @param corGrids is a vector of cor thresholds to be tested. Cor threshold defines how close two components are required to be, in terms of correlation, in order to match them.
#' @param matchGrids is a vector of match thresholds to be tested. Match threshold describes the proportion of sampling chains that need to contain a component in order to include it in the robust components.
#' @inheritParams matchFactors
#' @export

MSE.Grids <- function(Ymtx, comps, match.res, corGrids, matchGrids){
  R <- length(comps)

  K.grids <-
    matrix(NA, nrow=length(corGrids), ncol=length(matchGrids),
           dimnames=list(corThr=corGrids, matchThr=matchGrids))
  mse <- list(all=NULL, survived=NULL)
  mse$all <- mse$survived <- array(NA, c(length(corGrids), length(matchGrids), R))

  for (i in 1:length(corGrids)){
    for (j in 1:length(matchGrids)){

      if(ncol(match.res[[i]][[j]]$indices)>0){
        K.grids[i,j] <- match.res[[i]][[j]]$Krobust # ifelse(tmp$Krobust > max(K.by.rep), K.by.rep, tmp$Krobust)
        if (K.grids[i,j]>0){ # at least 1 factor
          for (r in 1:R){
            ## use all "matched" factors
            indx.use <- match.res[[i]][[j]]$indices[r,]
            yhat <- Reduce("+", comps[[r]][abs(indx.use)])
            mse$all[i,j,r] <- mean((Ymtx-yhat)^2, na.rm=T)
            ## use the "matched" factors surviving corThr
            indx.use <- match.res[[i]][[j]]$indices[r,][match.res[[i]][[j]]$indices[r,] >0]
            yhat <- Reduce("+", comps[[r]][indx.use])
            mse$survived[i,j,r] <- mean((Ymtx-yhat)^2, na.rm=T)
          }
        }
      }

    }
  }
  return(list(K.grids=K.grids, mse=mse))
}
