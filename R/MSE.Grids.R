#' A function to compute MSE.
#'
#' A function to compute MSE b/w observed and reconstructed data matrices for each replicate for
#' each pair of (corThr, matchThr).
#'
#' @param Ymtx No description.
#' @param corGrids No description.
#' @param matchGrids No description.
#' @inheritParams matchFactors
#' @export

MSE.Grids <- function(Ymtx, maxK, comps, corGrids, matchGrids){
  R = length(comps)

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
