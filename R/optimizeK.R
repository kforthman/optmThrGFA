# a function to identify optimal (corThr, matchThr)
optimizeK <- function(K.grids, mse.array){ # c(N.corGrids, N.matchGrids, N.reps)
  reps = dim(mse.array)[3]
  mse.m <- apply(mse.array, 1:2, mean, na.rm=T)
  mse.m <- matrix(as.vector(mse.m), nrow=nrow(K.grids), ncol=ncol(K.grids),
                  dimnames=list(corThr=rownames(K.grids),
                                matchThr=colnames(K.grids)))
  corGrids <- as.numeric(rownames(K.grids))
  matchGrids <- as.numeric(colnames(K.grids))

  # find K that minimize MSE
  mse.min <- min(mse.m, na.rm=T)
  ind.min <- which(mse.m==mse.min, arr.ind=T)
  Krobust.min <- K.grids[ind.min][1]

  # find K using 1-SE rule
  min.mse.se <- sd(mse.array[ind.min[1,1], ind.min[1,2], ], na.rm=T)
  mseThr <- mse.m[ind.min[1,1], ind.min[1,2]] + min.mse.se
  K.grids.tmp <- K.grids
  K.grids.tmp[mse.m > mseThr | K.grids.tmp==0] <- NA
  Krobust.1se <- min(K.grids.tmp, na.rm=T)
  ind.1se <- which(K.grids.tmp == Krobust.1se, arr.ind=T)

  par.min <- cbind(corGrids[ind.min[,1]], matchGrids[ind.min[,2]], K.grids[ind.min])
  par.1se <- cbind(corGrids[ind.1se[,1]], matchGrids[ind.1se[,2]], K.grids[ind.1se])
  colnames(par.min) <- colnames(par.1se) <- c('opt.corThr', 'opt.matchThr', 'optK')

  return(list(mse.m=mse.m, mse.min=mse.min,
              min.mse.se=min.mse.se, mseThr=mseThr,
              par.min=par.min, Krobust.min=Krobust.min,   ## min-MSE criterion
              par.1se=par.1se, Krobust.1se=Krobust.1se))  ## 1-SE criterion
}
