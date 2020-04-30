#' Flip signs of factor loadings to ensure factors in the same direction across replicates.
#'
#' @param models No description.
#' @param rob No description.
#' @param use.unmatched No description.
#' @export

w.signs <- function(models, rob, use.unmatched=F){
  indices <- abs(rob$indices)
  if (!use.unmatched){
    indices[rob$indices<0] <- 0
  }
  n.reps <- length(models)
  maxK <- ncol(indices)
  ref.rep <- which.max(apply(rob$indices>0, 1, sum))[1]
  rep.ind <- c(1:n.reps)[-ref.rep]
  for (r in rep.ind){
    sign.tmp <- diag(sign(cor(models[[ref.rep]]$W.Summ$p50[, indices[ref.rep,]],
                              models[[r]]$W.Summ$p50[, indices[r,]])))
    len.tmp <- length(sign.tmp) - maxK
    if(len.tmp<0){
      sign.tmp <- c(sign.tmp, rep(0, -len.tmp))
    } else if (len.tmp>0){
      sign.tmp <- sign.tmp[1:maxK]
    }
    indices[r, ] <- indices[r, ] * sign.tmp
    if (!use.unmatched){
      indices[indices==0] <- NA
    }
  }
  return(indices)
}
