#' A function to compute goodness-of-fit on W matrix elements
#'
#' A function to compute goodness-of-fit on W matrix elements. It requires matching factors between estimated and true W matriices
#'
#' @param est No description.
#' @param true No description.
#' @export

GOF <- function(est, true){
  nonZero.k <- apply(est, 2, sd)>0
  K = sum(nonZero.k)
  if(K==0){
    return(out=c(cor=NA, RMSE=NA))
  } else{
    matched <- est[, nonZero.k]
    if(K==1){
      r = cor(est[, nonZero.k], true[,1])
      matched <- matched*sign(r)
      out = c(
        cor(matched, true[, 1]),
        sqrt(mean((matched-true[, 1])^2))
      )
    } else if (K>1){
      idx <- factor.sign <- rep(NA, K)
      for (k in 1:K){
        r <- apply(est[, nonZero.k], 2, function(x){ cor(x, true[,k])})
        idx[k] <- which.max(abs(r))
        factor.sign[k] <- sign(r)[idx[k]]
      }
      for (k in 1:K){
        matched[,k] <- est[, nonZero.k][, idx[k]]*factor.sign[k]
      }
      out = c(
        cor(as.vector(matched), as.vector(true[, 1:K])),
        sqrt(mean((as.vector(matched)-as.vector(true[, 1:K]))^2))
      )
    }
    return(out)
  }
}
