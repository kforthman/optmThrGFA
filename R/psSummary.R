#' Functions to extract robust factors.
#'
#' Functions to extract robust factors
#' (0) Using the GFA::gfa() resulting object to compute
#' (1) posterior medians
#' (2) credible interval limits
#' (3) factor-wise reconstratced data (i.e. X * W)
#'
#' @param gfa.res No description.
#' @param credible.lv No description.
#' @export


psSummary <- function(gfa.res, credible.lv){
  p.vec <- c((1-credible.lv)/2, 0.5, (1+credible.lv)/2)
  # dimensions = 3 (quantiles) x D x K
  tmp <- apply(gfa.res$posterior$W, 2:3, function(x) quantile(x, p.vec) )
  gfa.res$W.Summ <- list(lo=tmp[1,,], p50=tmp[2,,], hi=tmp[3,,])

  rownames(gfa.res$W.Summ$lo) <- rownames(gfa.res$W.Summ$p50) <- rownames(gfa.res$W.Summ$hi) <- rownames(gfa.res$W)

  # dimensions = 3 (quantiles) x N x K
  tmp <- apply(gfa.res$posterior$X, 2:3, function(x) quantile(x, p.vec) )
  gfa.res$X.Summ <- list(lo=tmp[1,,], p50=tmp[2,,], hi=tmp[3,,])
  gfa.res$Yhat.p50 <- list()
  for (k in 1:dim(tmp)[3]){
    # if only one factor, need to redefine gfa.res$X.Summ$p50 as a matrix
    if (dim(tmp)[3] == 1) {
      dummyX <- as.matrix(gfa.res$X.Summ$p50);
      dummyW <- as.matrix(gfa.res$W.Summ$p50);
    } else {
      dummyX <- gfa.res$X.Summ$p50;
      dummyW <- gfa.res$W.Summ$p50;
    }
    gfa.res$Yhat.p50[[k]] <- matrix(dummyX[,k], ncol=1) %*%
      matrix(dummyW[,k], nrow=1) # N x D
  }
  names(gfa.res$Yhat.p50) <- paste0("K", 1:length(gfa.res$Yhat.p50))
  out <- gfa.res[names(gfa.res) %in% c("K", "W.Summ", "X.Summ", "Yhat.p50")]
  return(out)
}
