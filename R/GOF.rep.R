#' A function to compute goodness-of-fit on W matrix elements
#'
#' A function to compute goodness-of-fit on W matrix elements. It requires matching factors between estimated and true W matriices
#'
#' @param gfaList No description.
#' @param parm No description.
#' @param trueMtx No description.
#' @export

GOF.rep <- function(gfaList, parm, trueMtx){
  R = length(gfaList)
  K = ncol(trueMtx)
  out <- data.frame(Replicate = 1:R)
  out$RMSE <- out$cor <- NA
  for (rep in 1:R){
    if(parm=="W"){
      est <- gfaList[[rep]]$W.Summ$p50
    } else { est <- gfaList[[rep]]$X.Summ$p50 }
    # est <- ifelse(parm=="W", gfaList[[rep]]$W.Summ$p50, gfaList[[rep]]$X.Summ$p50)
    tmp <- GOF(est=est, true=trueMtx)
    out[rep, c("cor", "RMSE")] <- tmp
  }
  out$Replicate <- as.factor(out$Replicate)
  return(out)
}
