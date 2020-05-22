#' GFA Heatmap.
#'
#' A function to produce a heatmap for robust factor loadings.
#'
#'
#' @param robW No description.
#' @param block.names No description.
#' @param varIdx.by.block No description.
#' @param heatmap.rep No description.
#' @param factor.order No description.
#' @inheritParams w.plot
#' @export

gfa_heatmap <- function(robW, block.names, varIdx.by.block, conf.level,
                        sparse=F, heatmap.rep=FALSE, factor.order=NULL){
  n.rep <- max(robW$w.ci$Replicate)
  # heat maps
  gr <- varIdx.by.block; M <- length(gr)
  if (is.null(block.names)) {
    names(gr) <- paste("Source",1:M)
  } else { names(gr) <- block.names }
  gr1 <- c(0,cumsum(sapply(gr, length))); names(gr1) <- c(names(gr), "NA")
  if (n.rep > 1 & heatmap.rep){
    for (r in 1:n.rep){
      if(sparse){
        w.tmp <- matrix(robW$w.ci$Median[robW$w.ci$Replicate==r]*(1-robW$w.ci$contain.0[robW$w.ci$Replicate==r]),
                        nrow(robW$w.med), ncol(robW$w.med))
      } else {
        w.tmp <- matrix(robW$w.ci$Median[robW$w.ci$Replicate==r],
                        nrow(robW$w.med), ncol(robW$w.med))
      }
      colnames(w.tmp) <- 1:ncol(robW$w.med)
      rownames(w.tmp) <- rownames(robW$w.med)
      w.plot(w.tmp, D=nrow(w.tmp), K=ncol(w.tmp), gr1, conf.level, r)
    }
  }
  # print('Robust heat map')
  if(!is.null(factor.order)){
    w.plot(robW$w.med[, factor.order], D=nrow(robW$w.med), K=length(factor.order), gr1, conf.level, replicate=NULL)
  } else {

    dummy <- robW$w.med
    rowCounter <- 1
    for (LIST in 1:length(gr)) {
      for (VAR in 1:length(gr[[LIST]])) {
        rownames(dummy)[rowCounter] <- gr[[LIST]][VAR]
        rowCounter = rowCounter + 1;
      }
    }

    w.plot(robW$w.med, D=nrow(robW$w.med), K=ncol(robW$w.med), gr1, conf.level, replicate=NULL)
  }
}
