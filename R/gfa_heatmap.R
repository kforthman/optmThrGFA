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

gfa_heatmap <- function(robW, block.names, varIdx.by.block, conf.level, heatmap.rep=FALSE, factor.order=NULL){
  n.rep <- max(robW$w.ci$Replicate)

  # heat maps
  gr <- varIdx.by.block
  M <- length(gr)

  if (is.null(block.names)) {
    names(gr) <- paste("Source",1:M)
  } else { names(gr) <- block.names }

  gr1 <- c(0,cumsum(sapply(gr, length))); names(gr1) <- c(names(gr), "NA")

  if (n.rep > 1 & heatmap.rep){
    for (r in 1:n.rep){
      w.tmp <- matrix(robW$w.ci$Median[robW$w.ci$Replicate==r]*(1-robW$w.ci$contain.0[robW$w.ci$Replicate==r]),
                      nrow(robW$w.med), ncol(robW$w.med))
      colnames(w.tmp) <- 1:ncol(robW$w.med)
      rownames(w.tmp) <- rownames(robW$w.med)
      w.plot(
        w = w.tmp,
        D = nrow(w.tmp),
        K = ncol(w.tmp),
        gr1 = gr1,
        conf.level = conf.level,
        replicate = r
        )
    }
  }

  if(!is.null(factor.order)){
    w.plot(
      w = robW$w.med[, factor.order],
      D = nrow(robW$w.med),
      K = ncol(robW$w.med),
      gr1 = gr1,
      conf.level = conf.level,
      replicate = NULL
    )
  } else {
    w.plot(
      w = robW$w.med,
      D = nrow(robW$w.med),
      K = ncol(robW$w.med),
      gr1 = gr1,
      conf.level = conf.level,
      replicate = NULL
      )
  }
}
