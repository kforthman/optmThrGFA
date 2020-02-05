### 4b. a function to produce a heatmap for robust factor loadings
gfa_heatmap <- function(robW, block.names, varIdx.by.block, conf.level, heatmap.rep=FALSE, factor.order=NULL){
  n.rep <- max(robW$w.ci$Replicate)
  # heat maps
  gr <- varIdx.by.block; M <- length(gr)
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
      w_plot(w.tmp, D=nrow(w.tmp), K=ncol(w.tmp), gr1, conf.level, r)
    }
  }
  # print('Robust heat map')
  if(!is.null(factor.order)){
    w_plot(robW$w.med[, factor.order], D=nrow(robW$w.med), K=ncol(robW$w.med), gr1, conf.level, replicate=NULL)
  } else {
    w_plot(robW$w.med, D=nrow(robW$w.med), K=ncol(robW$w.med), gr1, conf.level, replicate=NULL)
  }
}
