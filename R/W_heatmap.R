#' For true W matrix parameters.
#'
#' @param W_DxK No description.
#' @param varIdx.by.block No description.
#' @param block.names No description.

W_heatmap <- function(W_DxK, varIdx.by.block, block.names){
  gr <- varIdx.by.block
  M <- length(gr)
  names(gr) <- block.names
  gr1 <- c(0, cumsum(sapply(gr, length)))
  names(gr1) <- c(names(gr), "NA")
  D <- nrow(W_DxK)
  K <- ncol(W_DxK)

  mar <- c(6,4,4,6)
  par(mar=mar)
  cols <- colorRampPalette(c("orange","red","white","blue","cyan"))(19)
  if(any(is.na(W_DxK))) cols <- colorRampPalette(c("orange","red","#DDDDDD","blue","cyan"))(19)
  M <- max(abs(W_DxK), na.rm=T)
  breaks <- seq(-M, M, length=20)

  title <- c("Matrix W^T", "Factors", "Features")

  if (K==1){
    image(as.matrix(W_DxK[,1]), col=cols, breaks=breaks, axes=F, main=title[1],
          xlab="", ylab="")
  } else {
    image(1:D, 1:K, W_DxK[,K:1], col=cols, breaks=breaks, axes=F, main=title[1],
          xlab="",ylab="")
  }
  title(xlab=title[3],line=mar[1]-1)
  title(ylab=title[2],line=mar[2]-1)
  box()
  par(las=2)
  if (K == 1){
    axis(1, (0:(D-1))/D, rownames(W_DxK), cex.axis=D^(-1/5))
    axis(2, K:1, colnames(W_DxK), cex.axis=K^(-1/5))
  } else {
    axis(1, 1:D, rownames(W_DxK), cex.axis=D^(-1/5))
    axis(2, K:1, colnames(W_DxK), cex.axis=K^(-1/5))
  }

  #Grouping
  par(xpd=T)
  mu <- gr1[-1]/2+gr1[-length(gr1)]/2
  N <- K
  for(i in 1:length(mu)) {
    if (K ==1){
      if(i!=length(mu)) lines(rep(gr1[i+1]-0.5,2)/D, c(-1, 1.05), lwd=2)
      text(mu[i]/D,1.065,names(gr1)[i])
    } else {
      if(i!=length(mu)) lines(rep(gr1[i+1]+1/2,2), c(.5, N*1.03+.5), lwd=2)
      text(mu[i],N*1.03+.5,names(gr1)[i])
    }
  }
  #Colorbar
  n <- length(cols)
  if (K==1){
    cba <- 1.1
    cbw <- 1/D
    for(i in 1:n){
      polygon(c(0,cbw,cbw,0)+cba, (c(0,0,N/n,N/n)+N*(i-1)/n+1/2)-1,
              col=cols[i], border=NA)
    }
    #Colorbar: axis
    lines(rep(cba+cbw,2),c(0,N)+1/2-1)
    m <- 10^floor(log10(M))
    m <- floor(M/m)*m
    for(l in c(-m,0,m)) {
      ly <- N*(l/M/2+.5)+1/2-1
      lines(cba+cbw-c(cbw,-cbw)/5, rep(ly,2))
      text(cba+cbw*2.5+0.02,ly,l)
    }
  } else {
    cba <- D + 1/2 + D/60
    cbw <- D/40
    for(i in 1:n){
      polygon(c(0,cbw,cbw,0)+cba, c(0,0,N/n,N/n)+N*(i-1)/n+1/2,
              col=cols[i], border=NA)
    }
    #Colorbar: axis
    lines(rep(cba+cbw,2),c(0,N)+1/2)
    m <- 10^floor(log10(M)); m <- floor(M/m)*m
    for(l in c(-m,0,m)) {
      ly <- N*(l/M/2+.5)+1/2
      lines(cba+cbw-c(cbw,-cbw)/5, rep(ly,2))
      text(cba+cbw*2.5,ly,l)
    }
  }
  par(xpd=F)
}
