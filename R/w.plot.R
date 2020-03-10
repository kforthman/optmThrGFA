#' Create a heatmap.
#'
#' A base function to create a heat map.
#'
#' @param w No description.
#' @param D No description.
#' @param K No description.
#' @param gr1 No description.
#' @param conf.level No description.
#' @param replicate No description.
#' @export

w.plot <- function(w, D, K, gr1, conf.level, replicate){
  mar <- c(6,4,4,6)
  par(mar=mar)
  cols <- colorRampPalette(c("orange","red","white","blue","cyan"))(19)
  if(any(is.na(w))) cols <- colorRampPalette(c("orange","red","#DDDDDD","blue","cyan"))(19)
  M <- max(abs(w),na.rm=T)
  breaks <- seq(-M,M,length=20)

  title <- c("Matrix W^T","Components","Features")
  if (!is.null(replicate)){
    title[1] <- paste0('Replicate ', replicate,': ', title[1])
  } else {
    if (is.null(conf.level)){
      title[1] <- paste0(title[1], ' (all components & ', round(sum(w!=0)), ' loadings)')
    } else if (!is.null(conf.level)){
      title[1] <- paste0(title[1], ' (', sum(w!=0),' non-zero loadings at ', conf.level*100, '% confidence)')
    }
  }

  if (K==1){
    image(as.matrix(w[,1]), col=cols, breaks=breaks, axes=F, main=title[1],
          xlab="", ylab="")
  } else {
    image(1:D, 1:K, w[,K:1], col=cols, breaks=breaks, axes=F, main=title[1],
          xlab="",ylab="")
  }
  title(xlab=title[3],line=mar[1]-1)
  title(ylab=title[2],line=mar[2]-1)
  box()
  par(las=2)
  if (K == 1){
    axis(1, (0:(D-1))/D, rownames(w), cex.axis=D^(-1/5))
    axis(2, K:1, colnames(w), cex.axis=K^(-1/5))
  } else {
    axis(1, 1:D, rownames(w), cex.axis=D^(-1/5))
    axis(2, K:1, colnames(w), cex.axis=K^(-1/5))
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
