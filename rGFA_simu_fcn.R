##### For true W matrix parameters #####
W_heatmap <- function(W_DxK, varIdx.by.block, block.names){
  gr <- varIdx.by.block; M <- length(gr)
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

##### for data simulation #####
## N: sample size
## tW.list: a list of matrices with K x D_{i} dimension
## sd.noise: SD of random Gaussian noise
data_simu <- function(N, W_DxK, varIdx.by.block, sd.noise){
  B = length(varIdx.by.block)
  K = ncol(W_DxK)
  D = nrow(W_DxK)
  # X <- matrix(rnorm(N*K), N, K) ## 0 score mean not guarantteed
  X <- matrix(NA, N, K)
  for (k in 1:K){
    X[, k] <- rnorm(N) ## each factor (column) has a mean closer to 0
  }
  colnames(X) <- paste0("K", 1:K)
  Y <- X %*% t(W_DxK) + sd.noise*matrix(rnorm(N*D), N, D)
  colnames(Y) <- paste0("V", 1:D)
  Y.list <- list()
  for (b in 1:B){
    Y.list[[b]] <- Y[, varIdx.by.block[[b]]]
  }
  return(list(Y.list=Y.list, X=X))
}

#### Revise robustComponents() to allow 0 matched factors
## Internal function for matching components
matchComponents <- function(comps, maxK, N, D, corThr, matchThr) {
  corList <- list()
  reps <- length(comps)
  rob <- list(Krobust=0,effect=c(),indices=matrix(NA,reps,0),cor=matrix(NA,reps,0))

  compStorage <- vector("list",length=reps) #Store the components that can still be used
  for(rep in 1:reps) compStorage[[rep]] <- which(sapply(comps[[rep]], function(x) {sum(abs(x)) > 0}))

  for (rep1 in 1:reps) {
    matching <- vector(length = maxK)
    sim <- array(NA, dim = c(maxK, reps))
    matchInd <- matrix(0, nrow = maxK, ncol = reps)

    for (k1 in compStorage[[rep1]]) {
      for (rep2 in which(sapply(compStorage, length)>0)) {
        #Correlation between the two components.
        #Note that we do not need to use the absolute values, as the correlation is computed in data space.
        d <- sapply(comps[[rep2]][compStorage[[rep2]]], function(x) cor(c(x), c(comps[[rep1]][[k1]])))
        sim[k1,rep2] <- max(d)
        matchInd[k1,rep2] <- compStorage[[rep2]][which.max(d)]
        if(sim[k1,rep2] < corThr) matchInd[k1,rep2] <- matchInd[k1,rep2]*-1 #Index saved for debugging purposes
      }

      if(sum(sim[k1,]>corThr, na.rm=T) >= matchThr*reps) { #Robust component found!
        # average over all similar components
        comp <- matrix(0, N, D)
        for(rep2 in which(matchInd[k1,] > 0)) {
          comp <- comp + comps[[rep2]][[matchInd[k1,rep2]]]
          compStorage[[rep2]] <- compStorage[[rep2]][!compStorage[[rep2]]==matchInd[k1,rep2]]
        }
        comp <- comp/sum(matchInd[k1,]>0)

        rob$Krobust <- rob$Krobust + 1
        rob$effect <- array(c(rob$effect, comp), dim=c(dim(comp),rob$K),
                            dimnames=list(NA, paste0("V", 1:D),
                                          paste0("K",1: rob$Krobust)))
        rob$indices <- cbind(rob$indices, matchInd[k1,])
        rob$cor <- cbind(rob$cor, sim[k1,])
      }
    }
  }
  if(rob$Krobust>0){
    rownames(rob$indices) <- rownames(rob$cor) <- paste0("rep",1:reps)
    colnames(rob$indices) <- colnames(rob$cor) <- paste0("K",1:rob$Krobust)
  }
  return(rob)
}

##### Compute factor-specific posterior mean of reconstructed data (X*t(W)) #####
pmXW_by_factor <- function(gfa.obj){
  comps <- list()
  for (k in 1:res$K) {
    comp <- crossprod(res$posterior$X[,,k], res$posterior$W[,,k])
    comps[[k]] <- comp / res$opts$iter.saved
  }
  return(comps)
}

##### a function adopted from GFA::matchComponents to match factors #####
matchFactors <- function(comps, maxK, corThr, matchThr){
  reps <- length(comps)
  rob <- list(Krobust=0, indices=matrix(NA,reps,0), cor=matrix(NA,reps,0))

  compStorage <- vector("list", length=reps) #Store the components that can still be used
  for(rep in 1:reps) compStorage[[rep]] <- which(sapply(comps[[rep]], function(x) {sum(abs(x)) > 0}))

  for (rep1 in 1:reps) {
    sim <- array(NA, dim = c(maxK, reps))
    matchInd <- matrix(0, nrow = maxK, ncol = reps)

    for (k1 in compStorage[[rep1]]) {
      for (rep2 in which(sapply(compStorage, length)>0)) {
        #Correlation between the two components.
        #Note that we do not need to use the absolute values, as the correlation is computed in data space.
        d <- sapply(comps[[rep2]][compStorage[[rep2]]], function(x) cor(c(x), c(comps[[rep1]][[k1]])))
        sim[k1,rep2] <- max(d)
        matchInd[k1,rep2] <- compStorage[[rep2]][which.max(d)]
        if(sim[k1,rep2] < corThr) matchInd[k1,rep2] <- matchInd[k1,rep2]*-1 #Index saved for debugging purposes
      }

      if(sum(sim[k1, ]>corThr, na.rm=T) - 1 >= matchThr*(reps-1)) { #Robust component found!
        rob$Krobust <- rob$Krobust + 1
        rob$indices <- cbind(rob$indices, matchInd[k1,])
        rob$cor <- cbind(rob$cor, sim[k1,])
        for(rep2 in which(matchInd[k1,] > 0)) {
          compStorage[[rep2]] <- compStorage[[rep2]][!compStorage[[rep2]]==matchInd[k1,rep2]]
        }
      }
    }
  }
  if(rob$Krobust>0){
    rownames(rob$indices) <- rownames(rob$cor) <- paste0("rep",1:reps)
    colnames(rob$indices) <- colnames(rob$cor) <- paste0("K", 1:rob$Krobust)
  }
  return(rob)
}

# a fcuntion to compute MSE b/w observed and reconstructed data matrices for each replicate for
#   each pair of (corThr, matchThr).
MSE_Grids <- function(Ymtx, maxK, comps, corGrids, matchGrids){
  R = length(comps)

  indices <- rep( list(list()), length(corGrids) )
  K.grids <-
    matrix(NA, nrow=length(corGrids), ncol=length(matchGrids),
           dimnames=list(corThr=corGrids, matchThr=matchGrids))
  mse <- list(all=NULL, survived=NULL)
  mse$all <- mse$survived <- array(NA, c(length(corGrids), length(matchGrids), R))

  for (i in 1:length(corGrids)){
    for (j in 1:length(matchGrids)){
      print(paste(Sys.time(),": ",i,",",j))

      tmp <- matchFactors(comps, maxK, corThr=corGrids[i], matchThr=matchGrids[j])
      if(ncol(tmp$indices)>0){
        indices[[i]][[j]] <- tmp$indices
        K.grids[i,j] <- tmp$Krobust # ifelse(tmp$Krobust > max(K.by.rep), K.by.rep, tmp$Krobust)
        if (K.grids[i,j]>0){ # at least 1 factor
          for (r in 1:R){
            ## use all "matched" factors
            indx.use <- tmp$indices[r,]
            yhat <- Reduce("+", comps[[r]][abs(indx.use)])
            mse$all[i,j,r] <- mean((Ymtx-yhat)^2, na.rm=T)
            ## use the "matched" factors surviving corThr
            indx.use <- tmp$indices[r,][tmp$indices[r,] >0]
            yhat <- Reduce("+", comps[[r]][indx.use])
            mse$survived[i,j,r] <- mean((Ymtx-yhat)^2, na.rm=T)
          }
        }
      }
      rm(tmp)
    }
  }
  return(list(K.grids=K.grids, indices=indices, mse=mse))
}
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

##### Functions to extract robust factors #####
### 0. Using the GFA::gfa() resulting object to compute
###    (1) posterior medians
###    (2) credible interval limits
###    (3) factor-wise reconstratced data (i.e. X * W)
psSummary <- function(gfa.res, credible.lv){
  p.vec <- c((1-credible.lv)/2, 0.5, (1+credible.lv)/2)
  # dimensions = 3 (quantiles) x D x K
  tmp <- apply(gfa.res$posterior$W, 2:3, function(x) quantile(x, p.vec) )
  gfa.res$W.Summ <- list(lo=tmp[1,,], p50=tmp[2,,], hi=tmp[3,,])
  # dimensions = 3 (quantiles) x N x K
  tmp <- apply(gfa.res$posterior$X, 2:3, function(x) quantile(x, p.vec) )
  gfa.res$X.Summ <- list(lo=tmp[1,,], p50=tmp[2,,], hi=tmp[3,,])
  gfa.res$Yhat.p50 <- list()
  for (k in 1:dim(tmp)[3]){
    gfa.res$Yhat.p50[[k]] <- matrix(gfa.res$X.Summ$p50[,k], ncol=1) %*%
      matrix(gfa.res$W.Summ$p50[,k], nrow=1) # N x D
  }
  names(gfa.res$Yhat.p50) <- paste0("K", 1:length(gfa.res$Yhat.p50))
  out <- gfa.res[names(gfa.res) %in% c("K", "W.Summ", "X.Summ", "Yhat.p50")]
  return(out)
}
### Flip signs of factor loadings to ensure factors in the same direction across replicates
w_signs <- function(models, rob, use.unmatched=F){
  indices <- abs(rob$indices)
  if (!use.unmatched){
    indices[rob$indices<0] <- 0
  }
  n.reps <- length(models)
  maxK <- ncol(indices)
  ref.rep <- which.min(apply(rob$indices<0, 1, sum))[1]
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

### a function to compute % variance explained
rob_var_exp <- function(models, indices, block.names, varIdx.by.block, use.unmatched=F, by.block=T){
  n.reps <- length(models)
  indices <- w_signs(models, indices, use.unmatched)
  K.rob <- ncol(indices)
  W.p50.rep <- array(NA, dim=c(nrow(models[[1]]$W.Summ$p50), K.rob, n.reps))
  for (r in 1:n.reps){
    idx_vec <- indices[r,]
    if (length(idx_vec)==1){
      W.p50.rep[,,r] <- sign(idx_vec)*models[[r]]$W.Summ$p50[, abs(idx_vec)]
    } else {
      W.p50.rep[,,r] <- sweep(models[[r]]$W.Summ$p50[, abs(idx_vec)],
                              MARGIN=2, sign(idx_vec), '*')
    }
  }
  ve.rep <- apply(W.p50.rep^2, 2:3, mean, na.rm=T)*100

  colnames(ve.rep) <- paste0('rep.', 1:n.reps)
  rownames(ve.rep) <- paste0('Component ', 1:K.rob)
  ve <- data.frame(matrix(NA, K.rob, 3))
  names(ve) <- c('Component', 'Mean', 'SE')
  ve$Component <- 1:K.rob
  ve$Mean <- apply(ve.rep, 1, mean, na.rm=T)
  ve$SE <- apply(ve.rep, 1, sd, na.rm=T)/sqrt(apply(!is.na(ve.rep), 1, sum))

  ve <- ve[order(-ve$Mean), ]
  ve <- within(ve, cum_var <- cumsum(ifelse(is.na(ve$Mean), 0, ve$Mean)))

  tot.ve.m <- mean(apply(W.p50.rep^2, 3, sum, na.rm=T))/dim(W.p50.rep)[1]*100
  tot.ve.s <- sd(apply(W.p50.rep^2, 3, sum, na.rm=T))/dim(W.p50.rep)[1]*100
  p <- ggplot(ve,
              aes(x=1:K.rob, y=Mean, ymin=Mean-SE, ymax=Mean+SE)) +
    geom_pointrange() +
    xlab('Robust components') + ylab('Percent variance explained') +
    theme_bw() # use a white background
  print(p + ggtitle(paste0('The ', K.rob, ' robust components explain ',
                           round(tot.ve.m,1), '+/-', round(tot.ve.s, 1),
                           '% variance of all variables')))

  # % variance explained by block
  if (by.block){
    nBlocks = length(block.names)
    ve.by.block <- ve.by.block.comp <- NULL # list()
    for (b in 1:nBlocks){
      ## variance explained by block by factor/component
      ve.tmp <- apply(W.p50.rep[varIdx.by.block[[b]],,]^2, 2:3, mean, na.rm=T)*100
      colnames(ve.tmp) <- paste0('rep.', 1:n.reps)
      rownames(ve.tmp) <- paste0('Component ', 1:K.rob)
      ve.comps <- data.frame(matrix(NA, K.rob, 3))
      names(ve.comps) <- c('Component', 'Mean', 'SE')
      ve.comps$Component <- 1:K.rob
      ve.comps$Mean <- apply(ve.tmp, 1, mean, na.rm=T)
      ve.comps$SE <- apply(ve.tmp, 1, sd, na.rm=T)/sqrt(apply(!is.na(ve.tmp), 1, sum))
      ve.comps$Block <- block.names[b]
      ve.by.block.comp <- rbind(ve.by.block.comp, ve.comps)

      ## variance explained by block
      ve.tmp <- apply(W.p50.rep[varIdx.by.block[[b]],,]^2, 3, sum, na.rm=T)/
        length(varIdx.by.block[[b]])
      m <- mean(ve.tmp)*100
      s <- sd(ve.tmp)/sqrt(n.reps)*100
      ve.by.block <- rbind(ve.by.block, c(m, s))
    }
    ve.by.block.comp$Block <- as.factor(ve.by.block.comp$Block)
    ve.by.block.comp$Block <- factor(ve.by.block.comp$Block, levels=block.names)
    ve.by.block.comp <- ve.by.block.comp[, c('Block', 'Component', 'Mean', 'SE')]

    ve.by.block <- data.frame(Block=block.names, ve.by.block)
    names(ve.by.block)[-1] <- c("Mean", "SE")

    return(list(
      indices=indices,
      ve=ve.rep,                 ## (Krobust x n.Reps): varaince explained (ve) per factor per replicate
      ve.summ=ve,                ## (Krobust rows): Mean and SE of ve per factor
      ve.by.block.comp=ve.by.block.comp, ## (Krobust x nBlocks rows): ve by block by factor
      ve.by.block=ve.by.block))  ## (nBlock rows): ve per block
  } else{
    return(list(indices=indices, ve=ve.rep, ve.summ=ve))
  }
}

### a function to compute robust loadings and robust scores
rob_wx <- function(models, indices, block.labs, var.labs=NULL){
  N <- nrow(models[[1]]$X.Summ$p50)
  if("W.Summ" %in% names(models[[1]])){
    D <- nrow(models[[1]]$W.Summ$p50)
  } else if ("W" %in% names(models[[1]])){
    D <- nrow(models[[1]]$W)
  }
  n.reps <- length(models)
  Krobust <- ncol(indices)

  # create a dataframe to store credible intervals of loadings across replicates
  # Henry edited on 2018-12-12: Add block names and preferred variable names
  df.base <- data.frame(Block=block.labs,
                        Variable=rownames(models[[1]]$W.Summ$p50),
                        var.lab=var.labs,
                        var.order=as.numeric(1:length(var.labs)))
  w.ci <- df.base[rep(seq_len(nrow(df.base)), n.reps*Krobust),]
  w.ci$Replicate <- rep(1:n.reps, each=D*Krobust)
  w.ci$Component <- rep(rep(1:Krobust, each=D), n.reps)
  w.ci <- w.ci[, c('Replicate', 'Component', 'Block', 'Variable', 'var.lab', 'var.order')]
  w.ci$Upper <- w.ci$Median <- w.ci$Lower <- 0

  # Extract posterior medians and credible intervals for loadings in each replicate
  for (r in 1:n.reps){
    for (k in 1:Krobust){
      w.ci[w.ci$Replicate==r & w.ci$Component==k, c('Lower', 'Median', 'Upper')] <-
        sign(indices[r,k]) *
        do.call(cbind, lapply(models[[r]]$W.Summ, function(x) x[, abs(indices[r,k])]))
      if(sign(indices[r,k])==-1){
        w.ci[w.ci$Replicate==r & w.ci$Component==k, c('Lower', 'Upper')] <-
          w.ci[w.ci$Replicate==r & w.ci$Component==k, c('Upper', 'Lower')]
      }
    }
  }
  # Henry edited 2018-08-27: change "w.ci$Lower==0 & w.ci$Upper==0" to "|"
  w.ci$contain.0 <- (w.ci$Lower * w.ci$Upper < 0 | w.ci$Lower==0 | w.ci$Upper==0)*1

  # compute medians of posterior medians and credible interval limits across replicates
  w.ci.med <- df.base[rep(seq_len(nrow(df.base)), Krobust),]
  w.ci.med$Component <- rep(1:Krobust, each=D)
  w.ci.med <- w.ci.med[, c('Component', 'Block', 'Variable', 'var.lab', 'var.order')]
  tmp <- aggregate(Lower  ~ Component + var.order, median, data=w.ci)
  w.ci.med <- merge(w.ci.med, tmp, by=c('Component', 'var.order'))
  tmp <- aggregate(Median ~ Component + var.order, median, data=w.ci)
  w.ci.med <- merge(w.ci.med, tmp, by=c('Component', 'var.order'))
  tmp <- aggregate(Upper  ~ Component + var.order, median, data=w.ci)
  w.ci.med <- merge(w.ci.med, tmp, by=c('Component', 'var.order'))

  w.ci.med <- w.ci.med[with(w.ci.med, order(w.ci.med$Component, w.ci.med$var.order)),
                       c('Component', 'Block', 'Variable', 'var.lab',
                         'var.order', 'Lower', 'Median', 'Upper')]
  # Henry edited 2018-08-27: change "w.ci$Lower==0 & w.ci$Upper==0" to "|"
  w.ci.med$contain.0 <- (w.ci.med$Lowe * w.ci.med$Upper<0 |
                           w.ci.med$Lower==0 | w.ci.med$Upper==0)
  w.ci.med$all.0 <- (w.ci.med$Lower==0 & w.ci.med$Upper==0)*1

  # compute median (across replicates) of poseterior medians if the credible interval excludes 0
  w.med <- matrix(w.ci.med$Median*(1-w.ci.med$contain.0), D, Krobust)
  colnames(w.med) <- 1:Krobust
  rownames(w.med) <- var.labs # rownames(models[[1]]$W)

  # compute the medians (across replicates) of posterior medians of factor scores
  x.rep <- array(NA, dim=c(N, Krobust, n.reps))
  for (r in 1:n.reps){
    x.rep[,,r] <- models[[r]]$X.Summ$p50[, abs(indices[r,])]
    x.rep[,,r] <- sweep(x.rep[,,r], MARGIN=2, sign(indices[r,]), '*')
  }
  x.rob <- apply(x.rep, 1:2, median, na.rm=T)
  colnames(x.rob) <- paste0('K', 1:Krobust)
  return(list(w.ci=w.ci, w.ci.med=w.ci.med, w.med=w.med, x.rob=x.rob))
}

### 4a. a base function to create a heat map
w_plot <- function(w, D, K, gr1, conf.level, replicate){
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
