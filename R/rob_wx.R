#' A function to compute robust loadings and robust scores.
#'
#' @param models No description.
#' @param indices No description.
#' @param block.labs No description.
#' @param var.labs No description.
#' @export

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
      if(!is.na(indices[r,k]) & indices[r,k]!=0){ ## Henry editted on 2020-03-25 to bypass unmatched factors
        w.ci[w.ci$Replicate==r & w.ci$Component==k, c('Lower', 'Median', 'Upper')] <-
          sign(indices[r,k]) *
          do.call(cbind, lapply(models[[r]]$W.Summ, function(x) x[, abs(indices[r,k])]))
        if(sign(indices[r,k])==-1){
          w.ci[w.ci$Replicate==r & w.ci$Component==k, c('Lower', 'Upper')] <-
            w.ci[w.ci$Replicate==r & w.ci$Component==k, c('Upper', 'Lower')]
        }
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
  # w.med <- matrix(w.ci.med$Median, D, Krobust)
  w.med <- matrix(w.ci.med$Median*(1-w.ci.med$contain.0), D, Krobust)
  colnames(w.med) <- 1:Krobust
  rownames(w.med) <- var.labs # rownames(models[[1]]$W)

  # compute the medians (across replicates) of posterior medians of factor scores
  x.rep <- array(NA, dim=c(N, Krobust, n.reps))
  for (r in 1:n.reps){
    x.rep[,,r][,indices[r,]!=0 & !is.na(indices[r,])] <- models[[r]]$X.Summ$p50[, abs(indices[r,])]
    x.rep[,,r] <- sweep(x.rep[,,r], MARGIN=2, sign(indices[r,]), '*')
  }
  x.rob <- apply(x.rep, 1:2, median, na.rm=T)
  colnames(x.rob) <- paste0('K', 1:Krobust)
  return(list(w.ci=w.ci, w.ci.med=w.ci.med, w.med=w.med, x.rob=x.rob))
}
