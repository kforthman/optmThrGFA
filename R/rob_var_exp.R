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
