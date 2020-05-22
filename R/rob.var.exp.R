#' A function to compute \% variance explained
#'
#' @param block.names No description.
#' @param varIdx.by.block No description.
#' @param by.block No description.
#' @inheritParams w.signs
#' @export

rob.var.exp <- function(models, indices, block.names, varIdx.by.block, use.unmatched=F, by.block=T){
  n.reps <- length(models)
  indices <- w.signs(models, indices, use.unmatched)

  if (use.unmatched){
    K.rob <- sum(colMeans(indices==0)<1)
  } else {
    K.rob <- sum(colMeans(is.na(indices))<1)
  }

  # hack for dealing with K.rob == 1
  if (K.rob != 1) {
    W.p50.rep <- array(NA, dim=c(nrow(models[[1]]$W.Summ$p50), K.rob, n.reps))
  } else {
    W.p50.rep <- array(NA, dim=c(length(models[[1]]$W.Summ$p50), K.rob, n.reps))
  }

  for (r in 1:n.reps){
    idx_vec <- indices[r,]
    if(use.unmatched){
      idx_vec <- idx_vec[idx_vec!=0]
    }else{
      idx_vec <- idx_vec[!is.na(idx_vec)]
    }
    if (length(idx_vec)==1){
      # hack for dealing with K.rob == 1
      if (!is.null(dim(a))) {
        W.p50.rep[,,r] <- sign(idx_vec)*models[[r]]$W.Summ$p50[, abs(idx_vec)]
      } else {
        W.p50.rep[,,r] <- sign(idx_vec)*models[[r]]$W.Summ$p50
      }
    } else {
      # hack to deal with cases in which idx_vec contains zeros or NAs
      signDummy <- sign(idx_vec);
      if (any(idx_vec == 0) | any(is.na(idx_vec) == TRUE))  {

        if (!is.na(any(idx_vec == 0))) {
          zeros <- which(idx_vec == 0)
          signDummy <-signDummy[-zeros]
        } else {
          zeros <- NULL;
        }
        if (any(is.na(idx_vec) == TRUE)) {
          NAs <- which(is.na(idx_vec) == TRUE)
          signDummy <-signDummy[-NAs]
        }

        if (is.null(dim(models[[r]]$W.Summ$p50))) {
          models[[r]]$W.Summ$p50 <- as.matrix(models[[r]]$W.Summ$p50);
        }

        sweepDummy <- sweep(as.matrix(models[[r]]$W.Summ$p50[, abs(idx_vec)]),
                            MARGIN=2, signDummy, '*')
        numRows <- dim(sweepDummy)[1];
        counter = 1;
        modData = matrix(NA, nrow = numRows, ncol = length(idx_vec))
        for (COL in 1:length(idx_vec)) {
          if (any(zeros == COL)) {
            modData[,COL] <- rep(0, numRows);
          } else {
            modData[,COL] = sweepDummy[,counter];
            counter = counter + 1;
          }
        }
        W.p50.rep[,,r] <- modData;
      } else {
        W.p50.rep[,,r] <- sweep(models[[r]]$W.Summ$p50[, abs(idx_vec)],
                                MARGIN=2, signDummy, '*')
      }
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
      # hack to deal with R automatically removing dimensions in "W.p50.rep" (which messes up the apply function) when there's either only one RF or only one variable
      dummy = dim(W.p50.rep);
      Krobust = dummy[2];
      numCurrVars <- length(varIdx.by.block[[b]])

      if (length(dim(W.p50.rep[varIdx.by.block[[b]],,])) < 3) {

        dummy1 = W.p50.rep[varIdx.by.block[[b]],,]^2; #dimensions will depend on which parameter (variable or RF) has only one instance
        if (numCurrVars == 1) {

          #ve.tmp dimensions = RFs x replicates (variables = 1 and its original dimension (1) gets removed -- unless there is only one RF, in which case you need to use as.matrix and transpose...)
          if (is.null(dim(dummy1))) { #if only one RF dummy will be a vector and the dim() function will return null (could also use is.vector)
            dummy1 <- t(as.matrix(dummy1))
          }
          ve.tmp <- dummy1; #dummy 1 is already in the robust factors x replicates format that we want

        } else

          #ve.tmp dimensions = variables x replicates (RF = 1 and its original dimension (2) gets removed)
          ve.tmp <- t(as.matrix(apply(dummy1, 2, mean, na.rm=T)*100)) #take the mean across variables and get a robust factors by replicates matrix
      } else {
        ve.tmp <- apply(W.p50.rep[varIdx.by.block[[b]],,]^2, 2:3, mean, na.rm=T)*100 #take the mean across variables and get a robust factors by replicates matrix
      }
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
      # similar edit to the one above to deal with dimensionality reduction messing up the apply function
      if (length(dim(W.p50.rep[varIdx.by.block[[b]],,])) < 3) {#jd hack
        ve.tmp <- apply(dummy1, 2, sum, na.rm=T)/
          length(varIdx.by.block[[b]])
      } else {
        ve.tmp <- apply(W.p50.rep[varIdx.by.block[[b]],,]^2, 3, sum, na.rm=T)/
          length(varIdx.by.block[[b]]) #sum of entire variable x robust factor matrix for each replicate
      }

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
