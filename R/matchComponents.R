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
