#' A function to match factors.
#'
#' @param comps No description.
#' @param maxK No description.
#' @param corThr No description.
#' @param matchThr No description.
#' @export
#'
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
