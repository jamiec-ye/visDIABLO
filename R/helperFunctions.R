#' Get correlation matrix
#
#' @param model A DIABLO model object.

getCorMat <- function(model, method = 'pearson') {
  object <- model
  keepA <- lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)
  # names(keepA$Proteomics) <- gsub('^.+_(.+)$', '\\1', names(keepA$Proteomics))
  # names(keepA$Transcriptomics) <- gsub('^(.+)?\nENSG.+$', '\\1', names(keepA$Transcriptomics))

  cord <- mapply(function(x, y, keep){
    cor(x[, keep], y, use = "pairwise", method = method)
  }, x = object$X, y = object$variates[-length(object$variates)], keep = keepA[-length(keepA)])

  simMatList <- vector("list", length(model$X))
  for(i in 1:length(cord)){
    for(j in 1:length(cord)){
      simMatList[[i]][[j]] <- cord[[i]] %*% t(cord[[j]])
    }
  }
  corMat <- do.call(rbind, lapply(simMatList, function(i) do.call(cbind, i)))
}
