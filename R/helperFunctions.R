#' Get correlation matrix
#
#' @param model A DIABLO model object.
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @importFrom stats cor

getCorMat <- function(model, method = 'pearson') {
  object <- model
  keepA <- lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)

  cord <- mapply(function(x, y, keep){
    stats::cor(x[, keep], y, use = "pairwise", method = method)
  }, x = object$X, y = object$variates[-length(object$variates)], keep = keepA[-length(keepA)])

  simMatList <- vector("list", length(model$X))
  for(i in 1:length(cord)){
    for(j in 1:length(cord)){
      simMatList[[i]][[j]] <- cord[[i]] %*% t(cord[[j]])
    }
  }
  corMat <- do.call(rbind, lapply(simMatList, function(i) do.call(cbind, i)))
}
