# HELPER cormat
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

# HELPER plotting
ggVarGrid <- function(plots, ...) {
  # extract the legend from one of the plots
  legend <- get_legend(plots[[1]])
  # add the legend to the row we made earlier
  pgrid <- plot_grid(plotlist = map(plots, ~ . + theme(legend.position = 'none')),
                     label_size = 20, labels = 'AUTO', ...)
  plot_grid(pgrid, legend, rel_widths = c(10, 1))
}

# HELPER pseudo ellipses
ggally_ellipse <- function (data, mapping, ...) {
  rangeX <- range(GGally:::eval_data_col(data, mapping$x), na.rm = TRUE)
  rangeY <- range(GGally:::eval_data_col(data, mapping$y), na.rm = TRUE)
  p <- ggplot(data = data) +
    geom_point(data = data.frame(rangeX = rangeX, rangeY = rangeY),
               mapping = aes(x = rangeX, y = rangeY), alpha = 0) +
    geom_point(mapping = mapping, ...) +
    stat_ellipse(mapping = mapping, ...)
  p
}

# HELPER unit circle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# HELPER rescale
rescale <- function(x, scalar = 1) {
  y <- sign(x)
  x <- abs(x)
  x <- scalar * (x - min(x))/(max(x) - min(x))
  x * y
}

#' create the plot given a standard input of 2 data.frames
#' 1. vars - a data.frame containing the variates
#' 2. vecs (optional) - a data.frame of the loadings
# ggProducePlot <- function(vars, labs, vecs = NULL) {
#   if(!is.null(vecs))
#     vars <- vars %>% mutate(`comp 1` = rescale(`comp 1`, 0.5), `comp 2` = rescale(`comp 2`, 0.5))
#
#   g <- ggplot(vars) +
#     scale_color_canva('Time', palette = 'Warm and cool') +
#     labs(x = sprintf('Component 1\n(%2.1f%% var. explained)', labs[1] * 100),
#          y = sprintf('Component 2\n(%2.1f%% var. explained)', labs[2] * 100))
#
#   if(!is.null(vecs)) {
#     vecs <- vecs %>%
#       mutate(sel = max(abs(x), abs(y))) %>%
#       arrange(desc(sel)) %>%
#       slice(1:10)
#
#     g <- g +
#       geom_vline(xintercept = 0, alpha = 0.2) +
#       geom_hline(yintercept = 0, alpha = 0.2) +
#       geom_path(data = circleFun(diameter = 1), aes(x, y), alpha = 0.2) +
#       geom_path(data = circleFun(diameter = 2), aes(x, y), alpha = 0.2) +
#       geom_text_repel(data = vecs, aes(x, y, label = names), size = 3) +
#       geom_segment(data = vecs, aes(x = 0, y = 0, xend = x, yend = y),
#                    alpha = 0.5,
#                    arrow = arrow(type = 'open', length = unit(0.25,"cm")))
#   }
#   g <- g +
#     geom_point(aes(`comp 1`, `comp 2`, colour = day), size = 3) +
#     stat_ellipse(aes(`comp 1`, `comp 2`, colour = day), level = 0.95) +
#     theme(axis.ticks = element_blank(), axis.text = element_blank())
# }

ggProducePlot <- function(vars, labs, vecs = NULL, topn = 10) {
  if(!is.null(vecs)) {
    vecs <- vecs %>%
      mutate(sel = max(abs(x), abs(y))) %>%
      arrange(desc(sel)) %>%
      slice(1:topn)
    list(component_plot = ggBindVars(vars, labs), variable_plot  = ggBindVecs(vecs, labs))
  }
  else
    ggBindVars(vars, labs)
}

ggBindVars <- function(vars, labs) {
  ggplot(vars) +
    geom_point(aes(`comp 1`, `comp 2`, colour = day), size = 3) +
    # stat_ellipse(aes(`comp 1`, `comp 2`, colour = day), level = 0.95) +
    scale_color_canva('Time', palette = 'Warm and cool') +
    labs(x = sprintf('Component 1\n(%2.1f%% var. explained)', labs[1] * 100),
         y = sprintf('Component 2\n(%2.1f%% var. explained)', labs[2] * 100))
}

ggBindVecs <- function(input, labs, vecs = NULL) UseMethod("ggBindVecs")

ggBindVecs.data.frame <- function(input, labs) {
  ggplot(input) +
    geom_vline(xintercept = 0, alpha = 0.2) +
    geom_hline(yintercept = 0, alpha = 0.2) +
    geom_path(data = circleFun(diameter = 1), aes(x, y), alpha = 0.2) +
    geom_path(data = circleFun(diameter = 2), aes(x, y), alpha = 0.2) +
    geom_text_repel(aes(x, y, label = names), size = 3) +
    geom_segment(aes(x = 0, y = 0, xend = x, yend = y),
                 alpha = 0.5,
                 arrow = arrow(type = 'open', length = unit(0.25,"cm"))) +
    labs(x = sprintf('Component 1\n(%2.1f%% var. explained)', labs[1] * 100),
         y = sprintf('Component 2\n(%2.1f%% var. explained)', labs[2] * 100))
}

ggBindVecs.ggplot <- function(input, labs, vecs) {
  input +
    geom_vline(xintercept = 0, alpha = 0.2) +
    geom_hline(yintercept = 0, alpha = 0.2) +
    geom_path(data = circleFun(diameter = 1), aes(x, y), alpha = 0.2) +
    geom_path(data = circleFun(diameter = 2), aes(x, y), alpha = 0.2) +
    geom_text_repel(data = vecs, aes(x, y, label = names), size = 3) +
    geom_segment(data = vecs, aes(x = 0, y = 0, xend = x, yend = y),
                 alpha = 0.5,
                 arrow = arrow(type = 'open', length = unit(0.25,"cm"))) +
    labs(x = sprintf('Component 1\n(%2.1f%% var. explained)', labs[1] * 100),
         y = sprintf('Component 2\n(%2.1f%% var. explained)', labs[2] * 100))
}






ggCompPlot <- function(mix, meta) UseMethod("ggCompPlot")

ggCompPlot.pca <- function(mix, meta) {
  vars <- cbind(mix$x[ , 1:2], meta)
  colnames(vars) <- c('comp 1', 'comp 2', colnames(meta))
  labs <- mix$explained_variance
  ggProducePlot(vars, labs)
}

ggCompPlot.spca <- function(mix, meta) {
  vars <- cbind(mix$x[ , 1:2], meta)
  colnames(vars) <- c('comp 1', 'comp 2', colnames(meta))
  labs <- mix$explained_variance
  vecs <- mixOmics::plotVar(mix, plot = F)
  ggProducePlot(vars, labs, vecs)
}

ggCompPlot.sipca <- function(mix, meta) {
  vars <- cbind(mix$x[ , 1:2], meta)
  colnames(vars) <- c('comp 1', 'comp 2', colnames(meta))
  labs <- mix$explained_variance
  vecs <- mixOmics::plotVar(mix, plot = F)
  ggProducePlot(vars, labs, vecs)
}

ggCompPlot.splsda <- function(mix, meta) {
  vars <- cbind(mix$variates$X[ , 1:2], meta)
  colnames(vars) <- c('comp 1', 'comp 2', colnames(meta))
  labs <- mix$explained_variance$X
  vecs <- mixOmics::plotVar(mix, plot = F)
  ggProducePlot(vars, labs, vecs)
}

ggCompPlot.block.splsda <- function(mix, meta) {
  pmap(.l = list(vars = head(mix$variates, -1),
                 labs = head(mix$explained_variance, -1),
                 vecs = mixOmics::plotVar(mix, plot = F) %>% split(.$Block)),
       function(vars, labs, vecs) {
         vars <- cbind(vars[ , 1:2], meta)
         colnames(vars) <- c('comp 1', 'comp 2', colnames(meta))
         ggProducePlot(vars, labs, vecs)
       })
}
