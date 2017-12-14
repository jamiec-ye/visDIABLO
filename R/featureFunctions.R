#' Gene set enrichment analysis
#
#' @param geneSet A list of gene symbols.
#' @param featureMapping A list of data frames containing 'Data.Names', 'Gene.Symbols', 'Display.Names' for each of the data blocks.
#' @importFrom dplyr arrange
#' @importFrom magrittr "%>%"

genesetEnrichment <- function(geneSet, featureMapping) {
  # Get dataset
  candidates <- geneSet
  if (!is.null(featureMapping)){
    for (name in candidates){
      warning(sprintf("%s", name))
    }
  }

  # Perform 'sear'
  # TODO: remove library()
  library(sear)
  result <- sear::sear(candidates) %>%
    dplyr::arrange(fdr)
  result <- result[1:20, c("collection", "geneset", "fdr")]
  result[, -c(1,2)] <- signif(result[, -c(1,2)], digits = 3)
  result
}


#' Match with existing Protein-Protein Interaction data
#
#' If there is a match, returns the index of the match in the PPI data in a vector of length 'edges'
#' 'edges' and 'data' both two column data frames
#
#' @param edges A data frame consisting of two proteins that have an edge.
#' @param data A data frame consisting of two columns of proteins that interact
#' @importFrom prodlim row.match

matchPPI <- function(edges, data){
  matches <- prodlim::row.match(edges, data)
}
