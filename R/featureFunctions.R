# Gene Set Enrichment Analysis using
# Simple (effin') Enrichment Analysis in R (sear)
genesetEnrichment <- function(geneSet) {
  candidates <- geneSet

  result <- sear::sear(candidates)
  profile <- result[order(result$fdr),]
}

# Rename features using HGNC symbols
# Must be in following formats:
# Transciptomics: XXXX\nENSGXXXXXXXXXXX
# Proteomics: XXXXXX_XXXXXXXXXX_XXXXXX_(HGNC or NA)
# Luminex Cytokine: Premade table
#
convertHGNC <- function(model) {
  M <- model
  # Rename transciptomics
  mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  transcriptomeEnsembl <- colnames(M$X$`Transcriptomics`) %>%
    map(strsplit, "\n") %>%
    map(1) %>%
    map(last) %>%
    unlist()
  transcriptomeEnsembl <- as.matrix(transcriptomeEnsembl)
  colnames(transcriptomeEnsembl) <- "ensembl_gene_id"

  humanToHGNC <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                                filters = 'ensembl_gene_id',
                                values = transcriptomeEnsembl,
                                mart = mart)
  humanToHGNC = humanToHGNC[!duplicated(humanToHGNC$ensembl_gene_id),]
  humanToHGNC <- merge(transcriptomeEnsembl, humanToHGNC, by = "ensembl_gene_id", all = TRUE)

  unmappedTranscriptomics <- which(humanToHGNC$hgnc_symbol == "")
  for(i in unmappedTranscriptomics){
    humanToHGNC$hgnc_symbol[i] <- as.character(humanToHGNC$ensembl_gene_id[i])
  }

  colnames(M$X$`Transcriptomics`) <- humanToHGNC$hgnc_symbol
  rownames(M$loadings$`Transcriptomics`) <- humanToHGNC$hgnc_symbol

  # Rename proteomics
  names <- colnames(M$X$`Proteomics`)

  # Parse by '_'
  allProtein <- names %>%
    map(strsplit, '_') %>%
    map(1)

  # att <- biomaRt::listAttributes(maRt)
  # att[grep('uniprot', ignore.case = T, biomaRt::listAttributes(maRt)[ , 1]), ]

  # Look up hgnc symbols for protein identifiers
  table <- biomaRt::getBM(attributes = c('uniprotswissprot', 'hgnc_symbol'),
                          filters = 'uniprotswissprot',
                          values = allProtein,
                          mart = maRt)

  # Remove last string
  lookup <- allProtein %>%
    map(head, -1)

  # Join lookup and table in a dataframe
  match <- lookup %>%
    lapply(function(x) {
      df <- data.frame(uniprotswissprot = x)
      df %>%
        left_join(table)
    })


  # Create hgncSymbol vector
  hgncSymbols <- allProtein %>%
    map(last) %>%
    unlist()
  missingNames <- which(hgncSymbols == 'NA')

  # Read from chart for missing names
  missingChart <- proteinMapping

  # Replace 'NA' with hgnc symbol from biomaRt
  for (i in missingNames){
    # Find first non-NA match with biomaRt
    foo <- as.data.frame(match[i])['hgnc_symbol']
    nonNAindex <- which(!is.na(foo$hgnc_symbol))
    firstNonNA <- min(nonNAindex)
    # If non-NA, use the name given by biomaRt
    if (!(is.na(foo$hgnc_symbol[firstNonNA]))){
      hgncSymbols[i] <- foo$hgnc_symbol[firstNonNA]
    }
    # If NA, then read from .csv file where names are manually searched and inputted
    else if (as.character(lookup[[i]][1]) %in% missingChart$Data.Names){
      rowIndex <- which(proteinMapping$Data.Names == as.character(lookup[[i]][1]))
      if (proteinMapping$Gene.Symbols[rowIndex] != 'N/A'){
        hgncSymbols[i] <- proteinMapping$Gene.Symbols[rowIndex]
      }
      # If cannot find any name, just return identifier
      else{
        hgncSymbols[i] <- as.character(lookup[i])[1]
      }
    }
  }

  colnames(M$X$`Proteomics`) <- hgncSymbols
  rownames(M$loadings$`Proteomics`) <- hgncSymbols

  ## List of renamed proteins in vector hgncSymbols ##
  # legacy

  # proteomicsNames1 <- read_csv("~/HLIPROOF/nameTable/proteinNames.csv")
  # proteomicsNames2 <- read_csv("~/HLIPROOF/nameTable/proteomics.csv")
  # proteomicsNames <- rbind(proteomicsNames1,proteomicsNames2)
  # unmappedProteomics <- which(proteomicsNames$`Gene Symbols` == "N/A")
  # for(i in unmappedProteomics){
  #   proteomicsNames$`Gene Symbols`[i] <- as.character(proteomicsNames$`Data Names`[i])
  # }
  #
  # proteinNames <- colnames(M$X$`Proteomics`)
  # for(i in 1:length(proteinNames)){
  #   if(proteinNames[i] %in% proteomicsNames$`Data Names`){
  #     index <- match(proteinNames[i], proteomicsNames$`Data Names`)
  #     proteinNames[i] <- proteomicsNames$`Gene Symbols`[index]
  #   }
  # }

  ## Bug w/ circos
  # colnames(M$X$`Proteomics`) <- proteinNames
  # rownames(M$loadings$`Proteomics`) <- proteinNames

  # Rename luminex cytokines
  luminexCytokineNames <- luminexcytokineMapping

  luminexNames <- colnames(M$X$`Luminex cytokine`)
  for(i in 1:length(luminexNames)){
    if(luminexNames[i] %in% luminexCytokineNames$`Data Names`){
      index <- match(luminexNames[i], luminexCytokineNames$`Data Names`)
      luminexNames[i] <- luminexCytokineNames$`Gene Symbols`[index]
    }
  }

  colnames(M$X$`Luminex cytokine`) <- luminexNames
  rownames(M$loadings$`Luminex cytokine`) <- luminexNames

  M
}

# If there is a match, returns the index of the match in the PPI data in a vector of length 'edges'
# 'edges' and 'data' both two column data frames
matchPPI <- function(edges, data){
  matches <- prodlim::row.match(edges, data)
}
