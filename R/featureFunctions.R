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
  proteomicsNames1 <- read_csv("~/HLIPROOF/nameTable/proteinNames.csv")
  proteomicsNames2 <- read_csv("~/HLIPROOF/nameTable/proteomics.csv")
  proteomicsNames <- rbind(proteomicsNames1,proteomicsNames2)
  unmappedProteomics <- which(proteomicsNames$`Gene Symbols` == "N/A")
  for(i in unmappedProteomics){
    proteomicsNames$`Gene Symbols`[i] <- as.character(proteomicsNames$`Data Names`[i])
  }

  proteinNames <- colnames(M$X$`Proteomics`)
  for(i in 1:length(proteinNames)){
    if(proteinNames[i] %in% proteomicsNames$`Data Names`){
      index <- match(proteinNames[i], proteomicsNames$`Data Names`)
      proteinNames[i] <- proteomicsNames$`Gene Symbols`[index]
    }
  }

  ## Bug w/ circos
  # colnames(M$X$`Proteomics`) <- proteinNames
  # rownames(M$loadings$`Proteomics`) <- proteinNames

  # Rename luminex cytokines
  luminexCytokineNames <- read_csv("~/HLIPROOF/nameTable/luminexCytokineNames.csv")

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
