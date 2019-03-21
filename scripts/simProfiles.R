simProfiles() <- function(nRibosomes, faFile, transAbundance) {
  ## simulate /pi_i and /rho_{ij} parameters and write to rawProfiles.txt file
  # nRibosomes: scalar; number of ribosomes per expt
  # faList: output list from readFAfile()
  # transAbundance: numeric vector; multinomial probilities for ribosomes per transcript
  
}

genTransAbundance <- function(faFile, exptLengths, exptAbundances, model="loess") {
  ## generate multinomial probabilities for ribosomes per transcript
  # faList: output list from readFAfile()
  # exptLengths: transcript lengths from model expt
  # exptAbundances: transcript abundances from model expt
  # model: model class for modelAbundances ~ modelLengths; one of c("loess", "lm")
  transLengths <- lengths(faFile) - 11 # extra 11 codons per transcript for padding
  exptModel <- do.call(model, exptAbundances ~ exptLengths) # model transcript abundance ~ length
  transAbundance <- predict(exptModel, transLengths) # use model to predict new abundances
  transAbundance <- transAbundance / sum(transAbundance) # scale to probabilities
  return(transAbundance)
}

readFAfile <- function(faFile) {
  ## read in .fa file, convert to list of vectors of codons per transcript
  # faFile: character; path to .fa file of transcript sequences
  rawFile <- readLines(faFile)
  transcriptStartLines <- grep(">", rawFile)
  transcriptNames <- sub(">", "", rawFile[transcriptStartLines])
  faList <- lapply(transcriptStartLines,
                                function(x) {
                                  sequence <- faList[x+1]
                                  nCodons <- nchar(sequence)/3
                                  codons <- substring(sequence, first=(3*(1:nCodons)-1), last=(3*(1:nCodons)+1))
                                  names(codons) <- as.character(seq.int(from=-6, length.out=nCodons))
                                })
  names(faList) <- transcriptNames
  return(faList)
}