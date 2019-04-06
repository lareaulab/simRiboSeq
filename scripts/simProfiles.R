simPi <- function(faList, exptLengths, exptAbundances) {
  ## generate per-transcript probabilities
  # faList: output list from readFAfile()
  # exptLengths: transcript lengths from a model experiment
  # exptAbundances: transcript abundances from a model experiment
  transLengths <- lengths(faList) - 11 # extra 11 codons per transcript for padding
  exptModel <- loess(exptAbundances ~ exptLengths) # model transcript abundance ~ length
  transAbundance <- predict(exptModel, transLengths) # use model to predict new abundances
  transAbundance[is.na(transAbundance)] <- min(transAbundance, na.rm=T)/2 # abundance for transcripts outside model domain
  transAbundance <- transAbundance + rpois(length(faList), 10) - 7 # add noise
  pis <- transAbundance / sum(transAbundance) # scale to probabilities, sum to 1
  names(pis) <- names(faList)
  return(pis)
}

simRho <- function(faList, pad5, pad3, codonTE) {
  ## generate per-codon probabilities for all transcripts
  # faList: output list from readFAfile()
  # pad5: scalar; number of codons padding 5' end of cds
  # pad3: scalar; number of codons padding 3' end of cds
  # codonTE: named numeric vector; translational efficiency of codon averaged over genome
  transCodons <- lapply(faList, 
                        function(transcript) {
                          transcript[(pad5+1):(length(transcript)-pad3)]
                        })
  rhos <- lapply(transCodons,
                 function(transcript) {
                   codonRhos <- codonTE[match(transcript, names(codonTE))] # match per-codon TE to genomic average
                   codonRhos <- codonRhos / sum(codonRhos) # scale to probabilities, sum to 1
                   return(codonRhos)
                 })
  names(rhos) <- names(faList)
  return(rhos)
}