simPi <- function(faList, exptLengths, exptAbundances, model="loess") {
  ## generate per-transcript probabilities and write to [[xx]]
  # faList: output list from readFAfile()
  # exptLengths: transcript lengths from a model experiment
  # exptAbundances: transcript abundances from a model experiment
  # model: model class for modelAbundances ~ modelLengths; one of c("loess", "lm")
  transLengths <- lengths(faList) - 11 # extra 11 codons per transcript for padding
  exptModel <- do.call(model, exptAbundances ~ exptLengths) # model transcript abundance ~ length
  transAbundance <- predict(exptModel, transLengths) # use model to predict new abundances
  pis <- transAbundance / sum(transAbundance) # scale to probabilities, sum to 1
  names(pis) <- names(faList)
  return(pis)
}

simRho <- function(faList, pad5, pad3, codonTE) {
  ## generate per-codon probabilities for all transcripts and write to [[xx]]
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

genRawProfiles() <- function(nRibosomes, pis, rhos) {
  ## generate ribosome counts per codon per transcript
  # nRibosomes: scalar; number of ribosomes per expt
  # pis: numeric vector; multinomial probilities for ribosomes per transcript [ output from simPi() ]
  # rhos: list of numeric vector; multinomial probabilities for ribosomes per codon per transcript [ output from simRho() ]
  riboPerTrans <- rmultinom(n=1, size=nRibosomes, prob=pis)
  rawProfile <- mapply(rmultinom, 
                       n=1,
                       size=riboPerTrans,
                       prob=rhos)
  names(rawProfile) <- names(pis)
  return(rawProfile)
}

simProfiles <- function() {
  ## simulate /pi_i and /rho_{ij} parameters and write to rawProfiles.txt file
  pis <- simPis()
  rhos <- simRhos()
  rawProfile <- genRawProfiles()
  # export parameters and rawProfile.txt
}

### yeast genome
source("scripts/helper.R")
yeastFAfile <- "refData/scer.transcripts.13cds10.fa"
yeastGeneCodons <- readFAfile(yeastFAfile, pad5=13, pad3=10)

