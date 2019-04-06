digest <- function(faList, riboCounts, delta5, delta3, minSize, maxSize) {
  ## generate footprint sequences for all transcripts (in faList and rawProfiles)
  # faList: output list from readFAfile()
  # riboCounts: output list from simProfiles()
  # delta5: named numeric vector; probabilities of 5' digest lengths
  # delta3: named numeric vector; probabilities of 3' digest lengths
  # minSize: scalar; minimum footprint size
  # maxSize: scalar: maximum footprint size
  footprints <- unlist(mapply(digest_transcript,
                       codonSequence=faList,
                       riboCounts=riboCounts,
                       delta5=delta5, delta3=delta3))
  footprints <- footprints[(nchar(footprints) >= minSize) & (nchar(footprints) <= maxSize)]
  return(footprints)
}

digest_transcript <- function(codonSequence, riboCounts, delta5, delta3) {
  ## generate footprint sequences for individual transcript
  # codonSequence: character vector; codons in transcript
  # riboCounts: numeric vector; ribosome counts per codon
  # delta5: named numeric vector; probabilities of 5' digest lengths
  # delta3: named numeric vector; probabilities of 3' digest lengths
  startCodon <- which(names(codonSequence)=="0")
  startPosition <- 3*(startCodon-1)+1
  # 1. convert vector of codons to vector of nucleotides
  ntSequence <- paste(codonSequence, collpase="")
  # 2. convert riboCounts to numeric vector of A site coordinates (in nucleotide units)
  AsitePositions <- unlist(mapply(rep, seq.int(length(riboCounts)), riboCounts))
  AsitePositions <- 3*(AsitePositions-1) + startPosition
  # 3. generate digest lengths
  digest5 <- sample(as.numeric(names(delta5)), size=length(AsitePositions), prob=delta5, replace=T)
  digest3 <- sample(as.numeric(names(delta3)), size=length(AsitePositions), prob=delta3, replace=T)
  # 4. generate footprint start and stop positions
  footprintStart <- AsitePositions - 3*5 - digest5 # 5 codons + 5' digest length
  footprintEnd <- AsitePositions + 2 + 3*3 + digest3 # 2 nt for Asite codon + 3 codons + 3' digest length
  # 5. extract footprint sequences
  footprints <- substring(ntSequence, first=footprintStart, last=footprintEnd)
  return(footprints)
}

getBiasRegion <- function(footprints, region, biasLength) {
  ## get nucleotide sequence for specified bias region
  # footprints: character vector; footprint sequences
  # region: scalar; 5 for 5' region, 3 for 3' region
  # biasLength: scalar; number of nucleotides in bias region
  footprintLengths <- sapply(footprints, nchar)
  if(region==5) {
    biasStart <- 1
    biasEnd <- biasLength
  }
  if(region==3) {
    biasEnd <- footprintLengths
    biasStart <- biasEnd - biasLength + 1
  }
  biasRegions <- unlist(mapply(substring,
                               text=footprints,
                               first=biasStart,
                               last=biasEnd))
  return(biasRegions)
}

ligate <- function(footprints, ligBias) {
  ## apply preferential 3' bias to digested footprints
  # footprints: character vector; sequences of digested footprints
  # ligBias: named numeric vector; probabilties of successful ligation for bias regions
  biasLength <- nchar(names(ligBias)[1])
  biasRegions <- getBiasRegion(footprints, region=3, biasLength=biasLength)
  ligateProbs <- ligBias[match(biasRegions, names(ligBias))]
  ligate_keep <- sapply(ligateProbs, function(p) rbinom(1, 1, p))
  footprints <- footprints[ligate_keep == 1]
  return(footprints)
}

circularize <- function(footprints, circBias) {
  ## apply preferential 5' bias to digested+ligated footprints
  # footprints: character vector; sequences of digested+ligated footprints
  # circBias: named numeric vector; probabilities of successful circularization for bias regions
  biasLength <- nchar(names(circBias)[1])
  biasRegions <- getBiasRegion(footprints, region=5, biasLength=biasLength)
  circProbs <- circBias[match(biasRegions, names(circBias))]
  circ_keep <- sapply(circProbs, function(p) rbinom(1, 1, p))
  footprints <- footprints[circ_keep == 1]
  return(footprints)
}


simFootprints <- function() {
  ## generate footprint sequences and pass through experimental procedures
  footprints <- digest()
  footprints <- ligate()
  footprints <- circularize()
  # lose footprints to sequencing
  # export footprint sequences
}
