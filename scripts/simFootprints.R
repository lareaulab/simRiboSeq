genRawProfiles <- function(nRibosomes, pis, rhos) {
  ## generate ribosome counts per codon per transcript
  # nRibosomes: scalar; number of ribosomes per expt
  # pis: numeric vector; multinomial probilities for ribosomes per transcript [ output from simPi() ]
  # rhos: list of numeric vector; multinomial probabilities for ribosomes per codon per transcript [ output from simRho() ]
  print("... generating raw profile ...")
  riboPerTrans <- rmultinom(n=1, size=nRibosomes, prob=pis)
  rawProfile <- mapply(rmultinom, 
                       n=1,
                       size=riboPerTrans,
                       prob=rhos)
  names(rawProfile) <- names(pis)
  return(rawProfile)
}

footprint <- setClass("footprint", slots=list(sequence="character",
                                              transcript="character",
                                              ASite="numeric"))

digest <- function(faList, riboCounts, delta5, delta3, minSize, maxSize) {
  ## generate footprint sequences for all transcripts (in faList and rawProfiles)
  # faList: output list from readFAfile()
  # riboCounts: output list from simProfiles()
  # delta5: named numeric vector; probabilities of 5' digest lengths
  # delta3: named numeric vector; probabilities of 3' digest lengths
  # minSize: scalar; minimum footprint size
  # maxSize: scalar: maximum footprint size
  print("... digesting footprints ...")
  footprints <- unlist(mapply(digest_transcript,
                              codonSequence=faList,
                              riboCounts=riboCounts,
                              transcriptName=names(faList),
                              MoreArgs=list(delta5=delta5, delta3=delta3)))
  footprintLengths <- sapply(footprints, function(x) nchar(x@sequence))
  footprints <- footprints[(footprintLengths >= minSize) & (footprintLengths <= maxSize)]
  return(footprints)
}

digest_transcript <- function(codonSequence, riboCounts, transcriptName, delta5, delta3) {
  ## generate footprint sequences for individual transcript
  # codonSequence: character vector; codons in transcript
  # riboCounts: numeric vector; ribosome counts per codon
  # transcriptName: character; name of transcript
  # delta5: named numeric vector; probabilities of 5' digest lengths
  # delta3: named numeric vector; probabilities of 3' digest lengths
  startCodon <- which(names(codonSequence)=="0")
  startPosition <- 3*(startCodon-1)+1
  # 1. convert vector of codons to vector of nucleotides
  ntSequence <- paste(codonSequence, collapse="")
  # 2. convert riboCounts to numeric vector of A site coordinates (in nucleotide units)
  Asites <- unlist(mapply(rep, seq.int(length(riboCounts)), riboCounts))
  AsitePositions <- 3*(Asites-1) + startPosition
  # 3. generate digest lengths
  digest5 <- sample(as.numeric(names(delta5)), size=length(AsitePositions), prob=delta5, replace=T)
  digest3 <- sample(as.numeric(names(delta3)), size=length(AsitePositions), prob=delta3, replace=T)
  # 4. generate footprint start and stop positions
  footprintStart <- AsitePositions - 3*5 - digest5 # 5 codons + 5' digest length
  footprintStart[footprintStart < 1] <- 1
  footprintEnd <- AsitePositions + 2 + 3*3 + digest3 # 2 nt for Asite codon + 3 codons + 3' digest length
  footprintEnd[footprintEnd > nchar(ntSequence)] <- nchar(ntSequence)
  # print(paste(transcriptName, nchar(ntSequence), "nt"))
  # print(range(footprintStart))
  # print(range(footprintEnd))
  # 5. extract footprint sequences
  footprints <- substring(ntSequence, first=footprintStart, last=footprintEnd)
  reads <- mapply(new, sequence=footprints, ASite=Asites, 
                  MoreArgs=list(Class="footprint", transcript=transcriptName))
  return(reads)
}

getBiasRegion <- function(footprints, region, biasLength) {
  ## get nucleotide sequence for specified bias region
  # footprints: list of "footprint" objects
  # region: scalar; 5 for 5' region, 3 for 3' region
  # biasLength: scalar; number of nucleotides in bias region
  footprintLengths <- sapply(footprints, function(x) nchar(x@sequence))
  if(region==5) {
    biasStart <- 1
    biasEnd <- biasLength
  }
  if(region==3) {
    biasEnd <- footprintLengths
    biasStart <- biasEnd - biasLength + 1
  }
  biasRegions <- unlist(mapply(substring,
                               text=sapply(footprints, function(x) x@sequence),
                               first=biasStart,
                               last=biasEnd))
  return(biasRegions)
}

biasScores <- read.table("~/simRiboSeq/refData/codon_scores.tsv")
rownames(biasScores) <- sort(codons)
p5bias <- biasScores[,1]
names(p5bias) <- sort(codons)
p5bias <- (p5bias+1)/(max(p5bias, na.rm=T)+1)
p5bias[is.na(p5bias)] <- 0
n3bias <- biasScores[,9]
names(n3bias) <- sort(codons)
n3bias <- (n3bias+1)/(max(n3bias, na.rm=T)+1)
n3bias[is.na(n3bias)] <- 0


ligate <- function(footprints, ligBias) {
  ## apply preferential 3' bias to digested footprints
  # footprints: character vector; sequences of digested footprints
  # ligBias: named numeric vector; probabilties of successful ligation for bias regions
  print("... ligating footprints ...")
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
  print("... circularizing footprints ...")
  biasLength <- nchar(names(circBias)[1])
  biasRegions <- getBiasRegion(footprints, region=5, biasLength=biasLength)
  circProbs <- circBias[match(biasRegions, names(circBias))]
  circ_keep <- sapply(circProbs, function(p) rbinom(1, 1, p))
  footprints <- footprints[circ_keep == 1]
  return(footprints)
}

simFootprints <- function(faList, nRibosomes, pis, rhos, 
                          delta5, delta3, minSize, maxSize, 
                          ligBias, circBias) {
  ## generate footprint sequences and pass through experimental procedures
  # faList: output list from readFAfile()
  # nRibosomes: scalar; number of ribosomes per expt
  # pis: numeric vector; multinomial probilities for ribosomes per transcript [ output from simPi() ]
  # rhos: list of numeric vector; multinomial probabilities for ribosomes per codon per transcript [ output from simRho() ]
  # delta5: named numeric vector; probabilities of 5' digest lengths
  # delta3: named numeric vector; probabilities of 3' digest lengths
  # minSize: scalar; minimum footprint size
  # maxSize: scalar: maximum footprint size
  # ligBias: named numeric vector; probabilties of successful ligation for bias regions
  # circBias: named numeric vector; probabilities of successful circularization for bias regions
  nRounds <- 1
  print(paste("Round", nRounds, "of simulating footprints"))
  codonCounts <- genRawProfiles(nRibosomes, pis, rhos)
  footprints <- digest(faList, codonCounts, delta5, delta3, minSize, maxSize)
  footprints <- ligate(footprints, ligBias)
  footprints <- circularize(footprints, circBias)
  print(paste(length(footprints), "of", nRibosomes, 
              "(", round(length(footprints)/nRibosomes*100, digits=2), 
              " % ) footprints simulated"))
  while(length(footprints) < nRibosomes) {
    nRounds <- nRounds + 1
    print(paste("Round", nRounds, "of simulating footprints"))
    tmpCodonCounts <- genRawProfiles(nRibosomes, pis, rhos)
    tmpFootprints <- digest(faList, tmpCodonCounts, delta5, delta3, minSize, maxSize)
    tmpFootprints <- ligate(footprints, ligBias)
    tmpFootprints <- circularize(footprints, circBias)
    if(length(tmpFootprints) > (nRibosomes-length(footprints))) {
      tmpFootprints <- sample(tmpFootprints, size=nRibosomes-length(footprints))
    }
    footprints <- append(footprints, tmpFootprints)
    print(paste(length(footprints), "of", nRibosomes, 
                "(", round(length(footprints)/nRibosomes*100, digits=2), 
                " % ) footprints simulated"))
  }
  return(footprints)
}