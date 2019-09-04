genRawProfiles <- function(nRibosomes, rhos, pis) {
  ## generate ribosome counts per codon per transcript
  # nRibosomes: scalar; number of ribosomes per expt
  # pis: list of numeric vector; multinomial probabilities for ribosomes per codon per transcript [ output from simPi() ]
  # rhos: numeric vector; multinomial probilities for ribosomes per transcript [ output from simRho() ]
  print("... ... generating raw profile")
  riboPerTrans <- rmultinom(n=1, size=nRibosomes, prob=rhos)
  rawProfile <- mapply(rmultinom, 
                       n=1,
                       size=riboPerTrans,
                       prob=pis)
  rawProfile <- lapply(rawProfile, function(x) x[,1])
  names(rawProfile) <- names(rhos)
  return(rawProfile)
}

footprint <- setClass("footprint", slots=list(sequence="character",
                                              transcript="character",
                                              ASite="numeric",
                                              digest5="numeric",
                                              digest3="numeric",
                                              id="numeric"))

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
  Asites <- unlist(mapply(rep, seq.int(length(riboCounts)), riboCounts)) - 1
  if(length(Asites)==0) {
    return(NULL)
  }
  AsitePositions <- 3*(Asites) + startPosition
  # 3. generate digest lengths
  digest5 <- sample(as.numeric(names(delta5)), size=length(AsitePositions), prob=delta5, replace=T)
  digest3 <- sample(as.numeric(names(delta3)), size=length(AsitePositions), prob=delta3, replace=T)
  # 4. generate footprint start and stop positions
  footprintStart <- AsitePositions - digest5 # 5' digest length
  footprintEnd <- AsitePositions + 2 + digest3 # 2 nt for Asite codon + 3' digest length
  # 5. extract footprint sequences
  footprints <- substring(text=ntSequence, first=footprintStart, last=footprintEnd)
  # 6. convert to "footprint" object
  ids <- unlist(mapply(seq.int, to=length(footprints)))
  reads <- mapply(footprint, sequence=footprints, ASite=Asites, transcript=transcriptName,
                  digest5=digest5, digest3=digest3, id=ids)
  return(reads)
}

digest <- function(faList, riboCounts, delta5, delta3, minSize=27, maxSize=31, 
                   digest_transcript) {
  ## generate footprint sequences for all transcripts (in faList and rawProfiles)
  # faList: output list from readFAfile()
  # riboCounts: output list from simProfiles()
  # delta5: named numeric vector; probabilities of 5' digest lengths
  # delta3: named numeric vector; probabilities of 3' digest lengths
  # minSize: scalar; minimum footprint size
  # maxSize: scalar; maximum footprint size
  # mc.cores: scalar; number of cores to use for parallelization
  # digest_transcript: function; function that digests individual transcript
  print("... ... digesting footprints")
  nTranscripts <- length(faList)
  # digest footprints
  footprints <- unlist(mcmapply(digest_transcript,
                                codonSequence=faList,
                                riboCounts=riboCounts,
                                transcriptName=names(faList),
                                MoreArgs=list(delta5=delta5, delta3=delta3)))
  footprintLengths <- sapply(footprints, function(x) nchar(x@sequence))
  footprints <- footprints[(footprintLengths >= minSize) & (footprintLengths <= maxSize)]
  return(footprints)
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
  biasRegions <- unlist(mcmapply(substring,
                                 text=sapply(footprints, function(x) x@sequence),
                                 first=biasStart,
                                 last=biasEnd))
  return(biasRegions)
}

ligate <- function(footprints, ligBias) {
  ## apply preferential 3' bias to digested footprints
  # footprints: list of "footprint" objects
  # ligBias: named numeric vector; probabilties of successful ligation for bias regions
  print("... ... ligating footprints ...")
  biasLength <- nchar(names(ligBias)[1])
  biasRegions <- getBiasRegion(footprints, region=3, biasLength=biasLength)
  ligateProbs <- ligBias[match(biasRegions, names(ligBias))]
  ligate_keep <- sapply(ligateProbs, function(p) rbinom(1, 1, p))
  footprints <- footprints[ligate_keep == 1]
  return(footprints)
}

RT <- function(footprints, RTBias) {
  ## add extra, non-templated base at 5' end of digested footprints
  # footprints: list of "footprint" objects
  # RTBias: named numeric vector, probabilities of adding additional base
  print("... ... reverse transcribing ...")
  bases <- names(RTBias)
  basesToAdd <- sample(bases, size=length(footprints), replace=T, prob=RTBias)
  for(i in (1:length(basesToAdd))[basesToAdd != ""]) {
    footprints[[i]]@sequence <- paste0(basesToAdd[i], footprints[[i]]@sequence, collapse="")
  }
  return(footprints)
}

circularize <- function(footprints, circBias) {
  ## apply preferential 5' bias to digested+ligated footprints
  # footprints: list of "footprint" objects
  # circBias: named numeric vector; probabilities of successful circularization for bias regions
  print("... ... circularizing footprints ...")
  biasLength <- nchar(names(circBias)[1])
  biasRegions <- getBiasRegion(footprints, region=5, biasLength=biasLength)
  circProbs <- circBias[match(biasRegions, names(circBias))]
  circ_keep <- sapply(circProbs, function(p) rbinom(1, 1, p))
  footprints <- footprints[circ_keep == 1]
  return(footprints)
}

simFootprints <- function(faList, nRibosomes, rhos, pis, 
                          delta5, delta3, minSize=27, maxSize=31, 
                          mc.cores=NULL, digest_transcript,
                          ligBias, RTBias, circBias) {
  ## generate footprint sequences and pass through experimental procedures
  # faList: output list from readFAfile()
  # nRibosomes: scalar; number of ribosomes per expt
  # pis: list of numeric vector; multinomial probabilities for ribosomes per codon per transcript [ output from simPi() ]
  # rhos: numeric vector; multinomial probilities for ribosomes per transcript [ output from simRho() ]
  # delta5: named numeric vector; probabilities of 5' digest lengths
  # delta3: named numeric vector; probabilities of 3' digest lengths
  # minSize: scalar; minimum footprint size
  # maxSize: scalar: maximum footprint size
  # mc.cores: scalar; number of cores to use for parallelization for digest()
  # digest_transcript: function; function that digests individual transcript
  # ligBias: named numeric vector; probabilties of successful ligation for bias regions
  # circBias: named numeric vector; probabilities of successful circularization for bias regions
  if(is.null(mc.cores)) {mc.cores <- parallel::detectCores()-10}
  print(paste("Number of cores:", mc.cores))
  footprintCluster <- makeCluster(mc.cores)
  nRounds <- 1
  print(paste("... Round", nRounds, "of simulating footprints"))
  codonCounts <- genRawProfiles(nRibosomes, rhos, pis)
  footprints <- digest(faList, codonCounts, delta5, delta3, minSize, maxSize,
                       digest_transcript=digest_transcript)
  footprints <- ligate(footprints, ligBias)
  footprints <- RT(footprints, RTBias)
  footprints <- circularize(footprints, circBias)
  print(paste(length(footprints), "of", nRibosomes, 
              "(", round(length(footprints)/nRibosomes*100, digits=2), 
              " % ) footprints simulated"))
  while(length(footprints) < nRibosomes) {
    nRounds <- nRounds + 1
    print(paste("... Round", nRounds, "of simulating footprints"))
    tmpCodonCounts <- genRawProfiles(nRibosomes, rhos, pis)
    tmpFootprints <- digest(faList, codonCounts, delta5, delta3, minSize, maxSize,
                            digest_transcript=digest_transcript)
    tmpFootprints <- ligate(footprints, ligBias)
    tmpFootprints <- circularize(footprints, circBias)
    if(length(tmpFootprints) > (nRibosomes-length(footprints))) {
      tmpFootprints <- sample(tmpFootprints, size=nRibosomes-length(footprints))
    }
    footprints <- append(footprints, tmpFootprints)
    print(paste("...", length(footprints), "of", nRibosomes, 
                "(", round(length(footprints)/nRibosomes*100, digits=2), 
                " % ) footprints simulated"))
  }
  stopCluster(footprintCluster)
  return(footprints)
}
