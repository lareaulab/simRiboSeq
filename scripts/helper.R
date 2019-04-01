readFAfile <- function(faFile, pad5, pad3) {
  ## read in .fa file, convert to list of vectors of codons per transcript
  # faFile: character; path to .fa file of transcript sequences
  # pad5: numeric; number of nt padding 5' end of cds
  # pad3: numeric; number of nt padding 3' end of cds
  rawFile <- readLines(faFile)
  transcriptStartLines <- grep(">", rawFile)
  transcriptNames <- sub(">", "", rawFile[transcriptStartLines])
  faList <- lapply(transcriptStartLines,
                   function(x) {
                     sequence <- rawFile[x+1]
                     nCodons <- floor((nchar(sequence) - (pad5 %% 3) - (pad3 %% 3))/3)
                     offset <- pad5 %% 3
                     codonSequence <- substring(sequence, first=(3*(1:nCodons-1)+1)+offset, last=(3*(1:nCodons))+offset)
                     names(codonSequence) <- as.character(seq.int(from=-floor(pad5/3), length.out=nCodons))
                     return(codonSequence)
                   })
  names(faList) <- transcriptNames
  return(faList)
}

readRawProfiles <- function(inputFile) {
  ## read in rawProfiles.txt file, convert to list of vectors of ribosome counts per codon per transcript
  
}