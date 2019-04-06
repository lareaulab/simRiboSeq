simTranscriptome <- function(transcriptLengths, codonDist) {
  ## simulate transcript sequences and write to .fa file
  # transcriptLengths: numeric vector; number of codons per transcript
  # codonDist: named numeric vector; genomic codon proportions (must add to 1)
  nTranscripts <- length(transcriptLengths)
  transcripts <- lapply(transcriptLengths+12, genSequence, codonDist) # extra codon padding for footprint + digestion
  transcripts <- sapply(transcripts, paste, collapse="")
  names(transcripts) <- paste0("gene", 1:nTranscripts)
  return(transcripts)
}

genSequence <- function(nCodons, codonDist) {
  ## generate artifical transcript sequence
  # nCodons: scalar; number of codons in transcript
  # codonDist: named numeric vector; genomic codon proportions (must add to 1)
  codonDist <- codonDist / sum(codonDist)
  sample(names(codonDist), size=nCodons, replace=T, prob=codonDist)
}