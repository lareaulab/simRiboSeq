simTranscriptome <- function(transcriptLengths, codonDist, outFile) {
  ## simulate transcript sequences and write to .fa file
  # transcriptLengths: numeric vector; number of codons per transcript
  # codonDist: named numeric vector; genomic codon proportions (must add to 1)
  # outFile: character; where to write output .fa file
  nTranscripts <- length(transcriptLengths)
  transcripts <- lapply(transcriptLengths+10, genSequence, codonDist) # extra codon padding for footprint + digestion
  names(transcripts) <- paste0("gene", 1:nTranscripts)
  outputFA <- rep(NULL, 2*nTranscripts)
  outputFA[2*(1:nTranscripts)-1] <- paste0(">", names(transcripts))
  outputFA[2*(1:nTranscripts)] <- sapply(transcripts, paste, collapse="")
  writeLines(outputFA, con=outFile, sep="\n")
  return(transcripts)
}

genSequence <- function(nCodons, codonDist) {
  ## generate artifical transcript sequence
  # nCodons: scalar; number of codons in transcript
  # codonDist: named numeric vector; genomic codon proportions (must add to 1)
  sample(names(codonDist), size=nCodons, replace=T, prob=codonDist)
}

### yeast genome
yeastGeneLengths <- read.table("refData/scer.transcripts.13cds10.lengths.txt", header=F)[,3]/3
yeastFAfile <- readLines("refData/scer.transcripts.13cds10.fa")
yeastGeneCodons <- lapply(grep(">", yeastFAfile),
                     function(x) {
                       sequence <- yeastFAfile[x+1]
                       nCodons <- floor(nchar(sequence)/3)
                       # account for 13cds10
                       codons <- substring(sequence, first=(3*(1:nCodons)-1), last=(3*(1:nCodons)+1))
                       # 0-based indexing
                       names(codons) <- as.character(seq.int(from=-4, length.out=nCodons))
                       return(codons)
                     })
codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")
yeastCodonCounts <- sapply(yeastGeneCodons,
                         function(gene) {
                           gene_trim <- gene[5:(length(gene)-4)]
                           sapply(codons, 
                                  function(codon) {
                                    sum(gene_trim==codon)
                                  })
                         })
yeastCodonDist <- rowSums(yeastCodonCounts)/sum(yeastCodonCounts)
randomCodonDist <- rep(1, length(codons))/length(codons)
names(randomCodonDist) <- codons
# simulate yeast transcriptome: yeast gene lengths, yeast codon usage
yeastGenes <- simTranscriptome(yeastGeneLengths, yeastCodonDist, 
                               file.path("outputs", "simYeast.fa"))
# simulate random transcriptome: yeast gene lengths, random codon usage
randomGenes <- simTranscriptome(yeastGeneLengths, randomCodonDist, 
                                file.path("outputs", "simRandom.fa"))
