readFAfile <- function(faFile, pad5, pad3) {
  ## read in .fa file where sequence broken up over multiple lines, convert to list of vectors of codons per transcript
  # faFile: character; path to .fa file of transcript sequences
  # pad5: numeric; number of nt padding 5' end of cds
  # pad3: numeric; number of nt padding 3' end of cds
  rawFile <- readLines(faFile)
  transcriptStartLines <- grep(">", rawFile)
  nTranscripts <- length(transcriptStartLines)
  transcriptNames <- sapply(transcriptStartLines,
                            function(x) {
                              gsub(">", "", strsplit(rawFile[x], split=" ")[[1]][1])
                            })
  transcriptStartLines <- c(transcriptStartLines, length(rawFile)+1)
  transcriptStartLines <- c(transcriptStartLines, length(rawFile)+1) # add extra line for bookkeeping
  faList <- sapply(1:nTranscripts,
                   function(x) {
                     startLine <- transcriptStartLines[x]+1
                     endLine <- transcriptStartLines[x+1]-1
                     transcriptSequence <- paste(rawFile[startLine:endLine], collapse="")
                     nCodons <- floor((nchar(transcriptSequence))/3)
                     sequenceOffset <- pad5 %% 3
                     codonSequence <- substring(transcriptSequence, 
                                                first=(3*(1:nCodons-1)+1)+sequenceOffset, 
                                                last=(3*(1:nCodons))+sequenceOffset)
                     names(codonSequence) <- as.character(seq.int(from=-floor(pad5/3), length.out=nCodons))
                     return(codonSequence)
                   })
  names(faList) <- transcriptNames
  return(faList)
}

countCodons <- function(faList, codons, pad5, pad3) {
  ## count up codons across transcriptome
  # faList: output list from readFAfile()
  # codons: character vector; codons in transcriptome
  # pad5: numeric; number of nt padding 5' end of cds
  # pad3: numeric; number of nt padding 3' end of cds
  pad5 <- floor(pad5/3)
  pad3 <- floor(pad3/3)
  codonCounts <- sapply(faList,
                        function(gene) {
                          gene_trim <- gene[(pad5+1):(length(gene)-pad3-1)]
                          sapply(codons,
                                 function(codon) {
                                   sum(gene_trim==codon)
                                 })
                        })
  return(codonCounts)
}

readRawProfiles <- function(inputFile) {
  ## read in rawProfiles.txt file, convert to list of vectors of ribosome counts per codon per transcript
  rawFile <- readLines(inputFile)
  transcriptNames <- sapply(rawFile,
                            function(x) {
                              strsplit(x, split="\t")[[1]][1]
                            })
  codonCounts <- lapply(rawFile,
                        function(x) {
                          counts <- strsplit(x, split="\t")[[1]]
                          nCodons <- length(counts)-1
                          counts <- as.numeric(counts[1:nCodons+1])
                          return(counts)
                        })
  names(codonCounts) <- transcriptNames
  return(codonCounts)
}

writeTranscriptomeFA <- function(transcripts, outFile) {
  ## write output from simTranscriptome() to .fa file
  # transcripts: character vector; transcript sequences
  # outFile: filename for output .fa file
  nTranscripts <- length(transcripts)
  transcripts <- sapply(transcripts, paste, collapse="")
  outputFA <- rep(NULL, 2*nTranscripts)
  outputFA[2*(1:nTranscripts)-1] <- paste0(">", names(transcripts))
  outputFA[2*(1:nTranscripts)] <- transcripts
  writeLines(outputFA, con=outFile, sep="\n")
}

writeFootprintsFA <- function(footprints, outFile) {
  ## write output from simFootprints() to .fa file
  # footprints: list of "footprints" objects
  # outFile: character; filename for output .fa file
  nFootprints <- length(footprints)
  outputFA <- rep(NULL, 2*nFootprints)
  footprintNames <- sapply(footprints, 
                           function(x) paste(x@transcript, x@ASite, x@digest5, x@digest3, 
                                             x@sequence, x@id, sep="_"))
  footprintSequences <- sapply(footprints, function(x) x@sequence)
  outputFA[2*(1:nFootprints)-1] <- paste0(">", footprintNames)
  outputFA[2*(1:nFootprints)] <- footprintSequences
  writeLines(outputFA, con=outFile, sep="\n")
}

readCtsByCodon <- function(filename) {
  ## read in cts_by_codon file (output by iXnos)
  # filename: character; path to cts_by_codon file
  rawFile <- readLines(filename)
  output <- lapply(rawFile,
                   function(x) {
                     txt <- strsplit(x, split="\t")[[1]]
                     return(txt[2:length(txt)]) # return only counts
                   })
  names(output) <- sapply(rawFile,
                          function(x) {
                            strsplit(x, split="\t")[[1]][1] # pull transcript names
                          })
  return(output)
}
