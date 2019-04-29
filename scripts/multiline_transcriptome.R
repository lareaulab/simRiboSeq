rm(list=ls())

faFile <- readLines("~/simRiboSeq/outputs/yeast_yeastCodons.fa")

nTranscripts <- length(faFile)/2
sequences <- faFile[2*(1:nTranscripts)]
names(sequences) <- faFile[2*(1:nTranscripts)-1]
nt_position <- 1
faLines <- lapply(1:length(sequences),
                    function(x) {
                      num_nt <- nchar(sequences[[x]])
                      starts <- seq(from=1, to = num_nt, by=60)
                      ends <- starts + 59
                      ends[length(ends)] <- num_nt
                      transcript <- sub(">", "", names(sequences)[x])
                      headerLine <- paste(names(sequences)[x],
                                          transcript,
                                          "utr5:", 3*6,
                                          "utr3:", 3*4,
                                          "Chr I", 
                                          paste0(nt_position, "-", nt_position+num_nt-3*6-3*4, ","))
                      faLines <- c(headerLine, 
                                   substring(sequences[x],
                                             first=starts,
                                             last=ends))
                      nt_position <- nt_position + num_nt
                      return(faLines)
                    })
faLines <- unlist(faLines)

writeLines(faLines, "~/iXnos/genome_data/yeast_yeastCodons_18cds12.fa")
