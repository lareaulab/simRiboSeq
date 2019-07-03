rm(list=ls())

outputDir <- "~/simRiboSeq/outputs/"

expts <- c("yeast_uniformCodons_uniformDelta_noBias",
           "yeast_uniformCodons_uniformDelta_withBias",
           "yeast_yeastCodons_uniformDelta_noBias",
           "yeast_yeastCodons_uniformDelta_withBias")

for(expt in expts) {
  faFiles <- grep("Rda", list.files(file.path(outputDir, expt)), value=T)
  nFootprints <- 0
  outFile <- file.path(outputDir, paste0(expt, ".fa"))
  for(part in faFiles) {
    print(paste("Loading", part))
    load(file=file.path(outputDir, expt, part))
    objName <- sub("Codons", "", sub("Delta", "", sub("80Mreads_", "", sub(".Rda", "", part))))
    nNewFootprints <- length(get(objName))
    outputFA <- rep(NULL, 2*nNewFootprints)
    footprintNames <- sapply(get(objName), 
                             function(x) {
                               paste(x@transcript, x@ASite, 
                                     "d5", x@digest5, "d3", x@digest3, 
                                     sep="_")
                             })
    footprintNames <- paste(footprintNames, 
                            seq.int(from=nFootprints+1, length.out=nNewFootprints),
                            sep="_")
    footprintSequences <- sapply(get(objName), function(x) x@sequence)
    outputFA[2*(1:nNewFootprints)-1] <- paste0(">", footprintNames)
    outputFA[2*(1:nNewFootprints)] <- footprintSequences
    cat(outputFA, file=outFile, sep="\n", append=T)
    rm(list=objName)
    nFootprints <- nFootprints + nNewFootprints
  }
}

quit(save="no")