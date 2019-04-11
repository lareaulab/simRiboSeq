setwd("~/simRiboSeq/outputs/")

faFiles <- grep(".*part.*Rda", list.files(), value=T)

for(x in 1:length(faFiles)) {
  print(paste("File", x, "of", length(faFiles)))
  load(file=faFiles[x])
  objName <- sub("Delta_80Mreads", "", sub("Codons", "", sub(".Rda", "", faFiles[x])))
  nFootprints <- length(get(objName))
  outputFA <- rep(NULL, 2*nFootprints)
  footprintNames <- sapply(get(objName), function(x) paste(x@transcript, x@ASite, sep="_"))
  footprintNames <- paste(get(objName), seq.int(nFootprints), sep="_")
  footprintSequences <- sapply(get(objName), function(x) x@sequence)
  outputFA[2*(1:nFootprints)-1] <- paste0(">", footprintNames)
  outputFA[2*(1:nFootprints)] <- footprintSequences
  writeLines(outputFA, con=sub(".Rda", ".fa", x))
  rm(list=objName)
}