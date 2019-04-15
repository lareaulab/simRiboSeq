library(foreach)
library(doParallel)

setwd("~/simRiboSeq/outputs/")

faFiles <- grep(".*part.*Rda", list.files(), value=T)
faFiles <- faFiles[-c(1:5)]

nCores <- parallel::detectCores()-10
cl <- makeCluster(nCores)
registerDoParallel(cl)

foreach(x=1:length(faFiles)) %dopar% {
  print(paste("File", x, "of", length(faFiles)))
  load(file=faFiles[x])
  objName <- sub("Delta_80Mreads", "", sub("Codons", "", sub(".Rda", "", faFiles[x])))
  nFootprints <- length(get(objName))
  outputFA <- rep(NULL, 2*nFootprints)
  footprintNames <- sapply(get(objName), function(x) paste(x@transcript, x@ASite, sep="_"))
  footprintNames <- paste(footprintNames, seq.int(nFootprints), sep="_")
  footprintSequences <- sapply(get(objName), function(x) x@sequence)
  outputFA[2*(1:nFootprints)-1] <- paste0(">", footprintNames)
  outputFA[2*(1:nFootprints)] <- footprintSequences
  outFile <- sub(".Rda", ".fa", faFiles[x])
  writeLines(outputFA, con=outFile, sep="\n")
  rm(list=objName)
}

stopCluster(cl)

yeastCodonFiles <- grep("yeastCodon.*part.*Rda", list.files(), value=T)

for(x in 1:length(yeastCodonFiles)) {
  print(paste("File", x, "of", length(yeastCodonFiles)))
  load(file=yeastCodonFiles[x])
  objName <- sub("Delta_80Mreads", "", sub("Codons", "", sub(".Rda", "", yeastCodonFiles[x])))
  fileName <- sub(".Rda", "", sub("yeastDelta", "uniformDelta", yeastCodonFiles[x]))
  new_objName <- sub("Delta_80Mreads", "", sub("Codons", "", sub(".Rda", "", fileName)))
  assign(new_objName, get(objName))
  nFootprints <- length(get(objName))
  outputFA <- rep(NULL, 2*nFootprints)
  footprintNames <- sapply(get(objName), function(x) paste(x@transcript, x@ASite, sep="_"))
  footprintNames <- paste(footprintNames, seq.int(nFootprints), sep="_")
  footprintSequences <- sapply(get(objName), function(x) x@sequence)
  outputFA[2*(1:nFootprints)-1] <- paste0(">", footprintNames)
  outputFA[2*(1:nFootprints)] <- footprintSequences
  writeLines(outputFA, con=paste0(fileName, ".fa"), sep="\n")
  save(list=new_objName, file=paste0(fileName, ".Rda"))
  rm(list=c(new_objName, objName))
}
