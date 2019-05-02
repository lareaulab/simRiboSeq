rm(list=ls())

### evaluate nucleotide frequency

# uniform genome
uniformGenome <- readLines("~/iXnos/genome_data/yeast_uniformCodons.fa")
uniformGenome_nGenes <- length(uniformGenome)/2
uniformGenome <- uniformGenome[2*(1:uniformGenome_nGenes)]
uniformGenome_nt <- strsplit(uniformGenome, split="")
uniformGenome_nt_count <- t(sapply(uniformGenome_nt,
                                   function(x) {
                                     c(sum(x=="A"),
                                       sum(x=="T"),
                                       sum(x=="C"),
                                       sum(x=="G"))
                                   }))
uniformGenome_nt_freq <- colSums(uniformGenome_nt_count)/sum(uniformGenome_nt_count)
names(uniformGenome_nt_freq) <- c("A", "T", "C", "G")

# yeast genome
yeastGenome <- readLines("~/iXnos/genome_data/yeast_yeastCodons.fa")
yeastGenome_nGenes <- length(yeastGenome)/2
yeastGenome <- yeastGenome[2*(1:yeastGenome_nGenes)]
yeastGenome_nt <- strsplit(yeastGenome, split="")
yeastGenome_nt_count <- t(sapply(yeastGenome_nt,
                                 function(x) {
                                   c(sum(x=="A"),
                                     sum(x=="T"),
                                     sum(x=="C"),
                                     sum(x=="G"))
                                 }))
yeastGenome_nt_freq <- colSums(yeastGenome_nt_count)/sum(yeastGenome_nt_count)
names(yeastGenome_nt_freq) <- c("A", "T", "C", "G")

# # plot
# par(mfrow=c(1,2))
# barplot(uniformGenome_nt_freq, 
#         main="uniform", ylab="frequency",
#         col=RColorBrewer::brewer.pal(4, "Set1"),
#         ylim=c(0, 0.35))
# barplot(yeastGenome_nt_freq,
#         main="yeast", ylab="frequency",
#         col=RColorBrewer::brewer.pal(4, "Set1"),
#         ylim=c(0, 0.35))

### evaluate codon distribution

source("~/simRiboSeq/scripts/helper.R")
codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")

# uniform genome
uniformGenome_codons <- readFAfile("~/iXnos/genome_data/yeast_uniformCodons.fa",
                                   pad5=3*6, pad3=3*34)
uniformGenome_codons_count <- t(sapply(uniformGenome_codons,
                                       function(gene) {
                                         sapply(codons,
                                                function(codon) {
                                                  sum(gene == codon)
                                                })
                                       }))
uniformGenome_codons_freq <- colSums(uniformGenome_codons_count)/sum(uniformGenome_codons_count)

# yeast genome
yeastGenome_codons <- readFAfile("~/iXnos/genome_data/yeast_yeastCodons.fa",
                                 pad5=3*6, pad3=3*34)
yeastGenome_codons_count <- t(sapply(yeastGenome_codons,
                                     function(gene) {
                                       sapply(codons,
                                              function(codon) {
                                                sum(gene == codon)
                                              })
                                     }))
yeastGenome_codons_freq <- colSums(yeastGenome_codons_count)/sum(yeastGenome_codons_count)

# # plot
# par(mfrow=c(2,1))
# barplot(uniformGenome_codons_freq,
#         main="uniform", ylab="frequency", las=2,
#         ylim=c(0,0.05))
# barplot(yeastGenome_codons_freq,
#         main="yeast", ylab="frequency", las=2,
#         ylim=c(0,0.05))

### evaluate fragment length by frame

processFile <- function(filepath) {
  ## filepath: character; path to *.wts.sam file
  print(paste("file:", filepath))
  counts <- data.frame(frame0=rep(0,5),
                       frame1=rep(0,5),
                       frame2=rep(0,5))
  con <- file(filepath, "r")
  lineCount <- 1
  while(T) {
    line <- readLines(con, n=1)
    if(length(line)==0) {
      break
    }
    if(grepl("^gene", line)) {
      line <- strsplit(line, split="\t")[[1]]
      readFrame <- as.numeric(line[4])%%3
      readLength <- nchar(line[10])
      counts[readLength-26, readFrame+1] <- counts[readLength-26, readFrame+1] + 1
    }
    lineCount <- lineCount + 1
    if((lineCount %% 100000)==0) {
      print(paste("line", lineCount))
    }
  }
  close(con)
  return(counts)
}

# uniform_noBias
uniform_noBias_samFile <- grep("wts.sam", 
                               list.files("~/iXnos/expts/yeast_uniformCodons_noBias/process/"),
                               value=T)
uniform_noBias_counts <- processFile(file.path("~/iXnos/expts/yeast_uniformCodons_noBias/process/",
                                               uniform_noBias_samFile))

# uniform
uniform_samFile <- grep("wts.sam", 
                        list.files("~/iXnos/expts/yeast_uniformCodons/process/"),
                        value=T)
uniform_counts <- processFile(file.path("~/iXnos/expts/yeast_uniformCodons/process/",
                                        uniform_samFile))

# uniform_corrected
uniform_corrected_samFile <- grep("wts.sam", 
                                  list.files("~/iXnos/expts/yeast_uniformCodons_corrected/process/"),
                                  value=T)
uniform_corrected_counts <- processFile(file.path("~/iXnos/expts/yeast_uniformCodons_corrected/process/",
                                                  uniform_samFile))

# yeast_noBias
yeast_noBias_samFile <- grep("wts.sam", 
                               list.files("~/iXnos/expts/yeast_yeastCodons_noBias/process/"),
                               value=T)
yeast_noBias_counts <- processFile(file.path("~/iXnos/expts/yeast_yeastCodons_noBias/process/",
                                               yeast_noBias_samFile))

# yeast
yeast_samFile <- grep("wts.sam", 
                        list.files("~/iXnos/expts/yeast_yeastCodons/process/"),
                        value=T)
yeast_counts <- processFile(file.path("~/iXnos/expts/yeast_yeastCodons/process/",
                                        yeast_samFile))

# yeast_corrected
yeast_corrected_samFile <- grep("wts.sam", 
                                  list.files("~/iXnos/expts/yeast_yeastCodons_corrected/process/"),
                                  value=T)
yeast_corrected_counts <- processFile(file.path("~/iXnos/expts/yeast_yeastCodons_corrected/process/",
                                                  yeast_samFile))

save.image(file=file.path("~/simRiboSeq/outputs/EDA.RData"))