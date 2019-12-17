rm(list=ls())

# load scripts ------------------------------------------------------------

library(parallel)

scriptDir <- "../scripts"
outputDir <- "../outputs"
refDir <- "../refData"

scripts <- c("helper.R", "simTranscriptome.R", "simRibosomeDist.R", "simFootprints.R")
for(x in scripts) {source(file.path(scriptDir, x))}

nCores <- 30

# reference data ----------------------------------------------------------

## general params:
codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")

# digest()
delta5_bias <- c(0.4, 0.3, 0.2)
names(delta5_bias) <- as.character(15:17)
delta3_bias <- c(0.35, 0.4, 0.25)
names(delta3_bias) <- as.character(9:11)
minSize <- 27
maxSize <- 31

## weinberg expt: codon TE scores for simRho()
weinberg_codonTEfile <- "tunney_supp_table_2_codon_scores.csv"
weinberg_codonScores <- read.csv(file.path(refDir, weinberg_codonTEfile),
                                 header=T, stringsAsFactors=F)
weinberg_codonTE <- weinberg_codonScores$X0
weinberg_codonTE <- weinberg_codonTE - min(weinberg_codonTE) # scale up so min = 0.1
weinberg_stopCodons <- codons[which(!(codons %in% weinberg_codonScores$codon))]
weinberg_codonTE <- c(weinberg_codonTE, rep(0, length(weinberg_stopCodons))) # add TE for stop codons
names(weinberg_codonTE) <- c(weinberg_codonScores$codon, weinberg_stopCodons)

## green expt: bias scores for ligate() & circularize()
green_biasFile <- "codon_scores.tsv"
green_biasScores <- read.table(file.path(refDir, green_biasFile))
colnames(green_biasScores) <- as.character(seq.int(from=-5, length.out=ncol(green_biasScores)))
green_p5bias <- exp(green_biasScores[,"-5"])
names(green_p5bias) <- sort(codons)
green_p5bias <- green_p5bias/max(green_p5bias, na.rm=T)
green_p5bias[is.na(green_p5bias)] <- 0
green_n3bias <- exp(green_biasScores[,"3"])
names(green_n3bias) <- sort(codons)
green_n3bias <- green_n3bias/max(green_n3bias, na.rm=T)
green_n3bias[is.na(green_n3bias)] <- 0
# par(mfrow=c(2,1))
# plot(density(green_biasScores[,"-5"], na.rm=T), main="p5 bias scores from green iXnos")
# plot(density(green_p5bias), main="scaled p5 probabilities")
# plot(density(green_biasScores[,"3"], na.rm=T), main="n3 bias scores from green iXnos")
# plot(density(green_n3bias), main="scaled n3 probabilities")

## RTBias
rtBias <- c(0.8, 0.08, 0.04, 0.03, 0.05)
names(rtBias) <- c("", "A", "T", "C", "G")
noRTbias <- c(1, 0, 0, 0, 0)
names(noRTbias) <- c("", "A", "T", "C", "G")

# simulation: simWeinberg -------------------------------------------------

exptName <- "riboViz_simWeinberg"

set.seed(42)

# 0. choose top 3 translated transcripts
weinberg_cts_by_codon <- readRawProfiles(file.path(refDir, "cts_by_codon.size.27.31.txt"))
weinberg_cts_by_transcript <- sapply(weinberg_cts_by_codon, sum)
transcripts <- sort(weinberg_cts_by_transcript, decreasing=T)[1:3]

# 1. specify transcriptome
# s. cerevisiae transcriptome
scer_transcriptome <- readFAfile(file.path(refDir, "scer.transcripts.20cds20.fa"),
                                 pad5=20, pad3=20)
transcripts_seq <- scer_transcriptome[match(names(transcripts), names(scer_transcriptome))]
writeTranscriptomeFA(transcripts_seq, 
                     file.path(outputDir, paste0(exptName, "_transcripts_18cds18.fa")))

# # check codon frequencies
# transcripts_codon_freq <- sapply(codons,
#                                  function(x) {
#                                    sapply(transcripts_seq,
#                                           function(y) {
#                                             sum(x == y)
#                                           })
#                                  })

# 2. specify ribosome distributions (probabilities)
transcripts_rho <- transcripts / sum(transcripts)
transcripts_pi <- weinberg_cts_by_codon[match(names(transcripts), names(weinberg_cts_by_codon))]
transcripts_pi <- lapply(transcripts_pi,
                         function(x) {
                           x / sum(x)
                         })
save(transcripts_rho, transcripts_pi,
     file=file.path(outputDir, paste0(exptName, "_rho_pi.Rda")))

# 3. simulate footprints
# biased delta
# biased recovery: green_p5bias, green_n3bias

partSize <- 1e6
nParts <- floor(sum(transcripts)/partSize)

if(!(exptName %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, paste0(exptName, "_footprints"))))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(exptName, "_footprints_part_", i)
  assign(partName, 
         value=simFootprints(transcripts_seq, nRibosomes=partSize, 
                             rhos=transcripts_rho, pis=transcripts_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             minSize=20, maxSize=40, mc.cores=nCores, 
                             ligBias=green_n3bias, circBias=green_p5bias,
                             RTBias=noRTbias, digest_transcript=digest_transcript))
  writeFootprintsFA(get(partName),
                    file.path(outputDir, exptName, paste0(partName, ".fa")))
  save(list=partName, 
       file=file.path(outputDir, exptName, paste0(partName, ".Rda")))
  rm(list=partName)
}

system(paste0("cat ", file.path(outputDir, exptName), "/*.fa > ", file.path(outputDir), "/", exptName, ".fa"))

# 4. reconstruct cts_by_codon
footprints <- readLines(file.path(outputDir, paste0(exptName, ".fa")))
footprints <- footprints[2*(1:(length(footprints)/2))-1]
footprints <- data.frame(transcript=sapply(footprints,
                                           function(x) {
                                             sub(">", "", strsplit(x, split="_")[[1]][1])
                                           }),
                         codon=sapply(footprints,
                                      function(x) {
                                        as.numeric(strsplit(x, split="_")[[1]][2])
                                      }))
footprint_cts <- data.frame(transcript=unlist(mapply(rep, names(transcripts_pi), lengths(transcripts_pi))),
                            codon=unlist(lapply(lengths(transcripts_pi), seq.int)))
rownames(footprint_cts) <- NULL
footprint_cts$rpfCount <- sapply(1:nrow(footprint_cts),
                                  function(x) {
                                    sum(footprints$transcript==footprint_cts$transcript[x] & 
                                          footprints$codon==footprint_cts$codon[x])
                                  })
write.table(footprint_cts, file=file.path(outputDir, paste0(exptName, "_rpfCounts.txt")),
            quote=F, row.names=F, col.names=T, sep="\t")

# simulation: simASite ----------------------------------------------------

exptName <- "riboViz_simASite"

set.seed(24)

# 0. choose top 3 translated transcripts
# (same as above)

# 1. specify transcriptome
# s. cerevisiae transcriptome
# (same as above)

# 2. specify ribosome distributions (probabilities)
transcripts_rho <- transcripts / sum(transcripts)
transcripts_pi <- simPi(transcripts_seq, pad5=6, pad3=7, codonTE=weinberg_codonTE)
save(transcripts_rho, transcripts_pi,
     file=file.path(outputDir, paste0(exptName, "_rho_pi.Rda")))

# 3. simulate footprints
# biased delta
# biased recovery: green_p5bias, green_n3bias

partSize <- 1e6
nParts <- floor(sum(transcripts)/partSize)

if(!(exptName %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, paste0(exptName, "_footprints"))))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(exptName, "_footprints_part_", i)
  assign(partName, 
         value=simFootprints(transcripts_seq, nRibosomes=partSize, 
                             rhos=transcripts_rho, pis=transcripts_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             minSize=20, maxSize=40, mc.cores=nCores, 
                             ligBias=green_n3bias, circBias=green_p5bias,
                             RTBias=noRTbias, digest_transcript=digest_transcript))
  writeFootprintsFA(get(partName),
                    file.path(outputDir, exptName, paste0(partName, ".fa")))
  save(list=partName, 
       file=file.path(outputDir, exptName, paste0(partName, ".Rda")))
  rm(list=partName)
}

system(paste0("cat ", file.path(outputDir, exptName), "/*.fa > ", file.path(outputDir), "/", exptName, ".fa"))

# 4. reconstruct cts_by_codon
footprints <- readLines(file.path(outputDir, paste0(exptName, ".fa")))
footprints <- footprints[2*(1:(length(footprints)/2))-1]
footprints <- data.frame(transcript=sapply(footprints,
                                           function(x) {
                                             sub(">", "", strsplit(x, split="_")[[1]][1])
                                           }),
                         codon=sapply(footprints,
                                      function(x) {
                                        as.numeric(strsplit(x, split="_")[[1]][2])
                                      }))
footprint_cts <- data.frame(transcript=unlist(mapply(rep, names(transcripts_pi), lengths(transcripts_pi))),
                            codon=unlist(lapply(lengths(transcripts_pi), seq.int)))
rownames(footprint_cts) <- NULL
footprint_cts$rpfCount <- sapply(1:nrow(footprint_cts),
                                 function(x) {
                                   sum(footprints$transcript==footprint_cts$transcript[x] & 
                                         footprints$codon==footprint_cts$codon[x])
                                 })
write.table(footprint_cts, file=file.path(outputDir, paste0(exptName, "_rpfCounts.txt")),
            quote=F, row.names=F, col.names=T, sep="\t")
