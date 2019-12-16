rm(list=ls())

setwd("~/simRiboSeq/scripts/")

# load scripts ------------------------------------------------------------

library(parallel)

scriptDir <- "."
outputDir <- "../outputs"
refDir <- "../refData"

scripts <- c("helper.R", "simTranscriptome.R", "simRibosomeDist.R", "simFootprints.R")
for(x in scripts) {source(file.path(scriptDir, x))}

# reference data ----------------------------------------------------------

## general params:
codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")
# simTranscriptome()
# uniformCodonDist <- rep(1, length(codons))/length(codons)
# names(uniformCodonDist) <- codons
# digest()
# delta5_uniform <- rep(1, 3)/3
# names(delta5_uniform) <- as.character(15:17)
# delta3_uniform <- rep(1, 3)/3
# names(delta3_uniform) <- as.character(9:11)
delta5_bias <- c(0.4, 0.3, 0.2)
names(delta5_bias) <- as.character(15:17)
delta3_bias <- c(0.35, 0.4, 0.25)
names(delta3_bias) <- as.character(9:11)
minSize <- 27
maxSize <- 31
# nReads_weinberg <- 72928017
# nReads_green <- 55018963
# nReads_lareau <- 85607480
nReads <- 8e7
partSize <- 1e6
nParts <- nReads/partSize

# ## yeast genome params for simTranscriptome()
yeastFAfile <- "scer.transcripts.13cds10.fa"
yeast_pad5 <- 13
yeast_pad3 <- 10
yeastFAlist <- readFAfile(file.path(refDir, yeastFAfile),
                          yeast_pad5, yeast_pad3)
# filter out transcripts w/ < 20 codons + padding
yeastFAlist <- yeastFAlist[lengths(yeastFAlist) > (yeast_pad5+yeast_pad3+20)]
yeastLengths <- lengths(yeastFAlist) - floor(yeast_pad5/3) - floor(yeast_pad3/3)
yeastCodonCounts <- countCodons(yeastFAlist, codons, yeast_pad5, yeast_pad3)
yeastCodonDist <- rowSums(yeastCodonCounts)/sum(yeastCodonCounts)
rm(yeastFAlist)

## weinberg expt: gene lengths and abundances for simPi()
weinberg_file <- "cts_by_codon.size.27.31.txt"
weinberg_data <- readRawProfiles(file.path(refDir, weinberg_file))
weinberg_lengths <- lengths(weinberg_data)
weinberg_abundances <- sapply(weinberg_data, sum)
weinberg_nRibosomes <- sum(weinberg_abundances)
rm(weinberg_data)

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
# /mnt/lareaulab/rtunney/iXnos/results/green/s28_cod_n5p4_nt_n15p14/epoch30/codon_scores.tsv
green_biasFile <- "codon_scores.tsv"
green_biasScores <- read.table(file.path(refDir, green_biasFile))
colnames(green_biasScores) <- as.character(seq.int(from=-5, length.out=ncol(green_biasScores)))
green_p5bias <- exp(green_biasScores[,"-5"])
names(green_p5bias) <- sort(codons)
green_p5bias <- green_p5bias/max(green_p5bias, na.rm=T)
# green_p5bias <- green_biasScores[,"-5"]
# green_p5bias <- (green_p5bias+1)/(max(green_p5bias, na.rm=T)+1)
green_p5bias[is.na(green_p5bias)] <- 0
green_n3bias <- exp(green_biasScores[,"3"])
names(green_n3bias) <- sort(codons)
green_n3bias <- green_n3bias/max(green_n3bias, na.rm=T)
# green_n3bias <- green_biasScores[,"3"]
# green_n3bias <- (green_n3bias+1)/(max(green_n3bias, na.rm=T)+1)
green_n3bias[is.na(green_n3bias)] <- 0
# par(mfrow=c(2,1))
# plot(density(green_biasScores[,"-5"], na.rm=T), main="p5 bias scores from green iXnos")
# plot(density(green_p5bias), main="scaled p5 probabilities")
# plot(density(green_biasScores[,"3"], na.rm=T), main="n3 bias scores from green iXnos")
# plot(density(green_n3bias), main="scaled n3 probabilities")

## bias scores for ligate() & circularize() --> uniform bias
p5bias <- rep(1, length(codons))
names(p5bias) <- codons
n3bias <- rep(1, length(codons))
names(n3bias) <- codons

## RTBias
rtBias <- c(0.8, 0.08, 0.04, 0.03, 0.05)
names(rtBias) <- c("", "A", "T", "C", "G")
noRTbias <- c(1, 0, 0, 0, 0)
names(noRTbias) <- c("", "A", "T", "C", "G")

# simulation 1 ------------------------------------------------------------

### scer transcriptome
### scer ribosome distributions
### delta5_bias, delta3_bias
### green p5 bias for ligBias (3')
### no p5 bias for circBias (5')
### no rtBias
### no footprint size selection

set.seed(68)

# 1. simulate transcriptome
scerFile <- "scer.transcripts.20cds20.fa"
scer <- readFAfile(file.path(refDir, scerFile), pad5=20, pad3=20)

# 2. simulate ribosome distributions
load(file.path(outputDir, "scer_rho_pi.Rda"))

# 3. simulate footprints
# biased delta5, delta3
# minSize=27, maxSize=31
# green data for ligBias (3') and circBias (5')
# no extra base: noRTbias

simulation_name <- "scer_biasDelta_allSizes_nop5bias"
if(!(simulation_name %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, simulation_name)))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(simulation_name, "_part", i)
  part_filename <- paste0(simulation_name, "_80Mreads_part", i)
  assign(partName,
         value=simFootprints(scer, nRibosomes=partSize,
                             rhos=scer_rho, pis=scer_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             ligBias=green_p5bias, RTBias=noRTbias, circBias=p5bias,
                             digest_transcript=digest_transcript,
                             minSize=10, maxSize=50))
  writeFootprintsFA(get(partName),
                    file.path(outputDir, simulation_name, paste0(part_filename, ".fa")))
  save(list=partName,
       file=file.path(outputDir, simulation_name, paste0(part_filename, ".Rda")))
  rm(list=partName)
}

# simulation 2 ------------------------------------------------------------

### scer transcriptome
### scer ribosome distributions
### delta5_bias, delta3_bias
### green p5 bias for ligBias (3')
### green p5 bias for circBias (5')
### no rtBias
### no footprint size selection

set.seed(13)

# 1. simulate transcriptome
scerFile <- "scer.transcripts.20cds20.fa"
scer <- readFAfile(file.path(refDir, scerFile), pad5=20, pad3=20)

# 2. simulate ribosome distributions
load(file.path(outputDir, "scer_rho_pi.Rda"))

# 3. simulate footprints
# biased delta5, delta3
# minSize=27, maxSize=31
# green data for ligBias (3') and circBias (5')
# no extra base: noRTbias

simulation_name <- "scer_biasDelta_allSizes_sameBias"
if(!(simulation_name %in% list.files(outputDir))) {
  system(paste("mkdir", file.path(outputDir, simulation_name)))
}

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0(simulation_name, "_part", i)
  part_filename <- paste0(simulation_name, "_80Mreads_part", i)
  assign(partName,
         value=simFootprints(scer, nRibosomes=partSize,
                             rhos=scer_rho, pis=scer_pi,
                             delta5=delta5_bias, delta3=delta3_bias,
                             ligBias=green_p5bias, RTBias=noRTbias, circBias=green_p5bias,
                             digest_transcript=digest_transcript))
  writeFootprintsFA(get(partName),
                    file.path(outputDir, simulation_name, paste0(part_filename, ".fa")))
  save(list=partName,
       file=file.path(outputDir, simulation_name, paste0(part_filename, ".Rda")))
  rm(list=partName)
}

# exit --------------------------------------------------------------------

quit(save="no")