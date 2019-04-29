rm(list=ls())

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
uniformCodonDist <- rep(1, length(codons))/length(codons)
names(uniformCodonDist) <- codons
# digest()
delta5_uniform <- rep(1, 3)/3
names(delta5_uniform) <- as.character(0:2)
delta3_uniform <- rep(1, 3)/3
names(delta3_uniform) <- as.character(0:2)
minSize <- 27
maxSize <- 31
# nReads_weinberg <- 72928017
# nReads_green <- 55018963
# nReads_lareau <- 85607480
nReads <- 8e7
partSize <- 1e6
nParts <- nReads/partSize

## yeast genome params for simTranscriptome()
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

# ## human genome params for simTranscriptome()
# humanFAfile <- "gencode.v22.transcript.13cds10.fa"
# human_pad5 <- 13
# human_pad3 <- 10
# humanFAlist <- readFAfile(file.path(refDir, humanFAfile),
#                           human_pad5, human_pad3)
# # filter out transcripts w/ < 20 codons + padding
# humanFAlist <- humanFAlist[lengths(humanFAlist) > (human_pad5+human_pad3+20)] 
# humanLengths <- lengths(humanFAlist) - floor(human_pad5/3) - floor(human_pad3/3)
# humanCodonCounts <- countCodons(humanFAlist, codons, human_pad5, human_pad3)
# humanCodonDist <- rowSums(humanCodonCounts)/sum(humanCodonCounts)

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

## bias scores for ligate() & circularize() --> uniform bias
p5bias <- rep(1, length(codons))
names(p5bias) <- codons
n3bias <- rep(1, length(codons))
names(n3bias) <- codons

# simulation 1 ------------------------------------------------------------

# set.seed(107)
# 
# # 1. simulate transcriptome
# # yeast genome
# # uniform codon distribution
# yeast_uniform <- readFAfile(file.path(outputDir, "yeast_uniformCodons.fa"), pad5=3*6, pad3=3*4)
# 
# # 2. simulate ribosome distributions
# # weinberg data for expt lengths, abundances, codon TE
# load(file.path(outputDir, "yeast_uniformCodons_rho_pi.Rda"))
# 
# # 3. simulate footprints
# # uniform delta5, delta3
# # no ligation or circularization bias
# # minSize=13, maxSize=36
# 
# system('mkdir ../outputs/yeast_uniformCodons_uniformDelta_noBias_80Mreads')
# 
# exptDir <- "yeast_uniformCodons_uniformDelta_noBias_80Mreads"
# 
# minSize <- 13
# maxSize <- 36
# 
# for(i in 1:nParts) {
#   print(paste("Part", i, "of", nParts))
#   partName <- paste0("yeast_uniform_uniform_part", i)
#   part_filename <- paste0("yeast_uniformCodons_uniformDelta_noBias_80Mreads_part", i)
#   assign(partName, 
#          value=simFootprints(yeast_uniform, nRibosomes=partSize, 
#                              rhos=yeast_uniform_rho, pis=yeast_uniform_pi,
#                              delta5=delta5_uniform, delta3=delta3_uniform,
#                              ligBias=n3bias, circBias=p5bias,
#                              digest_transcript=digest_transcript,
#                              minSize=minSize, maxSize=maxSize))
#   writeFootprintsFA(get(partName),
#                     file.path(outputDir, exptDir, paste0(part_filename, ".fa")))
#   save(list=partName, 
#        file=file.path(outputDir, exptDir, paste0(part_filename, ".Rda")))
#   rm(list=partName)
# }

# simulation 2 ------------------------------------------------------------

set.seed(29)

# 1. simulate transcriptome
# yeast genome
# uniform codon distribution
yeast_yeast <- readFAfile(file.path(outputDir, "yeast_yeastCodons.fa"), pad5=3*6, pad3=3*4)

# 2. simulate ribosome distributions
# weinberg data for expt lengths, abundances, codon TE
load(file.path(outputDir, "yeast_yeastCodons_rho_pi.Rda"))

# 3. simulate footprints
# uniform delta5, delta3
# no ligation or circularization bias
# minSize=13, maxSize=36

system('mkdir ../outputs/yeast_yeastCodons_uniformDelta_noBias_80Mreads')

exptDir <- "yeast_yeastCodons_uniformDelta_noBias_80Mreads"

minSize <- 13
maxSize <- 36

for(i in 1:nParts) {
  print(paste("Part", i, "of", nParts))
  partName <- paste0("yeast_uniform_uniform_part", i)
  part_filename <- paste0("yeast_yeastCodons_uniformDelta_noBias_80Mreads_part", i)
  assign(partName, 
         value=simFootprints(yeast_yeast, nRibosomes=partSize, 
                             rhos=yeast_yeast_rho, pis=yeast_yeast_pi,
                             delta5=delta5_uniform, delta3=delta3_uniform,
                             ligBias=n3bias, circBias=p5bias,
                             digest_transcript=digest_transcript,
                             minSize=minSize, maxSize=maxSize))
  writeFootprintsFA(get(partName),
                    file.path(outputDir, exptDir, paste0(part_filename, ".fa")))
  save(list=partName, 
       file=file.path(outputDir, exptDir, paste0(part_filename, ".Rda")))
  rm(list=partName)
}