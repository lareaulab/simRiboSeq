rm(list=ls())

library(ggplot2)
library(gridExtra)

### compare simulated parameters to choros learned parameters

codons <- apply(expand.grid(c("A", "T", "C", "G"),
                            c("A", "T", "C", "G"),
                            c("A", "T", "C", "G")),
                1, paste, collapse="")

delta5_uniform <- rep(1, 3)/3
names(delta5_uniform) <- as.character(0:2)
delta3_uniform <- rep(1, 3)/3
names(delta3_uniform) <- as.character(0:2)

## green expt: bias scores for ligate() & circularize()
# /mnt/lareaulab/rtunney/iXnos/results/green/s28_cod_n5p4_nt_n15p14/epoch30/codon_scores.tsv
green_biasFile <- "codon_scores.tsv"
green_biasScores <- read.table(file.path("~/simRiboSeq/refData/", green_biasFile))
colnames(green_biasScores) <- as.character(seq.int(from=-5, length.out=ncol(green_biasScores)))
green_p5bias <- green_biasScores[,"-5"]
names(green_p5bias) <- sort(codons)
green_p5bias <- (green_p5bias+1)/(max(green_p5bias, na.rm=T)+1)
green_p5bias[is.na(green_p5bias)] <- 0
green_n3bias <- green_biasScores[,"3"]
names(green_n3bias) <- sort(codons)
green_n3bias <- (green_n3bias+1)/(max(green_n3bias, na.rm=T)+1)
green_n3bias[is.na(green_n3bias)] <- 0

## uniform: /rho and /pi and /beta
load("~/simRiboSeq/outputs/yeast_uniformCodons_rho_pi.Rda")
uniform_choros_dir <- "~/choros/yeast_uniformCodons_uniformDelta_80Mreads/results/"
# rho
uniform_choros_rho <- as.numeric(readLines(file.path(uniform_choros_dir, "rho.corrected.txt")))
uniform_rho <- data.frame(simulation=yeast_uniform_rho,
                          choros=uniform_choros_rho)
uniform_rho_plot <- ggplot(uniform_rho, aes(simulation, choros)) + 
  geom_point(alpha=0.1, size=0.5) + 
  theme_classic() +
  ggtitle(expression(rho),
          subtitle=paste0("correlation = ", 
                          signif(cor(uniform_rho$simulation, uniform_rho$choros), digits=2)))
# pi
uniform_choros_pi <- readLines(file.path(uniform_choros_dir, "pi.corrected.txt"))
uniform_choros_pi <- lapply(strsplit(uniform_choros_pi, "\t"), as.numeric)
uniform_pi <- data.frame(simulation=unlist(yeast_uniform_pi),
                         choros=unlist(uniform_choros_pi))
uniform_pi$choros[is.infinite(uniform_pi$choros)] <- NA
uniform_pi_plot <- ggplot(uniform_pi, aes(simulation, choros)) + 
  geom_bin2d(bins=50) + 
  scale_fill_gradient(low="white", high="steelblue") +
  theme_classic() +
  ggtitle(expression(pi),
          subtitle=paste0("correlation = ", 
                          signif(cor(uniform_pi$simulation, uniform_pi$choros,
                                     use="complete.obs"), digits=2)))
uniform_choros_beta5 <- read.table(file.path(uniform_choros_dir, "beta5.100.txt"))
uniform_choros_beta3 <- read.table(file.path(uniform_choros_dir, "beta3.100.txt"))
uniform_beta <- data.frame(simulation=c(green_p5bias[match(codons, names(green_p5bias))],
                                        green_n3bias[match(codons, names(green_n3bias))]),
                           choros=c(uniform_choros_beta5$V2[match(codons, uniform_choros_beta5$V1)],
                                    uniform_choros_beta3$V2[match(codons, uniform_choros_beta3$V1)]))
uniform_beta$choros[is.infinite(uniform_beta$choros)] <- NA
uniform_beta_plot <- ggplot(uniform_beta, aes(simulation, choros)) +
  geom_point(alpha=0.5, size=0.7) +
  theme_classic() +
  ggtitle(expression(beta),
          subtitle=paste0("correlation = ", 
                          signif(cor(uniform_beta$simulation, uniform_beta$choros,
                                     use="complete.obs"), digits=2)))
grid.arrange(uniform_rho_plot, uniform_pi_plot, uniform_beta_plot)

## yeast: /rho and /pi and /beta
load("~/simRiboSeq/outputs/yeast_yeastCodons_rho_pi.Rda")
yeast_choros_dir <- "~/choros/yeast_yeastCodons_uniformDelta_80Mreads/results/"
# rho
yeast_choros_rho <- as.numeric(readLines(file.path(yeast_choros_dir, "rho.corrected.txt")))
yeast_rho <- data.frame(simulation=yeast_yeast_rho,
                          choros=yeast_choros_rho)
yeast_rho_plot <- ggplot(yeast_rho, aes(simulation, choros)) + 
  geom_point(alpha=0.1, size=0.5) + 
  theme_classic() +
  ggtitle(expression(rho),
          subtitle=paste0("correlation = ", 
                          signif(cor(yeast_rho$simulation, yeast_rho$choros), digits=2)))
# pi
yeast_choros_pi <- readLines(file.path(yeast_choros_dir, "pi.corrected.txt"))
yeast_choros_pi <- lapply(strsplit(yeast_choros_pi, "\t"), as.numeric)
yeast_pi <- data.frame(simulation=unlist(yeast_yeast_pi),
                         choros=unlist(yeast_choros_pi))
yeast_pi$choros[is.infinite(yeast_pi$choros)] <- NA
yeast_pi_plot <- ggplot(yeast_pi, aes(simulation, choros)) + 
  geom_bin2d(bins=50) + ylim(-17, -3) +
  scale_fill_gradient(low="white", high="steelblue") +
  theme_classic() +
  ggtitle(expression(pi),
          subtitle=paste0("correlation = ", 
                          signif(cor(yeast_pi$simulation, yeast_pi$choros,
                                     use="complete.obs"), digits=2)))
yeast_choros_beta5 <- read.table(file.path(yeast_choros_dir, "beta5.100.txt"))
yeast_choros_beta3 <- read.table(file.path(yeast_choros_dir, "beta3.100.txt"))
yeast_beta <- data.frame(simulation=c(green_p5bias[match(codons, names(green_p5bias))],
                                        green_n3bias[match(codons, names(green_n3bias))]),
                           choros=c(yeast_choros_beta5$V2[match(codons, yeast_choros_beta5$V1)],
                                    yeast_choros_beta3$V2[match(codons, yeast_choros_beta3$V1)]))
yeast_beta$choros[is.infinite(yeast_beta$choros)] <- NA
yeast_beta_plot <- ggplot(yeast_beta, aes(simulation, choros)) +
  geom_point(alpha=0.5, size=0.7) +
  theme_classic() +
  ggtitle(expression(beta),
          subtitle=paste0("correlation = ", 
                          signif(cor(yeast_beta$simulation, yeast_beta$choros,
                                     use="complete.obs"), digits=2)))
grid.arrange(yeast_rho_plot, yeast_pi_plot, yeast_beta_plot)
