rm(list=ls())

source(file.path("~/simRiboSeq/scripts/helper.R"))

# read in iXnos cts_by_codon output
expt <- "yeast_uniformCodons"
ixnos_cts_filename <- file.path("~/iXnos/expts", expt, "process/cts_by_codon.size.27.31.txt")
ixnos_cts <- readCtsByCodon(ixnos_cts_filename)

# set up empty counts list
gene_len <- read.table(file.path("~/iXnos/genome_data/", paste0(expt, "_lengths.txt")),
                       col.names=c("transcript", "utr5", "cds", "utr3"),
                       stringsAsFactors=F)
nGenes <- nrow(gene_len)
sim_cts <- lapply(gene_len$cds,
                  function(x) {
                    rep(0, times=x/3)
                  })
names(sim_cts) <- gene_len$transcript

faFiles_dir <- file.path("~/simRiboSeq/outputs/", paste0(expt, "_uniformDelta_80Mreads/"))
faFiles <- list.files(faFiles_dir)
faFiles <- grep(".fa", faFiles, value=T)

for(x in faFiles) {
  print(x)
  tmp <- readLines(file.path(faFiles_dir, x))
  nReads <- length(tmp)/2
  footprints <- data.frame(t(simplify2array(strsplit(tmp[2*(1:nReads)-1], split="_"))),
                           row.names=NULL,
                           stringsAsFactors=F)
  colnames(footprints) <- c("transcript", "Asite", "id")
  footprints$transcript <- sub(">", "", footprints$transcript)
  footprints$Asite <- as.numeric(footprints$Asite)
  for(x in 1:nReads) {
    sim_cts[[which(names(sim_cts)==footprints$transcript[x])]][footprints$Asite[x]] <- sim_cts[[which(names(sim_cts)==footprints$transcript[x])]][footprints$Asite[x]] + 1
  }
  rm(tmp)
}

save(sim_cts, ixnos_cts, expt, 
     file=file.path("~/simRiboSeq/outputs/", paste0(expt, "uniformDelta_80Mreads_counts.Rda")))