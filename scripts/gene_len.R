rm(list=ls())

source("~/simRiboSeq/scripts/helper.R")

yeast_uniform <- readFAfile("~/simRiboSeq/outputs/yeast_uniformCodons.fa", 
                            pad5=3*6, pad3=3*4)
yeast_uniform_lengths <- data.frame(gene_name = names(yeast_uniform),
                                    UTR5_length = 3*6,
                                    CDS_length = 3*(lengths(yeast_uniform)-6-4),
                                    UTR3_length=3*4)
write.table(yeast_uniform_lengths, file="~/iXnos/genome_data/yeast_uniformCodons_lengths.txt",
            sep="\t", row.names=F, col.names=F, quote=F)

yeast_yeast <- readFAfile("~/simRiboSeq/outputs/yeast_yeastCodons.fa", 
                            pad5=3*6, pad3=3*4)
yeast_yeast_lengths <- data.frame(gene_name = names(yeast_yeast),
                                    UTR5_length = 3*6,
                                    CDS_length = 3*(lengths(yeast_yeast)-6-4),
                                    UTR3_length=3*4)
write.table(yeast_yeast_lengths, file="~/iXnos/genome_data/yeast_yeastCodons_lengths.txt",
            sep="\t", row.names=F, col.names=F, quote=F)
