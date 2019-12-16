# leave-one-out analysis of codon position contributions

resultsDir <- "/mnt/lareaulab/amok/iXnos/results"

# expts <- c("scer_biasDelta_withBias_extraBase",
#            "scer_biasDelta_p5bias_extraBase")

expts <- c("scer_biasDelta_nop5bias",
           "scer_biasDelta_nop5bias_corrected")

# plots with expt as title
for(expt in expts) {
  leaveout_fname = file.path(resultsDir, expt, "leaveout_series/lr_leaveout_corrs.txt")
  full_fname = file.path(resultsDir, expt, "feat_neighborhood_series/linreg_corrs.txt")
  full.name = "lr_cod_n7p5_nt_n21p17"
  dt = read.delim(leaveout_fname, header = T, comment = "#")
  full = read.delim(full_fname, header = T, comment = "#")
  leaveout.mean = dt$test_corr
  full.row = which( full$model == full.name)
  full.mean = full$test_corr[full.row]
  label1 = c("all", NA, -7, NA, -5, NA, -3, NA,  "P", NA,  1,  NA, 3,  NA, 5)
  label2 = c(NA,     NA, NA, -6, NA, -4, NA, "E", NA,  "A", NA, 2,  NA, 4,  NA)
  par( mex = 0.65 )
  par( mar = c(6,5.5,4,3) )
  par( oma = c(0,1.5,1,0) )
  par( lwd = 0.75 )
  par( xpd = NA )
  if(grepl("yeast|scer", expt)) {
    # yMax <- ceiling(10*summary(full.mean-leaveout.mean)["Max."])/10
    yMax <- 0.6
  } else {
    yMax <- 0.2
  }
  plot( NA, 
        xlim = c( 1, 15 ),
        ylim = c( 0, yMax ),
        axes = F,
        xlab = "codon position",
        ylab = expression(paste(Delta, " correlation")),
        main = paste0(expt, ": linear regression"))
  rect( c(2.6:14.6), 0, c(3.4:15.4), full.mean - leaveout.mean,
        col = "mediumpurple3", border = NA )
  yLab <- as.character(seq(0,yMax,by=0.05))
  yLab[2*(1:(length(yLab)/2))] <- NA
  axis( 2, at = seq(0,yMax,by=0.05), labels = yLab, lwd = 0.75 )
  axis( 1, at = 1:15, padj = -1, labels = label1, tick = F, cex.axis = 0.7)
  axis( 1, at = 1:15, padj = -1, labels = label2, tick = F, cex.axis = 0.7)
}

