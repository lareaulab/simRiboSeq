# leave-one-out analysis of codon position contributions

resultsDir <- "/mnt/lareaulab/amok/iXnos/results"

expts <- grep("yeast_", list.files(resultsDir), value=T)

for(expt in expts) {
  
  leaveout_fname = file.path(resultsDir, expt, "leaveout_series/leaveout_corrs.txt")
  full_fname = file.path(resultsDir, expt, "feat_neighborhood_series/feat_neighborhood_corrs.txt")
  out_fname = file.path("/mnt/lareaulab/amok/simRiboSeq/outputs",paste0(expt, "_codoncorrs", ".pdf"))
  
  full.name = "full_cod_n7p5_nt_n21p17"
  
  dt = read.delim(leaveout_fname, header = T, comment = "#")
  full = read.delim(full_fname, header = T, comment = "#")
  
  leaveout.mean = apply(dt[1:13,2:11], 1, mean)
  
  full.row = which( full$rep_series == full.name)
  full.mean = apply(full[full.row,2:11], 1, mean)
  
  data = data.frame( cbind( t(full[full.row,2:11]), NA, t(dt[1:13,2:11]) ))
  
  label1 = c("all", NA, -7, NA, -5, NA, -3, NA,  "P", NA,  1,  NA, 3,  NA, 5)
  label2 = c(NA,     NA, NA, -6, NA, -4, NA, "E", NA,  "A", NA, 2,  NA, 4,  NA)
  
  # pdf( out_fname, width=2.3, height=1.67, pointsize = 7, useDingbats = F, bg = "white" )
  par( mex = 0.65 )
  par( mar = c(6,5.5,4,3) )
  par( oma = c(0,1.5,1,0) )
  par( lwd = 0.75 )
  par( xpd = NA )
  
  # summary(full.mean-leaveout.mean)
  yMax <- ceiling(10*summary(full.mean-leaveout.mean)["Max."])/10
  
  plot( NA, 
        xlim = c( 1, 15 ),
        ylim = c( 0, yMax ),
        axes = F,
        xlab = "codon position",
        ylab = expression(paste(Delta, " correlation")))
  
  rect( c(2.6:14.6), 0, c(3.4:15.4), full.mean - leaveout.mean,
        col = "mediumpurple3", border = NA )
  
  stripchart( full.mean - data, 
              vertical = T, 
              pch = 16,
              col = rgb(0.2,0.2,0.2,.3),
              cex = .4,
              add = T
  )
  yLab <- as.character(seq(0,yMax,by=0.05))
  yLab[2*(1:(length(yLab)/2))] <- NA
  
  axis( 2, at = seq(0,yMax,by=0.05), labels = yLab, lwd = 0.75 )
  axis( 1, at = 1:15, padj = -1, labels = label1, tick = F, cex.axis = 0.7)
  axis( 1, at = 1:15, padj = -1, labels = label2, tick = F, cex.axis = 0.7)
  
  # dev.off()
  
}
