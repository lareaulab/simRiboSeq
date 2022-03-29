#' Generate per-transcript probabilities from loess model trained on experimental data
#' 
#' @export
#' @param transcript_seq list of character vectors; output from read_fasta() or simulate_transcriptome()
#' @param expt_lengths numeric vector; transcript lengths from a model experiment
#' @param expt_abundances numeric vector; transcript abundances from a model experiment
#' @return numeric vector of per-transcript probabilities
simulate_rho <- function(transcript_seq, expt_lengths, expt_abundances) {
  transcript_lengths <- lengths(transcript_seq) - 10 # assume extra 10 codons for UTR padding
  rho_model <- loess(expt_abundances ~ expt_lengths) # model transcript abundance ~ length
  transcript_abundances <- predict(rho_model, transcript_lengths) # predict new abundances for transcript_seq
  transcript_abundances[is.na(transcript_abundances)] <- min(transcript_abundances, na.rm=T)/2 # for lengths outside model domain
  transcript_abundances <- transcript_abundances + rpois(length(transcript_seq), 10) - 7 # add noise
  rhos <- transcript_abundances / sum(transcript_abundances) # scale to probabilities, sum to 1
  names(rhos) <- names(transcript_seq)
  return(rhos)
}

#' Generate per-codon-position probabilities
#' 
#' @export
#' @param transcript_seq list of character vectors; output from read_fasta() or simulate_transcriptome()
#' @param utr5_length integer; length of 5' UTR in "codons"
#' @param utr3_length integer; length of 3' UTR in "codons
#' @param codon_TE named numeric vector; TE of codon averaged over transcriptome
simulate_pi <- function(transcript_seq, utr5_length=6, utr3_length=4, codon_TE) {
  transcript_codons <- lapply(transcript_seq,
                              function(transcript) {
                                transcript[seq.int(utr5_length+1,
                                                   length(transcript)-utr3_length)]
                              })
  pis <- lapply(transcript_codons,
                function(transcript) {
                  codon_pi <- codon_TE[match(transcript, names(codon_TE))] # grab codon TE
                  codon_pi <- codon_pi / sum(codon_pi) # scale to probabilities, sum to 1
                  return(codon_pi)
                })
  names(pis) <- names(transcript_seq)
  return(pis)
}
