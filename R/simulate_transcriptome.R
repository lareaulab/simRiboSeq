#' Generate a simulated transcriptome
#' 
#' @export
#' @param transcript_lengths numeric vector; transcript lengths in codons
#' @param codon_distribution named numeric vector; codon proportions in transcriptome
#' @return list of \code{length(transcript_lengths)} of character vectors corresponding to transcript sequences
simulate_transcriptome <- function(transcript_lengths, codon_distribution) {
  num_transcripts <- length(transcript_lengths)
  transcripts <- lapply(transcriptLengths+10, generate_transcript, codon_distribution) 
  # extra 10 codon padding for footprint + digestion
  names(transcripts) <- paste0("gene", 1:num_transcripts)
  return(transcripts)
}

#' Generate simulated mRNA transcript (5' UTR + CDS + 3' UTR)
#' 
#' @param num_codons integer; length of  transcript in codons
#' @param codon_distribution named numeric vector; codon proportions in transcriptome
#' @return character vector of \code{length(num_codons)} of mRNA transcript in trinucleotides
generate_transcript <- function(num_codons, codon_distribution) {
  # TODO: add test for valid codon_distribution for sample()
  codon_sequence <- sample(names(codon_distribution), 
                           size=num_codons, replace=T, 
                           prob=codon_distribution)
  names(codon_sequence) <- as.character(seq.int(from=-6, length.out=num_codons))
  # 6 codons of 5' UTR
  codon_sequence[7] <- "ATG" 
  # 4 codons of 3' UTR
  codon_sequence[length(codon_sequence)-4] <- sample(c("TAG", "TAA", "TGA"), size=1) 
  return(codon_sequence)
}
