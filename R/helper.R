#' Read in transcriptome fasta file, convert into vector of codons
#' 
#' @export
#' @param fasta_file character; path to transcriptome fasta file
#' @param utr5_length integer; length of 5' UTR
#' @param utr3_length integer; length of 3' UTR
#' @return list of character vectors of mRNA transcripts in trinucleotides
read_fasta <- function(fasta_file, utr5_length, utr3_length) {
  # TODO: allow flexible 5' UTR and 3' UTR lengths
  fasta_seq <- Biostrings::readDNAStringSet(fasta_file)
  transcript_names <- sapply(names(fasta_seq), 
                             function(x) strsplit(x, split=" ")[[1]][1])
  fasta_seq <- as.character(fasta_seq)
  fasta_seq <- lapply(fasta_seq,
                      function(x) {
                        num_codons <- floor(nchar(x)/3)
                        start_offset <- utr5_length %% 3
                        codon_sequence <- substring(x,
                                                    first=(3*(seq(num_codons)-1)+1)+start_offset,
                                                    last=(3*seq(num_codons))+start_offset)
                        names(codon_sequence) <- as.character(seq.int(from=-floor(utr5_length/3),
                                                                      length.out=num_codons))
                        return(codon_sequence)
                      })
  names(fasta_seq) <- transcript_names
  return(fasta_seq)
}

#' Compute codon proportions in transcriptome
#' 
#' @export
#' @param transcript_seq list of character vectors; output from read_fasta()
#' @param codons character vector; codons in transcriptome
#' @param utr5_length integer; length of 5' UTR
#' @param utr3_length integer; length of 3' UTR
#' @return matrix of codon appearances per transcript
count_codons <- function(transcript_seq, codons, utr5_length, utr3_length) {
  utr5_length <- floor(utr5_length / 3)
  utr3_length <- floor(utr3_length / 3)
  codon_counts <- sapply(transcript_seq,
                         function(transcript) {
                           trimmed_transcript <- transcript[seq.int(utr5_length+1,
                                                                    length(transcript)-utr3_length-1)]
                           sapply(codons, 
                                  function(codon) {
                                    sum(trimmed_transcript==codon)
                                  })
                         })
  return(codon_counts)
}

#' Count codon occurrences by transcript
#' 
#' @export
#' @param input_file character; path to input file
#' @return list of numeric vectors; ribosome counts per codon per transcript
read_raw_profiles <- function(input_file) {
  raw_file <- readLines(input_file)
  transcript_names <- sapply(raw_file,
                             function(x) {
                               strsplit(x, split="\t")[[1]][1]
                             })
  codon_counts <- lapply(raw_file,
                         function(x) {
                           counts <- strsplit(x, split="\t")[[1]]
                           counts <- as.numeric(counts[-1]) # first entry is transcript name
                           return(counts)
                         })
  names(codon_counts) <- transcript_names
  return(codon_counts)
}

#' Write output from simulate_transcriptome() to .fasta file
#' 
#' @export
#' @param transcripts list of character vectors; transcript sequences from simTranscriptome()
#' @param output_file character; path to output .fasta file
write_transcriptome <- function(transcripts, output_file) {
  transcripts_seq <- sapply(transcripts, paste, collapse="")
  transcripts_seq <- Biostrings::DNAStringSet(transcript_seq)
  Biostrings::writeXStringSet(transcripts_seq, output_file)
}

#' Write footprints to .fasta file
#' 
#' @export
#' @param footprints list of \code{footprint} objects
#' @param output_file character; path to output .fasta file
write_footprints_fasta <- function(footprints, output_file) {
  footprint_sequences <- sapply(footprints, function(x) x@sequence)
  names(footprint_sequences) <- sapply(footprints,
                                       function(x) {
                                         paste(x@transcript, x@Asite, x@digest5, 
                                               x@digest3, x@sequence, x@id, sep="_")
                                       })
  footprint_sequences <- Biostrings::DNAStringSet(footprint_sequences)
  Biostrings::writeXStringSet(footprint_sequences, output_file)
}

#' Write footprints to .fastq file
#' 
#' @export
#' @param footprints list of \code{footprint} objects
#' @param output_file character; path to output .fastq file
#' @param adapterr character; 3' adapter to ligate onto footprint sequence
write_footprints_fastq <- function(footprints, output_file, 
                                   adapter="CTGTAGGCACCATCAAT") {
  footprint_sequences <- sapply(footprints,
                                function(x) {
                                  paste0(x@sequence, adapter)
                                })
  names(footprint_sequences) <- sapply(footprints,
                                       function(x) {
                                         paste(x@transcript, x@A_site, x@digest_5, 
                                               x@digest_3, x@sequence, x@id, sep="_")
                                       })
  footprint_sequences <- Biostrings::DNAStringSet(footprint_sequences)
  Biostrings::writeXStringSet(footprint_sequences, output_file, format="fastq")
}

#' Read in cts_by_codon file from iXnos
#' 
#' @param input_file character; path to input file
#' @return list of numeric vectors; ribosome occupancies by codon position
read_cts_by_codon <- function(input_file) {
  raw_data <- readLines(input_file)
  cts_by_codon <- lapply(raw_data,
                         function(x) {
                           counts <- strsplit(x, split="\t")[[1]]
                           return(as.numeric(counts[-1]))
                         })
  names(cts_by_codon) <- sapply(raw_data,
                                function(x) {
                                  strsplit(x, split="\t")[[1]][1]
                                })
  return(cts_by_codon)
}