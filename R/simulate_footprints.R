#' Generate simulated profile of ribosome occupancy by codon position per transcript
#' 
#' @param num_ribosomes integer; number of footprints to simulate
#' @param rhos numeric vector; multinomial probabilities for ribosomes per transcript; output from simulate_rho()
#' @param pis list of numeric vectors; multinomial probabilities for ribosomes per codon position; output from simulate_pi()
#' @return list of numeric vector of ribosome occupancies by codon position
simulate_profiles <- function(nRibosomes, rhos, pis) {
  print("... ... generating raw profile")
  rpf_per_transcript <- rmultinom(n=1, size=nRibosomes, prob=rhos)
  raw_profile <- mapply(rmultinom, n=1, size=rpf_per_transcript, prob=pis)
  raw_profile <- lapply(raw_profile, function(x) x[,1])
  names(raw_profile) <- names(rhos)
  return(raw_profile)
}

#' Class for footprint objects
#' 
#' @export
footprint <- setClass("footprint", slots=list(sequence="character",
                                              transcript="character",
                                              A_site="numeric",
                                              digest_5="numeric",
                                              digest_3="numeric",
                                              id="numeric"))

#' Generate footprint sequences from ribosome profile
#' 
#' @param codon_sequence character vector; transcript sequence in codons
#' @param rpf_counts numeric vector; footprint counts per codon position
#' @param transcript_name character; name of transcript
#' @param delta_5 named numeric_vector; probabilities of 5' digest lengths
#' @param delta_3 named numeric vector; probabilities of 3' digest lengths
#' @return list of \code{footprint} objects
digest_transcript <- function(codon_sequence, rpf_counts, transcript_name, 
                              delta_5, delta_3) {
  start_codon <- which(names(codon_sequence) == "0")
  start_position <- 3*(start_codon-1) + 1
  # 1. convert codon_sequence to single character 
  nt_sequence <- paste(codon_sequence, collapse="")
  # 2. convert rpf_counts to numeric vector of A site coordinates (in nt)
  A_sites <- unlist(mapply(rep,
                           seq.int(length(rpf_counts)),
                           rpf_counts)) - 1
  if(length(A_sites) == 0) { return(NULL) }
  A_site_positions <- 3*A_sites + start_position
  # 3. generate digest lengths
  digest_5 <- sample(as.numeric(names(delta_5)), size=length(A_site_positions),
                     prob=delta_5, replace=T)
  digest_3 <- sample(as.numeric(names(delta_3)), size=length(A_site_positions),
                     prob=delta_3, replace=T)
  # 4. generate RPF start and stop positions
  rpf_start <- A_site_positions - digest_5 # 5' digest length
  rpf_end <- A_site_positions + 2 + digest_3 # 2nt for A site codon + 3' digest length
  # 5. extract RPF sequences
  rpf_seq <- substring(nt_sequence, rpf_start, rpf_end)
  # 6. convert to 'footprint' object
  ids <- unlist(mapply(seq.int, to=length(rpf_seq)))
  rpfs <- mapply(footprint, sequence=rpf_seq, A_site=A_sites, 
                 transcript=transcript_name, digest_5=digest_5, 
                 digest_3=digest_3, id=ids)
  return(rpfs)
}

#' Generate RPF sequences for all transcripts 
#' 
#' @param transcript_seq list of character vectors; output from read_fasta() or simulate_transcriptome()
#' @param rpf_counts list of numeric vectors; output from simulate_profiles()
#' @param delta_5 named numeric_vector; probabilities of 5' digest lengths
#' @param delta_3 named numeric vector; probabilities of 3' digest lengths
#' @param min_size integer; minimum RPF size
#' @param max_size integer; maximum RPF size
#' @param mc.cores integer; number of cores to use for parallelization
#' @param digest_transcript function; function that digests individual transcripts
#' @return list of \code{footprint} objects
digest <- function(transcript_seq, rpf_counts, delta_5, delta_3, 
                   min_size=27, max_size=31, digest_transcript) {
  print("... ... digesting footprints")
  num_transcripts <- length(transcript_seq)
  # digest footprints
  footprints <- unlist(parallel::mcmapply(digest_transcript,
                                          codon_sequence=transcript_seq,
                                          rpf_counts=rpf_counts,
                                          transcript_name=names(transcript_seq),
                                          MoreArgs=list(delta_5=delta_5, delta_3=delta_3)))
  rpf_lengths <- sapply(footprints, function(x) nchar(x@sequence))
  footprints <- footprints[(rpf_lengths >= min_size) & (rpf_lengths <= max_size)]
  return(footprints)
}

#' Get sequence for nucleotides at RPF edges
#' 
#' @param footprints list of \code{footprint} objects
#' @param region character; "f5" for 5' region or "f3" for 3' region
#' @param bias_length integer; number of nucleotides to constitute end region
#' @return character vector of sequences corresponding to bias region
get_bias_region <- function(footprints, region, bias_length) {
  rpf_lengths <- sapply(footprints, function(x) nchar(x@sequence))
  if(region=="f5") {
    bias_start <- 1
    bias_end <- bias_length
  }
  if(region=="f3") {
    bias_end <- rpf_lengths
    bias_start <- bias_end - bias_length + 1
  }
  bias_regions <- unlist(parallel::mcmapply(substring,
                                            text=sapply(footprints, function(x) x@sequence),
                                            first=bias_start,
                                            last=bias_end))
  return(bias_regions)
}

#' Apply 3' recovery bias to footprints
#' 
#' @param footprints list of \code{footprint} objects
#' @param lig_bias named numeric vector; probabilities for successful 3' ligation
#' @return subset of input \code{footprints} that were successfully "ligated"
ligate <- function(footprints, lig_bias) {
  print("... ... ligating footprints ...")
  bias_length <- nchar(names(lig_bias)[1])
  bias_regions <- get_bias_region(footprints, region="f3", bias_length=bias_length)
  ligation_prob <- lig_bias[match(bias_regions, names(lig_bias))]
  ligation_keep <- sapply(ligation_prob, function(p) rbinom(1, 1, p))
  footprints <- footprints[ligation_keep == 1]
  return(footprints)
}

#' Apply 5' recovery bias to footprints
#' 
#' @param footprints list of \code{footprint} objects
#' @param circ_bias named numeric vector; probabilities for successful 5' circularization
#' @return subset of input \code{footprints} that were successfully "circularized"
circularize <- function(footprints, circ_bias) {
  print("... ... circularizing footprints ...")
  bias_length <- nchar(names(circ_bias)[1])
  bias_regions <- get_bias_region(footprints, region="f5", bias_length=bias_length)
  circularization_prob <- circ_bias[match(bias_regions, names(circ_bias))]
  circularization_keep <- sapply(circularization_prob, function(p) rbinom(1, 1, p))
  footprints <- footprints[circularization_keep == 1]
  return(footprints)
}

#' Add extra non-templated base at 5' end of footprints
#' 
#' @param footprints list of \code{footprint} objects
#' @param rt_bias named numeric vector; probabilities of adding additional nt
#' @return input \code{footprints} with some footprints added NTA base
RT <- function(footprints, rt_bias) {
  print("... ... reverse transcribing ...")
  bases <- names(rt_bias)
  bases_to_add <- sample(bases, size=length(footprints), replace=T, prob=rt_bias)
  for(i in which(bases_to_add != "")) {
    footprints[[i]]@sequence <- paste0(bases_to_add[i], footprints[[i]]@sequence, collapse="")
  }
  return(footprints)
}

#' Simulate footprint sequences and pass through experimental procedures
#' 
#' @export
#' @param transcript_seq list of character vectors; output from read_fasta() or simulate_transcriptome()
#' @param num_ribosomes: integer; number of footprints to generate
#' @param rhos numeric vector; multinomial probabilities for ribosomes per transcript; output from simulate_rho()
#' @param pis list of numeric vectors; multinomial probabilities for ribosomes per codon position; output from simulate_pi()
#' @param delta_5 named numeric_vector; probabilities of 5' digest lengths
#' @param delta_3 named numeric vector; probabilities of 3' digest lengths
#' @param min_size integer; minimum RPF size
#' @param max_size integer; maximum RPF size
#' @param lig_bias named numeric vector; probabilities for successful 3' ligation
#' @param rt_bias named numeric vector; probabilities of adding additional nt
#' @param circ_bias named numeric vector; probabilities for successful 5' circularization
#' @param mc.cores integer; number of cores to use for parallelization
#' @return list of \code{footprint} objects
simulate_footprints <- function(transcript_seq, num_ribosomes, rhos, pis, 
                                delta_5, delta_3, min_size=27, max_size=31, 
                                lig_bias, rt_bias, circ_bias, mc.cores=NULL) {
  if(is.null(mc.cores)) {mc.cores <- parallel::detectCores()-10}
  print(paste("Number of cores:", mc.cores))
  rpf_cluster <- parallel::makeCluster(mc.cores)
  on.exit(parallel::stopCluster(rpf_cluster))
  num_rounds <- 1
  footprints <- c()
  while(length(footprints) < num_ribosomes) {
    print(paste("... Round", num_rounds, "of simulating footprints"))
    tmp_cts <- simulate_profiles(num_ribosomes, rhos, pis)
    tmp_rpf <- digest(transcript_seq, tmp_cts, delta_5, delta_3, min_size, max_size,
                      digest_transcript=digest_transcript)
    tmp_rpf <- ligate(tmp_rpf, lig_bias)
    tmp_rpf <- circularize(tmp_rpf, circ_bias)
    if(num_rounds == 1) {
      footprints <- tmp_rpf
    } else {
      if(length(tmp_rpf) > (num_ribosomes-length(footprints))) {
        tmp_rpf <- sample(tmp_rpf, size=num_ribosomes-length(footprints))
      }
      footprints <- append(footprints, tmp_rpf)
    }
    print(paste("...", length(footprints), "of", num_ribosomes, 
                "(", round(length(footprints)/num_ribosomes*100, digits=2), 
                " % ) footprints simulated"))
    num_rounds <- num_rounds + 1
    rm(tmp_cts, tmp_rpf)
  }
  return(footprints)
}
