% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_footprints.R
\name{simulate_footprints}
\alias{simulate_footprints}
\title{Simulate footprint sequences and pass through experimental procedures}
\usage{
simulate_footprints(
  transcript_seq,
  num_ribosomes,
  rhos,
  pis,
  delta_5,
  delta_3,
  min_size = 27,
  max_size = 31,
  lig_bias,
  rt_bias,
  circ_bias,
  mc.cores = NULL
)
}
\arguments{
\item{transcript_seq}{list of character vectors; output from read_fasta() or simulate_transcriptome()}

\item{rhos}{numeric vector; multinomial probabilities for ribosomes per transcript; output from simulate_rho()}

\item{pis}{list of numeric vectors; multinomial probabilities for ribosomes per codon position; output from simulate_pi()}

\item{delta_5}{named numeric_vector; probabilities of 5' digest lengths}

\item{delta_3}{named numeric vector; probabilities of 3' digest lengths}

\item{min_size}{integer; minimum RPF size}

\item{max_size}{integer; maximum RPF size}

\item{lig_bias}{named numeric vector; probabilities for successful 3' ligation}

\item{rt_bias}{named numeric vector; probabilities of adding additional nt}

\item{circ_bias}{named numeric vector; probabilities for successful 5' circularization}

\item{mc.cores}{integer; number of cores to use for parallelization}

\item{num_ribosomes:}{integer; number of footprints to generate}
}
\value{
list of \code{footprint} objects
}
\description{
Simulate footprint sequences and pass through experimental procedures
}
