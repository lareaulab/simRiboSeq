% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_footprints.R
\name{digest_transcript}
\alias{digest_transcript}
\title{Generate footprint sequences from ribosome profile}
\usage{
digest_transcript(
  codon_sequence,
  rpf_counts,
  transcript_name,
  delta_5,
  delta_3
)
}
\arguments{
\item{codon_sequence}{character vector; transcript sequence in codons}

\item{rpf_counts}{numeric vector; footprint counts per codon position}

\item{transcript_name}{character; name of transcript}

\item{delta_5}{named numeric_vector; probabilities of 5' digest lengths}

\item{delta_3}{named numeric vector; probabilities of 3' digest lengths}
}
\value{
list of \code{footprint} objects
}
\description{
Generate footprint sequences from ribosome profile
}
