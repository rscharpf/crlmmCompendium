\name{SNR}
\alias{SNR}
\docType{data}
\title{
Signal to noise ratio estimated from the CRLMM algorithm.
}
\description{
The signal to noise ratio (SNR) is an overall measure of the separation of the
genotype clusters for a particular sample.  The SNR can be useful as a
measure of quality control.  For Affymetrix 6.0, we generally suggest
dropping samples with a SNR below 5.
}
\usage{data(SNR)}
\format{
  The format is:
 num [1:1258] 6.6 7.05 7.83 7.35 6.47 ...
}
\source{
	HapMap phase III
}
\examples{
data(SNR)
\dontrun{
	##for an object of class \code{CNSet}, the SNR can be extracted by
	object$SNR[]
}
}
\keyword{datasets}
