\name{getHapMapIds}
\alias{getHapMapIds}
\title{
	Maps cel file names to HapMap ids
}
\description{
	Maps cel file names to HapMap ids 
}
\usage{
getHapMapIds(object)
}
\arguments{
  \item{object}{
  object of class \code{CNSet} in which the sample names are the CEL
  file names for HapMap phase 3
}
}
\value{
	character vector of HapMap identifiers (e.g., "NAXXXX")
}
\author{
R. Scharpf
}
\seealso{
	\code{\link{CNSet-class}}
}
\keyword{manip}
