\name{prePredictPanel}
\alias{prePredictPanel}
\title{
	Collect A and B intensities, etc from a CNSet object into a data.frame
}
\description{
	Creates a data.frame to be passed to the lattice function xyplot
}
\usage{
prePredictPanel(object, take.log = TRUE)
}
\arguments{
  \item{object}{
  object of class \code{CNSet}
}
  \item{take.log}{
  logical
}
}

\value{
object of class \code{data.frame}
}

\author{
R. Scharpf
}
\seealso{
	\code{\link{CNSet-class}}, \code{\link{xyplot}}
}
\examples{
library(crlmmCompendium)
library(crlmm)
data(genotypeSet)
df <- prePredictPanel(genotypeSet)
str(df)
}
\keyword{manip}

