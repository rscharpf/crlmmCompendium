\name{makeTransparent}
\alias{makeTransparent}
\title{
	utility function for setting transparency of colors
}
\description{
	utility function for setting transparency of colors
}
\usage{
makeTransparent(color, ...)
}
\arguments{
  \item{color}{
  character string indicating color
}
  \item{\dots}{
   passed to rgb function.  See in particular the alpha argument of
   rgb -- 0 means completely transparent and alpha =1 means not transparent.
}
}
\value{

	rbg output: "A character vector with elements of 7 or 9
     characters, "#" followed by the red, blue, green and optionally
     alpha values in hexadecimal (after rescaling to '0 ... 255')."

}

\author{
R. Scharpf
}

\seealso{
	\code{\link{rgb}}
}
\examples{
	makeTransparent("black", alpha=0.5)
}
\keyword{manip}

