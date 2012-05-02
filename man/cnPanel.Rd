\name{cnPanel}
\alias{cnPanel}
\title{
Lattice panel function for plotting raw copy number estimates versus
physical position 
}
\description{
Lattice panel function for plotting raw copy number estimates versus
physical position.  The segmentation results are overplotted as line
segments.  Hidden states from a hidden markov model are plotted just
beneath the 0 y-coordinate.
}
\usage{
cnPanel(x, y, ..., pch.cols, gt, cbs.segs, hmm.segs, shades, subscripts, add.ideogram = TRUE)
}
\arguments{
  \item{x}{
  physical position (integer vector) passed to panel.xyplot
}
  \item{y}{
  raw copy number estimates (numerical vector) passed to panel.xyplot 
}
  \item{\dots}{
  additional arguments to panel.xyplot
}
  \item{pch.cols}{
  color of plotting symbols to distinguish heterozygous genotypes (AB)
  from homozygous genotypes (AA or BB)
}
  \item{gt}{
  integer vector indicating the diallelic genotype calls (1=AA, 2=AB, 3=BB).
  Use 'NA' for nonpolymorphic loci.
}
  \item{cbs.segs}{
  RangedData object containing start and stop positions (unit is
  basepairs) and the segment mean labeled as 'seg.mean'
}
  \item{hmm.segs}{
  Results from HMM implemented in the VanillaICE package
}
  \item{shades}{
  Shades for coloring the hidden states
}
  \item{subscripts}{
  additional arguments for panel.xyplot
}
  \item{add.ideogram}{
  whether to add the ideogram at the top (y-coordinates are 5.6 - 5.9).
}
}
\value{
 Nothing
}
\author{
R. Scharpf
}

\seealso{
	\code{\link{panel.xyplot}}
}
\examples{
}
\keyword{plot}

