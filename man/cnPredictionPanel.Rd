\name{cnPredictionPanel}
\alias{cnPredictionPanel}
\title{
	Plot bivariate normal prediction regions on a lattice panel
}
\description{
	Bivariate normal prediction regions are added to a lattice panel.
}
\usage{
cnPredictionPanel(x, y, ..., object, copynumber, x.axis, line.col,
line.lwd, shades, subscripts, data.last =FALSE, highlight.index = NULL, scale.sd = rep(1, 2))
}
\arguments{
  \item{x}{
  passed to panel.xyplot
}
  \item{y}{
  passed to panel.xyplot
}
  \item{\dots}{
  passed to panel.xyplot
}
  \item{object}{
  object of class \code{CNSet}
}
 \item{copynumber}{
 Integer vector.  Which prediction regions to plot.
}
  \item{x.axis}{
  valid entries are 'A' or 'B', indicating whether the A or B allele
  is to be plotted on the x-axis.
}
  \item{line.col}{
  the color for the border of the prediction regions 
}
  \item{line.lwd}{
  currently ignored
}
  \item{shades}{
  the fill color for the prediction regions
}
  \item{subscripts}{
   passed to panel.xyplot
}
  \item{data.last}{
  whether to plot the data before or after plotting the prediction regions
}
  \item{highlight.index}{
  an integer vector. If not missing, this vector can be used to
  highlight the position of some samples in the A vs B scatterplot.
  These samples will be plotted with an 'X' instead of an 'o'.
}
  \item{scale.sd}{
	whether to scale the A and B standard deviations to provide a
	sample specific prediction region. This is experimental and
	defaults to 1 (no scaling)
}
}
\value{
	a lattice object
}
\author{
	R. Scharpf
}

\seealso{
	\code{\link{panel.xyplot}}, \code{\link{CNSet-class}}
}
\examples{
library("crlmm")
library("crlmmCompendium")


###################################################
### chunk number 27: loadObject-exampleData1
###################################################
#line 1164 "manuscript.Rnw"
data(exampleData1)


###################################################
### chunk number 30: defineLatticeObjects
###################################################
#line 1223 "manuscript.Rnw"
ldat <- prePredictPanel(exampleData1)
shades <- makeTransparent(brewer.pal(6, "BrBG"), alpha=0.6)[c(1,2,3,5,6)]
##replace the middle color (white) with something else
mykey <- simpleKey(as.character(0:4), points=FALSE, rectangles=TRUE, col="black", space="right", cex=0.7)
mykey$rectangles[["col"]] <- shades
fig2 <- xyplot(A~B|snp, ldat, cex=0.3, panel=cnPredictionPanel, object=exampleData1,
	       x.axis="B", copynumber=0:4, line.col=shades, line.lwd=1.5,
	       shades=shades, ylab=expression(log[2](I[A])), xlab=expression(log[2](I[B])),
	       par.strip.text=list(lines=0.9, cex=0.6),
	       key=mykey)


###################################################
### chunk number 31: ABscatterplots
###################################################
#line 1237 "manuscript.Rnw"
pars <- trellis.par.get()
pars$axis.text$cex <- 0.3
pars$xlab.text$cex <- 0.8
trellis.par.set("axis.text", pars$axis.text)
trellis.par.set("axis.text", pars$xlab.text)
print(fig2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

