\name{lmPanel}
\alias{lmPanel}

\title{
	For overlaying the line from the copy number model over
	boxplots of the genotype cluster intensities (intensity scale)
}
\description{
	For overlaying the line from the copy number model over
	boxplots of the genotype cluster intensities (intensity scale)
}
\usage{
lmPanel(..., line.col = "grey30", label.cex = 0.8, Ns, allele = "A", nu, ph, subscripts, ltext.y = 2500)
}
\arguments{
  \item{\dots}{
  passed to panel.bwplot
}
  \item{line.col}{
  color of line segment
}
  \item{label.cex}{
  size of font for plotting the number of samples for each genotype
}
  \item{Ns}{
  The number of samples for each genotype
}
  \item{allele}{
  which allele -- must be 'A' or 'B'
}
  \item{nu}{
  the intercept from the linear model
}
  \item{ph}{
  the slope from the linear model
}
  \item{subscripts}{
  ignored
}
  \item{ltext.y}{
  height to plot the genotype frequencies
}
}

\value{
	nothing
}

\author{
R. Scharpf
}

\seealso{
	\code{\link{panel.bwplot}}
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
### chunk number 28: dataCnFigs
###################################################
#line 1175 "manuscript.Rnw"
a <- t(as.matrix(A(exampleData1)))
gt <- t(as.matrix(calls(exampleData1)))
nuA <- as.numeric(nu(exampleData1, "A"))
phA <- as.numeric(phi(exampleData1, "A"))
col <- brewer.pal(7, "Accent")[c(1, 4, 7)]
NN <- Ns(exampleData1, i=1:16, j=1)[, , 1]
fns <- featureNames(exampleData1)
snpId <- matrix(fns, nrow(a), ncol(a), byrow=TRUE)
snpId <- factor(snpId, levels=fns, ordered=TRUE)  ##IMPORTANT
ldat <- data.frame(A=as.integer(a),
		   gt=as.factor(c("AA", "AB", "BB")[as.integer(gt)]),
		   snp=snpId)
(fig <- bwplot(A~gt|snp, ldat, cex=0.6, panel=lmPanel, nu=nuA, ph=phA,
	      fill="lightblue", Ns=NN,
	      par.strip.text=list(lines=0.9, cex=0.6), ylab=expression(I[A]), ltext.y=2500))


###################################################
### chunk number 29: boxplotA
###################################################
#line 1193 "manuscript.Rnw"
pars <- trellis.par.get()
pars$axis.text$cex <- 0.3
pars$xlab.text$cex <- 0.8
trellis.par.set("axis.text", pars$axis.text)
trellis.par.set("axis.text", pars$xlab.text)
print(fig)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

