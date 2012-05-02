###################################################
### code chunk number 5: compendium
###################################################
library("ff")
library("genefilter")
library("IRanges")
library("MASS")
library("VanillaICE")
library("crlmmCompendium")


###################################################
### code chunk number 36: loadObject-exampleData1
###################################################
if(!exists("exampleData1")) data(exampleData1)


###################################################
### code chunk number 39: defineLatticeObjects
###################################################
ldat <- prePredictPanel(exampleData1)
shades <- makeTransparent(brewer.pal(6, "BrBG"), alpha=0.6)[c(1,2,3,5,6)]
##replace the middle color (white) with something else
mykey <- simpleKey(as.character(0:4), points=FALSE, rectangles=TRUE, col="black", space="right", cex=0.7)
mykey$rectangles[["col"]] <- shades
(bvnfig <- xyplot(A~B|snp, ldat, cex=0.3, panel=cnPredictionPanel, object=exampleData1,
	       x.axis="B", copynumber=0:4, line.col=shades, line.lwd=1.5,
	       shades=shades, ylab=expression(log[2](I[A])), xlab=expression(log[2](I[B])),
	       par.strip.text=list(lines=0.9, cex=0.6),
		xlim=c(6,12.5), ylim=c(6,12.5),
	       key=mykey))


###################################################
### code chunk number 40: ABscatterplots
###################################################
pars <- trellis.par.get()
pars$axis.text$cex <- 0.5
pars$xlab.text$cex <- 0.7
pars$ylab.text$cex <- 0.7
trellis.par.set("axis.text", pars$axis.text)
trellis.par.set("xlab.text", pars$xlab.text)
trellis.par.set("ylab.text", pars$ylab.text)
print(bvnfig)


