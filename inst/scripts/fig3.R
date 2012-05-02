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
### code chunk number 37: dataCnFigs
###################################################
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
boxplotfig <- bwplot(A~gt|snp, ldat, cex=0.6, panel=lmPanel, nu=nuA,
		     ph=phA, fill="lightblue", Ns=NN,
		     par.strip.text=list(lines=0.9, cex=0.6),
		     ylab=expression(I[A]),
		     xlab=expression(I[B]), ltext.y=2500, label.cex=0.6)


###################################################
### code chunk number 38: boxplotA
###################################################
pars <- trellis.par.get()
pars$axis.text$cex <- 0.7
pars$xlab.text$cex <- 0.8
pars$ylab.text$cex <- 0.8
trellis.par.set("axis.text", pars$axis.text)
trellis.par.set("xlab.text", pars$xlab.text)
trellis.par.set("ylab.text", pars$ylab.text)
print(boxplotfig)


