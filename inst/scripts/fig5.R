###################################################
### code chunk number 5: compendium
###################################################
library("ff")
library("Biobase")
library("genefilter")
library("IRanges")
library("MASS")
library("VanillaICE")
library("crlmmCompendium")


###################################################
### code chunk number 47: loadObject-redonSet
###################################################
if(!exists("redonSet")) data(redonSet)


###################################################
### code chunk number 45: center
###################################################
copyNumber(redonSet) <- copyNumber(redonSet) - median(copyNumber(redonSet), na.rm=TRUE) + 200L


###################################################
### code chunk number 48: HMM
###################################################
cnStates <- c(0, 1, 2, 2, 3, 4)


###################################################
### code chunk number 49: viterbi
###################################################
fit.cn <- hmm(redonSet, TAUP=1e10, cnStates=cnStates, is.log=FALSE)
hmm.df <- as.data.frame(fit.cn)
print(hmm.df[, c(2:4,7, 10,11)])


###################################################
### code chunk number 51: dnacopy
###################################################
CNA.object <- CNA(genomdat=copyNumber(redonSet)/100,
		  chrom=chromosome(redonSet),
		  maploc=position(redonSet),
		  data.type="logratio",
		  sampleid=sampleNames(redonSet))
smu.object <- smooth.CNA(CNA.object)


###################################################
### code chunk number 52: cbs_segment
###################################################
cbs.segments <- segment(smu.object)
print(cbs.segments, showSegRows=TRUE)


###################################################
### code chunk number 53: coerce2Iranges
###################################################
cbs.out <- cbs.segments$output
cbs.segs1 <- RangedData(IRanges(cbs.out$loc.start, cbs.out$loc.end),
			numMarkers=cbs.out$num.mark,
			seg.mean=cbs.out$seg.mean,
			chrom=8L)


###################################################
### code chunk number 54: oligo.data.frame
###################################################
cbs.segs1 <- addCentromereBreaks(cbs.segs1)
cn <- as.numeric(copyNumber(redonSet))
gt <- as.integer(as.matrix(calls(redonSet)))
df <- data.frame(cn=cn, gt=gt, position=position(redonSet)/1e6)


###################################################
### code chunk number 55: redonFigSetup
###################################################
genotype.cols <- c("lightblue", "green3", "lightblue")
states <- unique(as.integer(factor(fit.cn$state, levels=c(1, 3, 4, 5))))
shades <- brewer.pal(10, "PRGn")
shades <- shades[c(2,4,1,8)]
shades[3] <- "white"
shades <- makeTransparent(shades, alpha=0.6)
mykey <- simpleKey(c("hom-del", "hem-del", "normal", "duplicated")[states[order(states)]], points=FALSE,
		   rectangles=TRUE, col="black", space="top", cex=0.7)
mykey$rectangles[["col"]] <- shades[states[order(states)]]


###################################################
### code chunk number 56: Fig5
###################################################
stdev <- mad(df$cn, na.rm=TRUE)/100
redonfig <- xyplot(cn/100~position, df, pch=".", panel=cnPanel,
		   ylim=c(-0.5,6), ylab="total copy number",
		   pch.cols=genotype.cols,
		   gt=df$gt,
		   hmm.segs=fit.cn,
		   cbs.segs=cbs.segs1,
		   scales=list(x=list(tick.number=12)),
		   lwd=1,
		   shades=shades, key=mykey, xlim=c(0,150), draw.key=TRUE,
		   xlab="physical position (Mb)",
		   add.ideogram=TRUE,
		   par.strip.text=list(lines=0.9, cex=0.6))
print(redonfig)
trellis.focus("panel", 1, 1)
ltext(median(df$position), 0, "HMM states", cex=0.9)


