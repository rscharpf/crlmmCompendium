redonsetUpdates <- function(){
	data(redonSet)
	featureData(redonSet) <- oligoClasses:::ADF2GDF(featureData(redonSet))
	copyNumber(redonSet) <- integerMatrix(copyNumber(redonSet), 100)
	save(redonSet, file="../../data/redonSet.rda")
}

test_redon <- function(){
	library(crlmmCompendium)
	library(Biobase)
	library(oligoClasses)
	data(redonSet)
	checkTrue(validObject(redonSet))
}

test_hmm1 <- function(){
	library(oligoClasses)
	library(crlmm)
	library(Biobase)
	library(VanillaICE)
	library(IRanges)
	load("../../data/redonSet.rda")
	cn <- copyNumber(redonSet)
	copyNumber(redonSet) <- integerMatrix(cn, 100)
	fit <- hmm(redonSet, is.log=FALSE, p.hom=0)
	fit <- fit[!state(fit) %in% c(3,4), ]
	ir <- RangedDataCNV(IRanges(3.6e6, 5.7e6),
			    chrom=8,
			    sampleId=sampleNames(redonSet))
	j <- subjectHits(findOverlaps(ir, fit))
	checkEquals(state(fit)[j], 5L)
}

test_hmm2 <- function(){
	load("../../data/cnSet.rda")
	library(VanillaICE)
	library(oligoClasses)
	library(Biobase)
	library(ff)
	outdir <- "/local_data/r00/rscharpf/crlmmCompendium"
	ldPath(outdir)
	invisible(open(cnSet))
	redonSet2 <- cnSet2oligoSet(cnSet, batch.name="SHELF", chrom=8)
	sampleNames(redonSet2) <- getHapMapIds(redonSet2)
	rSet <- redonSet2[[1]][, match("NA19007", sampleNames(redonSet2))]
	fit <- hmm(rSet, p.hom=0, prOutlierBAF=list(initial=1e-4, max=1e-3, maxROH=1e-3))
	rd <- fit[state(fit)==5, ]
	rd <- fit[!state(fit) %in% c(3,4), ]
	fig <- SNPchip:::xyplotLrrBaf(rd, rSet,
				      frame=200e3,
				      panel=SNPchip:::xypanelBaf,
				      scales=list(x="free"),
				      par.strip.text=list(cex=0.6),
				      cex=0.2,
				      state.col="black",
				      state.cex=0.8,
				      pch=21)
	checkEquals(state(rd)[1], 5)
}




