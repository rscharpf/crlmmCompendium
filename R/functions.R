getHapMapIds <- function(object){
	if(is(object$SKW, "ff")){
		open(object$SKW)
		open(object$SNR)
	}
	if(all(nchar(sampleNames(object)) == 7)) return(object)
	pathToSampleMap <- system.file("extdata", package="crlmmCompendium")
	sampleInfo <- read.table(file.path(pathToSampleMap, "passing_cels_sample_map.txt"),
				 sep="\t", stringsAsFactors=FALSE)
	colnames(sampleInfo) <- c("hapmapId", "celfile")
	pD <- pData(object)[,1:3]
	pD2 <- merge(pD, sampleInfo, by.y=2, by.x=0)
	rownames(pD2) <- pD2[[1]]
	pData(object) <- pD2[, -1]
	object$celFiles <- sampleNames(object)
	df <- data.frame(celFiles=sampleNames(object),
			 hapmapId=object$hapmapId,
			 stringsAsFactors=FALSE)

	hapmapId <- df$hapmapId
	hapmapId[duplicated(hapmapId)] <- paste(hapmapId[duplicated(hapmapId)], "dup", sep="_")
	hapmapId[duplicated(hapmapId)] <- paste(hapmapId[duplicated(hapmapId)], "01", sep="_")
	stopifnot(all(!duplicated(hapmapId)))
	return(hapmapId)
}

makeTransparent <- function(color, ...){
 	rgb.col <- col2rgb(color)/255
	rgb(rgb.col[1, ], rgb.col[2,], rgb.col[3,], ...)
}

prePredictPanel <- function(object, take.log=TRUE){
	b <- t(as.matrix(B(object)))
	a <- t(as.matrix(A(object)))
	gt <- t(as.matrix(calls(object)))
	gt.conf <- i2p(t(as.matrix(snpCallProbability(object))))
	fns <- featureNames(object)
	snpId <- matrix(fns, nrow(a), ncol(a), byrow=T)
	snpId <- factor(snpId, levels=fns, ordered=TRUE)  ##IMPORTANT
	if(take.log){
		df <- data.frame(A=log2(as.integer(a)), B=log2(as.integer(b)), gt=as.integer(gt), gt.conf=as.numeric(gt.conf), snp=snpId)
	} else 	df <- data.frame(A=as.integer(a), B=as.integer(b), gt=as.integer(gt),
				 gt.conf=as.numeric(gt.conf), snp=snpId)
	return(df)
}

cnPredictionPanel <- function(x, y, ...,
			      object,
			      copynumber, x.axis, line.col, line.lwd, shades, subscripts,
			      data.last=FALSE,
			      highlight.index=NULL,
			      scale.sd=rep(1,2)){
	if(length(scale.sd)==1) scale.sd <- rep(scale.sd,3)
	panel.grid(h=5, v=5)
	panel.xyplot(x, y, ...)
	if(!is.null(highlight.index)){
		##ii <- subscripts[highlight.index]
		ii <- highlight.index
		lpoints(x[ii], y[ii], pch="X", cex=1.5, col="black")
	}
	i <- panel.number()
	nuB <- as.numeric(nu(object, "B"))[i]
	phB <- as.numeric(phi(object, "B"))[i]
	nuA <- as.numeric(nu(object, "A"))[i]
	phA <- as.numeric(phi(object, "A"))[i]
	taus <- tau2(object)[, , , 1]
	cors <- corr(object)[, , 1]
	t2A <- taus[i, "A", "BB"]
	t2B <- taus[i, "B", "AA"]
	s2A <- taus[i, "A", "AA"]
	s2B <- taus[i, "B", "BB"]
	corrAB <- cors[i, "AB"]
	corrAA <- cors[i,"AA"]
	corrBB <- cors[i,"BB"]
	k <- 1
	for(CN in copynumber){
		for(CA in 0:CN){
			CB <- CN-CA
			A.scale <- sqrt(t2A*(CA==0) + s2A*(CA > 0))
			B.scale <- sqrt(t2B*(CB==0) + s2B*(CB > 0))
			if(CA == 0 | CB == 0){
				A.scale <- A.scale*scale.sd[1]
				B.scale <- B.scale*scale.sd[1]
			} else { ## both greater than zero
				A.scale <- A.scale*scale.sd[2]
				B.scale <- B.scale*scale.sd[2]
			}
			scale <- c(A.scale, B.scale)
			if(CA == 0 & CB > 0) rho <- corrBB
			if(CA > 0 & CB == 0) rho <- corrAA
			if(CA > 0 & CB > 0) rho <- corrAB
			if(CA == 0 & CB == 0) rho <- 0
			if(x.axis=="A"){
				dat.ellipse <- ellipse(x=rho, centre=c(log2(nuA+CA*phA), log2(nuB+CB*phB)), scale=scale)
			} else {
				dat.ellipse <- ellipse(x=rho, centre=c(log2(nuB+CB*phB), log2(nuA+CA*phA)), scale=rev(scale))
			}
			lpolygon(dat.ellipse[, 1], dat.ellipse[, 2], border=line.col[k], col=shades[k], ...)
		}
		k <- k+1
	}
	if(data.last) {
		panel.xyplot(x, y, ...)
		if(!is.null(highlight.index)){
			lpoints(x[ii], y[ii], pch="X", cex=1.5, col="black")
		}
	}
}

##insertCentromereBreak <- function(range.object, insertWhere=c("centromereStart", "centromereEnd")){
##	data(chromosomeAnnotation)
##	stopifnot(insertWhere %in% c("centromereStart", "centromereEnd"))
##	chrom <- unique(range.object$chrom)
##	stopifnot(all(is.integer(chrom)))
##	j <- which(colnames(chromosomeAnnotation) == insertWhere)
##	if(insertWhere=="centromereStart"){
##		index <- which(start(range.object) <= chromosomeAnnotation[chrom, 1] & end(range.object) >= chromosomeAnnotation[chrom, 1])
##		if(length(index) == 0) return(range.object)
##		##stopifnot(length(index) == 1)
##		end(range.object)[index] <- chromosomeAnnotation[chrom, 1]
##	}
##	if(insertWhere=="centromereEnd"){
##		index <- which(start(range.object) <= chromosomeAnnotation[chrom, 2] & end(range.object) >= chromosomeAnnotation[chrom, 2])
##		if(length(index) == 0) return(range.object)
##		##stopifnot(length(index) == 1)
##		start(range.object)[index] <- chromosomeAnnotation[chrom, 2]
##	}
##	range.object[index, ]
##}

my.rbind <- function(range1, range2){
	RangedData(IRanges(c(start(range1), start(range2)),
			   c(end(range1), end(range2))),
		   numMarkers=c(range1$numMarkers, range2$numMarkers),
		   seg.mean=c(range1$seg.mean, range2$seg.mean),
		   chrom=c(range1$chrom, range2$chrom))
}



addCentromereBreaks <- function(ranges.object){
##	range1 <- insertCentromereBreak(ranges.object, insertWhere="centromereStart")
##	range2 <- insertCentromereBreak(ranges.object, insertWhere="centromereEnd")
##	centromere.ranges <- my.rbind(range1, range2)
##	index <- which(start(ranges.object) == start(centromere.ranges)[1])
##	ranges.object <- ranges.object[-index, ]
##	rd <- my.rbind(ranges.object, centromere.ranges)
##	rd[sort(start(rd)), ]
	path.gap <- system.file("extdata", package="SNPchip")
	gaps <- readRDS(list.files(path.gap, pattern=paste("gap_", metadata(ranges.object)$genome, ".rda", sep=""), full.names=TRUE))
	tmp <- findOverlaps(gaps[, ], ranges.object)
	j <- subjectHits(tmp)
	i <- queryHits(tmp)
	gr <- gaps[i, ]
	gr <- ranges.object[j, ]
	start(gr) <- start(gaps)[i]
	end(gr) <- end(gaps)[i]
	cbs.segs2 <- c(ranges.object, gr)
	disjoin.gr <- disjoin(cbs.segs2)
	olaps <- findOverlaps(disjoin.gr, cbs.segs2)
	i <- queryHits(olaps); j <- subjectHits(olaps)
	elementMetadata(disjoin.gr)$seg.mean <- NA
	elementMetadata(disjoin.gr)$numberProbes <- NA
	elementMetadata(disjoin.gr)$seg.mean[i] <- values(cbs.segs2)$seg.mean[j]
	elementMetadata(disjoin.gr)$numberProbes[i] <- values(cbs.segs2)$numberProbes[j]
	## number of probes should be reestimated for the ranges that were broken by the centromere
	return(disjoin.gr)
}

lmPanel <- function(..., line.col="grey30", label.cex=0.8, Ns, allele="A", nu, ph, subscripts, ltext.y=2500){
	panel.grid(v=-1, h=0)
	panel.bwplot(...)
	i <- panel.number()
	if(allele=="A"){
		lsegments(x0=1, y0=nu[i]+2*ph[i],
			  x1=3,
			  y1=nu[i], col=line.col)
	}
	if(allele=="B"){
		lsegments(x0=1, y0=nu[i],
			  x1=3,
			  y1=nu[i]+2*ph[i], col=line.col)
	}
	ltext(1:3, y=ltext.y, labels=Ns[i, ], cex=label.cex,col="grey30")
}

##setMethod("posteriorMean", signature(object="oligoSnpSet"), function(object) assayDataElement(object, "posteriorMean"))
##setReplaceMethod("posteriorMean", signature(object="oligoSnpSet", value="matrix"), function(object, value) assayDataElementReplace(object, "posteriorMean", value))

cnPanel <- function(x, y, ..., pch.cols, gt, cbs.segs, hmm.segs, shades, subscripts, add.ideogram=TRUE){
	add.ideogram <- add.ideogram[[panel.number()]]
	##cbs.segs <- cbs.segs[[panel.number()]]
	draw.hmm.states <- ifelse(panel.number() <= length(hmm.segs), TRUE, FALSE)
	panel.grid(h=6, v=10)
	which.hom <- which(gt[subscripts] == 1 | gt[subscripts]==3)
	which.het <- which(gt[subscripts] == 2)
	panel.xyplot(x, y, col="grey60", ...)
	lpoints(x[which.hom], y[which.hom], col=pch.cols[1], ...)
	lpoints(x[which.het], y[which.het], col=pch.cols[2], ...)
	lsegments(x0=start(cbs.segs)/1e6, x1=end(cbs.segs)/1e6,
		  y0=values(cbs.segs)$seg.mean,
		  y1=values(cbs.segs)$seg.mean, ...)
	if(draw.hmm.states){
		hmm.segs <- hmm.segs[order(width(hmm.segs), decreasing=TRUE), ]
		stateindex <- as.integer(factor(state(hmm.segs),levels=c(1,2,3,5)))
		##stateindex <- as.integer(factor(fit.cn$state, levels=c(1, 5)))
		lrect(xleft=start(hmm.segs)/1e6,
		      xright=end(hmm.segs)/1e6,
		      ybottom=-0.4,
		      ytop=0,
		      border=shades[stateindex],
		      col=shades[stateindex])
	}
	ltext(130, 5, paste("MAD:", round(mad(y, na.rm=TRUE), 2)))
	if(add.ideogram){
		plotIdiogram(chromosome(hmm.segs)[1],
			     build=metadata(cbs.segs)$genome,
			     cytoband.ycoords=c(5.6, 5.9),
			     new=FALSE,
			     labels.cytoband=FALSE,
			     is.lattice=TRUE, unit="Mb")
	}
}


mySweave <- function(fname, ...){
	require(cacheSweave)
	require(tools)
	Sweave(paste(fname, ".Rnw", sep=""), ...)
	texi2dvi(paste(fname, ".tex", sep=""), pdf=TRUE)
}

run <- function() {
	require(cacheSweave)
	require(crlmmCompendium)
	mySweave("manuscript", driver=cacheSweaveDriver)
}

setMethod("show", signature(object="CNSet"), function(object){
	is.ff <- is(calls(object), "ff_matrix") | is(calls(object), "ffdf")
	if(is.ff){
		##to avoid warnings
		if("SKW" %in% varLabels(object)) {
			if(is(object$SKW, "ff"))
				open(object$SKW)
		}
		if("SNR" %in% varLabels(object)){
			if(is(object$SNR, "ff"))
				open(object$SNR)
		}
		if("gender" %in% varLabels(object)){
			if(is(object$gender, "ff"))
				open(object$gender)
		}
	}
	ad.class <- class(A(object))[1]
	cat("CNSet (assayData/batchStatistics elements: ", ad.class, ")\n", sep="")
	callNextMethod(object)
	bns <- batchNames(object)
	bns <- bns[-length(bns)]
	freq <- as.integer(table(batch(object)))
	L <- length(bns)
	if(L > 5){
		batch.freq <- paste(c(bns[1:5], "..."), c(freq[1:5], ""), collapse=", ")
		cat("batch:   ", batch.freq, "\n")
	} else cat("batch:   ", paste(bns, freq, sep=", "), "\n")
	adim <- list(nrow(object), length(batchNames(object)))
	cat("batchStatistics: ", length(ls(batchStatistics(object))), " elements, ", nrow(object), " features, ", length(unique(batch(object))), " batches\n")
})
