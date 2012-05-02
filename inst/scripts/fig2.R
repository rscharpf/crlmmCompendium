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
### code chunk number 23: loadObject-genotypeSet
###################################################
if(!exists("genotypeSet")) data(genotypeSet)


###################################################
### code chunk number 24: dataframeForClusterPlot
###################################################
df <- prePredictPanel(genotypeSet)


###################################################
### code chunk number 25: genotypeColor
###################################################
fill1 <- brewer.pal(3, "Set1")[df$gt]


###################################################
### code chunk number 26: confidenceColor
###################################################
gt.conf <- df$gt.conf
min.conf <- min(gt.conf)
max.conf <- max(gt.conf)
sc <- (gt.conf - min.conf)/(max.conf-min.conf)
fill2 <- sapply(sc, grey)


###################################################
### code chunk number 27: scandates
###################################################
dt <- strftime(protocolData(genotypeSet)$ScanDate, "%Y-%m-%d", usetz=FALSE)
range(dt)


###################################################
### code chunk number 28: dt.batch
###################################################
dt.batch <- split(dt, batch(genotypeSet))
sapply(dt.batch, range)


###################################################
### code chunk number 29: plateColor
###################################################
batch.scale <- which(batch(genotypeSet)=="SCALE")
batch.sloth <- which(batch(genotypeSet)=="SLOTH")
plate.cols <- brewer.pal(8, "Accent")[c(3, 8)]
fill3 <- rep("white", nrow(df))
fill3[batch.scale] <- plate.cols[1]
fill3[batch.sloth] <- plate.cols[2]


###################################################
### code chunk number 30: expandDataFrame
###################################################
df2 <- rbind(df, df, df)
df2$fill <- c(fill1, fill2, fill3)
colorby <- c("genotype", "confidence score", "plate")
df2$colorby <- factor(rep(colorby, each=nrow(df)), levels=colorby, ordered=TRUE)


###################################################
### code chunk number 31: clusterPlot
###################################################
(ABfig <- xyplot(A~B|colorby, df2,
	      panel=function(x, y, col, fill, plate.cols, ..., subscripts){
		      panel.grid(h=5,v=5)
		      panel.xyplot(x,y, col="grey60", fill=fill[subscripts], ...)
		      if(panel.number() == 3){
			      ## such that these plates are plotted last
			      lpoints(x[batch.scale], y[batch.scale], fill=plate.cols[1], ...)
			      lpoints(x[batch.sloth], y[batch.sloth], fill=plate.cols[2], ...)
		      }
	      },
	       aspect="iso",
	       fill=df2$fill, col=df2$col, cex=0.6, pch=21, plate.cols=plate.cols,
	       xlab=expression(log[2](I[B])),
	       ylab=expression(log[2](I[A])), main=featureNames(genotypeSet), layout=c(3,1),
	       par.strip.text=list(lines=0.9, cex=0.6)))


