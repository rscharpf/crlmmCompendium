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
### code chunk number 15: loadObject-snr
###################################################
if(!exists("SNR")) data(SNR)


###################################################
### code chunk number 12: snr
###################################################
(snrfig <- histogram(~SNR, breaks=100))


