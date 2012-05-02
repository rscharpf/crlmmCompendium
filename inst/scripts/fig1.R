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
### code chunk number 16: loadObject-snr
###################################################
if(!exists("SNR")) data(SNR)


###################################################
### code chunk number 13: snr
###################################################
(snrfig <- histogram(~SNR, breaks=100))


