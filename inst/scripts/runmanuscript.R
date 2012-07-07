library(cacheSweave)
if(FALSE){ ## release
	outdir <- "/nexsan2/disk3/rscharpf/crlmmCompendium/release"
} else { ## devel
	outdir <- "/local_data/r00/rscharpf/crlmmCompendium"
}
suppressWarnings(dir.create(outdir))
setCacheDir(outdir)
library(VanillaICE)
library(crlmmCompendium)
mySweave("~/Software/crlmmCompendium/inst/scripts/jss664", driver=cacheSweaveDriver)
crlmmCompendium:::run()##makes the pdf

##save(cnSet,file="~/cnSet_hapmapIII.rda")
