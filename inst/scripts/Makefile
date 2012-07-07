All: codechunks fig1 fig2 fig3 fig4 fig5 sync
## R code chunks for figures
codechunks:
	echo "Stangle(\"jss664.Rnw\", split=TRUE)" | R-devel --no-save --no-restore;
	cat Figures/manuscript-compendium.R > pkgs.R
	cat pkgs.R Figures/manuscript-loadObject-snr.R Figures/manuscript-snr.R > fig1.R
	cat pkgs.R Figures/manuscript-loadObject-genotypeSet.R Figures/manuscript-dataframeForClusterPlot.R Figures/manuscript-genotypeColor.R Figures/manuscript-confidenceColor.R Figures/manuscript-scandates.R Figures/manuscript-dt.batch.R Figures/manuscript-plateColor.R Figures/manuscript-expandDataFrame.R Figures/manuscript-clusterPlot.R > fig2.R
	cat pkgs.R Figures/manuscript-loadObject-exampleData1.R Figures/manuscript-dataCnFigs.R Figures/manuscript-boxplotA.R > fig3.R
	cat pkgs.R Figures/manuscript-loadObject-exampleData1.R Figures/manuscript-defineLatticeObjects.R Figures/manuscript-ABscatterplots.R > fig4.R
	cat pkgs.R Figures/manuscript-loadObject-redonSet.R Figures/manuscript-center.R Figures/manuscript-HMM.R Figures/manuscript-viterbi.R > hmm.R
	cat hmm.R Figures/manuscript-dnacopy.R Figures/manuscript-cbs_segment.R Figures/manuscript-coerce2Iranges.R Figures/manuscript-oligo.data.frame.R Figures/manuscript-redonFigSetup.R Figures/manuscript-Fig5.R > fig5.R
	cat pkgs.R Figures/manuscript-loadObject-snr.R Figures/manuscript-snr.R Figures/manuscript-loadObject-genotypeSet.R Figures/manuscript-dataframeForClusterPlot.R Figures/manuscript-genotypeColor.R Figures/manuscript-confidenceColor.R Figures/manuscript-loadObject-exampleData1.R Figures/manuscript-dataCnFigs.R Figures/manuscript-boxplotA.R Figures/manuscript-defineLatticeObjects.R Figures/manuscript-ABscatterplots.R Figures/manuscript-loadObject-redonSet.R Figures/manuscript-center.R Figures/manuscript-HMM.R Figures/manuscript-viterbi.R > standalone.R
	-rm hmm.R
	-rm Figures/manuscript-*.R


fig1: fig1.R
	cp Figures/manuscript-snrFig.pdf ../doc/org/img/fig1.pdf
	convert ../doc/org/img/fig1.pdf  ../doc/org/img/fig1.png
	convert -thumbnail 250 ../doc/org/img/fig1.pdf ../doc/org/img/thumb_fig1.png

fig2: Figures/manuscript-clusterColors.pdf
	cp Figures/manuscript-clusterColors.pdf ../doc/org/img/fig2.pdf
	convert ../doc/org/img/fig2.pdf  ../doc/org/img/fig2.png
	convert -thumbnail 250 ../doc/org/img/fig2.pdf ../doc/org/img/thumb_fig2.png

fig3: Figures/manuscript-boxPlotA.pdf
	cp Figures/manuscript-boxplotA.pdf ../doc/org/img/fig3.pdf
	convert ../doc/org/img/fig3.pdf  ../doc/org/img/fig3.png
	convert -thumbnail 250 ../doc/org/img/fig3.pdf ../doc/org/img/thumb_fig3.png

fig4: Figures/manuscript-ABscatterplots.pdf
	cp Figures/manuscript-ABscatterplots.pdf ../doc/org/img/fig4.pdf
	convert ../doc/org/img/fig4.pdf  ../doc/org/img/fig4.png
	convert -thumbnail 250 ../doc/org/img/fig4.pdf ../doc/org/img/thumb_fig4.png

fig5: Figures/manuscript-redonFigure.pdf
	cp Figures/manuscript-redonFigure.pdf ../doc/org/img/fig5.pdf
	convert ../doc/org/img/fig5.pdf  ../doc/org/img/fig5.png
	convert -thumbnail 250 ../doc/org/img/fig5.pdf ../doc/org/img/thumb_fig5.png

sync:
	rsync -avuzb ../doc/org/img/ ~/Public/public_html/crlmmCompendium/img
	cp fig[1-5].R ~/Public/public_html/crlmmCompendium/
	cp ~/Papers/2010/crlmmCompendium_1.1.0.tar.gz ~/Public/public_html/crlmmCompendium/
	cp jss664.Rnw ~/Public/public_html/crlmmCompendium

bib: jss664.aux
	bibtool -x jss664.aux -o jss664.bib -- 'expand.macros = ON'
	rm jss664.bbl jss664.blg
	## I had to go here
	##	http://blog.mikezhang.com/dcolumn/bibcapital.cgi
	## to capitalize the titles of the *.bib file

clean:
	-rm -rf manuscript.tex manuscript.bbl manuscript.blg manuscript.log manuscript.out manuscript.Rout *.sh.* .Rhistory .RData Rplots.pdf auto* *~ *.sh






