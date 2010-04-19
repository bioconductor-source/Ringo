all:	pdf

pdf:	Ringo.tex
	pdflatex Ringo
	bibtex Ringo
	pdflatex Ringo
	pdflatex Ringo
	make clean

clean:
	rm -f *.out *.bbl *.aux *.log *.blg *.brf *.toc  *~ *.eps
	rm -f *.nfs*
	rm -f Ringo-*.jpg Ringo-*.png Ringo-*.pdf
#	rm -f Rplots.pdf Rplots.ps # problem on build machine 'gewurz'
#	rm -f *.dvi *.ps Ringo.tex
