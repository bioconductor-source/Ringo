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
	rm -f Ringo-*.jpg Ringo-*.pdf Rplots.ps Rplots.pdf
#	rm -f *.dvi *.ps Ringo.tex
