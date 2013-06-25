
.SUFFIXES: .nw .tex .pdf .html .bib

.nw.tex:
	noweave -delay -index $< > $@
.tex.pdf:
	pdflatex $<
	# -bibtex $*
	# pdflatex $<
	# pdflatex $<
.nw.html:
	noweave -filter l2h -index -html $< > $@

prof.R: proftools.nw
	notangle -R$@ $< > $@
