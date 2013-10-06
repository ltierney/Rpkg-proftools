
.SUFFIXES: .nw .tex .pdf .html .bib

.nw.tex:
	noweave -delay -index $< > $@
.tex.pdf:
	pdflatex $<

.nw.html:
	noweave -filter l2h -index -html $< > $@

prof.R: proftools.nw
	notangle -R$@ $< > $@

clean:
	rm -f prof.R proftools.dvi proftools.pdf proftools.html \
	proftools.aux proftools.log proftools.out proftools.tex
