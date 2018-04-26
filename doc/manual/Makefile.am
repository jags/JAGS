SUBDIRS = figures

SOURCES = jags_user_manual.tex jags_user_manual.bib \
	jags_installation_manual.tex jags_developer_manual.tex \
	version.in 

FIGURES = figures/halft.pdf

LATEX = pdflatex
LATEXOPT = --shell-escape
NONSTOP = --interaction=nonstopmode

LATEXMK = latexmk
LATEXMKOPT = -pdf

docs: $(SOURCES) $(FIGURES) \
      jags_user_manual.pdf \
      jags_installation_manual.pdf \
      jags_developer_manual.pdf

.tex.pdf: 
	TEXINPUTS="$(srcdir):$$TEXINPUTS" \
	BIBINPUTS="$(srcdir):$$BIBINPUTS" \
        $(LATEXMK) $(LATEXMKOPT) $<

clean-local:
	rm -f jags_*_manual.{log,aux,dvi,bbl,blg,toc,pdf,fdb_latexmk,fls}

.PHONY: docs .tex.pdf clean-local distclean-local
