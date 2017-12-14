SUBDIRS = figures

EXTRA_DIST = jags_user_manual.tex jags_user_manual.bib	\
jags_installation_manual.tex jags_developer_manual.tex

docs: jags_user_manual.pdf jags_installation_manual.pdf jags_developer_manual.pdf

jags_user_manual.pdf: jags_user_manual.bbl jags_user_manual.tex
	TEXINPUTS="$(srcdir):$$TEXINPUTS" pdflatex jags_user_manual.tex

jags_user_manual.bbl: jags_user_manual.bib
	TEXINPUTS="$(srcdir):$$TEXINPUTS" pdflatex jags_user_manual.tex
	BIBINPUTS="$(srcdir):$$BIBINPUTS" bibtex jags_user_manual
	TEXINPUTS="$(srcdir):$$TEXINPUTS" pdflatex jags_user_manual.tex

jags_developer_manual.pdf: jags_developer_manual.tex
	TEXINPUTS="$(srcdir):$$TEXINPUTS" pdflatex jags_developer_manual.tex
	TEXINPUTS="$(srcdir):$$TEXINPUTS" pdflatex jags_developer_manual.tex

jags_installation_manual.pdf: jags_installation_manual.tex
	TEXINPUTS="$(srcdir):$$TEXINPUTS" pdflatex jags_installation_manual.tex
	TEXINPUTS="$(srcdir):$$TEXINPUTS" pdflatex jags_installation_manual.tex

clean-local:
	rm -f jags_*_manual.{log,aux,dvi,bbl,blg,toc,pdf}

