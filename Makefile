SRC = $(wildcard Session_*/Session_*.Rmd)

PDF   = $(SRC:.Rmd=.pdf)
HTML  = $(SRC:.Rmd=.html)
R     = $(SRC:.Rmd=.R)

RENDER_H = @Rscript -e "rmarkdown::render('$<', 'html_document', params=list(presentation=FALSE))"
RENDER_P = @Rscript -e "rmarkdown::render('$<', 'beamer_presentation', params=list(presentation=TRUE))"
RENDER_D = @Rscript -e "rmarkdown::render('$<', 'pdf_document', params=list(presentation=FALSE))"
RENDER_B = @Rscript -e "rmarkdown::render('$<', 'all')"
PURL = @Rscript -e "knitr::purl('$<', documentation = 2L, output = paste0(tools::file_path_sans_ext('$<'), '.R'))"

%.R:%.Rmd
	$(PURL)
%.html:%.Rmd
	$(RENDER_H)
	-rm -rf Session*/Session*.log
%Session_Preparation.pdf:%Session_Preparation.Rmd
	$(RENDER_D)
	-rm -rf Session*/Session*.log
%.pdf:%.Rmd
	$(RENDER_P)
	-rm -rf Session*/Session*.log

.PHONY: clean
.PHONY: tidy
.PHONY: r
.PHONY: pdf
.PHONY: html
.PHONY: all
	
all: 	$(PDF) $(HTML) $(R)
pdf:	$(PDF)
html:	$(HTML)
r: $(R)
clean:
	-rm -rf Session*/Session*.md
	-rm -rf Session*/Session*.tex
	-rm -rf Session*/Session*.pdf
	-rm -rf Session*/Session*.html
	-rm -rf Session*/Session*.R
	-rm -rf Session*/Session*.log
	-rm -rf Session*/Session*_files
tidy:
	-rm -rf Session*/Session*.md
	-rm -rf Session*/Session*.tex
	-rm -rf Session*/Session*.log
	-rm -rf Session*/Session*_files
