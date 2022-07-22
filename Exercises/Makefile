SRC = $(wildcard Exercise_*.Rmd)

PDF   = $(SRC:.Rmd=.pdf)
HTML  = $(SRC:.Rmd=_question.html)
SOL  = $(SRC:.Rmd=_solution.html)
R     = $(SRC:.Rmd=.R)

RENDER_H = @Rscript -e "rmarkdown::render('$<', 'html_document', gsub('.Rmd', '_question.html', gsub('__','Q_','$<')), params=list(solutions=FALSE, latex=FALSE))"
RENDER_S = @Rscript -e "rmarkdown::render('$<', 'html_document', gsub('.Rmd', '_solution.html', gsub('__','S_','$<')), params=list(solutions=TRUE, latex=FALSE))"
RENDER_D = @Rscript -e "rmarkdown::render('$<', 'pdf_document', gsub('.Rmd', '.pdf', gsub('__','Q_','$<')), params=list(solutions=FALSE, latex=TRUE))"
RENDER_B = @Rscript -e "rmarkdown::render('$<', 'all')"
PURL = @Rscript -e "knitr::purl('$<', documentation = 2L, output = paste0(tools::file_path_sans_ext('$<'), '.R'))"

%.R:%.Rmd
	$(PURL)
%_question.html:%.Rmd
	$(RENDER_H)
	-rm -rf Exercise_*.log
%_solution.html:%.Rmd
	$(RENDER_S)
	-rm -rf Exercise_*.log
%.pdf:%.Rmd
	$(RENDER_D)
	-rm -rf Exercise_*.log

.PHONY: clean
.PHONY: tidy
.PHONY: r
.PHONY: pdf
.PHONY: html
.PHONY: all
	
all: 	$(PDF) $(HTML) $(SOL)
pdf:	$(PDF)
html:	$(HTML) $(SOL)
r: $(R)
clean:
	-rm -rf Exercise_*.md
	-rm -rf Exercise_*.tex
	-rm -rf Exercise_*.pdf
	-rm -rf Exercise_*.html
	-rm -rf Exercise_*.R
	-rm -rf Exercise_*.log
	-rm -rf Exercise_*_files
tidy:
	-rm -rf Exercise_*.md
	-rm -rf Exercise_*.tex
	-rm -rf Exercise_*.log
	-rm -rf Exercise_*_files
