SRC = $(wildcard Exercise_*.Rmd)

PDF   = $(SRC:.Rmd=.pdf)
HTML  = $(SRC:.Rmd=.html)
SOL  = $(SRC:.Rmd=_solution.html)
R     = $(SRC:.Rmd=.R)

RENDER_H = @Rscript -e "f1 <- list.files(); rmarkdown::render('$<', 'html_document', gsub('.Rmd', '.html', gsub('__','Q_','$<')), params=list(solutions=FALSE, latex=FALSE)); f2 <- list.files(); file.copy(f2[!f2 %in% f1], '../', overwrite=TRUE)"
RENDER_S = @Rscript -e "f1 <- list.files(); rmarkdown::render('$<', 'html_document', gsub('.Rmd', '_solution.html', gsub('__','S_','$<')), params=list(solutions=TRUE, latex=FALSE)); f2 <- list.files(); file.copy(f2[!f2 %in% f1], '../Solutions/', overwrite=TRUE)"
RENDER_D = @Rscript -e "f1 <- list.files(); rmarkdown::render('$<', 'pdf_document', gsub('.Rmd', '.pdf', gsub('__','Q_','$<')), params=list(solutions=FALSE, latex=TRUE)); f2 <- list.files(); file.copy(f2[!f2 %in% f1], '../', overwrite=TRUE)"
RENDER_B = @Rscript -e "rmarkdown::render('$<', 'all')"
PURL = @Rscript -e "knitr::purl('$<', documentation = 2L, output = paste0(tools::file_path_sans_ext('$<'), '.R'))"

# TODO: new .PHONY compare which deletes any files under ignore_this_folder which are not identical to files in the appropriate place

%.R:%.Rmd
	$(PURL)
%.html:%.Rmd
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
	-rm -rf ../Exercise_*.pdf
	-rm -rf ../Exercise_*.html
	-rm -rf ../Solutions/Exercise_*.html
tidy:
	-rm -rf Exercise_*.md
	-rm -rf Exercise_*.tex
	-rm -rf Exercise_*.log
	-rm -rf Exercise_*_files
	-rm -rf ../Solutions/Exercise_0_Preparation_solution.html
