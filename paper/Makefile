TEX = pdflatex -interaction=nonstopmode -file-line-error
BIB = bibtex
OUTPUT = QuadIron_Library_companion_paper

.PHONY: all view

# create aux files -> create bbl file -> create pdf -> rename
all : aux bbl pdf rename

view :	# open paper
		open ${OUTPUT}.pdf

pdf :	# compile pdf
		$(TEX) quadiron.tex

bbl :	# make .bbl file for bib
		$(BIB) quadiron

aux :	# make aux file
		$(TEX) quadiron.tex

rename: # rename output file
		mv quadiron.pdf ${OUTPUT}.pdf

clean:  ## clean aux files
		rm *.aux *.log *.out
