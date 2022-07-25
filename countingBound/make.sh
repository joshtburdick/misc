#!/bin/bash
# Builds the document.

(cd R; R CMD BATCH --no-save --quiet figs1.R )
(cd py/figure/; ./zeroing.py)
# (cd py; ./approxRank.py)
pdflatex countingBound.tex
bibtex countingBound.aux
pdflatex countingBound.tex
pdflatex countingBound.tex
