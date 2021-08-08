#!/bin/bash
# Builds the document.

(cd R; R CMD BATCH --no-save --quiet < figs1.R )
pdflatex countingBound.tex
bibtex countingBound.aux
pdflatex countingBound.tex

