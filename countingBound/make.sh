#!/bin/bash
# Builds the document.

(cd R; R CMD BATCH --no-save --quiet figs1.R )
(cd py/figure/; ./zeroing.py)
(cd py/; ./lp_brute_force_1.py 6 3 1 > lp_brute_force_1_6_3.txt)
# (cd py; ./approxRank.py)
pdflatex countingBound.tex
bibtex countingBound.aux
pdflatex countingBound.tex
pdflatex countingBound.tex
