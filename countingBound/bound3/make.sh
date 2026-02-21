#!/bin/bash
# Builds the document.

# omitting this for now, since these keep re-running,
# and showing up in the git diffs
# (these may get rewritten in Python anyway)
# (cd R; R CMD BATCH --no-save --quiet figs1.R )

pdflatex countingBound3.tex
bibtex countingBound3.aux
pdflatex countingBound3.tex
pdflatex countingBound3.tex

