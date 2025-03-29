#!/bin/bash
# Builds the document.

# omitting this for now, since these keep re-running,
# and showing up in the git diffs
# (these may get rewritten in Python anyway)
# (cd R; R CMD BATCH --no-save --quiet figs1.R )

if [ $1 == "all" ]; then
    # XXX hideous directory-changing hack
    cd ../py/figure/
    # now writing figures in this directory
    # (for now, figures have been generated already)
    ./zeroing.py
    cd ../../bound2 
    mkdir -p ../py/fractions/bounds
    pwd
    (cd ../py/fractions; ./ip_bound_2.py 7 3 3 --result-file bounds/bounds_7_3_3.csv; ./bound_plot.py)
fi

pdflatex countingBound2.tex
bibtex countingBound2.aux
pdflatex countingBound2.tex
pdflatex countingBound2.tex

