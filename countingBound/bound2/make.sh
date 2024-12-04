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
#    mkdir -p py/bounds
    # XXX don't know which version of this to use
#    (cd py; ./lp_gate_bound_4.py 10 5 > bounds/bounds_10_5.csv)
    # (cd py; ./lp_gate_bound_6.py 7 3 5 > bounds/bounds_7_3.csv)
#    (cd py; ./boundsTable.py)
#    (cd py; ./bound_plot.py)
fi

pdflatex countingBound2.tex
bibtex countingBound2.aux
pdflatex countingBound2.tex
pdflatex countingBound2.tex

