#!/bin/bash
# Makes some plots.

mkdir -p bounds

./ip_layers.py 8 3 90 1,2,3,5,10,25 --result-file bounds/ip_layers_7_3_90.csv
./ip_layers_plot.py

