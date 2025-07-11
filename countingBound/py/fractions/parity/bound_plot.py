#!/usr/bin/env python3
# Plots the bounds.

import os
import pdb

import matplotlib.pyplot as plt
import numpy as np
import pandas
import seaborn as sns
from matplotlib.ticker import MaxNLocator

# output_dir = 'bound_plot/'
output_dir = '../../../bound2/'

def plot_bounds(bound_file, output_file):
    """Plots bounds from one file."""
    b = pandas.read_csv(bound_file)
    # b = b[b['Num. vertices']==6]
    # omitting "only step" bound
    b = b[ b['Constraints'] != 'Step' ]
    if False:
        # plot as separate lines
        g = sns.FacetGrid(b, col='Constraints',
            sharey=False, ylim=(0,25))
#            facet_kws=dict(sharey=False))
        g.map(sns.lineplot, 'Num. cliques', 'Min. gates',
            color='black')
    else:
        plt.figure(figsize=(6,3.5))
        # or "overlapping lines version" of this
        g = sns.lineplot(b, x='Num. cliques', y='Min. gates',
            hue='Constraints', alpha=0.5, lw=4)
#            palette=["#ff404050", "#40ff4050", "#202020a0"])
#        plt.xlim(0-0.2, 20+0.2)
#        plt.gca().set_xticks([0,5,10,15,20])
        # plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    #   other color possibilites:
    #        palette="YlGnBu")
    #        palette=["#a04040e0", "#40a000e0", "#20202080"])
    sns.move_legend(plt.gca(), "center")
    plt.tight_layout()
    plt.savefig(output_file)

os.makedirs(output_dir, exist_ok=True)
for s in ['ip_parity_nand_8_3_9']:
    plot_bounds(f'bounds/{s}.csv', f'{output_dir}/{s}.pdf')

