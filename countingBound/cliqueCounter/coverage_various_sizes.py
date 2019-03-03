#!/usr/bin/python
# Runs ./coverage for various parameter sizes.

import subprocess

def run_coverage(numVertices, edgeSize, numSamples=1000):
    """Runs coverage with various parameter settings.

    Saves results to a file named by parameters
    numVertices, edgeSize: problem size
    numSamples: how many (pairs of samples) to count
    Side effects: writes a table of results"""
    # hopefully this is searching for large enough cliques
    maxCliqueSize = min(7, numVertices)
    outputFile = ('coverage_' + str(numVertices) + '_' +
        str(edgeSize) + '_' + str(maxCliqueSize) + '_' +
        str(numSamples) + '.txt')
    print(outputFile)
    with open(outputFile, 'w') as output:
        subprocess.call(['./coverage', str(numVertices),
            str(edgeSize), str(maxCliqueSize),
            str(numSamples)], stdout=output)

for edgeSize in range(2,5):
    for numVertices in range(7,26):
        run_coverage(numVertices, edgeSize)

