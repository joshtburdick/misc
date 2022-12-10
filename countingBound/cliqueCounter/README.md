
This is code to count the sizes of cliques in random hypergraphs
(and how many hyperedges they cover).

The basic strategy is to use a bitmask to represent the hyperedges of
a graph, and precompute bitmasks for all the possible hypercliques.
A random bitstring corresponds to a random graph, which can be quickly
compared with the table of possible cliques.

One downside is that the table of cliques uses a lot of memory.
Also, the random cliques aren't expected to be all that large.
Therefore, it might have made more sense to just use, e.g.,
depth-first search.

This is not the cleanest C++ code ever. (It sometimes follows the
Google C++ guidelines, but not consistently.)

