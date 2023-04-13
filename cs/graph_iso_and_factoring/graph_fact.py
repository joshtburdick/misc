#!/usr/bin/env python3
# Checking a construction to reduce factoring to GI.
# (I'm guessing this won't work, but it seems simplest to just
# write a program to check some examples, rather than drawing
# graphs with pencil and paper.)

import pdb

import networkx
import networkx.algorithms.isomorphism.vf2pp

class GraphIsoFactor:
    """Constructs graphs to compare for factoring.

    """

    def __init__(self, primes):
        """Constructor.

        primes: list of primes to use
        """
        self.primes = primes
        # precompute graph for 1
        self.g1 = self.get_graph(1)

    def get_graph(self, x):
        """Makes a graph for x.

        x: the number to construct a graph for
        """
        # initialize the graph
        g = networkx.Graph()
        for i in range(len(self.primes)-1):
            for j in range(i+1, len(self.primes)):
                p1 = self.primes[i]
                p2 = self.primes[j]
                # add edges at each "end"
                self.add_end(g, p1, p2, 'A')
                self.add_end(g, p1, p2, 'B')
                # add "middle" edges
                self.add_middle(g, p1, p2, x)
        return g

    def add_end(self, g, p1, p2, end):
        """Adds on 'end-to-middle' edges.

        g: the graph to add to
        p1, p2: the primes to add an "end" for
        end: the "end" to add on (either 'A' or 'B')
        Side effects: adds 'end-to-middle' edges to g
        """
        m = p1 * p2
        for i in range(1, m):
            x1 = i % p1
            x2 = i % p2
            if (x1 > 0) and (x2 > 0):
                g.add_edge((end, p1, x1), (end, p1, p2, i))
                g.add_edge((end, p2, x2), (end, p1, p2, i))

    def add_end_CBIP(self, g, p1, p2, end):
        """Adds on a CBIP of end edges.

        g: the graph to add to
        p, p1, p2: the primes to add an "end" for
        end: the "end" to add on (either 'A' or 'B')
        Side effects: adds a CBIP of edges to g, of vertices
            which are equivalent mod p.
        """
        pass
        # FIXME


    def add_middle(self, g, p1, p2, x):
        """Adds the 'middle' connections.

        This adds connections between A and B, where
            a*b = x (mod p1*p2).
        g: the graph to add on to
        p1, p2: the primes giving the pair of edges to connect
            (with p1 < p2)
        x: the number to factor (_only_ mod p1*p2)
        Side effects: adds edges described above.
        """
        m = p1 * p2
        # we assume that x1 is relatively prime to p1 and p2
        x1 = x % m
        # loop through the possible factors (mod p1*p2)
        for a in range(1, m):
            for b in range(1, m):
                # XXX not worrying much about efficiency here
                # we assumed x1 is relatively prime to p1 and p2, so
                # I think we don't need to check that here
                if a*b % m == x1:
                    g.add_edge(('A', p1, p2, a), ('B', p1, p2, b))

    def print_factors_from_iso(self, iso):
        """Prints the factors implied by an isomorphism."""
        for end in ['A', 'B']:
            print(end + ':')
            for p in self.primes:
                y = iso[(end, p, 1)]
                print(f'{y} (mod {p})')
        print('------')

    def print_isos(self, x):
        """Print all the isomorphisms."""
        g = gif.get_graph(x)
        count = 0
        for iso in networkx.algorithms.isomorphism.vf2pp.vf2pp_all_isomorphisms(self.g1, g):
            # for now, just printing counts
            count += 1
#            if count % 1 == 0:
#                print(count)
        print(count)
            # self.print_factors_from_iso(iso)
            # print()

if __name__ == '__main__':
    gif = GraphIsoFactor([3,7])
    gif.print_isos(10)

# factors   # isos
# 3,5       96 = 2^5 * 3 = 4 * 4!
# 3,7       2880 = 2^6 * 3^2 * 14 = 4! * 5!
# 5,7       dunno, it's slooooow...
# summary: don't know what, if any, pattern this is

