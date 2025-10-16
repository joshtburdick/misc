#!/usr/bin/env python3
# Attempt at bounding rank, based on zeroing out edges.
# FIXME
# - use argparse?
# - add sharper upper bounds, based on number of vertices zeroed out?
#   (??? is this still needed? add_vertex_zeroing_constraints()
#   implements this, but doesn't seem to have any effect)

import argparse
import pdb
import sys

import fractions
import itertools
import numpy as np
import scipy.special
import scipy.stats

import pandas

sys.path.append("../")            # XXX
sys.path.append("../fractions/")  # XXX

import lp_helper
import pulp_helper

# Wrapper for comb(), with exact arithmetic.
def comb(n, k):
    return scipy.special.comb(n, k, exact=True)

# Hypergeometric distribution, with exact arithmetic.
def hyperg_frac(N, K, n, k):
    # based on https://en.wikipedia.org/wiki/Hypergeometric_distribution
    # note that we don't try to optimize this
    return fractions.Fraction(
        comb(K, k) * comb(N-K, n-k),
        comb(N, n))

class LpEdgeZeroing:
    """Attempt at bound by zeroing out vertices.
    """

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        # the number of possible cliques
        self.max_cliques = comb(n, k)
        # the number of cliques which could be "hit" by the zonked edge
        self.max_hit_by_edge = comb(n-2, k-2)
        # number of cliques, omitting those which contain the zonked edge
        self.max_cliques_low = self.max_cliques - self.max_hit_by_edge
        # the total number of functions
        self.num_functions = 2 ** self.max_cliques
        # variable names: the expected ("E") rank at each level,
        # and in the "high" and "low" sets at each level
        vars = ([("E", i) for i in range(self.max_cliques+1)
            + [("high", i) for i in range(self.max_cliques+1)
            + [("low", i) for i in range(self.max_cliques_low+1)])
        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(vars)



        # ??? omit all of these?
        # to count functions in "low" sets, omit cliques hit by an edge
        # (XXX currently we pad this with 0s at the top)
        self.num_functions_low = [comb(self.max_cliques_low, i) 
            for i in range(self.max_cliques+1)
        ]
        # to count functions in "high" sets, count all the
        # possible functions, and subtract off functions in "low" sets
        self.num_functions_high = [
            comb(self.max_cliques, i) - self.num_functions_low[i]
            for i in range(self.max_cliques+1)
        ]


    def add_expected_level_constraints(self):
        """Adds constraints for expected rank, at each level."""




    def add_average_rank_constraints(self):
        """Adds equality constraints on average rank of all of the functions.

        This is an average of the expected rank, for each set of functions,
        weighted by that set's size.
        """
        A = ([("high", i): self.num_functions_high[i]
            for i in range(self.max_cliques+1)]
            + [("low", i): self.num_functions_low[i]
            for i in range(self.max_cliques_low+1)])
        self.lp.add_constraint(A, "=", (self.num_functions-1) / 2)

    def add_edge_zeroing_constraints(self):





    def add_vertex_zeroing_constraints(self):
        """Adds constraints from zeroing out vertices.

        The sets of cliques are:
        C: the cliques in the larger set, using (up to) k+1 vertices
        B: the cliques zeroed by feeding in zeros to a vertex of C
        A: the cliques which are left over
        """
        # add case when v == self.k
        # FIXME: not sure this is sufficient...
        self.lp.add_constraint([((self.k, 1), 1)],
            '>=',
            comb(self.n, self.k)/2 + 1/2)

        # loop through number of vertices in larger graph
        for v in range(self.k+1, self.n+1):
            # add "zero cliques" constraint
            self.lp.add_constraint([((v, 0), 1)], "=", 0)
            # loop through number of cliques in that graph
            # (note that the trivial bound when C_size==0 is implied by the
            # overall "box" constraints on all the variables)
            for C_size in range(1, comb(v, k)+1):
                # Maximum number of cliques we might hit,
                # _assuming that only v vertices are present_
                # (If we pick a vertex randomly, and miss all of the
                # cliques, we get to "re-roll".)
                num_cliques_hitting_vertex = comb(v-1, k-1)
                # Bounds on number of cliques zeroed.
                # The min is at least 1 (assuming we can "re-roll"), and no
                # more than the difference between the number of cliques in the
                # larger set, and the number left over.
                min_zeroed = max(1, C_size - comb(v-1, k))
                # The max is limited by the current number of vertices (and how
                # many cliques hit one vertex), and the number of cliques.
                max_zeroed = min(num_cliques_hitting_vertex, C_size)
                # the range of possible number of cliques zeroed ...
                B_size = np.arange(min_zeroed, max_zeroed+1)
                # ... and the number left over
                A_size = C_size - B_size
                # the probability of some number of cliques being hit
                # (again, assuming that only v vertices are "in use" by
                # the hyperedges)
                # was:
                # p_zonk = scipy.stats.hypergeom(comb(v, k),
                #     C_size,
                #     num_cliques_hitting_vertex)
                def p_zonk(x):
                    print(f"v={v} k={k} C_size={C_size} num_={num_cliques_hitting_vertex} x={x}")
                    return hyperg_frac(comb(v, k),
                        C_size,
                        num_cliques_hitting_vertex,
                        x)
                # the probability of at least one clique being hit
                # (if this happens, we get to "re-roll", so we assume it doesn't)
                # Note that since this we assume at least one clique is hit,
                # we normalize by this probability.
                p_at_least_one_hit = 1 - p_zonk(0)
                # ??? does this give the same answer?
                # p_at_least_one_hit = np.array([p_zonk.pmf(z1) for z1 in z]).sum()

                # The constraint on the rank of the v-vertex set in C (which is "hit")
                # is a weighted sum of:
                # - the expected rank of what's left in A, after zonking,
                # - half the expected number of functions in A
                # - plus half the number of functions in B
                A = [((v, C_size), 1)]
                print(f"B_size={B_size}")
                A += [((v-1, C_size-j), -p_zonk(j) / p_at_least_one_hit)
                    for j in B_size]
                # note that this is the _expected_ number of functions in A
                # FIXME rename this?
                p_zonk_1 = np.array([p_zonk(b) for b in B_size])
                A_num_functions = ( p_zonk_1 * self.counts_max_vertices[v-1][A_size] ).sum() / p_at_least_one_hit
                # the number of functions (or, "sets of cliques") in B
                # ??? is this right?
                B_num_functions = self.counts_max_vertices[v][ C_size ]
                self.lp.add_constraint(A, '>=',
                    # Since we presumably "hit" a clique, the number of "new"
                    # functions is the number which include all the vertices.
                    (A_num_functions + B_num_functions) / 2.)

    def get_bounds_on_average(self, num_levels_below):
        """Gets bounds, on average rank for sets with lots of cliques.

        num_levels_below: minimize the weighted average of from
            [N-num_levels_below, N] cliques. (Thus, when num_levels_below==0,
            this minimizes the rank of finding _all_ the cliques.)
        Returns: Pandas DataFrame of results, and the value of the objective functions
        """


        # construct objective function, as weighted average of "high levels"
        N = self.num_cliques
        num_cliques = range(self.num_cliques-num_levels_below, self.num_cliques+1)
        num_functions = [comb(N, c) for c in num_cliques]
        total_functions = sum(num_functions)
        objective_function = {
            (self.n, num_cliques[i]): num_functions[i]/total_functions
            for i in range(len(num_cliques))
        }
        x = self.lp.solve_with_objective(objective_function)
        # get bounds for all the vertices, and any number of cliques
        bounds = [{
                "Num. cliques": c,
                "Expected rank": x[(self.n, c)]
            }
            for c in range(self.num_cliques+1)
        ]
        # XXX also include the value of the objective function
        bounds += [{"Num. cliques": -1, "Expected rank": x["__objective__"]
        }]

        return pandas.DataFrame(bounds)







if __name__ == '__main__':
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    num_levels_below = int(sys.argv[3])
    bound = LpEdgeZeroing(n, k)

    # FIXME add constraints based on options?
    # ??? do we need the "base case" of only one clique?
    bound.add_average_rank_constraints()
    bound.add_edge_zeroing_constraints()

    # get bound
    b = bound.get_bounds_on_average(num_levels_below)
    print(b)

