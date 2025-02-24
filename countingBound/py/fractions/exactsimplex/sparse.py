'''Sparse version of simplex LP solver.

A row will be represented by a Dict<Int, Fraction>, as it were.
The tableau will be represented by a Dict<Int, Dict<Int, Fraction>>.

This works in some cases, but not others.

It sometimes seems to fail on cases in which phase 1 isn't necessary.
Also, the objective function isn't always increasing (or decreasing),
which suggests a bug.
'''

import copy
import functools
import heapq
import pdb

'''
   Return a rectangular identity matrix with the specified diagonal entries, possibly
   starting in the middle.
'''
def identity(numRows, numCols, val=1, rowStart=0):
   return [[(val if i == j else 0) for j in range(numCols)]
               for i in range(rowStart, numRows)]


'''
   standardForm: [float], [[float]], [float], [[float]], [float], [[float]], [float] -> [float], [[float]], [float]
   Convert a linear program in general form to the standard form for the
   simplex algorithm.  The inputs are assumed to have the correct dimensions: cost
   is a length n list, greaterThans is an n-by-m matrix, gtThreshold is a vector
   of length m, with the same pattern holding for the remaining inputs. No
   dimension errors are caught, and we assume there are no unrestricted variables.
'''
def standardForm(cost, greaterThans=[], gtThreshold=[], lessThans=[], ltThreshold=[],
                equalities=[], eqThreshold=[], maximization=True):
   newVars = 0
   numRows = 0
   if gtThreshold != []:
      newVars += len(gtThreshold)
      numRows += len(gtThreshold)
   if ltThreshold != []:
      newVars += len(ltThreshold)
      numRows += len(ltThreshold)
   if eqThreshold != []:
      numRows += len(eqThreshold)

   if not maximization:
      cost = [-x for x in cost]

   if newVars == 0:
      return cost, equalities, eqThreshold

   newCost = list(cost) + [0] * newVars

   constraints = []
   threshold = []

   oldConstraints = [(greaterThans, gtThreshold, -1), (lessThans, ltThreshold, 1),
                     (equalities, eqThreshold, 0)]

   offset = 0
   for constraintList, oldThreshold, coefficient in oldConstraints:
      constraints += [c + r for c, r in zip(constraintList,
         identity(numRows, newVars, coefficient, offset))]

      threshold += oldThreshold
      offset += len(oldThreshold)

   return newCost, constraints, threshold

def num_entries(a):
    """Number of entries in a sparse matrix, represented as above.

    This is for tracking memory usage.
    """
    return sum([len(x) for x in a.values()])

def dot(a,b):
    columns = set(a.keys()).union(set(b.keys()))
    return sum(a[j] * b[j] for j in columns)

def column_orig(A, j):
    """Returns j'th column of A, as a list (with 0's for missing entries).
    FIXME this doesn't work so well for sparse matrices"""
    return [(row[j] if j in row else 0)
        for j in A]

def column(A, j):
    """Returns j'th column of A, as a list of entries.

    ??? should this return a dict?"""
    column = [(i, row[j])
        for (i, row) in A.items()
        if j in row and row[j]!=0]
    return column

# deprecated; try to omit this?
# def transpose(A):
#    return [column(A, j) for j in range(len(A[0]))]

def isPivotCol_orig(col):
    return (len([c for c in col if c == 0]) == len(col) - 1) and sum(col) == 1
def isPivotCol(col):
    """Pivot columns have only a 1 (in some row), and aren't column -1."""
    return len(col)==1 and col[0][1]==1 and col[0][0]!=-1

def variableValueForPivotColumn_1(tableau, column):
    pivotRow1 = [i for (i, x) in enumerate(column) if x == 1]
    assert(len(pivotRow1)==1)
    pivotRow = pivotRow1[0]
    return tableau[pivotRow][-1]

def variableValueForPivotColumn(tableau, j):
    """If we know j is a pivot column, gets
        the variable in the row which contains a 1.
    """
    col = column(tableau, j)
    if len(col)==1 and col[0][1]==1:
        i = col[0][0]
        return tableau[i][-1]
    else:
        return None

def initialTableau(c, A, b):
    """Constructs the initial tableau.

    c: the objective function, as a dict (indexed by column) of fractions
    A: the tableau of the matrix, as a dict of dict of fractions
    b: the constraints, as a dict (indexed by row) of fractions
    assume the last m columns of A are the slack variables; the initial basis is
    the set of slack variables
    """
    tableau = A
    tableau[-1] = c
    tableau[-1][-1] = 0
    for (i, x) in b.items():
        tableau[i][-1] = x
    return tableau

def primalSolution(tableau):
    # the pivot columns denote which variables are used
    # columns = transpose(tableau)
    # indices = [j for j, col in enumerate(columns[:-1]) if isPivotCol(col)]
    columns_list = [set(r.keys()) for r in tableau.values()]
    all_columns = functools.reduce(lambda a, b: a.union(b), columns_list)
    indices = [j for j in all_columns
        if isPivotCol(column(tableau, j))]
    return {colIndex: variableValueForPivotColumn(tableau, colIndex)
        for colIndex in indices}

def objectiveValue(tableau):
    return -(tableau[-1][-1])

def canImprove(tableau):
    lastRow = [a for (j,a) in tableau[-1].items()
        if j != -1]
    return any([x>0 for x in lastRow])

# this can be slightly faster
def moreThanOneMin(L):
    if len(L) <= 1:
        return False
    x,y = heapq.nsmallest(2, L, key=lambda x: x[1])
    return x == y

def sparsifyRows(tableau):
    """Removes 0's from rows (except for row or column -1)."""
    def sparsify1(row):
        return {j:x for (j,x) in row.items() if j==-1 or x!=0 }
    return {i: (sparsify1(row) if i!=-1 else row)
        for (i,row) in tableau.items()}

def findPivotIndex(tableau):
    # pick minimum positive index of the last row
    column_choices = [(j,x) for (j,x) in tableau[-1].items()
        if x>0 and j>=0]
    j = min(column_choices, key=lambda a: a[1])[0]

    # check if unbounded
    if all([x <= 0 for (i,x) in column(tableau, j)]):
        raise Exception('Linear program is unbounded.')

    # compute quotients
    quotients = [(i, r[-1] / r[j])
        for (i,r) in tableau.items()
        if i>=0 and j in r and r[j]>0]
    # check for degeneracy: more than one minimizer of the quotient
    if moreThanOneMin(quotients):
        raise Exception('Linear program is degenerate.')

    # pick row index minimizing the quotient
    i = min(quotients, key=lambda x: x[1])[0]

    return i, j


def pivotAbout(tableau, pivot):
    i,j = pivot

    # rescale row containing pivot, so that tableau[i][j] := 1
    pivotDenom = tableau[i][j]
    tableau[i] = {j: x / pivotDenom
        for (j,x) in tableau[i].items()}

    # subtract off the rescaled row, so that for the other rows,
    # the j'th column is 0
    # ??? can we just loop over rows in which column j is nonempty?
    for i1,row in tableau.items():
        if i1 != i and j in row:
            # was:
            # pivotRowMultiple = [y * tableau[k][j] for y in tableau[i]]
            # tableau[k] = [x - y for x,y in zip(tableau[k], pivotRowMultiple)]
            scale = tableau[i1][j]
            # note that we need to consider every column in which _either_
            # row i, or row i1, is non-zero
            for j1 in set(row).union(set(tableau[i])):
                x = tableau[i1][j1] if j1 in tableau[i1] else 0
                if j1 in tableau[i]:
                    x -= scale * tableau[i][j1]
                tableau[i1][j1] = x

'''
    Checks that a tableau is valid.
'''
def isValidTableau(tableau):
    # check that b >= 0
    for i,row in tableau.items():
        if i<0:
            continue
        if row[-1] < 0:
            return False
    # check number of basis columns
    # FIXME
    return True



'''
   The simplex algorithm, as a function from "tableau => tableau".

   This seems possibly convenient for the two-phase simplex algorithm.
'''
def tableauSimplex(tableau, verbosity=0):
   tableau = sparsifyRows(tableau)

   # make sure every non-objective row has b >= 0
   if True:
     for i, row in tableau.items():
        if i < 0:
            continue
        if row[-1] < 0:
            # this doesn't actually update the tableau!
            # row = {j: -x for j, x in row.items()} 
            # this works; but then the solver sometimes crashes...
            tableau[i] = {j: -x for j, x in row.items()} 
   # pdb.set_trace()

   if verbosity >= 1:
       print("Initial tableau:")
       for row in tableau.items():
          print(row)
       print()
   # if not isValidTableau(tableau):
   iter = 0
   print("iter\tobjective\tfloat(objective)\tn. tableau entries")
   print("\t".join([str(iter),
        str(-tableau[-1][-1]),
        str(float(-tableau[-1][-1])),
        str(num_entries(tableau))]))

   while canImprove(tableau):
      # copy tableau beforehand (for debugging)
      tableau_before = copy.deepcopy(tableau)
      if verbosity >= 3:
         print(f"tableau =")
         for x in tableau.items():
            print(x)
      pivot = findPivotIndex(tableau)
      if verbosity >= 2:
         print("Next pivot index is=%d,%d \n" % pivot)
      pivotAbout(tableau, pivot)
      tableau = sparsifyRows(tableau)

      # check tableau
      if not isValidTableau(tableau):
        pdb.set_trace()

      if verbosity >= 2:
          print("Tableau after pivot:")
          for row in tableau.items():
             print(row)
          print()

      # tableau = sparsifyRows(tableau)

      iter += 1
      print("\t".join([str(iter),
         str(objectiveValue(tableau)),
         str(float(objectiveValue(tableau))),
         str(num_entries(tableau))]))

   return tableau


'''
   simplex: [float], [[float]], [float] -> [float], float
   Solve the given standard-form linear program:

      max <c,x>
      s.t. Ax = b
           x >= 0

   providing the optimal solution x* and the value of the objective function

   Note that this only does one phase of the simplex algorithm.
'''
def simplex(c, A, b, verbosity=0):
   tableau = initialTableau(c, A, b)
   tableau = tableauSimplex(tableau, verbosity=verbosity)

   return tableau, primalSolution(tableau), objectiveValue(tableau)


'''
    Solve the given linear program, using the two-phase simplex algorithm.

    This is like simplex(), but uses the two-phase simplex algorithm,
    and so should be able to find an initial feasible solution.

    Attempting to follow description at
        https://sites.math.washington.edu/~burke/crs/407/notes/section3-18.pdf
'''
def simplex_two_phase(c, A, b, verbosity=0):
    # start with the initial simplex
    tableau = initialTableau(c, A, b)

    # first, renumber current objective from -1 to -2
    tableau[-2] = tableau.pop(-1)

    # first phase: add variable to optimize 
    # FIXME base this on the number of variables?
    phase_1_objective = 1000000
    tableau[-1] = { phase_1_objective: -1, -1: 0 }
    # artificial_vars.add(artificial_var_index)
    # artificial_var_index += 1
    for (i, row) in tableau.items():
        # -1 and -2 are objective functions; don't add artificial vars. for them
        if i < 0:
            continue
        row[phase_1_objective] = -1

    pdb.set_trace()
    tableau = tableauSimplex(tableau, verbosity=verbosity)
    pdb.set_trace()

    # check for feasibility
    if objectiveValue(tableau) > 0:
       raise ValueError("problem is infeasible")

    # second phase
    # first, remove phase 1 objective (which should no longer be needed)
    tableau.pop(-1)
    for (i, row) in tableau.items():
        try:
            row.pop(phase_1_objective)
        except KeyError:
            pass
    # restore the original objective
    tableau[-1] = tableau.pop(-2)

    tableau = tableauSimplex(tableau, verbosity=verbosity)

    return tableau, primalSolution(tableau), objectiveValue(tableau)

'''
    Solve the given linear program, using the two-phase simplex algorithm.

    This is like simplex(), but uses the two-phase simplex algorithm,
    and so should be able to find an initial feasible solution.

    Attempting to follow description at
        https://en.wikipedia.org/wiki/Simplex_algorithm
'''
def simplex_two_phase_v1(c, A, b, verbosity=0):
    # start with the initial simplex
    tableau = initialTableau(c, A, b)

    # first phase: add artificial variables
    # (first, renumber current objective from -1 to -2)
    # for (i, row) in tableau.items():
    #     row[-2] = row.pop(-1)
    tableau[-2] = tableau.pop(-1)

    artificial_vars = set()
    # FIXME base this on the number of variables?
    artificial_var_index = 1000000

    # this is the objective function for finding an initial feasible sol'n
    tableau[-1] = { -1: 0 }
    # artificial_vars.add(artificial_var_index)
    # artificial_var_index += 1
    for (i, row) in tableau.items():
        # -1 and -2 are objective functions; don't add artificial vars. for them
        if i < 0:
            continue
        tableau[-1][artificial_var_index] = -1
        row[artificial_var_index] = 1
        artificial_vars.add(artificial_var_index)
        artificial_var_index += 1

    # based on Wiki, it appears that the objective function needs to be tweaked first
    # (possibly this should be fused with the loop above?)
    for (i, row) in tableau.items():
        # again, skip objective functions
        if i < 0:
            continue
        # only consider rows which contain an artificial var
        if not artificial_vars.intersection(row):
            continue
        # at this point, the row must contain an artificial var (which presumably =1);
        # add it to row -1 (the first phase objective function) 
        for j in set(row).union(set(tableau[-1])):
            x = tableau[-1][j] if j in tableau[-1] else 0
            if j in row:
                x += row[j]
            tableau[-1][j] = x

    # pdb.set_trace()
    print(f"before phase1: objective = {objectiveValue(tableau)}")
    tableau = tableauSimplex(tableau, verbosity=verbosity)
    print(f"after phase1: objective = {objectiveValue(tableau)}")
    # pdb.set_trace()

    # FIXME check for feasibility

    # second phase
    # first, remove artificial variables (hopefully they're not needed)
    tableau.pop(-1)
    for (i, row) in tableau.items():
        for j in artificial_vars:
            if j in row:
                row.pop(j)
    # for (i, row) in tableau.items():
    #     row[-1] = row.pop(-2)
    tableau[-1] = tableau.pop(-2)

    tableau = tableauSimplex(tableau, verbosity=verbosity)

    return tableau, primalSolution(tableau), objectiveValue(tableau)


