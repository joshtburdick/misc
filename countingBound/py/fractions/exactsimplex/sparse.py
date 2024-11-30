'''Sparse version of simplex LP solver.

A row will be represented by a Dict<Int, Fraction>, as it were.
The tableau will be represented by a Dict<Int, Dict<Int, Fraction>>.
'''

import functools
import heapq
import pdb

'''
   Return a rectangular identity matrix with the specified diagonal entiries, possibly
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
        if j in row]
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
        if x>0 and j!=-1]
    j = min(column_choices, key=lambda a: a[1])[0]

    # check if unbounded
    if all([x <= 0 for (i,x) in column(tableau, j)]):
        raise Exception('Linear program is unbounded.')

    # check for degeneracy: more than one minimizer of the quotient
    quotients = [(i, r[-1] / r[j])
        for (i,r) in tableau.items()
        if i!=-1 and -1 in r and j in r and r[j]>0]
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
            for j1 in row:
                if j1 in tableau[i]:
                    tableau[i1][j1] -= scale * tableau[i][j1]


'''
   simplex: [float], [[float]], [float] -> [float], float
   Solve the given standard-form linear program:

      max <c,x>
      s.t. Ax = b
           x >= 0

   providing the optimal solution x* and the value of the objective function
'''
def simplex(c, A, b):
   tableau = initialTableau(c, A, b)
   print("Initial tableau:")
   for row in tableau.items():
      print(row)
   print()

   while canImprove(tableau):
      pivot = findPivotIndex(sparsifyRows(tableau))
      print("Next pivot index is=%d,%d \n" % pivot)
      pivotAbout(tableau, pivot)
      print("Tableau after pivot:")
      for row in tableau.items():
         print(row)
      print()

   tableau = sparsifyRows(tableau)

   return tableau, primalSolution(tableau), objectiveValue(tableau)
