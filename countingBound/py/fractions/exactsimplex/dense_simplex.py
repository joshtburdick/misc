# Experimental dense implementation of the simplex method.
# This is intended as a reference implementation.


import numpy as np




class DenseSimplex:
    """Dense simplex method.

    Rather than adding rows together, this maintains a square matrix
    which keeps track of which rows would've been multiplied together.
    """



    def __init__(self, A, b, c):
        """Constructor.

        A, b: the "Ax <= b" constraints
        c: the objective function (as something to be maximized)
        """
        m, n = A.shape

        # the tableau
        T1 = np.concat([A, np.eye(m), b.reshape((-1,1))], axis=1)
        T2 = np.concat([c, np.zeros(m+1)], axis=1)
        self.T = np.concat([T1, T2], axis=0)

        # matrix to do the row-math book-keeping
        self.M = np.eye(m.shape[0])

        # basis[j] is the column index in T corresponding to M[:,j]
        self.basis = range(n+1, n+m+1)

        # make sure all rows have b >= 0




