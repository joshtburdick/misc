#!/usr/bin/env python3
# Two-phase simplex algorithm from ChatGPT v. 3.5 .

from fractions import Fraction
from collections import defaultdict

class SparseSimplex:
    def __init__(self, c, A, b, senses):
        """
        c: list of objective function coefficients (maximize c^T x)
        A: list of constraint coefficient dictionaries (sparse: {col: value})
        b: list of RHS values
        senses: list of constraint types: '<=', '>=', or '='
        """
        self.num_vars = len(c)
        self.A = A
        self.b = [Fraction(x) for x in b]
        self.senses = senses
        self.c = [Fraction(x) for x in c]
        self.tableau = []
        self.basis = []
        self.num_constraints = len(A)
        self.var_index = self.num_vars  # next index for slacks/artificial vars
        self.artificials = []

        self.build_initial_tableau()

    def build_initial_tableau(self):
        for i in range(self.num_constraints):
            row = defaultdict(Fraction, self.A[i])
            sense = self.senses[i]

            if sense == '<=':
                row[self.var_index] = Fraction(1)
                self.basis.append(self.var_index)
                self.var_index += 1

            elif sense == '>=':
                row[self.var_index] = Fraction(-1)
                row[self.var_index + 1] = Fraction(1)
                self.basis.append(self.var_index + 1)
                self.artificials.append(self.var_index + 1)
                self.var_index += 2

            elif sense == '=':
                row[self.var_index] = Fraction(1)
                self.basis.append(self.var_index)
                self.artificials.append(self.var_index)
                self.var_index += 1

            else:
                raise ValueError("Invalid constraint sense")

            row['rhs'] = self.b[i]
            self.tableau.append(row)

        # Objective function (negated for maximization)
        self.objective_row = defaultdict(Fraction)
        for j in range(self.num_vars):
            self.objective_row[j] = -self.c[j]
        self.objective_row['rhs'] = Fraction(0)

    def pivot(self, pivot_row, pivot_col):
        pivot_element = self.tableau[pivot_row][pivot_col]
        new_row = defaultdict(Fraction)
        for k, v in self.tableau[pivot_row].items():
            new_row[k] = v / pivot_element
        self.tableau[pivot_row] = new_row

        for i, row in enumerate(self.tableau):
            if i != pivot_row and pivot_col in row:
                factor = row[pivot_col]
                for k, v in self.tableau[pivot_row].items():
                    row[k] -= factor * v

        # Update objective row
        if pivot_col in self.objective_row:
            factor = self.objective_row[pivot_col]
            for k, v in self.tableau[pivot_row].items():
                self.objective_row[k] -= factor * v

        self.basis[pivot_row] = pivot_col

    def find_pivot(self):
        # Bland's rule for stability (smallest index entering var)
        candidates = [j for j, v in self.objective_row.items() if j != 'rhs' and v < 0]
        if not candidates:
            return None, None  # optimal

        pivot_col = min(candidates)

        # Minimum ratio test
        ratios = []
        for i, row in enumerate(self.tableau):
            if pivot_col in row and row[pivot_col] > 0:
                ratio = row['rhs'] / row[pivot_col]
                ratios.append((ratio, i))

        if not ratios:
            raise Exception("Unbounded solution.")

        _, pivot_row = min(ratios)
        return pivot_row, pivot_col

    def phase_one(self):
        # Build artificial objective
        artificial_obj = defaultdict(Fraction)
        for a in self.artificials:
            artificial_obj[a] = Fraction(-1)
        artificial_obj['rhs'] = Fraction(0)

        # Adjust for initial infeasibility
        for i, row in enumerate(self.tableau):
            for a in self.artificials:
                if a in row:
                    factor = row[a]
                    for k, v in row.items():
                        artificial_obj[k] += factor * v

        # Iteratively pivot for phase one
        while True:
            candidates = [j for j, v in artificial_obj.items() if j != 'rhs' and v < 0]
            if not candidates:
                break  # feasible

            pivot_col = min(candidates)

            ratios = []
            for i, row in enumerate(self.tableau):
                if pivot_col in row and row[pivot_col] > 0:
                    ratio = row['rhs'] / row[pivot_col]
                    ratios.append((ratio, i))

            if not ratios:
                raise Exception("Infeasible problem.")

            _, pivot_row = min(ratios)

            # Pivot in tableau and artificial objective
            pivot_element = self.tableau[pivot_row][pivot_col]
            new_row = defaultdict(Fraction)
            for k, v in self.tableau[pivot_row].items():
                new_row[k] = v / pivot_element
            self.tableau[pivot_row] = new_row

            for i, row in enumerate(self.tableau):
                if i != pivot_row and pivot_col in row:
                    factor = row[pivot_col]
                    for k, v in new_row.items():
                        row[k] -= factor * v

            factor = artificial_obj[pivot_col]
            for k, v in new_row.items():
                artificial_obj[k] -= factor * v

            self.basis[pivot_row] = pivot_col

        if artificial_obj['rhs'] != 0:
            raise Exception("Infeasible problem.")

        # Remove artificial vars from tableau and objective
        for row in self.tableau:
            for a in self.artificials:
                if a in row:
                    del row[a]

        for a in self.artificials:
            if a in self.objective_row:
                del self.objective_row[a]

    def solve(self):
        if self.artificials:
            self.phase_one()

        while True:
            pivot_row, pivot_col = self.find_pivot()
            if pivot_row is None:
                break
            self.pivot(pivot_row, pivot_col)

        solution = [Fraction(0)] * self.num_vars
        for i, var in enumerate(self.basis):
            if var < self.num_vars:
                solution[var] = self.tableau[i]['rhs']

        optimal_value = self.objective_row['rhs']

        return optimal_value, solution


# Example usage
if __name__ == "__main__":
    c = [3, 2]
    A = [
        {0: 2, 1: 1},
        {0: 1, 1: 2},
        {0: 1, 1: -1},
        {0: 1, 1: 0}
    ]
    b = [18, 14, 4, 6]
    senses = ['<=', '<=', '<=', '=']

    simplex = SparseSimplex(c, A, b, senses)
    value, solution = simplex.solve()
    print("Optimal value:", float(value))
    print("Solution:", [float(x) for x in solution])

