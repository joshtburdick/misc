#!/usr/bin/env python4

import pulp

from fractions import Fraction

prob = pulp.LpProblem("sidsproblem", pulp.LpMaximize)

x1 = pulp.LpVariable("x_1",0,1)
x2 = pulp.LpVariable("x_2",0,1)

prob += x1 + x2
prob += (1 + 5*x1 + 4*x2 >= 1)
prob += (4*x1 - 11*x2 <= -1)
prob += (x1 <= 0.23)
prob += (x2 <= Fraction(2,9))
prob.writeLP("sidsproblem2.lp")

prob.solve(pulp.GLPK(options=['--exact']))

for v in prob.variables():
    print(v.name, "=", v.varValue)

