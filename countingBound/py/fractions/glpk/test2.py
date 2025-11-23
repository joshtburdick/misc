#!/usr/bin/env python3

import numpy as np
import pulp

from fractions import Fraction

prob = pulp.LpProblem("sidsproblem", pulp.LpMinimize)

x1 = pulp.LpVariable("x_1",5)
x2 = pulp.LpVariable("x_2",0,1)

prob += x1 + x2
prob += ((2*x1 + x2) >= 1)
prob += 0.5 * x1 + 0.3 * x2 == 1
# prob += (x1 + 2*x2 >= 1)
prob += (pulp.lpSum([1 * x1, Fraction(2, 1) * x2]) >= 1)

prob.writeLP("sidsproblem2.lp")

x = prob.solve(pulp.GLPK(options=['--exact']))
print(x)
for v in prob.variables():
    print(v.name, "=", v.varValue)

