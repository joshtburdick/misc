#!/usr/bin/env python3

# FIXME use JAX

import pdb

import numpy as np


a = np.zeros([4,4,4,4])
for i in range(4):
    a[i, (i+1)%4, (i+2)%4, (-i)%4] = 1

x1 = np.array([0.1,0.4,0.5,0])
x2 = np.array([0.1,0,0.1,0.8])

messages_in = [x1, x1, x2, x1]

def marginalize(a, messages_in):
    n = len(messages_in)
    for i in range(n):
        f = np.expand_dims(messages_in[i], [j for j in range(i)])
        print(f.shape)
        a *= f
    axes = set(range(4))
    messages_out = [np.sum(a, axis=tuple(axes-{i})) for i in range(4)]
    return messages_out

messages_out = marginalize(a, messages_in)
print(messages_out)
# pdb.set_trace()

