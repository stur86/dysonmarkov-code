"""
Helper code written for 

Simone Sturniolo
"A Dyson equation approach for averaging of classical and quantum observables on multiple realizations of Markov processes"

https://arxiv.org/abs/2004.01183

MIT License

Copyright (c) 2020 Simone Sturniolo

(full text in LICENSE file)
"""

import time
import numpy as np
from dysonavg import *

lie_group = np.array([
    [[0, 0,  0],
     [0, 0, -1],
     [0, 1,  0]],
    [[ 0, 0, 1],
     [ 0, 0, 0],
     [-1, 0, 0]],
    [[0, -1, 0],
     [1,  0, 0],
     [0,  0, 0]]
])

# Frequencies
axes = np.array([[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]])
w = 1.0j*np.sum(axes[:,:,None,None]*lie_group[None,:,:,:], axis=1)

N = len(w)
dt = 1e-3
steps = 10000

nurange = np.array([1e-2, 1e-1, 5e-1, 1, 10, 100])
times = []

# Force compilation of Numba functions
numbacompile()

for i, nu in enumerate(nurange):
    print('Testing for {0}'.format(nu))
    Q = expQ(N=N, nu=nu)
    Q[:,0] *= 0
    t0 = time.time()
    t, dy = dysonMarkovFID(w, Q, dt, steps=steps)
    t1 = time.time()
    c, mc = monteCarloMarkovFID(w, Q, dt, steps=steps, samples=500, 
                                use_gillespie=True, batch_size=20)
    t2 = time.time()

    mcv = np.tensordot(mc, [1,0,0], axes=(2,0))
    dyv = np.tensordot(dy, [1,0,0], axes=(2,0))

    data = np.real(np.concatenate([t[:,None],c[:,None], mcv, dyv], axis=1))

    times.append([t1-t0, t2-t1])

    np.savetxt('results/vector_{0}.dat'.format(i), data)

# Save times
np.savetxt('results/vector_times.dat', np.array(times)*1e3)

# # Extremes
# lf = np.average(np.cos(t[:,None]*w[None,:]), axis=1)
# np.savetxt('results/vector_static.dat', np.real([t, lf]).T)
# hf = np.cos(t*np.average(w))
# np.savetxt('results/scalar_fast.dat', np.real([t, hf]).T)

