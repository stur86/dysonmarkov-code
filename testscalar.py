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

# Frequencies
w = np.array([-0.5, 1, 2])
N = len(w)
dt = 1e-3
steps = 10000

nurange = np.array([1e-2, 1e-1, 1, 5, 10, 100])
times = []

# Force compilation of Numba functions
numbacompile()

for i, nu in enumerate(nurange):
    print('Testing for {0}'.format(nu))
    Q = expQ(N=N, nu=nu)
    t0 = time.time()
    t, dy = dysonMarkovFID(w, Q, dt, steps=steps)
    t1 = time.time()
    c, mc = monteCarloMarkovFID(w, Q, dt, steps=steps, samples=500, use_gillespie=True, batch_size=10)
    t2 = time.time()

    times.append([t1-t0, t2-t1])

    np.savetxt('results/scalar_{0}.dat'.format(i), np.real([t, c, mc, dy]).T)

# Save times
np.savetxt('results/scalar_times.dat', np.array(times)*1e3)

# Extremes
lf = np.average(np.cos(t[:,None]*w[None,:]), axis=1)
np.savetxt('results/scalar_static.dat', np.real([t, lf]).T)
hf = np.cos(t*np.average(w))
np.savetxt('results/scalar_fast.dat', np.real([t, hf]).T)

