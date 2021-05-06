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
from scipy.special import binom, factorial

# Frequencies
w = np.array([-0.5, 1, 2])
N = len(w)
dt0 = 1e-3
steps = 10000

dt_test = np.array([1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3])

nurange = np.array([1e-2, 1e-1, 1, 5, 10, 100])
times = []

numbacompile()

for i, nu in enumerate(nurange):
    print('Testing for {0}'.format(nu))

    Q = expQ(N=N, nu=nu)
    # Start with the reference
    samples_each = 1000
    samples_runs = 50

    mc_runs = []

    for j in range(samples_runs):
        _, y0 = monteCarloMarkovFID(w, Q, dt0, steps=steps,
                                    samples=samples_each, use_gillespie=True, 
                                    batch_size=20)
        mc_runs.append(y0)
    # Now get an overall reference as well as std
    mc_avg = np.average(mc_runs, axis=0)
    mc_std = np.std(mc_runs, axis=0)

    # Now the Dyson tests
    perf = []
    for dt1 in dt_test:
        st1 = int(steps*dt0/dt1)
        skip = int(steps/st1)
        t0 = time.time()
        t, y1 = dysonMarkovFID(w, Q, dt1, steps=st1)
        t1 = time.time()

        err = np.average(np.abs(np.real(y1-mc_avg[::skip])/mc_std[::skip])[1:])

        perf.append([t1-t0, err])

    np.savetxt('results/benchmark_{0}.dat'.format(i), perf)
