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

Sx = np.array([[0.0, 1.0], [1.0, 0.0j]])
Sy = np.array([[0.0, -1.0j], [1.0j, 0.0]])
Sz = np.array([[1.0, 0.0], [0.0j, -1.0]])

A = 1

II = np.kron(Sx, Sx)+np.kron(Sy, Sy)+np.kron(Sz, Sz)

H_list = [-II+A*np.kron(Sz, np.eye(2))+A*np.kron(np.eye(2), Sz),
          -II-A*np.kron(Sz, np.eye(2))-A*np.kron(np.eye(2), Sz),
          -II+A*np.kron(Sz, np.eye(2))-A*np.kron(np.eye(2), Sz),
          -II-A*np.kron(Sz, np.eye(2))+A*np.kron(np.eye(2), Sz)]
qsize = H_list[0].shape[0]
L_list = np.array([np.kron(H, np.eye(qsize)) -
                   np.kron(np.eye(qsize), H.T) for H in H_list])
rho0 = 0.25*np.ones((4, 4)).reshape((-1,))

N = len(H_list)
dt = 1e-3
steps = 10000

nurange = np.array([1e-2, 1e-1, 5e-1, 1, 10, 100])
times = []

# Force compilation of Numba functions
numbacompile()

for i, nu in enumerate(nurange):
    print('Testing for {0}'.format(nu))
    Q = expQ(N=N, nu=nu)
    t0 = time.time()
    t, dy = dysonMarkovFID(L_list, Q, dt, steps=steps)
    t1 = time.time()
    c, mc = monteCarloMarkovFID(L_list, Q, dt, steps=steps, samples=500, 
                                use_gillespie=True, batch_size=20)
    t2 = time.time()

    mcrho = np.tensordot(mc, rho0, axes=(2, 0))
    dyrho = np.tensordot(dy, rho0, axes=(2, 0))

    times.append([t1-t0, t2-t1])

    # Now expectation values
    dyII = np.sum(II.T.reshape((-1,))[None, :]*dyrho, axis=1)
    mcII = np.sum(II.T.reshape((-1,))[None, :]*mcrho, axis=1)

    # And quantum entropy
    dyS = []
    mcS = []
    for dyr, mcr in zip(dyrho, mcrho):
        evals, evecs = np.linalg.eigh(dyr.reshape((4,4)))
        logp = np.where(evals > 0, np.log(evals), 0)
        dyS.append(np.sum(-evals*logp))
        evals, evecs = np.linalg.eigh(mcr.reshape((4,4)))
        logp = np.where(evals > 0, np.log(evals), 0)
        mcS.append(np.sum(-evals*logp))

    np.savetxt('results/quantum_{0}.dat'.format(i),
               np.real(np.array([t, c, dyII, mcII, dyS, mcS]).T))

# Save times
np.savetxt('results/quantum_times.dat', np.array(times)*1e3)

# # Extremes
# lf = np.average(np.cos(t[:,None]*w[None,:]), axis=1)
# np.savetxt('results/vector_static.dat', np.real([t, lf]).T)
# hf = np.cos(t*np.average(w))
# np.savetxt('results/scalar_fast.dat', np.real([t, hf]).T)
