"""
Helper code written for 

Simone Sturniolo
"A Dyson equation approach for averaging of classical and quantum observables on multiple realizations of Markov processes"

https://arxiv.org/abs/2004.01183

MIT License

Copyright (c) 2020 Simone Sturniolo

(full text in LICENSE file)
"""

# Functions to be used for Dyson equation calculation of FID, as well as the Monte Carlo test benchmark
import numpy as np
from numba import jit, objmode
from scipy.integrate import cumtrapz
from scipy.linalg import expm
from scipy.linalg.lapack import ctrtri

# Exponential decorrelation Markov process matrix


def expM(N=None, weights=None, nu_t=0.1):

    if weights is None:
        if N is None:
            raise ValueError('One of N or weights must be specified')
        weights = np.ones(N)/N
    elif N is None:
        weights = np.array(weights)
        weights /= np.sum(weights)
        N = len(weights)

    d = np.exp(-nu_t)
    M = weights[:, None]*np.ones((N, N))*(1-d)
    M *= (1-np.eye(N))
    M += np.diag(1-np.sum(M, axis=0))

    return M

# Exponential decorrelation, Q-matrix


def expQ(N=None, weights=None, nu=0):

    if weights is None:
        if N is None:
            raise ValueError('One of N or weights must be specified')
        weights = np.ones(N)/N
    elif N is None:
        weights = np.array(weights)
        weights /= np.sum(weights)
        N = len(weights)

    Q = np.ones((N, N))*nu*weights[:, None]
    Q -= np.eye(N)*nu

    return Q

# Numba version of matrix exponential


@jit(nopython=True, cache=True)
def _nbexpm(M):
    e, ev = np.linalg.eig(M+0.j)
    return np.dot(np.dot(ev, np.diag(np.exp(e))), np.linalg.inv(ev))

@jit(nopython=True, cache=True)
def _objexpm(M):
    with objmode(answer='complex128[:,:]'):
        answer = expm(M)
    return answer

# Gillespie algorithm for simulation


@jit(nopython=True, cache=True)
def _gillMC(weights, Q, tmax, samples, tblock=100):

    N = len(weights)
    Qjump = Q*(1-np.eye(N))
    for i in range(N):
        Qjump[:, i] = np.cumsum(Qjump[:, i])
    Qsum = Qjump[-1, :].copy()
    for i in range(N):
        if Qsum[i] > 0:
            Qjump[:, i] /= Qsum[i]

    cumw = np.cumsum(weights)
    cumw /= cumw[-1]

    trajs = np.zeros((samples, tblock, 2))

    for s_i in range(samples):
        t = 0
        t_i = 0
        c = np.sum(np.random.random() > cumw)
        trajs[s_i, 0, 0] = 0
        trajs[s_i, 0, 1] = c

        while t < tmax:
            if Qsum[c] == 0:
                # This means literally no transition is possible
                dt = tmax
                c = c
            else:
                dt = -np.log(np.random.random())/Qsum[c]
                c = np.sum(np.random.random() > Qjump[:, c])

            t += dt
            t_i += 1

            trajs[s_i, t_i, 0] = t
            trajs[s_i, t_i, 1] = c

            if t_i >= trajs.shape[1]:
                trajs = np.concatenate((trajs,
                                        np.zeros((samples, tblock, 2))),
                                       axis=1)

    return trajs


@jit(nopython=True, cache=True)
def _gillMCfid(O_list, weights, Q, dt, steps, samples,
                   tblock=100):
    
    phi_trajs = np.zeros((samples, steps))
    conf_trajs = np.zeros((samples, steps))


    # First, compute trajectories
    trajs = _gillMC(weights, Q, dt*(steps-1), samples, tblock)

    # out_array has to be properly shaped already
    for i, tr in enumerate(trajs):
        # Starting configuration?
        c = int(tr[0, 1])
        # Times?
        t0 = 0
        event = 1
        conf_trajs[i, 0] = c
        for j in np.arange(1, steps):
            t = j*dt
            phi_trajs[i, j] = phi_trajs[i, j-1]
            # Is the next event before or after this?
            while tr[event, 0] < t:
                phi_trajs[i, j] += O_list[c]*(tr[event, 0]-t0)
                t0 = tr[event, 0]
                c = int(tr[event, 1])
                event += 1

            conf_trajs[i, j] = c
            phi_trajs[i, j] += O_list[c]*(t-t0)
            t0 = t

    fid = np.sum(np.exp(1.0j*phi_trajs), axis=0)/samples
    corr = np.zeros(steps)

    for i in range(steps):
        corr[i] = np.sum(conf_trajs[:,i] == conf_trajs[:,0])/samples

    return corr, fid


@jit(nopython=True, cache=True)
def _gillMCfidmat(O_list, weights, Q, dt, steps, samples,
                      tblock=100):
    
    N = len(O_list)
    D = O_list.shape[1]

    mat_trajs = np.zeros((samples, steps, D, D)) + 0.0j
    conf_trajs = np.zeros((samples, steps))

    # Diagonalise the matrices in O_list
    O_evals = np.zeros((N, D)) + 0.0j
    O_evecsL = np.zeros((N, D, D)) + 0.0j
    O_evecsR = np.zeros((N, D, D)) + 0.0j

    for i in range(N):
        evals, evecs = np.linalg.eig(O_list[i])
        O_evals[i] = evals
        O_evecsL[i] = evecs
        O_evecsR[i] = np.linalg.inv(evecs)


    # First, compute trajectories
    trajs = _gillMC(weights, Q, dt*(steps-1), samples, tblock)

    # out_array has to be properly shaped already
    for i, tr in enumerate(trajs):
        # Starting configuration?
        c = int(tr[0, 1])
        # Times?
        t0 = 0
        event = 1
        conf_trajs[i, 0] = c
        mat_trajs[i, 0] = np.eye(D)
        for j in np.arange(1, steps):
            t = j*dt
            mat_trajs[i, j] = mat_trajs[i, j-1].copy()
            # Is the next event before or after this?
            while tr[event, 0] < t:
                # Find the matrix exponential
                mdt = np.diag(np.exp(1.0j*O_evals[c]*(tr[event, 0]-t0)))
                mdt = np.dot(np.dot(O_evecsL[c], mdt), O_evecsR[c])

                mat_trajs[i, j] = np.dot(mdt, mat_trajs[i, j])
                t0 = tr[event, 0]
                c = int(tr[event, 1])
                event += 1

            conf_trajs[i, j] = c

            mdt = np.diag(np.exp(1.0j*O_evals[c]*(t-t0)))
            mdt = np.dot(np.dot(O_evecsL[c], mdt), O_evecsR[c])

            mat_trajs[i, j] = np.dot(mdt, mat_trajs[i, j])
            t0 = t

    fid = np.sum(mat_trajs, axis=0)/samples
    corr = np.zeros(steps)

    for i in range(steps):
        corr[i] = np.sum(conf_trajs[:,i] == conf_trajs[:,0])/samples

    return corr, fid

def numbacompile():
    # Force one usage of both Numba functions to compile them
    _gillMCfid(np.zeros(2), np.zeros(2), np.eye(2), 1e-2, 1, 1)
    _gillMCfidmat(np.zeros((2,3,3)), np.zeros(2), np.eye(2), 1e-2, 1, 1)

# Monte Carlo simulation of FID over a Markov process


def monteCarloMarkovFID(O_list, Q, dt, steps=100, samples=1000, weights=None,
                        batch_size=1000, use_gillespie=False):

    M = expm(Q*dt)
    O_list = np.array(O_list)
    N = O_list.shape[0]

    if weights is None:
        weights = np.ones(N)/N
    else:
        weights = np.array(weights)
        weights /= np.sum(weights)
        if weights.shape[0] != N:
            raise ValueError('Invalid weights array')

    csumw = np.cumsum(weights)
    csumM = np.cumsum(M, axis=0)

    isMat = (len(O_list.shape) == 3)

    if isMat:
        D = O_list.shape[1]
        fid = np.zeros((steps, D, D))+.0j
    else:
        if len(O_list.shape) != 1:
            raise RuntimeError('Invalid O_list')
        fid = np.zeros(steps)+.0j

    if not use_gillespie:

        for s0 in np.arange(0, samples, batch_size):
            sbatch = range(s0, min(s0+batch_size, samples))
            sn = len(sbatch)
            configs = np.zeros((sn, steps)).astype(int)
            configs[:, 0] = np.sum(np.random.random(
                sn)[:, None] > np.cumsum(weights)[None, :], axis=1)
            # Jumps?
            for i in range(1, steps):
                jumps = np.random.random(sn)
                configs[:, i] = np.sum(jumps[None, :] > np.cumsum(M, axis=0)[:,
                                                                             configs[:, i-1]],
                                       axis=0)

            # Self correlation function?
            corr = np.average(configs[:, 0, None] == configs, axis=0)
            if not isMat:
                phi = 1.0j*cumtrapz(O_list[configs], dx=dt, axis=1, initial=0.0)
                fidpart = np.exp(phi)
            else:
                eO_list = np.array([expm(1.0j*o*dt) for o in O_list])
                fidpart = np.zeros((sn, steps, D, D)).astype(np.complex)
                fidpart[:, 0] = np.eye(D)[None, :, :]
                for i in range(1, steps):
                    fidpart[:, i] = np.sum(eO_list[configs[:, i]][:, :, :, None]
                                           * fidpart[:, i-1][:, None, :, :], axis=2)
            fid += np.sum(fidpart, axis=0)
        fid /= samples

    else:
        corr = np.zeros(steps)
        for s0 in np.arange(0, samples, batch_size):
            sbatch = range(s0, min(s0+batch_size, samples))
            sn = len(sbatch)

            if not isMat:
                c, f = _gillMCfid(O_list, weights, Q, dt, steps, batch_size)
            else:
                c, f = _gillMCfidmat(O_list, weights, Q, dt, steps, batch_size)

            corr += np.real(c)*sn
            fid += f*sn

        corr /= samples
        fid /= samples

    return corr, fid


def dysonMarkovFID(O_list, Q, dt, steps=100, weights=None, use_trapz=True):

    M = expm(Q*dt)
    O_list = np.array(O_list)
    N = len(O_list)
    D = 1
    if len(O_list.shape) == 3:
        D = O_list.shape[1]

    if weights is None:
        weights = np.ones(N)/N
    else:
        weights = np.array(weights)
        weights /= np.sum(weights)
        if weights.shape[0] != N:
            raise ValueError('Invalid weights array')

    # Single step G0
    X = np.kron(M, np.eye(D))

    # Step matrix for iterative generation
    stepD = np.zeros((N*D, N*D)).astype(np.complex)
    stepF = np.zeros((N*D, N*D)).astype(np.complex)

    if use_trapz:
        for i, o in enumerate(O_list):
            stepD[i*D:(i+1)*D, i*D:(i+1)*D] = (np.eye(D)+0.5j*dt*o)
            stepF[i*D:(i+1)*D, i*D:(i+1) *
                  D] = np.linalg.inv(np.eye(D)-0.5j*dt*o)
        stepM = np.dot(stepF, np.dot(X, stepD))
    else:
        for i, o in enumerate(O_list):
            stepF[i*D:(i+1)*D, i*D:(i+1) *
                  D] = np.linalg.inv(np.eye(D)-1.0j*dt*o)
        stepM = np.dot(stepF, X)

    # Array of Gs
    Ginitial = np.kron(np.diag(weights), np.eye(D))
    G = Ginitial.copy()
    Glist = [np.eye(D)]
    for ti in range(1, steps):
        G = np.dot(stepM, G)
        Glist.append(np.sum(G.reshape((N, D, N, D)), axis=(0, 2)))
    Glist = np.array(Glist)

    if D == 1:
        Glist = Glist[:, 0, 0]

    return np.arange(0, steps)*dt, Glist
