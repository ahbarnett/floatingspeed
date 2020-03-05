# Python speed tests for 3D 1/r kernel between sources and targets.
#
# run from cmd line: python lap3dkernel.py

import numpy as np
from numpy import array,zeros,random
import numba
from time import perf_counter as tic

def lap3dcharge_native(y,q,x):
    """evaluate potentials & fields of 3D Laplace charges at non-self targets,
    via naive direct sum. Slow native-python reference implementation.

    Evaluate at the nt (non-self) targets x,

                1   ns-1 
    pot(x_i) = ---  sum  q_j / r_ij                 for i=0,..,nt
               4pi  j=0

    Here R_ij := x_i - y_j is src-targ displacement, where x_i, y_j in R^3.
    r_ij := |R_ij|, where |.| is the 2-norm in R^3.
    No targets coincident with any source are allowed, ie, no self-interaction.

    Inputs:
    y : ns*3 source locations
    q : ns   source charges
    x : nt*3 target locations

    Outputs:
    pot : ns potential values at targets

    Barnett 9/11/18. Cut down 3/3/20 to match language perf tests
    """
    y = np.atleast_2d(y)     # handle ns=1 case: make 1x3 not 3-vecs
    x = np.atleast_2d(x)
    ns = y.shape[0]
    nt = x.shape[0]
    pot = zeros(nt)
    for j in range(ns):       # loop over sources, ie vectorize over targs...
        R = x - y[j]                   # nt*3
        r2 = np.sum(R**2,axis=1)       # squared dists
        r = np.sqrt(r2)                # dists
        pot += (q[j]) / r     # contrib from this src
    return pot

#@numba.njit('void(f8[:,:],f8[:],f8[:,:],f8[:],f8[:,:],b1)',parallel=True,fastmath=True)   # explicit signature, makes it cache? but can't do optional args?
# Also note cache=True fails w/ parallel=True in numba 0.39
@numba.njit(parallel=True,fastmath=True)   # recompiles every run, slow
def lap3dcharge_numba(y,q,x,pot):
    """evaluate pot of 3D Laplace charges, non-self, naive sum,
    numba jit. Writes into pot.
    See lap3dcharge_native.
    pot passed in since njit fails with internal pot=zeros(nt)
    Barnett cut down from lap3dkernels in perilap3d, 3/3/20
    """
    y = np.atleast_2d(y)     # handle ns=1 case: make 1x3 not 3-vecs
    x = np.atleast_2d(x)
    ns = y.shape[0]
    nt = x.shape[0]
    assert(pot.shape==(nt,))
    for i in numba.prange(nt):    # loop over targs
        pot[i] = 0.0
        for j in range(ns):
            R0 = x[i,0]-y[j,0]
            R1 = x[i,1]-y[j,1]
            R2 = x[i,2]-y[j,2]
            r2 = R0**2+R1**2+R2**2
            r = np.sqrt(r2)
            pot[i] += q[j] / r

def perf_test():
    """ test speed of pot in lap3dcharge, eval speeds of slow & jit & self.
    Barnett 9/11/18. cut down 3/3/20
    """
    # perf tests...
    ns = 10000                    # sources
    nt = 10000                    # targs
    y = random.rand(ns,3)     # sources in [0,1]^3
    q = random.randn(ns)      # charges
    x = random.rand(nt,3)     # targs
    #y=np.asfortranarray(y); x=np.asfortranarray(x); q=np.asfortranarray(q)
    u = lap3dcharge_native(y,q,x)    # warm up
    t0=tic()
    u = lap3dcharge_native(y,q,x)    # native python
    t=tic()-t0
    print("native: %d src-targ pairs in %.3g s: %.3g Gpair/s" % (ns*nt,t,ns*nt/t/1e9))
    u2 = zeros(nt)    # numba version writes outputs to arguments
    lap3dcharge_numba(y,q,x,u2)    # warm up
    t0=tic()
    lap3dcharge_numba(y,q,x,u2)     # numba
    t =tic()-t0
    print("numba:  %d src-targ pairs in %.3g s: %.3g Gpair/s" % (ns*nt,t,ns*nt/t/1e9))
    print("pot err numba vs native:  %.3g"%(np.max(np.abs(u-u2))))

if __name__ == "__main__":
    perf_test()

