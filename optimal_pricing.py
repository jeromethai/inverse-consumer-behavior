'''
Created on Sep 5, 2014

@author: jeromethai
'''

from cvxopt import matrix
from utility_estimator import x_solver

"""
def x_solver(Q, r, p, H, z, soft):
    Optimization w.r.t. x-block
    
    Parameters
    ----------
    Q, r: parameters of the quadratic utility function
    p: price vector
    H: observation matrix
    z: observed demand vector
    soft: weight on the observation

def mis_solver(zs, H, ps, smooth, soft=1e6, max_iter=3):
    Solves the inverse optimization problem with missing values
    
    Parameters
    ----------
    zs: list of observed demand vectors
    H: observation matrix
    ps: list of price vectors
    smooth: regularization parameter on r
    soft: weight on the observation
    max_iter: number of iterations
    
    n = H.size[1]
    Q, r = -spdiag([10.]*n), matrix(10., (n,1))
    for k in range(max_iter):
        xs = [x_solver(Q, r, p, H, z, soft) for p,z in zip(ps,zs)]
        Q, r = inv_solver(xs, ps, compute_xmax(xs), smooth)
    return Q, r
"""

def solver(Q, r, i, d, p, soft=1e8, max_iter=3):
    """Pricing for product i that matches desired demand d
    given other prices p and quadratic utility function (Q,r)
    
    Parameters
    ----------
    (Q,r): quadratic utility function
    i: computes optimal price for product i, 0<=i<=n-1
    d: target demand for product i
    p: full price vector (price of product i, p^i, is arbitrary) 
    max_iter: number of iterations
    """
    n = len(r)
    H = matrix(0., (1,n)); H[i] = 1.
    q = p
    q[i] = 0.
    print q[i]
    for k in range(max_iter):
        x = x_solver(Q, r, q, H, d, soft)
        q[i] = max((H*(Q*x+r))[0], 0)
        print q[i]
    return q


if __name__ == '__main__':
    pass