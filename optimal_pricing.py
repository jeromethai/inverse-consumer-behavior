'''
Created on Sep 5, 2014

@author: jeromethai
'''

from cvxopt import matrix
from utility_estimator import x_solver


def solver(Q, r, i, d, p, soft=1e5, max_iter=5):
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
    q = matrix(p); q[i] = 100.
    print q[i]
    for k in range(max_iter):
        x = x_solver(Q, r, q, H, d, soft)
        q[i] = (H*(Q*x+r))[0]
        print q[i]
    return q


if __name__ == '__main__':
    pass