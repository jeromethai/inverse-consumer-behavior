'''
Created on Sep 3, 2014

@author: jeromethai
'''

import numpy as np
from cvxopt import matrix, spdiag, solvers
import rank_nullspace as rn
import Utility as U
from numpy import linalg as la


def compute_xmax(xs):
    """Compute xmax from the demands vectors xs"""
    n = len(xs[0])
    xmax = matrix(0., (n,1))
    for x in xs:
        for i in range(n):
            if x[i]>xmax[i]: xmax[i] = x[i]
    return xmax


def remove_values(data, K):
    """Remove entry K in each demand vector x
    returns the observation matrix
    
    Parameters
    ----------
    data: data points [(x^j, p^j)], x^j demand vector, p^j the price vector
    K: entry to be removed
    """
    xs = [d[0] for d in data] # list of demand vectors in equilibrium
    ps = [d[1] for d in data] # list of price vectors
    n, zs = len(xs[0]), []
    ind = range(K)+range(K+1,n)
    H, k = matrix(0.0, (n-1,n)), 0
    for i in ind: H[k, int(i)] = 1.0; k += 1
    for x in xs:
        z, k = matrix(0.0, (n-1,1)), 0
        for i in ind: z[k] = x[int(i)]; k += 1
        zs.append(z)
    return zip(zs, ps), H


def conic_constraint(n):
    """Construct conic constraint Q<=0 for the SDP"""
    G = [] # construct G
    for i in range(n):
        for j in range(i,n):
            zeros = [0.]*n**2
            zeros[i*n+j] = 1.
            G.append(zeros)
    for i in range(n): G.append([0.]*n**2)
    G = [matrix(G)]
    h = [matrix(0., (n,n))] # construct h
    return G, h


def linear_constraint(x):
    """Construct linear constraint matrix A for the SDP
    A*[vec(Q),r] = Q*x + r"""
    n = len(x)
    A = matrix(0., (n, (n*(n+1))/2))
    for i in range(n):
        for j in range(n):
            if j>=i: A[i,n*i+j-(i*(i+1))/2] = x[j]
            else: A[i,n*j+i-(j*(j+1))/2] = x[j]
    return matrix([[A],[spdiag([1.]*n)]])   


def linear_objective(xs):
    """Construct linear objective for the SDP
    f(Q,r) = sum_j{-(x^j)'*Q*x^j - r'*x^j}
    
    Parameters
    ----------
    xs: list of demand vectors
    """
    n, N = len(xs[0]), len(xs)
    c = matrix(0., ((n*(n+1))/2+n, 1))
    for i in range(n):
        for x in xs: c[(n*(n+1))/2+i] -= x[i]
        for j in range(i,n):
         if j == i:
             for x in xs: c[i*n+j-(i*(i+1))/2] -= x[i]**2
         else:
             for x in xs: c[i*n+j-(i*(i+1))/2] -= 2*x[i]*x[j]
    return c


def x_solver(Q, r, p, H, z, soft):
    """Optimization w.r.t. x-block
    
    Parameters
    ----------
    Q, r: parameters of the quadratic utility function
    p: price vector
    H: observation matrix
    z: observed demand vector
    soft: weight on the observation
    """
    n = len(r)
    P = soft*H.T*H - 2.*Q
    q = p - r - soft*H.T*z
    G = matrix([Q, -spdiag([1.]*n)])
    h = matrix([p-r, matrix(0., (n,1))])
    return solvers.qp(P,q,G,h)['x']


def inv_solver(xs, ps, xmax):
    """Optimization w.r.t. (Q,r)
    
    Parameters
    ----------
    xs: list of (optimal) demand vectors
    ps: list of price vectors
    xmax: maximum demand vector (for each product)
    """
    n, N = len(xs[0]), len(xs)
    # constraint Q <= 0
    G, h = conic_constraint(n)
    # constraint Q*xmax + r >= 0
    A = -linear_constraint(xmax)
    b = matrix(0., (n,1))
    # constraints Q*x^j + r <= p^j
    for j in range(N):
        A = matrix([A, linear_constraint(xs[j])])
        b = matrix([b, matrix(ps[j])])
    # constraint r >= 0
    tmp = matrix([[matrix(0., (n, (n*(n+1))/2))], [-spdiag([1.]*n)]])
    A = matrix([A, tmp])
    b = matrix([b, matrix(0., (n,1))])
    # linear objective
    c = linear_objective(xs)
    x = solvers.sdp(c, A, b, G, h)['x']
    Q, r = matrix(0., (n,n)), matrix(0., (n,1))
    for i in range(n):
        Q[i,i] = x[(n+1)*i-(i*(i+1))/2]
        r[i] = x[(n*(n+1))/2+i]
        for j in range(i+1,n):
            Q[i,j] = x[n*i+j-(i*(i+1))/2]; Q[j,i] = Q[i,j]
    return Q, r
    #return G, h, xmax, A, b, c
    


def mis_solver(zs, H, ps, soft=1e5, max_iter=3):
    """Solves the inverse optimization problem with missing values
    
    Parameters
    ----------
    zs: list of observed demand vectors
    H: observation matrix
    ps: list of price vectors
    soft: weight on the observation
    max_iter: number of iterations
    """
    n = H.size[1]
    Q, r = -spdiag([10.]*n), matrix(10., (n,1))
    for k in range(max_iter):
        xs = [x_solver(Q, r, p, H, z, soft) for p,z in zip(ps,zs)]
        Q, r = inv_solver(xs, ps, compute_xmax(xs))
    return Q, r


def main_solver(data, H=None, soft=1e5, max_iter=3):
    """Main solver that imputes the utility function
    with the best smoothing parameter
    
    Parameters
    ----------
    data: data points [(x^j, p^j)], x^j demand vector, p^j the price vector
    H: observation matrix
    soft: weight on the observation
    max_iter: maximum number of iterations
    """
    zs = [d[0] for d in data] # list of demand vectors in equilibrium
    ps = [d[1] for d in data] # list of price vectors
    if H is None: return inv_solver(zs, ps, compute_xmax(zs))
    else: return mis_solver(zs, H, ps, soft, max_iter)            
            

if __name__ == '__main__':
    pass