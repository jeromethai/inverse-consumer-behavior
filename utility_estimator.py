'''
Created on Sep 3, 2014

@author: jeromethai
'''

import numpy as np
from cvxopt import matrix, spdiag, solvers


def compute_xmax(xs):
    """Compute xmax from the demands vectors xs"""
    n = len(xs[0])
    xmax = matrix(0., (n,1))
    for x in xs:
        for i in range(n):
            if x[i]>xmax[i]: xmax[i] = x[i]
    return xmax


def remove_values(xs, K):
    """Remove K entries uniformly in each demand vector x"""
    n, zs, Hs = len(xs[0]), [], []
    for x in xs:
        ind = np.sort(np.random.permutation(n)[:n-K])
        z, H = matrix(0.0, (n-K,1)), matrix(0.0, (n-K,n))
        k = 0
        for i in ind:
            z[k] = x[int(i)]
            H[k, int(i)] = 1.0
            k += 1
        zs.append(z)
        Hs.append(H)
    return zs, Hs


def impute_values(zs, Hs):
    """Impute the missing demands for a product
    with the mean of observed values
    from the list on observed demands zs
    and the list of observation matrices Hs"""
    n = Hs[0].size[1]
    xmean, counts = [0.]*n, [0]*n
    for z, H in zip(zs, Hs):
        k = 0
        for i in range(n):
            if sum(H[:,i]) > 0:
                xmean[i] += z[k]
                counts[i] += 1
                k += 1
    for i in range(n): xmean[i] /= counts[i]
    xs = []
    for z, H in zip(zs, Hs):
        x = matrix(0.0, (n,1))
        k = 0
        for i in range(n):
            if sum(H[:,i]) > 0:
                x[i] = z[k]
                k += 1
            else:
                x[i] = xmean[i]
        xs.append(x)
    return xs 


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


def linear_objective(xs, smooth=0.0):
    """Construct linear objective for the SDP
    f(Q,r) = sum_j{-(x^j)'*Q*x^j - r'*x^j} + smooth*sum_i{r_i}"""
    n, N = len(xs[0]), len(xs)
    c = matrix(0., ((n*(n+1))/2+n, 1))
    for i in range(n):
        for x in xs: c[(n*(n+1))/2+i] -= x[i] + smooth
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
    G = matrix([Q, -spdiag([0.]*n)])
    h = matrix([p-r, matrix(0., (n,1))])
    return solvers.qp(P,q,G,h)['x']


def inv_solver(xs, ps, xmax, smooth=0.0):
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
    c = linear_objective(xs, smooth)
    x = solvers.sdp(c, A, b, G, h)['x']
    Q, r = matrix(0., (n,n)), matrix(0., (n,1))
    for i in range(n):
        Q[i,i] = x[(n+1)*i-(i*(i+1))/2]
        r[i] = x[(n*(n+1))/2+i]
        for j in range(i+1,n):
            Q[i,j] = x[n*i+j-(i*(i+1))/2]; Q[j,i] = Q[i,j]
    return Q, r
    #return G, h, xmax, A, b, c
    


def mis_solver(zs, Hs, ps, smooth, soft, xs0=None, max_iter=3):
    if xs0 is None: xs0 = est.impute_values(zs, Hs)
    pass


if __name__ == '__main__':
    pass