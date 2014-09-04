'''
Created on Sep 3, 2014

@author: jeromethai
'''

from cvxopt import matrix, spdiag, solvers


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
    sum_j{-(x^j)'*Q*x^j - r'*x^j}"""
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


def x_solver(Q, r, p, soft, obs, x_obs):
    """Optimization w.r.t. x-block
    
    Parameters
    ----------
    Q: 
    r:
    p:
    soft:
    obs: list of observed indinces vectors
    x_obs: list of observed demand vectors
    """
    


def inv_solver(data):
    """Optimization w.r.t. (Q,r)
    data = [(x^j,p^j)] list of tuples of demand vector and associated price vector
    """
    xs = [d[0] for d in data] # list of demand vectors in equilibrium
    ps = [d[1] for d in data] # list of price vectors
    n, N = len(xs[0]), len(xs)
    # constraint Q <= 0
    G, h = conic_constraint(n)
    # constraint Q*xmax + r >= 0
    xmax = matrix(0., (n,1)) # compute xmax
    for x in xs:
        for i in range(n):
            if x[i]>xmax[i]: xmax[i] = x[i]
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
    


def mis_solver():
    pass


if __name__ == '__main__':
    pass