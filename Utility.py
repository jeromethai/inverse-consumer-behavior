'''
Created on Sep 3, 2014

@author: jeromethai
'''


import numpy as np
import numpy.random as ra
from cvxopt import matrix, spdiag, solvers


class Utility:
    """Utility class containing type of the utility function"""
    def __init__(self, data, type):
        self.data = data #parameters of the utility function
        self.type = type #type of the utility function
        if type == 'sqrt' or type == 'quad': self.n = len(data[1])


    def utility(self, x):
        """Computes the utility for a supply x"""
        if self.type == 'sqrt': return sqrt_utility(self.data[0], self.data[1], x)
        if self.type == 'quad': return quad_utility(self.data[0], self.data[1], x)
        if self.type == 'other': pass
        
        
    def compute_demand(self, p):
        """Computes demand x given price p and utility function U
        such that x is solution of min p'x-U(x) s.t. x >=0 """
        
        G, h = spdiag([-1.0]*self.n), matrix(0.0, (self.n, 1))
        
        if self.type == 'quad':
            Q, r = self.data
            return solvers.qp(-Q, p-r, G, h)['x']

        if self.type == 'sqrt':
            def F(x=None, z=None):
                if x is None: return 0, matrix(1.0, (self.n, 1))
                u, Du, H = self.utility(x)
                f, Df  = p.T*x - u, p.T - Du
                if z is None: return f, Df
                return f, Df, -z[0]*H
            return solvers.cp(F, G, h)['x']
        
        
    def sample_data(self, N, pmin=8., pmax=12.):
        """Sample N demands from N prices uniformly sampled in |pmin, pmax|"""
        ps = [ra.uniform(pmin, pmax, self.n) for k in range(N)]
        print ps
        
        

def sqrt_utility(A, b, x):
    """Square root utility function: U(x) = 1'sqrt(Ax+b)"""
    tmp = np.sqrt(A*x + b)
    f = sum(tmp)[0]
    tmp = matrix([.5/x[0] for x in tmp])
    Df = tmp.T*A
    tmp = spdiag([-2.*x**3 for x in tmp])
    return f, Df, A.T*tmp*A
    
    
def quad_utility(Q, r, x):
    """Quadratic utility function: U(x) = (1/2)x'Qx + r'x"""
    return (.5*(x.T*Q*x) + r.T*x)[0]


def sample_utility(n, model, alpha, bmax):
    """Sample a sqrt-type utility function following a model
    b_i ~ U[0,bmax]
    'model 1': A = alpha*(I + 0.5*B) with B_ij ~ U[0,1]
    'model 2': A_ij = alpha*(0.7^|j-i|)
    'model 3': A_ij = alpha*1{i=j} + alpha*0.4*1{|i-j|=1} + alpha*0.2*1{|i-j|=2} + alpha*0.1*1{|i-j|=3}
    'model 4': A = alpha*(I + 0.5*B) with B_ij ~ Bernoulli(0.1)
    """
    
    A, b = matrix(0.0, (n,n)), matrix(ra.uniform(0,bmax,(n,1)))
    
    if model == 'model 1': A = matrix(ra.uniform(0,0.5,(n,n)))
    
    if model == 'model 2':
        for i in range(n):
            for j in range(i+1,n):
                A[i,j] = .7**abs(j-i); A[j,i] = A[i,j]

    if model == 'model 3':
        for i in range(n):
            for j in range(i+1,n):
                if abs(j-i)<4: A[i,j] = .8*.5**abs(j-i); A[j,i] = A[i,j]
                
    if model == 'model 4': A = 0.5*matrix(ra.binomial(1,0.1,(n,n)))
                
    for i in range(n): A[i,i] = 1.0
                
    print A
    return Utility((alpha*A,b), 'sqrt')


if __name__ == '__main__':
    pass