'''
Created on Sep 3, 2014

@author: jeromethai
'''


import numpy as np
import numpy.random as ra
from cvxopt import matrix, spdiag, solvers
import matplotlib.pylab as plt


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
        
        
    def visualize_AQ(self):
        """Visualize the matrix A if the utility function is sqrt
        Visualize the matrix Q if the utility function is  quad"""
        M = np.matrix(self.data[0])
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_aspect('equal')
        plt.imshow(M, interpolation='nearest', cmap=plt.cm.YlOrRd)
        plt.colorbar()
        plt.show()
        
        
    def sample_data(self, N, pmin=8., pmax=12.):
        """Sample N demands from N prices uniformly sampled in |pmin, pmax|"""
        ps = [matrix(ra.uniform(pmin, pmax, (self.n, 1))) for k in range(N)]
        xs = [self.compute_demand(p) for p in ps]
        return zip(xs,ps)
        

def sqrt_utility(A, b, x):
    """Square root utility function: U(x) = 1'sqrt(Ax+b)"""
    tmp1 = np.sqrt(A*x + b)
    f = sum(tmp1)[0]
    tmp2 = matrix([.5/a[0] for a in tmp1])
    Df = tmp2.T*A
    tmp3 = spdiag([-2.*a**3 for a in tmp2])
    return f, Df, A.T*tmp3*A
    
    
def quad_utility(Q, r, x):
    """Quadratic utility function: U(x) = (1/2)x'Qx + r'x"""
    return (.5*(x.T*Q*x) + r.T*x)[0]


def sample_utility(n, model, alpha, beta, bmax):
    """Sample a sqrt-type utility function following a model
    b_i ~ U[0,bmax]
    'model 1': A = alpha*(I + B) with B_ij ~ U[0,beta]
    'model 2': A_ij = alpha*(beta^|j-i|)
    'model 3': A = alpha*(I + 0.5*B) with B_ij ~ Bernoulli(beta)
    """
    A, b = matrix(0.0, (n,n)), matrix(ra.uniform(0,bmax,(n,1)))
    
    if model == 1: A = matrix(ra.uniform(0,beta,(n,n)))
    
    if model == 2:
        for i in range(n):
            for j in range(n/2):
                A[i, int(np.mod(i+j+1,n))] = beta**(j+1)
                A[i, int(np.mod(i-(j+1),n))] = beta**(j+1)
                
    if model == 3: A = 0.5*matrix(ra.binomial(1,beta,(n,n)))
                
    for i in range(n): A[i,i] = 1.0
    
    return Utility((alpha*A,b), 'sqrt')


if __name__ == '__main__':
    pass