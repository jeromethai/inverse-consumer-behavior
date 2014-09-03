'''
Created on Sep 3, 2014

@author: jeromethai
'''

import Utility as U
from cvxopt import matrix
import numpy.random as ra


def test1(n, alpha, bmax, model):
    #Q = matrix([[-1., 0.],[0., -1.]])
    #r = matrix([1., 1.])
    #u1 = U.Utility((Q,r), 'quad')
    #print u1.type, u1.n, u1.utility(matrix([1., 1.]))
    u2 = U.sample_utility(n, model, alpha, bmax)
    #f, Df, H = u2.utility(matrix([1.]*n))
    #print f
    #print Df
    #print H
    #print u1.compute_demand(matrix([.5, .5]))
    x = u2.compute_demand(matrix(ra.uniform(8.,12.,(n,1))))
    print x
    print sum(x)
    #print u2.sample_data(10)


def test2():
    pass


def main():
    #test1(5, 30, 1, 'model 1')
    #test1(5, 25, 1, 'model 2')
    #test1(5, 35, 1, 'model 3')
    test1(5, 55, 1, 'model 4')


if __name__ == '__main__':
    main()