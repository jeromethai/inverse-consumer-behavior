'''
Created on Sep 3, 2014

@author: jeromethai
'''


import numpy as np
from cvxopt import matrix
import utility_estimator as est
import generate_samples as sample
import Utility as U


def test1():
    G, h, xmax, A, b, c = est.inv_solver([[1,1,1,1,1],[1,2,3,4,5]],[[0,0,0,0,0],[1,2,3,4,5]])
    print G[0]
    print h[0]
    print xmax
    print A.size
    print A
    print b
    print c
    #for g in G: print g


def test2(N, model):
    u, data = sample.sample_utility_data(N, model)
    Q, r = est.inv_solver(data)
    u_est, x_est, x_true = U.Utility((Q,r), 'quad'), [], []
    for d in data:
        x_true.append(d[1])
        x_est.append(u_est.compute_demand(d[1]))
    x_true, x_est = matrix(x_true), matrix(x_est)
    print np.linalg.norm(x_true-x_est, 1)/np.linalg.norm(x_true, 1)


def main():
    test2(200, 4)


if __name__ == '__main__':
    main()