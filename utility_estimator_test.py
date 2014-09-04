'''
Created on Sep 3, 2014

@author: jeromethai
'''


from numpy import linalg as la
from cvxopt import matrix
import utility_estimator as est
import generate_samples as sample
import Utility as U


def test1(N, model):
    u, data = sample.sample_utility_data(N, model)
    xs = [d[0] for d in data] # list of demand vectors in equilibrium
    ps = [d[1] for d in data] # list of price vectors
    xmax = est.compute_xmax(xs)
    error = []
    smooth = [0.0, 0.01, 1.0, 100.0, 10000.0]
    for k in range(len(smooth)):
        Q, r = est.inv_solver(xs, ps, xmax, smooth[k])
        u_est, x_est, x_true = U.Utility((Q,r), 'quad'), [], []
        for d in data:
            x_true.append(d[0])
            x_est.append(u_est.compute_demand(d[1]))
        x_true, x_est = matrix(x_true), matrix(x_est)
        error.append(la.norm(x_true-x_est, 1)/la.norm(x_true, 1))
    print error


def test2(N, model, K):
    u, data = sample.sample_utility_data(N, model)
    xs = [d[0] for d in data] # list of demand vectors in equilibrium
    ps = [d[1] for d in data] # list of price vectors
    zs, Hs = est.remove_values(xs, K)
    xs_i = est.impute_values(zs, Hs)
    for x,z,H,xi in zip(xs,zs,Hs,xs_i):
        print x
        print z
        print H
        print xi


def main():
    #test1(200, 4)
    test2(10, 4, 2)


if __name__ == '__main__':
    main()