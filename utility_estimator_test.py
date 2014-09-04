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
    error, Qs, rs, ws = [], [], [], []
    smooth = [0.0, 0.01, 1.0, 100.0]
    for k in range(len(smooth)):
        Q, r = est.inv_solver(xs, ps, xmax, smooth[k])
        Qs.append(Q)
        rs.append(r)
        ws.append(la.eig(Q)[0])
        u_est, x_est, x_true = U.Utility((Q,r), 'quad'), [], []
        for d in data:
            x_true.append(d[0])
            x_est.append(u_est.compute_demand(d[1]))
        x_true, x_est = matrix(x_true), matrix(x_est)
        error.append(la.norm(x_true-x_est, 1)/la.norm(x_true, 1))
    for Q in Qs: print Q
    for r in rs: print r
    for w in ws: print w
    print error


def test2(N, model, K):
    u, data = sample.sample_utility_data(N, model)
    xs = [d[0] for d in data] # list of demand vectors in equilibrium
    ps = [d[1] for d in data] # list of price vectors
    zs, Hs = est.remove_values(xs, K)
    #xs_i = est.impute_values(zs, Hs)
    smooth = [0.0, 0.01, 1.0, 100.0]
    soft = [10000.0]
    error = []
    for i in range(len(soft)):
        e = []
        for j in range(len(smooth)):
            Q, r = est.mis_solver(zs, Hs, ps, smooth[j], soft[i])
            u_est, x_est, x_true = U.Utility((Q,r), 'quad'), [], []
            for d in data:
                x_true.append(d[0])
                x_est.append(u_est.compute_demand(d[1]))
            x_true, x_est = matrix(x_true), matrix(x_est)
            e.append(la.norm(x_true-x_est, 1)/la.norm(x_true, 1))
        error.append(e)
    print error
            

def main():
    #test1(200, 3)
    test2(200, 2, 1)


if __name__ == '__main__':
    main()