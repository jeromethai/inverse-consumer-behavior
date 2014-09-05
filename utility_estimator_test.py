'''
Created on Sep 3, 2014

@author: jeromethai
'''


from numpy import linalg as la
from cvxopt import matrix
import utility_estimator as est
import generate_samples as sample
import Utility as U


def test(N, model, K):
    """Test utility_estimator"""
    u, data = sample.sample_utility_data(N, model)
    data_mis, H = est.remove_values(data, K)
    Q1, r1, smooth1 = est.main_solver(data)
    Q2, r2, smooth2 = est.main_solver(data_mis, H)
    u1 = U.Utility((Q1,r1), 'quad')
    u2 = U.Utility((Q2,r2), 'quad')
    x_true = matrix([d[0] for d in data])
    x1 = matrix([u1.compute_demand(d[1]) for d in data])
    x2 = matrix([u2.compute_demand(d[1]) for d in data])
    print smooth1, smooth2
    print la.norm(x1-x_true, 1)/la.norm(x_true, 1)
    print la.norm(x2-x_true, 1)/la.norm(x_true, 1)
    print H
            

def main():
    test(200, 2, 0)


if __name__ == '__main__':
    main()