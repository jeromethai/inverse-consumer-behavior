'''
Created on Sep 3, 2014

@author: jeromethai
'''


import numpy as np
from numpy import linalg as la
from cvxopt import matrix
import utility_estimator as est
import generate_samples as sample
import Utility as U
import matplotlib.pyplot as plt


def display_results(u1, u2, u_true, error1, error2, xmax):
    """Display results
    
    Parameters
    ----------
    u1: estimated utility function 1
    u2: estimated utility function 2
    u_true: true utility function
    error1: relative error in estimated demand from u1
    error2: relative error in estimated demand from u2
    xmax: maximum prices vector
    """
    index = [0,2]
    for i in index:
        xdata = np.linspace(0.0, xmax[i], num=100)
        e = matrix([0.]*5); e[i] = 1.0
        vals1 = [u1.utility(matrix(x*e)) for x in xdata]
        vals2 = [u2.utility(matrix(x*e)) for x in xdata]
        true_vals = [u_true.utility(matrix(x*e))[0] for x in xdata]
        offset1 = np.mean(true_vals) - np.mean(vals1)
        offset2 = np.mean(true_vals) - np.mean(vals2)
        plt.plot( xdata, [offset1+v for v in vals1], 'b', label='full')
        plt.plot( xdata, [offset2+v for v in vals2], 'g', label='missing')
        plt.plot( xdata, true_vals, 'r', label='true')
        plt.xlabel('Demand in product {}'.format(i+1))
        plt.ylabel('Utility')
        plt.title('Estimated utility function, err1={:.0f}%, err2={:.0f}%'.format(100*error1, 100*error2))
        plt.legend(loc=0)
        plt.show()


def test(N, model, K):
    """Test utility_estimator
    
    Parameters:
    -----------
    N: number of samples
    model: utility function model
    K: missing value is demand in product K
    """
    u, data = sample.sample_utility_data(N, model)
    data_mis, H = est.remove_values(data, K)
    Q1, r1, smooth1 = est.main_solver(data)
    Q2, r2, smooth2 = est.main_solver(data_mis, H)
    u1 = U.Utility((Q1,r1), 'quad')
    u2 = U.Utility((Q2,r2), 'quad')
    xs = [d[0] for d in data]
    x_true = matrix(xs)
    x1 = matrix([u1.compute_demand(d[1]) for d in data])
    x2 = matrix([u2.compute_demand(d[1]) for d in data])
    e1 = la.norm(x1-x_true, 1)/la.norm(x_true, 1)
    e2 = la.norm(x2-x_true, 1)/la.norm(x_true, 1)
    print e1,e2
    xmax = est.compute_xmax(xs)
    display_results(u1, u2, u, e1, e2, xmax)            


def main():
    test(200, 2, 0)


if __name__ == '__main__':
    main()