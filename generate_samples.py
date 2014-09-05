'''
Created on Sep 3, 2014

@author: jeromethai
'''

import Utility as U
from cvxopt import matrix
import numpy.random as ra


def sample_utility_data(N, model, n=5, bmax=1):
    """Sample a utility function
    then sample N data points (x^j,p^j)
    with 5 products (n=5)
    """
    alpha = [30.,25.,35.,55.]
    u = U.sample_utility(n, model, alpha[model-1], bmax)
    return u, u.sample_data(N)


def test(N, model):
    u, data = sample_utility_data(N, model)
    for d in data:
        print sum(d[0])
        print matrix([[d[0]],[d[1]]])


def main():
    test(5, 4)


if __name__ == '__main__':
    main()