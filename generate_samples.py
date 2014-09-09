'''
Created on Sep 3, 2014

@author: jeromethai
'''

import Utility as U
from cvxopt import matrix
import numpy.random as ra


def sample_utility_data(N, model, n=5, alpha=50., beta=0.3, bmax=1):
    """Sample a utility function
    then sample N data points (x^j,p^j)
    with 5 products (n=5)
    """
    u = U.sample_utility(n, model, alpha, beta, bmax)
    return u, u.sample_data(N)


def test(N, model):
    u, data = sample_utility_data(N, model)
    for d in data:
        print sum(d[0])
        print matrix([[d[0]],[d[1]]])
    print u.data[0]


def main():
    test(5, 1)


if __name__ == '__main__':
    main()