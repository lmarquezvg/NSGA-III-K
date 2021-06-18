"""
Reference set generated with the systematic approach proposed by Das and Dennis.

I. Das and J. E. Dennis, "Normal-Boundary Intersection: A New Method for 
Generating the Pareto Surface in Nonlinear Multicriteria Optimization Problems," 
in SIAM Journal on Optimization, vol. 8, no. 3, pp. 631-657, 1998.

Y. Tian, R. Cheng, X. Zhang, and Y. Jin, "PlatEMO: A MATLAB platform for 
evolutionary multi-objective optimization," in IEEE Computational Intelligence 
Magazine, vol. 12, no. 4, pp. 73-87, 2017.
"""

import numpy as np
import numpy.matlib
import itertools
import math

def nchoosek(v, k):
    """Returns a matrix with all combinations of v taken k at a time"""
    return np.array(list(itertools.combinations(v, k)))

def uniformPoints(H1, H2, m):
    """Generates a two layers of reference points with H1 and H2 divisions respectively"""
    W1 = nchoosek(list(range(1, H1+m)), m-1)-np.matlib.repmat(np.arange(0, m-1), math.comb(H1+m-1, m-1), 1)-1
    W1 = (np.hstack((W1, np.zeros([len(W1), 1])+H1))-np.hstack((np.zeros([len(W1), 1]), W1)))/H1
    if (H2 > 0):
        W2 = nchoosek(list(range(1, H2+m)), m-1)-np.matlib.repmat(np.arange(0, m-1), math.comb(H2+m-1, m-1), 1)-1
        W2 = (np.hstack((W2, np.zeros([len(W2), 1])+H2))-np.hstack((np.zeros([len(W2), 1]), W2)))/H2
        W2 = W2/2+1/(2*m)
        W = np.vstack((W1, W2))
    else:
        W = W1
    W[W<1e-6] = 1e-6
    return W
