"""
DTLZ4 test problem.

K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, "Scalable Test Problems for 
Evolutionary Multiobjective Optimization," in Evolutionary Multiobjective 
Optimization. Theoretical Advances and Applications, pp. 105-145, 2005.
"""

import numpy as np

def parameters(m):
    """Returns number of decision variables, lower bounds, and upper bounds"""
    k = 10
    n = m-1+k
    lb = np.zeros(n)
    ub = np.ones(n)
    return n,lb,ub

def evaluate(individual,m):
    """Evaluates an individual for the DTLZ4 test problem"""
    evaluation = np.zeros(m)
    n = len(individual)
    s = 0
    for i in range(m-1,n):
        s += (individual[i]-0.5)**2
    g = s
    for i in range(0,m):
        mult = 1
        for j in range(0,m-i-1):
            mult *= np.cos((individual[j]**100)*np.pi/2)
        if (i == 0):
            evaluation[i] = (1+g)*mult
        else:
            evaluation[i] = (1+g)*mult*np.sin((individual[m-i-1]**100)*np.pi/2)
    return evaluation
