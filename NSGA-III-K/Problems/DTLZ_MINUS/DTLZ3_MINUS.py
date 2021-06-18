"""
DTLZ3_MINUS test problem.

K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, "Scalable Test Problems for 
Evolutionary Multiobjective Optimization," in Evolutionary Multiobjective 
Optimization. Theoretical Advances and Applications, pp. 105-145, 2005.

H. Ishibuchi, Y. Setoguchi, H. Masuda, and Y. Nojima, "Performance of 
Decomposition-Based Many-Objective Algorithms Strongly Depends on Pareto Front 
Shapes," in IEEE Transactions on Evolutionary Computation, vol. 21, no. 2, 
pp. 169-190, 2017.
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
    """Evaluates an individual for the DTLZ3_MINUS test problem"""
    evaluation = np.zeros(m)
    n = len(individual)
    s = 0
    for i in range(m-1,n):
        s += ((individual[i]-0.5)**2-np.cos(20*np.pi*(individual[i]-0.5)))
    g = 100*(n-m+1+s)
    for i in range(0,m):
        mult = 1
        for j in range(0,m-i-1):
            mult *= np.cos(individual[j]*np.pi/2)
        if (i == 0):
            evaluation[i] = (1+g)*mult
        else:
            evaluation[i] = (1+g)*mult*np.sin(individual[m-i-1]*np.pi/2)
    return -evaluation
