"""
DTLZ7_MINUS test problem.

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
    k = 20
    n = m-1+k
    lb = np.zeros(n)
    ub = np.ones(n)
    return n,lb,ub

def evaluate(individual,m):
    """Evaluates an individual for the DTLZ7_MINUS test problem"""
    evaluation = np.zeros(m)
    n = len(individual)
    s = 0
    for i in range(m-1,n):
        s += individual[i]
    g = 1+9/(n-m+1)*s
    s = 0
    for i in range(0,m-1):
        evaluation[i] = individual[i]
        s += evaluation[i]/(1+g)*(1+np.sin(3*np.pi*evaluation[i]))
    h = m-s
    evaluation[m-1] = (1+g)*h
    return -evaluation
