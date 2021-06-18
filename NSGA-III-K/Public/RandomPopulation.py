"""
Random population.
"""

import numpy as np

class population:
    """Creates a class for population"""
    def __init__(Pop, dec, obj):
        Pop.dec = dec
        Pop.obj = obj

def randomPopulation(N,n,m,lb,ub,evaluate):
    """Generates a random population"""
    decision = lb+(np.random.rand(N,n)*(ub-lb))
    objective = np.zeros([N,m])
    for i in range(0,N):
        objective[i] = evaluate(decision[i],m)
    return population(decision,objective)
